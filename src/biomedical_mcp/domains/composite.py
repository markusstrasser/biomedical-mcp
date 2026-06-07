"""Composite domain — multi-API compound queries for variant review and gene investigation.

These tools call multiple underlying APIs and return structured, provenance-tagged results.
Per-enrichment timeout pattern from BioMCP: each section has an independent timeout
so a slow/failed API doesn't block the response.
"""

import logging
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeout

from fastmcp import FastMCP

from biomedical_mcp.cache import Cache
from biomedical_mcp.clingen import ClinGen
from biomedical_mcp.gnomad import GnomAD
from biomedical_mcp.panelapp import PanelApp
from biomedical_mcp.hpo import HPO
from biomedical_mcp.gtex import GTEx
from biomedical_mcp.hgnc import HGNC
from biomedical_mcp.myvariant import MyVariant
from biomedical_mcp.litvar import LitVar
from biomedical_mcp.orphanet import Orphanet
from biomedical_mcp.gwas_catalog import GWASCatalog
from biomedical_mcp.provenance import make_provenance
from biomedical_mcp.variant_input import normalize_variant
from biomedical_mcp.composite_core import run_sections, describe
from biomedical_mcp import composite_sections as cs

log = logging.getLogger(__name__)

# Thread pool for concurrent API calls (APIs use sync httpx)
_executor = ThreadPoolExecutor(max_workers=8)

# Per-enrichment timeout (seconds) — optional sections that miss this deadline
# return a warning instead of blocking the response. Pattern from BioMCP.
ENRICHMENT_TIMEOUT = 8


def _gather_sync(*callables, timeout: int = ENRICHMENT_TIMEOUT):
    """Run multiple sync functions concurrently via threads, return results in order.

    Each callable gets its own timeout — a slow API doesn't block the rest.
    Failed/timed-out calls return an error dict instead of raising.
    """
    futures = [_executor.submit(fn) for fn in callables]
    results = []
    for future in futures:
        try:
            results.append(future.result(timeout=timeout))
        except FuturesTimeout:
            results.append({"error": "timeout", "detail": f"Enrichment timed out after {timeout}s"})
        except Exception as exc:
            results.append({"error": f"{type(exc).__name__}: {exc}"})
    return results


def create_server(cache: Cache) -> FastMCP:
    clingen = ClinGen(cache)
    gnomad = GnomAD(cache)
    panelapp = PanelApp(cache)
    hpo = HPO(cache)
    gtex = GTEx(cache)
    hgnc = HGNC(cache)
    myvariant = MyVariant(cache)
    litvar = LitVar(cache)
    orphanet = Orphanet(cache)
    gwas = GWASCatalog(cache)

    server = FastMCP("composite")

    @server.tool(tags={"variant-review"})
    def variant_context(variant_id: str, gene_symbol: str | None = None) -> dict:
        """Comprehensive variant context in one call: population frequency, ClinVar,
        predictions, CGI drug associations, COSMIC, GWAS, gene constraint, expression,
        panels, literature.

        Accepts multiple formats: rsID ("rs1801133"), HGVS ("chr7:g.140453136A>T"),
        or gene+protein change ("BRAF V600E", "BRAF p.Val600Glu").

        Args:
            variant_id: Variant identifier in any supported format.
            gene_symbol: Gene symbol (optional, auto-detected from results if omitted).
        """
        # Normalize variant input — handles "BRAF V600E", "p.Val600Glu", etc.
        parsed = normalize_variant(variant_id)
        if parsed.format == "search_hint":
            return {
                "variant_id": variant_id,
                "status": "needs_search",
                "parsed_format": parsed.format,
                "suggestion": parsed.suggestion,
                "gene": parsed.gene,
                "protein_change": parsed.protein_change,
            }

        lookup_id = parsed.normalized
        if parsed.gene and not gene_symbol:
            gene_symbol = parsed.gene

        # Phase 1: variant-level lookups (concurrent, per-enrichment timeout)
        results = _gather_sync(
            lambda: myvariant.lookup(lookup_id),
            lambda: litvar.variant_publications(lookup_id, limit=5),
            lambda: gnomad.variant_frequency(lookup_id),
            lambda: gwas.variant_associations(lookup_id),
        )
        mv_result, lit_result, gn_result, gwas_result = results

        # MyVariant returns a list for multi-allelic rsIDs — take first entry
        if isinstance(mv_result, list) and mv_result:
            mv_result = mv_result[0]

        # Extract gene from MyVariant if not provided
        if not gene_symbol and isinstance(mv_result, dict):
            gene_symbol = (mv_result.get("predictions", {}) or {}).get("gene")
            if not gene_symbol:
                cv = mv_result.get("clinvar", {})
                if isinstance(cv, dict):
                    gene_symbol = cv.get("gene")

        # Phase 2: gene-level lookups (concurrent, if gene known)
        gene_data = {}
        if gene_symbol:
            gene_results = _gather_sync(
                lambda: gnomad.gene_constraint(gene_symbol),
                lambda: panelapp.gene_panels(gene_symbol, confidence="3"),
                lambda: gtex.top_tissues(gene_symbol, limit=5),
            )
            constraint_result, panels_result, expr_result = gene_results
            gene_data = {
                "symbol": gene_symbol,
                "constraint": constraint_result,
                "panels": panels_result,
                "top_tissues": expr_result,
            }

        result = {
            "variant_id": variant_id,
            "normalized_id": lookup_id,
            "myvariant": mv_result,
            "gnomad": gn_result,
            "literature": lit_result,
            "gwas": gwas_result,
            "gene": gene_data,
            "provenance": {
                "myvariant": make_provenance("myvariant", evidence_grade="C4"),
                "gnomad": make_provenance("gnomad", data_version="gnomad_r4", evidence_grade="E5"),
                "litvar": make_provenance("litvar", evidence_grade="D4"),
                "gwas_catalog": make_provenance("gwas_catalog", evidence_grade="C3"),
                "panelapp": make_provenance("panelapp"),
                "gtex": make_provenance("gtex"),
            },
        }

        # Include CGI/COSMIC if present in MyVariant results
        if isinstance(mv_result, dict):
            if "cgi" in mv_result:
                result["cgi_drug_associations"] = mv_result["cgi"]
            if "cosmic" in mv_result:
                result["cosmic"] = mv_result["cosmic"]

        return result

    # Client registry passed to the section-fetcher builders. The client classes
    # ARE the adapter layer; the registry binds sections to their methods.
    clients = {
        "hgnc": hgnc, "gnomad": gnomad, "panelapp": panelapp, "gtex": gtex,
        "hpo": hpo, "clingen": clingen, "orphanet": orphanet,
        "myvariant": myvariant, "litvar": litvar, "gwas_catalog": gwas,
    }

    @server.tool(tags={"gene-lookup"})
    def gene_dossier(gene_symbol: str, sections: list[str] | None = None) -> dict:
        """Comprehensive gene report via a standardized partial-failure envelope.

        Sections (default: all): nomenclature, constraint, panels, expression,
        phenotypes, gene_disease_validity, dosage_sensitivity, rare_diseases.
        Call describe_sections("gene") to see sources and descriptions.

        One dead upstream API degrades only its section (status="error"), never the
        whole dossier. Pass `sections` to fetch a subset.

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "CYP2D6", "SCN1A").
            sections: Subset of section names, or omit for all default sections.
        """
        fetchers = cs.gene_fetchers(gene_symbol, clients)
        return run_sections("gene", gene_symbol, sections, fetchers, cs.SECTIONS["gene"])

    @server.tool(tags={"discovery"})
    def describe_sections(entity: str) -> dict:
        """List the available sections for a composite entity, with sources and
        descriptions. Use this to discover valid `sections` values before calling
        gene_dossier / variant_context / drug_profile / protein_profile / disease_profile.

        Args:
            entity: One of "gene", "variant", "drug", "protein", "disease".
        """
        meta = cs.SECTIONS.get(entity)
        if meta is None:
            return {
                "error": f"unknown entity '{entity}'",
                "valid_entities": sorted(cs.SECTIONS.keys()),
            }
        return describe(entity, meta)

    @server.tool(tags={"variant-review"})
    def batch_variant_lookup(variant_ids: str) -> dict:
        """Look up multiple variants in one call. Comma-separated IDs (max 10).

        Accepts any supported format per variant: rsIDs, HGVS, "BRAF V600E".

        Args:
            variant_ids: Comma-separated variant IDs (e.g. "rs1801133,rs1801131,BRAF V600E").
        """
        ids = [v.strip() for v in variant_ids.split(",") if v.strip()][:10]
        parsed = [normalize_variant(v) for v in ids]

        # Separate exact lookups from search hints
        exact_ids = [p.normalized for p in parsed if p.format != "search_hint"]
        hints = [p for p in parsed if p.format == "search_hint"]

        results = {}
        if exact_ids:
            batch = myvariant.batch(exact_ids)
            for item in batch:
                vid = item.get("variant_id", "unknown")
                results[vid] = item

        if hints:
            results["needs_search"] = [
                {"input": h.raw, "suggestion": h.suggestion} for h in hints
            ]

        return {
            "count": len(results),
            "results": results,
            "provenance": make_provenance("myvariant", evidence_grade="C4"),
        }

    @server.tool(tags={"admin"})
    def health_check() -> dict:
        """Check connectivity to all upstream APIs. Returns status per source."""
        sources = {
            "myvariant": myvariant,
            "gnomad": gnomad,
            "gtex": gtex,
            "panelapp": panelapp,
            "hpo": hpo,
            "hgnc": hgnc,
            "clingen": clingen,
            "litvar": litvar,
            "orphanet": orphanet,
            "gwas_catalog": gwas,
        }
        statuses = {}
        for name, client in sources.items():
            try:
                ok = client.validate()
                statuses[name] = "ok" if ok else "error"
            except Exception as exc:
                statuses[name] = f"error: {type(exc).__name__}"

        ok_count = sum(1 for v in statuses.values() if v == "ok")
        return {
            "status": f"{ok_count}/{len(statuses)} ok",
            "sources": statuses,
        }

    return server
