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

    @server.tool(tags={"gene-lookup"})
    def gene_dossier(gene_symbol: str) -> dict:
        """Comprehensive gene report: names, constraint, pathways, diseases, expression,
        panels, gene-disease validity, dosage sensitivity, rare disease associations.

        Calls HGNC + gnomAD + PanelApp + GTEx + HPO + ClinGen + Orphanet concurrently.

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "CYP2D6", "SCN1A").
        """
        results = _gather_sync(
            lambda: hgnc.gene_names(gene_symbol),
            lambda: gnomad.gene_constraint(gene_symbol),
            lambda: panelapp.gene_panels(gene_symbol, confidence="all"),
            lambda: gtex.top_tissues(gene_symbol, limit=10),
            lambda: hpo.gene_phenotypes(gene_symbol),
            lambda: clingen.gene_validity(gene_symbol),
            lambda: clingen.gene_dosage(gene_symbol),
            lambda: orphanet.gene_diseases(gene_symbol),
        )
        (names_result, constraint_result, panels_result, expr_result,
         pheno_result, validity_result, dosage_result, orphanet_result) = results

        def safe(result):
            return result if not isinstance(result, Exception) else {"error": str(result)}

        return {
            "gene_symbol": gene_symbol,
            "nomenclature": safe(names_result),
            "constraint": safe(constraint_result),
            "panels": safe(panels_result),
            "expression": safe(expr_result),
            "phenotypes": safe(pheno_result),
            "gene_disease_validity": safe(validity_result),
            "dosage_sensitivity": safe(dosage_result),
            "rare_diseases": safe(orphanet_result),
            "provenance": {
                "hgnc": make_provenance("hgnc"),
                "gnomad": make_provenance("gnomad", data_version="gnomad_r4"),
                "panelapp": make_provenance("panelapp"),
                "gtex": make_provenance("gtex"),
                "hpo": make_provenance("hpo"),
                "clingen": make_provenance("clingen_gene_validity", evidence_grade="A1"),
                "orphanet": make_provenance("orphanet", evidence_grade="B2"),
            },
        }

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
