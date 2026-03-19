"""Composite domain — multi-API compound queries for variant review and gene investigation.

These tools call multiple underlying APIs and return structured, provenance-tagged results.
"""

import asyncio
import logging
from concurrent.futures import ThreadPoolExecutor

from fastmcp import FastMCP

from biomedical_mcp.cache import Cache
from biomedical_mcp.gnomad import GnomAD
from biomedical_mcp.panelapp import PanelApp
from biomedical_mcp.hpo import HPO
from biomedical_mcp.gtex import GTEx
from biomedical_mcp.hgnc import HGNC
from biomedical_mcp.myvariant import MyVariant
from biomedical_mcp.litvar import LitVar
from biomedical_mcp.provenance import make_provenance

log = logging.getLogger(__name__)

# Thread pool for concurrent API calls (APIs use sync httpx)
_executor = ThreadPoolExecutor(max_workers=6)


def _gather_sync(*callables):
    """Run multiple sync functions concurrently via threads, return results in order."""
    loop = asyncio.new_event_loop()
    try:
        futures = [loop.run_in_executor(_executor, fn) for fn in callables]
        results = loop.run_until_complete(asyncio.gather(*futures, return_exceptions=True))
    finally:
        loop.close()
    return results


def create_server(cache: Cache) -> FastMCP:
    gnomad = GnomAD(cache)
    panelapp = PanelApp(cache)
    hpo = HPO(cache)
    gtex = GTEx(cache)
    hgnc = HGNC(cache)
    myvariant = MyVariant(cache)
    litvar = LitVar(cache)

    server = FastMCP("composite")

    @server.tool(tags={"variant-review"})
    def variant_context(variant_id: str, gene_symbol: str | None = None) -> dict:
        """Comprehensive variant context in one call: population frequency, ClinVar,
        predictions, gene constraint, expression, panels, literature.

        Calls gnomAD + MyVariant + GTEx + PanelApp + LitVar2 concurrently.

        Args:
            variant_id: rsID (e.g. "rs1801133"), gnomAD ID (e.g. "1-11796321-G-A"), or HGVS.
            gene_symbol: Gene symbol (optional, enriches with constraint/expression/panels).
        """
        # Phase 1: variant-level lookups (concurrent)
        results = _gather_sync(
            lambda: myvariant.lookup(variant_id),
            lambda: litvar.variant_publications(variant_id, limit=5),
            lambda: gnomad.variant_frequency(variant_id) if not variant_id.startswith("rs") else {"note": "Use gnomAD-style ID for frequency"},
        )
        mv_result, lit_result, gn_result = results

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
                "constraint": constraint_result if not isinstance(constraint_result, Exception) else {"error": str(constraint_result)},
                "panels": panels_result if not isinstance(panels_result, Exception) else {"error": str(panels_result)},
                "top_tissues": expr_result if not isinstance(expr_result, Exception) else {"error": str(expr_result)},
            }

        return {
            "variant_id": variant_id,
            "myvariant": mv_result if not isinstance(mv_result, Exception) else {"error": str(mv_result)},
            "gnomad": gn_result if not isinstance(gn_result, Exception) else {"error": str(gn_result)},
            "literature": lit_result if not isinstance(lit_result, Exception) else {"error": str(lit_result)},
            "gene": gene_data,
            "provenance": {
                "myvariant": make_provenance("myvariant", evidence_grade="C4"),
                "gnomad": make_provenance("gnomad", data_version="gnomad_r4", evidence_grade="E5"),
                "litvar": make_provenance("litvar", evidence_grade="D4"),
                "panelapp": make_provenance("panelapp"),
                "gtex": make_provenance("gtex"),
            },
        }

    @server.tool(tags={"gene-lookup"})
    def gene_dossier(gene_symbol: str) -> dict:
        """Comprehensive gene report: names, constraint, pathways, diseases, expression, panels.

        Calls HGNC + gnomAD + PanelApp + GTEx + HPO concurrently.

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "CYP2D6", "SCN1A").
        """
        results = _gather_sync(
            lambda: hgnc.gene_names(gene_symbol),
            lambda: gnomad.gene_constraint(gene_symbol),
            lambda: panelapp.gene_panels(gene_symbol, confidence="all"),
            lambda: gtex.top_tissues(gene_symbol, limit=10),
            lambda: hpo.gene_phenotypes(gene_symbol),
        )
        names_result, constraint_result, panels_result, expr_result, pheno_result = results

        def safe(result):
            return result if not isinstance(result, Exception) else {"error": str(result)}

        return {
            "gene_symbol": gene_symbol,
            "nomenclature": safe(names_result),
            "constraint": safe(constraint_result),
            "panels": safe(panels_result),
            "expression": safe(expr_result),
            "phenotypes": safe(pheno_result),
            "provenance": {
                "hgnc": make_provenance("hgnc"),
                "gnomad": make_provenance("gnomad", data_version="gnomad_r4"),
                "panelapp": make_provenance("panelapp"),
                "gtex": make_provenance("gtex"),
                "hpo": make_provenance("hpo"),
            },
        }

    return server
