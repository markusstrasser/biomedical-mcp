"""Panels domain — curated gene panels by condition (PanelApp, Genomics England)."""

from fastmcp import FastMCP

from biomedical_mcp.panelapp import PanelApp
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    panelapp = PanelApp(cache)

    server = FastMCP("panels")

    @server.tool(tags={"variant-review", "gene-lookup"})
    def gene_panels(gene_symbol: str, confidence: str = "3") -> dict:
        """Which clinical panels include this gene? Default: Green (high confidence) only.

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "KCNQ1").
            confidence: Minimum confidence: "3"=Green (default), "2"=Amber, "1"=Red, "all"=any.
        """
        return panelapp.gene_panels(gene_symbol, confidence=confidence)

    @server.tool(tags={"gene-lookup"})
    def panel_genes(panel_id: str, signed_off: bool = True, confidence: str = "3") -> dict:
        """All genes in a panel. Default: signed-off version, Green confidence only.

        Args:
            panel_id: PanelApp panel ID (numeric string, e.g. "55" for Cardiomyopathy).
            signed_off: Use signed-off (reviewed) version (default True).
            confidence: Minimum confidence: "3"=Green, "2"=Amber, "1"=Red, "all"=any.
        """
        return panelapp.panel_genes(panel_id, signed_off=signed_off, confidence=confidence)

    @server.tool(tags={"gene-lookup"})
    def panel_search(query: str, limit: int = 20) -> dict:
        """Search panels by disease keyword.

        Args:
            query: Disease/condition keyword (e.g. "cardiomyopathy", "epilepsy").
            limit: Max results (default 20).
        """
        return panelapp.panel_search(query, limit=limit)

    return server
