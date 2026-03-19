"""Pathways domain — metabolic, signaling, biological pathways."""

from fastmcp import FastMCP

from biomedical_mcp.kegg import KEGG
from biomedical_mcp.reactome import Reactome
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    kegg = KEGG(cache)
    reactome = Reactome(cache)

    server = FastMCP("pathways")

    @server.tool()
    def kegg_gene(gene_symbol: str, organism: str = "hsa") -> dict:
        """Get KEGG pathways for a gene.

        Args:
            gene_symbol: Gene symbol (e.g. "TP53", "EGFR").
            organism: KEGG organism code (default "hsa" = human).
        """
        return kegg.gene_pathways(gene_symbol, organism=organism)

    @server.tool()
    def kegg_detail(pathway_id: str) -> dict:
        """Get detailed info for a KEGG pathway: genes, compounds, diseases.

        Args:
            pathway_id: KEGG pathway ID (e.g. "hsa04110" for cell cycle).
        """
        return kegg.pathway_info(pathway_id)

    @server.tool()
    def kegg_search(query: str, database: str = "pathway") -> dict:
        """Search KEGG databases by keyword.

        Args:
            query: Search term (e.g. "apoptosis", "insulin").
            database: KEGG database: "pathway", "genes", "compound", "disease", "drug".
        """
        return kegg.find(query, database=database)

    @server.tool()
    def reactome_gene(gene_symbol: str) -> dict:
        """Get Reactome biological pathways for a gene.

        Args:
            gene_symbol: Gene symbol (e.g. "TP53", "BRCA1").
        """
        return reactome.pathways_for_gene(gene_symbol)

    @server.tool()
    def reactome_detail(pathway_id: str) -> dict:
        """Get Reactome pathway details: summary, compartments, diagram link.

        Args:
            pathway_id: Reactome stable ID (e.g. "R-HSA-109582" for hemostasis).
        """
        return reactome.pathway_detail(pathway_id)

    @server.tool()
    def reactome_enrichment(gene_list: list[str]) -> dict:
        """Run Reactome pathway enrichment analysis on a gene list.

        Args:
            gene_list: List of gene symbols (e.g. ["TP53", "BRCA1", "ATM", "CHEK2"]).
        """
        return reactome.enrichment(gene_list)

    return server
