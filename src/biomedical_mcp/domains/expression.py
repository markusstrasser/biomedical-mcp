"""Expression domain — tissue expression, eQTLs (GTEx Portal v2)."""

from fastmcp import FastMCP

from biomedical_mcp.gtex import GTEx
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    gtex = GTEx(cache)

    server = FastMCP("expression")

    @server.tool(tags={"gene-lookup", "expression"})
    def gene_expression(gene_symbol: str | None = None, ensembl_id: str | None = None) -> dict:
        """Median TPM by tissue for a gene. Tissues below TPM 1.0 are filtered out.

        Args:
            gene_symbol: Gene symbol (e.g. "CYP2D6", "BRCA1").
            ensembl_id: Ensembl gene ID (e.g. "ENSG00000100197").
        """
        return gtex.gene_expression(gene_symbol=gene_symbol, ensembl_id=ensembl_id)

    @server.tool(tags={"gene-lookup", "expression"})
    def top_tissues(gene_symbol: str, limit: int = 10) -> dict:
        """Top N tissues by expression level for a gene.

        Args:
            gene_symbol: Gene symbol.
            limit: Number of top tissues (default 10).
        """
        return gtex.top_tissues(gene_symbol, limit=limit)

    @server.tool(tags={"variant-review", "expression"})
    def eqtl(gene_symbol: str, tissue: str | None = None) -> dict:
        """Significant eQTLs for a gene, optionally filtered by tissue.

        Args:
            gene_symbol: Gene symbol.
            tissue: GTEx tissue ID (e.g. "Liver", "Whole_Blood"). None = all tissues.
        """
        return gtex.eqtl(gene_symbol, tissue=tissue)

    return server
