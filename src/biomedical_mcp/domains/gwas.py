"""GWAS domain — variant-trait associations (GWAS Catalog, EBI)."""

from fastmcp import FastMCP

from biomedical_mcp.gwas_catalog import GWASCatalog
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    gwas = GWASCatalog(cache)

    server = FastMCP("gwas")

    @server.tool(tags={"variant-review"})
    def variant_associations(rsid: str, limit: int = 25) -> dict:
        """GWAS associations for a variant: traits, p-values, effect sizes, genes.

        Args:
            rsid: dbSNP rsID (e.g. "rs1801133", "rs7903146").
            limit: Max associations (default 25).
        """
        return gwas.variant_associations(rsid, limit=limit)

    @server.tool(tags={"gene-lookup"})
    def gene_associations(gene_symbol: str, limit: int = 25) -> dict:
        """GWAS associations for a gene: traits, p-values, risk alleles.

        Args:
            gene_symbol: Gene symbol (e.g. "FTO", "APOE").
            limit: Max associations (default 25).
        """
        return gwas.gene_associations(gene_symbol, limit=limit)

    @server.tool()
    def trait_search(query: str, limit: int = 10) -> dict:
        """Search GWAS traits/phenotypes by keyword.

        Args:
            query: Trait keyword (e.g. "type 2 diabetes", "height").
            limit: Max results (default 10).
        """
        return gwas.trait_search(query, limit=limit)

    @server.tool(tags={"variant-review", "surveillance"})
    def new_for_variants(rsids: list[str], since_date: str = "") -> dict:
        """Check GWAS Catalog for new associations for a list of rsIDs.

        Useful for surveillance: "have any of my variants gained new GWAS hits?"
        Caps at 50 rsIDs per call.

        Args:
            rsids: List of dbSNP rsIDs (e.g. ["rs1801133", "rs7903146"]). Max 50.
            since_date: Optional ISO date filter (e.g. "2024-01-01"). Only return
                        associations published after this date.
        """
        return gwas.new_for_variants(rsids, since_date=since_date)

    return server
