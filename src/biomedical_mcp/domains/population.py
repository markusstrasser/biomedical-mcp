"""Population domain — allele frequencies, constraint metrics (gnomAD)."""

from fastmcp import FastMCP

from biomedical_mcp.gnomad import GnomAD
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    gnomad = GnomAD(cache)

    server = FastMCP("population")

    @server.tool(tags={"variant-review", "population"})
    def variant_frequency(variant_id: str) -> dict:
        """Get population allele frequencies from gnomAD (genome + exome) by population.

        Args:
            variant_id: gnomAD-style ID (e.g. "1-55516888-G-A"), rsID, or HGVS.
        """
        return gnomad.variant_frequency(variant_id)

    @server.tool(tags={"gene-lookup", "population"})
    def gene_constraint(gene_symbol: str) -> dict:
        """Get gnomAD gene constraint: pLI, LOEUF (oe_lof_upper), mis_z, syn_z, lof_z.

        High pLI (>0.9) or low LOEUF (<0.35) = gene is intolerant to loss-of-function.

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "TP53").
        """
        return gnomad.gene_constraint(gene_symbol)

    return server
