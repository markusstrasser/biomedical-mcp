"""Variants domain — genetic variant annotations, ClinVar, batch lookup."""

from fastmcp import FastMCP

from biomedical_mcp.myvariant import MyVariant
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    myvariant = MyVariant(cache)

    server = FastMCP("variants")

    @server.tool()
    def lookup(variant_id: str) -> dict:
        """Look up variant annotations: ClinVar significance, gnomAD allele frequency, CADD score, dbNSFP predictions, dbSNP, CIViC, SNPedia.

        Args:
            variant_id: rsID (e.g. "rs1805007"), HGVS (e.g. "chr16:g.89919683C>T"), or dbSNP-style ID.
        """
        return myvariant.lookup(variant_id)

    @server.tool()
    def clinvar(gene_symbol: str, significance: str | None = None, limit: int = 25) -> dict:
        """Search ClinVar variants for a gene, optionally filtered by clinical significance.

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "TP53").
            significance: Filter (e.g. "Pathogenic", "Likely pathogenic", "Benign").
            limit: Max results (default 25, max 100).
        """
        results = myvariant.clinvar_search(gene_symbol, significance=significance, limit=limit)
        return {"gene": gene_symbol, "significance_filter": significance, "count": len(results), "variants": results}

    @server.tool()
    def batch(variant_ids: list[str]) -> dict:
        """Batch lookup up to 100 variants at once. Returns same annotations as variants_lookup.

        Args:
            variant_ids: List of rsIDs or HGVS IDs (max 100).
        """
        results = myvariant.batch(variant_ids)
        return {"count": len(results), "variants": results}

    return server
