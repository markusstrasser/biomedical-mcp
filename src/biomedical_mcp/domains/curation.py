"""Curation domain — ClinGen gene-disease validity and dosage sensitivity."""

from fastmcp import FastMCP

from biomedical_mcp.clingen import ClinGen
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    clingen = ClinGen(cache)

    server = FastMCP("curation")

    @server.tool(tags={"variant-review", "gene-lookup", "curation"})
    def gene_validity(gene_symbol: str) -> dict:
        """ClinGen gene-disease validity: classification (Definitive/Strong/Moderate/Limited/Disputed/Refuted), mode of inheritance, expert panel.

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "TP53", "SCN5A").
        """
        return clingen.gene_validity(gene_symbol)

    @server.tool(tags={"curation"})
    def disease_validity(disease_label: str) -> dict:
        """Search ClinGen gene-disease validity curations by disease name.

        Args:
            disease_label: Disease name or partial match (e.g. "breast cancer", "cardiomyopathy").
        """
        return clingen.disease_validity(disease_label)

    @server.tool(tags={"curation"})
    def validity_by_classification(classification: str = "Definitive") -> dict:
        """List all gene-disease pairs with a specific ClinGen classification level.

        Args:
            classification: One of "Definitive", "Strong", "Moderate", "Limited", "Disputed", "Refuted" (default: "Definitive").
        """
        return clingen.classification_summary(classification)

    @server.tool(tags={"variant-review", "gene-lookup", "curation"})
    def gene_dosage(gene_symbol: str) -> dict:
        """ClinGen dosage sensitivity: haploinsufficiency and triplosensitivity scores.

        Scores: "Sufficient evidence" / "Some evidence" / "Little evidence" / "No evidence" / "Autosomal recessive" / "Gene associated with AR condition" / "Dosage sensitivity unlikely".

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "TP53").
        """
        return clingen.gene_dosage(gene_symbol)

    return server
