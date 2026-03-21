"""Rare disease domain — Orphanet (Orphadata) rare disease database."""

from fastmcp import FastMCP

from biomedical_mcp.orphanet import Orphanet
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    orphanet = Orphanet(cache)

    server = FastMCP("rare_disease")

    @server.tool(tags={"gene-lookup", "rare-disease"})
    def gene_rare_diseases(gene_symbol: str) -> dict:
        """Rare diseases associated with a gene (Orphanet).

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "CFTR", "DMD").
        """
        return orphanet.gene_diseases(gene_symbol)

    @server.tool(tags={"rare-disease"})
    def disease_genes(orphacode: int) -> dict:
        """Genes associated with a rare disease — includes association type and external refs.

        Args:
            orphacode: Orphanet disease code (e.g. 166024).
        """
        return orphanet.disease_genes(orphacode)

    @server.tool(tags={"rare-disease"})
    def disease_natural_history(orphacode: int) -> dict:
        """Inheritance mode, age of onset, and age of death for a rare disease.

        Args:
            orphacode: Orphanet disease code (e.g. 166024).
        """
        return orphanet.natural_history(orphacode)

    @server.tool(tags={"rare-disease"})
    def disease_epidemiology(orphacode: int) -> dict:
        """Prevalence and incidence data for a rare disease.

        Args:
            orphacode: Orphanet disease code (e.g. 166024).
        """
        return orphanet.epidemiology(orphacode)

    @server.tool(tags={"rare-disease"})
    def disease_phenotypes(orphacode: int) -> dict:
        """HPO phenotypes associated with a rare disease, with frequency.

        Args:
            orphacode: Orphanet disease code (e.g. 166024).
        """
        return orphanet.phenotypes(orphacode)

    @server.tool(tags={"rare-disease"})
    def disease_cross_references(orphacode: int) -> dict:
        """Cross-references to ICD-10, ICD-11, OMIM for a rare disease.

        Args:
            orphacode: Orphanet disease code (e.g. 166024).
        """
        return orphanet.cross_reference(orphacode)

    @server.tool(tags={"rare-disease"})
    def omim_to_orphanet(omim_code: str) -> dict:
        """Find Orphanet rare diseases by OMIM code.

        Args:
            omim_code: OMIM code (e.g. "113705" for BRCA1).
        """
        return orphanet.search_by_omim(omim_code)

    return server
