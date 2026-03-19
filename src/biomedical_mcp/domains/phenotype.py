"""Phenotype domain — HPO terms, gene-phenotype associations, disease profiles.

Advisory lookups only — phenotype concordance validation stays in selve.
"""

from fastmcp import FastMCP

from biomedical_mcp.hpo import HPO
from biomedical_mcp.monarch import Monarch
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    hpo = HPO(cache)
    monarch = Monarch(cache)

    server = FastMCP("phenotype")

    @server.tool(tags={"phenotype"})
    def hpo_search(query: str, limit: int = 10) -> dict:
        """Search HPO terms by keyword (e.g. "seizure", "cardiomyopathy").

        Args:
            query: Search term.
            limit: Max results (default 10).
        """
        return hpo.search(query, limit=limit)

    @server.tool(tags={"phenotype"})
    def hpo_term(hpo_id: str) -> dict:
        """Get HPO term details: definition, synonyms, parent/child terms.

        Args:
            hpo_id: HPO identifier (e.g. "HP:0001250" for seizure).
        """
        return hpo.term(hpo_id)

    @server.tool(tags={"phenotype", "gene-lookup"})
    def hpo_genes(hpo_id: str) -> dict:
        """Get genes associated with an HPO phenotype term.

        Args:
            hpo_id: HPO identifier (e.g. "HP:0001250").
        """
        return hpo.genes(hpo_id)

    @server.tool(tags={"phenotype", "gene-lookup"})
    def gene_phenotypes(gene_symbol: str) -> dict:
        """Get HPO phenotypes associated with a gene. Advisory — not for concordance validation.

        Args:
            gene_symbol: Gene symbol (e.g. "SCN1A", "EBF3").
        """
        return hpo.gene_phenotypes(gene_symbol)

    @server.tool(tags={"phenotype"})
    def disease_phenotypes(disease_id: str, limit: int = 50) -> dict:
        """Get phenotype profile for a disease from Monarch (cross-species).

        Args:
            disease_id: OMIM, MONDO, or Orphanet ID (e.g. "OMIM:601419", "MONDO:0007254").
            limit: Max phenotypes (default 50).
        """
        return monarch.disease_phenotypes(disease_id, limit=limit)

    return server
