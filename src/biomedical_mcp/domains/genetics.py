"""Genetics domain — gene annotation, IDs, cross-references, sequences."""

from fastmcp import FastMCP

from biomedical_mcp.ensembl import Ensembl
from biomedical_mcp.mygene import MyGene
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    mygene = MyGene(cache)
    ensembl = Ensembl(cache)

    server = FastMCP("genetics")

    @server.tool()
    def gene_info(gene_id: str) -> dict:
        """Get gene annotations: Entrez/Ensembl/UniProt IDs, KEGG pathways, GO terms, genomic position, aliases.

        Args:
            gene_id: Gene symbol (e.g. "CYP2D6"), Entrez ID (e.g. "1565"), or Ensembl ID (e.g. "ENSG00000100197").
        """
        result = mygene.gene_info(gene_id)
        if not result:
            return {"error": f"Gene not found: {gene_id}"}
        return result

    @server.tool()
    def gene_search(query: str, species: str = "human", limit: int = 10) -> dict:
        """Search genes by keyword.

        Args:
            query: Search term (e.g. "cytochrome P450", "breast cancer").
            species: Species filter (default "human").
            limit: Max results (default 10, max 50).
        """
        results = mygene.search(query, species=species, limit=limit)
        return {"query": query, "count": len(results), "results": results}

    @server.tool()
    def ensembl_gene(symbol: str | None = None, ensembl_id: str | None = None, species: str = "homo_sapiens") -> dict:
        """Get gene annotation from Ensembl: ID, biotype, chromosomal location, description.

        Args:
            symbol: Gene symbol (e.g. "BRCA1").
            ensembl_id: Ensembl gene ID (e.g. "ENSG00000012048"). Takes priority.
            species: Species (default "homo_sapiens").
        """
        if ensembl_id:
            return ensembl.gene_by_id(ensembl_id)
        if symbol:
            return ensembl.gene_by_symbol(symbol, species=species)
        return {"error": "Provide either symbol or ensembl_id"}

    @server.tool()
    def ensembl_xrefs(ensembl_id: str) -> dict:
        """Get cross-references from Ensembl ID to UniProt, HGNC, RefSeq, CCDS, etc.

        Args:
            ensembl_id: Ensembl gene/transcript ID (e.g. "ENSG00000012048").
        """
        return ensembl.xrefs(ensembl_id)

    @server.tool()
    def ensembl_sequence(ensembl_id: str, seq_type: str = "genomic") -> dict:
        """Get nucleotide or protein sequence for an Ensembl ID.

        Args:
            ensembl_id: Ensembl gene/transcript ID.
            seq_type: Sequence type: "genomic", "cds", "cdna", or "protein".
        """
        return ensembl.sequence(ensembl_id, seq_type=seq_type)

    return server
