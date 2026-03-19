"""Proteins domain — structure, function, interactions, variants."""

from fastmcp import FastMCP

from biomedical_mcp.uniprot import UniProt
from biomedical_mcp.alphafold import AlphaFold
from biomedical_mcp.stringdb import StringDB
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    uniprot = UniProt(cache)
    alphafold = AlphaFold(cache)
    stringdb = StringDB(cache)

    server = FastMCP("proteins")

    @server.tool()
    def protein(accession: str | None = None, gene_symbol: str | None = None) -> dict:
        """Get protein details: function, subcellular location, domains, disease associations, GO terms.

        Args:
            accession: UniProt accession (e.g. "P38398" for BRCA1).
            gene_symbol: Gene symbol (e.g. "BRCA1"). Resolved to human reviewed entry.
        """
        result = uniprot.protein(accession=accession, gene_symbol=gene_symbol)
        if not result:
            return {"error": f"Protein not found: accession={accession}, gene_symbol={gene_symbol}"}
        return result

    @server.tool()
    def variants(accession: str | None = None, gene_symbol: str | None = None, limit: int = 50) -> dict:
        """Get known natural variants and mutagenesis data for a protein.

        Args:
            accession: UniProt accession.
            gene_symbol: Gene symbol (e.g. "TP53").
            limit: Max variants (default 50).
        """
        result = uniprot.variants(accession=accession, gene_symbol=gene_symbol, limit=limit)
        if not result:
            return {"error": f"Protein not found: accession={accession}, gene_symbol={gene_symbol}"}
        return result

    @server.tool()
    def search(query: str, organism: str = "human", limit: int = 10) -> dict:
        """Search UniProt proteins by keyword.

        Args:
            query: Search term (e.g. "kinase", "DNA repair").
            organism: Organism filter (default "human").
            limit: Max results (default 10, max 50).
        """
        results = uniprot.search(query, organism=organism, limit=limit)
        return {"query": query, "count": len(results), "results": results}

    @server.tool()
    def structure(uniprot_id: str) -> dict:
        """Get AlphaFold predicted structure info: model URL, confidence (pLDDT), PAE image.

        Args:
            uniprot_id: UniProt accession (e.g. "P38398" for BRCA1, "P04637" for TP53).
        """
        return alphafold.prediction(uniprot_id)

    @server.tool()
    def interactions(protein: str, species: str = "9606", min_score: int = 700, limit: int = 25) -> dict:
        """Get protein-protein interactions from STRING DB.

        Args:
            protein: Protein/gene name (e.g. "TP53", "EGFR").
            species: NCBI taxonomy ID (default "9606" = human).
            min_score: Minimum combined score 0-1000 (default 700 = high confidence).
            limit: Max interaction partners (default 25).
        """
        return stringdb.interactions(protein, species=species, min_score=min_score, limit=limit)

    @server.tool()
    def enrichment(proteins: list[str], species: str = "9606") -> dict:
        """Get functional enrichment (GO, KEGG, Reactome terms) for a set of proteins.

        Args:
            proteins: List of protein/gene names (e.g. ["TP53", "BRCA1", "ATM"]).
            species: NCBI taxonomy ID (default "9606" = human).
        """
        return stringdb.functional_enrichment(proteins, species=species)

    return server
