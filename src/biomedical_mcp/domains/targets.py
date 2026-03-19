"""Targets domain — disease-gene associations, pharmacogenetics, search."""

from fastmcp import FastMCP

from biomedical_mcp.opentargets import OpenTargets
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    ot = OpenTargets(cache)

    server = FastMCP("targets")

    @server.tool()
    def search(query: str, entity_type: str | None = None, limit: int = 10) -> dict:
        """Search Open Targets for targets, diseases, or drugs.

        Args:
            query: Search term (gene name, disease, drug).
            entity_type: Filter to "target", "disease", or "drug". None searches all.
            limit: Max results (default 10, max 50).
        """
        hits = ot.search(query, entity_type=entity_type, limit=limit)
        return {"query": query, "total": len(hits), "hits": hits}

    @server.tool()
    def target_info(gene_symbol: str | None = None, ensembl_id: str | None = None) -> dict:
        """Get target details: description, tractability, safety, known drugs.

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "CYP2D6").
            ensembl_id: Ensembl gene ID (e.g. "ENSG00000012048"). Takes priority.
        """
        eid = ot.resolve_target(gene_symbol=gene_symbol, ensembl_id=ensembl_id)
        if not eid:
            return {"error": f"Could not resolve target: gene_symbol={gene_symbol}, ensembl_id={ensembl_id}"}
        info = ot.target_info(eid)
        return info if info else {"error": f"Target not found: {eid}"}

    @server.tool()
    def disease_associations(gene_symbol: str | None = None, ensembl_id: str | None = None, limit: int = 25) -> dict:
        """Get diseases associated with a gene, scored by evidence strength.

        Args:
            gene_symbol: Gene symbol (e.g. "TP53").
            ensembl_id: Ensembl gene ID.
            limit: Max associations (default 25, max 100).
        """
        eid = ot.resolve_target(gene_symbol=gene_symbol, ensembl_id=ensembl_id)
        if not eid:
            return {"error": f"Could not resolve target: gene_symbol={gene_symbol}, ensembl_id={ensembl_id}"}
        return ot.disease_associations(eid, limit=limit)

    @server.tool()
    def disease_targets(disease_id: str, limit: int = 25) -> dict:
        """Get top gene targets for a disease by association score.

        Args:
            disease_id: EFO disease ID (e.g. "EFO_0000305" for breast carcinoma). Use targets_search to find IDs.
            limit: Max targets (default 25, max 100).
        """
        return ot.disease_targets(disease_id, limit=limit)

    @server.tool()
    def pharmacogenetics(gene_symbol: str | None = None, ensembl_id: str | None = None) -> dict:
        """Get pharmacogenetic interactions for a gene: drug-variant-phenotype mappings.

        Args:
            gene_symbol: Gene symbol (e.g. "CYP2D6", "CYP3A4").
            ensembl_id: Ensembl gene ID.
        """
        eid = ot.resolve_target(gene_symbol=gene_symbol, ensembl_id=ensembl_id)
        if not eid:
            return {"error": f"Could not resolve target: gene_symbol={gene_symbol}, ensembl_id={ensembl_id}"}
        return ot.pharmacogenetics(eid)

    @server.tool()
    def drug_info(drug_name: str | None = None, chembl_id: str | None = None) -> dict:
        """Get drug details from Open Targets: type, phase, mechanisms of action, indications.

        Args:
            drug_name: Drug name (e.g. "ibuprofen").
            chembl_id: ChEMBL compound ID (e.g. "CHEMBL521"). Takes priority.
        """
        cid = ot.resolve_drug(drug_name=drug_name, chembl_id=chembl_id)
        if not cid:
            return {"error": f"Could not resolve drug: drug_name={drug_name}, chembl_id={chembl_id}"}
        info = ot.drug_info(cid)
        return info if info else {"error": f"Drug not found: {cid}"}

    return server
