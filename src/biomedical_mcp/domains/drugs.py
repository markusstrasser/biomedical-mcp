"""Drugs domain — compounds, mechanisms, bioactivity, labels, safety, PGx alleles."""

from fastmcp import FastMCP

from biomedical_mcp.chembl import ChEMBL
from biomedical_mcp.openfda import OpenFDA
from biomedical_mcp.pharmvar import PharmVar
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    chembl = ChEMBL(cache)
    openfda = OpenFDA(cache)
    pharmvar = PharmVar(cache)

    server = FastMCP("drugs")

    @server.tool()
    def compound(name: str | None = None, chembl_id: str | None = None) -> dict:
        """Get compound/molecule details: molecular properties, max phase, ATC, synonyms.

        Args:
            name: Drug/compound name (e.g. "ibuprofen").
            chembl_id: ChEMBL molecule ID (e.g. "CHEMBL521").
        """
        cid = chembl.resolve_compound(name=name, chembl_id=chembl_id)
        if not cid:
            return {"error": f"Could not resolve compound: name={name}, chembl_id={chembl_id}"}
        return chembl.compound(cid)

    @server.tool()
    def target(gene_symbol: str | None = None, chembl_id: str | None = None, uniprot_id: str | None = None) -> dict:
        """Get ChEMBL target details: type, organism, cross-references.

        Args:
            gene_symbol: Gene symbol (e.g. "EGFR").
            chembl_id: ChEMBL target ID (e.g. "CHEMBL203").
            uniprot_id: UniProt accession (e.g. "P00533").
        """
        result = chembl.target(gene_symbol=gene_symbol, chembl_id=chembl_id, uniprot_id=uniprot_id)
        if not result:
            return {"error": f"Target not found: gene_symbol={gene_symbol}, chembl_id={chembl_id}, uniprot_id={uniprot_id}"}
        return result

    @server.tool()
    def mechanism(compound_name: str | None = None, chembl_id: str | None = None) -> dict:
        """Get mechanism of action for a compound, with target details.

        Args:
            compound_name: Drug name (e.g. "ibuprofen").
            chembl_id: ChEMBL molecule ID.
        """
        cid = chembl.resolve_compound(name=compound_name, chembl_id=chembl_id)
        if not cid:
            return {"error": f"Could not resolve compound: name={compound_name}, chembl_id={chembl_id}"}
        mechanisms = chembl.mechanism(cid)
        return {"chembl_id": cid, "mechanisms": mechanisms}

    @server.tool()
    def bioactivity(target_chembl_id: str, limit: int = 25) -> dict:
        """Get bioactivity data (IC50/Ki/EC50) for tested compounds against a target.

        Args:
            target_chembl_id: ChEMBL target ID (e.g. "CHEMBL203"). Use drugs_target to find IDs.
            limit: Max results (default 25, max 100).
        """
        activities = chembl.bioactivity(target_chembl_id, limit=limit)
        return {"target_chembl_id": target_chembl_id, "activities": activities}

    @server.tool()
    def indications(compound_name: str | None = None, chembl_id: str | None = None, limit: int = 25) -> dict:
        """Get approved and clinical indications for a compound.

        Args:
            compound_name: Drug name.
            chembl_id: ChEMBL molecule ID.
            limit: Max results (default 25, max 100).
        """
        cid = chembl.resolve_compound(name=compound_name, chembl_id=chembl_id)
        if not cid:
            return {"error": f"Could not resolve compound: name={compound_name}, chembl_id={chembl_id}"}
        result = chembl.drug_indications(cid, limit=limit)
        return {"chembl_id": cid, "indications": result}

    @server.tool()
    def adverse_events(drug_name: str, reaction: str | None = None, serious: bool | None = None, limit: int = 25) -> dict:
        """Get FDA adverse event reports for a drug.

        Args:
            drug_name: Generic drug name (e.g. "ibuprofen", "metformin").
            reaction: Filter by reaction (e.g. "NAUSEA", "HEADACHE").
            serious: Filter by seriousness (True = serious only).
            limit: Max results (default 25, max 100).
        """
        return openfda.adverse_events(drug_name, reaction=reaction, serious=serious, limit=limit)

    @server.tool()
    def label(drug_name: str, sections: list[str] | None = None) -> dict:
        """Get FDA drug labeling: boxed warnings, dosing, interactions, contraindications.

        Args:
            drug_name: Generic drug name (e.g. "metformin", "warfarin").
            sections: Specific sections to return. Options: boxed_warning, warnings, dosage_and_administration, drug_interactions, adverse_reactions, indications_and_usage, contraindications, clinical_pharmacology, pharmacokinetics, use_in_specific_populations.
        """
        return openfda.drug_label(drug_name, sections=sections)

    @server.tool()
    def recalls(drug_name: str, classification: str | None = None, limit: int = 25) -> dict:
        """Get FDA drug recall enforcement reports.

        Args:
            drug_name: Generic drug name.
            classification: Filter by class ("Class I" = most serious, "Class II", "Class III").
            limit: Max results (default 25, max 100).
        """
        return openfda.recalls(drug_name, classification=classification, limit=limit)

    @server.tool(tags={"pgx"})
    def star_alleles(gene: str, limit: int = 50) -> dict:
        """Get star allele definitions from PharmVar for a PGx gene.

        Requires PHARMVAR_API_KEY env var (free at pharmvar.org).

        Args:
            gene: PGx gene (e.g. "CYP2D6", "CYP2C19", "CYP3A4").
            limit: Max alleles (default 50).
        """
        return pharmvar.star_alleles(gene, limit=limit)

    return server
