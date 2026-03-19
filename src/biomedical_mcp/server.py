"""Biomedical MCP server — 9 APIs, 29 tools for genes, drugs, variants, proteins, trials."""

import logging
import os
from contextlib import asynccontextmanager
from pathlib import Path

from fastmcp import FastMCP, Context

from biomedical_mcp.cache import Cache
from biomedical_mcp.opentargets import OpenTargets
from biomedical_mcp.chembl import ChEMBL
from biomedical_mcp.clinicaltrials import ClinicalTrials
from biomedical_mcp.icd10 import ICD10
from biomedical_mcp.npi import NPI
from biomedical_mcp.myvariant import MyVariant
from biomedical_mcp.openfda import OpenFDA
from biomedical_mcp.uniprot import UniProt
from biomedical_mcp.mygene import MyGene

log = logging.getLogger(__name__)

DEFAULT_DATA_DIR = Path.home() / ".local" / "share" / "biomedical-mcp"


def create_mcp(data_dir: Path | None = None) -> FastMCP:
    data_dir = data_dir or Path(os.environ.get("BIOMEDICAL_MCP_DATA", DEFAULT_DATA_DIR))

    @asynccontextmanager
    async def lifespan(server):
        data_dir.mkdir(parents=True, exist_ok=True)
        cache = Cache(data_dir / "cache.db")
        ot = OpenTargets(cache)
        chembl = ChEMBL(cache)
        ct = ClinicalTrials(cache)
        icd10 = ICD10(cache)
        npi = NPI(cache)
        myvariant = MyVariant(cache)
        openfda = OpenFDA(cache)
        uniprot = UniProt(cache)
        mygene = MyGene(cache)
        log.info("biomedical-mcp started (cache: %s)", data_dir / "cache.db")
        yield {
            "cache": cache, "ot": ot, "chembl": chembl, "ct": ct,
            "icd10": icd10, "npi": npi, "myvariant": myvariant,
            "openfda": openfda, "uniprot": uniprot, "mygene": mygene,
        }

    mcp = FastMCP(
        "biomedical",
        instructions=(
            "Biomedical data lookup via 9 APIs: Open Targets, ChEMBL, ClinicalTrials.gov, ICD-10, NPI, "
            "MyVariant.info, OpenFDA, UniProt, MyGene.info.\n\n"
            "Gene/target tools: ot_target_info, ot_disease_associations, ot_pharmacogenetics, chembl_target, chembl_bioactivity, gene_info, gene_search\n"
            "Drug tools: ot_drug_info, chembl_compound, chembl_mechanism, chembl_drug_indications, openfda_adverse_events, openfda_drug_label, openfda_recalls\n"
            "Variant tools: variant_lookup, variant_clinvar, variant_batch\n"
            "Protein tools: uniprot_protein, uniprot_variants, uniprot_search\n"
            "Disease tools: ot_disease_targets, ct_search, ct_trial_detail, ct_stats\n"
            "Code/provider tools: icd10_search, icd10_lookup, npi_search, npi_lookup\n"
            "General search: ot_search"
        ),
        lifespan=lifespan,
    )

    # ── Open Targets ──────────────────────────────────────────────

    @mcp.tool()
    def ot_search(ctx: Context, query: str, entity_type: str | None = None, limit: int = 10) -> dict:
        """Search Open Targets for targets, diseases, or drugs.

        Args:
            query: Search term (gene name, disease, drug).
            entity_type: Filter to "target", "disease", or "drug". None searches all.
            limit: Max results (default 10, max 50).
        """
        ot: OpenTargets = ctx.request_context["ot"]
        hits = ot.search(query, entity_type=entity_type, limit=limit)
        return {"query": query, "total": len(hits), "hits": hits}

    @mcp.tool()
    def ot_target_info(ctx: Context, gene_symbol: str | None = None, ensembl_id: str | None = None) -> dict:
        """Get target details: description, tractability, safety, known drugs.

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "CYP2D6"). Resolved via OT search.
            ensembl_id: Ensembl gene ID (e.g. "ENSG00000012048"). Takes priority.
        """
        ot: OpenTargets = ctx.request_context["ot"]
        eid = ot.resolve_target(gene_symbol=gene_symbol, ensembl_id=ensembl_id)
        if not eid:
            return {"error": f"Could not resolve target: gene_symbol={gene_symbol}, ensembl_id={ensembl_id}"}
        info = ot.target_info(eid)
        if not info:
            return {"error": f"Target not found: {eid}"}
        return info

    @mcp.tool()
    def ot_disease_associations(ctx: Context, gene_symbol: str | None = None, ensembl_id: str | None = None, limit: int = 25) -> dict:
        """Get diseases associated with a gene, scored by evidence strength.

        Args:
            gene_symbol: Gene symbol (e.g. "TP53").
            ensembl_id: Ensembl gene ID.
            limit: Max associations (default 25, max 100).
        """
        ot: OpenTargets = ctx.request_context["ot"]
        eid = ot.resolve_target(gene_symbol=gene_symbol, ensembl_id=ensembl_id)
        if not eid:
            return {"error": f"Could not resolve target: gene_symbol={gene_symbol}, ensembl_id={ensembl_id}"}
        return ot.disease_associations(eid, limit=limit)

    @mcp.tool()
    def ot_disease_targets(ctx: Context, disease_id: str, limit: int = 25) -> dict:
        """Get top gene targets for a disease by association score.

        Args:
            disease_id: EFO disease ID (e.g. "EFO_0000305" for breast carcinoma). Use ot_search to find IDs.
            limit: Max targets (default 25, max 100).
        """
        ot: OpenTargets = ctx.request_context["ot"]
        return ot.disease_targets(disease_id, limit=limit)

    @mcp.tool()
    def ot_pharmacogenetics(ctx: Context, gene_symbol: str | None = None, ensembl_id: str | None = None) -> dict:
        """Get pharmacogenetic interactions for a gene: drug-variant-phenotype mappings.

        Args:
            gene_symbol: Gene symbol (e.g. "CYP2D6", "CYP3A4").
            ensembl_id: Ensembl gene ID.
        """
        ot: OpenTargets = ctx.request_context["ot"]
        eid = ot.resolve_target(gene_symbol=gene_symbol, ensembl_id=ensembl_id)
        if not eid:
            return {"error": f"Could not resolve target: gene_symbol={gene_symbol}, ensembl_id={ensembl_id}"}
        return ot.pharmacogenetics(eid)

    @mcp.tool()
    def ot_drug_info(ctx: Context, drug_name: str | None = None, chembl_id: str | None = None) -> dict:
        """Get drug details: type, phase, mechanisms of action, indications.

        Args:
            drug_name: Drug name (e.g. "ibuprofen"). Resolved via OT search.
            chembl_id: ChEMBL compound ID (e.g. "CHEMBL521"). Takes priority.
        """
        ot: OpenTargets = ctx.request_context["ot"]
        cid = ot.resolve_drug(drug_name=drug_name, chembl_id=chembl_id)
        if not cid:
            return {"error": f"Could not resolve drug: drug_name={drug_name}, chembl_id={chembl_id}"}
        info = ot.drug_info(cid)
        if not info:
            return {"error": f"Drug not found: {cid}"}
        return info

    # ── ChEMBL ────────────────────────────────────────────────────

    @mcp.tool()
    def chembl_compound(ctx: Context, name: str | None = None, chembl_id: str | None = None) -> dict:
        """Get compound/molecule details: molecular properties, max phase, ATC, synonyms.

        Args:
            name: Drug/compound name (e.g. "ibuprofen").
            chembl_id: ChEMBL molecule ID (e.g. "CHEMBL521").
        """
        chembl_client: ChEMBL = ctx.request_context["chembl"]
        cid = chembl_client.resolve_compound(name=name, chembl_id=chembl_id)
        if not cid:
            return {"error": f"Could not resolve compound: name={name}, chembl_id={chembl_id}"}
        return chembl_client.compound(cid)

    @mcp.tool()
    def chembl_target(ctx: Context, gene_symbol: str | None = None, chembl_id: str | None = None, uniprot_id: str | None = None) -> dict:
        """Get ChEMBL target details: type, organism, cross-references.

        Args:
            gene_symbol: Gene symbol to search (e.g. "EGFR").
            chembl_id: ChEMBL target ID (e.g. "CHEMBL203").
            uniprot_id: UniProt accession (e.g. "P00533").
        """
        chembl_client: ChEMBL = ctx.request_context["chembl"]
        result = chembl_client.target(gene_symbol=gene_symbol, chembl_id=chembl_id, uniprot_id=uniprot_id)
        if not result:
            return {"error": f"Target not found: gene_symbol={gene_symbol}, chembl_id={chembl_id}, uniprot_id={uniprot_id}"}
        return result

    @mcp.tool()
    def chembl_mechanism(ctx: Context, compound_name: str | None = None, chembl_id: str | None = None) -> dict:
        """Get mechanism of action for a compound, with target details.

        Args:
            compound_name: Drug name (e.g. "ibuprofen").
            chembl_id: ChEMBL molecule ID.
        """
        chembl_client: ChEMBL = ctx.request_context["chembl"]
        cid = chembl_client.resolve_compound(name=compound_name, chembl_id=chembl_id)
        if not cid:
            return {"error": f"Could not resolve compound: name={compound_name}, chembl_id={chembl_id}"}
        mechanisms = chembl_client.mechanism(cid)
        return {"chembl_id": cid, "mechanisms": mechanisms}

    @mcp.tool()
    def chembl_bioactivity(ctx: Context, target_chembl_id: str, limit: int = 25) -> dict:
        """Get bioactivity data (IC50/Ki/EC50) for tested compounds against a target.

        Args:
            target_chembl_id: ChEMBL target ID (e.g. "CHEMBL203"). Use chembl_target to find IDs.
            limit: Max results (default 25, max 100).
        """
        chembl_client: ChEMBL = ctx.request_context["chembl"]
        activities = chembl_client.bioactivity(target_chembl_id, limit=limit)
        return {"target_chembl_id": target_chembl_id, "activities": activities}

    @mcp.tool()
    def chembl_drug_indications(ctx: Context, compound_name: str | None = None, chembl_id: str | None = None, limit: int = 25) -> dict:
        """Get approved and clinical indications for a compound.

        Args:
            compound_name: Drug name.
            chembl_id: ChEMBL molecule ID.
            limit: Max results (default 25, max 100).
        """
        chembl_client: ChEMBL = ctx.request_context["chembl"]
        cid = chembl_client.resolve_compound(name=compound_name, chembl_id=chembl_id)
        if not cid:
            return {"error": f"Could not resolve compound: name={compound_name}, chembl_id={chembl_id}"}
        indications = chembl_client.drug_indications(cid, limit=limit)
        return {"chembl_id": cid, "indications": indications}

    # ── ClinicalTrials.gov ────────────────────────────────────────

    @mcp.tool()
    def ct_search(
        ctx: Context,
        condition: str | None = None,
        intervention: str | None = None,
        gene: str | None = None,
        status: str | None = None,
        limit: int = 20,
    ) -> dict:
        """Search clinical trials by condition, intervention, gene, or status.

        Args:
            condition: Disease/condition (e.g. "breast cancer").
            intervention: Drug or therapy (e.g. "pembrolizumab").
            gene: Gene name for gene-related trials (e.g. "BRCA1").
            status: Filter by status: RECRUITING, COMPLETED, ACTIVE_NOT_RECRUITING, etc.
            limit: Max results (default 20, max 100).
        """
        ct_client: ClinicalTrials = ctx.request_context["ct"]
        return ct_client.search(condition=condition, intervention=intervention, gene=gene, status=status, limit=limit)

    @mcp.tool()
    def ct_trial_detail(ctx: Context, nct_id: str) -> dict:
        """Get full trial details: eligibility, outcomes, sponsor, design.

        Args:
            nct_id: ClinicalTrials.gov identifier (e.g. "NCT04379596").
        """
        ct_client: ClinicalTrials = ctx.request_context["ct"]
        return ct_client.trial_detail(nct_id)

    @mcp.tool()
    def ct_stats(ctx: Context, condition: str) -> dict:
        """Get trial counts by phase and status for a condition.

        Args:
            condition: Disease/condition (e.g. "lung cancer").
        """
        ct_client: ClinicalTrials = ctx.request_context["ct"]
        return ct_client.stats(condition)

    # ── ICD-10 ────────────────────────────────────────────────────

    @mcp.tool()
    def icd10_search(ctx: Context, query: str, limit: int = 10) -> dict:
        """Search ICD-10-CM codes by description.

        Args:
            query: Search term (e.g. "diabetes", "hypertension").
            limit: Max results (default 10, max 50).
        """
        icd10_client: ICD10 = ctx.request_context["icd10"]
        results = icd10_client.search(query, limit=limit)
        return {"query": query, "results": results}

    @mcp.tool()
    def icd10_lookup(ctx: Context, code: str) -> dict:
        """Look up an ICD-10-CM code for its description and subcodes.

        Args:
            code: ICD-10-CM code (e.g. "E11", "E11.65").
        """
        icd10_client: ICD10 = ctx.request_context["icd10"]
        result = icd10_client.lookup(code)
        if not result:
            return {"error": f"ICD-10 code not found: {code}"}
        return result

    # ── NPI ────────────────────────────────────────────────────────

    @mcp.tool()
    def npi_search(
        ctx: Context,
        name: str | None = None,
        specialty: str | None = None,
        city: str | None = None,
        state: str | None = None,
        limit: int = 10,
    ) -> dict:
        """Search the NPI registry for healthcare providers.

        Args:
            name: Provider name ("John Smith" or just "Smith").
            specialty: Taxonomy description (e.g. "Cardiology", "Internal Medicine").
            city: City name.
            state: Two-letter state code (e.g. "CA", "NY").
            limit: Max results (default 10, max 200).
        """
        npi_client: NPI = ctx.request_context["npi"]
        results = npi_client.search(name=name, specialty=specialty, city=city, state=state, limit=limit)
        return {"results": results, "count": len(results)}

    @mcp.tool()
    def npi_lookup(ctx: Context, npi: str) -> dict:
        """Look up a specific NPI number for full provider details.

        Args:
            npi: 10-digit NPI number.
        """
        npi_client: NPI = ctx.request_context["npi"]
        result = npi_client.lookup(npi)
        if not result:
            return {"error": f"NPI not found: {npi}"}
        return result

    # ── MyVariant.info ─────────────────────────────────────────────

    @mcp.tool()
    def variant_lookup(ctx: Context, variant_id: str) -> dict:
        """Look up variant annotations: ClinVar significance, gnomAD allele frequency, CADD score, dbNSFP predictions, dbSNP, CIViC, SNPedia.

        Args:
            variant_id: rsID (e.g. "rs1805007"), HGVS (e.g. "chr16:g.89919683C>T"), or dbSNP-style ID.
        """
        mv: MyVariant = ctx.request_context["myvariant"]
        return mv.lookup(variant_id)

    @mcp.tool()
    def variant_clinvar(ctx: Context, gene_symbol: str, significance: str | None = None, limit: int = 25) -> dict:
        """Search ClinVar variants for a gene, optionally filtered by clinical significance.

        Args:
            gene_symbol: Gene symbol (e.g. "BRCA1", "TP53").
            significance: Filter by significance (e.g. "Pathogenic", "Likely pathogenic", "Benign").
            limit: Max results (default 25, max 100).
        """
        mv: MyVariant = ctx.request_context["myvariant"]
        variants = mv.clinvar_search(gene_symbol, significance=significance, limit=limit)
        return {"gene": gene_symbol, "significance_filter": significance, "count": len(variants), "variants": variants}

    @mcp.tool()
    def variant_batch(ctx: Context, variant_ids: list[str]) -> dict:
        """Batch lookup up to 100 variants at once. Returns same annotations as variant_lookup.

        Args:
            variant_ids: List of rsIDs or HGVS IDs (max 100).
        """
        mv: MyVariant = ctx.request_context["myvariant"]
        results = mv.batch(variant_ids)
        return {"count": len(results), "variants": results}

    # ── OpenFDA ────────────────────────────────────────────────────

    @mcp.tool()
    def openfda_adverse_events(ctx: Context, drug_name: str, reaction: str | None = None, serious: bool | None = None, limit: int = 25) -> dict:
        """Get FDA adverse event reports for a drug.

        Args:
            drug_name: Generic drug name (e.g. "ibuprofen", "metformin").
            reaction: Filter by reaction (e.g. "NAUSEA", "HEADACHE").
            serious: Filter by seriousness (True = serious only, False = non-serious only).
            limit: Max results (default 25, max 100).
        """
        fda: OpenFDA = ctx.request_context["openfda"]
        return fda.adverse_events(drug_name, reaction=reaction, serious=serious, limit=limit)

    @mcp.tool()
    def openfda_drug_label(ctx: Context, drug_name: str, sections: list[str] | None = None) -> dict:
        """Get FDA drug labeling: boxed warnings, dosing, interactions, contraindications.

        Args:
            drug_name: Generic drug name (e.g. "metformin", "warfarin").
            sections: Specific sections to return. Options: boxed_warning, warnings, dosage_and_administration, drug_interactions, adverse_reactions, indications_and_usage, contraindications, clinical_pharmacology, pharmacokinetics, use_in_specific_populations. Default: all.
        """
        fda: OpenFDA = ctx.request_context["openfda"]
        return fda.drug_label(drug_name, sections=sections)

    @mcp.tool()
    def openfda_recalls(ctx: Context, drug_name: str, classification: str | None = None, limit: int = 25) -> dict:
        """Get FDA drug recall enforcement reports.

        Args:
            drug_name: Generic drug name.
            classification: Filter by class ("Class I" = most serious, "Class II", "Class III").
            limit: Max results (default 25, max 100).
        """
        fda: OpenFDA = ctx.request_context["openfda"]
        return fda.recalls(drug_name, classification=classification, limit=limit)

    # ── UniProt ────────────────────────────────────────────────────

    @mcp.tool()
    def uniprot_protein(ctx: Context, accession: str | None = None, gene_symbol: str | None = None) -> dict:
        """Get protein details: function, subcellular location, domains, disease associations, GO terms, cross-references.

        Args:
            accession: UniProt accession (e.g. "P38398" for BRCA1).
            gene_symbol: Gene symbol (e.g. "BRCA1"). Resolved to human reviewed entry.
        """
        up: UniProt = ctx.request_context["uniprot"]
        result = up.protein(accession=accession, gene_symbol=gene_symbol)
        if not result:
            return {"error": f"Protein not found: accession={accession}, gene_symbol={gene_symbol}"}
        return result

    @mcp.tool()
    def uniprot_variants(ctx: Context, accession: str | None = None, gene_symbol: str | None = None, limit: int = 50) -> dict:
        """Get known natural variants and mutagenesis data for a protein.

        Args:
            accession: UniProt accession.
            gene_symbol: Gene symbol (e.g. "TP53").
            limit: Max variants (default 50).
        """
        up: UniProt = ctx.request_context["uniprot"]
        result = up.variants(accession=accession, gene_symbol=gene_symbol, limit=limit)
        if not result:
            return {"error": f"Protein not found: accession={accession}, gene_symbol={gene_symbol}"}
        return result

    @mcp.tool()
    def uniprot_search(ctx: Context, query: str, organism: str = "human", limit: int = 10) -> dict:
        """Search UniProt proteins by keyword.

        Args:
            query: Search term (e.g. "kinase", "DNA repair").
            organism: Organism filter (default "human"). Use NCBI taxonomy ID for non-human.
            limit: Max results (default 10, max 50).
        """
        up: UniProt = ctx.request_context["uniprot"]
        results = up.search(query, organism=organism, limit=limit)
        return {"query": query, "count": len(results), "results": results}

    # ── MyGene.info ────────────────────────────────────────────────

    @mcp.tool()
    def gene_info(ctx: Context, gene_id: str) -> dict:
        """Get gene annotations: Entrez/Ensembl/UniProt IDs, KEGG pathways, GO terms, genomic position, aliases.

        Args:
            gene_id: Gene symbol (e.g. "CYP2D6"), Entrez ID (e.g. "1565"), or Ensembl ID (e.g. "ENSG00000100197").
        """
        mg: MyGene = ctx.request_context["mygene"]
        result = mg.gene_info(gene_id)
        if not result:
            return {"error": f"Gene not found: {gene_id}"}
        return result

    @mcp.tool()
    def gene_search(ctx: Context, query: str, species: str = "human", limit: int = 10) -> dict:
        """Search genes by keyword.

        Args:
            query: Search term (e.g. "cytochrome P450", "breast cancer").
            species: Species filter (default "human").
            limit: Max results (default 10, max 50).
        """
        mg: MyGene = ctx.request_context["mygene"]
        results = mg.search(query, species=species, limit=limit)
        return {"query": query, "count": len(results), "results": results}

    return mcp


def main():
    mcp = create_mcp()
    mcp.run()
