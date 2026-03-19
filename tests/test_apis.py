"""Live API smoke tests — one request per API to verify connectivity and response shape."""

import tempfile
from pathlib import Path

import pytest

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


@pytest.fixture
def cache(tmp_path):
    return Cache(tmp_path / "test_cache.db")


# ── Open Targets ──────────────────────────────────────────────


class TestOpenTargets:
    def test_search(self, cache):
        ot = OpenTargets(cache)
        hits = ot.search("BRCA1", entity_type="target", limit=3)
        assert len(hits) > 0
        assert hits[0]["id"].startswith("ENSG")

    def test_target_info(self, cache):
        ot = OpenTargets(cache)
        info = ot.target_info("ENSG00000012048")  # BRCA1
        assert info is not None
        assert info["approvedSymbol"] == "BRCA1"

    def test_disease_associations(self, cache):
        ot = OpenTargets(cache)
        result = ot.disease_associations("ENSG00000012048", limit=5)
        assert result["target"] == "BRCA1"
        assert len(result["associations"]) > 0

    def test_resolve_target(self, cache):
        ot = OpenTargets(cache)
        eid = ot.resolve_target(gene_symbol="TP53")
        assert eid is not None
        assert eid.startswith("ENSG")

    def test_drug_info(self, cache):
        ot = OpenTargets(cache)
        info = ot.drug_info("CHEMBL25")  # aspirin
        assert info is not None
        assert info["name"].lower() == "aspirin"


# ── ChEMBL ────────────────────────────────────────────────────


class TestChEMBL:
    def test_resolve_compound(self, cache):
        ch = ChEMBL(cache)
        cid = ch.resolve_compound(name="ibuprofen")
        assert cid is not None
        assert cid.startswith("CHEMBL")

    def test_compound(self, cache):
        ch = ChEMBL(cache)
        result = ch.compound("CHEMBL521")  # ibuprofen
        assert result["name"] == "IBUPROFEN"
        assert float(result["max_phase"]) == 4.0

    def test_mechanism(self, cache):
        ch = ChEMBL(cache)
        mechs = ch.mechanism("CHEMBL521")
        assert len(mechs) > 0
        assert mechs[0]["mechanism"] is not None

    def test_target_search(self, cache):
        ch = ChEMBL(cache)
        result = ch.target(gene_symbol="EGFR")
        assert result is not None
        assert result["chembl_id"].startswith("CHEMBL")


# ── ClinicalTrials.gov ────────────────────────────────────────


class TestClinicalTrials:
    def test_search(self, cache):
        ct = ClinicalTrials(cache)
        result = ct.search(condition="breast cancer", limit=5)
        assert result["total"] > 0
        assert len(result["studies"]) > 0
        assert result["studies"][0]["nct_id"].startswith("NCT")

    def test_trial_detail(self, cache):
        ct = ClinicalTrials(cache)
        # Search for a real trial first
        search = ct.search(condition="diabetes", limit=1)
        nct_id = search["studies"][0]["nct_id"]
        detail = ct.trial_detail(nct_id)
        assert detail["nct_id"] == nct_id
        assert "eligibility" in detail

    def test_stats(self, cache):
        ct = ClinicalTrials(cache)
        result = ct.stats("lung cancer")
        assert result["total"] > 0
        assert len(result["by_phase"]) > 0


# ── ICD-10 ────────────────────────────────────────────────────


class TestICD10:
    def test_search(self, cache):
        icd = ICD10(cache)
        results = icd.search("diabetes", limit=5)
        assert len(results) > 0
        assert "code" in results[0]
        assert "description" in results[0]

    def test_lookup(self, cache):
        icd = ICD10(cache)
        result = icd.lookup("E11")
        assert result is not None
        assert "diabetes" in result["description"].lower() or "E11" in result["code"]


# ── NPI ────────────────────────────────────────────────────────


class TestNPI:
    def test_search(self, cache):
        npi_client = NPI(cache)
        results = npi_client.search(specialty="Cardiology", state="CA", limit=3)
        assert len(results) > 0
        assert results[0]["npi"] is not None

    def test_lookup(self, cache):
        npi_client = NPI(cache)
        # Search first to get a real NPI
        results = npi_client.search(specialty="Internal Medicine", state="NY", limit=1)
        if results:
            npi_num = results[0]["npi"]
            detail = npi_client.lookup(npi_num)
            assert detail is not None
            assert detail["npi"] == npi_num


# ── MyVariant.info ─────────────────────────────────────────────


class TestMyVariant:
    def test_lookup(self, cache):
        mv = MyVariant(cache)
        result = mv.lookup("rs1805007")  # MC1R red hair variant — multi-allelic
        # rsIDs with multiple alt alleles return a list
        if isinstance(result, list):
            assert len(result) > 0
            first = result[0]
        else:
            first = result
        assert first["variant_id"] is not None
        assert any(k in first for k in ("clinvar", "gnomad_exome", "cadd", "dbsnp"))

    def test_clinvar_search(self, cache):
        mv = MyVariant(cache)
        results = mv.clinvar_search("BRCA1", significance="Pathogenic", limit=5)
        assert len(results) > 0
        assert results[0]["variant_id"] is not None

    def test_batch(self, cache):
        mv = MyVariant(cache)
        results = mv.batch(["rs1805007", "rs1800497"])
        # Multi-allelic rsIDs expand to multiple entries
        assert len(results) >= 2
        assert all(r["variant_id"] is not None for r in results)


# ── OpenFDA ────────────────────────────────────────────────────


class TestOpenFDA:
    def test_adverse_events(self, cache):
        fda = OpenFDA(cache)
        result = fda.adverse_events("ibuprofen", limit=3)
        assert result["drug"] == "ibuprofen"
        assert result["total"] > 0
        assert len(result["events"]) > 0

    def test_drug_label(self, cache):
        fda = OpenFDA(cache)
        result = fda.drug_label("metformin")
        assert "error" not in result
        assert result["drug"] == "metformin"
        assert len(result["sections"]) > 0

    def test_recalls(self, cache):
        fda = OpenFDA(cache)
        # Metformin has had NDMA-related recalls
        result = fda.recalls("metformin", limit=3)
        assert result["drug"] == "metformin"
        # May or may not have recalls, just check shape
        assert "recalls" in result


# ── UniProt ────────────────────────────────────────────────────


class TestUniProt:
    def test_protein_by_gene(self, cache):
        up = UniProt(cache)
        result = up.protein(gene_symbol="BRCA1")
        assert result is not None
        assert result["gene"] == "BRCA1"
        assert result["function"] is not None
        assert len(result["domains"]) > 0

    def test_variants(self, cache):
        up = UniProt(cache)
        result = up.variants(gene_symbol="TP53")
        assert result is not None
        assert result["variant_count"] > 0

    def test_search(self, cache):
        up = UniProt(cache)
        results = up.search("kinase", limit=5)
        assert len(results) > 0
        assert results[0]["accession"] is not None


# ── MyGene.info ────────────────────────────────────────────────


class TestMyGene:
    def test_gene_info_by_symbol(self, cache):
        mg = MyGene(cache)
        result = mg.gene_info("CYP2D6")
        assert result is not None
        assert result["symbol"] == "CYP2D6"
        assert result["entrez_id"] is not None
        assert result["ensembl_gene"] is not None

    def test_gene_info_kegg(self, cache):
        mg = MyGene(cache)
        result = mg.gene_info("CYP2D6")
        assert result is not None
        # CYP2D6 should have KEGG pathways
        assert len(result["kegg_pathways"]) > 0

    def test_search(self, cache):
        mg = MyGene(cache)
        results = mg.search("cytochrome P450", limit=5)
        assert len(results) > 0
        assert results[0]["symbol"] is not None
