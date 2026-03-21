"""Live API smoke tests for new domains — one request per API to verify connectivity."""

import pytest

from biomedical_mcp.cache import Cache
from biomedical_mcp.gnomad import GnomAD
from biomedical_mcp.panelapp import PanelApp
from biomedical_mcp.hpo import HPO
from biomedical_mcp.monarch import Monarch
from biomedical_mcp.gtex import GTEx
from biomedical_mcp.hgnc import HGNC
from biomedical_mcp.gwas_catalog import GWASCatalog
from biomedical_mcp.litvar import LitVar
from biomedical_mcp.interpro import InterPro
from biomedical_mcp.pdb import PDB


@pytest.fixture
def cache(tmp_path):
    return Cache(tmp_path / "test_cache.db")


# ── gnomAD ────────────────────────────────────────────────────


class TestGnomAD:
    def test_variant_frequency(self, cache):
        gn = GnomAD(cache)
        result = gn.variant_frequency("1-11796321-G-A")  # rs1801133, MTHFR C677T — common
        assert "error" not in result
        assert "genome" in result or "exome" in result

    def test_gene_constraint(self, cache):
        gn = GnomAD(cache)
        result = gn.gene_constraint("BRCA1")
        assert "error" not in result
        # Should have constraint metrics
        assert result.get("pLI") is not None or result.get("oe_lof") is not None

    def test_validate(self, cache):
        gn = GnomAD(cache)
        assert gn.validate() is True


# ── PanelApp ──────────────────────────────────────────────────


class TestPanelApp:
    def test_gene_panels(self, cache):
        pa = PanelApp(cache)
        result = pa.gene_panels("BRCA1", confidence="all")
        assert result["gene_symbol"] == "BRCA1"
        assert result["count"] > 0

    def test_panel_search(self, cache):
        pa = PanelApp(cache)
        result = pa.panel_search("epilepsy")
        assert result["count"] > 0

    def test_validate(self, cache):
        pa = PanelApp(cache)
        assert pa.validate() is True


# ── HPO ───────────────────────────────────────────────────────


class TestHPO:
    def test_search(self, cache):
        hpo = HPO(cache)
        result = hpo.search("seizure", limit=5)
        assert result["count"] > 0

    def test_term(self, cache):
        hpo = HPO(cache)
        result = hpo.term("HP:0001250")  # Seizure
        assert result["hpo_id"] == "HP:0001250"
        assert result["name"] is not None

    def test_validate(self, cache):
        hpo = HPO(cache)
        assert hpo.validate() is True


# ── Monarch ───────────────────────────────────────────────────


class TestMonarch:
    def test_search(self, cache):
        mon = Monarch(cache)
        result = mon.search("BRCA1", limit=3)
        assert result["count"] > 0

    def test_validate(self, cache):
        mon = Monarch(cache)
        assert mon.validate() is True


# ── GTEx ──────────────────────────────────────────────────────


class TestGTEx:
    def test_gene_expression(self, cache):
        gt = GTEx(cache)
        result = gt.gene_expression(gene_symbol="CYP2D6")
        # GTEx API may return different shapes; just verify no crash
        assert "gene" in result

    def test_validate(self, cache):
        gt = GTEx(cache)
        assert gt.validate() is True


# ── HGNC ──────────────────────────────────────────────────────


class TestHGNC:
    def test_gene_names(self, cache):
        hg = HGNC(cache)
        result = hg.gene_names("BRCA1")
        assert result["symbol"] == "BRCA1"
        assert result["hgnc_id"] is not None
        assert result["name"] is not None

    def test_search(self, cache):
        hg = HGNC(cache)
        result = hg.search("cytochrome P450")
        assert result["count"] > 0

    def test_validate(self, cache):
        hg = HGNC(cache)
        assert hg.validate() is True


# ── GWAS Catalog ──────────────────────────────────────────────


class TestGWASCatalog:
    def test_trait_search(self, cache):
        gc = GWASCatalog(cache)
        result = gc.trait_search("diabetes", limit=5)
        assert result["count"] > 0

    def test_validate(self, cache):
        gc = GWASCatalog(cache)
        assert gc.validate() is True


# ── LitVar2 ──────────────────────────────────────────────────


class TestLitVar:
    def test_variant_publications(self, cache):
        lv = LitVar(cache)
        result = lv.variant_publications("rs1801133")
        assert result["variant_id"] == "rs1801133"
        # Common variant should have publications
        assert result["publication_count"] > 0

    def test_validate(self, cache):
        lv = LitVar(cache)
        assert lv.validate() is True


# ── InterPro ──────────────────────────────────────────────────


class TestInterPro:
    def test_validate(self, cache):
        ip = InterPro(cache)
        assert ip.validate() is True


# ── PDB ───────────────────────────────────────────────────────


class TestPDB:
    def test_structure_detail(self, cache):
        p = PDB(cache)
        result = p.structure_detail("1TUP")  # TP53 DNA-binding domain
        assert result["pdb_id"] == "1TUP"
        assert result.get("title") is not None

    def test_validate(self, cache):
        p = PDB(cache)
        assert p.validate() is True


# ── ISBT Blood Groups ────────────────────────────────────────


class TestISBT:
    def test_systems(self, cache):
        from biomedical_mcp.isbt import ISBT
        isbt = ISBT(cache)
        result = isbt.systems()
        assert result["count"] >= 48  # 48 blood group systems minimum
        assert any(s["symbol"] == "ABO" for s in result["systems"])

    def test_alleles(self, cache):
        from biomedical_mcp.isbt import ISBT
        isbt = ISBT(cache)
        result = isbt.alleles("ABO")
        assert result["count"] > 100  # ABO has 200+ alleles
        assert any(a["isbt_allele"].startswith("ABO*") for a in result["alleles"])

    def test_search_alleles(self, cache):
        from biomedical_mcp.isbt import ISBT
        isbt = ISBT(cache)
        result = isbt.search_alleles(system_symbol="FY")
        assert result["count"] > 10  # FY has 30+ alleles

    def test_variant_lookup_exonic(self, cache):
        from biomedical_mcp.isbt import ISBT
        isbt = ISBT(cache)
        result = isbt.variant_lookup(rsid="rs8176719")  # ABO O allele 261delG
        assert result["count"] >= 1

    def test_validate(self, cache):
        from biomedical_mcp.isbt import ISBT
        isbt = ISBT(cache)
        assert isbt.validate() is True


# ── Orphanet ─────────────────────────────────────────────────


class TestOrphanet:
    def test_gene_diseases(self, cache):
        from biomedical_mcp.orphanet import Orphanet
        orp = Orphanet(cache)
        result = orp.gene_diseases("BRCA1")
        assert result["count"] > 0
        assert any("susceptibility" in (d.get("association_type") or "").lower()
                    or "disease-causing" in (d.get("association_type") or "").lower()
                    for d in result["diseases"])

    def test_natural_history(self, cache):
        from biomedical_mcp.orphanet import Orphanet
        orp = Orphanet(cache)
        result = orp.natural_history(558)  # Marfan syndrome
        assert result["name"] == "Marfan syndrome"
        assert "Autosomal dominant" in result["inheritance"]

    def test_omim_lookup(self, cache):
        from biomedical_mcp.orphanet import Orphanet
        orp = Orphanet(cache)
        result = orp.search_by_omim("154700")  # Marfan
        assert result["count"] > 0

    def test_omim_404_graceful(self, cache):
        from biomedical_mcp.orphanet import Orphanet
        orp = Orphanet(cache)
        result = orp.search_by_omim("999999999")  # Nonexistent
        assert result["count"] == 0
        assert "note" in result

    def test_validate(self, cache):
        from biomedical_mcp.orphanet import Orphanet
        orp = Orphanet(cache)
        assert orp.validate() is True


# ── ClinGen ──────────────────────────────────────────────────


class TestClinGen:
    def test_gene_validity(self, cache):
        from biomedical_mcp.clingen import ClinGen
        cg = ClinGen(cache)
        result = cg.gene_validity("BRCA1")
        assert result["count"] >= 1
        assert any(c["classification"] == "Definitive" for c in result["curations"])

    def test_gene_dosage(self, cache):
        from biomedical_mcp.clingen import ClinGen
        cg = ClinGen(cache)
        result = cg.gene_dosage("BRCA1")
        assert result["count"] >= 1

    def test_disease_search(self, cache):
        from biomedical_mcp.clingen import ClinGen
        cg = ClinGen(cache)
        result = cg.disease_validity("cardiomyopathy")
        assert result["count"] > 10  # Many cardiomyopathy gene curations


# ── Composite ─────────────────────────────────────────────────


class TestComposite:
    def test_gene_dossier(self, cache):
        """Test gene_dossier compound tool."""
        from biomedical_mcp.domains.composite import create_server as create_composite
        # Just verify the module imports and server creates without error
        server = create_composite(cache)
        assert server is not None
