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

    def test_variant_associations(self, cache):
        gc = GWASCatalog(cache)
        result = gc.variant_associations("rs7903146", limit=5)
        assert result["association_count"] > 0
        assert result["associations"][0]["study"].startswith("GCST")

    def test_gene_associations(self, cache):
        gc = GWASCatalog(cache)
        result = gc.gene_associations("FTO", limit=3)
        assert result["association_count"] > 0
        assert "FTO" in result["associations"][0]["genes"]


V2_ASSOC_PAGE = {
    "_embedded": {"associations": [{
        "association_id": 226328771, "pvalue_mantissa": 1, "pvalue_exponent": -14,
        "or_value": "1.12", "beta": "-", "accession_id": "GCST90476468",
        "pubmed_id": "39024449", "mapped_genes": ["FTO"],
        "snp_effect_allele": ["rs1421085-T"],
        "efo_traits": [{"efo_id": "EFO_0004338", "efo_trait": "body weight"}],
    }]},
    "page": {"size": 1, "totalElements": 1, "totalPages": 1, "number": 0},
}
V2_TRAIT_PAGE = {"_embedded": {"efo_traits": [{
    "efo_trait": "asthma", "efo_id": "MONDO_0004979",
    "uri": "http://purl.obolibrary.org/obo/MONDO_0004979",
}]}}


class TestGWASCatalogV2Shapes:
    """Offline: v2 response shapes map onto the stable output keys."""

    def _patch(self, monkeypatch, gc, pub_date="2024-07-19"):
        calls = []

        def fake_get(path, params=None, headers=None):
            calls.append((path, params))
            if path == "/associations":
                return V2_ASSOC_PAGE
            if path == "/efo-traits":
                return V2_TRAIT_PAGE
            if path.startswith("/publications/"):
                return {"pubmed_id": "39024449", "publication_date": pub_date}
            raise AssertionError(path)

        monkeypatch.setattr(gc, "_get", fake_get)
        return calls

    def test_variant_associations(self, cache, monkeypatch):
        gc = GWASCatalog(cache)
        calls = self._patch(monkeypatch, gc)
        result = gc.variant_associations("rs1421085", limit=5)
        assert calls[0] == ("/associations", {"rs_id": "rs1421085", "size": 5})
        a = result["associations"][0]
        assert a == {
            "pvalue": 1, "pvalue_exponent": -14, "risk_allele": "rs1421085-T",
            "or_beta": "1.12", "genes": ["FTO"], "trait": "body weight",
            "traits": ["body weight"], "study": "GCST90476468", "pubmed_id": "39024449",
        }

    def test_gene_associations(self, cache, monkeypatch):
        gc = GWASCatalog(cache)
        calls = self._patch(monkeypatch, gc)
        result = gc.gene_associations("FTO", limit=5)
        assert calls[0][1]["mapped_gene"] == "FTO"
        assert result["association_count"] == 1

    def test_empty_page_is_not_found(self, cache, monkeypatch):
        gc = GWASCatalog(cache)
        monkeypatch.setattr(gc, "_get", lambda path, params=None, headers=None: {"page": {"totalElements": 0}})
        result = gc.variant_associations("rs0", limit=5)
        assert result["association_count"] == 0 and "not found" in result["note"]

    def test_trait_search(self, cache, monkeypatch):
        gc = GWASCatalog(cache)
        self._patch(monkeypatch, gc)
        result = gc.trait_search("asthma", limit=5)
        assert result["traits"] == [{"trait": "asthma", "short_form": "MONDO_0004979",
                                     "uri": "http://purl.obolibrary.org/obo/MONDO_0004979"}]

    def test_new_for_variants_date_filter(self, cache, monkeypatch):
        gc = GWASCatalog(cache)
        self._patch(monkeypatch, gc, pub_date="2024-07-19")
        kept = gc.new_for_variants(["rs1421085"], since_date="2024-01-01")
        assert kept["variants"]["rs1421085"]["associations"][0]["publication_date"] == "2024-07-19"

    def test_new_for_variants_drops_older(self, cache, monkeypatch):
        gc = GWASCatalog(cache)
        self._patch(monkeypatch, gc, pub_date="2020-01-01")
        dropped = gc.new_for_variants(["rs1421085"], since_date="2024-01-01")
        assert dropped["variants"]["rs1421085"]["association_count"] == 0


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
