"""Live API smoke tests for new domains — one request per API to verify connectivity."""

import tempfile
from pathlib import Path

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
