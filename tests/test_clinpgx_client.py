"""Offline tests for the ClinPGx client.

Monkeypatches _cached_get to return canned PostgREST rows.
No network calls are made.
"""

from __future__ import annotations

import pytest

from biomedical_mcp.cache import Cache
from biomedical_mcp.clinpgx import ClinPGx, _level_rank
from biomedical_mcp.errors import NotFoundError, SourceUnavailableError


# ── Canned fixture data (PostgREST shapes from live API probing) ──────────────

PAIR_VIEW_CODEINE = [
    {
        "pairid": 1001,
        "drugid": "RxNorm:2670",
        "drugname": "codeine",
        "genesymbol": "CYP2D6",
        "guidelinename": "CYP2D6 and Codeine",
        "guidelineurl": "https://www.clinpgx.org/guideline/PA166251464",
        "cpiclevel": "A",
        "clinpgxlevel": "1A",
        "pgxtesting": "Actionable PGx",
        "pmids": ["28407481"],
        "usedforrecommendation": "Yes",
        "provisional": False,
    },
    {
        "pairid": 1002,
        "drugid": "RxNorm:2670",
        "drugname": "codeine",
        "genesymbol": "CYP3A4",
        "guidelinename": None,
        "guidelineurl": None,
        "cpiclevel": "C",
        "clinpgxlevel": None,
        "pgxtesting": "Informative PGx",
        "pmids": None,
        "usedforrecommendation": "No",
        "provisional": True,
    },
]

PAIR_VIEW_CYP2D6 = [
    {
        "pairid": 2001,
        "drugid": "RxNorm:7531",
        "drugname": "nortriptyline",
        "genesymbol": "CYP2D6",
        "guidelinename": "CYP2D6, CYP2C19 and Tricyclic Antidepressants",
        "guidelineurl": "https://www.clinpgx.org/guideline/PA166251445",
        "cpiclevel": "A",
        "clinpgxlevel": "1A",
        "pgxtesting": "Informative PGx",
        "pmids": ["23486447"],
        "usedforrecommendation": "Yes",
        "provisional": False,
    },
    {
        "pairid": 2002,
        "drugid": "RxNorm:6813",
        "drugname": "methadone",
        "genesymbol": "CYP2D6",
        "guidelinename": "CYP2D6, OPRM1, COMT, and Opioids",
        "guidelineurl": "https://www.clinpgx.org/guideline/PA166251454",
        "cpiclevel": "C",
        "clinpgxlevel": "3",
        "pgxtesting": None,
        "pmids": ["33387367"],
        "usedforrecommendation": "No",
        "provisional": False,
    },
    {
        "pairid": 2003,
        "drugid": "RxNorm:77492",
        "drugname": "tamsulosin",
        "genesymbol": "CYP2D6",
        "guidelinename": None,
        "guidelineurl": None,
        "cpiclevel": "B/C",
        "clinpgxlevel": None,
        "pgxtesting": "Actionable PGx",
        "pmids": None,
        "usedforrecommendation": "n/a",
        "provisional": True,
    },
]

RECS_CODEINE = [
    {
        "recommendationid": 5001,
        "lookupkey": {"CYP2D6": "0.0"},
        "drugname": "codeine",
        "guidelinename": "CYP2D6 and Codeine",
        "guidelineurl": "https://www.clinpgx.org/guideline/PA166251464",
        "implications": {"CYP2D6": "Reduced morphine formation; reduced analgesia"},
        "drugrecommendation": "Avoid codeine use. Use an alternative analgesic.",
        "classification": "Strong",
        "phenotypes": {"CYP2D6": "Poor Metabolizer"},
        "activityscore": {"CYP2D6": "0.0"},
        "population": "general",
        "comments": "n/a",
    },
    {
        "recommendationid": 5002,
        "lookupkey": {"CYP2D6": "3.0"},
        "drugname": "codeine",
        "guidelinename": "CYP2D6 and Codeine",
        "guidelineurl": "https://www.clinpgx.org/guideline/PA166251464",
        "implications": {"CYP2D6": "Increased morphine formation; risk of toxicity"},
        "drugrecommendation": "Avoid codeine. Use alternative opioid with caution.",
        "classification": "Strong",
        "phenotypes": {"CYP2D6": "Ultrarapid Metabolizer"},
        "activityscore": {"CYP2D6": "3.0"},
        "population": "general",
        "comments": "Especially dangerous in nursing mothers or neonates.",
    },
]

RECS_NORTRIPTYLINE: list[dict] = []
RECS_METHADONE: list[dict] = []
RECS_TAMSULOSIN: list[dict] = []


# ── Fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture
def cache(tmp_path):
    return Cache(tmp_path / "test_cache.db")


@pytest.fixture
def client(cache):
    return ClinPGx(cache)


# ── _level_rank unit tests ────────────────────────────────────────────────────


class TestLevelRank:
    def test_order(self):
        assert _level_rank("A") < _level_rank("B")
        assert _level_rank("B") < _level_rank("B/C")
        assert _level_rank("B/C") < _level_rank("C")
        assert _level_rank("C") < _level_rank("D")

    def test_null_sorts_last(self):
        assert _level_rank(None) > _level_rank("D")

    def test_unknown_sorts_last(self):
        assert _level_rank("X") > _level_rank("D")

    def test_whitespace_tolerant(self):
        assert _level_rank(" A ") == _level_rank("A")


# ── pgx_for_drug ─────────────────────────────────────────────────────────────


class TestPgxForDrug:
    def _patch(self, monkeypatch, client):
        """Patch _cached_get to return canned data based on (prefix, path, params)."""
        def fake_cached_get(prefix, path, params=None, **kw):
            if prefix == "pairs" and params and "codeine" in params.get("drugname", ""):
                return PAIR_VIEW_CODEINE, False
            if prefix == "recs" and params and "codeine" in params.get("drugname", ""):
                return RECS_CODEINE, False
            return [], False

        monkeypatch.setattr(client, "_cached_get", fake_cached_get)

    def test_returns_pairs(self, client, monkeypatch):
        self._patch(monkeypatch, client)
        result = client.pgx_for_drug("codeine")
        assert result["drug"] == "codeine"
        assert result["pair_count"] == 2
        assert result["source"] == "cpic"

    def test_sorted_by_cpic_level(self, client, monkeypatch):
        self._patch(monkeypatch, client)
        result = client.pgx_for_drug("codeine")
        levels = [p["cpic_level"] for p in result["pairs"]]
        # A should come before C
        assert levels[0] == "A"
        assert levels[1] == "C"

    def test_recommendations_attached(self, client, monkeypatch):
        self._patch(monkeypatch, client)
        result = client.pgx_for_drug("codeine")
        # CYP2D6 pair (level A) should have recommendations
        cyp2d6_pair = next(p for p in result["pairs"] if p["gene"] == "CYP2D6")
        assert len(cyp2d6_pair["recommendations"]) == 2

    def test_recommendation_fields(self, client, monkeypatch):
        self._patch(monkeypatch, client)
        result = client.pgx_for_drug("codeine")
        cyp2d6_pair = next(p for p in result["pairs"] if p["gene"] == "CYP2D6")
        rec = cyp2d6_pair["recommendations"][0]
        assert "phenotypes" in rec
        assert "recommendation" in rec
        assert "classification" in rec
        assert "implications" in rec
        assert "activity_score" in rec

    def test_comments_na_becomes_none(self, client, monkeypatch):
        self._patch(monkeypatch, client)
        result = client.pgx_for_drug("codeine")
        cyp2d6_pair = next(p for p in result["pairs"] if p["gene"] == "CYP2D6")
        # First rec has comments="n/a" → should be None
        assert cyp2d6_pair["recommendations"][0]["comments"] is None
        # Second rec has real comment → should be present
        assert cyp2d6_pair["recommendations"][1]["comments"] is not None

    def test_guideline_url_present(self, client, monkeypatch):
        self._patch(monkeypatch, client)
        result = client.pgx_for_drug("codeine")
        cyp2d6_pair = next(p for p in result["pairs"] if p["gene"] == "CYP2D6")
        assert cyp2d6_pair["guideline_url"].startswith("https://")

    def test_not_found_raises(self, client, monkeypatch):
        monkeypatch.setattr(client, "_cached_get", lambda *a, **kw: ([], False))
        with pytest.raises(NotFoundError):
            client.pgx_for_drug("nonexistent_drug_xyz")

    def test_network_error_raises_source_unavailable(self, client, monkeypatch):
        import httpx

        def raise_network(*a, **kw):
            raise httpx.ConnectError("connect failed")

        monkeypatch.setattr(client, "_cached_get", raise_network)
        with pytest.raises(SourceUnavailableError):
            client.pgx_for_drug("codeine")


# ── pgx_for_gene ──────────────────────────────────────────────────────────────


class TestPgxForGene:
    def _patch(self, monkeypatch, client):
        recs_map = {
            "nortriptyline": RECS_NORTRIPTYLINE,
            "methadone": RECS_METHADONE,
            "tamsulosin": RECS_TAMSULOSIN,
        }

        def fake_cached_get(prefix, path, params=None, **kw):
            if prefix == "pairs" and params and "CYP2D6" in params.get("genesymbol", ""):
                return PAIR_VIEW_CYP2D6, False
            if prefix == "recs":
                dn = (params or {}).get("drugname", "eq.").replace("eq.", "")
                return recs_map.get(dn, []), False
            return [], False

        monkeypatch.setattr(client, "_cached_get", fake_cached_get)

    def test_returns_pairs(self, client, monkeypatch):
        self._patch(monkeypatch, client)
        result = client.pgx_for_gene("CYP2D6")
        assert result["gene"] == "CYP2D6"
        assert result["pair_count"] == 3
        assert result["source"] == "cpic"

    def test_sorted_by_cpic_level(self, client, monkeypatch):
        self._patch(monkeypatch, client)
        result = client.pgx_for_gene("CYP2D6")
        levels = [p["cpic_level"] for p in result["pairs"]]
        # A first, then B/C, then C
        assert levels[0] == "A"
        assert levels[1] == "B/C"
        assert levels[2] == "C"

    def test_gene_normalised_to_uppercase(self, client, monkeypatch):
        """Lowercase input should be normalised before hitting the API."""
        calls = []

        def fake_cached_get(prefix, path, params=None, **kw):
            calls.append(params)
            if prefix == "pairs":
                return PAIR_VIEW_CYP2D6, False
            return [], False

        monkeypatch.setattr(client, "_cached_get", fake_cached_get)
        client.pgx_for_gene("cyp2d6")
        pair_call = next(c for c in calls if c and "genesymbol" in c)
        assert pair_call["genesymbol"] == "eq.CYP2D6"

    def test_not_found_raises(self, client, monkeypatch):
        monkeypatch.setattr(client, "_cached_get", lambda *a, **kw: ([], False))
        with pytest.raises(NotFoundError):
            client.pgx_for_gene("NOTAREAL_GENE")

    def test_network_error_raises_source_unavailable(self, client, monkeypatch):
        import httpx

        def raise_network(*a, **kw):
            raise httpx.ConnectError("connect failed")

        monkeypatch.setattr(client, "_cached_get", raise_network)
        with pytest.raises(SourceUnavailableError):
            client.pgx_for_gene("CYP2D6")

    def test_pair_fields_complete(self, client, monkeypatch):
        self._patch(monkeypatch, client)
        result = client.pgx_for_gene("CYP2D6")
        pair = result["pairs"][0]  # nortriptyline, level A
        assert pair["gene"] == "CYP2D6"
        assert pair["drug"] == "nortriptyline"
        assert pair["cpic_level"] == "A"
        assert pair["guideline_url"].startswith("https://")
        assert isinstance(pair["pmids"], list)
        assert isinstance(pair["provisional"], bool)
        assert isinstance(pair["recommendations"], list)
