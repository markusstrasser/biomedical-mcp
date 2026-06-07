"""Offline unit tests for the DDInter client.

No network — all HTTP calls are monkeypatched. Tests verify:
  1. Drug name resolution via search-source endpoint.
  2. Interaction parsing: severity level mapping, mechanism label extraction.
  3. Graceful empty return when drug is not found in DDInter.
  4. Cache hit avoids a second POST call.
  5. validate() returns True on a successful probe, False on failure.
"""

from __future__ import annotations

import pytest

from biomedical_mcp.cache import Cache
from biomedical_mcp.ddinter import DDInter, _LEVEL_MAP, _mechanism_label


# ── Fixtures ─────────────────────────────────────────────────────────────────


@pytest.fixture
def cache(tmp_path):
    return Cache(tmp_path / "ddinter_test.db")


@pytest.fixture
def client(cache):
    return DDInter(cache)


# ── Canned server responses ───────────────────────────────────────────────────

SEARCH_ASPIRIN = {
    "draw": 1,
    "recordsTotal": 1,
    "recordsFiltered": 1,
    "data": [
        {
            "internal_id": "DDInter20",
            "drug": "Acetylsalicylic acid",
            "drugbank_id": "DB00945",
            "version": "1.0",
            "type": "small molecule",
            "indication": "Pain relief...",
        }
    ],
}

SEARCH_EMPTY = {"draw": 1, "recordsTotal": 0, "recordsFiltered": 0, "data": []}

INTERACTIONS_ASPIRIN = {
    "draw": 1,
    "recordsTotal": 694,
    "recordsFiltered": 694,
    "data": [
        {
            "drug_id": "DDInter900",
            "drug_name": "Ibuprofen",
            "interaction_id": 10165,
            "level": 3,
            "excretion": "0",
            "absorption": "0",
            "others": "0",
            "antagonistic_effect": "1",
            "distribution": "0",
            "metabolism": "0",
            "synergistic_effect": "1",
        },
        {
            "drug_id": "DDInter9",
            "drug_name": "Acalabrutinib",
            "interaction_id": 2728,
            "level": 2,
            "excretion": "0",
            "absorption": "0",
            "others": "0",
            "antagonistic_effect": "0",
            "distribution": "0",
            "metabolism": "1",
            "synergistic_effect": "0",
        },
        {
            "drug_id": "DDInter999",
            "drug_name": "SomeDrug",
            "interaction_id": 9999,
            "level": 1,
            "excretion": "1",
            "absorption": "0",
            "others": "0",
            "antagonistic_effect": "0",
            "distribution": "0",
            "metabolism": "0",
            "synergistic_effect": "0",
        },
        {
            "drug_id": "DDInter888",
            "drug_name": "UnknownDrug",
            "interaction_id": 8888,
            "level": 0,
            "excretion": "0",
            "absorption": "0",
            "others": "0",
            "antagonistic_effect": "0",
            "distribution": "0",
            "metabolism": "0",
            "synergistic_effect": "0",
        },
    ],
}


# ── Helper ────────────────────────────────────────────────────────────────────


def _make_post_mock(responses: list[dict]):
    """Return a mock for _post_datatable that yields responses in order."""
    responses_iter = iter(responses)

    def mock_post(path, *, keyword=None, length=25):
        return next(responses_iter)

    return mock_post


# ── Tests: level map & mechanism label (pure functions) ───────────────────────


class TestPureFunctions:
    def test_level_map_major(self):
        assert _LEVEL_MAP[3] == "Major"

    def test_level_map_moderate(self):
        assert _LEVEL_MAP[2] == "Moderate"

    def test_level_map_minor(self):
        assert _LEVEL_MAP[1] == "Minor"

    def test_level_map_unknown(self):
        assert _LEVEL_MAP[0] == "Unknown"

    def test_mechanism_label_single(self):
        row = {"metabolism": "1", "absorption": "0", "antagonistic_effect": "0",
               "distribution": "0", "excretion": "0", "synergistic_effect": "0", "others": "0"}
        assert _mechanism_label(row) == "metabolism"

    def test_mechanism_label_multiple(self):
        row = {"absorption": "1", "metabolism": "1", "distribution": "0",
               "excretion": "0", "synergistic_effect": "0", "antagonistic_effect": "0", "others": "0"}
        label = _mechanism_label(row)
        assert "absorption" in label
        assert "metabolism" in label

    def test_mechanism_label_none_when_empty(self):
        row = {k: "0" for k in ("absorption", "distribution", "metabolism",
                                "excretion", "synergistic_effect", "antagonistic_effect", "others")}
        assert _mechanism_label(row) is None


# ── Tests: interactions() ─────────────────────────────────────────────────────


class TestInteractions:
    def test_interactions_returns_list(self, client, monkeypatch):
        monkeypatch.setattr(
            client, "_post_datatable",
            _make_post_mock([SEARCH_ASPIRIN, INTERACTIONS_ASPIRIN])
        )
        result = client.interactions("aspirin", limit=10)
        assert isinstance(result, list)
        assert len(result) == 4

    def test_first_result_is_ibuprofen_major(self, client, monkeypatch):
        monkeypatch.setattr(
            client, "_post_datatable",
            _make_post_mock([SEARCH_ASPIRIN, INTERACTIONS_ASPIRIN])
        )
        result = client.interactions("aspirin", limit=10)
        first = result[0]
        assert first["interacting_drug"] == "Ibuprofen"
        assert first["severity"] == "Major"
        assert first["ddinter_id"] == "DDInter900"
        assert first["interaction_id"] == 10165

    def test_severity_moderate(self, client, monkeypatch):
        monkeypatch.setattr(
            client, "_post_datatable",
            _make_post_mock([SEARCH_ASPIRIN, INTERACTIONS_ASPIRIN])
        )
        result = client.interactions("aspirin", limit=10)
        second = result[1]
        assert second["severity"] == "Moderate"
        assert second["interacting_drug"] == "Acalabrutinib"

    def test_severity_minor(self, client, monkeypatch):
        monkeypatch.setattr(
            client, "_post_datatable",
            _make_post_mock([SEARCH_ASPIRIN, INTERACTIONS_ASPIRIN])
        )
        result = client.interactions("aspirin", limit=10)
        third = result[2]
        assert third["severity"] == "Minor"

    def test_severity_unknown(self, client, monkeypatch):
        monkeypatch.setattr(
            client, "_post_datatable",
            _make_post_mock([SEARCH_ASPIRIN, INTERACTIONS_ASPIRIN])
        )
        result = client.interactions("aspirin", limit=10)
        fourth = result[3]
        assert fourth["severity"] == "Unknown"

    def test_mechanism_parsed(self, client, monkeypatch):
        monkeypatch.setattr(
            client, "_post_datatable",
            _make_post_mock([SEARCH_ASPIRIN, INTERACTIONS_ASPIRIN])
        )
        result = client.interactions("aspirin", limit=10)
        # Ibuprofen row has antagonistic_effect + synergistic_effect
        assert result[0]["mechanism"] is not None
        assert "antagonistic" in result[0]["mechanism"]
        # Acalabrutinib row has only metabolism
        assert result[1]["mechanism"] == "metabolism"
        # Unknown row has no mechanism flags
        assert result[3]["mechanism"] is None

    def test_not_found_returns_empty(self, client, monkeypatch):
        monkeypatch.setattr(
            client, "_post_datatable",
            _make_post_mock([SEARCH_EMPTY])
        )
        result = client.interactions("nonexistent_drug_xyz")
        assert result == []

    def test_cache_avoids_second_call(self, client, monkeypatch):
        call_count = {"n": 0}

        def counting_post(path, *, keyword=None, length=25):
            call_count["n"] += 1
            if call_count["n"] == 1:
                return SEARCH_ASPIRIN
            return INTERACTIONS_ASPIRIN

        monkeypatch.setattr(client, "_post_datatable", counting_post)

        # First call: 2 POSTs (search + interactions)
        client.interactions("aspirin", limit=10)
        assert call_count["n"] == 2

        # Second call: 0 POSTs (both search resolve and interactions are cached)
        client.interactions("aspirin", limit=10)
        assert call_count["n"] == 2  # no additional calls


# ── Tests: validate() ─────────────────────────────────────────────────────────


class TestValidate:
    def test_validate_true_when_data_present(self, client, monkeypatch):
        monkeypatch.setattr(
            client, "_post_datatable",
            _make_post_mock([SEARCH_ASPIRIN])
        )
        assert client.validate() is True

    def test_validate_false_on_empty(self, client, monkeypatch):
        monkeypatch.setattr(
            client, "_post_datatable",
            _make_post_mock([SEARCH_EMPTY])
        )
        assert client.validate() is False

    def test_validate_false_on_exception(self, client, monkeypatch):
        def raise_exc(**_):
            raise ConnectionError("network down")

        monkeypatch.setattr(client, "_post_datatable", raise_exc)
        assert client.validate() is False


# ── caught-red-handed (close review) ─────────────────────────────────────────

def test_interactions_handles_null_or_empty_level(tmp_path):
    """A null/empty `level` from the undocumented endpoint must not crash parsing."""
    from biomedical_mcp.cache import Cache
    from biomedical_mcp.ddinter import DDInter
    c = DDInter(Cache(tmp_path / "c.db"))
    c._resolve_drug_id = lambda name: "DDInter1"  # type: ignore
    c._post_datatable = lambda path, length=500: {  # type: ignore
        "data": [
            {"drug_name": "warfarin", "level": None},
            {"drug_name": "aspirin", "level": ""},
            {"drug_name": "ibuprofen", "level": 3},
        ]
    }
    rows = c.interactions("testdrug")
    sev = {r["interacting_drug"]: r["severity"] for r in rows}
    assert sev["warfarin"] == "Unknown"
    assert sev["aspirin"] == "Unknown"
    assert sev["ibuprofen"] == "Major"
