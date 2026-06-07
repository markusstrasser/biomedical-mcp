"""Offline tests for the Somatic client (CIViC + OncoKB).

No network: _post (CIViC GraphQL) and _get (OncoKB REST) are monkeypatched.
All assertions use canned response payloads that mirror real API shapes.
"""

from __future__ import annotations

import os
import pytest

from biomedical_mcp.cache import Cache
from biomedical_mcp.somatic import Somatic, _parse_combined, _parse_civic_evidence, _parse_oncokb_response


# ── _parse_combined ────────────────────────────────────────────────────────────

def test_parse_combined_two_args():
    gene, change = _parse_combined("BRAF", "V600E")
    assert gene == "BRAF"
    assert change == "V600E"


def test_parse_combined_single_string():
    gene, change = _parse_combined("BRAF V600E", None)
    assert gene == "BRAF"
    assert change == "V600E"


def test_parse_combined_strips_whitespace():
    gene, change = _parse_combined("  EGFR  ", "  L858R  ")
    assert gene == "EGFR"
    assert change == "L858R"


def test_parse_combined_invalid_raises():
    with pytest.raises(ValueError, match="Cannot parse"):
        _parse_combined("BRAF", None)   # single word, no space


# ── _parse_civic_evidence ──────────────────────────────────────────────────────

_CIVIC_NODES = [
    {
        "id": 1409,
        "evidenceType": "PREDICTIVE",
        "evidenceLevel": "A",
        "evidenceDirection": "SUPPORTS",
        "significance": "SENSITIVITYRESPONSE",
        "variantOrigin": "SOMATIC",
        "status": "ACCEPTED",
        "description": "Chapman et al. phase III trial...",
        "disease": {"name": "Skin Melanoma"},
        "therapies": [{"name": "Vemurafenib"}],
        "source": {
            "citation": "Chapman et al., 2011",
            "citationId": "21639808",
            "sourceType": "PUBMED",
            "link": "https://pubmed.ncbi.nlm.nih.gov/21639808/",
        },
    },
    {
        "id": 9851,
        "evidenceType": "PREDICTIVE",
        "evidenceLevel": "A",
        "evidenceDirection": "SUPPORTS",
        "significance": "SENSITIVITYRESPONSE",
        "variantOrigin": "SOMATIC",
        "status": "ACCEPTED",
        "description": "Kopetz et al. colon cancer...",
        "disease": {"name": "Colorectal Cancer"},
        "therapies": [{"name": "Encorafenib"}, {"name": "Cetuximab"}],
        "source": {
            "citation": "Kopetz et al., 2019",
            "citationId": "31566309",
            "sourceType": "PUBMED",
            "link": "https://pubmed.ncbi.nlm.nih.gov/31566309/",
        },
    },
]


def test_parse_civic_evidence_basic():
    items = _parse_civic_evidence(_CIVIC_NODES)
    assert len(items) == 2
    first = items[0]
    assert first["id"] == 1409
    assert first["evidence_type"] == "PREDICTIVE"
    assert first["evidence_level"] == "A"
    assert first["significance"] == "SENSITIVITYRESPONSE"
    assert first["disease"] == "Skin Melanoma"
    assert first["therapies"] == ["Vemurafenib"]
    assert first["citation"] == "Chapman et al., 2011"
    assert first["pmid"] == "21639808"


def test_parse_civic_evidence_multi_therapy():
    items = _parse_civic_evidence(_CIVIC_NODES)
    second = items[1]
    assert second["therapies"] == ["Encorafenib", "Cetuximab"]


def test_parse_civic_evidence_empty():
    assert _parse_civic_evidence([]) == []


def test_parse_civic_evidence_missing_fields():
    """Nodes with missing optional fields should not crash."""
    nodes = [{"id": 999, "evidenceType": "DIAGNOSTIC"}]
    items = _parse_civic_evidence(nodes)
    assert len(items) == 1
    assert items[0]["disease"] is None
    assert items[0]["therapies"] == []


# ── _parse_oncokb_response ─────────────────────────────────────────────────────

_ONCOKB_RAW = {
    "oncogenic": "Oncogenic",
    "geneExist": True,
    "variantExist": True,
    "mutationEffect": {
        "knownEffect": "Gain-of-function",
        "description": "BRAF V600E activates the MAPK pathway...",
        "citations": {"pmids": ["15035987"], "abstracts": []},
    },
    "highestSensitiveLevel": "LEVEL_1",
    "highestResistanceLevel": None,
    "highestDiagnosticImplicationLevel": "LEVEL_Dx2",
    "highestPrognosticImplicationLevel": None,
    "geneSummary": "BRAF is frequently mutated in melanoma.",
    "variantSummary": "The BRAF V600E mutation is known to be oncogenic.",
    "treatments": [
        {
            "level": "LEVEL_1",
            "drugs": [{"drugName": "Vemurafenib"}, {"drugName": "Cobimetinib"}],
            "levelAssociatedCancerType": {"name": "Melanoma"},
            "pmids": ["26091043"],
        }
    ],
    "diagnosticImplications": [
        {
            "levelOfEvidence": "LEVEL_Dx2",
            "tumorType": {"name": "Hairy Cell Leukemia"},
            "alterations": ["V600E"],
        }
    ],
}


def test_parse_oncokb_response_fields():
    result = _parse_oncokb_response(_ONCOKB_RAW)
    assert result["oncogenic"] == "Oncogenic"
    assert result["mutation_effect"] == "Gain-of-function"
    assert result["highest_sensitive_level"] == "LEVEL_1"
    assert result["highest_resistance_level"] is None
    assert result["gene_summary"].startswith("BRAF")
    assert result["variant_summary"].startswith("The BRAF")


def test_parse_oncokb_response_treatments():
    result = _parse_oncokb_response(_ONCOKB_RAW)
    treatments = result["treatments"]
    assert len(treatments) == 1
    t = treatments[0]
    assert t["level"] == "LEVEL_1"
    assert "Vemurafenib" in t["drugs"]
    assert t["cancer_type"] == "Melanoma"
    assert "26091043" in t["pmids"]


def test_parse_oncokb_response_diagnostic():
    result = _parse_oncokb_response(_ONCOKB_RAW)
    diag = result["diagnostic_implications"]
    assert len(diag) == 1
    assert diag[0]["level"] == "LEVEL_Dx2"
    assert diag[0]["cancer_type"] == "Hairy Cell Leukemia"


def test_parse_oncokb_response_empty():
    result = _parse_oncokb_response({})
    assert result["oncogenic"] is None
    assert result["treatments"] == []
    assert result["diagnostic_implications"] == []


# ── Somatic.civic_evidence — monkeypatched _post ──────────────────────────────

_CIVIC_GRAPHQL_RESPONSE = {
    "data": {
        "evidenceItems": {
            "totalCount": 109,
            "nodes": _CIVIC_NODES,
        }
    }
}


@pytest.fixture
def somatic(tmp_path):
    cache = Cache(tmp_path / "cache.db")
    return Somatic(cache)


def test_civic_evidence_parses_response(somatic, monkeypatch):
    monkeypatch.setattr(somatic, "_post", lambda path, json_data=None, data=None: _CIVIC_GRAPHQL_RESPONSE)
    items = somatic.civic_evidence("BRAF", "V600E")
    assert len(items) == 2
    assert items[0]["evidence_type"] == "PREDICTIVE"
    assert items[0]["disease"] == "Skin Melanoma"
    assert items[1]["therapies"] == ["Encorafenib", "Cetuximab"]


def test_civic_evidence_caches(somatic, monkeypatch):
    call_count = {"n": 0}

    def fake_post(path, json_data=None, data=None):
        call_count["n"] += 1
        return _CIVIC_GRAPHQL_RESPONSE

    monkeypatch.setattr(somatic, "_post", fake_post)
    somatic.civic_evidence("BRAF", "V600E")
    somatic.civic_evidence("BRAF", "V600E")   # second call — should hit cache
    assert call_count["n"] == 1


def test_civic_evidence_graphql_error_returns_empty(somatic, monkeypatch):
    """GraphQL errors are logged; partial/empty nodes still return a list."""
    monkeypatch.setattr(somatic, "_post", lambda path, json_data=None, data=None: {
        "errors": [{"message": "some error"}],
        "data": {"evidenceItems": {"totalCount": 0, "nodes": []}},
    })
    items = somatic.civic_evidence("UNKNOWN", "X999X")
    assert items == []


# ── Somatic.oncokb_annotation — no token ──────────────────────────────────────

def test_oncokb_no_token_returns_error_dict_not_exception(somatic, monkeypatch):
    """Without ONCOKB_TOKEN the oncokb section returns an error dict, never raises."""
    monkeypatch.delenv("ONCOKB_TOKEN", raising=False)
    result = somatic.oncokb_annotation("BRAF", "V600E")
    assert isinstance(result, dict)
    assert result.get("type") == "ApiKeyRequiredError"
    assert result.get("env_var") == "ONCOKB_TOKEN"
    assert "error" in result
    assert "oncokb" in result.get("error", "").lower() or result.get("env_var") == "ONCOKB_TOKEN"


# ── Somatic.oncokb_annotation — with token ────────────────────────────────────

def test_oncokb_with_token_parses_response(somatic, monkeypatch):
    """With a token set, _get result is parsed into annotation dict."""
    monkeypatch.setenv("ONCOKB_TOKEN", "test-token-abc")

    # Monkeypatch the internal _oncokb_client.get directly
    class FakeResponse:
        def raise_for_status(self): pass
        def json(self): return _ONCOKB_RAW

    monkeypatch.setattr(somatic._oncokb_client, "get", lambda path, params=None, headers=None: FakeResponse())
    monkeypatch.setattr(somatic, "_rate_wait", lambda: None)

    result = somatic.oncokb_annotation("BRAF", "V600E")
    assert result["oncogenic"] == "Oncogenic"
    assert result["mutation_effect"] == "Gain-of-function"
    assert result["highest_sensitive_level"] == "LEVEL_1"
    assert len(result["treatments"]) == 1


def test_oncokb_with_token_caches(somatic, monkeypatch):
    monkeypatch.setenv("ONCOKB_TOKEN", "test-token-abc")
    call_count = {"n": 0}

    class CountingResponse:
        def raise_for_status(self): pass
        def json(self):
            call_count["n"] += 1
            return _ONCOKB_RAW

    monkeypatch.setattr(somatic._oncokb_client, "get",
                        lambda path, params=None, headers=None: CountingResponse())
    monkeypatch.setattr(somatic, "_rate_wait", lambda: None)

    somatic.oncokb_annotation("BRAF", "V600E")
    somatic.oncokb_annotation("BRAF", "V600E")
    assert call_count["n"] == 1


# ── Somatic.annotate — combined method ────────────────────────────────────────

def test_annotate_combined_string_no_token(somatic, monkeypatch):
    """annotate('BRAF V600E') with no token: civic runs, oncokb returns error dict."""
    monkeypatch.delenv("ONCOKB_TOKEN", raising=False)
    monkeypatch.setattr(somatic, "_post", lambda path, json_data=None, data=None: _CIVIC_GRAPHQL_RESPONSE)

    result = somatic.annotate("BRAF V600E")
    assert result["gene"] == "BRAF"
    assert result["protein_change"] == "V600E"
    assert isinstance(result["civic"], list)
    assert len(result["civic"]) == 2
    assert result["civic_total"] == 2
    # oncokb should contain error dict, not raise
    assert result["oncokb"]["type"] == "ApiKeyRequiredError"


def test_annotate_two_args_with_token(somatic, monkeypatch):
    """annotate('BRAF', 'V600E') with token: both sections succeed."""
    monkeypatch.setenv("ONCOKB_TOKEN", "test-token")
    monkeypatch.setattr(somatic, "_post", lambda path, json_data=None, data=None: _CIVIC_GRAPHQL_RESPONSE)

    class FakeResponse:
        def raise_for_status(self): pass
        def json(self): return _ONCOKB_RAW

    monkeypatch.setattr(somatic._oncokb_client, "get", lambda path, params=None, headers=None: FakeResponse())
    monkeypatch.setattr(somatic, "_rate_wait", lambda: None)

    result = somatic.annotate("BRAF", "V600E")
    assert result["gene"] == "BRAF"
    assert result["protein_change"] == "V600E"
    assert len(result["civic"]) == 2
    assert result["oncokb"]["oncogenic"] == "Oncogenic"


def test_annotate_civic_fails_oncokb_still_errors_gracefully(somatic, monkeypatch):
    """Even if CIViC returns empty, oncokb error dict is still returned cleanly."""
    monkeypatch.delenv("ONCOKB_TOKEN", raising=False)
    monkeypatch.setattr(somatic, "_post", lambda path, json_data=None, data=None: {
        "data": {"evidenceItems": {"totalCount": 0, "nodes": []}}
    })
    result = somatic.annotate("UNKNOWN V999X")
    assert result["civic"] == []
    assert result["oncokb"]["type"] == "ApiKeyRequiredError"


# ── validate ──────────────────────────────────────────────────────────────────

def test_validate_returns_true_on_success(somatic, monkeypatch):
    monkeypatch.setattr(somatic, "_post", lambda path, json_data=None, data=None: {"data": {"__typename": "Query"}})
    assert somatic.validate() is True


def test_validate_returns_false_on_exception(somatic, monkeypatch):
    def fail(*a, **kw): raise ConnectionError("down")
    monkeypatch.setattr(somatic, "_post", fail)
    assert somatic.validate() is False
