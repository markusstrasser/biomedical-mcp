"""Offline tests for the PharmVar client.

Monkeypatches _get with canned /genes/{symbol} records (shape from live probing 2026-10-08).
No network calls are made.
"""

from __future__ import annotations

import httpx
import pytest

from biomedical_mcp.cache import Cache
from biomedical_mcp.pharmvar import PharmVar


def _allele(name: str, kind: str, function: str = "normal function") -> dict:
    return {"geneSymbol": "CYP2C19", "alleleName": name, "pvId": f"PV{len(name)}", "alleleType": kind,
            "function": function, "activityScore": None, "evidenceLevel": "Definitive"}


GENE_CYP2C19 = {
    "geneSymbol": "CYP2C19",
    "alleles": [
        _allele("CYP2C19*2.018", "Sub", "no function"),
        _allele("CYP2C19*10", "Core", "decreased function"),
        _allele("CYP2C19*2", "Core", "no function"),
        _allele("CYP2C19*1", "Core"),
        _allele("CYP2C19*2.001", "Sub", "no function"),
    ],
}


def _status_error(code: int) -> httpx.HTTPStatusError:
    req = httpx.Request("GET", "https://www.pharmvar.org/api-service/genes/X")
    return httpx.HTTPStatusError("err", request=req, response=httpx.Response(code, request=req))


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("PHARMVAR_API_KEY", "test-key")
    return PharmVar(Cache(tmp_path / "cache.db"))


def test_key_sent_as_api_key_header(client):
    assert client.client.headers["API-Key"] == "test-key"
    assert "Authorization" not in client.client.headers


def test_uses_gene_endpoint_and_orders_core_first(client, monkeypatch):
    calls = []
    monkeypatch.setattr(client, "_get", lambda path, params=None, headers=None: calls.append(path) or GENE_CYP2C19)
    r = client.star_alleles("cyp2c19", limit=4)
    assert calls == ["/genes/CYP2C19"]
    assert [a["name"] for a in r["alleles"]] == ["CYP2C19*1", "CYP2C19*2", "CYP2C19*10", "CYP2C19*2.001"]
    assert (r["allele_count"], r["core_count"]) == (5, 3)


def test_second_call_is_cached_and_honours_limit(client, monkeypatch):
    monkeypatch.setattr(client, "_get", lambda *a, **kw: GENE_CYP2C19)
    client.star_alleles("CYP2C19")
    monkeypatch.setattr(client, "_get", lambda *a, **kw: pytest.fail("should hit cache"))
    r = client.star_alleles("CYP2C19", limit=2)
    assert r["_cache_hit"] is True and len(r["alleles"]) == 2 and r["allele_count"] == 5


def test_rejected_key_returns_actionable_error(client, monkeypatch):
    def raise_401(*a, **kw):
        raise _status_error(401)
    monkeypatch.setattr(client, "_get", raise_401)
    r = client.star_alleles("CYP2C19")
    assert r["type"] == "ApiKeyRequiredError" and r["env_var"] == "PHARMVAR_API_KEY"


def test_unknown_gene_lists_supported_genes(client, monkeypatch):
    def raise_404(*a, **kw):
        raise _status_error(404)
    monkeypatch.setattr(client, "_get", raise_404)
    monkeypatch.setattr(client, "_cached_get", lambda *a, **kw: (["CYP2D6", "CYP2C19"], False))
    r = client.star_alleles("TPMT")
    assert r["type"] == "NotFoundError" and "CYP2C19, CYP2D6" in r["suggestion"]


def test_gene_without_alleles_is_not_found(client, monkeypatch):
    monkeypatch.setattr(client, "_get", lambda *a, **kw: {"geneSymbol": "TPMT", "alleles": []})
    monkeypatch.setattr(client, "_cached_get", lambda *a, **kw: (["CYP2D6"], False))
    assert client.star_alleles("TPMT")["type"] == "NotFoundError"
