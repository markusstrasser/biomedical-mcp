"""LitVar2 (NCBI) REST API client — variant-level literature mining."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class LitVar(BaseClient):
    domain_name = "litvar"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def variant_publications(self, variant_id: str, limit: int = 20) -> dict:
        """Papers mentioning a variant (rsID or HGVS).

        Returns publication count and PMIDs from LitVar2.
        Falls back to gene-level search note if variant returns 0 hits.
        """
        key = self._cache_key("pubs", variant_id)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        # LitVar2 autocomplete to get variant metadata + publication count
        autocomplete = self._get("/variant/autocomplete/", params={"query": variant_id})
        if not autocomplete:
            return {"variant_id": variant_id, "publication_count": 0, "pmids": [],
                    "note": "No LitVar2 match. Try gene-level search."}

        first = autocomplete[0] if isinstance(autocomplete, list) else autocomplete
        litvar_id = first.get("_id") or variant_id
        pub_count = first.get("pmids_count", 0)
        genes = first.get("gene", [])
        hgvs = first.get("hgvs") or first.get("name")

        # Fetch actual PMIDs
        pmids = []
        if pub_count > 0:
            try:
                import urllib.parse
                encoded_id = urllib.parse.quote(litvar_id, safe="")
                pub_data = self._get(f"/variant/get/{encoded_id}/publications")
                raw_pmids = pub_data.get("pmids", []) if isinstance(pub_data, dict) else pub_data if isinstance(pub_data, list) else []
                pmids = raw_pmids[:limit]
            except Exception:
                pass

        result = {
            "variant_id": variant_id,
            "litvar_id": litvar_id,
            "hgvs": hgvs,
            "genes": genes,
            "publication_count": pub_count,
            "pmids": pmids,
        }
        self.cache.set(key, result)
        return result

    def validate(self) -> bool:
        try:
            self._get("/variant/autocomplete/", params={"query": "rs1801133"})
            return True
        except Exception:
            return False
