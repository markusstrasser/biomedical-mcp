"""Monarch Initiative v3 API client — cross-species phenotype-gene-disease associations."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class Monarch(BaseClient):
    domain_name = "monarch"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def disease_phenotypes(self, disease_id: str, limit: int = 50) -> dict:
        """Get phenotype profile for a disease (OMIM, MONDO, Orphanet ID)."""
        key = self._cache_key("disease_pheno", disease_id)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        data = self._get(f"/association", params={
            "subject": disease_id,
            "predicate": "biolink:has_phenotype",
            "limit": min(limit, 200),
        })
        items = data.get("items", []) if isinstance(data, dict) else []
        result = {
            "disease_id": disease_id,
            "phenotype_count": len(items),
            "phenotypes": [
                {
                    "hpo_id": item.get("object"),
                    "name": item.get("object_label"),
                    "frequency": item.get("frequency_qualifier_label"),
                }
                for item in items
            ],
        }
        self.cache.set(key, result)
        return result

    def search(self, query: str, category: str | None = None, limit: int = 10) -> dict:
        """Search Monarch for genes, diseases, or phenotypes."""
        params = {"q": query, "limit": min(limit, 50)}
        if category:
            params["category"] = category
        data, _ = self._cached_get("search", "/search", params)
        items = data.get("items", []) if isinstance(data, dict) else []
        return {
            "query": query,
            "count": len(items),
            "results": [
                {
                    "id": item.get("id"),
                    "name": item.get("name"),
                    "category": item.get("category"),
                    "description": (item.get("description") or "")[:300],
                }
                for item in items[:limit]
            ],
        }

    def validate(self) -> bool:
        try:
            self._get("/search", params={"q": "BRCA1", "limit": 1})
            return True
        except Exception:
            return False
