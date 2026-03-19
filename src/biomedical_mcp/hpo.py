"""HPO (Human Phenotype Ontology) REST API client — phenotype terms and gene associations."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class HPO(BaseClient):
    domain_name = "hpo"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def search(self, query: str, limit: int = 10) -> dict:
        """Search HPO terms by keyword."""
        data, _ = self._cached_get(
            "search", "/hp/search",
            params={"q": query, "max": min(limit, 50)},
            ttl_days=self.ttl_days,
        )
        terms = data.get("terms", []) if isinstance(data, dict) else data if isinstance(data, list) else []
        return {
            "query": query,
            "count": len(terms),
            "terms": [
                {
                    "hpo_id": t.get("id"),
                    "name": t.get("name"),
                    "definition": t.get("definition"),
                }
                for t in terms[:limit]
            ],
        }

    def term(self, hpo_id: str) -> dict:
        """Get HPO term details: definition, synonyms, parents, children."""
        key = self._cache_key("term", hpo_id)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        data = self._get(f"/hp/terms/{hpo_id}")
        result = {
            "hpo_id": data.get("id"),
            "name": data.get("name"),
            "definition": data.get("definition"),
            "synonyms": data.get("synonyms", []),
            "parents": [{"id": p.get("id"), "name": p.get("name")} for p in data.get("parents", [])],
            "children": [{"id": c.get("id"), "name": c.get("name")} for c in data.get("children", [])],
        }
        self.cache.set(key, result)
        return result

    def genes(self, hpo_id: str) -> dict:
        """Get genes associated with an HPO term."""
        key = self._cache_key("genes", hpo_id)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        data = self._get(f"/hp/terms/{hpo_id}/genes")
        genes = data.get("genes", []) if isinstance(data, dict) else data if isinstance(data, list) else []
        result = {
            "hpo_id": hpo_id,
            "gene_count": len(genes),
            "genes": [
                {
                    "gene_symbol": g.get("symbol") or g.get("name"),
                    "entrez_id": g.get("entrezGeneId"),
                }
                for g in genes
            ],
        }
        self.cache.set(key, result)
        return result

    def gene_phenotypes(self, gene_symbol: str) -> dict:
        """Get HPO phenotypes associated with a gene."""
        key = self._cache_key("gene_pheno", gene_symbol)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        data = self._get(f"/hp/genes/{gene_symbol}/terms")
        terms = data.get("terms", []) if isinstance(data, dict) else data if isinstance(data, list) else []
        result = {
            "gene_symbol": gene_symbol,
            "phenotype_count": len(terms),
            "phenotypes": [
                {"hpo_id": t.get("id"), "name": t.get("name")}
                for t in terms
            ],
        }
        self.cache.set(key, result)
        return result

    def validate(self) -> bool:
        try:
            self._get("/hp/search", params={"q": "seizure", "max": 1})
            return True
        except Exception:
            return False
