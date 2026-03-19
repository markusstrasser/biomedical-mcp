"""HGNC (HUGO Gene Nomenclature Committee) REST API client — official gene names."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class HGNC(BaseClient):
    domain_name = "hgnc"

    def __init__(self, cache: Cache):
        super().__init__(cache)
        # HGNC requires Accept: application/json header
        self.client.headers["Accept"] = "application/json"

    def gene_names(self, symbol: str) -> dict:
        """Official HGNC nomenclature: symbol, name, aliases, previous symbols, locus group."""
        key = self._cache_key("symbol", symbol)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        data = self._get(f"/fetch/symbol/{symbol}")
        docs = data.get("response", {}).get("docs", [])
        if not docs:
            return {"symbol": symbol, "error": "not_found"}

        doc = docs[0]
        result = {
            "hgnc_id": doc.get("hgnc_id"),
            "symbol": doc.get("symbol"),
            "name": doc.get("name"),
            "locus_group": doc.get("locus_group"),
            "locus_type": doc.get("locus_type"),
            "aliases": doc.get("alias_symbol", []),
            "previous_symbols": doc.get("prev_symbol", []),
            "ensembl_id": doc.get("ensembl_gene_id"),
            "entrez_id": doc.get("entrez_id"),
            "uniprot_ids": doc.get("uniprot_ids", []),
            "gene_family": doc.get("gene_group", []),
            "chromosome": doc.get("location"),
        }
        self.cache.set(key, result)
        return result

    def search(self, query: str, limit: int = 10) -> dict:
        """Search HGNC by keyword."""
        data, _ = self._cached_get("search", f"/search/{query}")
        docs = data.get("response", {}).get("docs", [])
        return {
            "query": query,
            "count": len(docs),
            "results": [
                {
                    "hgnc_id": d.get("hgnc_id"),
                    "symbol": d.get("symbol"),
                    "name": d.get("name"),
                    "locus_group": d.get("locus_group"),
                }
                for d in docs[:limit]
            ],
        }

    def validate(self) -> bool:
        try:
            self._get("/fetch/symbol/BRCA1")
            return True
        except Exception:
            return False
