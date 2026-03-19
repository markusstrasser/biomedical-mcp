"""InterPro (EBI) REST API client — protein domains, families, sites."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class InterPro(BaseClient):
    domain_name = "interpro"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def protein_domains(self, uniprot_id: str | None = None,
                        gene_symbol: str | None = None) -> dict:
        """Get domain/family/site annotations for a protein."""
        query = uniprot_id or gene_symbol
        if not query:
            return {"error": "Provide uniprot_id or gene_symbol"}

        key = self._cache_key("domains", query)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        # InterPro API: /protein/uniprot/{accession} returns entries matching
        path = f"/protein/uniprot/{query}"
        try:
            data = self._get(path, params={"format": "json"})
        except Exception:
            # Try as gene symbol search
            if gene_symbol and not uniprot_id:
                return {"gene_symbol": gene_symbol, "error": "not_found",
                        "note": "InterPro requires UniProt accession. Resolve via genetics_gene_info first."}
            return {"query": query, "error": "not_found"}

        # Extract domain entries
        entries = []
        results = data.get("results", []) if isinstance(data, dict) else []
        for result in results:
            metadata = result.get("metadata", {})
            entries.append({
                "accession": metadata.get("accession"),
                "name": metadata.get("name"),
                "type": metadata.get("type"),  # domain, family, homologous_superfamily, etc.
                "source_database": metadata.get("source_database"),
                "go_terms": [g.get("identifier") for g in metadata.get("go_terms", [])],
            })

        result = {
            "query": query,
            "entry_count": len(entries),
            "entries": entries,
        }
        self.cache.set(key, result)
        return result

    def validate(self) -> bool:
        try:
            self._get("/entry/interpro", params={"page_size": 1, "format": "json"})
            return True
        except Exception:
            return False
