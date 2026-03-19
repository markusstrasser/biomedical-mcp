"""RCSB PDB REST API client — experimental protein structures."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class PDB(BaseClient):
    domain_name = "pdb"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def structures_for_gene(self, uniprot_id: str) -> dict:
        """Get experimental PDB structures for a UniProt accession."""
        key = self._cache_key("gene_structs", uniprot_id)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        # RCSB search API
        search_payload = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_id,
                },
            },
            "return_type": "entry",
            "request_options": {
                "results_content_type": ["experimental"],
                "paginate": {"start": 0, "rows": 25},
                "sort": [{"sort_by": "rcsb_accession_info.deposit_date", "direction": "desc"}],
            },
        }

        try:
            resp = self.client.post(
                "https://search.rcsb.org/rcsbsearch/v2/query",
                json=search_payload,
            )
            resp.raise_for_status()
            search_data = resp.json()
        except Exception as exc:
            log.warning("PDB search failed for %s: %s", uniprot_id, exc)
            return {"uniprot_id": uniprot_id, "error": "search_failed", "structures": []}

        pdb_ids = [r.get("identifier") for r in search_data.get("result_set", [])]

        structures = []
        for pdb_id in pdb_ids[:10]:  # detail for top 10
            detail = self._structure_detail(pdb_id)
            if detail:
                structures.append(detail)

        result = {
            "uniprot_id": uniprot_id,
            "total_structures": search_data.get("total_count", len(pdb_ids)),
            "structures": structures,
        }
        self.cache.set(key, result)
        return result

    def _structure_detail(self, pdb_id: str) -> dict | None:
        """Get detail for a single PDB entry."""
        key = self._cache_key("detail", pdb_id)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return cached

        try:
            data = self._get(f"/rest/v1/core/entry/{pdb_id}")
        except Exception:
            return None

        entry = data.get("rcsb_entry_info", {}) if isinstance(data, dict) else {}
        citation = data.get("rcsb_primary_citation", {}) if isinstance(data, dict) else {}
        result = {
            "pdb_id": data.get("rcsb_id", pdb_id),
            "title": (data.get("struct", {}) or {}).get("title"),
            "method": entry.get("experimental_method"),
            "resolution": entry.get("resolution_combined", [None])[0] if entry.get("resolution_combined") else None,
            "deposit_date": (data.get("rcsb_accession_info", {}) or {}).get("deposit_date"),
            "citation_doi": citation.get("pdbx_database_id_doi"),
        }
        self.cache.set(key, result)
        return result

    def structure_detail(self, pdb_id: str) -> dict:
        """Get full detail for a PDB entry."""
        result = self._structure_detail(pdb_id)
        if not result:
            return {"pdb_id": pdb_id, "error": "not_found"}
        return result

    def validate(self) -> bool:
        try:
            self._get("/rest/v1/core/entry/1TUP")
            return True
        except Exception:
            return False
