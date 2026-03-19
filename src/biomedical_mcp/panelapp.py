"""PanelApp (Genomics England) REST API client — curated gene panels."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class PanelApp(BaseClient):
    domain_name = "panelapp"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def gene_panels(self, gene_symbol: str, confidence: str = "3") -> dict:
        """Which panels include a gene. confidence: '3'=Green, '2'=Amber, '1'=Red, 'all'=any."""
        key = self._cache_key("gene", gene_symbol, confidence)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        data, _ = self._cached_get("gene_raw", f"/genes/{gene_symbol}/")

        results = data.get("results", []) if isinstance(data, dict) else []
        panels = []
        for entry in results:
            panel = entry.get("panel") or {}
            conf_level = str(entry.get("confidence_level", ""))
            if confidence != "all" and conf_level < confidence:
                continue
            panels.append({
                "panel_id": panel.get("id"),
                "panel_name": panel.get("name"),
                "panel_version": panel.get("version"),
                "confidence_level": conf_level,
                "confidence_label": _confidence_label(conf_level),
                "mode_of_inheritance": entry.get("mode_of_inheritance"),
                "phenotypes": entry.get("phenotypes", []),
            })

        result = {
            "gene_symbol": gene_symbol,
            "confidence_filter": confidence,
            "count": len(panels),
            "panels": panels,
        }
        self.cache.set(key, result)
        return result

    def panel_genes(self, panel_id: int | str, signed_off: bool = True,
                    confidence: str = "3") -> dict:
        """All genes in a panel. Default: signed-off version, Green confidence only."""
        endpoint = f"/panels/{panel_id}/"
        if signed_off:
            endpoint = f"/panels/signedoff/{panel_id}/"

        key = self._cache_key("panel", str(panel_id), str(signed_off), confidence)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        try:
            data, _ = self._cached_get("panel_raw", endpoint)
        except Exception:
            # Fallback to non-signedoff if signedoff not available
            if signed_off:
                data, _ = self._cached_get("panel_raw", f"/panels/{panel_id}/")
            else:
                raise

        genes_raw = data.get("genes", []) if isinstance(data, dict) else []
        genes = []
        for gene in genes_raw:
            conf_level = str(gene.get("confidence_level", ""))
            if confidence != "all" and conf_level < confidence:
                continue
            genes.append({
                "gene_symbol": gene.get("gene_data", {}).get("gene_symbol"),
                "confidence_level": conf_level,
                "confidence_label": _confidence_label(conf_level),
                "mode_of_inheritance": gene.get("mode_of_inheritance"),
                "phenotypes": gene.get("phenotypes", []),
                "evidence": gene.get("evidence", []),
            })

        result = {
            "panel_id": data.get("id"),
            "panel_name": data.get("name"),
            "panel_version": data.get("version"),
            "signed_off": signed_off,
            "confidence_filter": confidence,
            "gene_count": len(genes),
            "genes": genes,
        }
        self.cache.set(key, result)
        return result

    def panel_search(self, query: str, limit: int = 20) -> dict:
        """Search panels by disease/keyword."""
        params = {"search": query, "page_size": min(limit, 100)}
        data, cache_hit = self._cached_get("search", "/panels/", params)
        results = data.get("results", []) if isinstance(data, dict) else []
        return {
            "query": query,
            "count": len(results),
            "panels": [
                {
                    "panel_id": p.get("id"),
                    "name": p.get("name"),
                    "version": p.get("version"),
                    "gene_count": p.get("stats", {}).get("number_of_genes", 0),
                    "relevant_disorders": p.get("relevant_disorders", []),
                }
                for p in results[:limit]
            ],
        }

    def validate(self) -> bool:
        try:
            self._get("/panels/", params={"page_size": 1})
            return True
        except Exception:
            return False


def _confidence_label(level: str) -> str:
    return {"3": "Green", "2": "Amber", "1": "Red", "0": "None"}.get(level, f"Level {level}")
