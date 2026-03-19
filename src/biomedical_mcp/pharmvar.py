"""PharmVar REST API client — star allele definitions for pharmacogenomics genes.

Requires PHARMVAR_API_KEY environment variable (free, register at pharmvar.org).
"""

from __future__ import annotations

import logging
import os

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

PHARMVAR_BASE = "https://www.pharmvar.org/api-service"


class PharmVar(BaseClient):
    domain_name = "pharmvar"

    def __init__(self, cache: Cache):
        super().__init__(cache, base_url=PHARMVAR_BASE)
        api_key = os.environ.get("PHARMVAR_API_KEY", "")
        if api_key:
            self.client.headers["Authorization"] = api_key
        else:
            log.warning("PHARMVAR_API_KEY not set — PharmVar tools will return auth errors")

    def star_alleles(self, gene: str, limit: int = 50) -> dict:
        """Get star allele definitions for a PGx gene (e.g., CYP2D6, CYP2C19)."""
        key = self._cache_key("alleles", gene)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        try:
            data = self._get("/alleles", params={"gene": gene})
        except Exception as exc:
            if "401" in str(exc):
                return {"gene": gene, "error": "auth_required",
                        "note": "Set PHARMVAR_API_KEY env var (free at pharmvar.org)"}
            raise

        alleles = data if isinstance(data, list) else data.get("alleles", []) if isinstance(data, dict) else []
        result = {
            "gene": gene,
            "allele_count": len(alleles),
            "alleles": [
                {
                    "name": a.get("alleleName") or a.get("name"),
                    "function": a.get("function"),
                    "activity_score": a.get("activityScore"),
                }
                for a in alleles[:limit]
            ],
        }
        self.cache.set(key, result)
        return result

    def validate(self) -> bool:
        if not os.environ.get("PHARMVAR_API_KEY"):
            return False
        try:
            self._get("/genes")
            return True
        except Exception:
            return False
