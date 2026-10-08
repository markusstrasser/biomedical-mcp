"""PharmVar REST API client — star allele definitions for pharmacogenomics genes.

Requires PHARMVAR_API_KEY environment variable (free, register at pharmvar.org).

API shape (probed 2026-10-08):
- The key goes in the `API-Key` header; `Authorization` returns 401.
- `/genes/{symbol}` returns one gene with all its alleles (~1 MB, ~20 s server time).
- `/alleles` ignores `gene`/`geneSymbol` and returns every gene's alleles (~25 MB, ~50 s).
- `/genes/list` returns the supported gene symbols.
"""

from __future__ import annotations

import logging
import os
import re

import httpx

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache
from biomedical_mcp.errors import ApiKeyRequiredError, NotFoundError

log = logging.getLogger(__name__)

PHARMVAR_BASE = "https://www.pharmvar.org/api-service"
DOCS_URL = "https://www.pharmvar.org"
# Gene records take ~20 s to serve; the base client's 30 s default leaves no margin.
TIMEOUT = httpx.Timeout(120.0, connect=15.0)


def _star_order(name: str) -> tuple:
    """Sort key: CYP2C19*2 < CYP2C19*2.001 < CYP2C19*10 (numeric, not lexical)."""
    star = name.split("*", 1)[1] if "*" in name else name
    return tuple(int(n) for n in re.findall(r"\d+", star)), star


class PharmVar(BaseClient):
    domain_name = "pharmvar"

    def __init__(self, cache: Cache):
        super().__init__(cache, base_url=PHARMVAR_BASE)
        self.client.timeout = TIMEOUT
        api_key = os.environ.get("PHARMVAR_API_KEY", "")
        if api_key:
            self.client.headers["API-Key"] = api_key
        else:
            log.warning("PHARMVAR_API_KEY not set — PharmVar tools will return auth errors")

    def star_alleles(self, gene: str, limit: int = 50) -> dict:
        """Get star allele definitions for a PGx gene (e.g., CYP2D6, CYP2C19).

        Core alleles come first in star order, then suballeles.
        """
        gene = gene.strip().upper()
        key = self._cache_key("gene", gene)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is None:
            try:
                data = self._get(f"/genes/{gene}")
            except httpx.HTTPStatusError as exc:
                status = exc.response.status_code
                if status == 401:
                    err = ApiKeyRequiredError("PharmVar", "PHARMVAR_API_KEY", DOCS_URL).to_dict()
                    err["note"] = "Key missing or rejected."
                    return {"gene": gene, **err}
                if status == 404:
                    return self._not_curated(gene)
                raise
            # Genes PharmVar does not curate (e.g. TPMT) come back 200 with no alleles.
            if not data.get("alleles"):
                return self._not_curated(gene)
            alleles = sorted(
                (
                    {
                        "name": a.get("alleleName"),
                        "core": a.get("alleleType") == "Core",
                        "function": a.get("function"),
                        "activity_score": a.get("activityScore"),
                        "evidence_level": a.get("evidenceLevel"),
                        "pv_id": a.get("pvId"),
                    }
                    for a in data.get("alleles", [])
                ),
                key=lambda a: (not a["core"], _star_order(a["name"] or "")),
            )
            cached = {"gene": data.get("geneSymbol", gene), "alleles": alleles}
            self.cache.set(key, cached)
            hit = False
        else:
            hit = True
        alleles = cached["alleles"]
        return {
            "gene": cached["gene"],
            "allele_count": len(alleles),
            "core_count": sum(a["core"] for a in alleles),
            "alleles": alleles[:limit],
            **({"_cache_hit": True} if hit else {}),
        }

    def _not_curated(self, gene: str) -> dict:
        return {"gene": gene, **NotFoundError(
            "PharmVar allele set", gene, f"PharmVar curates: {', '.join(self.genes())}").to_dict()}

    def genes(self) -> list[str]:
        """Gene symbols PharmVar curates."""
        data, _ = self._cached_get("genes", "/genes/list")
        return sorted(data) if isinstance(data, list) else []

    def validate(self) -> bool:
        if not os.environ.get("PHARMVAR_API_KEY"):
            return False
        try:
            return bool(self.genes())
        except Exception:
            return False
