"""DDInter 2.0 client — drug-drug interaction database.

API notes (verified 2026-06-07):
  DDInter has NO documented public REST API. It exposes Django DataTables
  AJAX endpoints that power its web UI. These are stable, return JSON, and
  work without authentication or cookies.

  Endpoint 1 — drug search (POST):
    POST /server/search-source/
    Body: keyword=<name>&draw=1&start=0&length=<n>
    Returns: {"data": [{"internal_id": "DDInter20", "drug": "Acetylsalicylic acid", ...}]}
    The `internal_id` field is the opaque string ID used in all subsequent calls.

  Endpoint 2 — interactions for a drug (POST):
    POST /server/interact-with/<internal_id>/
    Body: draw=1&start=0&length=<n>&search[value]=&search[regex]=false
    Returns: {"recordsTotal": N, "data": [{"drug_name": ..., "level": 1|2|3, ...}]}

  Severity levels (from explanation page):
    1 → Minor    (limited clinical effect, usually no therapy change)
    2 → Moderate (may exacerbate disease or require therapy change)
    3 → Major    (life-threatening, requires medical intervention)
    0 → Unknown  (no severity annotation; from Sci Transl Med corpus)

  Mechanism flags are boolean strings "0"/"1" on these keys:
    absorption, distribution, metabolism, excretion,
    synergistic_effect, antagonistic_effect, others
"""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

BASE_URL = "https://ddinter2.scbdd.com"

_LEVEL_MAP = {1: "Minor", 2: "Moderate", 3: "Major", 0: "Unknown"}

_MECHANISM_KEYS = (
    "absorption",
    "distribution",
    "metabolism",
    "excretion",
    "synergistic_effect",
    "antagonistic_effect",
    "others",
)


def _mechanism_label(row: dict) -> str | None:
    """Return a comma-joined string of active mechanism flags, or None."""
    active = [k for k in _MECHANISM_KEYS if row.get(k) == "1"]
    if not active:
        return None
    return ", ".join(k.replace("_", " ") for k in active)


class DDInter(BaseClient):
    """Client for DDInter 2.0 drug-drug interaction database.

    Uses undocumented DataTables AJAX endpoints (POST, JSON responses).
    No API key required. Rate-limit: 1 req/s (conservative; server has no
    documented limit but is a small academic resource).
    """

    domain_name = "ddinter"

    def __init__(self, cache: Cache):
        super().__init__(cache, base_url=BASE_URL)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _post_datatable(
        self,
        path: str,
        *,
        keyword: str | None = None,
        length: int = 25,
    ) -> dict:
        """POST a DataTables AJAX request and return the parsed JSON."""
        body: dict[str, str | int] = {
            "draw": 1,
            "start": 0,
            "length": length,
            "search[value]": keyword or "",
            "search[regex]": "false",
        }
        if keyword:
            body["keyword"] = keyword
        self._rate_wait()
        resp = self.client.post(
            path,
            data=body,
            headers={"X-Requested-With": "XMLHttpRequest"},
        )
        resp.raise_for_status()
        return resp.json()

    def _resolve_drug_id(self, drug_name: str) -> str | None:
        """Resolve a drug name to a DDInter internal_id string (e.g. 'DDInter20').

        Tries an exact case-insensitive match first, then falls back to the
        first result from the server's own search ranking.
        """
        cache_key = self._cache_key("resolve", drug_name.lower())
        cached = self.cache.get(cache_key, max_age_days=self.ttl_days)
        if cached is not None:
            return cached.get("internal_id")

        data = self._post_datatable("/server/search-source/", keyword=drug_name, length=10)
        rows = data.get("data", [])
        if not rows:
            return None

        name_lower = drug_name.lower()
        # Prefer exact name match
        for row in rows:
            if row.get("drug", "").lower() == name_lower:
                iid = row["internal_id"]
                self.cache.set(cache_key, {"internal_id": iid})
                return iid

        # Fall back to first result
        iid = rows[0]["internal_id"]
        self.cache.set(cache_key, {"internal_id": iid})
        return iid

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def interactions(self, drug_name: str, limit: int = 25) -> list[dict]:
        """Return drug-drug interactions for *drug_name* from DDInter 2.0.

        Each result dict has:
          interacting_drug  str   — name of the interacting drug
          ddinter_id        str   — DDInter internal ID of interacting drug
          severity          str   — "Major" | "Moderate" | "Minor" | "Unknown"
          mechanism         str | None  — comma-joined mechanism categories
          interaction_id    int   — DDInter interaction record ID

        Returns an empty list if the drug is not found in DDInter.
        """
        internal_id = self._resolve_drug_id(drug_name)
        if not internal_id:
            log.debug("DDInter: no match for %r", drug_name)
            return []

        cache_key = self._cache_key("interactions", internal_id, str(limit))
        cached = self.cache.get(cache_key, max_age_days=self.ttl_days)
        if cached is not None:
            return cached

        path = f"/server/interact-with/{internal_id}/"
        data = self._post_datatable(path, length=min(limit, 500))

        results = []
        for row in data.get("data", [])[:limit]:
            level_int = int(row.get("level", 0))
            results.append(
                {
                    "interacting_drug": row.get("drug_name", ""),
                    "ddinter_id": row.get("drug_id", ""),
                    "severity": _LEVEL_MAP.get(level_int, "Unknown"),
                    "mechanism": _mechanism_label(row),
                    "interaction_id": row.get("interaction_id"),
                }
            )

        self.cache.set(cache_key, results)
        return results

    def validate(self) -> bool:
        """Health check — verify DDInter search endpoint is reachable."""
        try:
            data = self._post_datatable("/server/search-source/", keyword="aspirin", length=1)
            return bool(data.get("data"))
        except Exception as exc:
            log.warning("DDInter health check failed: %s", exc)
            return False
