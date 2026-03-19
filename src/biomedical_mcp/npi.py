"""NPI/NPPES Registry API client."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

NPI_BASE = "https://npiregistry.cms.hhs.gov/api/"


def _is_retryable(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 500, 502, 503, 504)
    return isinstance(exc, (httpx.ConnectError, httpx.ReadTimeout))


class NPI:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(
            timeout=15,
            headers={"User-Agent": "biomedical-mcp/0.1"},
        )

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10), retry=retry_if_exception(_is_retryable))
    def _get(self, params: dict) -> dict:
        resp = self.client.get(NPI_BASE, params={"version": "2.1", **params})
        resp.raise_for_status()
        return resp.json()

    def _cached(self, prefix: str, params: dict, max_age_days: int = 7) -> dict:
        raw = f"{prefix}:{sorted(params.items())}"
        key = f"npi:{hashlib.md5(raw.encode()).hexdigest()}"
        cached = self.cache.get(key, max_age_days=max_age_days)
        if cached is not None:
            return cached
        result = self._get(params)
        self.cache.set(key, result)
        return result

    def _extract_provider(self, result: dict) -> dict:
        basic = result.get("basic", {})
        taxonomies = result.get("taxonomies", [])
        addresses = result.get("addresses", [])
        practice_addr = next((a for a in addresses if a.get("address_purpose") == "LOCATION"), addresses[0] if addresses else {})
        return {
            "npi": result.get("number"),
            "type": "individual" if result.get("enumeration_type") == "NPI-1" else "organization",
            "name": f"{basic.get('first_name', '')} {basic.get('last_name', '')}".strip() or basic.get("organization_name"),
            "credential": basic.get("credential"),
            "gender": basic.get("gender"),
            "status": basic.get("status"),
            "specialties": [
                {"taxonomy": t.get("code"), "description": t.get("desc"), "primary": t.get("primary")}
                for t in taxonomies
            ],
            "address": {
                "line1": practice_addr.get("address_1"),
                "city": practice_addr.get("city"),
                "state": practice_addr.get("state"),
                "postal_code": practice_addr.get("postal_code"),
            } if practice_addr else None,
        }

    def search(
        self,
        name: str | None = None,
        specialty: str | None = None,
        city: str | None = None,
        state: str | None = None,
        limit: int = 10,
    ) -> list[dict]:
        params = {"limit": min(limit, 200)}
        if name:
            parts = name.split(maxsplit=1)
            if len(parts) == 2:
                params["first_name"] = parts[0]
                params["last_name"] = parts[1]
            else:
                params["last_name"] = parts[0]
        if specialty:
            params["taxonomy_description"] = specialty
        if city:
            params["city"] = city
        if state:
            params["state"] = state

        data = self._cached("search", params)
        results = data.get("results", [])
        return [self._extract_provider(r) for r in results]

    def lookup(self, npi: str) -> dict | None:
        data = self._cached("lookup", {"number": npi})
        results = data.get("results", [])
        if not results:
            return None
        return self._extract_provider(results[0])
