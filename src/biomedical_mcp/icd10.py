"""NLM ICD-10-CM API client."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

ICD10_BASE = "https://clinicaltables.nlm.nih.gov/api/icd10cm/v3"


def _is_retryable(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 500, 502, 503, 504)
    return isinstance(exc, (httpx.ConnectError, httpx.ReadTimeout))


class ICD10:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(
            base_url=ICD10_BASE,
            timeout=15,
            headers={"User-Agent": "biomedical-mcp/0.1"},
        )

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10), retry=retry_if_exception(_is_retryable))
    def _get(self, path: str, params: dict) -> list:
        resp = self.client.get(path, params=params)
        resp.raise_for_status()
        return resp.json()

    def _cached(self, prefix: str, params: dict, max_age_days: int = 90) -> list:
        raw = f"{prefix}:{sorted(params.items())}"
        key = f"icd10:{hashlib.md5(raw.encode()).hexdigest()}"
        cached = self.cache.get(key, max_age_days=max_age_days)
        if cached is not None:
            return cached
        result = self._get("/search", params)
        self.cache.set(key, result)
        return result

    def search(self, query: str, limit: int = 10) -> list[dict]:
        # NLM API returns: [totalCount, [codes], null, [[code, description], ...]]
        data = self._cached("search", {"sf": "code,name", "terms": query, "maxList": min(limit, 50)})
        if not data or len(data) < 4:
            return []
        pairs = data[3] if data[3] else []
        return [{"code": p[0], "description": p[1]} for p in pairs]

    def lookup(self, code: str) -> dict | None:
        # Search for exact code
        data = self._cached("lookup", {"sf": "code,name", "terms": code, "maxList": 20})
        if not data or len(data) < 4:
            return None
        pairs = data[3] if data[3] else []
        # Find exact match
        for p in pairs:
            if p[0] == code:
                return {"code": p[0], "description": p[1]}
        # Find prefix matches (subcodes)
        matches = [{"code": p[0], "description": p[1]} for p in pairs if p[0].startswith(code)]
        if matches:
            return {"code": code, "description": matches[0]["description"], "subcodes": matches}
        return None
