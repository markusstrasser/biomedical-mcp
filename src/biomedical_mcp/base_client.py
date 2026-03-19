"""Base API client with caching, rate limiting, retry, and provenance."""

from __future__ import annotations

import hashlib
import logging
import threading
import time

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache
from biomedical_mcp.config import get_domain
from biomedical_mcp.provenance import make_provenance

log = logging.getLogger(__name__)

# Global per-host rate limiter (thread-safe)
_rate_lock = threading.Lock()
_last_call: dict[str, float] = {}


def _wait_rate_limit(host: str, max_rate: float) -> None:
    """Block until the next request to `host` is allowed."""
    min_interval = 1.0 / max_rate
    with _rate_lock:
        last = _last_call.get(host, 0.0)
        now = time.monotonic()
        wait_time = min_interval - (now - last)
        if wait_time > 0:
            time.sleep(wait_time)
        _last_call[host] = time.monotonic()


def _is_retryable(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 500, 502, 503, 504)
    return isinstance(exc, (httpx.ConnectError, httpx.ReadTimeout))


class BaseClient:
    """Base for all API clients — provides caching, rate limiting, retry."""

    domain_name: str = ""  # override in subclass

    def __init__(self, cache: Cache, base_url: str | None = None):
        self.cache = cache
        cfg = get_domain(self.domain_name)
        self.base_url = base_url or cfg.get("base_url", "")
        self.ttl_days = cfg.get("ttl_days", 7)
        self.rate_limit = cfg.get("rate_limit", 1.0)
        self.client = httpx.Client(
            base_url=self.base_url,
            timeout=30.0,
            headers={"User-Agent": "biomedical-mcp/0.3"},
        )

    def _rate_wait(self) -> None:
        """Wait for rate limit on this client's host."""
        host = httpx.URL(self.base_url).host
        _wait_rate_limit(host, self.rate_limit)

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=2, max=15),
           retry=retry_if_exception(_is_retryable))
    def _get(self, path: str, params: dict | None = None, headers: dict | None = None) -> dict | list:
        self._rate_wait()
        resp = self.client.get(path, params=params, headers=headers or {})
        resp.raise_for_status()
        return resp.json()

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=2, max=15),
           retry=retry_if_exception(_is_retryable))
    def _post(self, path: str, json_data: dict | None = None, data: dict | None = None) -> dict | list:
        self._rate_wait()
        resp = self.client.post(path, json=json_data, data=data)
        resp.raise_for_status()
        return resp.json()

    def _cache_key(self, prefix: str, *parts: str) -> str:
        raw = f"{self.domain_name}:{prefix}:{':'.join(str(p) for p in parts)}"
        return hashlib.md5(raw.encode()).hexdigest()

    def _cached_get(self, prefix: str, path: str, params: dict | None = None,
                    ttl_days: int | None = None, headers: dict | None = None) -> tuple[dict | list, bool]:
        """Cached GET. Returns (data, cache_hit)."""
        raw = f"{path}:{sorted(params.items()) if params else ''}"
        key = self._cache_key(prefix, raw)
        cached = self.cache.get(key, max_age_days=ttl_days or self.ttl_days)
        if cached is not None:
            return cached, True
        result = self._get(path, params, headers=headers)
        self.cache.set(key, result)
        return result, False

    def _provenance(self, *, cache_hit: bool = False, data_version: str | None = None,
                    evidence_grade: str | None = None) -> dict:
        return make_provenance(
            self.domain_name,
            data_version=data_version,
            evidence_grade=evidence_grade,
            cache_hit=cache_hit,
        )

    def validate(self) -> bool:
        """Health check — override in subclass. Returns True if API is reachable."""
        return True
