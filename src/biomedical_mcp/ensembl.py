"""Ensembl REST API — gene, transcript, and regulatory annotation."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

ENSEMBL_URL = "https://rest.ensembl.org"


def _transient(exc: BaseException) -> bool:
    return isinstance(exc, (httpx.TimeoutException, httpx.HTTPStatusError))


class Ensembl:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(timeout=30.0, headers={"Content-Type": "application/json"})

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10),
           retry=retry_if_exception(_transient))
    def _get(self, path: str, params: dict | None = None) -> dict | list:
        resp = self.client.get(f"{ENSEMBL_URL}{path}", params=params or {})
        resp.raise_for_status()
        return resp.json()

    def gene_by_symbol(self, symbol: str, species: str = "homo_sapiens") -> dict:
        """Look up gene by symbol — returns Ensembl ID, biotype, location, description."""
        key = hashlib.md5(f"ensembl:symbol:{symbol}:{species}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=30)
        if cached is not None:
            return cached

        try:
            data = self._get(f"/lookup/symbol/{species}/{symbol}", {"expand": "1"})
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 400:
                return {"error": f"Gene not found: {symbol}"}
            raise

        result = {
            "ensembl_id": data.get("id"),
            "symbol": data.get("display_name"),
            "description": data.get("description"),
            "biotype": data.get("biotype"),
            "chromosome": data.get("seq_region_name"),
            "start": data.get("start"),
            "end": data.get("end"),
            "strand": data.get("strand"),
            "assembly": data.get("assembly_name"),
            "source": data.get("source"),
        }
        self.cache.set(key, result)
        return result

    def gene_by_id(self, ensembl_id: str) -> dict:
        """Look up gene by Ensembl ID."""
        key = hashlib.md5(f"ensembl:id:{ensembl_id}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=30)
        if cached is not None:
            return cached

        try:
            data = self._get(f"/lookup/id/{ensembl_id}", {"expand": "1"})
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 400:
                return {"error": f"ID not found: {ensembl_id}"}
            raise

        result = {
            "ensembl_id": data.get("id"),
            "symbol": data.get("display_name"),
            "description": data.get("description"),
            "biotype": data.get("biotype"),
            "chromosome": data.get("seq_region_name"),
            "start": data.get("start"),
            "end": data.get("end"),
            "strand": data.get("strand"),
            "assembly": data.get("assembly_name"),
        }
        self.cache.set(key, result)
        return result

    def xrefs(self, ensembl_id: str) -> dict:
        """Get cross-references to external databases (UniProt, HGNC, RefSeq, etc.)."""
        key = hashlib.md5(f"ensembl:xrefs:{ensembl_id}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=30)
        if cached is not None:
            return cached

        data = self._get(f"/xrefs/id/{ensembl_id}")
        xrefs = {}
        for x in data:
            db = x.get("dbname", "unknown")
            if db not in xrefs:
                xrefs[db] = []
            xrefs[db].append({
                "primary_id": x.get("primary_id"),
                "display_id": x.get("display_id"),
                "description": x.get("description"),
            })

        result = {"ensembl_id": ensembl_id, "xrefs": xrefs}
        self.cache.set(key, result)
        return result

    def sequence(self, ensembl_id: str, seq_type: str = "genomic") -> dict:
        """Get sequence for a gene/transcript. seq_type: genomic, cds, cdna, protein."""
        key = hashlib.md5(f"ensembl:seq:{ensembl_id}:{seq_type}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=30)
        if cached is not None:
            return cached

        data = self._get(f"/sequence/id/{ensembl_id}", {"type": seq_type})
        result = {
            "ensembl_id": data.get("id"),
            "molecule": data.get("molecule"),
            "seq_type": seq_type,
            "length": len(data.get("seq", "")),
            "sequence": data.get("seq", "")[:5000],  # cap to avoid huge responses
        }
        self.cache.set(key, result)
        return result

    def validate(self) -> bool:
        """Health check — ping Ensembl REST API."""
        try:
            resp = self.client.get(f"{ENSEMBL_URL}/info/ping")
            return resp.status_code == 200
        except Exception:
            return False
