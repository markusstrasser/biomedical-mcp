"""AlphaFold DB — predicted protein structure and confidence data."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/api"


def _transient(exc: BaseException) -> bool:
    return isinstance(exc, (httpx.TimeoutException, httpx.HTTPStatusError))


class AlphaFold:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(timeout=30.0)

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10),
           retry=retry_if_exception(_transient))
    def _get(self, path: str) -> dict | list:
        resp = self.client.get(f"{ALPHAFOLD_URL}{path}")
        resp.raise_for_status()
        return resp.json()

    def prediction(self, uniprot_id: str) -> dict:
        """Get AlphaFold structure prediction for a UniProt accession."""
        key = hashlib.md5(f"alphafold:prediction:{uniprot_id}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=30)
        if cached is not None:
            return cached

        try:
            data = self._get(f"/prediction/{uniprot_id}")
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return {"error": f"No AlphaFold prediction for: {uniprot_id}"}
            raise

        # API returns a list, take first entry
        entry = data[0] if isinstance(data, list) and data else data
        result = {
            "uniprot_id": entry.get("uniprotAccession"),
            "gene": entry.get("gene"),
            "organism": entry.get("organismScientificName"),
            "model_url": entry.get("pdbUrl"),
            "cif_url": entry.get("cifUrl"),
            "pae_image_url": entry.get("paeImageUrl"),
            "model_confidence_url": entry.get("bcifUrl"),
            "mean_plddt": entry.get("globalMetricValue"),
            "model_version": entry.get("latestVersion"),
            "sequence_length": entry.get("uniprotEnd"),
        }
        self.cache.set(key, result)
        return result

    def pae(self, uniprot_id: str) -> dict:
        """Get predicted aligned error (PAE) summary for a protein structure."""
        key = hashlib.md5(f"alphafold:pae:{uniprot_id}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=30)
        if cached is not None:
            return cached

        try:
            data = self._get(f"/prediction/{uniprot_id}")
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return {"error": f"No AlphaFold prediction for: {uniprot_id}"}
            raise

        entry = data[0] if isinstance(data, list) and data else data
        result = {
            "uniprot_id": entry.get("uniprotAccession"),
            "pae_doc_url": entry.get("paeDocUrl"),
            "pae_image_url": entry.get("paeImageUrl"),
            "mean_plddt": entry.get("globalMetricValue"),
            "confidence_version": entry.get("confidenceVersion"),
        }
        self.cache.set(key, result)
        return result
