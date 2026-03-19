"""STRING DB — protein-protein interaction networks."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

STRING_URL = "https://string-db.org/api/json"
HUMAN_TAXID = "9606"


def _transient(exc: BaseException) -> bool:
    return isinstance(exc, (httpx.TimeoutException, httpx.HTTPStatusError))


class StringDB:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(timeout=30.0)

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10),
           retry=retry_if_exception(_transient))
    def _get(self, endpoint: str, params: dict) -> list | dict:
        url = f"{STRING_URL}/{endpoint}"
        resp = self.client.get(url, params=params)
        resp.raise_for_status()
        return resp.json()

    def resolve_id(self, protein: str, species: str = HUMAN_TAXID) -> str | None:
        """Resolve protein name to STRING identifier."""
        key = hashlib.md5(f"string:resolve:{protein}:{species}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=30)
        if cached is not None:
            return cached.get("string_id")

        results = self._get("get_string_ids", {
            "identifiers": protein, "species": species, "limit": 1,
        })
        if results and isinstance(results, list):
            string_id = results[0].get("stringId")
            self.cache.set(key, {"string_id": string_id, "results": results})
            return string_id
        return None

    def interactions(self, protein: str, species: str = HUMAN_TAXID,
                     min_score: int = 700, limit: int = 25) -> dict:
        """Get protein-protein interactions."""
        key = hashlib.md5(f"string:interactions:{protein}:{species}:{min_score}:{limit}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=7)
        if cached is not None:
            return cached

        string_id = self.resolve_id(protein, species)
        if not string_id:
            return {"error": f"Could not resolve protein: {protein}"}

        results = self._get("interaction_partners", {
            "identifiers": string_id, "species": species,
            "required_score": min_score, "limit": limit,
        })
        interactions = []
        for r in results:
            interactions.append({
                "partner_a": r.get("preferredName_A", ""),
                "partner_b": r.get("preferredName_B", ""),
                "combined_score": r.get("score", 0),
                "experimental": r.get("escore", 0),
                "coexpression": r.get("ascore", 0),
                "database": r.get("dscore", 0),
                "textmining": r.get("tscore", 0),
            })

        result = {"protein": protein, "species": species, "count": len(interactions),
                  "interactions": interactions}
        self.cache.set(key, result)
        return result

    def functional_enrichment(self, proteins: list[str], species: str = HUMAN_TAXID) -> dict:
        """Get functional enrichment (GO, KEGG, Reactome) for a set of proteins."""
        proteins_key = ",".join(sorted(proteins))
        key = hashlib.md5(f"string:enrichment:{proteins_key}:{species}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=7)
        if cached is not None:
            return cached

        results = self._get("enrichment", {
            "identifiers": "\r".join(proteins), "species": species,
        })
        enrichment = {}
        for r in results:
            category = r.get("category", "unknown")
            if category not in enrichment:
                enrichment[category] = []
            enrichment[category].append({
                "term": r.get("term", ""),
                "description": r.get("description", ""),
                "p_value": r.get("p_value", 1.0),
                "fdr": r.get("fdr", 1.0),
                "gene_count": r.get("number_of_genes", 0),
                "genes": r.get("preferredNames", ""),
            })

        result = {"proteins": proteins, "species": species, "categories": enrichment}
        self.cache.set(key, result)
        return result
