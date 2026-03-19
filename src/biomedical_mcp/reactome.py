"""Reactome — biological pathway analysis and lookup."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

REACTOME_URL = "https://reactome.org"
CONTENT_URL = f"{REACTOME_URL}/ContentService"
ANALYSIS_URL = f"{REACTOME_URL}/AnalysisService"


def _transient(exc: BaseException) -> bool:
    return isinstance(exc, (httpx.TimeoutException, httpx.HTTPStatusError))


class Reactome:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(timeout=30.0, headers={"Accept": "application/json"})

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10),
           retry=retry_if_exception(_transient))
    def _get(self, url: str, params: dict | None = None) -> dict | list:
        resp = self.client.get(url, params=params or {})
        resp.raise_for_status()
        return resp.json()

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10),
           retry=retry_if_exception(_transient))
    def _post(self, url: str, data: str, content_type: str = "text/plain") -> dict:
        resp = self.client.post(url, content=data, headers={"Content-Type": content_type})
        resp.raise_for_status()
        return resp.json()

    def pathways_for_gene(self, gene_symbol: str, species: str = "Homo sapiens") -> dict:
        """Get Reactome pathways containing a gene."""
        key = hashlib.md5(f"reactome:gene_pathways:{gene_symbol}:{species}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=14)
        if cached is not None:
            return cached

        try:
            data = self._get(f"{CONTENT_URL}/data/mapping/UniProt/{gene_symbol}/pathways")
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                # Try via gene symbol directly
                try:
                    data = self._get(f"{CONTENT_URL}/search/query",
                                     {"query": gene_symbol, "types": "Pathway", "species": species})
                    entries = data.get("results", [{}])[0].get("entries", []) if data.get("results") else []
                    pathways = [{
                        "stable_id": e.get("stId", ""),
                        "name": e.get("name", ""),
                        "species": species,
                    } for e in entries[:25]]
                    result = {"gene": gene_symbol, "count": len(pathways), "pathways": pathways}
                    self.cache.set(key, result)
                    return result
                except Exception:
                    return {"error": f"No pathways found for: {gene_symbol}"}
            raise

        pathways = [{
            "stable_id": p.get("stId", ""),
            "name": p.get("displayName", ""),
            "species": p.get("speciesName", ""),
            "is_disease": p.get("isDisease", False),
        } for p in data if isinstance(p, dict)]

        result = {"gene": gene_symbol, "count": len(pathways), "pathways": pathways[:50]}
        self.cache.set(key, result)
        return result

    def pathway_detail(self, pathway_id: str) -> dict:
        """Get detailed info for a Reactome pathway."""
        key = hashlib.md5(f"reactome:pathway:{pathway_id}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=14)
        if cached is not None:
            return cached

        data = self._get(f"{CONTENT_URL}/data/query/{pathway_id}")
        result = {
            "stable_id": data.get("stId"),
            "name": data.get("displayName"),
            "species": data.get("speciesName"),
            "is_disease": data.get("isDisease", False),
            "summation": [s.get("text", "") for s in data.get("summation", [])],
            "compartment": [c.get("displayName", "") for c in data.get("compartment", [])],
            "has_diagram": data.get("hasDiagram", False),
            "diagram_url": f"{REACTOME_URL}/PathwayBrowser/#/{pathway_id}" if data.get("hasDiagram") else None,
        }
        self.cache.set(key, result)
        return result

    def enrichment(self, gene_list: list[str], species: str = "Homo sapiens") -> dict:
        """Run pathway enrichment analysis on a gene list."""
        genes_str = "\n".join(gene_list)
        key = hashlib.md5(f"reactome:enrichment:{genes_str}:{species}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=3)
        if cached is not None:
            return cached

        data = self._post(
            f"{ANALYSIS_URL}/identifiers/projection?pageSize=25&page=1",
            genes_str,
        )
        pathways = []
        for p in data.get("pathways", []):
            entities = p.get("entities", {})
            pathways.append({
                "stable_id": p.get("stId", ""),
                "name": p.get("name", ""),
                "p_value": entities.get("pValue"),
                "fdr": entities.get("fdr"),
                "found": entities.get("found", 0),
                "total": entities.get("total", 0),
                "ratio": entities.get("ratio", 0),
            })

        result = {
            "input_genes": len(gene_list),
            "found_identifiers": data.get("identifiersNotFound", 0),
            "pathways_found": len(pathways),
            "pathways": pathways,
        }
        self.cache.set(key, result)
        return result
