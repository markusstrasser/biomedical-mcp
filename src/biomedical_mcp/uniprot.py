"""UniProt REST API client — protein function, domains, variants, and cross-references."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

UNIPROT_BASE = "https://rest.uniprot.org"

PRIORITY_FEATURE_TYPES = {
    "Domain", "Natural variant", "Mutagenesis", "Active site",
    "Binding site", "Modified residue",
}


def _is_retryable(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 500, 502, 503, 504)
    return isinstance(exc, (httpx.ConnectError, httpx.ReadTimeout))


class UniProt:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(
            base_url=UNIPROT_BASE,
            timeout=30,
            headers={"Accept": "application/json", "User-Agent": "biomedical-mcp/0.1"},
        )

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=2, max=15), retry=retry_if_exception(_is_retryable))
    def _get(self, path: str, params: dict | None = None) -> dict | list:
        resp = self.client.get(path, params=params)
        resp.raise_for_status()
        return resp.json()

    def _cached(self, prefix: str, path: str, params: dict | None = None, max_age_days: int = 30) -> dict | list:
        raw = f"{path}:{sorted(params.items()) if params else ''}"
        key = f"uniprot:{prefix}:{hashlib.md5(raw.encode()).hexdigest()}"
        cached = self.cache.get(key, max_age_days=max_age_days)
        if cached is not None:
            return cached
        result = self._get(path, params)
        self.cache.set(key, result)
        return result

    def _resolve_accession(self, accession: str | None = None, gene_symbol: str | None = None) -> str | None:
        """Resolve a gene symbol to a UniProt accession (human, reviewed)."""
        if accession:
            return accession
        if not gene_symbol:
            return None
        data = self._cached(
            "resolve", "/uniprotkb/search",
            {"query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true", "format": "json", "size": "1"},
        )
        results = data.get("results", [])
        if not results:
            return None
        return results[0].get("primaryAccession")

    def protein(self, accession: str | None = None, gene_symbol: str | None = None) -> dict | None:
        """Get protein details: function, domains, disease associations, cross-refs."""
        acc = self._resolve_accession(accession, gene_symbol)
        if not acc:
            return None
        data = self._cached("protein", f"/uniprotkb/{acc}", {"format": "json"})
        return self._extract_protein(data)

    def variants(self, accession: str | None = None, gene_symbol: str | None = None, limit: int = 50) -> dict | None:
        """Get known natural variants and mutagenesis data."""
        acc = self._resolve_accession(accession, gene_symbol)
        if not acc:
            return None
        data = self._cached("protein", f"/uniprotkb/{acc}", {"format": "json"})
        features = data.get("features", [])
        variant_features = [
            f for f in features
            if f.get("type") in ("Natural variant", "Mutagenesis")
        ][:limit]
        return {
            "accession": acc,
            "gene": self._get_gene_name(data),
            "variant_count": len(variant_features),
            "variants": [self._extract_feature(f) for f in variant_features],
        }

    def search(self, query: str, organism: str = "human", limit: int = 10) -> list[dict]:
        """Search proteins by keyword."""
        org_id = "9606" if organism.lower() == "human" else organism
        data = self._cached(
            "search", "/uniprotkb/search",
            {"query": f"{query} AND organism_id:{org_id}", "format": "json", "size": str(min(limit, 50))},
        )
        results = data.get("results", [])
        return [
            {
                "accession": r.get("primaryAccession"),
                "gene": self._get_gene_name(r),
                "protein_name": self._get_protein_name(r),
                "organism": r.get("organism", {}).get("scientificName"),
                "length": r.get("sequence", {}).get("length"),
            }
            for r in results
        ]

    def _extract_protein(self, data: dict) -> dict:
        """Extract structured protein summary from raw UniProt response."""
        result = {
            "accession": data.get("primaryAccession"),
            "gene": self._get_gene_name(data),
            "protein_name": self._get_protein_name(data),
            "organism": data.get("organism", {}).get("scientificName"),
            "length": data.get("sequence", {}).get("length"),
        }

        # Function from comments
        comments = data.get("comments", [])
        for c in comments:
            ctype = c.get("commentType")
            if ctype == "FUNCTION":
                texts = c.get("texts", [])
                result["function"] = texts[0].get("value", "")[:2000] if texts else None
            elif ctype == "SUBCELLULAR LOCATION":
                locs = c.get("subcellularLocations", [])
                result["subcellular_location"] = [
                    loc.get("location", {}).get("value")
                    for loc in locs if loc.get("location")
                ][:10]
            elif ctype == "DISEASE":
                disease = c.get("disease", {})
                if disease:
                    result.setdefault("diseases", []).append({
                        "name": disease.get("diseaseId"),
                        "description": disease.get("description", "")[:500],
                        "acronym": disease.get("acronym"),
                    })

        # Priority features
        features = data.get("features", [])
        priority = [f for f in features if f.get("type") in PRIORITY_FEATURE_TYPES]
        result["domains"] = [
            self._extract_feature(f) for f in priority
            if f.get("type") == "Domain"
        ][:20]
        result["active_sites"] = [
            self._extract_feature(f) for f in priority
            if f.get("type") in ("Active site", "Binding site")
        ][:10]

        # GO terms
        dbrefs = data.get("uniProtKBCrossReferences", [])
        go_terms = [
            {"id": ref.get("id"), "term": self._go_term_from_props(ref.get("properties", []))}
            for ref in dbrefs if ref.get("database") == "GO"
        ][:20]
        if go_terms:
            result["go_terms"] = go_terms

        # Cross-references summary (count by database)
        db_counts = {}
        for ref in dbrefs:
            db = ref.get("database")
            if db:
                db_counts[db] = db_counts.get(db, 0) + 1
        result["cross_ref_databases"] = dict(sorted(db_counts.items(), key=lambda x: -x[1])[:20])

        return result

    def _extract_feature(self, feature: dict) -> dict:
        loc = feature.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        return {
            "type": feature.get("type"),
            "description": feature.get("description", "")[:300],
            "start": start,
            "end": end,
            "alternativeSequence": feature.get("alternativeSequence", {}).get("originalSequence"),
        }

    @staticmethod
    def _get_gene_name(data: dict) -> str | None:
        genes = data.get("genes", [])
        if genes:
            return genes[0].get("geneName", {}).get("value")
        return None

    @staticmethod
    def _get_protein_name(data: dict) -> str | None:
        desc = data.get("proteinDescription", {})
        rec = desc.get("recommendedName", {})
        if rec:
            return rec.get("fullName", {}).get("value")
        sub = desc.get("submissionNames", [])
        if sub:
            return sub[0].get("fullName", {}).get("value")
        return None

    @staticmethod
    def _go_term_from_props(props: list) -> str | None:
        for p in props:
            if p.get("key") == "GoTerm":
                return p.get("value")
        return None
