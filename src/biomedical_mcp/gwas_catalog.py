"""GWAS Catalog (EBI) REST API client — variant-trait associations."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class GWASCatalog(BaseClient):
    domain_name = "gwas_catalog"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def variant_associations(self, rsid: str, limit: int = 25) -> dict:
        """GWAS associations for a variant (rsID)."""
        key = self._cache_key("var_assoc", rsid)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        try:
            data = self._get(
                f"/singleNucleotidePolymorphisms/{rsid}/associations",
                params={"projection": "associationBySnp"},
            )
        except Exception as exc:
            import httpx as _httpx
            if isinstance(exc, _httpx.HTTPStatusError) and exc.response.status_code == 404:
                return {"rsid": rsid, "association_count": 0, "associations": [],
                        "note": "Variant not found in GWAS Catalog"}
            log.warning("GWAS Catalog API error for %s: %s", rsid, exc)
            return {"rsid": rsid, "association_count": 0, "associations": [],
                    "error": "api_error", "note": f"GWAS Catalog API error: {type(exc).__name__}"}
        associations = self._extract_associations(data, limit)
        result = {
            "rsid": rsid,
            "association_count": len(associations),
            "associations": associations,
        }
        self.cache.set(key, result)
        return result

    def gene_associations(self, gene_symbol: str, limit: int = 25) -> dict:
        """GWAS associations for a gene."""
        key = self._cache_key("gene_assoc", gene_symbol)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        try:
            data = self._get(f"/genes/{gene_symbol}/associations")
        except Exception:
            return {"gene_symbol": gene_symbol, "association_count": 0, "associations": [],
                    "note": "No GWAS associations found (gene may not be in GWAS Catalog)"}
        associations = self._extract_associations(data, limit)
        result = {
            "gene_symbol": gene_symbol,
            "association_count": len(associations),
            "associations": associations,
        }
        self.cache.set(key, result)
        return result

    def trait_search(self, query: str, limit: int = 10) -> dict:
        """Search EFO traits by keyword."""
        data, _ = self._cached_get(
            "trait_search",
            "/efoTraits",
            params={"search": query, "size": min(limit, 50)},
        )
        traits = []
        for item in (data.get("_embedded", {}).get("efoTraits", []) if isinstance(data, dict) else []):
            traits.append({
                "trait": item.get("trait"),
                "short_form": item.get("shortForm"),
                "uri": item.get("uri"),
            })
        return {"query": query, "count": len(traits), "traits": traits[:limit]}

    def _extract_associations(self, data: dict | list, limit: int) -> list[dict]:
        """Extract associations from GWAS Catalog response."""
        raw = []
        if isinstance(data, dict):
            raw = data.get("_embedded", {}).get("associations", [])
        elif isinstance(data, list):
            raw = data

        associations = []
        for assoc in raw[:limit]:
            strongest = assoc.get("strongestRiskAlleles", [{}])
            risk_allele = strongest[0].get("riskAlleleName") if strongest else None
            loci = assoc.get("loci", [{}])
            genes = []
            for locus in loci:
                for gene_entry in locus.get("authorReportedGenes", []):
                    genes.append(gene_entry.get("geneName"))

            # Traits from efoTraits array (HAL response format)
            efo_traits = assoc.get("efoTraits", [])
            trait_names = [t.get("trait") for t in efo_traits if t.get("trait")]

            # Study accession from _links (not available as top-level dict in projection)
            study_accession = None
            study_obj = assoc.get("study")
            if isinstance(study_obj, dict):
                study_accession = study_obj.get("accessionId")
            elif isinstance(assoc.get("_links", {}), dict):
                study_href = assoc.get("_links", {}).get("study", {}).get("href", "")
                if "/studies/" in study_href:
                    study_accession = study_href.rsplit("/studies/", 1)[-1]

            associations.append({
                "pvalue": assoc.get("pvalueMantissa"),
                "pvalue_exponent": assoc.get("pvalueExponent"),
                "risk_allele": risk_allele,
                "or_beta": assoc.get("orPerCopyNum"),
                "genes": genes,
                "trait": ", ".join(trait_names) if trait_names else None,
                "traits": trait_names,
                "study": study_accession,
            })
        return associations

    def validate(self) -> bool:
        try:
            self._get("/metadata")
            return True
        except Exception:
            return False
