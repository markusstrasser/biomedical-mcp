"""Orphanet (Orphadata) REST API client — rare disease data."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


def _extract_results(data: dict | list) -> dict | list:
    """Orphadata wraps responses in {data: {results: ...}} — unwrap."""
    if isinstance(data, dict):
        inner = data.get("data", data)
        if isinstance(inner, dict):
            return inner.get("results", inner)
    return data


class Orphanet(BaseClient):
    domain_name = "orphanet"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def gene_diseases(self, gene_symbol: str) -> dict:
        """Get rare diseases associated with a gene."""
        path = f"/rd-associated-genes/genes/names/{gene_symbol}"
        data, cache_hit = self._cached_get("gene", path)
        results = _extract_results(data)

        diseases = []
        if isinstance(results, list):
            for entry in results:
                for assoc in entry.get("DisorderGeneAssociation", []):
                    disorder = assoc.get("Disorder", entry)
                    diseases.append({
                        "orphacode": disorder.get("ORPHAcode") or entry.get("ORPHAcode"),
                        "name": disorder.get("Preferred term") or disorder.get("Name"),
                        "association_type": assoc.get("DisorderGeneAssociationType"),
                        "association_status": assoc.get("DisorderGeneAssociationStatus"),
                    })
        elif isinstance(results, dict):
            for assoc in results.get("DisorderGeneAssociation", []):
                disorder = assoc.get("Disorder", {})
                diseases.append({
                    "orphacode": disorder.get("ORPHAcode") or results.get("ORPHAcode"),
                    "name": disorder.get("Preferred term") or disorder.get("Name"),
                    "association_type": assoc.get("DisorderGeneAssociationType"),
                    "association_status": assoc.get("DisorderGeneAssociationStatus"),
                })

        return {
            "gene_symbol": gene_symbol,
            "count": len(diseases),
            "diseases": diseases,
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="B2"),
        }

    def disease_genes(self, orphacode: int) -> dict:
        """Get genes associated with a rare disease by ORPHAcode."""
        path = f"/rd-associated-genes/orphacodes/{orphacode}"
        data, cache_hit = self._cached_get("disease_genes", path)
        results = _extract_results(data)

        genes = []
        if isinstance(results, dict):
            for assoc in results.get("DisorderGeneAssociation", []):
                gene = assoc.get("Gene", {})
                genes.append({
                    "symbol": gene.get("Symbol"),
                    "name": gene.get("name"),
                    "gene_type": gene.get("GeneType"),
                    "locus": [l.get("GeneLocus") for l in gene.get("Locus", [])],
                    "association_type": assoc.get("DisorderGeneAssociationType"),
                    "association_status": assoc.get("DisorderGeneAssociationStatus"),
                    "external_refs": {
                        ref.get("Source"): ref.get("Reference")
                        for ref in gene.get("ExternalReference", [])
                    },
                })

        disease_name = results.get("Preferred term") if isinstance(results, dict) else None

        return {
            "orphacode": orphacode,
            "disease_name": disease_name,
            "count": len(genes),
            "genes": genes,
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="B2"),
        }

    def natural_history(self, orphacode: int) -> dict:
        """Inheritance, age of onset, age of death for a rare disease."""
        path = f"/rd-natural_history/orphacodes/{orphacode}"
        data, cache_hit = self._cached_get("nathistory", path)
        results = _extract_results(data)

        if not isinstance(results, dict):
            return {"error": f"No natural history for ORPHAcode {orphacode}"}

        return {
            "orphacode": orphacode,
            "name": results.get("Preferred term"),
            "disorder_type": results.get("Typology"),
            "disorder_group": results.get("DisorderGroup"),
            "inheritance": results.get("TypeOfInheritance", []),
            "age_of_onset": results.get("AverageAgeOfOnset", []),
            "age_of_death": results.get("AverageAgeOfDeath", []),
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="B2"),
        }

    def epidemiology(self, orphacode: int) -> dict:
        """Prevalence and incidence data for a rare disease."""
        path = f"/rd-epidemiology/orphacodes/{orphacode}"
        data, cache_hit = self._cached_get("epidemiology", path)
        results = _extract_results(data)

        if not isinstance(results, dict):
            return {"error": f"No epidemiology for ORPHAcode {orphacode}"}

        return {
            "orphacode": orphacode,
            "name": results.get("Preferred term"),
            "prevalence": results.get("Prevalence", []),
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="B2"),
        }

    def phenotypes(self, orphacode: int) -> dict:
        """HPO phenotypes associated with a rare disease."""
        path = f"/rd-phenotypes/orphacodes/{orphacode}"
        data, cache_hit = self._cached_get("phenotypes", path)
        results = _extract_results(data)

        if not isinstance(results, dict):
            return {"error": f"No phenotypes for ORPHAcode {orphacode}"}

        hpo_list = results.get("HPODisorderAssociation", [])
        phenotypes = [
            {
                "hpo_id": p.get("HPO", {}).get("HPOId"),
                "hpo_term": p.get("HPO", {}).get("HPOTerm"),
                "frequency": p.get("HPOFrequency"),
            }
            for p in hpo_list
        ]

        return {
            "orphacode": orphacode,
            "name": results.get("Preferred term"),
            "count": len(phenotypes),
            "phenotypes": phenotypes,
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="B2"),
        }

    def cross_reference(self, orphacode: int) -> dict:
        """Cross-references to ICD-10, ICD-11, OMIM for a rare disease."""
        path = f"/rd-cross-referencing/orphacodes/{orphacode}"
        data, cache_hit = self._cached_get("xref", path)
        results = _extract_results(data)

        if not isinstance(results, dict):
            return {"error": f"No cross-references for ORPHAcode {orphacode}"}

        return {
            "orphacode": orphacode,
            "name": results.get("Preferred term"),
            "external_references": results.get("ExternalReference", []),
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="B2"),
        }

    def search_by_omim(self, omim_code: str) -> dict:
        """Find rare diseases by OMIM code (disease OMIM, not gene OMIM)."""
        path = f"/rd-cross-referencing/omims/{omim_code}"
        try:
            data, cache_hit = self._cached_get("omim", path)
        except Exception as e:
            if "404" in str(e):
                return {
                    "omim_code": omim_code,
                    "count": 0,
                    "diseases": [],
                    "note": "OMIM code not found in Orphanet. Use disease OMIM codes (not gene OMIM).",
                    "provenance": self._provenance(evidence_grade="B2"),
                }
            raise
        results = _extract_results(data)

        diseases = []
        if isinstance(results, list):
            for r in results:
                diseases.append({
                    "orphacode": r.get("ORPHAcode"),
                    "name": r.get("Preferred term"),
                })
        elif isinstance(results, dict):
            diseases.append({
                "orphacode": results.get("ORPHAcode"),
                "name": results.get("Preferred term"),
            })

        return {
            "omim_code": omim_code,
            "count": len(diseases),
            "diseases": diseases,
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="B2"),
        }

    def validate(self) -> bool:
        try:
            self._get("/rd-associated-genes/genes", params={"page": 1})
            return True
        except Exception:
            return False
