"""ISBT Blood Group Database REST API client — alleles, antigens, systems."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class ISBT(BaseClient):
    domain_name = "isbt"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def systems(self) -> dict:
        """Get all blood group systems with gene/antigen/allele counts."""
        data, cache_hit = self._cached_get("systems", "/system/summary")
        systems = []
        for s in (data if isinstance(data, list) else []):
            genes = s.get("genes", [])
            antigens = s.get("antigens", [])
            alleles = s.get("alleles", [])
            systems.append({
                "id": s.get("id"),
                "isbt_number": s.get("isbt_number"),
                "name": s.get("name"),
                "symbol": s.get("symbol"),
                "category": s.get("category"),
                "gene_count": len(genes),
                "antigen_count": len(antigens),
                "allele_count": len(alleles),
                "genes": [g.get("name") for g in genes],
            })
        return {
            "count": len(systems),
            "systems": systems,
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="A1"),
        }

    def system_detail(self, system_id: int) -> dict:
        """Full detail for a blood group system including genes, antigens, alleles."""
        data, cache_hit = self._cached_get("system", f"/system/{system_id}")
        if isinstance(data, dict) and "error" in data:
            return {"error": data.get("message", "System not found")}

        genes = data.get("genes", [])
        antigens = data.get("antigens", [])

        return {
            "id": data.get("id"),
            "isbt_number": data.get("isbt_number"),
            "name": data.get("name"),
            "symbol": data.get("symbol"),
            "category": data.get("category"),
            "description": data.get("description"),
            "genes": [
                {
                    "name": g.get("name"),
                    "isbt_gene": g.get("isbt_gene"),
                    "isbt_transcript": g.get("isbt_transcript"),
                    "isbt_protein": g.get("isbt_protein"),
                    "chromosome": g.get("chromosome"),
                    "exon_count": g.get("exon_count"),
                    "hg38_start": g.get("hg38_start"),
                    "hg38_end": g.get("hg38_end"),
                }
                for g in genes
            ],
            "antigens": [
                {
                    "isbt_code": a.get("isbt_code"),
                    "name": a.get("name"),
                    "display_name": a.get("display_name"),
                    "isbt_number": a.get("isbt_number"),
                }
                for a in antigens
            ],
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="A1"),
        }

    def alleles(self, system_symbol: str) -> dict:
        """Get all alleles for a blood group system by symbol (e.g. 'ABO', 'RHD')."""
        # First get all alleles, then filter by system
        data, cache_hit = self._cached_get("all_alleles", "/allele")
        all_alleles = data if isinstance(data, list) else []

        filtered = [
            a for a in all_alleles
            if a.get("system", {}).get("symbol", "").upper() == system_symbol.upper()
        ]

        alleles = []
        for a in filtered:
            alleles.append({
                "isbt_allele": a.get("isbt_allele"),
                "phenotype": a.get("isbt_phenotype"),
                "snp": a.get("isbt_snp"),
                "null_allele": a.get("null_allele", False),
                "weak_allele": a.get("weak_allele", False),
                "partial_allele": a.get("partial_allele", False),
                "group": a.get("table_row_group"),
                "gene": a.get("gene", {}).get("name"),
                "notes": a.get("notes"),
                "approved": a.get("approved"),
            })

        return {
            "system": system_symbol.upper(),
            "count": len(alleles),
            "alleles": alleles,
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="A1"),
        }

    def allele_detail(self, allele_id: int) -> dict:
        """Get full detail for a specific allele by ID."""
        data, cache_hit = self._cached_get("allele", f"/allele/{allele_id}")
        if isinstance(data, dict) and "error" in data:
            return {"error": data.get("message", "Allele not found")}

        return {
            "id": data.get("id"),
            "isbt_allele": data.get("isbt_allele"),
            "system": data.get("system", {}).get("symbol"),
            "gene": data.get("gene", {}).get("name"),
            "phenotype": data.get("isbt_phenotype"),
            "snp": data.get("isbt_snp"),
            "null_allele": data.get("null_allele"),
            "weak_allele": data.get("weak_allele"),
            "partial_allele": data.get("partial_allele"),
            "group": data.get("table_row_group"),
            "alternate_names": data.get("alternate_names"),
            "notes": data.get("notes"),
            "comment": data.get("comment"),
            "approved": data.get("approved"),
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="A1"),
        }

    def antigens(self, system_symbol: str | None = None) -> dict:
        """Get antigens, optionally filtered by system symbol."""
        data, cache_hit = self._cached_get("all_antigens", "/antigen")
        all_antigens = data if isinstance(data, list) else []

        if system_symbol:
            all_antigens = [
                a for a in all_antigens
                if a.get("system", {}).get("symbol", "").upper() == system_symbol.upper()
            ]

        antigens = [
            {
                "isbt_code": a.get("isbt_code"),
                "name": a.get("name"),
                "display_name": a.get("display_name"),
                "isbt_number": a.get("isbt_number"),
                "system": a.get("system", {}).get("symbol"),
            }
            for a in all_antigens
        ]

        return {
            "filter": system_symbol or "all",
            "count": len(antigens),
            "antigens": antigens,
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="A1"),
        }

    def search_alleles(self, isbt_allele: str | None = None,
                       system_symbol: str | None = None,
                       phenotype: str | None = None) -> dict:
        """Search alleles by ISBT allele name, system, or phenotype.

        ISBT search is field-specific (no free-text). Provide at least one filter.
        """
        params = {}
        if isbt_allele:
            params["isbt_allele"] = isbt_allele
        if system_symbol:
            params["system_symbol"] = system_symbol
        if phenotype:
            params["isbt_phenotype"] = phenotype
        if not params:
            return {"error": "Provide at least one of: isbt_allele, system_symbol, phenotype"}

        data, cache_hit = self._cached_get("search", "/allele/search", params=params)
        results = data if isinstance(data, list) else data.get("results", []) if isinstance(data, dict) else []

        alleles = [
            {
                "id": a.get("id"),
                "isbt_allele": a.get("isbt_allele"),
                "system": a.get("system", {}).get("symbol"),
                "phenotype": a.get("isbt_phenotype"),
                "snp": a.get("isbt_snp"),
                "group": a.get("table_row_group"),
                "approved": a.get("approved"),
            }
            for a in results
        ]

        return {
            "query": params,
            "count": len(alleles),
            "alleles": alleles,
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="A1"),
        }

    def variant_lookup(self, rsid: str | None = None,
                       chromosome: str | None = None, position: int | None = None,
                       ref: str | None = None, alt: str | None = None) -> dict:
        """Look up blood group variants by rsID or GRCh38 coordinates."""
        params = {}
        if rsid:
            params["rsid"] = rsid
        if chromosome:
            params["grch38_chr"] = chromosome
        if position is not None:
            params["grch38_pos"] = str(position)
        if ref:
            params["grch38_ref"] = ref
        if alt:
            params["grch38_alt"] = alt
        if not params:
            return {"error": "Provide rsid or chromosome+position"}

        data, cache_hit = self._cached_get("variant", "/variant/search", params=params)
        results = data if isinstance(data, list) else []

        variants = [
            {
                "id": v.get("id"),
                "rsid": v.get("rsid"),
                "hgvs_genomic_grch38": v.get("hgvs_genomic_grch38"),
                "hgvs_transcript": v.get("hgvs_transcript"),
                "hgvs_predicted_protein": v.get("hgvs_predicted_protein"),
                "exon": v.get("exon"),
                "gene": v.get("gene_name"),
                "grch38_chr": v.get("grch38_chr"),
                "grch38_pos": v.get("grch38_pos"),
            }
            for v in results
        ]

        return {
            "query": params,
            "count": len(variants),
            "variants": variants,
            "provenance": self._provenance(cache_hit=cache_hit, evidence_grade="A1"),
        }

    def validate(self) -> bool:
        try:
            self._get("/system")
            return True
        except Exception:
            return False
