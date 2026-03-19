"""GTEx Portal v2 REST API client — tissue expression, eQTLs."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

# Minimum TPM to report a tissue as "expressed" (noise filter per review)
MIN_TPM_THRESHOLD = 1.0


class GTEx(BaseClient):
    domain_name = "gtex"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def _resolve_gencode_id(self, gene_symbol: str) -> str | None:
        """Resolve gene symbol to GTEx gencodeId via reference API."""
        key = self._cache_key("resolve", gene_symbol)
        cached = self.cache.get(key, max_age_days=30)
        if cached is not None:
            return cached.get("gencodeId")

        data = self._get("/reference/gene", params={"geneId": gene_symbol})
        genes = data.get("data", []) if isinstance(data, dict) else []
        if not genes:
            return None
        gencode_id = genes[0].get("gencodeId")
        self.cache.set(key, {"gencodeId": gencode_id})
        return gencode_id

    def gene_expression(self, gene_symbol: str | None = None,
                        ensembl_id: str | None = None) -> dict:
        """Median TPM by tissue for a gene. Filters tissues below TPM 1.0 by default."""
        if not gene_symbol and not ensembl_id:
            return {"error": "Provide gene_symbol or ensembl_id"}

        # GTEx v2 requires gencodeId (versioned Ensembl ID)
        gencode_id = ensembl_id
        if gene_symbol and not gencode_id:
            gencode_id = self._resolve_gencode_id(gene_symbol)
        if not gencode_id:
            return {"error": f"Could not resolve gene: {gene_symbol or ensembl_id}"}

        key = self._cache_key("expr", gencode_id)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        data = self._get("/expression/medianGeneExpression",
                         params={"gencodeId": gencode_id, "datasetId": "gtex_v8"})

        raw_items = data.get("data", []) if isinstance(data, dict) else data if isinstance(data, list) else []
        tissues = []
        for entry in raw_items:
            tpm = entry.get("median", 0)
            if tpm >= MIN_TPM_THRESHOLD:
                tissues.append({
                    "tissue_id": entry.get("tissueSiteDetailId"),
                    "median_tpm": round(tpm, 2),
                })

        tissues.sort(key=lambda t: t["median_tpm"], reverse=True)
        result = {
            "gene": gene_symbol or ensembl_id,
            "gencode_id": gencode_id,
            "expressed_tissue_count": len(tissues),
            "min_tpm_threshold": MIN_TPM_THRESHOLD,
            "tissues": tissues,
        }
        self.cache.set(key, result)
        return result

    def top_tissues(self, gene_symbol: str, limit: int = 10) -> dict:
        """Top N tissues by expression level for a gene."""
        expr = self.gene_expression(gene_symbol=gene_symbol)
        if "error" in expr:
            return expr
        tissues = expr.get("tissues", [])[:limit]
        return {
            "gene": gene_symbol,
            "top_tissues": tissues,
        }

    def eqtl(self, gene_symbol: str, tissue: str | None = None) -> dict:
        """Significant single-tissue eQTLs for a gene."""
        gencode_id = self._resolve_gencode_id(gene_symbol)
        if not gencode_id:
            return {"error": f"Could not resolve gene: {gene_symbol}"}

        params = {"gencodeId": gencode_id}
        if tissue:
            params["tissueSiteDetailId"] = tissue

        key = self._cache_key("eqtl", gene_symbol, tissue or "all")
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        data = self._get("/association/singleTissueEqtl", params=params)
        eqtls = data.get("singleTissueEqtl") or (data if isinstance(data, list) else [])

        result = {
            "gene": gene_symbol,
            "tissue_filter": tissue,
            "eqtl_count": len(eqtls),
            "eqtls": [
                {
                    "variant_id": e.get("variantId"),
                    "tissue": e.get("tissueSiteDetailId"),
                    "pvalue": e.get("pValue"),
                    "nes": e.get("nes"),
                }
                for e in eqtls[:50]
            ],
        }
        self.cache.set(key, result)
        return result

    def validate(self) -> bool:
        try:
            self._get("/", params={})
            return True
        except Exception:
            return False
