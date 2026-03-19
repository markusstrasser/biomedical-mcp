"""gnomAD GraphQL API client — population frequencies, constraint, coverage."""

from __future__ import annotations

import logging

import httpx

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

DATASET = "gnomad_r4"


def _variant_query(variant_id: str, dataset: str) -> str:
    return f"""{{
  variant(variantId: "{variant_id}", dataset: {dataset}) {{
    variant_id rsids
    genome {{ ac an af populations {{ id ac an }} }}
    exome {{ ac an af populations {{ id ac an }} }}
  }}
}}"""


def _constraint_query(gene_symbol: str) -> str:
    return f"""{{
  gene(gene_symbol: "{gene_symbol}", reference_genome: GRCh38) {{
    gene_id symbol
    gnomad_constraint {{
      pLI oe_lof oe_lof_lower oe_lof_upper
      oe_mis oe_mis_lower oe_mis_upper
      mis_z syn_z lof_z
    }}
  }}
}}"""


class GnomAD(BaseClient):
    domain_name = "gnomad"

    def __init__(self, cache: Cache):
        super().__init__(cache)
        # gnomAD uses POST to a single GraphQL endpoint
        self.client = httpx.Client(
            timeout=30.0,
            headers={"User-Agent": "biomedical-mcp/0.3", "Content-Type": "application/json"},
        )

    def _graphql(self, query: str, variables: dict | None = None) -> dict:
        """Execute a gnomAD GraphQL query."""
        self._rate_wait()
        payload: dict = {"query": query}
        if variables:
            payload["variables"] = variables
        resp = self.client.post(
            f"{self.base_url}",
            json=payload,
        )
        resp.raise_for_status()
        data = resp.json()
        if "errors" in data:
            log.warning("gnomAD GraphQL errors: %s", data["errors"])
        return data.get("data", {})

    def variant_frequency(self, variant_id: str, dataset: str = DATASET) -> dict:
        """Population allele frequencies for a variant.

        Args:
            variant_id: gnomAD-style variant ID (e.g. "1-55516888-G-A") or rsID.
        """
        key = self._cache_key("var", variant_id, dataset)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        data = self._graphql(_variant_query(variant_id, dataset), {})
        variant = data.get("variant")
        if not variant:
            return {"variant_id": variant_id, "error": "not_found"}

        result = {
            "variant_id": variant.get("variant_id"),
            "rsids": variant.get("rsids", []),
            "dataset": dataset,
        }

        for source in ("genome", "exome"):
            src_data = variant.get(source)
            if src_data:
                pops = {}
                for pop in (src_data.get("populations") or []):
                    pop_ac = pop.get("ac", 0)
                    pop_an = pop.get("an", 0)
                    pops[pop["id"]] = round(pop_ac / pop_an, 6) if pop_an > 0 else 0
                result[source] = {
                    "af": src_data.get("af"),
                    "ac": src_data.get("ac"),
                    "an": src_data.get("an"),
                    "populations": pops,
                }

        self.cache.set(key, result)
        return result

    def gene_constraint(self, gene_symbol: str) -> dict:
        """Constraint metrics for a gene: pLI, LOEUF, mis_z, syn_z, lof_z."""
        key = self._cache_key("constraint", gene_symbol)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        data = self._graphql(_constraint_query(gene_symbol), {})
        gene = data.get("gene")
        if not gene:
            return {"gene_symbol": gene_symbol, "error": "not_found"}

        constraint = gene.get("gnomad_constraint") or {}
        result = {
            "gene_symbol": gene.get("symbol"),
            "ensembl_id": gene.get("gene_id"),
            "pLI": constraint.get("pLI"),
            "oe_lof": constraint.get("oe_lof"),
            "oe_lof_upper": constraint.get("oe_lof_upper"),
            "oe_mis": constraint.get("oe_mis"),
            "oe_mis_upper": constraint.get("oe_mis_upper"),
            "mis_z": constraint.get("mis_z"),
            "syn_z": constraint.get("syn_z"),
            "lof_z": constraint.get("lof_z"),
        }
        self.cache.set(key, result)
        return result

    def validate(self) -> bool:
        try:
            data = self._graphql("{ meta { clinvar_release_date } }", None)
            return "meta" in data
        except Exception:
            return False
