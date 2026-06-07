"""Somatic variant annotation — CIViC clinical evidence + OncoKB oncogenicity.

CIViC: GraphQL API at https://civicdb.org/api/graphql — open, no auth.
       Queries evidenceItems filtered by molecularProfileName (e.g. "BRAF V600E").

OncoKB: REST API at https://www.oncokb.org/api/v1 — requires ONCOKB_TOKEN env var
        (Bearer auth). Returns oncogenicity, mutation effect, treatment levels.
        Missing token → returns ApiKeyRequiredError.to_dict(), never raises.

Usage:
    from biomedical_mcp.somatic import Somatic
    client = Somatic(cache)
    result = client.annotate("BRAF", "V600E")
    result = client.annotate("BRAF V600E")   # combined form auto-parsed
"""

from __future__ import annotations

import hashlib
import logging
import os

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache
from biomedical_mcp.errors import ApiKeyRequiredError

log = logging.getLogger(__name__)

# Base URLs are hardcoded — these are public/semi-public APIs with stable URLs
# that don't share a config-domain (CIViC is oncology-specific, OncoKB requires
# auth). We register a single "somatic" domain in config.py for TTL + rate
# defaults; both hostnames are encoded here, matching how OpenTargets does it.
CIVIC_URL = "https://civicdb.org/api/graphql"
ONCOKB_URL = "https://www.oncokb.org/api/v1"

ONCOKB_ENV_VAR = "ONCOKB_TOKEN"
ONCOKB_DOCS = "https://www.oncokb.org/apiAccess"

# CIViC GraphQL query — filter by molecularProfileName ("BRAF V600E").
# Verified fields via live introspection 2026-06-07:
#   EvidenceItem: id, evidenceType, evidenceLevel, evidenceDirection,
#                 significance, status, description, disease, therapies, source
#   Disease: name
#   Therapy: name
#   Source: citation, citationId (citationId = PMID for PubMed sources)
_CIVIC_QUERY = """
query($name: String!, $first: Int!, $status: EvidenceStatusFilter) {
  evidenceItems(
    molecularProfileName: $name
    status: $status
    first: $first
  ) {
    totalCount
    nodes {
      id
      evidenceType
      evidenceLevel
      evidenceDirection
      significance
      status
      description
      variantOrigin
      disease { name }
      therapies { name }
      source {
        citation
        citationId
        sourceType
        link
      }
    }
  }
}
"""


def _parse_combined(gene_or_combined: str, protein_change: str | None) -> tuple[str, str]:
    """Parse gene + protein_change from either separate args or combined 'BRAF V600E'."""
    if protein_change is not None:
        return gene_or_combined.strip(), protein_change.strip()
    parts = gene_or_combined.strip().split(None, 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    raise ValueError(
        f"Cannot parse gene+change from '{gene_or_combined}'. "
        "Pass gene and protein_change separately or as 'GENE CHANGE' (e.g. 'BRAF V600E')."
    )


def _parse_civic_evidence(nodes: list[dict]) -> list[dict]:
    """Flatten CIViC evidenceItems nodes into a clean list."""
    out = []
    for node in nodes:
        disease = node.get("disease") or {}
        source = node.get("source") or {}
        therapies = [t.get("name") for t in (node.get("therapies") or []) if t.get("name")]
        out.append({
            "id": node.get("id"),
            "evidence_type": node.get("evidenceType"),
            "evidence_level": node.get("evidenceLevel"),
            "evidence_direction": node.get("evidenceDirection"),
            "significance": node.get("significance"),
            "variant_origin": node.get("variantOrigin"),
            "disease": disease.get("name"),
            "therapies": therapies,
            "description": node.get("description"),
            "citation": source.get("citation"),
            "pmid": source.get("citationId"),
            "source_type": source.get("sourceType"),
            "source_link": source.get("link"),
        })
    return out


def _parse_oncokb_response(raw: dict) -> dict:
    """Extract the clinically actionable fields from an OncoKB annotation response."""
    mut = raw.get("mutationEffect") or {}
    return {
        "oncogenic": raw.get("oncogenic"),
        "gene_exist": raw.get("geneExist"),
        "variant_exist": raw.get("variantExist"),
        "mutation_effect": mut.get("knownEffect"),
        "mutation_effect_description": mut.get("description"),
        "highest_sensitive_level": raw.get("highestSensitiveLevel"),
        "highest_resistance_level": raw.get("highestResistanceLevel"),
        "highest_diagnostic_level": raw.get("highestDiagnosticImplicationLevel"),
        "highest_prognostic_level": raw.get("highestPrognosticImplicationLevel"),
        "gene_summary": raw.get("geneSummary"),
        "variant_summary": raw.get("variantSummary"),
        "treatments": _parse_oncokb_treatments(raw.get("treatments") or []),
        "diagnostic_implications": _parse_oncokb_implications(raw.get("diagnosticImplications") or []),
    }


def _parse_oncokb_treatments(treatments: list[dict]) -> list[dict]:
    out = []
    for t in treatments:
        tumor = t.get("levelAssociatedCancerType") or {}
        drugs = [d.get("drugName") for d in (t.get("drugs") or []) if d.get("drugName")]
        out.append({
            "level": t.get("level"),
            "drugs": drugs,
            "cancer_type": tumor.get("name"),
            "pmids": t.get("pmids") or [],
        })
    return out


def _parse_oncokb_implications(implications: list[dict]) -> list[dict]:
    out = []
    for imp in implications:
        tumor = imp.get("tumorType") or {}
        out.append({
            "level": imp.get("levelOfEvidence"),
            "cancer_type": tumor.get("name"),
            "alterations": imp.get("alterations") or [],
        })
    return out


class Somatic(BaseClient):
    """Somatic variant annotation — CIViC (open) + OncoKB (token-gated).

    domain_name = "somatic" → reads TTL/rate from config.py DOMAINS["somatic"].
    CIViC requests go to CIVIC_URL via _post (GraphQL POST).
    OncoKB requests go to ONCOKB_URL via a separate httpx call (different base_url).
    """

    domain_name = "somatic"

    def __init__(self, cache: Cache):
        # base_url points to CIViC (used by self.client in _post).
        super().__init__(cache, base_url=CIVIC_URL)
        # Separate client for OncoKB (different base URL, may have auth header).
        import httpx
        self._oncokb_client = httpx.Client(
            base_url=ONCOKB_URL,
            timeout=30.0,
            headers={"User-Agent": "biomedical-mcp/0.5", "Accept": "application/json"},
        )

    # ── CIViC ─────────────────────────────────────────────────────────────────

    def civic_evidence(
        self,
        gene: str,
        protein_change: str,
        *,
        limit: int = 25,
        status: str = "ACCEPTED",
    ) -> list[dict]:
        """Fetch CIViC evidence items for gene + protein_change.

        Args:
            gene: Hugo gene symbol (e.g. "BRAF").
            protein_change: Single-letter protein change (e.g. "V600E").
            limit: Max evidence items to return (default 25).
            status: CIViC status filter — "ACCEPTED" (default), "SUBMITTED", or "ALL".

        Returns:
            List of flattened evidence item dicts.
        """
        mp_name = f"{gene.upper()} {protein_change}"
        cache_key = self._cache_key("civic", mp_name, str(limit), status)
        cached = self.cache.get(cache_key, max_age_days=self.ttl_days)
        if cached is not None:
            return cached

        payload = {
            "query": _CIVIC_QUERY,
            "variables": {
                "name": mp_name,
                "first": min(limit, 100),
                "status": status if status != "ALL" else None,
            },
        }
        # Remove None variables (GraphQL treats absent vs null differently here)
        payload["variables"] = {k: v for k, v in payload["variables"].items() if v is not None}

        raw = self._post("", json_data=payload)
        data = raw.get("data", {}) if isinstance(raw, dict) else {}
        errors = raw.get("errors") if isinstance(raw, dict) else None
        if errors:
            log.warning("CIViC GraphQL errors for %s: %s", mp_name, errors)

        nodes = (data.get("evidenceItems") or {}).get("nodes") or []
        result = _parse_civic_evidence(nodes)
        self.cache.set(cache_key, result)
        return result

    # ── OncoKB ────────────────────────────────────────────────────────────────

    def oncokb_annotation(
        self,
        gene: str,
        protein_change: str,
        *,
        tumor_type: str | None = None,
    ) -> dict:
        """Fetch OncoKB oncogenicity + treatment level annotation.

        Requires ONCOKB_TOKEN environment variable (Bearer token).
        Returns ApiKeyRequiredError.to_dict() if token is absent — never raises.

        Args:
            gene: Hugo gene symbol (e.g. "BRAF").
            protein_change: Protein change (e.g. "V600E").
            tumor_type: Optional tumor type for treatment-level context.

        Returns:
            Parsed annotation dict, or error dict if token is missing.
        """
        token = os.environ.get(ONCOKB_ENV_VAR)
        if not token:
            return ApiKeyRequiredError(
                api="OncoKB",
                env_var=ONCOKB_ENV_VAR,
                docs_url=ONCOKB_DOCS,
            ).to_dict()

        cache_key = self._cache_key("oncokb", gene.upper(), protein_change, tumor_type or "")
        cached = self.cache.get(cache_key, max_age_days=self.ttl_days)
        if cached is not None:
            return cached

        params: dict = {"hugoSymbol": gene.upper(), "alteration": protein_change}
        if tumor_type:
            params["tumorType"] = tumor_type

        self._rate_wait()
        resp = self._oncokb_client.get(
            "/annotate/mutations/byProteinChange",
            params=params,
            headers={"Authorization": f"Bearer {token}"},
        )
        resp.raise_for_status()
        raw = resp.json()
        result = _parse_oncokb_response(raw)
        self.cache.set(cache_key, result)
        return result

    # ── Combined annotate ─────────────────────────────────────────────────────

    def annotate(
        self,
        gene_or_combined: str,
        protein_change: str | None = None,
        *,
        tumor_type: str | None = None,
        civic_limit: int = 25,
        civic_status: str = "ACCEPTED",
    ) -> dict:
        """Annotate a somatic variant with CIViC evidence and OncoKB classification.

        Accepts either:
          annotate("BRAF", "V600E")
          annotate("BRAF V600E")

        CIViC always runs (open API). OncoKB runs if ONCOKB_TOKEN is set;
        otherwise the oncokb key contains an error dict, not an exception.

        Returns:
            {
              "gene": "BRAF",
              "protein_change": "V600E",
              "civic": [... evidence items ...],
              "civic_total": <int>,
              "oncokb": {... annotation ...} | {"error": ..., "type": "ApiKeyRequiredError", ...}
            }
        """
        gene, change = _parse_combined(gene_or_combined, protein_change)

        civic_items = self.civic_evidence(
            gene, change, limit=civic_limit, status=civic_status
        )
        oncokb_result = self.oncokb_annotation(gene, change, tumor_type=tumor_type)

        return {
            "gene": gene,
            "protein_change": change,
            "civic": civic_items,
            "civic_total": len(civic_items),
            "oncokb": oncokb_result,
        }

    def validate(self) -> bool:
        """Health check — verify CIViC GraphQL endpoint is reachable."""
        try:
            raw = self._post("", json_data={"query": "{ __typename }"})
            return isinstance(raw, dict) and "data" in raw
        except Exception:
            return False
