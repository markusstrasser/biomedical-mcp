"""Provenance envelope for tool outputs — tracks source, version, evidence grade."""

from __future__ import annotations

from datetime import datetime, timezone


def make_provenance(
    source_api: str,
    *,
    data_version: str | None = None,
    evidence_grade: str | None = None,
    cache_hit: bool = False,
) -> dict:
    """Build a provenance dict for tool responses.

    Evidence grades follow the constitutional hierarchy:
      A1 = systematic review, B2 = functional assay, C3 = ClinVar 3+ stars,
      C4 = ClinVar 1 star, D4 = case report, E5 = in silico / population freq,
      F6 = LLM-generated
    """
    return {
        "source_api": source_api,
        "data_version": data_version,
        "evidence_grade": evidence_grade,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "cache_hit": cache_hit,
    }


def envelope(data: dict, provenance: dict) -> dict:
    """Wrap tool output in a provenance envelope."""
    return {"data": data, "provenance": provenance}
