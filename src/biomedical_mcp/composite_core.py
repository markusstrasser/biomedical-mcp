"""Composite core — the standardized partial-failure envelope + section runner.

This is the CONTRACT every entity composite conforms to (gene_dossier,
variant_context, drug_profile, protein_profile, disease_profile). It exists so
the consolidated tool surface (a handful of entity composites instead of ~85 raw
tools) degrades gracefully: one dead upstream API fails ONE section, never the
whole dossier, and never returns a silent null that the model hallucinates over.

Design notes (decision: agent-infra/decisions/2026-06-07-biomedical-mcp-tool-consolidation.md):
- `Section` is pure static metadata — no client references — so describe_sections()
  can enumerate sections with zero per-call cost.
- `run_sections()` runs the requested (or default) sections concurrently with a
  per-section timeout, wrapping each in a SectionResult with explicit status.
- The 32 client classes ARE the adapter layer; composites call them via fetchers.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeout
from dataclasses import dataclass
from typing import Any, Callable

# Shared thread pool — clients use sync httpx, so sections run on threads.
_executor = ThreadPoolExecutor(max_workers=8)

# Per-section deadline. A slow/down API yields a timeout SectionResult instead
# of blocking the whole composite. Pattern inherited from the prior composite.py.
SECTION_TIMEOUT = 8


@dataclass(frozen=True)
class Section:
    """Static metadata for one composite section. Drives `describe_sections`.

    No client references — pure data, defined once in composite_sections.SECTIONS.
    """

    name: str
    sources: tuple[str, ...]          # source_api keys, e.g. ("gnomad",)
    default: bool                     # included when the caller omits `sections`
    description: str
    evidence_grade: str | None = None  # constitutional grade (A1..F6), if applicable
    params: tuple[str, ...] = ()       # optional extra params this section honors


# A fetcher returns the section's raw data (any JSON-able value). It may raise;
# run_sections() catches and classifies the failure.
Fetcher = Callable[[], Any]


def _has_data(data: Any) -> bool:
    """Did the fetcher return 'found' data vs empty/None/error-shaped?"""
    if data is None:
        return False
    if isinstance(data, dict) and data.get("error"):
        return False
    if isinstance(data, (list, dict, str)) and len(data) == 0:
        return False
    return True


def run_sections(
    entity: str,
    identifier: str,
    requested: list[str] | None,
    fetchers: dict[str, Fetcher],
    sections_meta: tuple[Section, ...],
    *,
    timeout: int = SECTION_TIMEOUT,
    extra: dict | None = None,
) -> dict:
    """Run sections concurrently and assemble the standardized envelope.

    Args:
        entity: entity type ("gene", "variant", "drug", ...).
        identifier: the resolved id/symbol the sections were fetched for.
        requested: explicit section names, or None to run all `default=True`.
        fetchers: section name -> zero-arg callable returning that section's data.
        sections_meta: the entity's Section tuple (for metadata + validation).
        timeout: per-section deadline in seconds.
        extra: optional top-level keys to merge into the result (e.g. normalized id).

    Returns an envelope:
        {
          "entity", "id", "overall_status": "ok|partial|error",
          "sections": { name: {"status", "data", "provenance"?, "error"?} }
        }
    Section status is one of: ok, empty, error, not_applicable.
    """
    meta_by_name = {s.name: s for s in sections_meta}
    if requested is None:
        names = [s.name for s in sections_meta if s.default]
    else:
        names = [n.strip() for n in requested if n and n.strip()]

    sections_out: dict[str, dict] = {}
    runnable: list[str] = []

    for n in names:
        if n not in meta_by_name:
            sections_out[n] = {
                "status": "not_applicable",
                "data": None,
                "error": {
                    "type": "unknown_section",
                    "retryable": False,
                    "message": f"'{n}' is not a section of '{entity}'. "
                               f"Call describe_sections('{entity}') for valid sections.",
                },
            }
        elif n not in fetchers:
            sections_out[n] = {
                "status": "not_applicable",
                "data": None,
                "error": {
                    "type": "unavailable",
                    "retryable": False,
                    "message": f"section '{n}' has no fetcher bound for this request",
                },
            }
        else:
            runnable.append(n)

    futures = {n: _executor.submit(fetchers[n]) for n in runnable}
    for n in runnable:
        meta = meta_by_name[n]
        prov = {"sources": list(meta.sources), "evidence_grade": meta.evidence_grade}
        try:
            data = futures[n].result(timeout=timeout)
            status = "ok" if _has_data(data) else "empty"
            sections_out[n] = {"status": status, "data": data, "provenance": prov}
        except FuturesTimeout:
            sections_out[n] = {
                "status": "error", "data": None, "provenance": prov,
                "error": {"type": "timeout", "retryable": True,
                          "message": f"{n} timed out after {timeout}s"},
            }
        except Exception as exc:  # noqa: BLE001 — surface any upstream failure per-section
            sections_out[n] = {
                "status": "error", "data": None, "provenance": prov,
                "error": {"type": type(exc).__name__, "retryable": False,
                          "message": str(exc)},
            }

    statuses = {s["status"] for s in sections_out.values()}
    has_good = bool(statuses & {"ok", "empty"})
    has_bad = bool(statuses & {"error", "not_applicable"})
    if has_good and not has_bad:
        overall = "ok"
    elif has_good and has_bad:
        overall = "partial"
    else:
        overall = "error"

    out = {
        "entity": entity,
        "id": identifier,
        "overall_status": overall,
        "sections": sections_out,
    }
    if extra:
        out.update(extra)
    return out


def describe(entity: str, sections_meta: tuple[Section, ...]) -> dict:
    """Serialize a entity's sections for the describe_sections meta-tool."""
    return {
        "entity": entity,
        "sections": [
            {
                "name": s.name,
                "sources": list(s.sources),
                "default": s.default,
                "description": s.description,
                "evidence_grade": s.evidence_grade,
                "params": list(s.params),
            }
            for s in sections_meta
        ],
    }
