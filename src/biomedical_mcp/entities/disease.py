"""Disease composite entity.

Each entities/<entity>.py exports this exact interface so domains/composite.py can
register composites generically (one tool per module, zero per-entity wiring):

    ENTITY       str    — entity key ("disease")
    TOOL_NAME    str    — MCP tool name ("disease_profile")
    DESCRIPTION  str    — tool docstring shown to the agent
    NEEDS        tuple  — client keys this entity's fetchers require
    SECTIONS     tuple[Section, ...]  — static section metadata (describe_sections)
    resolve(identifier, clients) -> (resolved_id, short_circuit | None)
    fetchers(identifier, clients, **opts) -> {section_name: zero-arg callable}

ID-RESOLUTION STRATEGY (the core complexity):
============================================================
Disease identifiers come in four forms:

  1. Plain name: "cystic fibrosis" — no canonical structure.
  2. ORPHA:<int> — Orphanet numeric code, e.g. ORPHA:586.
  3. MONDO:<int> — MONDO ontology CURIE, e.g. MONDO:0009861.
  4. OMIM:<int> — OMIM disease code, e.g. OMIM:219700.
  5. EFO:<str>  — EBI EFO code, e.g. EFO:0000270. (OpenTargets native)

resolve() does lightweight prefix detection and returns (identifier, None) in
all cases — no expensive network round-trip. The per-section fetchers handle
the rest:

  Orphanet sections (natural_history, epidemiology, disease_genes):
    • ORPHA:<int>  → strip prefix, pass int ORPHAcode directly ✓
    • OMIM:<int>   → call orphanet.search_by_omim() to find ORPHAcodes ✓
    • MONDO/EFO    → cannot cheaply resolve to ORPHAcode; skip with error-dict
    • Plain name   → orphanet has no free-text search endpoint; skip with error-dict

  Monarch (phenotypes):
    • MONDO / OMIM / ORPHA CURIE form accepted by Monarch as disease_id ✓
    • Plain name → fallback to monarch.search() to find best CURIE, then disease_phenotypes
    • ORPHA:<int> → format as "ORPHA:<int>" CURIE before calling Monarch

  OpenTargets (associations):
    • MONDO / EFO  → pass directly as disease_id ✓
    • OMIM          → pass as "OMIM:<int>" — OT accepts this ✓
    • ORPHA:<int>  → pass as "Orphanet_<int>" (OT EFO mapping convention) ✓
    • Plain name   → resolve via OT search(), take top hit id

  ICD-10 (coding):
    • Always by name search — works for any input form ✓
    • For plain names, use the identifier directly
    • For IDs, use the original identifier string as a fallback name

LIMITATIONS:
  - Orphanet sections silently error for MONDO/EFO/plain-name inputs — these
    identifiers cannot be cheaply mapped to ORPHAcodes without a lookup that
    requires the Orphanet API (which has no CURIE→ORPHAcode search endpoint).
    Workaround: pass ORPHA: or OMIM: prefix.
  - Monarch search-based fallback for plain names returns the first hit, which
    may not be the intended disease if the name is ambiguous.
  - OpenTargets requires an EFO/MONDO-style id. For ORPHA inputs we use the
    Orphanet_<int> convention but OT may not have every rare disease.
  - ICD-10 search is name-based; structured lookup (icd10.lookup) is only run
    when the input looks like an ICD-10 code (letter + digits).
"""

from __future__ import annotations

import re
from typing import Any

from biomedical_mcp.composite_core import Section, Fetcher

ENTITY = "disease"
TOOL_NAME = "disease_profile"
DESCRIPTION = (
    "Comprehensive disease report via a standardized partial-failure envelope. "
    "Sections (default: all): associations, phenotypes, natural_history, "
    "epidemiology, disease_genes, coding. "
    "Call describe_sections('disease') for sources/descriptions. "
    "One dead upstream API degrades only its section; the profile still returns. "
    "Pass `identifier` = disease name (e.g. 'cystic fibrosis') or a structured id "
    "(ORPHA:586, MONDO:0009861, OMIM:219700, EFO:0000270). "
    "Orphanet sections (natural_history, epidemiology, disease_genes) work best "
    "with ORPHA: or OMIM: prefix. Monarch/OpenTargets work with any form."
)

NEEDS = ("orphanet", "monarch", "opentargets", "icd10")

SECTIONS: tuple[Section, ...] = (
    Section(
        "associations",
        ("opentargets",),
        True,
        "Target-disease associations from Open Targets (gene targets, scores, evidence sources). "
        "Requires a resolvable disease id (MONDO/EFO/OMIM/ORPHA); plain names use OT search.",
    ),
    Section(
        "phenotypes",
        ("monarch",),
        True,
        "HPO phenotype profile from Monarch Initiative. Accepts MONDO/OMIM/ORPHA CURIEs "
        "or falls back to Monarch search for plain names.",
    ),
    Section(
        "natural_history",
        ("orphanet",),
        True,
        "Inheritance patterns, age of onset, age of death from Orphanet. "
        "Requires ORPHA: or OMIM: input; returns error-shaped dict for other id forms.",
        "B2",
    ),
    Section(
        "epidemiology",
        ("orphanet",),
        True,
        "Prevalence and incidence data from Orphanet. "
        "Requires ORPHA: or OMIM: input; returns error-shaped dict for other id forms.",
        "B2",
    ),
    Section(
        "disease_genes",
        ("orphanet",),
        True,
        "Genes associated with a rare disease from Orphanet. "
        "Requires ORPHA: or OMIM: input; returns error-shaped dict for other id forms.",
        "B2",
    ),
    Section(
        "coding",
        ("icd10",),
        True,
        "ICD-10-CM codes via NLM search. Works for all input forms (name-based search).",
    ),
)


# ── ID form detection helpers ────────────────────────────────────────────────

_ORPHA_RE = re.compile(r"^ORPHA:(\d+)$", re.IGNORECASE)
_MONDO_RE = re.compile(r"^MONDO:\d+$", re.IGNORECASE)
_OMIM_RE = re.compile(r"^OMIM:(\d+)$", re.IGNORECASE)
_EFO_RE = re.compile(r"^EFO:[A-Z0-9_]+$", re.IGNORECASE)
_ICD10_RE = re.compile(r"^[A-Z]\d{2}(\.\d+)?$")  # e.g. E84, E84.0


def _parse_orpha_code(identifier: str) -> int | None:
    """Return integer ORPHAcode if identifier is ORPHA:<int>, else None."""
    m = _ORPHA_RE.match(identifier)
    return int(m.group(1)) if m else None


def _parse_omim_code(identifier: str) -> str | None:
    """Return bare OMIM code string if identifier is OMIM:<int>, else None."""
    m = _OMIM_RE.match(identifier)
    return m.group(1) if m else None


def _is_curie(identifier: str) -> bool:
    """True if identifier is a structured CURIE (MONDO/OMIM/ORPHA/EFO)."""
    return bool(
        _ORPHA_RE.match(identifier)
        or _MONDO_RE.match(identifier)
        or _OMIM_RE.match(identifier)
        or _EFO_RE.match(identifier)
    )


def _opentargets_id(identifier: str) -> str | None:
    """Derive an OT-compatible disease id from a structured identifier.

    Returns:
        - MONDO/EFO as-is (OT accepts these natively).
        - OMIM:XXXXXX as "OMIM_XXXXXX" (OT EFO convention).
        - ORPHA:123   as "Orphanet_123" (OT EFO convention).
        - Plain name  → None (caller must use OT search).
    """
    if _MONDO_RE.match(identifier) or _EFO_RE.match(identifier):
        return identifier
    m_omim = _OMIM_RE.match(identifier)
    if m_omim:
        return f"OMIM_{m_omim.group(1)}"
    m_orpha = _ORPHA_RE.match(identifier)
    if m_orpha:
        return f"Orphanet_{m_orpha.group(1)}"
    return None


def _monarch_id(identifier: str) -> str | None:
    """Return a Monarch-compatible disease CURIE, or None if plain name."""
    # MONDO, OMIM, EFO accepted directly by Monarch
    if _MONDO_RE.match(identifier) or _OMIM_RE.match(identifier) or _EFO_RE.match(identifier):
        return identifier
    # ORPHA:123 → Monarch expects "ORPHA:123" (same format)
    if _ORPHA_RE.match(identifier):
        return identifier
    return None


def _search_name_for_display(identifier: str) -> str:
    """Extract a human-readable name from any identifier form for ICD-10/search.

    For CURIEs, strip the prefix for a better search term. For plain names, use as-is.
    """
    for prefix in ("ORPHA:", "MONDO:", "OMIM:", "EFO:"):
        if identifier.upper().startswith(prefix):
            # Return the bare code — ICD-10 search works better with a disease name,
            # so callers should prefer the identifier as-is for CURIE inputs.
            return identifier
    return identifier


# ── resolve ──────────────────────────────────────────────────────────────────

def resolve(identifier: str, clients: dict[str, Any]) -> tuple[str, dict | None]:
    """Detect and validate identifier form. Returns (identifier, None).

    No network call here — ID form detection is pure regex. Each fetcher
    handles its own resolution needs (see module docstring strategy).
    A non-None second element would short-circuit run_sections; we never
    do this so every section always gets a chance to run (partial is acceptable).
    """
    return identifier, None


# ── fetchers ─────────────────────────────────────────────────────────────────

def fetchers(identifier: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:  # noqa: C901
    """Bind per-section callables. Each fetcher handles its own ID resolution."""
    limit = opts.get("limit", 25)
    c = clients

    # ── associations (OpenTargets) ────────────────────────────────────────────
    def _associations() -> dict:
        ot_id = _opentargets_id(identifier)
        if ot_id is None:
            # Plain name: resolve via OT search
            hits = c["opentargets"].search(identifier, entity_type="disease", limit=5)
            if not hits:
                return {"error": f"disease '{identifier}' not found in OpenTargets search"}
            ot_id = hits[0]["id"]
        return c["opentargets"].disease_targets(ot_id, limit=limit)

    # ── phenotypes (Monarch) ──────────────────────────────────────────────────
    def _phenotypes() -> dict:
        curie = _monarch_id(identifier)
        if curie is None:
            # Plain name: search Monarch for top hit, then fetch phenotypes
            search_result = c["monarch"].search(identifier, category="biolink:Disease", limit=3)
            results = search_result.get("results", [])
            if not results:
                return {"error": f"disease '{identifier}' not found in Monarch search"}
            curie = results[0]["id"]
        return c["monarch"].disease_phenotypes(curie, limit=limit)

    # ── natural_history (Orphanet, requires ORPHAcode) ────────────────────────
    def _natural_history() -> dict:
        orpha_code = _parse_orpha_code(identifier)
        if orpha_code is not None:
            return c["orphanet"].natural_history(orpha_code)
        omim_code = _parse_omim_code(identifier)
        if omim_code is not None:
            # Resolve OMIM → ORPHAcode via cross-reference lookup
            xref = c["orphanet"].search_by_omim(omim_code)
            diseases = xref.get("diseases", [])
            if not diseases or diseases[0].get("orphacode") is None:
                return {"error": f"OMIM:{omim_code} not found in Orphanet cross-reference"}
            return c["orphanet"].natural_history(diseases[0]["orphacode"])
        # MONDO, EFO, or plain name — no cheap mapping available
        return {
            "error": (
                f"Cannot resolve '{identifier}' to an Orphanet ORPHAcode. "
                "Use ORPHA:<int> or OMIM:<int> for Orphanet sections."
            )
        }

    # ── epidemiology (Orphanet, requires ORPHAcode) ───────────────────────────
    def _epidemiology() -> dict:
        orpha_code = _parse_orpha_code(identifier)
        if orpha_code is not None:
            return c["orphanet"].epidemiology(orpha_code)
        omim_code = _parse_omim_code(identifier)
        if omim_code is not None:
            xref = c["orphanet"].search_by_omim(omim_code)
            diseases = xref.get("diseases", [])
            if not diseases or diseases[0].get("orphacode") is None:
                return {"error": f"OMIM:{omim_code} not found in Orphanet cross-reference"}
            return c["orphanet"].epidemiology(diseases[0]["orphacode"])
        return {
            "error": (
                f"Cannot resolve '{identifier}' to an Orphanet ORPHAcode. "
                "Use ORPHA:<int> or OMIM:<int> for Orphanet sections."
            )
        }

    # ── disease_genes (Orphanet, requires ORPHAcode) ─────────────────────────
    def _disease_genes() -> dict:
        orpha_code = _parse_orpha_code(identifier)
        if orpha_code is not None:
            return c["orphanet"].disease_genes(orpha_code)
        omim_code = _parse_omim_code(identifier)
        if omim_code is not None:
            xref = c["orphanet"].search_by_omim(omim_code)
            diseases = xref.get("diseases", [])
            if not diseases or diseases[0].get("orphacode") is None:
                return {"error": f"OMIM:{omim_code} not found in Orphanet cross-reference"}
            return c["orphanet"].disease_genes(diseases[0]["orphacode"])
        return {
            "error": (
                f"Cannot resolve '{identifier}' to an Orphanet ORPHAcode. "
                "Use ORPHA:<int> or OMIM:<int> for Orphanet sections."
            )
        }

    # ── coding (ICD-10) ───────────────────────────────────────────────────────
    def _coding() -> dict:
        # If identifier looks like an ICD-10 code already, try direct lookup
        if _ICD10_RE.match(identifier):
            result = c["icd10"].lookup(identifier)
            if result:
                return {"input": identifier, "matched": result, "search_results": []}
        # Always run search — works for names and CURIE-style inputs
        search_term = identifier  # ICD-10 search handles partial names reasonably
        results = c["icd10"].search(search_term, limit=10)
        return {
            "input": identifier,
            "count": len(results),
            "codes": results,
        }

    return {
        "associations": _associations,
        "phenotypes": _phenotypes,
        "natural_history": _natural_history,
        "epidemiology": _epidemiology,
        "disease_genes": _disease_genes,
        "coding": _coding,
    }
