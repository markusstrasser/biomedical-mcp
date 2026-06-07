"""Drug composite entity — ChEMBL + OpenFDA.

Conforms to the entities/<entity>.py interface established by gene.py:

    ENTITY       str    — "drug"
    TOOL_NAME    str    — "drug_profile"
    DESCRIPTION  str    — tool docstring shown to the agent
    NEEDS        tuple  — client keys this entity's fetchers require
    SECTIONS     tuple[Section, ...]  — static section metadata
    resolve(identifier, clients) -> (resolved_id, short_circuit | None)
    fetchers(identifier, clients, **opts) -> {section_name: zero-arg callable}

ChEMBL-id resolution strategy
──────────────────────────────
resolve() passes the drug name through unchanged (returns (drug_name, None)).
Each ChEMBL-backed fetcher calls chembl.resolve_compound(name=identifier)
itself.  This is cheap because ChEMBL's client caches the lookup; a cold
resolve is one HTTP call, subsequent calls hit the in-process cache.  The
alternative — resolving once in resolve() and stashing the id — would require
widening the resolve() contract (it currently returns the identifier string,
not an arbitrary payload).  Relying on client-side caching keeps the interface
clean and matches how gene.py avoids resolve() side-effects.

If resolve_compound returns None (unknown drug), each ChEMBL fetcher raises
ValueError so run_sections() classifies only those sections as "error".
OpenFDA fetchers use the raw name string and are unaffected.

Deferred (Phase 4)
──────────────────
- Pharmacogenomics: CPIC / PharmGKB (CYP2D6 phenotype → dose recommendations).
  Will be added as a new "pharmacogenomics" Section here.
- Drug-drug interactions: DDInter API.
  Will be added as a new "interactions" Section here.
"""

from __future__ import annotations

from typing import Any

from biomedical_mcp.composite_core import Section, Fetcher

ENTITY = "drug"
TOOL_NAME = "drug_profile"
DESCRIPTION = (
    "Comprehensive drug profile via a standardized partial-failure envelope. "
    "Sections (default on: compound, mechanism, label, adverse_events; "
    "default off: indications, recalls). "
    "Call describe_sections('drug') for sources/descriptions. "
    "One dead upstream API degrades only its section; the profile still returns. "
    "Pass `identifier` = generic or brand drug name (e.g. 'aspirin', 'ibuprofen', 'atorvastatin')."
)

NEEDS = ("chembl", "openfda")

SECTIONS: tuple[Section, ...] = (
    Section(
        "compound",
        ("chembl",),
        True,
        "ChEMBL compound metadata: type, max clinical phase, physicochemical "
        "properties (MW, ALogP, HBD/HBA, PSA, Ro5 violations), ATC codes, "
        "synonyms, first approval year.",
    ),
    Section(
        "mechanism",
        ("chembl",),
        True,
        "Mechanisms of action from ChEMBL: action type, target name / ChEMBL id, "
        "supporting references.",
    ),
    Section(
        "label",
        ("openfda",),
        True,
        "FDA drug labeling: boxed warnings, contraindications, dosage, "
        "drug interactions, adverse reactions, pharmacokinetics, "
        "and use in specific populations.",
    ),
    Section(
        "adverse_events",
        ("openfda",),
        True,
        "FDA FAERS adverse event reports: total report count + sampled events "
        "with reaction terms, seriousness flags, and co-administered drugs.",
    ),
    Section(
        "indications",
        ("chembl",),
        False,  # non-default — can be expensive / verbose
        "ChEMBL drug indications: MeSH/EFO terms and max clinical phase per indication.",
    ),
    Section(
        "recalls",
        ("openfda",),
        False,  # non-default
        "FDA drug recall enforcement reports: reason, classification (Class I/II/III), "
        "status, and product description.",
    ),
)


# ── ChEMBL-id resolution helper ──────────────────────────────────────────────

def _require_chembl_id(identifier: str, chembl_client: Any) -> str:
    """Resolve drug name → ChEMBL molecule id, or raise ValueError.

    The ChEMBL client caches resolve_compound() by query string, so calling
    this from multiple fetchers for the same identifier costs at most one
    network round-trip.
    """
    chembl_id = chembl_client.resolve_compound(name=identifier)
    if not chembl_id:
        raise ValueError(
            f"Could not resolve '{identifier}' to a ChEMBL molecule id. "
            "Check spelling or try an alternative name."
        )
    return chembl_id


# ── Public interface ──────────────────────────────────────────────────────────

def resolve(identifier: str, clients: dict[str, Any]) -> tuple[str, dict | None]:
    """Pass the drug name through unchanged.

    ChEMBL id resolution is deferred to the individual fetchers so they can
    fail independently (a resolve failure degrades only ChEMBL sections, not
    OpenFDA ones).  No short-circuit is returned here.
    """
    return identifier, None


def fetchers(identifier: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    """Return per-section zero-arg callables for identifier (drug name).

    ChEMBL fetchers call _require_chembl_id() at execution time; if resolution
    fails they raise and run_sections() marks only those sections as "error".
    OpenFDA fetchers use the name string directly and are not affected.
    """
    c = clients

    def _compound() -> dict:
        cid = _require_chembl_id(identifier, c["chembl"])
        return c["chembl"].compound(cid)

    def _mechanism() -> list[dict]:
        cid = _require_chembl_id(identifier, c["chembl"])
        return c["chembl"].mechanism(cid)

    def _indications() -> list[dict]:
        cid = _require_chembl_id(identifier, c["chembl"])
        return c["chembl"].drug_indications(cid, limit=25)

    return {
        "compound": _compound,
        "mechanism": _mechanism,
        "label": lambda: c["openfda"].drug_label(identifier),
        "adverse_events": lambda: c["openfda"].adverse_events(identifier),
        "indications": _indications,
        "recalls": lambda: c["openfda"].recalls(identifier),
    }
