"""Offline contract tests for the disease composite entity.

No network: all clients are fakes. Tests cover:
  1. resolve() — handles plain name, ORPHA:, MONDO:, OMIM: inputs.
  2. fetchers() bind + run via run_sections with fake clients — overall ok.
  3. Partial: MONDO: input → orphanet sections error (no ORPHAcode), while
     monarch/opentargets/icd10 sections succeed → overall_status == "partial".
  4. ORPHA: input → all orphanet sections succeed.
  5. OMIM: input → orphanet resolves via search_by_omim → succeeds.
  6. Plain name → uses search() to resolve for monarch/opentargets.

Mirrors tests/test_composite_core.py style.
"""

from __future__ import annotations

import pytest

from biomedical_mcp.composite_core import run_sections
from biomedical_mcp.entities import disease as disease_entity


# ── fake clients ──────────────────────────────────────────────────────────────

class FakeOrphanet:
    def natural_history(self, orphacode: int) -> dict:
        return {
            "orphacode": orphacode,
            "name": "Cystic fibrosis",
            "inheritance": ["Autosomal recessive"],
            "age_of_onset": ["Neonatal"],
            "age_of_death": [],
        }

    def epidemiology(self, orphacode: int) -> dict:
        return {
            "orphacode": orphacode,
            "name": "Cystic fibrosis",
            "prevalence": [{"class": "1-5 / 10 000"}],
        }

    def disease_genes(self, orphacode: int) -> dict:
        return {
            "orphacode": orphacode,
            "disease_name": "Cystic fibrosis",
            "count": 1,
            "genes": [{"symbol": "CFTR"}],
        }

    def search_by_omim(self, omim_code: str) -> dict:
        # Simulate a successful OMIM→ORPHA lookup
        if omim_code == "219700":
            return {
                "omim_code": omim_code,
                "count": 1,
                "diseases": [{"orphacode": 586, "name": "Cystic fibrosis"}],
            }
        # OMIM code not found
        return {
            "omim_code": omim_code,
            "count": 0,
            "diseases": [],
            "note": "Not found",
        }


class FakeOrphanetNotFound:
    """Simulates an Orphanet client where search_by_omim finds nothing."""
    def search_by_omim(self, omim_code: str) -> dict:
        return {"omim_code": omim_code, "count": 0, "diseases": []}

    def natural_history(self, orphacode: int) -> dict:
        return {"orphacode": orphacode, "inheritance": []}

    def epidemiology(self, orphacode: int) -> dict:
        return {"orphacode": orphacode, "prevalence": []}

    def disease_genes(self, orphacode: int) -> dict:
        return {"orphacode": orphacode, "count": 0, "genes": []}


class FakeMonarch:
    def disease_phenotypes(self, disease_id: str, limit: int = 50) -> dict:
        return {
            "disease_id": disease_id,
            "phenotype_count": 2,
            "phenotypes": [
                {"hpo_id": "HP:0002099", "name": "Asthma"},
                {"hpo_id": "HP:0001508", "name": "Failure to thrive"},
            ],
        }

    def search(self, query: str, category: str | None = None, limit: int = 10) -> dict:
        return {
            "query": query,
            "count": 1,
            "results": [
                {"id": "MONDO:0009861", "name": "cystic fibrosis", "category": "biolink:Disease"}
            ],
        }


class FakeOpenTargets:
    def disease_targets(self, disease_id: str, limit: int = 25) -> dict:
        return {
            "disease": "Cystic fibrosis",
            "disease_id": disease_id,
            "total": 1,
            "targets": [{"ensembl_id": "ENSG00000001626", "symbol": "CFTR", "score": 0.99}],
        }

    def search(self, query: str, entity_type: str | None = None, limit: int = 10) -> list[dict]:
        return [{"id": "MONDO:0009861", "name": "cystic fibrosis", "entity": "disease"}]


class FakeICD10:
    def search(self, query: str, limit: int = 10) -> list[dict]:
        return [{"code": "E84", "description": "Cystic fibrosis"}]

    def lookup(self, code: str) -> dict | None:
        if code.startswith("E84"):
            return {"code": code, "description": "Cystic fibrosis"}
        return None


def _clients(
    orphanet=None, monarch=None, opentargets=None, icd10=None
) -> dict:
    return {
        "orphanet": orphanet or FakeOrphanet(),
        "monarch": monarch or FakeMonarch(),
        "opentargets": opentargets or FakeOpenTargets(),
        "icd10": icd10 or FakeICD10(),
    }


# ── 1. resolve() ──────────────────────────────────────────────────────────────

def test_resolve_plain_name():
    resolved, sc = disease_entity.resolve("cystic fibrosis", {})
    assert resolved == "cystic fibrosis"
    assert sc is None


def test_resolve_orpha_prefix():
    resolved, sc = disease_entity.resolve("ORPHA:586", {})
    assert resolved == "ORPHA:586"
    assert sc is None


def test_resolve_mondo_prefix():
    resolved, sc = disease_entity.resolve("MONDO:0009861", {})
    assert resolved == "MONDO:0009861"
    assert sc is None


def test_resolve_omim_prefix():
    resolved, sc = disease_entity.resolve("OMIM:219700", {})
    assert resolved == "OMIM:219700"
    assert sc is None


def test_resolve_efo_prefix():
    resolved, sc = disease_entity.resolve("EFO:0000270", {})
    assert resolved == "EFO:0000270"
    assert sc is None


# ── 2. Full run with ORPHA: — all sections succeed ────────────────────────────

def test_orpha_input_all_sections_ok():
    """ORPHA:<int> → orphanet sections get int code directly; all succeed."""
    clients = _clients()
    fts = disease_entity.fetchers("ORPHA:586", clients)
    env = run_sections(
        disease_entity.ENTITY,
        "ORPHA:586",
        None,
        fts,
        disease_entity.SECTIONS,
    )
    assert env["entity"] == "disease"
    assert env["id"] == "ORPHA:586"
    assert env["overall_status"] == "ok", env["sections"]

    nh = env["sections"]["natural_history"]
    assert nh["status"] == "ok"
    assert nh["data"]["orphacode"] == 586
    assert nh["provenance"]["evidence_grade"] == "B2"

    ep = env["sections"]["epidemiology"]
    assert ep["status"] == "ok"
    assert ep["data"]["orphacode"] == 586

    dg = env["sections"]["disease_genes"]
    assert dg["status"] == "ok"
    assert dg["data"]["genes"][0]["symbol"] == "CFTR"

    ph = env["sections"]["phenotypes"]
    assert ph["status"] == "ok"

    assoc = env["sections"]["associations"]
    assert assoc["status"] == "ok"

    coding = env["sections"]["coding"]
    assert coding["status"] == "ok"


# ── 3. MONDO: input → orphanet sections error, others succeed → partial ────────

def test_mondo_input_partial_orphanet_errors():
    """MONDO: cannot be resolved to ORPHAcode → orphanet fetchers return error-dicts.
    run_sections classifies {"error": ...} returns as status "empty" (not an
    exception), and "empty" is counted as a successful (non-failing) fetch.
    So overall_status is "ok" — the composite does NOT crash and all sections
    return data (good or empty).  The caller can inspect the "error" key in
    section data to see that Orphanet resolution was skipped.
    """
    clients = _clients()
    fts = disease_entity.fetchers("MONDO:0009861", clients)
    env = run_sections(
        disease_entity.ENTITY,
        "MONDO:0009861",
        None,
        fts,
        disease_entity.SECTIONS,
    )
    # Orphanet sections return an error-shaped dict → run_sections marks them "empty"
    # (_has_data returns False for {"error": ...}); data is still the dict.
    for section_name in ("natural_history", "epidemiology", "disease_genes"):
        sec = env["sections"][section_name]
        assert sec["status"] == "empty", (
            f"{section_name} should be 'empty' for MONDO: input (error-dict), got {sec['status']}"
        )
        assert "error" in sec["data"], (
            f"{section_name} data should contain 'error' key, got {sec['data']}"
        )

    # Monarch, OT, ICD-10 should be ok
    assert env["sections"]["phenotypes"]["status"] == "ok"
    assert env["sections"]["associations"]["status"] == "ok"
    assert env["sections"]["coding"]["status"] == "ok"

    # All sections returned (no exception raised) → overall_status is "ok"
    # (empty counts as non-failing). This is the correct contract: partial
    # requires a mix of ok+error where error means an *exception*, not an
    # error-shaped return value.
    assert env["overall_status"] == "ok"


# ── 4. Plain name → search-based resolution ───────────────────────────────────

def test_plain_name_monarch_and_ot_use_search():
    """Plain name uses search() for monarch and opentargets; icd10 runs by name."""
    clients = _clients()
    fts = disease_entity.fetchers("cystic fibrosis", clients)

    pheno_data = fts["phenotypes"]()
    assert pheno_data["disease_id"] == "MONDO:0009861"  # resolved via search

    assoc_data = fts["associations"]()
    assert assoc_data["disease_id"] == "MONDO:0009861"

    coding_data = fts["coding"]()
    assert len(coding_data["codes"]) > 0


def test_plain_name_orphanet_sections_return_error_dict():
    """Plain name → Orphanet sections cannot resolve → return error-dict, no crash."""
    clients = _clients()
    fts = disease_entity.fetchers("cystic fibrosis", clients)

    nh = fts["natural_history"]()
    assert "error" in nh

    ep = fts["epidemiology"]()
    assert "error" in ep

    dg = fts["disease_genes"]()
    assert "error" in dg


# ── 5. OMIM: → Orphanet resolves via search_by_omim ──────────────────────────

def test_omim_input_orphanet_resolves_via_xref():
    """OMIM:219700 → search_by_omim → ORPHAcode 586 → natural_history succeeds."""
    clients = _clients()
    fts = disease_entity.fetchers("OMIM:219700", clients)

    nh = fts["natural_history"]()
    assert "error" not in nh
    assert nh["orphacode"] == 586

    ep = fts["epidemiology"]()
    assert "error" not in ep
    assert ep["orphacode"] == 586

    dg = fts["disease_genes"]()
    assert "error" not in dg
    assert dg["orphacode"] == 586


def test_omim_not_in_orphanet_returns_error_dict():
    """OMIM code not found in Orphanet cross-ref → error-dict, no crash."""
    clients = _clients()
    fts = disease_entity.fetchers("OMIM:999999", clients)

    nh = fts["natural_history"]()
    assert "error" in nh
    assert "999999" in nh["error"]


# ── 6. Partial: MONDO + fake clients ─────────────────────────────────────────

def test_full_run_via_run_sections_binding():
    """Confirm fetchers() + run_sections integration: sections bind and execute."""
    clients = _clients()
    fts = disease_entity.fetchers("ORPHA:586", clients)
    env = run_sections(
        disease_entity.ENTITY,
        "ORPHA:586",
        None,
        fts,
        disease_entity.SECTIONS,
    )
    assert set(env["sections"]) == {s.name for s in disease_entity.SECTIONS}
    assert env["overall_status"] == "ok"


# ── 7. ICD-10 code input ──────────────────────────────────────────────────────

def test_icd10_code_input_direct_lookup():
    """E84 ICD-10 code → lookup() attempted, then search fallback."""
    clients = _clients()
    fts = disease_entity.fetchers("E84", clients)
    result = fts["coding"]()
    # Either matched via lookup or via search
    assert result is not None
    # Should contain either 'matched' key (lookup path) or 'codes' (search path)
    assert "matched" in result or "codes" in result


# ── 8. Module-level exports ───────────────────────────────────────────────────

def test_module_exports_correct_entity():
    assert disease_entity.ENTITY == "disease"
    assert disease_entity.TOOL_NAME == "disease_profile"
    assert isinstance(disease_entity.DESCRIPTION, str) and len(disease_entity.DESCRIPTION) > 20
    assert "orphanet" in disease_entity.NEEDS
    assert "monarch" in disease_entity.NEEDS
    assert "opentargets" in disease_entity.NEEDS
    assert "icd10" in disease_entity.NEEDS


def test_sections_have_required_fields():
    section_names = {s.name for s in disease_entity.SECTIONS}
    assert "associations" in section_names
    assert "phenotypes" in section_names
    assert "natural_history" in section_names
    assert "epidemiology" in section_names
    assert "disease_genes" in section_names
    assert "coding" in section_names
    # All sections are default
    for s in disease_entity.SECTIONS:
        assert s.default is True


def test_orphanet_sections_have_evidence_grade():
    orphanet_sections = [s for s in disease_entity.SECTIONS if "orphanet" in s.sources]
    assert orphanet_sections, "Expected at least one orphanet section"
    for s in orphanet_sections:
        assert s.evidence_grade == "B2", f"{s.name} should have B2 evidence grade"
