"""Offline contract tests for the drug composite entity.

No network — all clients are fakes. Tests verify:
  1. Module exports (ENTITY, TOOL_NAME, DESCRIPTION, NEEDS, SECTIONS).
  2. resolve() contract.
  3. fetchers() bind correctly and all six sections run via run_sections().
  4. A failed ChEMBL resolve degrades only ChEMBL sections (compound, mechanism,
     indications) while OpenFDA sections (label, adverse_events, recalls) still
     succeed — overall_status is "partial", not "error".
  5. Default sections (compound, mechanism, label, adverse_events) run; non-default
     sections (indications, recalls) do not run when sections=None.

Style mirrors tests/test_composite_core.py.
"""

import pytest

from biomedical_mcp.composite_core import run_sections, describe
from biomedical_mcp.entities import drug as drug_entity


# ── Fake clients ─────────────────────────────────────────────────────────────

CHEMBL_ID = "CHEMBL25"

FAKE_COMPOUND = {
    "chembl_id": CHEMBL_ID,
    "name": "aspirin",
    "type": "Small molecule",
    "max_phase": 4,
    "molecular_weight": "180.16",
    "alogp": "1.31",
    "hba": 3,
    "hbd": 1,
    "psa": "63.60",
    "ro5_violations": 0,
    "atc_classifications": ["N02BA01"],
    "synonyms": ["acetylsalicylic acid"],
    "first_approval": 1899,
}

FAKE_MECHANISM = [
    {
        "mechanism": "Cyclooxygenase inhibitor",
        "action_type": "INHIBITOR",
        "target_name": "Prostaglandin G/H synthase 1",
        "target_chembl_id": "CHEMBL612545",
        "references": [{"type": "PubMed", "id": "12345"}],
    }
]

FAKE_INDICATIONS = [
    {
        "mesh_id": "D010146",
        "mesh_heading": "Pain",
        "efo_id": "EFO_0000546",
        "efo_term": "pain",
        "max_phase": 4,
        "references": [],
    }
]

FAKE_LABEL = {
    "drug": "aspirin",
    "brand_names": ["Bayer Aspirin"],
    "manufacturer": ["Bayer HealthCare LLC"],
    "route": ["oral"],
    "product_type": ["HUMAN OTC DRUG"],
    "sections": {"indications_and_usage": "For the temporary relief of headache..."},
}

FAKE_ADVERSE_EVENTS = {
    "drug": "aspirin",
    "total": 12345,
    "events": [
        {
            "serious": 1,
            "seriousness_death": None,
            "reactions": ["Gastrointestinal haemorrhage"],
            "drugs": [{"name": "ASPIRIN", "indication": None, "characterization": "1"}],
            "receive_date": "20230101",
        }
    ],
}

FAKE_RECALLS = {
    "drug": "aspirin",
    "total": 3,
    "recalls": [
        {
            "recall_number": "D-1234-2023",
            "reason": "Labeling",
            "status": "Terminated",
            "classification": "Class III",
            "product_description": "Aspirin 81mg tablets",
            "recall_initiation_date": "20230101",
            "voluntary": "Voluntary: Firm Initiated",
        }
    ],
}


class FakeChEMBL:
    """Fake ChEMBL client — resolve_compound returns the molecule_chembl_id string."""

    def resolve_compound(self, name: str | None = None, chembl_id: str | None = None) -> str | None:
        # The real resolve_compound returns the chembl_id string (e.g. "CHEMBL25"),
        # NOT a dict — confirmed from chembl.py line 60: `return molecules[0]["molecule_chembl_id"]`
        if name and name.lower() == "aspirin":
            return CHEMBL_ID
        return None  # simulate unknown drug

    def compound(self, chembl_id: str) -> dict:
        assert chembl_id == CHEMBL_ID
        return FAKE_COMPOUND

    def mechanism(self, chembl_id: str) -> list[dict]:
        assert chembl_id == CHEMBL_ID
        return FAKE_MECHANISM

    def drug_indications(self, chembl_id: str, limit: int = 25) -> list[dict]:
        assert chembl_id == CHEMBL_ID
        return FAKE_INDICATIONS


class FakeChEMBLUnresolvable:
    """Fake ChEMBL that always fails to resolve."""

    def resolve_compound(self, name: str | None = None, chembl_id: str | None = None) -> str | None:
        return None

    def compound(self, chembl_id: str) -> dict:
        raise AssertionError("compound() should not be called when resolve fails")

    def mechanism(self, chembl_id: str) -> list[dict]:
        raise AssertionError("mechanism() should not be called when resolve fails")

    def drug_indications(self, chembl_id: str, limit: int = 25) -> list[dict]:
        raise AssertionError("drug_indications() should not be called when resolve fails")


class FakeOpenFDA:
    def drug_label(self, drug_name: str, sections: list[str] | None = None) -> dict:
        return FAKE_LABEL

    def adverse_events(self, drug_name: str, reaction: str | None = None,
                       serious: bool | None = None, limit: int = 25) -> dict:
        return FAKE_ADVERSE_EVENTS

    def recalls(self, drug_name: str, classification: str | None = None, limit: int = 25) -> dict:
        return FAKE_RECALLS


class FakeClinPGx:
    def pgx_for_drug(self, drug_name: str) -> list[dict]:
        return [{"gene": "CYP2C9", "drug": drug_name, "cpic_level": "A",
                 "guideline_url": "https://www.clinpgx.org/x"}]


class FakeDDInter:
    def interactions(self, drug_name: str, limit: int = 25) -> list[dict]:
        return [{"interacting_drug": "warfarin", "severity": "Major", "mechanism": "PD"}]


GOOD_CLIENTS = {"chembl": FakeChEMBL(), "openfda": FakeOpenFDA(),
                "clinpgx": FakeClinPGx(), "ddinter": FakeDDInter()}
BAD_CHEMBL_CLIENTS = {"chembl": FakeChEMBLUnresolvable(), "openfda": FakeOpenFDA(),
                      "clinpgx": FakeClinPGx(), "ddinter": FakeDDInter()}


# ── Module export tests ───────────────────────────────────────────────────────

def test_module_exports_entity():
    assert drug_entity.ENTITY == "drug"


def test_module_exports_tool_name():
    assert drug_entity.TOOL_NAME == "drug_profile"


def test_module_exports_description():
    assert drug_entity.DESCRIPTION
    assert "drug" in drug_entity.DESCRIPTION.lower()


def test_module_exports_needs():
    assert "chembl" in drug_entity.NEEDS
    assert "openfda" in drug_entity.NEEDS


def test_module_exports_six_sections():
    names = {s.name for s in drug_entity.SECTIONS}
    assert names == {"compound", "mechanism", "label", "adverse_events", "indications", "recalls", "pharmacogenomics", "interactions"}


def test_default_sections_are_compound_mechanism_label_adverse_events():
    defaults = {s.name for s in drug_entity.SECTIONS if s.default}
    assert defaults == {"compound", "mechanism", "label", "adverse_events", "pharmacogenomics"}


def test_non_default_sections_are_indications_and_recalls():
    non_defaults = {s.name for s in drug_entity.SECTIONS if not s.default}
    assert non_defaults == {"indications", "recalls", "interactions"}


def test_sources_correct():
    by_name = {s.name: s for s in drug_entity.SECTIONS}
    assert by_name["compound"].sources == ("chembl",)
    assert by_name["mechanism"].sources == ("chembl",)
    assert by_name["indications"].sources == ("chembl",)
    assert by_name["label"].sources == ("openfda",)
    assert by_name["adverse_events"].sources == ("openfda",)
    assert by_name["recalls"].sources == ("openfda",)


# ── resolve() tests ───────────────────────────────────────────────────────────

def test_resolve_returns_name_unchanged():
    resolved, short_circuit = drug_entity.resolve("aspirin", GOOD_CLIENTS)
    assert resolved == "aspirin"
    assert short_circuit is None


def test_resolve_does_not_short_circuit():
    _, short_circuit = drug_entity.resolve("unknowndrug", GOOD_CLIENTS)
    assert short_circuit is None


# ── fetchers() + run_sections() happy path ────────────────────────────────────

def test_default_sections_run_and_overall_ok():
    fs = drug_entity.fetchers("aspirin", GOOD_CLIENTS)
    env = run_sections("drug", "aspirin", None, fs, drug_entity.SECTIONS)
    assert env["entity"] == "drug"
    assert env["id"] == "aspirin"
    assert env["overall_status"] == "ok"
    # Only default sections should run when sections=None
    assert set(env["sections"]) == {"compound", "mechanism", "label", "adverse_events", "pharmacogenomics"}


def test_compound_section_data():
    fs = drug_entity.fetchers("aspirin", GOOD_CLIENTS)
    env = run_sections("drug", "aspirin", ["compound"], fs, drug_entity.SECTIONS)
    assert env["sections"]["compound"]["status"] == "ok"
    data = env["sections"]["compound"]["data"]
    assert data["chembl_id"] == CHEMBL_ID
    assert data["name"] == "aspirin"
    assert data["max_phase"] == 4


def test_mechanism_section_data():
    fs = drug_entity.fetchers("aspirin", GOOD_CLIENTS)
    env = run_sections("drug", "aspirin", ["mechanism"], fs, drug_entity.SECTIONS)
    assert env["sections"]["mechanism"]["status"] == "ok"
    mechs = env["sections"]["mechanism"]["data"]
    assert isinstance(mechs, list)
    assert mechs[0]["action_type"] == "INHIBITOR"


def test_label_section_data():
    fs = drug_entity.fetchers("aspirin", GOOD_CLIENTS)
    env = run_sections("drug", "aspirin", ["label"], fs, drug_entity.SECTIONS)
    assert env["sections"]["label"]["status"] == "ok"
    data = env["sections"]["label"]["data"]
    assert data["drug"] == "aspirin"
    assert "brand_names" in data


def test_adverse_events_section_data():
    fs = drug_entity.fetchers("aspirin", GOOD_CLIENTS)
    env = run_sections("drug", "aspirin", ["adverse_events"], fs, drug_entity.SECTIONS)
    assert env["sections"]["adverse_events"]["status"] == "ok"
    data = env["sections"]["adverse_events"]["data"]
    assert data["total"] == 12345


def test_indications_section_non_default_but_runnable():
    fs = drug_entity.fetchers("aspirin", GOOD_CLIENTS)
    env = run_sections("drug", "aspirin", ["indications"], fs, drug_entity.SECTIONS)
    assert env["sections"]["indications"]["status"] == "ok"
    data = env["sections"]["indications"]["data"]
    assert isinstance(data, list)
    assert data[0]["mesh_heading"] == "Pain"


def test_recalls_section_non_default_but_runnable():
    fs = drug_entity.fetchers("aspirin", GOOD_CLIENTS)
    env = run_sections("drug", "aspirin", ["recalls"], fs, drug_entity.SECTIONS)
    assert env["sections"]["recalls"]["status"] == "ok"
    data = env["sections"]["recalls"]["data"]
    assert data["total"] == 3


def test_all_six_sections_explicit():
    fs = drug_entity.fetchers("aspirin", GOOD_CLIENTS)
    all_sections = [s.name for s in drug_entity.SECTIONS]
    env = run_sections("drug", "aspirin", all_sections, fs, drug_entity.SECTIONS)
    assert env["overall_status"] == "ok"
    assert set(env["sections"]) == {"compound", "mechanism", "label", "adverse_events", "indications", "recalls", "pharmacogenomics", "interactions"}


# ── ChEMBL resolve failure → partial degradation ─────────────────────────────

def test_failed_chembl_resolve_degrades_only_chembl_sections():
    """When ChEMBL can't resolve the drug name, the three ChEMBL-backed sections
    (compound, mechanism, indications) should fail with status=error, while the
    two default OpenFDA sections (label, adverse_events) should still succeed.
    overall_status must be "partial", not "error".
    """
    fs = drug_entity.fetchers("unknowndrug", BAD_CHEMBL_CLIENTS)
    # Run all default sections (compound, mechanism, label, adverse_events)
    env = run_sections("drug", "unknowndrug", None, fs, drug_entity.SECTIONS)

    assert env["overall_status"] == "partial"

    # ChEMBL sections fail
    assert env["sections"]["compound"]["status"] == "error"
    assert env["sections"]["mechanism"]["status"] == "error"
    assert "Could not resolve" in env["sections"]["compound"]["error"]["message"]

    # OpenFDA sections succeed
    assert env["sections"]["label"]["status"] == "ok"
    assert env["sections"]["adverse_events"]["status"] == "ok"


def test_failed_chembl_resolve_indications_also_error_when_requested():
    """indications is non-default, but if explicitly requested it should also
    fail (status=error) when ChEMBL resolution fails."""
    fs = drug_entity.fetchers("unknowndrug", BAD_CHEMBL_CLIENTS)
    env = run_sections("drug", "unknowndrug", ["indications"], fs, drug_entity.SECTIONS)
    assert env["sections"]["indications"]["status"] == "error"
    assert "Could not resolve" in env["sections"]["indications"]["error"]["message"]


def test_failed_chembl_recall_still_ok():
    """recalls is OpenFDA-backed; ChEMBL failure does not affect it."""
    fs = drug_entity.fetchers("unknowndrug", BAD_CHEMBL_CLIENTS)
    env = run_sections("drug", "unknowndrug", ["recalls"], fs, drug_entity.SECTIONS)
    assert env["sections"]["recalls"]["status"] == "ok"


# ── describe() meta-tool ──────────────────────────────────────────────────────

def test_describe_drug_lists_all_sections():
    d = describe("drug", drug_entity.SECTIONS)
    assert d["entity"] == "drug"
    names = {s["name"] for s in d["sections"]}
    assert names == {"compound", "mechanism", "label", "adverse_events", "indications", "recalls", "pharmacogenomics", "interactions"}


def test_describe_carries_default_flags():
    d = describe("drug", drug_entity.SECTIONS)
    by_name = {s["name"]: s for s in d["sections"]}
    assert by_name["compound"]["default"] is True
    assert by_name["indications"]["default"] is False
    assert by_name["recalls"]["default"] is False


def test_describe_carries_source_keys():
    d = describe("drug", drug_entity.SECTIONS)
    by_name = {s["name"]: s for s in d["sections"]}
    assert by_name["compound"]["sources"] == ["chembl"]
    assert by_name["label"]["sources"] == ["openfda"]
