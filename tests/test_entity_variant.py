"""Offline contract tests for the variant composite entity.

No network: all clients are replaced with lightweight fake classes.
Mirrors the style of tests/test_composite_core.py.
"""

from biomedical_mcp.composite_core import run_sections
from biomedical_mcp.entities import variant as variant_entity


# ── resolve() ────────────────────────────────────────────────────────────────

def test_resolve_rsid_passes_through():
    resolved_id, sc = variant_entity.resolve("rs113488022", {})
    assert resolved_id == "rs113488022"
    assert sc is None


def test_resolve_hgvs_passes_through():
    resolved_id, sc = variant_entity.resolve("chr7:g.140453136A>T", {})
    assert resolved_id == "chr7:g.140453136A>T"
    assert sc is None


def test_resolve_gene_protein_passes_through():
    """'BRAF V600E' parses as gene_protein (not search_hint) → no short-circuit.

    normalize_variant("BRAF V600E") returns format='gene_protein', normalized='BRAF V600E'.
    The resolve() contract only short-circuits on format=='search_hint'.
    """
    resolved_id, sc = variant_entity.resolve("BRAF V600E", {})
    assert sc is None
    # normalized echoes the gene+change string, ready for downstream lookup
    assert resolved_id == "BRAF V600E"


def test_resolve_bare_protein_change_short_circuits():
    """'p.Val600Glu' (no gene context) → search_hint short-circuit."""
    resolved_id, sc = variant_entity.resolve("p.Val600Glu", {})
    assert sc is not None
    assert sc["status"] == "needs_search"
    assert resolved_id == "p.Val600Glu"


def test_resolve_residue_shorthand_short_circuits():
    """'PTPN22 620W' (residue shorthand) → search_hint short-circuit."""
    resolved_id, sc = variant_entity.resolve("PTPN22 620W", {})
    assert sc is not None
    assert sc["status"] == "needs_search"
    assert sc["gene"] == "PTPN22"
    assert sc["suggestion"] is not None


def test_resolve_freetext_short_circuits():
    """Unrecognized text → search_hint short-circuit."""
    _, sc = variant_entity.resolve("hello world", {})
    assert sc is not None
    assert sc["status"] == "needs_search"


# ── module-level constants ────────────────────────────────────────────────────

def test_entity_constants():
    assert variant_entity.ENTITY == "variant"
    assert variant_entity.TOOL_NAME == "variant_context"
    assert isinstance(variant_entity.DESCRIPTION, str) and len(variant_entity.DESCRIPTION) > 10
    assert "myvariant" in variant_entity.NEEDS
    assert "gnomad" in variant_entity.NEEDS
    assert "litvar" in variant_entity.NEEDS
    assert "gwas_catalog" in variant_entity.NEEDS


def test_sections_names_and_defaults():
    names = [s.name for s in variant_entity.SECTIONS]
    assert "annotation" in names
    assert "population_frequency" in names
    assert "literature" in names
    assert "gwas" in names
    # all four are default
    for s in variant_entity.SECTIONS:
        assert s.default is True


def test_sections_evidence_grades():
    by_name = {s.name: s for s in variant_entity.SECTIONS}
    assert by_name["annotation"].evidence_grade == "C4"
    assert by_name["population_frequency"].evidence_grade == "E5"
    assert by_name["literature"].evidence_grade == "D4"
    assert by_name["gwas"].evidence_grade == "C3"


# ── fetchers + run_sections (fake clients, no network) ───────────────────────

class FakeMyVariant:
    def lookup(self, variant_id):
        return {
            "variant_id": variant_id,
            "clinvar": {"significance": ["Pathogenic"], "gene": "BRCA1"},
            "predictions": {"gene": "BRCA1", "sift_pred": "D"},
        }


class FakeGnomAD:
    def variant_frequency(self, variant_id):
        return {"variant_id": variant_id, "genome": {"af": 0.0001}}


class FakeLitVar:
    def variant_publications(self, variant_id, limit=20):
        return {"variant_id": variant_id, "publication_count": 42, "pmids": ["12345"]}


class FakeGWASCatalog:
    def variant_associations(self, rsid, limit=25):
        return {"rsid": rsid, "association_count": 2, "associations": [{"trait": "Breast cancer"}]}


def _make_clients():
    return {
        "myvariant": FakeMyVariant(),
        "gnomad": FakeGnomAD(),
        "litvar": FakeLitVar(),
        "gwas_catalog": FakeGWASCatalog(),
    }


def test_fetchers_bind_and_run_via_run_sections():
    clients = _make_clients()
    bound = variant_entity.fetchers("rs113488022", clients)
    env = run_sections("variant", "rs113488022", None, bound, variant_entity.SECTIONS)

    assert env["entity"] == "variant"
    assert env["id"] == "rs113488022"
    assert env["overall_status"] == "ok"
    assert set(env["sections"]) == {"annotation", "population_frequency", "literature", "gwas"}

    assert env["sections"]["annotation"]["status"] == "ok"
    assert env["sections"]["population_frequency"]["status"] == "ok"
    assert env["sections"]["literature"]["status"] == "ok"
    assert env["sections"]["gwas"]["status"] == "ok"


def test_annotation_data_fields():
    clients = _make_clients()
    bound = variant_entity.fetchers("rs113488022", clients)
    env = run_sections("variant", "rs113488022", None, bound, variant_entity.SECTIONS)
    ann = env["sections"]["annotation"]["data"]
    assert ann["clinvar"]["gene"] == "BRCA1"
    assert ann["predictions"]["sift_pred"] == "D"


def test_literature_limit_opt():
    """literature_limit kwarg should be forwarded to the litvar client."""
    calls = []

    class CaptureLitVar:
        def variant_publications(self, variant_id, limit=20):
            calls.append(limit)
            return {"publication_count": 1, "pmids": []}

    clients = _make_clients()
    clients["litvar"] = CaptureLitVar()
    bound = variant_entity.fetchers("rs113488022", clients, literature_limit=3)
    bound["literature"]()
    assert calls == [3]


# ── list→first multi-allelic handling ────────────────────────────────────────

def test_annotation_list_to_first_element():
    """MyVariant sometimes returns a list for multi-allelic rsIDs — take first."""

    class MultiAllelicMyVariant:
        def lookup(self, variant_id):
            return [
                {"variant_id": "rs1_a", "clinvar": {"significance": ["Pathogenic"]}},
                {"variant_id": "rs1_b", "clinvar": {"significance": ["Benign"]}},
            ]

    clients = _make_clients()
    clients["myvariant"] = MultiAllelicMyVariant()
    bound = variant_entity.fetchers("rs1", clients)
    data = bound["annotation"]()
    assert data["variant_id"] == "rs1_a"
    assert data["clinvar"]["significance"] == ["Pathogenic"]


def test_annotation_empty_list_returns_empty_dict():
    """Empty list from MyVariant → empty dict (not a crash)."""

    class EmptyMyVariant:
        def lookup(self, variant_id):
            return []

    clients = _make_clients()
    clients["myvariant"] = EmptyMyVariant()
    bound = variant_entity.fetchers("rs_empty", clients)
    data = bound["annotation"]()
    assert data == {}


# ── partial-failure degradation ───────────────────────────────────────────────

def test_one_dead_client_partial_not_total_failure():
    """If gnomad is down, the other three sections still succeed."""

    class DeadGnomAD:
        def variant_frequency(self, variant_id):
            raise RuntimeError("gnomAD unavailable")

    clients = _make_clients()
    clients["gnomad"] = DeadGnomAD()
    bound = variant_entity.fetchers("rs113488022", clients)
    env = run_sections("variant", "rs113488022", None, bound, variant_entity.SECTIONS)

    assert env["overall_status"] == "partial"
    assert env["sections"]["population_frequency"]["status"] == "error"
    assert env["sections"]["annotation"]["status"] == "ok"
    assert env["sections"]["literature"]["status"] == "ok"
    assert env["sections"]["gwas"]["status"] == "ok"


# ── provenance carried through ─────────────────────────────────────────────

def test_provenance_sources_present():
    clients = _make_clients()
    bound = variant_entity.fetchers("rs113488022", clients)
    env = run_sections("variant", "rs113488022", None, bound, variant_entity.SECTIONS)
    prov = env["sections"]["annotation"]["provenance"]
    assert prov["sources"] == ["myvariant"]
    assert prov["evidence_grade"] == "C4"

    prov_gn = env["sections"]["population_frequency"]["provenance"]
    assert prov_gn["sources"] == ["gnomad"]
    assert prov_gn["evidence_grade"] == "E5"
