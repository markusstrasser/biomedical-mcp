"""Offline contract tests for the composite core (envelope + section runner + registry).

No network: fetchers are plain callables / fake clients. Proves the consolidation
contract that all entity composites conform to.
"""

from biomedical_mcp.composite_core import Section, run_sections, describe
from biomedical_mcp.entities import gene as gene_entity


# ── envelope + status logic ──────────────────────────────────────────────────

META = (
    Section("a", ("src_a",), True, "section a"),
    Section("b", ("src_b",), True, "section b"),
    Section("c", ("src_c",), False, "non-default section c"),
)


def test_all_ok_overall_ok():
    fetchers = {"a": lambda: {"x": 1}, "b": lambda: [1, 2]}
    env = run_sections("thing", "id1", None, fetchers, META)
    assert env["entity"] == "thing"
    assert env["id"] == "id1"
    assert env["overall_status"] == "ok"
    # defaults only (a, b) — c is not default and not requested
    assert set(env["sections"]) == {"a", "b"}
    assert env["sections"]["a"]["status"] == "ok"
    assert env["sections"]["a"]["data"] == {"x": 1}
    assert env["sections"]["a"]["provenance"]["sources"] == ["src_a"]


def test_one_dead_section_is_partial_not_total_failure():
    def boom():
        raise RuntimeError("upstream 500")

    fetchers = {"a": lambda: {"x": 1}, "b": boom}
    env = run_sections("thing", "id1", None, fetchers, META)
    assert env["overall_status"] == "partial"
    assert env["sections"]["a"]["status"] == "ok"
    assert env["sections"]["b"]["status"] == "error"
    assert env["sections"]["b"]["error"]["type"] == "RuntimeError"
    assert env["sections"]["b"]["error"]["message"] == "upstream 500"
    assert env["sections"]["b"]["data"] is None


def test_all_dead_is_error():
    def boom():
        raise ValueError("nope")

    fetchers = {"a": boom, "b": boom}
    env = run_sections("thing", "id1", None, fetchers, META)
    assert env["overall_status"] == "error"


def test_empty_result_classified_empty_not_error():
    fetchers = {"a": lambda: None, "b": lambda: {}}
    env = run_sections("thing", "id1", None, fetchers, META)
    # empty still counts as a successful (non-failing) fetch → overall ok
    assert env["overall_status"] == "ok"
    assert env["sections"]["a"]["status"] == "empty"
    assert env["sections"]["b"]["status"] == "empty"


def test_error_dict_from_client_classified_empty():
    # clients sometimes return {"error": "..."} instead of raising
    fetchers = {"a": lambda: {"error": "not found"}, "b": lambda: {"ok": 1}}
    env = run_sections("thing", "id1", None, fetchers, META)
    assert env["sections"]["a"]["status"] == "empty"


def test_requested_subset_only_runs_those():
    fetchers = {"a": lambda: {"x": 1}, "b": lambda: {"y": 2}, "c": lambda: {"z": 3}}
    env = run_sections("thing", "id1", ["c"], fetchers, META)
    assert set(env["sections"]) == {"c"}
    assert env["sections"]["c"]["status"] == "ok"


def test_unknown_section_is_not_applicable():
    fetchers = {"a": lambda: {"x": 1}}
    env = run_sections("thing", "id1", ["a", "bogus"], fetchers, META)
    assert env["sections"]["bogus"]["status"] == "not_applicable"
    assert env["sections"]["bogus"]["error"]["type"] == "unknown_section"
    # a ok + bogus not_applicable → partial
    assert env["overall_status"] == "partial"


# ── describe + registry ──────────────────────────────────────────────────────

def test_describe_gene_lists_sections():
    d = describe("gene", gene_entity.SECTIONS)
    assert d["entity"] == "gene"
    names = {s["name"] for s in d["sections"]}
    assert "nomenclature" in names
    assert "constraint" in names
    assert "rare_diseases" in names
    # metadata carried
    constraint = next(s for s in d["sections"] if s["name"] == "constraint")
    assert constraint["sources"] == ["gnomad"]
    assert constraint["evidence_grade"] == "E5"


def test_gene_fetchers_bind_to_clients_and_run():
    # Fake clients exercise the real gene_fetchers + run_sections path, no network.
    class FakeHGNC:
        def gene_names(self, s): return {"symbol": s, "name": "BRCA1 DNA repair"}

    class FakeGnomad:
        def gene_constraint(self, s): return {"pLI": 0.99, "loeuf": 0.2}

    class FakePanelApp:
        def gene_panels(self, s, confidence): return [{"panel": "Hereditary cancer"}]

    class FakeGTEx:
        def top_tissues(self, s, limit): return [{"tissue": "breast"}]

    class FakeHPO:
        def gene_phenotypes(self, s): return [{"hp": "HP:0000001"}]

    class FakeClinGen:
        def gene_validity(self, s): return {"classification": "Definitive"}
        def gene_dosage(self, s): return {"haploinsufficiency": 3}

    class FakeOrphanet:
        def gene_diseases(self, s): return [{"orpha": "ORPHA:145"}]

    clients = {
        "hgnc": FakeHGNC(), "gnomad": FakeGnomad(), "panelapp": FakePanelApp(),
        "gtex": FakeGTEx(), "hpo": FakeHPO(), "clingen": FakeClinGen(),
        "orphanet": FakeOrphanet(),
    }
    fetchers = gene_entity.fetchers("BRCA1", clients)
    env = run_sections("gene", "BRCA1", None, fetchers, gene_entity.SECTIONS)
    assert env["overall_status"] == "ok"
    assert env["sections"]["nomenclature"]["data"]["symbol"] == "BRCA1"
    assert env["sections"]["constraint"]["data"]["pLI"] == 0.99
    assert env["sections"]["dosage_sensitivity"]["data"]["haploinsufficiency"] == 3
    # every default gene section ran
    assert set(env["sections"]) == {s.name for s in gene_entity.SECTIONS}
