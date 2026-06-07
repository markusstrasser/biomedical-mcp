"""Offline contract tests for the protein composite entity.

No network: all clients are fakes. Validates:
  - resolve() passes accessions through unchanged
  - resolve() turns a gene symbol into a UniProt accession
  - resolve() short-circuits cleanly on resolution failure
  - fetchers() binds correctly and runs through run_sections with fakes
  - Missing uniprot client during resolution returns a clean error dict
"""

from biomedical_mcp.composite_core import run_sections
from biomedical_mcp.entities import protein as protein_entity

# ── fake clients ─────────────────────────────────────────────────────────────

FAKE_ACCESSION = "P38398"
FAKE_GENE = "BRCA1"


class FakeUniProt:
    """Returns deterministic protein data; supports both accession= and gene_symbol=."""

    def protein(self, accession=None, gene_symbol=None):
        # Resolution path: called with gene_symbol only
        if gene_symbol and not accession:
            if gene_symbol.upper() == FAKE_GENE:
                return {"accession": FAKE_ACCESSION, "gene": FAKE_GENE,
                        "protein_name": "Breast cancer type 1 susceptibility protein",
                        "function": "DNA repair via homologous recombination."}
            return None  # simulate not-found for unknown symbols
        # Direct accession path
        return {"accession": accession or FAKE_ACCESSION, "gene": FAKE_GENE,
                "protein_name": "Breast cancer type 1 susceptibility protein",
                "function": "DNA repair via homologous recombination."}

    def variants(self, accession=None, gene_symbol=None, limit=50):
        return {"accession": accession or FAKE_ACCESSION, "gene": FAKE_GENE,
                "variant_count": 2,
                "variants": [{"type": "Natural variant", "description": "V→M at 1313"}]}


class FakeAlphaFold:
    def prediction(self, uniprot_id):
        return {"uniprot_id": uniprot_id, "mean_plddt": 88.4, "model_url": "https://example.com/af.pdb"}


class FakePDB:
    def structures_for_gene(self, uniprot_id):
        return {"uniprot_id": uniprot_id, "total_structures": 3,
                "structures": [{"pdb_id": "1JNX", "method": "X-RAY DIFFRACTION"}]}


class FakeStringDB:
    def resolve_id(self, protein, species="9606"):
        return f"9606.ENSP00000{protein[:4].upper()}"

    def functional_enrichment(self, proteins, species="9606"):
        return {"proteins": proteins, "species": species,
                "categories": {"Process": [{"term": "GO:0006281", "description": "DNA repair",
                                            "fdr": 0.001, "gene_count": 1}]}}


class FakeInterPro:
    def protein_domains(self, uniprot_id=None, gene_symbol=None):
        return {"query": uniprot_id or gene_symbol, "entry_count": 2,
                "entries": [{"accession": "IPR002093", "name": "BRCA2/BRCA1 domain",
                             "type": "domain", "source_database": "SMART"}]}


def _make_clients():
    return {
        "uniprot": FakeUniProt(),
        "alphafold": FakeAlphaFold(),
        "pdb": FakePDB(),
        "stringdb": FakeStringDB(),
        "interpro": FakeInterPro(),
    }


# ── resolve() tests ──────────────────────────────────────────────────────────

def test_resolve_accession_passthrough():
    """A bare UniProt accession should pass through without calling any client."""
    acc, sc = protein_entity.resolve("P38398", {})
    assert acc == "P38398"
    assert sc is None


def test_resolve_accession_passthrough_O_prefix():
    acc, sc = protein_entity.resolve("O60260", {})
    assert acc == "O60260"
    assert sc is None


def test_resolve_gene_symbol_finds_accession():
    """Gene symbol should be resolved to accession via uniprot.protein(gene_symbol=...)."""
    clients = _make_clients()
    acc, sc = protein_entity.resolve("BRCA1", clients)
    assert sc is None
    assert acc == FAKE_ACCESSION


def test_resolve_gene_symbol_unknown_short_circuits():
    """Unknown gene symbol → short-circuit error dict, not an exception."""
    clients = _make_clients()
    identifier, sc = protein_entity.resolve("UNKNOWNGENE999", clients)
    assert sc is not None
    assert sc["error"] == "resolution_failed"
    assert "UNKNOWNGENE999" in sc["message"] or sc["identifier"] == "UNKNOWNGENE999"


def test_resolve_missing_uniprot_client_short_circuits():
    """Missing uniprot client → clean short-circuit, not KeyError."""
    identifier, sc = protein_entity.resolve("BRCA1", {})
    assert sc is not None
    assert sc["error"] == "resolution_failed"
    assert "uniprot client" in sc["message"]


def test_resolve_client_raises_short_circuits():
    """If uniprot.protein() raises, resolve() returns a clean short-circuit."""
    class BoomUniProt:
        def protein(self, accession=None, gene_symbol=None):
            raise RuntimeError("timeout")

    clients = {**_make_clients(), "uniprot": BoomUniProt()}
    identifier, sc = protein_entity.resolve("BRCA1", clients)
    assert sc is not None
    assert sc["error"] == "resolution_failed"
    assert "timeout" in sc["message"]


# ── fetchers + run_sections tests ────────────────────────────────────────────

def test_fetchers_bind_and_run_default_sections():
    """All default sections run and return ok data through the composite envelope."""
    clients = _make_clients()
    bound = protein_entity.fetchers(FAKE_ACCESSION, clients)
    env = run_sections(
        protein_entity.ENTITY,
        FAKE_ACCESSION,
        None,  # run defaults
        bound,
        protein_entity.SECTIONS,
    )
    assert env["entity"] == "protein"
    assert env["id"] == FAKE_ACCESSION
    assert env["overall_status"] == "ok"

    default_names = {s.name for s in protein_entity.SECTIONS if s.default}
    assert set(env["sections"]) == default_names

    assert env["sections"]["function"]["status"] == "ok"
    assert env["sections"]["function"]["data"]["gene"] == FAKE_GENE
    assert env["sections"]["variants"]["data"]["variant_count"] == 2
    assert env["sections"]["structure"]["data"]["mean_plddt"] == 88.4
    assert env["sections"]["experimental_structures"]["data"]["total_structures"] == 3
    assert env["sections"]["domains"]["data"]["entry_count"] == 2


def test_fetchers_interactions_opt_in():
    """interactions section is default=False; runs when explicitly requested."""
    clients = _make_clients()
    bound = protein_entity.fetchers(FAKE_ACCESSION, clients)
    env = run_sections(
        protein_entity.ENTITY,
        FAKE_ACCESSION,
        ["interactions"],
        bound,
        protein_entity.SECTIONS,
    )
    assert env["sections"]["interactions"]["status"] == "ok"
    cats = env["sections"]["interactions"]["data"]["categories"]
    assert "Process" in cats


def test_fetchers_one_dead_section_is_partial():
    """One failing fetcher degrades only its section; envelope stays partial not error."""
    class DeadAlphaFold:
        def prediction(self, uid):
            raise ConnectionError("AlphaFold down")

    clients = {**_make_clients(), "alphafold": DeadAlphaFold()}
    bound = protein_entity.fetchers(FAKE_ACCESSION, clients)
    env = run_sections(
        protein_entity.ENTITY,
        FAKE_ACCESSION,
        None,
        bound,
        protein_entity.SECTIONS,
    )
    assert env["overall_status"] == "partial"
    assert env["sections"]["structure"]["status"] == "error"
    assert env["sections"]["function"]["status"] == "ok"


def test_interactions_resolve_failure_returns_error_dict_not_exception():
    """If STRING-DB resolve_id returns None, interactions returns error dict (not raises)."""
    class UnresolvableStringDB:
        def resolve_id(self, protein, species="9606"):
            return None

        def functional_enrichment(self, proteins, species="9606"):
            raise AssertionError("should not be called")

    clients = {**_make_clients(), "stringdb": UnresolvableStringDB()}
    bound = protein_entity.fetchers(FAKE_ACCESSION, clients)
    env = run_sections(
        protein_entity.ENTITY,
        FAKE_ACCESSION,
        ["interactions"],
        bound,
        protein_entity.SECTIONS,
    )
    # error dict from client → classified as "empty" by composite core
    assert env["sections"]["interactions"]["status"] == "empty"


# ── metadata contract ────────────────────────────────────────────────────────

def test_protein_metadata_exports():
    assert protein_entity.ENTITY == "protein"
    assert protein_entity.TOOL_NAME == "protein_profile"
    assert isinstance(protein_entity.DESCRIPTION, str) and len(protein_entity.DESCRIPTION) > 20
    assert "uniprot" in protein_entity.NEEDS
    assert "alphafold" in protein_entity.NEEDS
    assert "pdb" in protein_entity.NEEDS
    assert "stringdb" in protein_entity.NEEDS
    assert "interpro" in protein_entity.NEEDS


def test_sections_have_expected_names():
    names = {s.name for s in protein_entity.SECTIONS}
    assert "function" in names
    assert "variants" in names
    assert "structure" in names
    assert "experimental_structures" in names
    assert "domains" in names
    assert "interactions" in names


def test_interactions_is_not_default():
    sec = next(s for s in protein_entity.SECTIONS if s.name == "interactions")
    assert sec.default is False


def test_default_sections_are_default():
    for name in ("function", "variants", "structure", "experimental_structures", "domains"):
        sec = next(s for s in protein_entity.SECTIONS if s.name == name)
        assert sec.default is True, f"expected {name} to be default"
