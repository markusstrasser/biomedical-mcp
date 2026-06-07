"""Section registry — the single source of truth for composite entities.

Two halves, intentionally separated:
  - SECTIONS: static metadata per entity (drives describe_sections, no clients).
  - *_fetchers(identifier, clients, **opts): bind each section to a client call.

Adding a new source = add a Section to SECTIONS + a line to the entity's
*_fetchers builder. No new MCP tool, no server.py edit. This is how CPIC/PharmGKB
(PGx), DDInter (interactions), and CIViC/OncoKB (somatic) will land.

Entities: gene (done), variant, drug, protein, disease (filled per-entity).
"""

from __future__ import annotations

from typing import Any

from biomedical_mcp.composite_core import Section, Fetcher

# ── Static section metadata ──────────────────────────────────────────────────
# Pure data. `describe_sections(entity)` returns this; it carries no per-call cost.

SECTIONS: dict[str, tuple[Section, ...]] = {
    "gene": (
        Section("nomenclature", ("hgnc",), True,
                "Official HGNC symbol, full name, aliases, previous symbols, gene family."),
        Section("constraint", ("gnomad",), True,
                "gnomAD loss-of-function constraint (pLI, LOEUF, missense Z).", "E5"),
        Section("panels", ("panelapp",), True,
                "Curated clinical gene panels the gene appears on (PanelApp)."),
        Section("expression", ("gtex",), True,
                "Top tissues by expression (GTEx).", params=("limit",)),
        Section("phenotypes", ("hpo",), True,
                "Gene→phenotype associations (HPO)."),
        Section("gene_disease_validity", ("clingen",), True,
                "ClinGen gene-disease validity classification.", "A1"),
        Section("dosage_sensitivity", ("clingen",), True,
                "ClinGen haploinsufficiency / triplosensitivity dosage scores."),
        Section("rare_diseases", ("orphanet",), True,
                "Orphanet rare-disease associations, inheritance, epidemiology.", "B2"),
    ),
    # Filled by per-entity workers against the same contract (Phase 2):
    "variant": (),
    "drug": (),
    "protein": (),
    "disease": (),
}


# ── Fetcher builders ─────────────────────────────────────────────────────────
# Each returns {section_name: zero-arg callable}. `clients` is a dict of
# instantiated client objects (built once in domains/composite.py).

def gene_fetchers(symbol: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    limit = opts.get("expression_limit", 10)
    c = clients
    return {
        "nomenclature": lambda: c["hgnc"].gene_names(symbol),
        "constraint": lambda: c["gnomad"].gene_constraint(symbol),
        "panels": lambda: c["panelapp"].gene_panels(symbol, confidence="all"),
        "expression": lambda: c["gtex"].top_tissues(symbol, limit=limit),
        "phenotypes": lambda: c["hpo"].gene_phenotypes(symbol),
        "gene_disease_validity": lambda: c["clingen"].gene_validity(symbol),
        "dosage_sensitivity": lambda: c["clingen"].gene_dosage(symbol),
        "rare_diseases": lambda: c["orphanet"].gene_diseases(symbol),
    }


# Stubs filled in Phase 2 (variant_context already has logic to port):
def variant_fetchers(variant_id: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    raise NotImplementedError("variant sections — Phase 2")


def drug_fetchers(drug: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    raise NotImplementedError("drug sections — Phase 2")


def protein_fetchers(protein_id: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    raise NotImplementedError("protein sections — Phase 2")


def disease_fetchers(disease_id: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    raise NotImplementedError("disease sections — Phase 2")
