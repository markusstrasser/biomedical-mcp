"""Gene composite entity — the template every other entity module follows.

Each entities/<entity>.py exports this exact interface so domains/composite.py can
register composites generically (one tool per module, zero per-entity wiring):

    ENTITY       str    — entity key ("gene")
    TOOL_NAME    str    — MCP tool name ("gene_dossier")
    DESCRIPTION  str    — tool docstring shown to the agent
    NEEDS        tuple  — client keys this entity's fetchers require
    SECTIONS     tuple[Section, ...]  — static section metadata (describe_sections)
    resolve(identifier, clients) -> (resolved_id, short_circuit | None)
    fetchers(identifier, clients, **opts) -> {section_name: zero-arg callable}

Adding a source to an entity = one Section + one fetcher line here. No new MCP tool.
"""

from __future__ import annotations

from typing import Any

from biomedical_mcp.composite_core import Section, Fetcher

ENTITY = "gene"
TOOL_NAME = "gene_dossier"
DESCRIPTION = (
    "Comprehensive gene report via a standardized partial-failure envelope. "
    "Sections (default: all): nomenclature, constraint, panels, expression, "
    "phenotypes, gene_disease_validity, dosage_sensitivity, rare_diseases. "
    "Call describe_sections('gene') for sources/descriptions. One dead upstream API "
    "degrades only its section; the dossier still returns. "
    "Pass `identifier` = gene symbol (e.g. 'BRCA1', 'CYP2D6')."
)

NEEDS = ("hgnc", "gnomad", "panelapp", "gtex", "hpo", "clingen", "orphanet", "clinpgx")

SECTIONS: tuple[Section, ...] = (
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
    Section("pharmacogenomics", ("clinpgx",), True,
            "CPIC/ClinPGx gene-drug pairs (all drugs whose metabolism this gene "
            "affects), CPIC levels, guideline URLs, phenotype→recommendation rows."),
)


def resolve(identifier: str, clients: dict[str, Any]) -> tuple[str, dict | None]:
    """Resolve/validate the identifier. Returns (resolved_id, short_circuit_or_None).

    Gene symbols need no resolution. Entities that do (variant normalize, drug→chembl_id)
    override this; returning a non-None dict short-circuits the composite (e.g. needs_search).
    """
    return identifier, None


def fetchers(identifier: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    limit = opts.get("expression_limit", 10)
    c = clients
    return {
        "nomenclature": lambda: c["hgnc"].gene_names(identifier),
        "constraint": lambda: c["gnomad"].gene_constraint(identifier),
        "panels": lambda: c["panelapp"].gene_panels(identifier, confidence="all"),
        "expression": lambda: c["gtex"].top_tissues(identifier, limit=limit),
        "phenotypes": lambda: c["hpo"].gene_phenotypes(identifier),
        "gene_disease_validity": lambda: c["clingen"].gene_validity(identifier),
        "dosage_sensitivity": lambda: c["clingen"].gene_dosage(identifier),
        "rare_diseases": lambda: c["orphanet"].gene_diseases(identifier),
        "pharmacogenomics": lambda: c["clinpgx"].pgx_for_gene(identifier),
    }
