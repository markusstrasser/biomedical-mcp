"""Protein composite entity.

Each entities/<entity>.py exports this exact interface so domains/composite.py can
register composites generically (one tool per module, zero per-entity wiring):

    ENTITY       str    — entity key ("protein")
    TOOL_NAME    str    — MCP tool name ("protein_profile")
    DESCRIPTION  str    — tool docstring shown to the agent
    NEEDS        tuple  — client keys this entity's fetchers require
    SECTIONS     tuple[Section, ...]  — static section metadata (describe_sections)
    resolve(identifier, clients) -> (resolved_id, short_circuit | None)
    fetchers(identifier, clients, **opts) -> {section_name: zero-arg callable}

identifier = UniProt accession (e.g. "P38398") or gene symbol (e.g. "BRCA1").
resolve() detects by pattern and resolves gene symbols to accessions via UniProt.
"""

from __future__ import annotations

import re
from typing import Any

from biomedical_mcp.composite_core import Section, Fetcher

ENTITY = "protein"
TOOL_NAME = "protein_profile"
DESCRIPTION = (
    "Comprehensive protein report via a standardized partial-failure envelope. "
    "identifier = UniProt accession (e.g. 'P38398', 'O60260') or gene symbol (e.g. 'BRCA1'). "
    "Gene symbols are resolved to a UniProt accession automatically. "
    "Sections (default on): function, variants, structure, experimental_structures, domains. "
    "Sections (default off): interactions (STRING functional enrichment, heavier). "
    "Call describe_sections('protein') for per-section sources and descriptions. "
    "One dead upstream API degrades only its section; the profile still returns."
)

NEEDS = ("uniprot", "alphafold", "pdb", "stringdb", "interpro")

SECTIONS: tuple[Section, ...] = (
    Section("function", ("uniprot",), True,
            "UniProt protein summary: function, subcellular location, disease links, "
            "GO terms, active sites, and cross-reference database counts."),
    Section("variants", ("uniprot",), True,
            "Known natural variants and mutagenesis data from UniProt features."),
    Section("structure", ("alphafold",), True,
            "AlphaFold predicted structure: mean pLDDT confidence, model/CIF/PAE URLs."),
    Section("experimental_structures", ("pdb",), True,
            "Experimental PDB structures (X-ray, cryo-EM) linked to the UniProt accession."),
    Section("domains", ("interpro",), True,
            "InterPro domain, family, and site annotations for the protein accession."),
    Section("interactions", ("stringdb",), False,
            "STRING-DB functional enrichment (GO, KEGG, Reactome) for the resolved protein. "
            "Default=False — heavier enrichment query; enable explicitly when pathway context is needed."),
)

# Regex for a bare UniProt accession: letter + digits [+ more segments].
# Covers canonical format like P38398, O60260, Q9UER7, A0A000.
_ACCESSION_RE = re.compile(r"^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$",
                           re.IGNORECASE)


def _looks_like_accession(identifier: str) -> bool:
    """Return True if identifier matches a UniProt accession pattern."""
    return bool(_ACCESSION_RE.match(identifier.strip()))


def resolve(identifier: str, clients: dict[str, Any]) -> tuple[str, dict | None]:
    """Resolve identifier to a UniProt accession.

    - If identifier is already a UniProt accession (by pattern), return it as-is.
    - Otherwise treat as gene symbol: call uniprot.protein(gene_symbol=...) to find
      the canonical accession, then return it.
    - If resolution fails, return (identifier, short_circuit_error) so the composite
      reports cleanly rather than silently propagating an unresolved id.
    """
    identifier = identifier.strip()
    if _looks_like_accession(identifier):
        return identifier, None

    # Resolve gene symbol → accession via UniProt.
    uniprot = clients.get("uniprot")
    if uniprot is None:
        return identifier, {
            "error": "resolution_failed",
            "message": "uniprot client not available for gene symbol resolution",
            "identifier": identifier,
        }

    try:
        protein_data = uniprot.protein(gene_symbol=identifier)
    except Exception as exc:
        return identifier, {
            "error": "resolution_failed",
            "message": f"UniProt resolution raised: {exc}",
            "identifier": identifier,
        }

    if not protein_data:
        return identifier, {
            "error": "resolution_failed",
            "message": f"No reviewed UniProt entry found for gene symbol '{identifier}'",
            "identifier": identifier,
        }

    accession = protein_data.get("accession")
    if not accession:
        return identifier, {
            "error": "resolution_failed",
            "message": f"UniProt returned data without an accession for '{identifier}'",
            "identifier": identifier,
        }

    return accession, None


def fetchers(identifier: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    """Build section fetchers for the resolved UniProt accession.

    identifier: resolved UniProt accession (from resolve()).
    opts:
        species (str): STRING-DB taxon id, default "9606".
    """
    c = clients
    accession = identifier
    species = opts.get("species", "9606")

    def _interactions() -> dict:
        string_id = c["stringdb"].resolve_id(accession, species)
        if not string_id:
            return {"error": f"Could not resolve {accession!r} in STRING-DB"}
        return c["stringdb"].functional_enrichment([string_id], species)

    return {
        "function": lambda: c["uniprot"].protein(accession=accession),
        "variants": lambda: c["uniprot"].variants(accession=accession),
        "structure": lambda: c["alphafold"].prediction(accession),
        "experimental_structures": lambda: c["pdb"].structures_for_gene(accession),
        "domains": lambda: c["interpro"].protein_domains(uniprot_id=accession),
        "interactions": _interactions,
    }
