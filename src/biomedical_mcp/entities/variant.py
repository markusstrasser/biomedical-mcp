"""Variant composite entity — conforms to the entities/gene.py interface.

Sections:
  annotation          myvariant.lookup → ClinVar/predictions/CGI/COSMIC   (C4)
  population_frequency gnomad.variant_frequency                            (E5)
  literature          litvar.variant_publications                          (D4)
  gwas                gwas_catalog.variant_associations                    (C3)

resolve() calls normalize_variant() and short-circuits for gene+protein inputs
(format == "search_hint") that cannot be resolved to a direct API identifier
without a prior search step.
"""

from __future__ import annotations

from typing import Any

from biomedical_mcp.composite_core import Section, Fetcher
from biomedical_mcp.variant_input import normalize_variant

ENTITY = "variant"
TOOL_NAME = "variant_context"
DESCRIPTION = (
    "Comprehensive variant report via a standardized partial-failure envelope. "
    "Sections (default: all): annotation, population_frequency, literature, gwas. "
    "Call describe_sections('variant') for sources/descriptions. "
    "Accepts: rsID (rs1801133), HGVS (chr7:g.140453136A>T), gene+change (BRAF V600E). "
    "Gene+protein inputs that cannot be directly resolved return status='needs_search' "
    "with a suggestion to call variants_search first."
)

NEEDS = ("myvariant", "gnomad", "litvar", "gwas_catalog", "somatic")

SECTIONS: tuple[Section, ...] = (
    Section(
        "annotation",
        ("myvariant",),
        True,
        "ClinVar significance/conditions, pathogenicity predictions (SIFT, PolyPhen2, CADD), "
        "CGI drug associations, COSMIC somatic context.",
        "C4",
    ),
    Section(
        "population_frequency",
        ("gnomad",),
        True,
        "gnomAD r4 allele frequencies across global populations (genome + exome).",
        "E5",
    ),
    Section(
        "literature",
        ("litvar",),
        True,
        "Variant-level publication count and PMIDs from LitVar2 (NCBI).",
        "D4",
        params=("literature_limit",),
    ),
    Section(
        "gwas",
        ("gwas_catalog",),
        True,
        "GWAS Catalog trait associations for the variant (rsID-based).",
        "C3",
    ),
    Section(
        "somatic",
        ("somatic",),
        False,  # oncology-specific — opt-in, not run by default
        "Somatic/oncology annotation: CIViC clinical evidence (type, level, "
        "significance, disease, therapies) + OncoKB oncogenicity & treatment levels. "
        "Needs gene+protein-change input (e.g. 'BRAF V600E'); OncoKB needs ONCOKB_TOKEN.",
    ),
)


def resolve(identifier: str, clients: dict[str, Any]) -> tuple[str, dict | None]:
    """Normalize the variant identifier.

    Returns (normalized_id, None) for directly queryable formats (rsid, hgvs_genomic).
    Returns (raw, short_circuit_dict) for gene+protein inputs that need a search step.
    """
    parsed = normalize_variant(identifier)
    if parsed.format == "search_hint":
        return identifier, {
            "status": "needs_search",
            "parsed_format": parsed.format,
            "suggestion": parsed.suggestion,
            "gene": parsed.gene,
            "protein_change": parsed.protein_change,
            "raw": parsed.raw,
        }
    return parsed.normalized, None


def _annotation_fetcher(identifier: str, myvariant_client: Any) -> dict:
    """Call myvariant.lookup; collapse multi-allelic list to first entry."""
    result = myvariant_client.lookup(identifier)
    if isinstance(result, list):
        result = result[0] if result else {}
    return result


def _somatic_fetcher(identifier: str, somatic_client: Any) -> dict:
    """Somatic annotation for gene+protein-change inputs ('BRAF V600E').

    Re-parses the identifier; for rsID/HGVS (no gene+change) returns an
    informative empty result rather than a hard failure.
    """
    parsed = normalize_variant(identifier)
    if parsed.gene and parsed.protein_change:
        return somatic_client.annotate(parsed.gene, parsed.protein_change)
    return {
        "gene": parsed.gene,
        "protein_change": parsed.protein_change,
        "civic": [],
        "oncokb": {"error": "Somatic annotation needs gene+protein-change input "
                            "(e.g. 'BRAF V600E')."},
        "note": f"Input format '{parsed.format}' lacks gene+protein_change for CIViC/OncoKB.",
    }


def fetchers(identifier: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    lit_limit = opts.get("literature_limit", 5)
    c = clients
    return {
        "annotation": lambda: _annotation_fetcher(identifier, c["myvariant"]),
        "population_frequency": lambda: c["gnomad"].variant_frequency(identifier),
        "literature": lambda: c["litvar"].variant_publications(identifier, limit=lit_limit),
        "gwas": lambda: c["gwas_catalog"].variant_associations(identifier),
        "somatic": lambda: _somatic_fetcher(identifier, c["somatic"]),
    }
