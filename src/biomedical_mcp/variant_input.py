"""Variant input normalization — handle multiple notation formats.

Ported from BioMCP's classify_variant_input / parse_variant_id pattern.
Supports: rsIDs, HGVS genomic, gene+protein change ("BRAF V600E"),
long-form ("BRAF p.Val600Glu"), prefixed short ("BRAF p.V600E").
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Three-letter → one-letter amino acid codes
_AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

# Regex patterns
_RSID_RE = re.compile(r"^(rs\d+)$", re.IGNORECASE)
_HGVS_GENOMIC_RE = re.compile(r"^(chr\w+:g\.\d+\w+>\w+)$", re.IGNORECASE)
# "BRAF V600E" or "EGFR L858R" — gene + single-letter protein change
_GENE_PROTEIN_RE = re.compile(r"^([A-Z][A-Z0-9]+)\s+([A-Z]\d+[A-Z*])$")
# "BRAF p.V600E" — gene + prefixed short protein change
_GENE_PROT_PREFIX_RE = re.compile(r"^([A-Z][A-Z0-9]+)\s+p\.([A-Z]\d+[A-Z*])$", re.IGNORECASE)
# "BRAF p.Val600Glu" — gene + long-form three-letter protein change
_GENE_PROT_LONG_RE = re.compile(
    r"^([A-Z][A-Z0-9]+)\s+p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)$"
)
# Standalone "p.Val600Glu" or "p.V600E"
_PROT_ONLY_LONG_RE = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)$")
_PROT_ONLY_SHORT_RE = re.compile(r"^p\.([A-Z]\d+[A-Z*])$", re.IGNORECASE)
# Gene + residue shorthand "PTPN22 620W" — number + single letter
_GENE_RESIDUE_RE = re.compile(r"^([A-Z][A-Z0-9]+)\s+(\d+[A-Z])$")
# Standalone protein change "R620W"
_BARE_PROT_RE = re.compile(r"^([A-Z])(\d+)([A-Z*])$")


@dataclass
class VariantId:
    """Parsed variant identifier."""
    format: str  # "rsid", "hgvs_genomic", "gene_protein", "search_hint"
    raw: str
    normalized: str  # the ID to pass to APIs
    gene: str | None = None
    protein_change: str | None = None
    suggestion: str | None = None  # actionable hint if format is "search_hint"


def _three_to_one(aa3: str) -> str:
    """Convert three-letter amino acid code to one-letter."""
    return _AA3TO1.get(aa3, aa3)


def normalize_variant(raw_input: str) -> VariantId:
    """Parse and normalize a variant identifier from various formats.

    Returns a VariantId with the normalized ID ready for API lookup,
    or a search_hint with suggestion if the input needs a search step.
    """
    s = raw_input.strip()
    if not s:
        return VariantId("search_hint", s, s, suggestion="Empty variant ID. Provide an rsID (rs1234), HGVS (chr7:g.140453136A>T), or gene+change (BRAF V600E).")

    # rsID: rs113488022
    m = _RSID_RE.match(s)
    if m:
        return VariantId("rsid", s, m.group(1).lower())

    # HGVS genomic: chr7:g.140453136A>T
    m = _HGVS_GENOMIC_RE.match(s)
    if m:
        return VariantId("hgvs_genomic", s, m.group(1))

    # Gene + short protein change: BRAF V600E
    m = _GENE_PROTEIN_RE.match(s)
    if m:
        gene, change = m.group(1), m.group(2)
        return VariantId("gene_protein", s, s, gene=gene, protein_change=change)

    # Gene + prefixed short: BRAF p.V600E
    m = _GENE_PROT_PREFIX_RE.match(s)
    if m:
        gene, change = m.group(1), m.group(2)
        return VariantId("gene_protein", s, f"{gene} {change}", gene=gene, protein_change=change)

    # Gene + long-form: BRAF p.Val600Glu
    m = _GENE_PROT_LONG_RE.match(s)
    if m:
        gene = m.group(1)
        aa_from = _three_to_one(m.group(2))
        pos = m.group(3)
        aa_to = _three_to_one(m.group(4))
        change = f"{aa_from}{pos}{aa_to}"
        return VariantId("gene_protein", s, f"{gene} {change}", gene=gene, protein_change=change)

    # Standalone p.Val600Glu
    m = _PROT_ONLY_LONG_RE.match(s)
    if m:
        aa_from = _three_to_one(m.group(1))
        pos = m.group(2)
        aa_to = _three_to_one(m.group(3))
        change = f"{aa_from}{pos}{aa_to}"
        return VariantId("search_hint", s, change, protein_change=change,
                         suggestion=f"Protein change without gene. Try: variants_search with hgvsp={change}")

    # Standalone p.V600E
    m = _PROT_ONLY_SHORT_RE.match(s)
    if m:
        change = m.group(1)
        return VariantId("search_hint", s, change, protein_change=change,
                         suggestion=f"Protein change without gene. Try: variants_search with hgvsp={change}")

    # Gene + residue alias: PTPN22 620W
    m = _GENE_RESIDUE_RE.match(s)
    if m:
        gene, alias = m.group(1), m.group(2)
        return VariantId("search_hint", s, s, gene=gene,
                         suggestion=f"Residue shorthand — try: variants_search with gene={gene} to find the variant.")

    # Bare protein change: R620W
    m = _BARE_PROT_RE.match(s)
    if m:
        change = s
        return VariantId("search_hint", s, change, protein_change=change,
                         suggestion=f"Bare protein change — try: variants_search with hgvsp={change}")

    # Anything else — treat as free text search
    return VariantId("search_hint", s, s,
                     suggestion=f"Unrecognized format. Try: variants_search with query=\"{s}\"")
