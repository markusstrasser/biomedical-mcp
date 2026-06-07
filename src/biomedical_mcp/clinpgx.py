"""ClinPGx client — CPIC pharmacogenomics guidelines via api.cpicpgx.org (PostgREST).

CPIC was rebranded to ClinPGx in 2026; the public API remains at api.cpicpgx.org/v1
(PostgREST). No authentication required. No API key.

PharmGKB note: PharmGKB's legacy REST API (api.pharmgkb.org/v1) returned 404 on all
probed endpoints as of 2026-06-07. CPIC's pair_view already carries the PharmGKB
curated evidence level via the `clinpgxlevel` field (e.g. 1A, 2A, 3), so downstream
consumers get integrated CPIC + PharmGKB signal from this single source.

Views used:
  pair_view       — gene-drug pairs with CPIC level, guideline URL, PGx testing flag
  recommendation_view — phenotype→recommendation rows keyed by drugname

CPIC level sort order (A best → D worst, then B/C, then unknowns):
  A > B > B/C > C > D (B/C is a real level in the data, sorted between B and C)
"""

from __future__ import annotations

import logging

import httpx

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache
from biomedical_mcp.errors import NotFoundError, SourceUnavailableError

log = logging.getLogger(__name__)

CPIC_BASE = "https://api.cpicpgx.org/v1"

# CPIC level ordering — lower index = higher clinical priority
_LEVEL_ORDER: dict[str, int] = {
    "A": 0,
    "B": 1,
    "B/C": 2,
    "C": 3,
    "D": 4,
}


def _level_rank(level: str | None) -> int:
    """Return sort key for a CPIC level string; unknown/null levels sort last."""
    if not level:
        return 99
    return _LEVEL_ORDER.get(level.strip().upper().replace(" ", ""), 98)


def _trim_rec(rec: dict) -> dict:
    """Flatten a recommendation_view row to the fields consumers care about."""
    return {
        "phenotypes": rec.get("phenotypes") or {},
        "activity_score": rec.get("activityscore") or {},
        "recommendation": rec.get("drugrecommendation"),
        "classification": rec.get("classification"),
        "implications": rec.get("implications") or {},
        "population": rec.get("population"),
        "comments": rec.get("comments") if rec.get("comments") not in (None, "n/a") else None,
    }


class ClinPGx(BaseClient):
    """CPIC/ClinPGx pharmacogenomics guidelines client.

    Provides gene-drug pair data (CPIC levels) and phenotype→recommendation
    lookups via the public PostgREST API at api.cpicpgx.org/v1.
    No authentication required.
    """

    domain_name = "clinpgx"

    def __init__(self, cache: Cache) -> None:
        super().__init__(cache, base_url=CPIC_BASE)
        # PostgREST uses Accept header for content negotiation
        self.client.headers["Accept"] = "application/json"

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _pairs(self, *, drug: str | None = None, gene: str | None = None) -> list[dict]:
        """Fetch pair_view rows for a drug name or gene symbol."""
        params: dict = {}
        if drug:
            params["drugname"] = f"eq.{drug}"
        elif gene:
            params["genesymbol"] = f"eq.{gene}"
        else:
            raise ValueError("Provide drug or gene")

        data, cache_hit = self._cached_get("pairs", "/pair_view", params)
        rows: list[dict] = data if isinstance(data, list) else []
        return rows

    def _recommendations_for_drug(self, drug_name: str) -> list[dict]:
        """Fetch recommendation_view rows for a specific drug."""
        params = {"drugname": f"eq.{drug_name}"}
        data, _ = self._cached_get("recs", "/recommendation_view", params)
        return data if isinstance(data, list) else []

    def _build_pair_entry(self, pair: dict, recs: list[dict]) -> dict:
        """Build a structured result dict from a pair_view row + its recommendations."""
        return {
            "gene": pair.get("genesymbol"),
            "drug": pair.get("drugname"),
            "cpic_level": pair.get("cpiclevel"),
            "clinpgx_level": pair.get("clinpgxlevel"),
            "pgx_testing": pair.get("pgxtesting"),
            "guideline_name": pair.get("guidelinename"),
            "guideline_url": pair.get("guidelineurl"),
            "used_for_recommendation": pair.get("usedforrecommendation"),
            "provisional": pair.get("provisional", False),
            "pmids": pair.get("pmids") or [],
            "recommendations": [_trim_rec(r) for r in recs],
        }

    # ── Public API ────────────────────────────────────────────────────────────

    def pgx_for_drug(self, drug_name: str) -> dict:
        """Gene-drug pairs + CPIC level + recommendations for a drug.

        Returns all genes implicated in the pharmacogenomics of the named drug,
        sorted by CPIC level (A first). Each entry includes guideline URL,
        PGx testing classification, and phenotype→recommendation rows.

        Args:
            drug_name: Generic drug name, e.g. "codeine", "clopidogrel", "warfarin".

        Returns:
            {
                "drug": str,
                "pairs": [
                    {
                        "gene": str,
                        "cpic_level": str,           # "A", "B", "B/C", "C", "D"
                        "clinpgx_level": str | None,  # "1A", "2A", "3", etc.
                        "pgx_testing": str | None,
                        "guideline_name": str | None,
                        "guideline_url": str | None,
                        "used_for_recommendation": str | None,
                        "provisional": bool,
                        "pmids": list[str],
                        "recommendations": [
                            {
                                "phenotypes": dict,
                                "activity_score": dict,
                                "recommendation": str,
                                "classification": str,
                                "implications": dict,
                                "population": str,
                                "comments": str | None,
                            }
                        ],
                    }
                ],
                "pair_count": int,
                "source": "cpic",
            }

        Raises:
            NotFoundError: if no pairs found for this drug.
            SourceUnavailableError: on network failure.
        """
        drug_lower = drug_name.strip().lower()
        try:
            pairs = self._pairs(drug=drug_lower)
        except httpx.HTTPError as exc:
            raise SourceUnavailableError(
                "ClinPGx", str(exc), "Try again or check api.cpicpgx.org status"
            ) from exc

        if not pairs:
            raise NotFoundError(
                "Drug",
                drug_name,
                "Check spelling or try the generic name. "
                "Not all drugs have CPIC guidelines.",
            )

        # Build results — fetch recs per unique drug name (usually just the one)
        drug_names_in_pairs = {p["drugname"] for p in pairs if p.get("drugname")}
        recs_by_drug: dict[str, list[dict]] = {}
        for dn in drug_names_in_pairs:
            try:
                recs_by_drug[dn] = self._recommendations_for_drug(dn)
            except Exception:
                recs_by_drug[dn] = []

        entries = [
            self._build_pair_entry(p, recs_by_drug.get(p.get("drugname", ""), []))
            for p in pairs
        ]
        entries.sort(key=lambda e: _level_rank(e["cpic_level"]))

        return {
            "drug": drug_lower,
            "pairs": entries,
            "pair_count": len(entries),
            "source": "cpic",
        }

    def pgx_for_gene(self, gene_symbol: str) -> dict:
        """Gene-drug pairs + CPIC level + recommendations for a gene.

        Returns all drugs for which this gene has pharmacogenomics implications,
        sorted by CPIC level (A first). Each entry includes guideline URL and
        phenotype→recommendation rows fetched per drug.

        Args:
            gene_symbol: HGNC gene symbol, e.g. "CYP2D6", "CYP2C19", "TPMT".

        Returns:
            {
                "gene": str,
                "pairs": [...],   # same shape as pgx_for_drug pairs list
                "pair_count": int,
                "source": "cpic",
            }

        Raises:
            NotFoundError: if no pairs found for this gene.
            SourceUnavailableError: on network failure.
        """
        gene_upper = gene_symbol.strip().upper()
        try:
            pairs = self._pairs(gene=gene_upper)
        except httpx.HTTPError as exc:
            raise SourceUnavailableError(
                "ClinPGx", str(exc), "Try again or check api.cpicpgx.org status"
            ) from exc

        if not pairs:
            raise NotFoundError(
                "Gene",
                gene_symbol,
                "Check the gene symbol (use HGNC format, e.g. CYP2D6). "
                "Not all genes have CPIC guidelines.",
            )

        # Fetch recommendations per drug; cache means repeated calls are free
        recs_by_drug: dict[str, list[dict]] = {}
        for p in pairs:
            dn = p.get("drugname")
            if dn and dn not in recs_by_drug:
                try:
                    recs_by_drug[dn] = self._recommendations_for_drug(dn)
                except Exception:
                    recs_by_drug[dn] = []

        entries = [
            self._build_pair_entry(p, recs_by_drug.get(p.get("drugname", ""), []))
            for p in pairs
        ]
        entries.sort(key=lambda e: _level_rank(e["cpic_level"]))

        return {
            "gene": gene_upper,
            "pairs": entries,
            "pair_count": len(entries),
            "source": "cpic",
        }

    def validate(self) -> bool:
        """Health check — probe pair_view for one row."""
        try:
            self._get("/pair_view", params={"limit": "1"})
            return True
        except Exception:
            return False
