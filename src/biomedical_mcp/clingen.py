"""ClinGen gene-disease validity and dosage sensitivity — CSV bulk download client.

ClinGen doesn't expose a REST API for gene-disease validity curations.
Data is downloaded as CSV from their web interface and cached locally in SQLite.
"""

from __future__ import annotations

import csv
import io
import logging
import time

import httpx

from biomedical_mcp.cache import Cache
from biomedical_mcp.provenance import make_provenance

log = logging.getLogger(__name__)

# CSV download endpoints (require X-Requested-With header)
VALIDITY_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
DOSAGE_URL = "https://search.clinicalgenome.org/kb/gene-dosage/download"

DOWNLOAD_HEADERS = {
    "X-Requested-With": "XMLHttpRequest",
    "User-Agent": "biomedical-mcp/0.4",
}

# Refresh CSVs at most once per day
REFRESH_INTERVAL_SECONDS = 86400


class ClinGen:
    """ClinGen data client backed by cached CSV downloads."""

    def __init__(self, cache: Cache):
        self.cache = cache
        self._validity_data: list[dict] | None = None
        self._dosage_data: list[dict] | None = None
        self._last_refresh: float = 0
        self._http = httpx.Client(timeout=60.0, headers=DOWNLOAD_HEADERS)

    def _ensure_loaded(self) -> None:
        """Load or refresh CSV data from cache or remote."""
        now = time.time()
        if self._validity_data is not None and (now - self._last_refresh) < REFRESH_INTERVAL_SECONDS:
            return

        # Try cache first
        cached_v = self.cache.get("clingen:validity_csv", max_age_days=1)
        cached_d = self.cache.get("clingen:dosage_csv", max_age_days=1)

        if cached_v and cached_d:
            self._validity_data = cached_v
            self._dosage_data = cached_d
            self._last_refresh = now
            return

        # Download fresh
        self._validity_data = self._download_validity()
        self._dosage_data = self._download_dosage()
        self._last_refresh = now

        if self._validity_data:
            self.cache.set("clingen:validity_csv", self._validity_data)
        if self._dosage_data:
            self.cache.set("clingen:dosage_csv", self._dosage_data)

    def _download_validity(self) -> list[dict]:
        """Download and parse gene-disease validity CSV."""
        try:
            resp = self._http.get(VALIDITY_URL)
            resp.raise_for_status()
        except Exception as e:
            log.warning("ClinGen validity download failed: %s", e)
            return []

        return _parse_validity_csv(resp.text)

    def _download_dosage(self) -> list[dict]:
        """Download and parse dosage sensitivity CSV."""
        try:
            resp = self._http.get(DOSAGE_URL)
            resp.raise_for_status()
        except Exception as e:
            log.warning("ClinGen dosage download failed: %s", e)
            return []

        return _parse_dosage_csv(resp.text)

    def gene_validity(self, gene_symbol: str) -> dict:
        """Get gene-disease validity curations for a gene."""
        self._ensure_loaded()
        matches = [
            r for r in (self._validity_data or [])
            if r["gene_symbol"].upper() == gene_symbol.upper()
        ]

        return {
            "gene_symbol": gene_symbol,
            "count": len(matches),
            "curations": matches,
            "provenance": make_provenance(
                "clingen_gene_validity",
                evidence_grade="A1",
                cache_hit=bool(self._validity_data),
            ),
        }

    def disease_validity(self, disease_label: str) -> dict:
        """Search gene-disease validity curations by disease name (substring match)."""
        self._ensure_loaded()
        query = disease_label.lower()
        matches = [
            r for r in (self._validity_data or [])
            if query in r["disease_label"].lower()
        ]

        return {
            "disease_query": disease_label,
            "count": len(matches),
            "curations": matches,
            "provenance": make_provenance(
                "clingen_gene_validity",
                evidence_grade="A1",
                cache_hit=bool(self._validity_data),
            ),
        }

    def classification_summary(self, classification: str = "Definitive") -> dict:
        """Get all gene-disease pairs with a specific classification level."""
        self._ensure_loaded()
        matches = [
            r for r in (self._validity_data or [])
            if r["classification"].lower() == classification.lower()
        ]

        return {
            "classification": classification,
            "count": len(matches),
            "curations": matches,
            "provenance": make_provenance("clingen_gene_validity", evidence_grade="A1"),
        }

    def gene_dosage(self, gene_symbol: str) -> dict:
        """Get dosage sensitivity (haploinsufficiency/triplosensitivity) for a gene."""
        self._ensure_loaded()
        matches = [
            r for r in (self._dosage_data or [])
            if r["gene_symbol"].upper() == gene_symbol.upper()
        ]

        return {
            "gene_symbol": gene_symbol,
            "count": len(matches),
            "dosage": matches,
            "provenance": make_provenance(
                "clingen_dosage",
                evidence_grade="A1",
                cache_hit=bool(self._dosage_data),
            ),
        }


def _parse_validity_csv(text: str) -> list[dict]:
    """Parse ClinGen gene-disease validity CSV (has header rows to skip)."""
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)

    # Find the header row with "GENE SYMBOL"
    header_idx = None
    for i, row in enumerate(rows):
        if row and "GENE SYMBOL" in row:
            header_idx = i
            break

    if header_idx is None:
        log.warning("Could not find header in ClinGen validity CSV")
        return []

    # Skip separator row (+++++...)
    data_start = header_idx + 1
    if data_start < len(rows) and rows[data_start][0].startswith("+"):
        data_start += 1

    results = []
    for row in rows[data_start:]:
        if len(row) < 7 or not row[0].strip():
            continue
        results.append({
            "gene_symbol": row[0].strip(),
            "gene_id": row[1].strip(),
            "disease_label": row[2].strip(),
            "disease_id": row[3].strip(),
            "moi": row[4].strip(),
            "sop": row[5].strip(),
            "classification": row[6].strip(),
            "report_url": row[7].strip() if len(row) > 7 else "",
            "classification_date": row[8].strip() if len(row) > 8 else "",
            "expert_panel": row[9].strip() if len(row) > 9 else "",
        })

    return results


def _parse_dosage_csv(text: str) -> list[dict]:
    """Parse ClinGen dosage sensitivity CSV."""
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)

    # Find the header row with "GENE SYMBOL"
    header_idx = None
    for i, row in enumerate(rows):
        if row and "GENE SYMBOL" in row:
            header_idx = i
            break

    if header_idx is None:
        log.warning("Could not find header in ClinGen dosage CSV")
        return []

    # Skip separator row if present
    data_start = header_idx + 1
    if data_start < len(rows) and rows[data_start] and rows[data_start][0].startswith("+"):
        data_start += 1

    results = []
    for row in rows[data_start:]:
        if len(row) < 4 or not row[0].strip():
            continue
        results.append({
            "gene_symbol": row[0].strip(),
            "hgnc_id": row[1].strip(),
            "haploinsufficiency": row[2].strip(),
            "triplosensitivity": row[3].strip(),
            "report_url": row[4].strip() if len(row) > 4 else "",
            "date": row[5].strip() if len(row) > 5 else "",
        })

    return results
