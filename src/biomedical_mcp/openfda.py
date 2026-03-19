"""OpenFDA REST API client — adverse events, drug labels, and recalls."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

OPENFDA_BASE = "https://api.fda.gov"

LABEL_SECTIONS = [
    "boxed_warning", "warnings", "dosage_and_administration",
    "drug_interactions", "adverse_reactions", "indications_and_usage",
    "contraindications", "clinical_pharmacology", "pharmacokinetics",
    "use_in_specific_populations",
]


def _is_retryable(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 500, 502, 503, 504)
    return isinstance(exc, (httpx.ConnectError, httpx.ReadTimeout))


class OpenFDA:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(
            base_url=OPENFDA_BASE,
            timeout=30,
            headers={"User-Agent": "biomedical-mcp/0.1"},
        )

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=2, max=15), retry=retry_if_exception(_is_retryable))
    def _get(self, path: str, params: dict | None = None) -> dict:
        resp = self.client.get(path, params=params)
        resp.raise_for_status()
        return resp.json()

    def _cached(self, prefix: str, path: str, params: dict | None = None, max_age_days: int = 7) -> dict:
        raw = f"{path}:{sorted(params.items()) if params else ''}"
        key = f"openfda:{prefix}:{hashlib.md5(raw.encode()).hexdigest()}"
        cached = self.cache.get(key, max_age_days=max_age_days)
        if cached is not None:
            return cached
        result = self._get(path, params)
        self.cache.set(key, result)
        return result

    def adverse_events(self, drug_name: str, reaction: str | None = None, serious: bool | None = None, limit: int = 25) -> dict:
        """Get adverse event reports for a drug."""
        search = f'patient.drug.openfda.generic_name:"{drug_name}"'
        if reaction:
            search += f' AND patient.reaction.reactionmeddrapt:"{reaction}"'
        if serious is not None:
            search += f" AND serious:{1 if serious else 2}"
        params = {"search": search, "limit": min(limit, 100)}
        data = self._cached("ae", "/drug/event.json", params)
        results = data.get("results", [])
        return {
            "drug": drug_name,
            "total": data.get("meta", {}).get("results", {}).get("total", 0),
            "events": [self._extract_event(e) for e in results],
        }

    def drug_label(self, drug_name: str, sections: list[str] | None = None) -> dict:
        """Get drug labeling information."""
        params = {"search": f'openfda.generic_name:"{drug_name}"', "limit": 1}
        data = self._cached("label", "/drug/label.json", params)
        results = data.get("results", [])
        if not results:
            return {"error": f"No label found for: {drug_name}"}
        label = results[0]
        requested = sections or LABEL_SECTIONS
        extracted = {}
        for s in requested:
            val = label.get(s)
            if val:
                # Label sections are arrays of strings
                extracted[s] = val[0][:3000] if isinstance(val, list) and val else str(val)[:3000]
        openfda = label.get("openfda", {})
        return {
            "drug": drug_name,
            "brand_names": openfda.get("brand_name", []),
            "manufacturer": openfda.get("manufacturer_name", []),
            "route": openfda.get("route", []),
            "product_type": openfda.get("product_type", []),
            "sections": extracted,
        }

    def recalls(self, drug_name: str, classification: str | None = None, limit: int = 25) -> dict:
        """Get drug recall enforcement reports."""
        search = f'openfda.generic_name:"{drug_name}"'
        if classification:
            search += f' AND classification:"{classification}"'
        params = {"search": search, "limit": min(limit, 100)}
        data = self._cached("recall", "/drug/enforcement.json", params)
        results = data.get("results", [])
        return {
            "drug": drug_name,
            "total": data.get("meta", {}).get("results", {}).get("total", 0),
            "recalls": [
                {
                    "recall_number": r.get("recall_number"),
                    "reason": r.get("reason_for_recall"),
                    "status": r.get("status"),
                    "classification": r.get("classification"),
                    "product_description": r.get("product_description", "")[:500],
                    "recall_initiation_date": r.get("recall_initiation_date"),
                    "voluntary": r.get("voluntary_mandated"),
                }
                for r in results
            ],
        }

    def _extract_event(self, event: dict) -> dict:
        """Extract key fields from an adverse event report."""
        reactions = [r.get("reactionmeddrapt") for r in (event.get("patient", {}).get("reaction") or [])]
        drugs = []
        for d in (event.get("patient", {}).get("drug") or [])[:5]:
            drugs.append({
                "name": d.get("medicinalproduct"),
                "indication": d.get("drugindication"),
                "characterization": d.get("drugcharacterization"),
            })
        return {
            "serious": event.get("serious"),
            "seriousness_death": event.get("seriousnessdeath"),
            "seriousness_hospitalization": event.get("seriousnesshospitalization"),
            "reactions": reactions[:10],
            "drugs": drugs,
            "outcome": event.get("patient", {}).get("patientonsetage"),
            "receive_date": event.get("receivedate"),
        }
