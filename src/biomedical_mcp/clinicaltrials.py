"""ClinicalTrials.gov v2 API client."""

import hashlib
import logging
from collections import Counter

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

CT_BASE = "https://clinicaltrials.gov/api/v2"


def _is_retryable(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 500, 502, 503, 504)
    return isinstance(exc, (httpx.ConnectError, httpx.ReadTimeout))


class ClinicalTrials:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(
            base_url=CT_BASE,
            timeout=30,
            headers={"User-Agent": "biomedical-mcp/0.1"},
        )

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=2, max=15), retry=retry_if_exception(_is_retryable))
    def _get(self, path: str, params: dict | None = None) -> dict:
        resp = self.client.get(path, params=params)
        resp.raise_for_status()
        return resp.json()

    def _cached(self, prefix: str, path: str, params: dict, max_age_days: int = 7) -> dict:
        raw = f"{path}:{sorted(params.items())}"
        key = f"ct:{prefix}:{hashlib.md5(raw.encode()).hexdigest()}"
        cached = self.cache.get(key, max_age_days=max_age_days)
        if cached is not None:
            return cached
        result = self._get(path, params)
        self.cache.set(key, result)
        return result

    def _extract_study(self, study: dict) -> dict:
        ps = study.get("protocolSection", {})
        ident = ps.get("identificationModule", {})
        status = ps.get("statusModule", {})
        design = ps.get("designModule", {})
        sponsor = ps.get("sponsorCollaboratorsModule", {})
        conditions = ps.get("conditionsModule", {})
        interventions = ps.get("armsInterventionsModule", {})
        enrollment = design.get("enrollmentInfo", {})
        return {
            "nct_id": ident.get("nctId"),
            "title": ident.get("briefTitle"),
            "status": status.get("overallStatus"),
            "phases": design.get("phases", []),
            "study_type": design.get("studyType"),
            "enrollment": enrollment.get("count"),
            "enrollment_type": enrollment.get("type"),
            "conditions": conditions.get("conditions", []),
            "interventions": [
                {"type": i.get("type"), "name": i.get("name")}
                for i in interventions.get("interventions", [])
            ],
            "sponsor": (sponsor.get("leadSponsor") or {}).get("name"),
            "start_date": (status.get("startDateStruct") or {}).get("date"),
            "completion_date": (status.get("primaryCompletionDateStruct") or {}).get("date"),
        }

    def search(
        self,
        condition: str | None = None,
        intervention: str | None = None,
        gene: str | None = None,
        status: str | None = None,
        limit: int = 20,
    ) -> dict:
        params = {"pageSize": min(limit, 100)}
        if condition:
            params["query.cond"] = condition
        if intervention:
            params["query.intr"] = intervention
        if gene:
            params["query.term"] = gene
        if status:
            params["filter.overallStatus"] = status

        data = self._cached("search", "/studies", params)
        studies = [self._extract_study(s) for s in data.get("studies", [])]
        return {
            "total": data.get("totalCount", len(studies)),
            "returned": len(studies),
            "studies": studies,
        }

    def trial_detail(self, nct_id: str) -> dict:
        data = self._cached("detail", f"/studies/{nct_id}", {})
        ps = data.get("protocolSection", {})
        base = self._extract_study(data)
        # Add detailed fields
        eligibility = ps.get("eligibilityModule", {})
        outcomes = ps.get("outcomesModule", {})
        design = ps.get("designModule", {})
        base.update({
            "official_title": ps.get("identificationModule", {}).get("officialTitle"),
            "description": ps.get("descriptionModule", {}).get("briefSummary"),
            "eligibility": {
                "criteria": eligibility.get("eligibilityCriteria"),
                "sex": eligibility.get("sex"),
                "min_age": eligibility.get("minimumAge"),
                "max_age": eligibility.get("maximumAge"),
                "healthy_volunteers": eligibility.get("healthyVolunteers"),
            },
            "design_info": design.get("designInfo", {}),
            "primary_outcomes": [
                {"measure": o.get("measure"), "timeframe": o.get("timeFrame")}
                for o in outcomes.get("primaryOutcomes", [])
            ],
            "secondary_outcomes": [
                {"measure": o.get("measure"), "timeframe": o.get("timeFrame")}
                for o in outcomes.get("secondaryOutcomes", [])
            ],
        })
        return base

    def stats(self, condition: str) -> dict:
        data = self._cached("stats", "/studies", {"query.cond": condition, "pageSize": 100})
        studies = data.get("studies", [])
        phase_counts: Counter = Counter()
        status_counts: Counter = Counter()
        for s in studies:
            ps = s.get("protocolSection", {})
            for phase in ps.get("designModule", {}).get("phases", ["N/A"]):
                phase_counts[phase] += 1
            status_counts[ps.get("statusModule", {}).get("overallStatus", "UNKNOWN")] += 1
        return {
            "condition": condition,
            "total": data.get("totalCount", len(studies)),
            "sampled": len(studies),
            "by_phase": dict(phase_counts.most_common()),
            "by_status": dict(status_counts.most_common()),
        }
