"""Open Targets Platform GraphQL client."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# --- GraphQL queries ---

SEARCH_QUERY = """
query($q: String!, $size: Int!, $entities: [String!]) {
  search(queryString: $q, entityNames: $entities, page: {size: $size, index: 0}) {
    total
    hits { id name entity description }
  }
}
"""

TARGET_INFO_QUERY = """
query($id: String!) {
  target(ensemblId: $id) {
    id
    approvedSymbol
    approvedName
    biotype
    functionDescriptions
    tractability { label modality value }
    safetyLiabilities { event datasource }
    knownDrugs: drugAndClinicalCandidates {
      count
      rows {
        drug { id name drugType maximumClinicalTrialPhase: maximumClinicalStage }
        diseases { disease { id name } }
        phase: maxClinicalStage
      }
    }
  }
}
"""

DISEASE_ASSOCIATIONS_QUERY = """
query($id: String!, $size: Int!) {
  target(ensemblId: $id) {
    id
    approvedSymbol
    associatedDiseases(page: {size: $size, index: 0}) {
      count
      rows {
        disease { id name }
        score
        datasourceScores { id score }
      }
    }
  }
}
"""

DISEASE_TARGETS_QUERY = """
query($id: String!, $size: Int!) {
  disease(efoId: $id) {
    id
    name
    associatedTargets(page: {size: $size, index: 0}) {
      count
      rows {
        target { id approvedSymbol approvedName }
        score
        datasourceScores { id score }
      }
    }
  }
}
"""

PHARMACOGENETICS_QUERY = """
query($id: String!) {
  target(ensemblId: $id) {
    id
    approvedSymbol
    pharmacogenomics {
      variantRsId
      genotype
      genotypeAnnotationText
      drugs { drugFromSource drugId }
      phenotypeText
      pgxCategory
      evidenceLevel
      isDirectTarget
      datasourceId
    }
  }
}
"""

DRUG_INFO_QUERY = """
query($id: String!) {
  drug(chemblId: $id) {
    id
    name
    drugType
    maximumClinicalTrialPhase: maximumClinicalStage
    drugWarnings { warningType }
    description
    synonyms { label }
    mechanismsOfAction {
      rows {
        mechanismOfAction
        targets { id approvedSymbol }
      }
    }
    indications {
      count
      rows {
        disease { id name }
        maxPhaseForIndication: maxClinicalStage
      }
    }
  }
}
"""


# Live maxClinicalStage values (2026-10-08) are strings, most advanced first; unlisted values sort last.
STAGE_ORDER = ["APPROVAL", "PHASE_4", "PHASE_3", "PHASE_2_3", "PHASE_2", "PHASE_1_2", "PHASE_1", "WITHDRAWAL", "UNKNOWN"]


def _stage_sort_key(row: dict) -> tuple:
    stage = row.get("phase")
    rank = STAGE_ORDER.index(stage) if stage in STAGE_ORDER else len(STAGE_ORDER)
    return (rank, ((row.get("drug") or {}).get("name") or "").upper())


def _is_retryable(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 500, 502, 503, 504)
    return isinstance(exc, (httpx.ConnectError, httpx.ReadTimeout))


class OpenTargets:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(
            timeout=30,
            headers={"Content-Type": "application/json", "User-Agent": "biomedical-mcp/0.1"},
        )

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=2, max=15), retry=retry_if_exception(_is_retryable))
    def _query(self, query: str, variables: dict) -> dict:
        resp = self.client.post(OT_URL, json={"query": query, "variables": variables})
        resp.raise_for_status()
        data = resp.json()
        if "errors" in data:
            raise ValueError(f"GraphQL errors: {data['errors']}")
        return data["data"]

    def _cached(self, prefix: str, query: str, variables: dict, max_age_days: int = 30) -> dict:
        key = f"ot:{prefix}:{hashlib.md5(str(sorted(variables.items())).encode()).hexdigest()}"
        cached = self.cache.get(key, max_age_days=max_age_days)
        if cached is not None:
            return cached
        result = self._query(query, variables)
        self.cache.set(key, result)
        return result

    def resolve_target(self, gene_symbol: str | None = None, ensembl_id: str | None = None) -> str | None:
        if ensembl_id and ensembl_id.startswith("ENSG"):
            return ensembl_id
        if not gene_symbol:
            return None
        hits = self.search(gene_symbol, entity_type="target", limit=5)
        # Prefer exact symbol match
        for h in hits:
            if h.get("name", "").upper() == gene_symbol.upper():
                return h["id"]
        return hits[0]["id"] if hits else None

    def resolve_drug(self, drug_name: str | None = None, chembl_id: str | None = None) -> str | None:
        if chembl_id and chembl_id.startswith("CHEMBL"):
            return chembl_id
        if not drug_name:
            return None
        hits = self.search(drug_name, entity_type="drug", limit=5)
        return hits[0]["id"] if hits else None

    def search(self, query: str, entity_type: str | None = None, limit: int = 10) -> list[dict]:
        entities = [entity_type] if entity_type else None
        data = self._cached("search", SEARCH_QUERY, {"q": query, "size": min(limit, 50), "entities": entities})
        search_data = data.get("search", {})
        return search_data.get("hits", [])

    def target_info(self, ensembl_id: str) -> dict:
        """Target details with known drugs, most advanced stage first (count is the full total).

        Keeps every APPROVAL row, then fills with the next stages up to 10 rows in total.

        knownDrugs now comes from drugAndClinicalCandidates: one row per drug with a
        `diseases` list (was one row per drug-disease pair with `disease`), and
        `phase`/`maximumClinicalTrialPhase` are clinical-stage strings (e.g. "APPROVAL").
        The per-row trial `status` has no equivalent and is dropped.
        """
        data = self._cached("target", TARGET_INFO_QUERY, {"id": ensembl_id})
        target = data.get("target")
        known = (target or {}).get("knownDrugs")
        if known:
            for row in known["rows"]:
                row["diseases"] = [d for d in row.get("diseases") or [] if d.get("disease")]
            ranked = sorted(known["rows"], key=_stage_sort_key)
            approved = sum(1 for r in ranked if r.get("phase") == "APPROVAL")
            known["rows"] = ranked[: max(10, approved)]
        return target

    def disease_associations(self, ensembl_id: str, limit: int = 25) -> dict:
        data = self._cached("assoc", DISEASE_ASSOCIATIONS_QUERY, {"id": ensembl_id, "size": min(limit, 100)})
        target = data.get("target") or {}
        assoc = target.get("associatedDiseases") or {}
        return {
            "target": target.get("approvedSymbol"),
            "total": assoc.get("count", 0),
            "associations": [
                {
                    "disease_id": r["disease"]["id"],
                    "disease_name": r["disease"]["name"],
                    "score": r["score"],
                    "sources": {s["id"]: s["score"] for s in r.get("datasourceScores", []) if s["score"] > 0},
                }
                for r in assoc.get("rows", [])
            ],
        }

    def disease_targets(self, disease_id: str, limit: int = 25) -> dict:
        data = self._cached("dtargets", DISEASE_TARGETS_QUERY, {"id": disease_id, "size": min(limit, 100)})
        disease = data.get("disease") or {}
        assoc = disease.get("associatedTargets") or {}
        return {
            "disease": disease.get("name"),
            "disease_id": disease.get("id"),
            "total": assoc.get("count", 0),
            "targets": [
                {
                    "ensembl_id": r["target"]["id"],
                    "symbol": r["target"]["approvedSymbol"],
                    "name": r["target"]["approvedName"],
                    "score": r["score"],
                }
                for r in assoc.get("rows", [])
            ],
        }

    def pharmacogenetics(self, ensembl_id: str, limit: int = 50) -> dict:
        data = self._cached("pgx", PHARMACOGENETICS_QUERY, {"id": ensembl_id})
        target = data.get("target") or {}
        all_entries = target.get("pharmacogenomics") or []
        entries = all_entries[:limit]
        return {
            "target": target.get("approvedSymbol"),
            "entries": entries,
            "total_available": len(all_entries),
            "truncated": len(all_entries) > limit,
        }

    def drug_info(self, chembl_id: str) -> dict:
        """Drug details from Open Targets.

        `maximumClinicalTrialPhase` and indication `maxPhaseForIndication` are now
        clinical-stage strings (e.g. "APPROVAL"); `hasBeenWithdrawn` is derived from a
        "Withdrawn" drug warning; `synonyms` are flattened to labels.
        `linkedDiseases` and `linkedTargets` have no equivalent and are dropped.
        """
        data = self._cached("drug", DRUG_INFO_QUERY, {"id": chembl_id})
        drug = data.get("drug")
        if drug:
            warnings = drug.pop("drugWarnings", None) or []
            drug["hasBeenWithdrawn"] = any(w.get("warningType") == "Withdrawn" for w in warnings)
            drug["synonyms"] = sorted({s["label"] for s in drug.get("synonyms") or []})
        return drug
