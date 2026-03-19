"""ChEMBL REST API client."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


def _is_retryable(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 500, 502, 503, 504)
    return isinstance(exc, (httpx.ConnectError, httpx.ReadTimeout))


class ChEMBL:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(
            base_url=CHEMBL_BASE,
            timeout=30,
            headers={"Accept": "application/json", "User-Agent": "biomedical-mcp/0.1"},
        )

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=2, max=15), retry=retry_if_exception(_is_retryable))
    def _get(self, path: str, params: dict | None = None) -> dict:
        resp = self.client.get(path, params=params)
        resp.raise_for_status()
        return resp.json()

    def _cached(self, prefix: str, path: str, params: dict | None = None, max_age_days: int = 30) -> dict:
        raw = f"{path}:{sorted(params.items()) if params else ''}"
        key = f"chembl:{prefix}:{hashlib.md5(raw.encode()).hexdigest()}"
        cached = self.cache.get(key, max_age_days=max_age_days)
        if cached is not None:
            return cached
        result = self._get(path, params)
        self.cache.set(key, result)
        return result

    def resolve_compound(self, name: str | None = None, chembl_id: str | None = None) -> str | None:
        if chembl_id and chembl_id.startswith("CHEMBL"):
            return chembl_id
        if not name:
            return None
        data = self._cached("msearch", f"/molecule/search.json", {"q": name})
        molecules = data.get("molecules", [])
        if not molecules:
            return None
        # Prefer exact name match
        for m in molecules:
            if m.get("pref_name", "").upper() == name.upper():
                return m["molecule_chembl_id"]
        return molecules[0]["molecule_chembl_id"]

    def compound(self, chembl_id: str) -> dict:
        data = self._cached("mol", f"/molecule/{chembl_id}.json")
        props = data.get("molecule_properties", {}) or {}
        return {
            "chembl_id": data.get("molecule_chembl_id"),
            "name": data.get("pref_name"),
            "type": data.get("molecule_type"),
            "max_phase": data.get("max_phase"),
            "molecular_weight": props.get("full_mwt"),
            "alogp": props.get("alogp"),
            "hba": props.get("hba"),
            "hbd": props.get("hbd"),
            "psa": props.get("psa"),
            "ro5_violations": props.get("num_ro5_violations"),
            "atc_classifications": data.get("atc_classifications", []),
            "synonyms": [s["molecule_synonym"] for s in (data.get("molecule_synonyms") or [])],
            "first_approval": data.get("first_approval"),
        }

    def target(self, gene_symbol: str | None = None, chembl_id: str | None = None, uniprot_id: str | None = None) -> dict | None:
        if chembl_id and chembl_id.startswith("CHEMBL"):
            data = self._cached("tgt", f"/target/{chembl_id}.json")
        elif uniprot_id:
            data = self._cached("tgt_uni", "/target.json", {"target_components__accession": uniprot_id})
            targets = data.get("targets", [])
            if not targets:
                return None
            data = targets[0]
        elif gene_symbol:
            search = self._cached("tgt_search", "/target/search.json", {"q": gene_symbol})
            targets = search.get("targets", [])
            if not targets:
                return None
            data = targets[0]
        else:
            return None
        components = data.get("target_components", []) or []
        return {
            "chembl_id": data.get("target_chembl_id"),
            "name": data.get("pref_name"),
            "type": data.get("target_type"),
            "organism": data.get("organism"),
            "accessions": [c.get("accession") for c in components if c.get("accession")],
            "cross_references": [
                {"source": x.get("xref_src"), "id": x.get("xref_id")}
                for c in components
                for x in (c.get("target_component_xrefs") or [])
            ][:20],
        }

    def mechanism(self, chembl_id: str) -> list[dict]:
        data = self._cached("mech", "/mechanism.json", {"molecule_chembl_id": chembl_id})
        return [
            {
                "mechanism": m.get("mechanism_of_action"),
                "action_type": m.get("action_type"),
                "target_name": m.get("target_pref_name"),
                "target_chembl_id": m.get("target_chembl_id"),
                "references": [{"type": r.get("ref_type"), "id": r.get("ref_id")} for r in (m.get("mechanism_refs") or [])],
            }
            for m in data.get("mechanisms", [])
        ]

    def bioactivity(self, target_chembl_id: str, limit: int = 25) -> list[dict]:
        data = self._cached("act", "/activity.json", {"target_chembl_id": target_chembl_id, "limit": min(limit, 100)})
        return [
            {
                "molecule_chembl_id": a.get("molecule_chembl_id"),
                "molecule_name": a.get("molecule_pref_name"),
                "type": a.get("standard_type"),
                "value": a.get("standard_value"),
                "units": a.get("standard_units"),
                "relation": a.get("standard_relation"),
                "assay_type": a.get("assay_type"),
            }
            for a in data.get("activities", [])
        ]

    def drug_indications(self, chembl_id: str, limit: int = 25) -> list[dict]:
        data = self._cached("ind", "/drug_indication.json", {"molecule_chembl_id": chembl_id, "limit": min(limit, 100)})
        return [
            {
                "mesh_id": i.get("mesh_id"),
                "mesh_heading": i.get("mesh_heading"),
                "efo_id": i.get("efo_id"),
                "efo_term": i.get("efo_term"),
                "max_phase": i.get("max_phase_for_ind"),
                "references": [{"type": r.get("ref_type"), "id": r.get("ref_id")} for r in (i.get("indication_refs") or [])],
            }
            for i in data.get("drug_indications", [])
        ]
