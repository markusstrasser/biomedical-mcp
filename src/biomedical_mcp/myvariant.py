"""MyVariant.info REST API client — variant annotations from 26+ field groups."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

MYVARIANT_BASE = "https://myvariant.info/v1"

DEFAULT_FIELDS = (
    "clinvar,gnomad_exome,gnomad_genome,cadd,dbnsfp.genename,"
    "dbnsfp.sift,dbnsfp.polyphen2,dbsnp,civic,snpedia"
)


def _is_retryable(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 500, 502, 503, 504)
    return isinstance(exc, (httpx.ConnectError, httpx.ReadTimeout))


class MyVariant:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(
            base_url=MYVARIANT_BASE,
            timeout=30,
            headers={"User-Agent": "biomedical-mcp/0.1"},
        )

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=2, max=15), retry=retry_if_exception(_is_retryable))
    def _get(self, path: str, params: dict | None = None) -> dict:
        resp = self.client.get(path, params=params)
        resp.raise_for_status()
        return resp.json()

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=2, max=15), retry=retry_if_exception(_is_retryable))
    def _post(self, path: str, data: dict) -> list | dict:
        resp = self.client.post(path, data=data)
        resp.raise_for_status()
        return resp.json()

    def _cached(self, prefix: str, path: str, params: dict | None = None, max_age_days: int = 30) -> dict:
        raw = f"{path}:{sorted(params.items()) if params else ''}"
        key = f"myvariant:{prefix}:{hashlib.md5(raw.encode()).hexdigest()}"
        cached = self.cache.get(key, max_age_days=max_age_days)
        if cached is not None:
            return cached
        result = self._get(path, params)
        self.cache.set(key, result)
        return result

    @staticmethod
    def _normalize_variant_id(variant_id: str) -> str:
        """Normalize variant ID to MyVariant.info HGVS format.

        Accepts: rsID, HGVS (chr:g.posRef>Alt), or VCF-style (chr:pos:ref>alt).
        """
        if variant_id.startswith("rs"):
            return variant_id
        if ":g." in variant_id:
            return variant_id
        # VCF-style: chr10:129867906:G>A → chr10:g.129867906G>A
        parts = variant_id.split(":")
        if len(parts) == 3:
            chrom, pos, change = parts
            return f"{chrom}:g.{pos}{change}"
        return variant_id

    def lookup(self, variant_id: str, fields: str | None = None) -> dict | list[dict]:
        """Look up a single variant by rsID, HGVS, or chrX:g.posA>B.

        Also accepts VCF-style chr:pos:ref>alt (auto-converted to HGVS).
        rsIDs with multiple alt alleles return a list of annotations.
        """
        normalized = self._normalize_variant_id(variant_id)
        f = fields or DEFAULT_FIELDS
        try:
            data = self._cached("var", f"/variant/{normalized}", {"fields": f})
        except httpx.HTTPStatusError as exc:
            if exc.response.status_code == 404:
                return {"variant_id": normalized, "error": "not_found", "query": variant_id}
            raise
        if isinstance(data, list):
            return [self._extract_annotations(d) for d in data]
        return self._extract_annotations(data)

    def clinvar_search(self, gene_symbol: str, significance: str | None = None, limit: int = 25) -> list[dict]:
        """Search ClinVar variants for a gene, optionally filtered by significance."""
        q = f"clinvar.gene.symbol:{gene_symbol}"
        if significance:
            q += f" AND clinvar.rcv.clinical_significance:{significance}"
        params = {
            "q": q,
            "fields": "clinvar,gnomad_exome,cadd,_id",
            "size": min(limit, 100),
        }
        data = self._cached("cvar", "/query", params)
        hits = data.get("hits", [])
        return [self._extract_clinvar_hit(h) for h in hits]

    def batch(self, variant_ids: list[str], fields: str | None = None) -> list[dict]:
        """Batch lookup up to 100 variants via POST."""
        normalized = [self._normalize_variant_id(vid) for vid in variant_ids[:100]]
        f = fields or DEFAULT_FIELDS
        # Batch results aren't cached individually — cache the whole batch
        raw = f"batch:{','.join(sorted(normalized))}:{f}"
        key = f"myvariant:batch:{hashlib.md5(raw.encode()).hexdigest()}"
        cached = self.cache.get(key, max_age_days=30)
        if cached is not None:
            return cached
        result = self._post("/variant", data={"ids": ",".join(normalized), "fields": f})
        if isinstance(result, dict):
            result = [result]
        extracted = []
        for record in result:
            if record.get("notfound"):
                extracted.append({"variant_id": record.get("_id", record.get("query")), "error": "not_found"})
            else:
                extracted.append(self._extract_annotations(record))
        self.cache.set(key, extracted)
        return extracted

    def _extract_annotations(self, data: dict) -> dict:
        """Extract structured annotations from raw MyVariant response."""
        result = {"variant_id": data.get("_id")}

        # ClinVar
        cv = data.get("clinvar")
        if cv:
            rcv = cv.get("rcv", [])
            if isinstance(rcv, dict):
                rcv = [rcv]
            result["clinvar"] = {
                "variant_id": cv.get("variant_id"),
                "gene": cv.get("gene", {}).get("symbol") if isinstance(cv.get("gene"), dict) else None,
                "significance": [r.get("clinical_significance") for r in rcv[:5]],
                "conditions": [r.get("conditions", {}).get("name") if isinstance(r.get("conditions"), dict) else None for r in rcv[:5]],
            }

        # gnomAD
        for src in ("gnomad_exome", "gnomad_genome"):
            gn = data.get(src)
            if gn:
                af = gn.get("af", {})
                result[src] = {
                    "allele_freq": af if isinstance(af, (int, float)) else af.get("af") if isinstance(af, dict) else None,
                    "filter": gn.get("filter"),
                }

        # CADD
        cadd = data.get("cadd")
        if cadd:
            result["cadd"] = {
                "phred": cadd.get("phred"),
                "raw_score": cadd.get("rawscore"),
            }

        # dbNSFP predictions
        dbnsfp = data.get("dbnsfp")
        if dbnsfp:
            sift = dbnsfp.get("sift", {})
            pp2 = dbnsfp.get("polyphen2", {})
            result["predictions"] = {
                "gene": dbnsfp.get("genename"),
                "sift_pred": sift.get("pred") if isinstance(sift, dict) else None,
                "sift_score": sift.get("score") if isinstance(sift, dict) else None,
                "polyphen2_pred": pp2.get("hdiv", {}).get("pred") if isinstance(pp2, dict) else None,
                "polyphen2_score": pp2.get("hdiv", {}).get("score") if isinstance(pp2, dict) else None,
            }

        # dbSNP
        dbsnp = data.get("dbsnp")
        if dbsnp:
            result["dbsnp"] = {
                "rsid": dbsnp.get("rsid"),
                "ref": dbsnp.get("ref"),
                "alt": dbsnp.get("alt"),
                "chrom": dbsnp.get("chrom"),
                "position": dbsnp.get("hg19", {}).get("start") if isinstance(dbsnp.get("hg19"), dict) else None,
            }

        # CIViC
        civic = data.get("civic")
        if civic:
            result["civic"] = {
                "name": civic.get("name"),
                "description": civic.get("description"),
                "evidence_items": civic.get("evidence_items") if isinstance(civic.get("evidence_items"), int) else None,
            }

        # SNPedia
        snpedia = data.get("snpedia")
        if snpedia:
            result["snpedia"] = {"text": snpedia.get("text", "")[:500]}

        return result

    def _extract_clinvar_hit(self, hit: dict) -> dict:
        """Extract a ClinVar search result."""
        result = {"variant_id": hit.get("_id")}
        cv = hit.get("clinvar", {})
        if cv:
            rcv = cv.get("rcv", [])
            if isinstance(rcv, dict):
                rcv = [rcv]
            result["significance"] = [r.get("clinical_significance") for r in rcv[:3]]
            gene = cv.get("gene", {})
            result["gene"] = gene.get("symbol") if isinstance(gene, dict) else None

        gn = hit.get("gnomad_exome", {})
        if gn:
            af = gn.get("af", {})
            result["gnomad_af"] = af if isinstance(af, (int, float)) else af.get("af") if isinstance(af, dict) else None

        cadd = hit.get("cadd", {})
        if cadd:
            result["cadd_phred"] = cadd.get("phred")

        return result
