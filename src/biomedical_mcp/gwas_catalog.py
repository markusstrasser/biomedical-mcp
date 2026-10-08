"""GWAS Catalog (EBI) REST API client — variant-trait associations."""

from __future__ import annotations

import logging

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class GWASCatalog(BaseClient):
    domain_name = "gwas_catalog"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    def variant_associations(self, rsid: str, limit: int = 25) -> dict:
        """GWAS associations for a variant (rsID), via v2 ``/associations?rs_id=``.

        Rows come strongest first (server-side ``sort=p_value&direction=asc``);
        ``association_count`` is the total match count, the list is capped at ``limit``.
        An unknown rsID returns an empty page in v2 (no 404).
        """
        key = self._cache_key("var_assoc_v2", rsid)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        try:
            data = self._get("/associations", params={"rs_id": rsid, "size": min(limit, 200),
                                     "sort": "p_value", "direction": "asc"})
        except Exception as exc:
            log.warning("GWAS Catalog API error for %s: %s", rsid, exc)
            return {"rsid": rsid, "association_count": 0, "associations": [],
                    "error": "api_error", "note": f"GWAS Catalog API error: {type(exc).__name__}"}
        associations = self._extract_associations(data, limit)
        result = {
            "rsid": rsid,
            "association_count": self._total(data, associations),
            "associations": associations,
        }
        if not associations:
            result["note"] = "Variant not found in GWAS Catalog"
        self.cache.set(key, result)
        return result

    def gene_associations(self, gene_symbol: str, limit: int = 25) -> dict:
        """GWAS associations for a gene, via v2 ``/associations?mapped_gene=``.

        v2 filters on Ensembl-mapped genes; v1 used author-reported genes. Strongest
        first; ``association_count`` is the total, the list is capped at ``limit``.
        """
        key = self._cache_key("gene_assoc_v2", gene_symbol)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        try:
            data = self._get("/associations",
                             params={"mapped_gene": gene_symbol, "size": min(limit, 200),
                                     "sort": "p_value", "direction": "asc"})
        except Exception as exc:
            log.warning("GWAS Catalog API error for %s: %s", gene_symbol, exc)
            return {"gene_symbol": gene_symbol, "association_count": 0, "associations": [],
                    "error": "api_error", "note": f"GWAS Catalog API error: {type(exc).__name__}"}
        associations = self._extract_associations(data, limit)
        result = {
            "gene_symbol": gene_symbol,
            "association_count": self._total(data, associations),
            "associations": associations,
        }
        if not associations:
            result["note"] = "No GWAS associations found (gene may not be in GWAS Catalog)"
        self.cache.set(key, result)
        return result

    def trait_search(self, query: str, limit: int = 10) -> dict:
        """Search EFO traits by keyword (v2 ``/efo-traits?efo_trait=``, substring match)."""
        data, _ = self._cached_get(
            "trait_search_v2",
            "/efo-traits",
            params={"efo_trait": query, "size": min(limit, 50)},
        )
        traits = []
        for item in (data.get("_embedded", {}).get("efo_traits", []) if isinstance(data, dict) else []):
            traits.append({
                "trait": item.get("efo_trait"),
                "short_form": item.get("efo_id"),
                "uri": item.get("uri"),
            })
        return {"query": query, "count": len(traits), "traits": traits[:limit]}

    @staticmethod
    def _total(data: dict | list, associations: list) -> int:
        """All matching associations (v2 ``page.totalElements``), not just the returned page."""
        if isinstance(data, dict):
            total = (data.get("page") or {}).get("totalElements")
            if isinstance(total, int):
                return total
        return len(associations)

    def _extract_associations(self, data: dict | list, limit: int) -> list[dict]:
        """Extract associations from a v2 ``/associations`` page.

        ``genes`` are v2 ``mapped_genes`` (v1 returned author-reported genes, which
        the v2 list view does not carry). ``or_beta`` is ``or_value`` (v1
        ``orPerCopyNum``); beta-only associations give None, as in v1.
        """
        raw = []
        if isinstance(data, dict):
            raw = data.get("_embedded", {}).get("associations", [])
        elif isinstance(data, list):
            raw = data

        associations = []
        for assoc in raw[:limit]:
            alleles = assoc.get("snp_effect_allele") or []
            trait_names = [t.get("efo_trait") for t in assoc.get("efo_traits", []) if t.get("efo_trait")]
            or_beta = assoc.get("or_value")
            if or_beta is None:
                or_beta = assoc.get("or_per_copy_num")
            associations.append({
                "pvalue": assoc.get("pvalue_mantissa"),
                "pvalue_exponent": assoc.get("pvalue_exponent"),
                "risk_allele": alleles[0] if alleles else None,
                "or_beta": or_beta,
                "genes": list(assoc.get("mapped_genes") or []),
                "trait": ", ".join(trait_names) if trait_names else None,
                "traits": trait_names,
                "study": assoc.get("accession_id"),
                "pubmed_id": assoc.get("pubmed_id"),
            })
        return associations

    def new_for_variants(self, rsids: list[str], since_date: str = "") -> dict:
        """Check for new GWAS associations for a list of rsIDs.

        Args:
            rsids: List of rsIDs (max 50).
            since_date: ISO date filter (e.g. "2024-01-01"). Only return associations
                        with lastUpdateDate after this date.
        """
        if len(rsids) > 50:
            return {"error": "Too many rsIDs — max 50 per call", "count": len(rsids)}

        per_variant = {}
        for rsid in rsids:
            assoc_data = self.variant_associations(rsid, limit=100)
            associations = assoc_data.get("associations", [])

            if since_date and associations:
                # Filter by study date if available in the association data
                filtered = []
                for a in associations:
                    pmid = a.get("pubmed_id")
                    # v2 studies carry no date; the publication record does
                    if pmid and since_date:
                        try:
                            pub_data, _ = self._cached_get(
                                "publication_detail", f"/publications/{pmid}",
                            )
                            pub_date = ""
                            if isinstance(pub_data, dict):
                                pub_date = pub_data.get("publication_date") or ""
                            if pub_date >= since_date:
                                a["publication_date"] = pub_date
                                filtered.append(a)
                        except Exception:
                            # Include association if we can't verify the date
                            filtered.append(a)
                    else:
                        filtered.append(a)
                associations = filtered

            per_variant[rsid] = {
                "association_count": len(associations),
                "associations": associations,
            }

        return {
            "variant_count": len(rsids),
            "since_date": since_date or None,
            "variants": per_variant,
        }

    def validate(self) -> bool:
        try:
            self._get("/metadata")
            return True
        except Exception:
            return False
