"""Supplement safety APIs — OpenFDA CAERS, NIH DSLD, PharmGKB."""

from __future__ import annotations

import logging
from collections import Counter

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class Supplements(BaseClient):
    domain_name = "supplements"

    def __init__(self, cache: Cache):
        super().__init__(cache)

    # --- OpenFDA CAERS (dietary supplement adverse events) ---

    def adverse_events(self, product: str, limit: int = 20) -> dict:
        """Query OpenFDA CAERS for dietary supplement adverse event reports."""
        key = self._cache_key("caers", product, str(limit))
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        try:
            data = self._get(
                "/food/event.json",
                params={
                    "search": f'products.name_brand:"{product}"',
                    "limit": min(limit, 100),
                },
            )
        except Exception as exc:
            log.warning("OpenFDA CAERS error for %s: %s", product, exc)
            return {
                "product": product,
                "total_reports": 0,
                "error": "api_error",
                "note": f"OpenFDA CAERS error: {type(exc).__name__}",
            }

        meta_total = data.get("meta", {}).get("results", {}).get("total", 0)
        results = data.get("results", [])

        # Aggregate reactions and outcomes across all returned reports
        all_reactions: list[str] = []
        outcome_counts: Counter[str] = Counter()
        serious_count = 0
        sample_reports = []

        for report in results:
            reactions = report.get("reactions", [])
            all_reactions.extend(reactions)
            outcomes = report.get("outcomes", [])
            for o in outcomes:
                outcome_counts[o] += 1
            if any(o in ("SERIOUS INJURIES", "DEATH", "LIFE THREATENING")
                   for o in outcomes):
                serious_count += 1
            if len(sample_reports) < 5:
                sample_reports.append({
                    "date": report.get("date_started"),
                    "reactions": reactions[:5],
                    "outcomes": outcomes,
                    "products": [
                        p.get("name_brand")
                        for p in (report.get("products") or [])[:3]
                    ],
                })

        reaction_counts = Counter(all_reactions)
        top_reactions = [
            {"reaction": r, "count": c}
            for r, c in reaction_counts.most_common(10)
        ]

        result = {
            "product": product,
            "total_reports": meta_total,
            "serious_count": serious_count,
            "top_reactions": top_reactions,
            "outcome_distribution": dict(outcome_counts),
            "sample_reports": sample_reports,
        }
        self.cache.set(key, result)
        return result

    # --- NIH DSLD (Dietary Supplement Label Database) ---

    def label_verify(self, product: str, size: int = 5) -> dict:
        """Query NIH DSLD for supplement label data."""
        key = self._cache_key("dsld", product, str(size))
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        # DSLD uses a different base URL — override per-call
        dsld_base = "https://api.ods.od.nih.gov/dsld/v9"
        try:
            self._rate_wait()
            resp = self.client.get(
                f"{dsld_base}/search-filter",
                params={"q": product, "from": 0, "size": min(size, 25)},
            )
            resp.raise_for_status()
            search_data = resp.json()
        except Exception as exc:
            log.warning("NIH DSLD search error for %s: %s", product, exc)
            return {
                "product": product,
                "matches": [],
                "error": "api_error",
                "note": f"NIH DSLD error: {type(exc).__name__}",
            }

        hits = search_data.get("hits", search_data.get("results", []))
        if isinstance(hits, dict):
            hits = hits.get("hits", [])

        matches = []
        for hit in hits[:size]:
            source = hit.get("_source", hit) if isinstance(hit, dict) else {}
            dsld_id = source.get("dsld_id") or hit.get("_id")
            match = {
                "dsld_id": dsld_id,
                "name": source.get("product_name") or source.get("LanguaL Product Name"),
                "brand": source.get("brand_name"),
                "on_market": source.get("on_market") or source.get("Market Status"),
            }

            # Fetch label detail if we have an ID
            if dsld_id:
                try:
                    self._rate_wait()
                    label_resp = self.client.get(f"{dsld_base}/label/{dsld_id}")
                    label_resp.raise_for_status()
                    label_data = label_resp.json()

                    ingredients = []
                    for ing in (label_data.get("ingredients") or
                                label_data.get("Ingredients") or []):
                        ingredients.append({
                            "name": ing.get("ingredient_name") or ing.get("Name"),
                            "amount": ing.get("amount") or ing.get("Amount"),
                            "unit": ing.get("unit") or ing.get("Unit"),
                            "dv_percent": ing.get("dv_percent") or ing.get("% Daily Value"),
                            "form": ing.get("ingredient_form") or ing.get("Form"),
                        })
                    match["ingredients"] = ingredients
                except Exception as exc:
                    log.debug("DSLD label fetch failed for %s: %s", dsld_id, exc)
                    match["ingredients"] = []

            matches.append(match)

        result = {
            "product": product,
            "match_count": len(matches),
            "matches": matches,
        }
        self.cache.set(key, result)
        return result

    # --- PharmGKB (pharmacogenomics) ---

    def pharmgkb_lookup(self, query: str) -> dict:
        """Query PharmGKB for pharmacogenomics data (gene-drug associations)."""
        key = self._cache_key("pgkb", query)
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        pgkb_base = "https://api.pharmgkb.org/v1/data"

        # Try as gene first, then as drug
        annotations = []
        drug_labels = []

        # Clinical annotations by gene
        try:
            self._rate_wait()
            resp = self.client.get(
                f"{pgkb_base}/clinicalAnnotation",
                params={
                    "view": "base",
                    "location.genes.symbol": query,
                },
            )
            resp.raise_for_status()
            data = resp.json()
            for ann in (data.get("data") or []):
                annotations.append({
                    "id": ann.get("id"),
                    "gene": query,
                    "drugs": [c.get("name") for c in (ann.get("relatedChemicals") or [])],
                    "phenotypes": [p.get("name") for p in (ann.get("relatedDiseases") or [])],
                    "level_of_evidence": ann.get("evidenceLevel") or ann.get("level"),
                })
        except Exception as exc:
            log.debug("PharmGKB gene lookup for %s: %s", query, exc)

        # Drug labels
        try:
            self._rate_wait()
            resp = self.client.get(
                f"{pgkb_base}/drugLabel",
                params={
                    "view": "base",
                    "relatedChemicals.name": query,
                },
            )
            resp.raise_for_status()
            data = resp.json()
            for label in (data.get("data") or []):
                drug_labels.append({
                    "id": label.get("id"),
                    "name": label.get("name"),
                    "source": label.get("source"),
                    "genes": [g.get("symbol") for g in (label.get("relatedGenes") or [])],
                    "chemicals": [c.get("name") for c in (label.get("relatedChemicals") or [])],
                })
        except Exception as exc:
            log.debug("PharmGKB drug label lookup for %s: %s", query, exc)

        if not annotations and not drug_labels:
            return {
                "query": query,
                "clinical_annotations": [],
                "drug_labels": [],
                "note": "No PharmGKB data found for query",
            }

        result = {
            "query": query,
            "clinical_annotation_count": len(annotations),
            "clinical_annotations": annotations,
            "drug_label_count": len(drug_labels),
            "drug_labels": drug_labels,
        }
        self.cache.set(key, result)
        return result

    def validate(self) -> bool:
        """Check OpenFDA CAERS endpoint reachability."""
        try:
            self._get("/food/event.json", params={"limit": 1})
            return True
        except Exception:
            return False
