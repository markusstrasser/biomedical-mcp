"""USDA FoodData Central API client — nutrient composition lookup."""

from __future__ import annotations

import logging
import os

from biomedical_mcp.base_client import BaseClient
from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)


class Nutrition(BaseClient):
    domain_name = "nutrition"

    def __init__(self, cache: Cache):
        super().__init__(cache)
        self.api_key = os.environ.get("FDC_API_KEY", "DEMO_KEY")

    def food_nutrients(self, query: str, limit: int = 5) -> dict:
        """Search USDA FoodData Central for nutrient composition."""
        key = self._cache_key("fdc_search", query, str(limit))
        cached = self.cache.get(key, max_age_days=self.ttl_days)
        if cached is not None:
            return {**cached, "_cache_hit": True}

        try:
            search_data = self._get(
                "/fdc/v1/foods/search",
                params={
                    "query": query,
                    "pageSize": min(limit, 25),
                    "api_key": self.api_key,
                },
            )
        except Exception as exc:
            log.warning("USDA FDC search error for %s: %s", query, exc)
            return {
                "query": query,
                "foods": [],
                "error": "api_error",
                "note": f"USDA FDC error: {type(exc).__name__}",
            }

        foods = []
        for item in (search_data.get("foods") or [])[:limit]:
            fdc_id = item.get("fdcId")
            food = {
                "fdc_id": fdc_id,
                "description": item.get("description"),
                "brand_owner": item.get("brandOwner"),
                "data_type": item.get("dataType"),
                "serving_size": item.get("servingSize"),
                "serving_size_unit": item.get("servingSizeUnit"),
            }

            # Extract nutrients from search results (already included)
            nutrients = {}
            for n in (item.get("foodNutrients") or []):
                name = n.get("nutrientName")
                if name:
                    nutrients[name] = {
                        "value": n.get("value"),
                        "unit": n.get("unitName"),
                    }
            food["nutrients"] = nutrients

            # If search results lack nutrients, fetch detail
            if not nutrients and fdc_id:
                try:
                    detail = self._get(
                        f"/fdc/v1/food/{fdc_id}",
                        params={"api_key": self.api_key},
                    )
                    for n in (detail.get("foodNutrients") or []):
                        nutrient = n.get("nutrient", {})
                        name = nutrient.get("name")
                        if name:
                            nutrients[name] = {
                                "value": n.get("amount"),
                                "unit": nutrient.get("unitName"),
                            }
                    food["nutrients"] = nutrients
                except Exception as exc:
                    log.debug("FDC detail fetch failed for %s: %s", fdc_id, exc)

            foods.append(food)

        result = {
            "query": query,
            "total_hits": search_data.get("totalHits", len(foods)),
            "food_count": len(foods),
            "foods": foods,
        }
        self.cache.set(key, result)
        return result

    def validate(self) -> bool:
        """Check USDA FDC endpoint reachability."""
        try:
            self._get(
                "/fdc/v1/foods/search",
                params={"query": "apple", "pageSize": 1, "api_key": self.api_key},
            )
            return True
        except Exception:
            return False
