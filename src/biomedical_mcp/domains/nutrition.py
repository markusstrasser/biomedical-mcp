"""Nutrition domain — food nutrient composition (USDA FoodData Central)."""

from fastmcp import FastMCP

from biomedical_mcp.nutrition import Nutrition
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    client = Nutrition(cache)

    server = FastMCP("nutrition")

    @server.tool(tags={"food"})
    def food_nutrients(query: str, limit: int = 5) -> dict:
        """Search USDA FoodData Central for nutrient composition per 100g.

        Returns matching foods with full nutrient profiles (calories, macros,
        vitamins, minerals). Uses DEMO_KEY by default; set FDC_API_KEY env var
        for higher rate limits.

        Args:
            query: Food search query (e.g. "chicken breast", "brown rice", "kale").
            limit: Max food matches (default 5, max 25).
        """
        return client.food_nutrients(query, limit=limit)

    return server
