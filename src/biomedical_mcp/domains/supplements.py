"""Supplements domain — adverse events, label verification, pharmacogenomics."""

from fastmcp import FastMCP

from biomedical_mcp.supplements import Supplements
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    client = Supplements(cache)

    server = FastMCP("supplements")

    @server.tool(tags={"safety"})
    def adverse_events(product: str, limit: int = 20) -> dict:
        """Dietary supplement adverse event reports from FDA CAERS.

        Returns total reports, serious count, top 10 reactions, outcome distribution,
        and up to 5 sample reports.

        Args:
            product: Supplement product name (e.g. "melatonin", "ashwagandha", "fish oil").
            limit: Max reports to analyze (default 20, max 100).
        """
        return client.adverse_events(product, limit=limit)

    @server.tool(tags={"labels"})
    def label_verify(product: str, size: int = 5) -> dict:
        """Verify supplement label claims via NIH Dietary Supplement Label Database (DSLD).

        Returns matching products with ingredients, amounts, units, daily value %,
        and market status.

        Args:
            product: Supplement name to search (e.g. "vitamin D", "magnesium glycinate").
            size: Max product matches (default 5, max 25).
        """
        return client.label_verify(product, size=size)

    @server.tool(tags={"pgx"})
    def pharmgkb_lookup(query: str) -> dict:
        """PharmGKB pharmacogenomics lookup — gene-drug interactions and clinical annotations.

        Searches by gene symbol (e.g. "CYP2D6") or drug name (e.g. "warfarin").
        Returns clinical annotations with evidence levels and drug label information.

        Args:
            query: Gene symbol or drug name (e.g. "CYP2D6", "warfarin", "codeine").
        """
        return client.pharmgkb_lookup(query)

    return server
