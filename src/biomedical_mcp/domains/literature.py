"""Literature domain — variant-level literature mining (LitVar2, NCBI)."""

from fastmcp import FastMCP

from biomedical_mcp.litvar import LitVar
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    litvar = LitVar(cache)

    server = FastMCP("literature")

    @server.tool(tags={"variant-review", "literature"})
    def variant_publications(variant_id: str, limit: int = 20) -> dict:
        """Papers mentioning a genetic variant. Uses NCBI LitVar2.

        Note: Works best with rsIDs. Rare variants without rsIDs may return 0 hits.

        Args:
            variant_id: rsID (e.g. "rs1801133") or HGVS notation.
            limit: Max publications (default 20).
        """
        return litvar.variant_publications(variant_id, limit=limit)

    return server
