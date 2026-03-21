"""Blood groups domain — ISBT blood group systems, alleles, antigens."""

from fastmcp import FastMCP

from biomedical_mcp.isbt import ISBT
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    isbt = ISBT(cache)

    server = FastMCP("bloodgroups")

    @server.tool(tags={"blood-group"})
    def systems() -> dict:
        """List all ISBT blood group systems with gene/antigen/allele counts.

        Returns: 48 blood group systems (ABO, RH, MNS, FY, JK, etc.) plus
        transcription factors, collections, and series.
        """
        return isbt.systems()

    @server.tool(tags={"blood-group"})
    def system_detail(system_id: int) -> dict:
        """Full detail for a blood group system: genes (with coordinates), antigens.

        Args:
            system_id: ISBT system ID (1=ABO, 4=RH, 8=FY, 9=JK, etc.).
        """
        return isbt.system_detail(system_id)

    @server.tool(tags={"blood-group"})
    def alleles(system_symbol: str) -> dict:
        """All alleles for a blood group system — names, phenotypes, SNPs, classifications.

        Args:
            system_symbol: ISBT system symbol (e.g. "ABO", "RHD", "RHCE", "FY", "JK", "KEL").
        """
        return isbt.alleles(system_symbol)

    @server.tool(tags={"blood-group"})
    def search_alleles(isbt_allele: str = "", system_symbol: str = "", phenotype: str = "") -> dict:
        """Search blood group alleles by allele name, system symbol, or phenotype.

        Field-specific search — provide at least one filter.

        Args:
            isbt_allele: ISBT allele name (e.g. "FY*02", "ABO*A1.01").
            system_symbol: ISBT system symbol (e.g. "FY", "ABO", "RHD").
            phenotype: ISBT phenotype (e.g. "O", "Fy(a-b-)", "Kel:1").
        """
        return isbt.search_alleles(
            isbt_allele=isbt_allele or None,
            system_symbol=system_symbol or None,
            phenotype=phenotype or None,
        )

    @server.tool(tags={"blood-group"})
    def antigens(system_symbol: str = "") -> dict:
        """Get blood group antigens, optionally filtered by system.

        Args:
            system_symbol: ISBT system symbol to filter (e.g. "ABO", "RH"). Empty = all antigens.
        """
        return isbt.antigens(system_symbol or None)

    @server.tool(tags={"blood-group", "variant-review"})
    def variant_lookup(rsid: str = "", chromosome: str = "", position: int = 0,
                       ref: str = "", alt: str = "") -> dict:
        """Look up blood group variants by rsID or GRCh38 genomic coordinates.

        Args:
            rsid: dbSNP rsID (e.g. "rs505922").
            chromosome: GRCh38 chromosome (e.g. "9"). Use with position.
            position: GRCh38 position. Use with chromosome.
            ref: Reference allele.
            alt: Alternative allele.
        """
        return isbt.variant_lookup(
            rsid=rsid or None,
            chromosome=chromosome or None,
            position=position if position else None,
            ref=ref or None,
            alt=alt or None,
        )

    return server
