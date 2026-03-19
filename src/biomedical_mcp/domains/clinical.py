"""Clinical domain — trials, ICD-10 codes, provider registry."""

from fastmcp import FastMCP

from biomedical_mcp.clinicaltrials import ClinicalTrials
from biomedical_mcp.icd10 import ICD10
from biomedical_mcp.npi import NPI
from biomedical_mcp.cache import Cache


def create_server(cache: Cache) -> FastMCP:
    ct = ClinicalTrials(cache)
    icd10 = ICD10(cache)
    npi = NPI(cache)

    server = FastMCP("clinical")

    @server.tool()
    def trial_search(
        condition: str | None = None,
        intervention: str | None = None,
        gene: str | None = None,
        status: str | None = None,
        limit: int = 20,
    ) -> dict:
        """Search clinical trials by condition, intervention, gene, or status.

        Args:
            condition: Disease/condition (e.g. "breast cancer").
            intervention: Drug or therapy (e.g. "pembrolizumab").
            gene: Gene name for gene-related trials (e.g. "BRCA1").
            status: Filter: RECRUITING, COMPLETED, ACTIVE_NOT_RECRUITING, etc.
            limit: Max results (default 20, max 100).
        """
        return ct.search(condition=condition, intervention=intervention, gene=gene, status=status, limit=limit)

    @server.tool()
    def trial_detail(nct_id: str) -> dict:
        """Get full trial details: eligibility, outcomes, sponsor, design.

        Args:
            nct_id: ClinicalTrials.gov identifier (e.g. "NCT04379596").
        """
        return ct.trial_detail(nct_id)

    @server.tool()
    def trial_stats(condition: str) -> dict:
        """Get trial counts by phase and status for a condition.

        Args:
            condition: Disease/condition (e.g. "lung cancer").
        """
        return ct.stats(condition)

    @server.tool()
    def icd10_search(query: str, limit: int = 10) -> dict:
        """Search ICD-10-CM codes by description.

        Args:
            query: Search term (e.g. "diabetes", "hypertension").
            limit: Max results (default 10, max 50).
        """
        results = icd10.search(query, limit=limit)
        return {"query": query, "results": results}

    @server.tool()
    def icd10_lookup(code: str) -> dict:
        """Look up an ICD-10-CM code for its description and subcodes.

        Args:
            code: ICD-10-CM code (e.g. "E11", "E11.65").
        """
        result = icd10.lookup(code)
        if not result:
            return {"error": f"ICD-10 code not found: {code}"}
        return result

    @server.tool()
    def provider_search(
        name: str | None = None,
        specialty: str | None = None,
        city: str | None = None,
        state: str | None = None,
        limit: int = 10,
    ) -> dict:
        """Search the NPI registry for healthcare providers.

        Args:
            name: Provider name ("John Smith" or just "Smith").
            specialty: Taxonomy description (e.g. "Cardiology", "Internal Medicine").
            city: City name.
            state: Two-letter state code (e.g. "CA", "NY").
            limit: Max results (default 10, max 200).
        """
        results = npi.search(name=name, specialty=specialty, city=city, state=state, limit=limit)
        return {"results": results, "count": len(results)}

    @server.tool()
    def provider_lookup(npi_number: str) -> dict:
        """Look up a specific NPI number for full provider details.

        Args:
            npi_number: 10-digit NPI number.
        """
        result = npi.lookup(npi_number)
        if not result:
            return {"error": f"NPI not found: {npi_number}"}
        return result

    return server
