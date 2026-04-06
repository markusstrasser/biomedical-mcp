"""Centralized domain registry — URLs, TTLs, rate limits, enable flags."""

from __future__ import annotations

# TTL tiers (days)
TTL_STABLE = 30      # gene names, constraint, panels, protein domains, allele definitions
TTL_CURATED = 7      # disease associations, expression, GWAS, ClinVar
TTL_VOLATILE = 1     # clinical trials, literature

DOMAINS: dict[str, dict] = {
    # Existing APIs
    "mygene":          {"base_url": "https://myvariant.info/v1",                        "ttl_days": TTL_CURATED,  "rate_limit": 3.0},
    "ensembl":         {"base_url": "https://rest.ensembl.org",                         "ttl_days": TTL_STABLE,   "rate_limit": 10.0},
    "opentargets":     {"base_url": "https://api.platform.opentargets.org/api/v4/graphql", "ttl_days": TTL_CURATED, "rate_limit": 1.0},
    "chembl":          {"base_url": "https://www.ebi.ac.uk/chembl/api/data",            "ttl_days": TTL_STABLE,   "rate_limit": 3.0},
    "openfda":         {"base_url": "https://api.fda.gov",                              "ttl_days": TTL_CURATED,  "rate_limit": 0.6},
    "uniprot":         {"base_url": "https://rest.uniprot.org",                         "ttl_days": TTL_STABLE,   "rate_limit": 3.0},
    "alphafold":       {"base_url": "https://alphafold.ebi.ac.uk/api",                  "ttl_days": TTL_STABLE,   "rate_limit": 3.0},
    "string":          {"base_url": "https://string-db.org/api",                        "ttl_days": TTL_CURATED,  "rate_limit": 1.0},
    "kegg":            {"base_url": "https://rest.kegg.jp",                             "ttl_days": TTL_STABLE,   "rate_limit": 3.0},
    "reactome":        {"base_url": "https://reactome.org/ContentService",              "ttl_days": TTL_STABLE,   "rate_limit": 3.0},
    "myvariant":       {"base_url": "https://myvariant.info/v1",                        "ttl_days": TTL_CURATED,  "rate_limit": 3.0},
    "clinicaltrials":  {"base_url": "https://clinicaltrials.gov/api/v2",                "ttl_days": TTL_VOLATILE, "rate_limit": 3.0},
    "icd10":           {"base_url": "https://clinicaltables.nlm.nih.gov/api",           "ttl_days": TTL_STABLE,   "rate_limit": 5.0},
    "npi":             {"base_url": "https://npiregistry.cms.hhs.gov/api",              "ttl_days": TTL_CURATED,  "rate_limit": 3.0},
    # New Tier 1 APIs
    "gnomad":          {"base_url": "https://gnomad.broadinstitute.org/api",            "ttl_days": TTL_STABLE,   "rate_limit": 1.0},
    "panelapp":        {"base_url": "https://panelapp.genomicsengland.co.uk/api/v1",    "ttl_days": TTL_STABLE,   "rate_limit": 3.0},
    "hpo":             {"base_url": "https://ontology.jax.org/api",                     "ttl_days": TTL_STABLE,   "rate_limit": 3.0},
    "monarch":         {"base_url": "https://api.monarchinitiative.org/v3/api",         "ttl_days": TTL_CURATED,  "rate_limit": 3.0},
    # New Tier 2 APIs
    "gtex":            {"base_url": "https://gtexportal.org/api/v2",                    "ttl_days": TTL_CURATED,  "rate_limit": 3.0},
    "hgnc":            {"base_url": "https://rest.genenames.org",                       "ttl_days": TTL_STABLE,   "rate_limit": 5.0},
    "gwas_catalog":    {"base_url": "https://www.ebi.ac.uk/gwas/rest/api",              "ttl_days": TTL_CURATED,  "rate_limit": 1.0},
    # New Tier 3 APIs
    "litvar":          {"base_url": "https://www.ncbi.nlm.nih.gov/research/litvar2-api", "ttl_days": TTL_VOLATILE, "rate_limit": 3.0},
    "interpro":        {"base_url": "https://www.ebi.ac.uk/interpro/api",               "ttl_days": TTL_STABLE,   "rate_limit": 3.0},
    "pdb":             {"base_url": "https://data.rcsb.org",                            "ttl_days": TTL_STABLE,   "rate_limit": 5.0},
    # New: ISBT, Orphanet, ClinGen
    "isbt":            {"base_url": "https://api-blooddatabase.isbtweb.org",            "ttl_days": TTL_STABLE,   "rate_limit": 3.0},
    "orphanet":        {"base_url": "https://api.orphadata.com",                        "ttl_days": TTL_CURATED,  "rate_limit": 3.0},
    # Supplements & nutrition
    "supplements":     {"base_url": "https://api.fda.gov",                              "ttl_days": TTL_CURATED,  "rate_limit": 4.0},
    "nutrition":       {"base_url": "https://api.nal.usda.gov",                         "ttl_days": TTL_STABLE,   "rate_limit": 3.0},
}


def get_domain(name: str) -> dict:
    """Get config for a domain, or sensible defaults."""
    return DOMAINS.get(name, {"ttl_days": TTL_CURATED, "rate_limit": 1.0})
