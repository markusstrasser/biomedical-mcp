"""Bio MCP server — consolidated entity-composite surface over 35 biomedical APIs.

Default profile ("composite"): 5 entity composites (gene_dossier, variant_context,
drug_profile, protein_profile, disease_profile) + describe_sections + bio_search +
health_check — 8 tools, ~1.2k at-rest tokens. Each composite fans out across many
APIs and returns a standardized partial-failure envelope.

`full` profile (BIOMEDICAL_MCP_PROFILE=full): additionally mounts ~85 raw per-source
tools across 20 domains (genetics, targets, drugs, proteins, pathways, variants,
clinical, population, panels, phenotype, expression, gwas, literature, bloodgroups,
rare_disease, curation, supplements, nutrition) for long-tail capability + maintenance.

Consolidation rationale: decisions/2026-06-07-biomedical-mcp-tool-consolidation.md (agent-infra).
"""

import logging
import os
from pathlib import Path

from fastmcp import FastMCP

from biomedical_mcp.cache import Cache
from biomedical_mcp.middleware import TelemetryMiddleware
from biomedical_mcp.domains import (
    genetics, targets, drugs, proteins, pathways, variants, clinical,
    population, panels, phenotype, expression, gwas, literature, composite,
    bloodgroups, rare_disease, curation, supplements, nutrition,
)

log = logging.getLogger(__name__)

DEFAULT_DATA_DIR = Path.home() / ".local" / "share" / "biomedical-mcp"

INSTRUCTIONS = """\
Biomedical data lookup via 32 APIs across 20 domains. Tools are namespace-prefixed.

DOMAINS:
  genetics_*      — gene annotation, IDs, sequences, nomenclature (MyGene, Ensembl, HGNC)
  targets_*       — disease-gene associations, pharmacogenetics (Open Targets)
  drugs_*         — compounds, mechanisms, labels, safety (ChEMBL, OpenFDA)
  proteins_*      — structure, function, interactions, domains (UniProt, AlphaFold, STRING, InterPro, PDB)
  pathways_*      — metabolic & signaling pathways (KEGG, Reactome)
  variants_*      — genetic variant annotations, ClinVar (MyVariant.info)
  clinical_*      — trials, ICD-10 codes, providers (ClinicalTrials.gov, NPI)
  population_*    — allele frequencies, constraint (gnomAD)
  expression_*    — tissue expression, eQTLs (GTEx)
  panels_*        — curated gene panels by condition (PanelApp)
  gwas_*          — GWAS variant-trait associations, variant surveillance (GWAS Catalog)
  phenotype_*     — phenotype terms, gene-phenotype (HPO, Monarch)
  literature_*    — variant-level literature mining (LitVar2)
  composite_*     — multi-API compound queries (variant_context, gene_dossier)
  bloodgroups_*   — blood group systems, alleles, antigens (ISBT)
  rare_disease_*  — rare disease data, inheritance, epidemiology (Orphanet)
  curation_*      — gene-disease validity, dosage sensitivity (ClinGen)
  supplements_*   — adverse events (FDA CAERS), label verification (NIH DSLD), PGx (PharmGKB)
  nutrition_*     — food nutrient composition (USDA FoodData Central)

QUICK REFERENCE:
  Gene lookup:        genetics_gene_info, genetics_ensembl_gene, genetics_gene_names
  Gene→disease:       targets_disease_associations, curation_gene_validity
  Gene→pathways:      pathways_kegg_gene, pathways_reactome_gene
  Gene→interactions:  proteins_interactions
  Variant annotation: variants_lookup, variants_clinvar
  Drug info:          drugs_compound, drugs_mechanism, drugs_label
  Protein structure:  proteins_structure (AlphaFold)
  Clinical trials:    clinical_trial_search
  Enrichment:         proteins_enrichment (STRING), pathways_reactome_enrichment
  Population freq:    population_variant_frequency
  Gene constraint:    population_gene_constraint
  Gene panels:        panels_gene_panels
  Tissue expression:  expression_gene_expression
  GWAS hits:          gwas_variant_associations, gwas_new_for_variants
  Phenotypes:         phenotype_gene_phenotypes, phenotype_hpo_search
  Literature:         literature_variant_publications
  Protein domains:    proteins_domains (InterPro)
  PDB structures:     proteins_experimental_structures, proteins_pdb_detail
  Blood groups:       bloodgroups_systems, bloodgroups_alleles, bloodgroups_search_alleles
  Rare diseases:      rare_disease_gene_rare_diseases, rare_disease_disease_natural_history
  Gene-disease proof: curation_gene_validity, curation_gene_dosage
  Supplement safety:  supplements_adverse_events, supplements_label_verify
  PGx interactions:   supplements_pharmgkb_lookup
  Food nutrients:     nutrition_food_nutrients
  EVERYTHING about a variant: composite_variant_context
  EVERYTHING about a gene:    composite_gene_dossier
"""


COMPOSITE_INSTRUCTIONS = """\
Biomedical data lookup via 5 entity composites (default profile). Each composite fans
out across many upstream APIs and returns a standardized partial-failure envelope:
one dead source degrades only its section (status="error"), never the whole result.

ENTITY COMPOSITES (pass `identifier`; optional `sections` to fetch a subset):
  composite_gene_dossier(identifier)      — gene: nomenclature, constraint, panels,
                                            expression, phenotypes, validity, dosage, rare diseases
  composite_variant_context(identifier)   — variant (rsID/HGVS/"BRAF V600E"): annotation
                                            (ClinVar/predictions/CGI/COSMIC), pop freq, literature, GWAS
  composite_drug_profile(identifier)      — drug: compound, mechanism, label, adverse events (+indications, recalls)
  composite_protein_profile(identifier)   — protein (gene symbol/UniProt acc): function, variants,
                                            structure, domains, experimental structures (+interactions)
  composite_disease_profile(identifier)   — disease (name/MONDO/ORPHA/OMIM/EFO): associations,
                                            phenotypes, natural history, epidemiology, genes, coding

DISCOVERY:
  composite_describe_sections(entity)     — list a composite's sections, sources, descriptions
  composite_bio_search(entity, query)     — resolve a name/keyword to identifiers, then call the composite
  composite_health_check()                — upstream API connectivity

Set BIOMEDICAL_MCP_PROFILE=full to additionally expose ~85 raw per-source tools
(sequence retrieval, clinical trials, pathways, supplements, nutrition, blood groups,
HPO/gene search, eQTLs, and other long-tail capabilities not yet in a composite section).
"""


def create_mcp(data_dir: Path | None = None) -> FastMCP:
    data_dir = data_dir or Path(os.environ.get("BIOMEDICAL_MCP_DATA", DEFAULT_DATA_DIR))
    data_dir.mkdir(parents=True, exist_ok=True)

    cache = Cache(data_dir / "cache.db")
    evicted = cache.cleanup(max_age_days=90)
    if evicted:
        log.info("bio-mcp cache cleanup: evicted %d expired entries", evicted)
    log.info("bio-mcp started (cache: %s)", data_dir / "cache.db")

    profile = os.environ.get("BIOMEDICAL_MCP_PROFILE", "composite").lower()

    main = FastMCP(
        "bio",
        instructions=INSTRUCTIONS if profile == "full" else COMPOSITE_INSTRUCTIONS,
        middleware=[TelemetryMiddleware()],
    )

    # Default ("composite") profile: the consolidated entity-composite surface only
    # (gene_dossier, variant_context, drug_profile, protein_profile, disease_profile,
    # describe_sections, bio_search, health_check). ~85 raw per-source tools are NOT
    # loaded — this is the at-rest-token win. See
    # agent-infra/decisions/2026-06-07-biomedical-mcp-tool-consolidation.md.
    main.mount(composite.create_server(cache), namespace="composite")

    if profile == "full":
        # Escape hatch: every raw per-source tool, for maintenance / expert use /
        # long-tail capability not yet promoted into a composite section.
        log.info("bio-mcp profile=full — mounting all raw domains")
        main.mount(genetics.create_server(cache), namespace="genetics")
        main.mount(targets.create_server(cache), namespace="targets")
        main.mount(drugs.create_server(cache), namespace="drugs")
        main.mount(proteins.create_server(cache), namespace="proteins")
        main.mount(pathways.create_server(cache), namespace="pathways")
        main.mount(variants.create_server(cache), namespace="variants")
        main.mount(clinical.create_server(cache), namespace="clinical")
        main.mount(population.create_server(cache), namespace="population")
        main.mount(panels.create_server(cache), namespace="panels")
        main.mount(phenotype.create_server(cache), namespace="phenotype")
        main.mount(expression.create_server(cache), namespace="expression")
        main.mount(gwas.create_server(cache), namespace="gwas")
        main.mount(literature.create_server(cache), namespace="literature")
        main.mount(bloodgroups.create_server(cache), namespace="bloodgroups")
        main.mount(rare_disease.create_server(cache), namespace="rare_disease")
        main.mount(curation.create_server(cache), namespace="curation")
        main.mount(supplements.create_server(cache), namespace="supplements")
        main.mount(nutrition.create_server(cache), namespace="nutrition")

    # Prompt templates — structured review workflows
    @main.prompt()
    def variant_review(variant_id: str, gene: str = "") -> str:
        """Structured variant review checklist using bio MCP tools."""
        gene_part = f" in {gene}" if gene else ""
        return f"""Review variant {variant_id}{gene_part}. Use bio MCP tools for each step:

1. **Population frequency** → population_variant_frequency("{variant_id}")
2. **ClinVar classification** → variants_lookup("{variant_id}")
3. **In silico predictions** → (included in variants_lookup: CADD, SIFT, PolyPhen2)
4. **Gene constraint** → population_gene_constraint("{gene}") [pLI, LOEUF]
5. **Tissue expression** → expression_gene_expression(gene_symbol="{gene}")
6. **Gene panels** → panels_gene_panels("{gene}") [clinical relevance]
7. **Literature** → literature_variant_publications("{variant_id}")
8. **Protein domain** → proteins_domains(uniprot_id=...) [is variant in critical domain?]

Or use composite_variant_context("{variant_id}", gene_symbol="{gene}") for steps 1-7 in one call.

Assess evidence tier per constitutional hierarchy (C3+ for clinical action)."""

    @main.prompt()
    def gene_dossier_prompt(gene: str) -> str:
        """Comprehensive gene investigation workflow."""
        return f"""Build a complete dossier for {gene}. Use bio MCP tools:

1. **Official names** → genetics_gene_names("{gene}") [aliases, previous symbols]
2. **Gene function** → genetics_gene_info("{gene}") [GO terms, pathways]
3. **Constraint** → population_gene_constraint("{gene}") [LoF intolerance]
4. **Disease associations** → targets_disease_associations(gene_symbol="{gene}")
5. **Clinical panels** → panels_gene_panels("{gene}") [curated relevance]
6. **Expression profile** → expression_gene_expression(gene_symbol="{gene}")
7. **Phenotypes** → phenotype_gene_phenotypes("{gene}") [HPO terms]
8. **GWAS associations** → gwas_gene_associations("{gene}")
9. **Drug targets** → targets_target_info(gene_symbol="{gene}")
10. **Protein structure** → proteins_structure(uniprot_id=...)

Or use composite_gene_dossier("{gene}") for steps 1-7 in one call."""

    @main.prompt()
    def pgx_review(gene: str) -> str:
        """Pharmacogenomics gene review workflow."""
        return f"""Review pharmacogenomics for {gene}:

1. **Star alleles** → drugs_star_alleles("{gene}") [PharmVar definitions]
2. **PGx interactions** → targets_pharmacogenetics(gene_symbol="{gene}") [OT PGx data]
3. **Drug mechanisms** → drugs_mechanism(compound_name=...) [for each interacting drug]
4. **FDA labels** → drugs_label(drug_name=..., sections=["boxed_warning","clinical_pharmacology"])
5. **Adverse events** → drugs_adverse_events(drug_name=...)

Cross-reference with ClinPGx guidelines if available."""

    return main


def main():
    mcp = create_mcp()
    mcp.run()
