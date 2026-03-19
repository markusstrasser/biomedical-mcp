"""Bio MCP server — 25 APIs, ~80 tools across 14 domains.

Uses FastMCP 3 mount() composition: each domain is a sub-server with
namespace-prefixed tools for discoverability.

Domains: genetics, targets, drugs, proteins, pathways, variants, clinical,
         population, panels, phenotype, expression, gwas, literature.
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
)

log = logging.getLogger(__name__)

DEFAULT_DATA_DIR = Path.home() / ".local" / "share" / "biomedical-mcp"

INSTRUCTIONS = """\
Biomedical data lookup via 27 APIs across 15 domains. Tools are namespace-prefixed.

DOMAINS:
  genetics_*    — gene annotation, IDs, sequences, nomenclature (MyGene, Ensembl, HGNC)
  targets_*     — disease-gene associations, pharmacogenetics (Open Targets)
  drugs_*       — compounds, mechanisms, labels, safety (ChEMBL, OpenFDA)
  proteins_*    — structure, function, interactions, domains (UniProt, AlphaFold, STRING, InterPro, PDB)
  pathways_*    — metabolic & signaling pathways (KEGG, Reactome)
  variants_*    — genetic variant annotations, ClinVar (MyVariant.info)
  clinical_*    — trials, ICD-10 codes, providers (ClinicalTrials.gov, NPI)
  population_*  — allele frequencies, constraint (gnomAD)
  expression_*  — tissue expression, eQTLs (GTEx)
  panels_*      — curated gene panels by condition (PanelApp)
  gwas_*        — GWAS variant-trait associations (GWAS Catalog)
  phenotype_*   — phenotype terms, gene-phenotype (HPO, Monarch)
  literature_*  — variant-level literature mining (LitVar2)
  composite_*   — multi-API compound queries (variant_context, gene_dossier)

QUICK REFERENCE:
  Gene lookup:        genetics_gene_info, genetics_ensembl_gene, genetics_gene_names
  Gene→disease:       targets_disease_associations
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
  GWAS hits:          gwas_variant_associations
  Phenotypes:         phenotype_gene_phenotypes, phenotype_hpo_search
  Literature:         literature_variant_publications
  Protein domains:    proteins_domains (InterPro)
  PDB structures:     proteins_experimental_structures, proteins_pdb_detail
  EVERYTHING about a variant: composite_variant_context
  EVERYTHING about a gene:    composite_gene_dossier
"""


def create_mcp(data_dir: Path | None = None) -> FastMCP:
    data_dir = data_dir or Path(os.environ.get("BIOMEDICAL_MCP_DATA", DEFAULT_DATA_DIR))
    data_dir.mkdir(parents=True, exist_ok=True)

    cache = Cache(data_dir / "cache.db")
    log.info("bio-mcp started (cache: %s)", data_dir / "cache.db")

    main = FastMCP(
        "bio",
        instructions=INSTRUCTIONS,
        middleware=[TelemetryMiddleware()],
    )

    # Existing domains (7)
    main.mount(genetics.create_server(cache), namespace="genetics")
    main.mount(targets.create_server(cache), namespace="targets")
    main.mount(drugs.create_server(cache), namespace="drugs")
    main.mount(proteins.create_server(cache), namespace="proteins")
    main.mount(pathways.create_server(cache), namespace="pathways")
    main.mount(variants.create_server(cache), namespace="variants")
    main.mount(clinical.create_server(cache), namespace="clinical")

    # New domains (6)
    main.mount(population.create_server(cache), namespace="population")
    main.mount(panels.create_server(cache), namespace="panels")
    main.mount(phenotype.create_server(cache), namespace="phenotype")
    main.mount(expression.create_server(cache), namespace="expression")
    main.mount(gwas.create_server(cache), namespace="gwas")
    main.mount(literature.create_server(cache), namespace="literature")
    main.mount(composite.create_server(cache), namespace="composite")

    return main


def main():
    mcp = create_mcp()
    mcp.run()
