"""Biomedical MCP server — 14 APIs, 41 tools across 7 domains.

Uses FastMCP 3 mount() composition: each domain is a sub-server with
namespace-prefixed tools for discoverability.

Domains: genetics, targets, drugs, proteins, pathways, variants, clinical.
"""

import logging
import os
from pathlib import Path

from fastmcp import FastMCP

from biomedical_mcp.cache import Cache
from biomedical_mcp.middleware import TelemetryMiddleware
from biomedical_mcp.domains import genetics, targets, drugs, proteins, pathways, variants, clinical

log = logging.getLogger(__name__)

DEFAULT_DATA_DIR = Path.home() / ".local" / "share" / "biomedical-mcp"

INSTRUCTIONS = """\
Biomedical data lookup via 14 APIs across 7 domains. Tools are namespace-prefixed.

DOMAINS:
  genetics_*  — gene annotation, IDs, sequences (MyGene, Ensembl)
  targets_*   — disease-gene associations, pharmacogenetics (Open Targets)
  drugs_*     — compounds, mechanisms, labels, safety (ChEMBL, OpenFDA)
  proteins_*  — structure, function, interactions (UniProt, AlphaFold, STRING)
  pathways_*  — metabolic & signaling pathways (KEGG, Reactome)
  variants_*  — genetic variant annotations, ClinVar (MyVariant.info)
  clinical_*  — trials, ICD-10 codes, providers (ClinicalTrials.gov, NPI)

QUICK REFERENCE:
  Gene lookup:        genetics_gene_info, genetics_ensembl_gene
  Gene→disease:       targets_disease_associations
  Gene→pathways:      pathways_kegg_gene, pathways_reactome_gene
  Gene→interactions:  proteins_interactions
  Variant annotation: variants_lookup, variants_clinvar
  Drug info:          drugs_compound, drugs_mechanism, drugs_label
  Protein structure:  proteins_structure (AlphaFold)
  Clinical trials:    clinical_trial_search
  Enrichment:         proteins_enrichment (STRING), pathways_reactome_enrichment
"""


def create_mcp(data_dir: Path | None = None) -> FastMCP:
    data_dir = data_dir or Path(os.environ.get("BIOMEDICAL_MCP_DATA", DEFAULT_DATA_DIR))
    data_dir.mkdir(parents=True, exist_ok=True)

    cache = Cache(data_dir / "cache.db")
    log.info("biomedical-mcp started (cache: %s)", data_dir / "cache.db")

    main = FastMCP(
        "biomedical",
        instructions=INSTRUCTIONS,
        middleware=[TelemetryMiddleware()],
    )

    main.mount(genetics.create_server(cache), namespace="genetics")
    main.mount(targets.create_server(cache), namespace="targets")
    main.mount(drugs.create_server(cache), namespace="drugs")
    main.mount(proteins.create_server(cache), namespace="proteins")
    main.mount(pathways.create_server(cache), namespace="pathways")
    main.mount(variants.create_server(cache), namespace="variants")
    main.mount(clinical.create_server(cache), namespace="clinical")

    return main


def main():
    mcp = create_mcp()
    mcp.run()
