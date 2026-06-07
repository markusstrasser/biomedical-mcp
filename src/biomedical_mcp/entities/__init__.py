"""Entity composite registry — assembles the consolidated tool surface.

ENTITY_MODULES is the list of entity composites exposed by default. Each module
conforms to the interface documented in entities/gene.py. domains/composite.py
registers one tool per module and builds the union of needed clients via CLIENT_FACTORY.

Add an entity: write entities/<name>.py to the interface, import it here, append to
ENTITY_MODULES. Add a source to an existing entity: edit that module only.
"""

from __future__ import annotations

from typing import Any, Callable

from biomedical_mcp.cache import Cache

# client key -> zero-arg-after-cache constructor. Lazy imports keep startup light
# and let the union-of-needs build only the clients actually used.
from biomedical_mcp.hgnc import HGNC
from biomedical_mcp.gnomad import GnomAD
from biomedical_mcp.panelapp import PanelApp
from biomedical_mcp.gtex import GTEx
from biomedical_mcp.hpo import HPO
from biomedical_mcp.clingen import ClinGen
from biomedical_mcp.orphanet import Orphanet
from biomedical_mcp.myvariant import MyVariant
from biomedical_mcp.litvar import LitVar
from biomedical_mcp.gwas_catalog import GWASCatalog
from biomedical_mcp.chembl import ChEMBL
from biomedical_mcp.openfda import OpenFDA
from biomedical_mcp.opentargets import OpenTargets
from biomedical_mcp.pharmvar import PharmVar
from biomedical_mcp.uniprot import UniProt
from biomedical_mcp.alphafold import AlphaFold
from biomedical_mcp.stringdb import StringDB
from biomedical_mcp.interpro import InterPro
from biomedical_mcp.pdb import PDB
from biomedical_mcp.monarch import Monarch
from biomedical_mcp.icd10 import ICD10
from biomedical_mcp.clinpgx import ClinPGx
from biomedical_mcp.ddinter import DDInter
from biomedical_mcp.somatic import Somatic

from biomedical_mcp.entities import gene, variant, drug, protein, disease

CLIENT_FACTORY: dict[str, Callable[[Cache], Any]] = {
    "hgnc": HGNC, "gnomad": GnomAD, "panelapp": PanelApp, "gtex": GTEx,
    "hpo": HPO, "clingen": ClinGen, "orphanet": Orphanet, "myvariant": MyVariant,
    "litvar": LitVar, "gwas_catalog": GWASCatalog, "chembl": ChEMBL,
    "openfda": OpenFDA, "opentargets": OpenTargets, "pharmvar": PharmVar,
    "uniprot": UniProt, "alphafold": AlphaFold, "stringdb": StringDB,
    "interpro": InterPro, "pdb": PDB, "monarch": Monarch, "icd10": ICD10,
    "clinpgx": ClinPGx, "ddinter": DDInter, "somatic": Somatic,
}

ENTITY_MODULES = [gene, variant, drug, protein, disease]


def build_clients(cache: Cache) -> dict[str, Any]:
    """Instantiate exactly the clients the active entity modules need (deduped)."""
    needed: set[str] = set()
    for mod in ENTITY_MODULES:
        needed.update(mod.NEEDS)
    return {key: CLIENT_FACTORY[key](cache) for key in needed}


def sections_index() -> dict[str, tuple]:
    """entity -> SECTIONS, for describe_sections."""
    return {mod.ENTITY: mod.SECTIONS for mod in ENTITY_MODULES}
