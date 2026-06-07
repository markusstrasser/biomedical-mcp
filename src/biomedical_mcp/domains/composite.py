"""Composite domain — the consolidated entity-composite tool surface.

Registers ONE tool per entity module in biomedical_mcp.entities (gene_dossier,
variant_context, drug_profile, protein_profile, disease_profile), plus the
describe_sections discovery tool and an admin health_check. Each composite returns
the standardized partial-failure envelope from composite_core.

This replaces ~85 raw per-source tools as the default surface. The raw domain tools
remain available behind the `full` server profile (see server.py). Adding an entity
or a source touches only entities/<entity>.py — never this file.
"""

from __future__ import annotations

import logging

from fastmcp import FastMCP

from biomedical_mcp.cache import Cache
from biomedical_mcp.composite_core import run_sections, describe
from biomedical_mcp import entities
from biomedical_mcp.mygene import MyGene
from biomedical_mcp.opentargets import OpenTargets
from biomedical_mcp.uniprot import UniProt
from biomedical_mcp.monarch import Monarch

log = logging.getLogger(__name__)


def create_server(cache: Cache) -> FastMCP:
    clients = entities.build_clients(cache)
    # Search clients for bio_search (name/keyword → identifier resolution).
    search_clients = {
        "mygene": MyGene(cache), "opentargets": OpenTargets(cache),
        "uniprot": UniProt(cache), "monarch": Monarch(cache),
    }
    server = FastMCP("composite")

    def _make(mod):
        def composite_tool(identifier: str, sections: list[str] | None = None) -> dict:
            resolved, short_circuit = mod.resolve(identifier, clients)
            if short_circuit is not None:
                return short_circuit
            fetch = mod.fetchers(resolved, clients)
            extra = {"input": identifier} if resolved != identifier else None
            return run_sections(mod.ENTITY, resolved, sections, fetch, mod.SECTIONS, extra=extra)

        composite_tool.__name__ = mod.TOOL_NAME
        composite_tool.__doc__ = mod.DESCRIPTION
        return composite_tool

    for mod in entities.ENTITY_MODULES:
        if not mod.SECTIONS:
            log.debug("skipping unfilled entity composite: %s", mod.ENTITY)
            continue
        server.tool(name=mod.TOOL_NAME, tags={mod.ENTITY})(_make(mod))

    @server.tool(tags={"discovery"})
    def describe_sections(entity: str) -> dict:
        """List available sections for a composite entity, with sources and
        descriptions. Use before calling an entity composite to discover valid
        `sections` values.

        Args:
            entity: One of the registered entities (gene, variant, drug, protein, disease).
        """
        idx = entities.sections_index()
        meta = idx.get(entity)
        if not meta:
            return {
                "error": f"unknown or unfilled entity '{entity}'",
                "valid_entities": [e for e, m in idx.items() if m],
            }
        return describe(entity, meta)

    @server.tool(tags={"discovery"})
    def bio_search(entity: str, query: str, limit: int = 10) -> dict:
        """Resolve a name/keyword to candidate identifiers, then call the matching
        composite (gene_dossier, drug_profile, etc.) with the returned id.

        Args:
            entity: One of "gene", "drug", "protein", "disease".
            query: Free-text name or keyword (e.g. "cytochrome P450", "imatinib").
            limit: Max results (default 10).
        """
        e = entity.lower()
        try:
            if e == "gene":
                hits = search_clients["mygene"].search(query, limit=limit)
            elif e == "drug":
                hits = search_clients["opentargets"].search(query, entity_type="drug", limit=limit)
            elif e == "protein":
                hits = search_clients["uniprot"].search(query, limit=limit)
            elif e == "disease":
                hits = search_clients["monarch"].search(query, category="disease", limit=limit)
            else:
                return {"error": f"search not supported for entity '{entity}'",
                        "supported": ["gene", "drug", "protein", "disease"]}
            return {"entity": e, "query": query, "results": hits,
                    "next": f"call composite_{e}_dossier / composite_{e}_profile with a result id"}
        except Exception as exc:  # noqa: BLE001
            return {"error": f"{type(exc).__name__}: {exc}", "entity": e, "query": query}

    @server.tool(tags={"admin"})
    def health_check() -> dict:
        """Check connectivity to all upstream APIs backing the composites."""
        statuses: dict[str, str] = {}
        for name, client in clients.items():
            try:
                statuses[name] = "ok" if client.validate() else "error"
            except Exception as exc:  # noqa: BLE001
                statuses[name] = f"error: {type(exc).__name__}"
        ok_count = sum(1 for v in statuses.values() if v == "ok")
        return {"status": f"{ok_count}/{len(statuses)} ok", "sources": statuses}

    return server
