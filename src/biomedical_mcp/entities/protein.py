"""Protein composite entity — STUB (filled in Phase 2). Conforms to entities/gene.py."""

from __future__ import annotations

from typing import Any

from biomedical_mcp.composite_core import Section, Fetcher

ENTITY = "protein"
TOOL_NAME = "protein_profile"
DESCRIPTION = "Protein profile (stub). identifier = gene symbol or UniProt accession."
NEEDS: tuple[str, ...] = ()
SECTIONS: tuple[Section, ...] = ()


def resolve(identifier: str, clients: dict[str, Any]) -> tuple[str, dict | None]:
    return identifier, None


def fetchers(identifier: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    return {}
