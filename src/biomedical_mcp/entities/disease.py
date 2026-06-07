"""Disease composite entity — STUB (filled in Phase 2). Conforms to entities/gene.py."""

from __future__ import annotations

from typing import Any

from biomedical_mcp.composite_core import Section, Fetcher

ENTITY = "disease"
TOOL_NAME = "disease_profile"
DESCRIPTION = "Disease profile (stub). identifier = disease name, MONDO/ORPHA/OMIM id."
NEEDS: tuple[str, ...] = ()
SECTIONS: tuple[Section, ...] = ()


def resolve(identifier: str, clients: dict[str, Any]) -> tuple[str, dict | None]:
    return identifier, None


def fetchers(identifier: str, clients: dict[str, Any], **opts: Any) -> dict[str, Fetcher]:
    return {}
