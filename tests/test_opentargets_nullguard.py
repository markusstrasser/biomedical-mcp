"""Regression: Open Targets GraphQL returns {"disease": null} / {"target": null}
for ids it doesn't know. data.get("disease", {}) returns None (default only applies
to a MISSING key), so the next .get() crashed with AttributeError — caught by the
live smoke on disease_profile("ORPHA:586") (2026-06-07).
"""

from biomedical_mcp.cache import Cache
from biomedical_mcp.opentargets import OpenTargets


def _ot(tmp_path, payload):
    ot = OpenTargets(Cache(tmp_path / "c.db"))
    ot._cached = lambda *a, **k: payload  # type: ignore
    return ot


def test_disease_targets_null_disease(tmp_path):
    ot = _ot(tmp_path, {"disease": None})
    out = ot.disease_targets("Orphanet_586")
    assert out["total"] == 0
    assert out["targets"] == []


def test_disease_associations_null_target(tmp_path):
    ot = _ot(tmp_path, {"target": None})
    out = ot.disease_associations("ENSG00000000000")
    assert out["total"] == 0
    assert out["associations"] == []


def test_pharmacogenetics_null_target(tmp_path):
    ot = _ot(tmp_path, {"target": None})
    out = ot.pharmacogenetics("ENSG00000000000")
    assert out  # does not raise; returns a dict
