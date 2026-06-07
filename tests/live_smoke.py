"""Live end-to-end smoke for the entity composites — runs against REAL APIs.

WHY THIS EXISTS (read before deleting):
The offline suite uses fake clients and is structurally blind to integration,
concurrency, and real response-schema drift. On 2026-06-08 the offline suite (145
tests) AND a cross-model `/critique close` both passed while THIS smoke immediately
surfaced 3 real bugs: a sqlite-cache concurrency crash under the section ThreadPool,
an OpenTargets `{"disease": null}` AttributeError, and ClinPGx raising on non-PGx
genes. Offline-green + model-clean is NECESSARY, not SUFFICIENT, for this server.

WHEN TO RUN: after touching any client (src/biomedical_mcp/*.py), the section runner
(composite_core.py), or an entity module (entities/*.py). Before claiming integration
code works.

    uv run python3 tests/live_smoke.py            # default profile (composite)
    ONCOKB_TOKEN=... uv run python3 tests/live_smoke.py   # also exercises OncoKB

This is a DIAGNOSTIC, not a gate: it always exits 0 and never hard-fails on a transient
upstream outage (that would be noise). It prints a per-section status table; a human/agent
reads it. `error` rows that aren't transient are the signal.
"""

import os
import pathlib
import tempfile
import time

from biomedical_mcp.cache import Cache
from biomedical_mcp import entities
from biomedical_mcp.composite_core import run_sections


def _run(entity, ident, clients, mods, sections=None):
    mod = mods[entity]
    t0 = time.monotonic()
    resolved, short = mod.resolve(ident, clients)
    if short is not None:
        print(f"\n[{entity}] {ident!r} -> SHORT-CIRCUIT: {short.get('status')}")
        return
    env = run_sections(entity, resolved, sections, mod.fetchers(resolved, clients),
                       mod.SECTIONS, timeout=25)
    dt = time.monotonic() - t0
    print(f"\n[{entity}] {ident!r} -> overall={env['overall_status']} ({dt:.1f}s)")
    for name, s in env["sections"].items():
        st = s["status"]
        extra = ""
        if st == "ok":
            d = s["data"]
            extra = f"({len(d) if isinstance(d, (list, dict)) else 1} keys/items)"
        elif st == "error":
            extra = f"-> {s['error']['type']}: {s['error']['message'][:70]}"
        print(f"    {st:13s} {name:24s} {extra}")


def main() -> None:
    cache = Cache(pathlib.Path(tempfile.mkdtemp()) / "live_smoke.db")
    clients = entities.build_clients(cache)
    mods = {m.ENTITY: m for m in entities.ENTITY_MODULES}

    print("=" * 70)
    print("biomedical-mcp LIVE smoke — real APIs. 'error' rows (non-transient) = signal.")
    print("=" * 70)
    _run("gene", "BRCA1", clients, mods)
    _run("drug", "warfarin", clients, mods,
         sections=["compound", "label", "pharmacogenomics", "interactions"])
    _run("variant", "rs1801133", clients, mods)
    _run("variant", "BRAF V600E", clients, mods, sections=["annotation", "somatic"])
    _run("protein", "BRCA1", clients, mods)
    _run("disease", "cystic fibrosis", clients, mods)
    _run("disease", "ORPHA:586", clients, mods)

    if not os.environ.get("ONCOKB_TOKEN"):
        print("\n[note] ONCOKB_TOKEN unset — somatic OncoKB half returns its api-key "
              "error dict by design (CIViC half still runs).")
    print("\nDone. Review 'error' rows above; transient upstream outages are not bugs.")


if __name__ == "__main__":
    main()
