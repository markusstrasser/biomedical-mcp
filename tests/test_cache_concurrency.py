"""Regression: Cache must survive concurrent access from the composite thread pool.

Pre-fix (one sqlite connection, no lock) this raised
'sqlite3.InterfaceError: bad parameter or other API misuse' under thread races —
caught by the live smoke on gene_dossier's concurrent clingen calls (2026-06-07).
"""

import threading

from biomedical_mcp.cache import Cache


def test_concurrent_get_set_no_interface_error(tmp_path):
    cache = Cache(tmp_path / "c.db")
    errors = []

    def hammer(tid: int):
        try:
            for i in range(50):
                cache.set(f"k{tid}-{i}", {"v": i})
                cache.get(f"k{tid}-{i}", max_age_days=1)
        except Exception as exc:  # noqa: BLE001
            errors.append(exc)

    threads = [threading.Thread(target=hammer, args=(t,)) for t in range(16)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert not errors, f"concurrent cache access raised: {errors[:3]}"
    # spot-check a value round-tripped
    assert cache.get("k0-0", max_age_days=1) == {"v": 0}
