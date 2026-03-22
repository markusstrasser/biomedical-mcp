"""SQLite response cache with TTL."""

import json
import sqlite3
from pathlib import Path

SCHEMA = """
CREATE TABLE IF NOT EXISTS cache (
    cache_key TEXT PRIMARY KEY,
    response TEXT NOT NULL,
    cached_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
"""


class Cache:
    def __init__(self, db_path: Path):
        self.db_path = db_path
        db_path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(str(db_path), check_same_thread=False)
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA busy_timeout=5000")
        self.conn.executescript(SCHEMA)

    def get(self, key: str, max_age_days: int = 7) -> dict | list | None:
        row = self.conn.execute(
            """SELECT response FROM cache
               WHERE cache_key = ? AND cached_at > datetime('now', ?)""",
            (key, f"-{max_age_days} days"),
        ).fetchone()
        return json.loads(row[0]) if row else None

    def set(self, key: str, value) -> None:
        self.conn.execute(
            "INSERT OR REPLACE INTO cache (cache_key, response) VALUES (?, ?)",
            (key, json.dumps(value)),
        )
        self.conn.commit()

    def cleanup(self, max_age_days: int = 90) -> int:
        """Delete expired cache entries older than max_age_days. Returns count deleted."""
        cursor = self.conn.execute(
            "DELETE FROM cache WHERE cached_at < datetime('now', ?)",
            (f"-{max_age_days} days",),
        )
        deleted = cursor.rowcount
        if deleted > 0:
            self.conn.execute("PRAGMA optimize")
        self.conn.commit()
        return deleted
