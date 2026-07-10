"""
ARIA Memory System
------------------
Hierarchical persistent memory for omics experiments.

Structure:
  Wing   → Experiment (one per project)
  Hall   → Modality (RNA, ATAC, HiC, ChIP, CUT&RUN, CUT&TAG)
  Room   → Analysis type (QC, differential, integration, etc.)
  Closet → Compressed findings (~170 tokens on startup)
  Tunnel → Cross-modal connections (RNA <-> ATAC <-> HiC links)

Storage: SQLite (local, no cloud, no external services)
Loading: Progressive (L0 -> L3), only what is needed per query
"""

from __future__ import annotations
import os
import sqlite3
import json
import tempfile
import threading
from datetime import datetime
from pathlib import Path
from typing import Optional


SCHEMA = """
CREATE TABLE IF NOT EXISTS wings (
    id         TEXT PRIMARY KEY,
    name       TEXT NOT NULL,
    organism   TEXT,
    genome     TEXT,
    created_at TEXT,
    updated_at TEXT,
    summary    TEXT
);

CREATE TABLE IF NOT EXISTS halls (
    id       TEXT PRIMARY KEY,
    wing_id  TEXT NOT NULL,
    modality TEXT NOT NULL,
    status   TEXT DEFAULT 'pending',
    FOREIGN KEY (wing_id) REFERENCES wings(id)
);

CREATE TABLE IF NOT EXISTS rooms (
    id         TEXT PRIMARY KEY,
    hall_id    TEXT NOT NULL,
    analysis   TEXT NOT NULL,
    params     TEXT,
    tool       TEXT,
    created_at TEXT,
    FOREIGN KEY (hall_id) REFERENCES halls(id)
);

CREATE TABLE IF NOT EXISTS findings (
    id          TEXT PRIMARY KEY,
    room_id     TEXT NOT NULL,
    content     TEXT NOT NULL,
    confidence  TEXT NOT NULL,
    is_stale    INTEGER DEFAULT 0,
    valid_from  TEXT,
    valid_until TEXT,
    FOREIGN KEY (room_id) REFERENCES rooms(id)
);

CREATE TABLE IF NOT EXISTS tunnels (
    id          TEXT PRIMARY KEY,
    wing_id     TEXT NOT NULL,
    from_hall   TEXT NOT NULL,
    to_hall     TEXT NOT NULL,
    entity      TEXT NOT NULL,
    description TEXT,
    confidence  TEXT,
    created_at  TEXT,
    FOREIGN KEY (wing_id) REFERENCES wings(id)
);

CREATE TABLE IF NOT EXISTS decisions (
    id         TEXT PRIMARY KEY,
    wing_id    TEXT NOT NULL,
    checkpoint TEXT,
    question   TEXT,
    decision   TEXT,
    rationale  TEXT,
    made_at    TEXT,
    made_by    TEXT DEFAULT 'user',
    FOREIGN KEY (wing_id) REFERENCES wings(id)
);
"""


class ARIAMemory:
    """
    Persistent memory for a single ARIA installation.
    One SQLite file per lab, multiple wings (experiments).
    """

    def __init__(self, db_path: str = "~/.aria/memory.db"):
        if db_path == ":memory:":
            self.db_path = db_path
        else:
            self.db_path = str(Path(db_path).expanduser())
            Path(self.db_path).parent.mkdir(parents=True, exist_ok=True)

        # Agents run on background threads (Orchestrator._dispatch_agents);
        # share one connection with check_same_thread=False and serialize
        # every read/write through a single lock so concurrent commits cannot
        # interleave and corrupt the WAL.
        self._conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self._lock = threading.RLock()

        # WAL allows concurrent readers and a single writer without "database
        # is locked" errors; safe to enable on shared, file-backed DBs.
        if db_path != ":memory:":
            try:
                self._conn.execute("PRAGMA journal_mode=WAL")
                self._conn.execute("PRAGMA synchronous=NORMAL")
            except sqlite3.DatabaseError:
                pass

        self._conn.row_factory = sqlite3.Row
        with self._lock:
            self._conn.executescript(SCHEMA)
            self._migrate_decisions_checkpoint_to_text()
            self._conn.commit()

    # ── Internal helpers (all DB access funnels through these) ──────────

    def _migrate_decisions_checkpoint_to_text(self) -> None:
        columns = self._conn.execute("PRAGMA table_info(decisions)").fetchall()
        checkpoint = next((col for col in columns if col[1] == "checkpoint"),
                          None)
        if checkpoint is None or str(checkpoint[2]).upper() == "TEXT":
            return
        self._conn.executescript(
            """
            ALTER TABLE decisions RENAME TO decisions_old;
            CREATE TABLE decisions (
                id         TEXT PRIMARY KEY,
                wing_id    TEXT NOT NULL,
                checkpoint TEXT,
                question   TEXT,
                decision   TEXT,
                rationale  TEXT,
                made_at    TEXT,
                made_by    TEXT DEFAULT 'user',
                FOREIGN KEY (wing_id) REFERENCES wings(id)
            );
            INSERT OR REPLACE INTO decisions
            SELECT id, wing_id, CAST(checkpoint AS TEXT), question, decision,
                   rationale, made_at, made_by
            FROM decisions_old;
            DROP TABLE decisions_old;
            """
        )

    def _write(self, sql: str, params: tuple = ()):
        with self._lock:
            self._conn.execute(sql, params)
            self._conn.commit()

    def _fetchone(self, sql: str, params: tuple = ()):
        with self._lock:
            return self._conn.execute(sql, params).fetchone()

    def _fetchall(self, sql: str, params: tuple = ()):
        with self._lock:
            return self._conn.execute(sql, params).fetchall()

    # ── Minimized per-experiment export ──────────────────────────────────

    def _copy_scoped_rows(self, dest, table: str,
                          select_sql: str, params: tuple) -> int:
        """Copy the rows returned by ``select_sql`` from the live DB into the
        same-named table in ``dest``, matching columns by NAME (robust to any
        column-order drift). Read stays inside the caller's ``self._lock``."""
        cols = [c[1] for c in dest.execute(f"PRAGMA table_info({table})").fetchall()]
        rows = self._conn.execute(select_sql, params).fetchall()
        if not rows:
            return 0
        data = [tuple(r[c] for c in cols) for r in rows]
        placeholders = ",".join("?" * len(cols))
        dest.executemany(
            f"INSERT INTO {table} ({','.join(cols)}) VALUES ({placeholders})", data)
        return len(data)

    def export_experiment_snapshot(self, experiment_id: str, dest_path) -> dict:
        """Write a MINIMIZED SQLite holding ONLY ``experiment_id``'s wing subtree
        (the wing + its halls/rooms/findings + its tunnels/decisions) to
        ``dest_path``. It NEVER copies the whole lab DB, so a capsule/report cannot
        leak other experiments' state. The read runs inside the shared lock so a
        pending WAL write cannot tear the snapshot. The completed SQLite is checked
        and atomically published, so a failed export cannot replace a prior valid
        snapshot. Returns per-table row counts."""
        dest = Path(dest_path)
        dest.parent.mkdir(parents=True, exist_ok=True)
        handle = tempfile.NamedTemporaryFile(
            prefix=f".{dest.name}.", suffix=".tmp", dir=dest.parent, delete=False
        )
        tmp = Path(handle.name)
        handle.close()
        try:
            out = sqlite3.connect(str(tmp))
            try:
                out.execute("PRAGMA foreign_keys=ON")
                out.executescript(SCHEMA)
                counts: dict[str, int] = {}
                with self._lock:
                    hall_ids = [r["id"] for r in self._conn.execute(
                        "SELECT id FROM halls WHERE wing_id=?", (experiment_id,))]
                    room_ids: list = []
                    if hall_ids:
                        ph = ",".join("?" * len(hall_ids))
                        room_ids = [r["id"] for r in self._conn.execute(
                            f"SELECT id FROM rooms WHERE hall_id IN ({ph})",
                            tuple(hall_ids))]
                    counts["wings"] = self._copy_scoped_rows(
                        out, "wings", "SELECT * FROM wings WHERE id=?", (experiment_id,))
                    counts["halls"] = self._copy_scoped_rows(
                        out, "halls", "SELECT * FROM halls WHERE wing_id=?",
                        (experiment_id,))
                    counts["tunnels"] = self._copy_scoped_rows(
                        out, "tunnels", "SELECT * FROM tunnels WHERE wing_id=?",
                        (experiment_id,))
                    counts["decisions"] = self._copy_scoped_rows(
                        out, "decisions", "SELECT * FROM decisions WHERE wing_id=?",
                        (experiment_id,))
                    counts["rooms"] = (self._copy_scoped_rows(
                        out, "rooms",
                        f"SELECT * FROM rooms WHERE hall_id IN "
                        f"({','.join('?' * len(hall_ids))})", tuple(hall_ids))
                        if hall_ids else 0)
                    counts["findings"] = (self._copy_scoped_rows(
                        out, "findings",
                        f"SELECT * FROM findings WHERE room_id IN "
                        f"({','.join('?' * len(room_ids))})", tuple(room_ids))
                        if room_ids else 0)
                out.commit()
                integrity = out.execute("PRAGMA integrity_check").fetchone()[0]
                foreign_key_errors = out.execute("PRAGMA foreign_key_check").fetchall()
                if integrity != "ok" or foreign_key_errors:
                    raise sqlite3.DatabaseError(
                        "scoped snapshot failed SQLite integrity validation"
                    )
                wing_ids = {
                    row[0] for row in out.execute("SELECT id FROM wings").fetchall()
                }
                if wing_ids not in ({experiment_id}, set()):
                    raise sqlite3.DatabaseError(
                        "scoped snapshot contains an unexpected experiment"
                    )
            finally:
                out.close()
            with tmp.open("rb") as stream:
                os.fsync(stream.fileno())
            os.replace(tmp, dest)
            try:
                dir_fd = os.open(dest.parent, os.O_RDONLY)
                try:
                    os.fsync(dir_fd)
                finally:
                    os.close(dir_fd)
            except OSError:
                pass
        except Exception:
            tmp.unlink(missing_ok=True)
            raise
        return {
            "scope": experiment_id,
            "tables": counts,
            "integrity_check": integrity,
        }

    # ── Wings ────────────────────────────────────────────────────────────

    def create_wing(self, experiment_id: str, name: str,
                    organism: str = None, genome: str = None) -> str:
        now = datetime.now().isoformat()
        self._write(
            "INSERT OR REPLACE INTO wings VALUES (?,?,?,?,?,?,?)",
            (experiment_id, name, organism, genome, now, now, None),
        )
        return experiment_id

    def update_wing_summary(self, experiment_id: str, summary: str):
        self._write(
            "UPDATE wings SET summary=?, updated_at=? WHERE id=?",
            (summary, datetime.now().isoformat(), experiment_id),
        )

    def get_wing(self, experiment_id: str) -> Optional[dict]:
        row = self._fetchone("SELECT * FROM wings WHERE id=?", (experiment_id,))
        return dict(row) if row else None

    def list_wings(self) -> list[dict]:
        rows = self._fetchall(
            "SELECT id, name, organism, genome, updated_at, summary FROM wings"
        )
        return [dict(r) for r in rows]

    # ── Halls ────────────────────────────────────────────────────────────

    def create_hall(self, hall_id: str, wing_id: str, modality: str) -> str:
        self._write(
            "INSERT OR REPLACE INTO halls VALUES (?,?,?,?)",
            (hall_id, wing_id, modality, "pending"),
        )
        return hall_id

    def update_hall_status(self, hall_id: str, status: str):
        self._write("UPDATE halls SET status=? WHERE id=?", (status, hall_id))

    def get_halls(self, wing_id: str) -> list[dict]:
        rows = self._fetchall("SELECT * FROM halls WHERE wing_id=?", (wing_id,))
        return [dict(r) for r in rows]

    # ── Rooms ────────────────────────────────────────────────────────────

    def create_room(self, room_id: str, hall_id: str, analysis: str,
                    params: dict = None, tool: str = None) -> str:
        self._write(
            "INSERT OR REPLACE INTO rooms VALUES (?,?,?,?,?,?)",
            (room_id, hall_id, analysis,
             json.dumps(params or {}), tool,
             datetime.now().isoformat()),
        )
        return room_id

    # ── Findings ─────────────────────────────────────────────────────────

    def store_finding(self, finding_id: str, room_id: str,
                      content: str, confidence: str) -> str:
        self._write(
            "INSERT OR REPLACE INTO findings VALUES (?,?,?,?,?,?,?)",
            (finding_id, room_id, content, confidence, 0,
             datetime.now().isoformat(), None),
        )
        return finding_id

    def invalidate_finding(self, finding_id: str):
        self._write(
            "UPDATE findings SET is_stale=1, valid_until=? WHERE id=?",
            (datetime.now().isoformat(), finding_id),
        )

    def get_active_findings(self, room_id: str) -> list[dict]:
        rows = self._fetchall(
            "SELECT * FROM findings WHERE room_id=? AND is_stale=0",
            (room_id,),
        )
        return [dict(r) for r in rows]

    # ── Tunnels ──────────────────────────────────────────────────────────

    def create_tunnel(self, tunnel_id: str, wing_id: str,
                      from_hall: str, to_hall: str,
                      entity: str, description: str,
                      confidence: str = "medium"):
        self._write(
            "INSERT OR REPLACE INTO tunnels VALUES (?,?,?,?,?,?,?,?)",
            (tunnel_id, wing_id, from_hall, to_hall,
             entity, description, confidence,
             datetime.now().isoformat()),
        )

    def get_tunnels(self, wing_id: str, entity: str = None) -> list[dict]:
        if entity:
            rows = self._fetchall(
                "SELECT * FROM tunnels WHERE wing_id=? AND entity=?",
                (wing_id, entity),
            )
        else:
            rows = self._fetchall(
                "SELECT * FROM tunnels WHERE wing_id=?", (wing_id,)
            )
        return [dict(r) for r in rows]

    # ── Decisions ────────────────────────────────────────────────────────

    def store_decision(self, decision_id: str, wing_id: str,
                       checkpoint: str | int | float, question: str,
                       decision: str, rationale: str = "",
                       made_by: str = "user"):
        self._write(
            "INSERT OR REPLACE INTO decisions VALUES (?,?,?,?,?,?,?,?)",
            (decision_id, wing_id, str(checkpoint), question,
             decision, rationale, datetime.now().isoformat(), made_by),
        )

    def get_decisions(self, wing_id: str) -> list[dict]:
        rows = self._fetchall(
            "SELECT * FROM decisions WHERE wing_id=? ORDER BY made_at",
            (wing_id,),
        )
        return [dict(r) for r in rows]

    # ── L0 Startup context ───────────────────────────────────────────────

    def startup_context(self) -> str:
        """
        Minimal context loaded on every ARIA startup.
        Target: ~170 tokens. Just enough to orient the system.
        """
        wings = self.list_wings()
        if not wings:
            return "No experiments found. ARIA ready for first analysis."
        lines = ["ARIA Memory -- Active Experiments:"]
        for w in wings[-5:]:
            summary = w.get("summary") or "No summary yet"
            lines.append(
                f"  [{w['id'][:8]}] {w['name']} "
                f"({w.get('organism','?')}/{w.get('genome','?')}) "
                f"-- {summary}"
            )
        return "\n".join(lines)

    def close(self):
        with self._lock:
            self._conn.close()
