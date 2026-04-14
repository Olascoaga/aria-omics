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
import sqlite3
import json
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
    checkpoint INTEGER,
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
        self._conn = sqlite3.connect(self.db_path)
        self._conn.row_factory = sqlite3.Row
        self._conn.executescript(SCHEMA)
        self._conn.commit()

    # ── Wings ────────────────────────────────────────────────────────────

    def create_wing(self, experiment_id: str, name: str,
                    organism: str = None, genome: str = None) -> str:
        now = datetime.now().isoformat()
        self._conn.execute(
            "INSERT OR REPLACE INTO wings VALUES (?,?,?,?,?,?,?)",
            (experiment_id, name, organism, genome, now, now, None)
        )
        self._conn.commit()
        return experiment_id

    def update_wing_summary(self, experiment_id: str, summary: str):
        self._conn.execute(
            "UPDATE wings SET summary=?, updated_at=? WHERE id=?",
            (summary, datetime.now().isoformat(), experiment_id)
        )
        self._conn.commit()

    def get_wing(self, experiment_id: str) -> Optional[dict]:
        row = self._conn.execute(
            "SELECT * FROM wings WHERE id=?", (experiment_id,)
        ).fetchone()
        return dict(row) if row else None

    def list_wings(self) -> list[dict]:
        rows = self._conn.execute(
            "SELECT id, name, organism, genome, updated_at, summary FROM wings"
        ).fetchall()
        return [dict(r) for r in rows]

    # ── Halls ────────────────────────────────────────────────────────────

    def create_hall(self, hall_id: str, wing_id: str, modality: str) -> str:
        self._conn.execute(
            "INSERT OR REPLACE INTO halls VALUES (?,?,?,?)",
            (hall_id, wing_id, modality, "pending")
        )
        self._conn.commit()
        return hall_id

    def update_hall_status(self, hall_id: str, status: str):
        self._conn.execute(
            "UPDATE halls SET status=? WHERE id=?", (status, hall_id)
        )
        self._conn.commit()

    def get_halls(self, wing_id: str) -> list[dict]:
        rows = self._conn.execute(
            "SELECT * FROM halls WHERE wing_id=?", (wing_id,)
        ).fetchall()
        return [dict(r) for r in rows]

    # ── Rooms ────────────────────────────────────────────────────────────

    def create_room(self, room_id: str, hall_id: str, analysis: str,
                    params: dict = None, tool: str = None) -> str:
        self._conn.execute(
            "INSERT OR REPLACE INTO rooms VALUES (?,?,?,?,?,?)",
            (room_id, hall_id, analysis,
             json.dumps(params or {}), tool,
             datetime.now().isoformat())
        )
        self._conn.commit()
        return room_id

    # ── Findings ─────────────────────────────────────────────────────────

    def store_finding(self, finding_id: str, room_id: str,
                      content: str, confidence: str) -> str:
        self._conn.execute(
            "INSERT OR REPLACE INTO findings VALUES (?,?,?,?,?,?,?)",
            (finding_id, room_id, content, confidence, 0,
             datetime.now().isoformat(), None)
        )
        self._conn.commit()
        return finding_id

    def invalidate_finding(self, finding_id: str):
        self._conn.execute(
            "UPDATE findings SET is_stale=1, valid_until=? WHERE id=?",
            (datetime.now().isoformat(), finding_id)
        )
        self._conn.commit()

    def get_active_findings(self, room_id: str) -> list[dict]:
        rows = self._conn.execute(
            "SELECT * FROM findings WHERE room_id=? AND is_stale=0",
            (room_id,)
        ).fetchall()
        return [dict(r) for r in rows]

    # ── Tunnels ──────────────────────────────────────────────────────────

    def create_tunnel(self, tunnel_id: str, wing_id: str,
                      from_hall: str, to_hall: str,
                      entity: str, description: str,
                      confidence: str = "medium"):
        self._conn.execute(
            "INSERT OR REPLACE INTO tunnels VALUES (?,?,?,?,?,?,?,?)",
            (tunnel_id, wing_id, from_hall, to_hall,
             entity, description, confidence,
             datetime.now().isoformat())
        )
        self._conn.commit()

    def get_tunnels(self, wing_id: str, entity: str = None) -> list[dict]:
        if entity:
            rows = self._conn.execute(
                "SELECT * FROM tunnels WHERE wing_id=? AND entity=?",
                (wing_id, entity)
            ).fetchall()
        else:
            rows = self._conn.execute(
                "SELECT * FROM tunnels WHERE wing_id=?", (wing_id,)
            ).fetchall()
        return [dict(r) for r in rows]

    # ── Decisions ────────────────────────────────────────────────────────

    def store_decision(self, decision_id: str, wing_id: str,
                       checkpoint: int, question: str,
                       decision: str, rationale: str = "",
                       made_by: str = "user"):
        self._conn.execute(
            "INSERT OR REPLACE INTO decisions VALUES (?,?,?,?,?,?,?,?)",
            (decision_id, wing_id, checkpoint, question,
             decision, rationale, datetime.now().isoformat(), made_by)
        )
        self._conn.commit()

    def get_decisions(self, wing_id: str) -> list[dict]:
        rows = self._conn.execute(
            "SELECT * FROM decisions WHERE wing_id=? ORDER BY made_at",
            (wing_id,)
        ).fetchall()
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
        self._conn.close()
