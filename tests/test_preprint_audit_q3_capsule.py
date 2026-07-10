"""Preprint-readiness audit — Q3 (interim of blockers A2 + A6).

A2: the lab keeps ONE SQLite for all experiments (wings), so copying the whole DB
into the report/capsule leaked every other experiment's state. The memory snapshot
must now export ONLY the current experiment's wing subtree.

A6: the reproducible report dir name is deterministic, so a rerun reused the same
directory and left stale artifacts that the capsule then bundled. Each reproducible
build must start from a clean directory.

Tracker: memory/audit/ARIA_PLAN_AUDITORIA_preprint_journal_2026-07-09.md
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

from aria.memory.memory import ARIAMemory


def _seed(mem: ARIAMemory, exp: str) -> None:
    mem.create_wing(exp, name=f"{exp}-name", organism="human", genome="hg38")
    hall = f"{exp}-hall"
    mem.create_hall(hall, exp, "bulk_RNA")
    room = f"{exp}-room"
    mem.create_room(room, hall, "differential_expression")
    mem.store_finding(f"{exp}-find", room, content="DE result", confidence="high")
    mem.create_tunnel(f"{exp}-tun", exp, hall, hall, entity="GENE1", description="d")
    mem.store_decision(f"{exp}-dec", exp, checkpoint=1, question="q?", decision="yes")


def test_export_experiment_snapshot_is_scoped_to_one_experiment(tmp_path):
    mem = ARIAMemory(db_path=str(tmp_path / "lab.db"))
    _seed(mem, "expA")
    _seed(mem, "expB")

    dest = tmp_path / "snapshot.sqlite"
    result = mem.export_experiment_snapshot("expA", dest)
    assert result["scope"] == "expA"
    assert result["tables"]["wings"] == 1  # only expA's wing

    conn = sqlite3.connect(str(dest))
    try:
        # Valid SQLite, not a torn copy.
        assert conn.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
        # Only expA rows exist across every scoped table; expB never leaks.
        for table, col in (("wings", "id"), ("halls", "wing_id"),
                           ("tunnels", "wing_id"), ("decisions", "wing_id")):
            vals = {r[0] for r in conn.execute(f"SELECT {col} FROM {table}")}
            assert vals == {"expA"}, f"{table} leaked: {vals}"
        # Subtree (rooms/findings) came along and is non-empty.
        assert conn.execute("SELECT COUNT(*) FROM rooms").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM findings").fetchone()[0] == 1
    finally:
        conn.close()


def test_snapshot_of_unknown_experiment_is_empty_not_full_db(tmp_path):
    mem = ARIAMemory(db_path=str(tmp_path / "lab.db"))
    _seed(mem, "expA")
    _seed(mem, "expB")

    dest = tmp_path / "snapshot.sqlite"
    result = mem.export_experiment_snapshot("does-not-exist", dest)
    assert all(v == 0 for v in result["tables"].values())
    conn = sqlite3.connect(str(dest))
    try:
        # No other experiment's rows are present.
        assert conn.execute("SELECT COUNT(*) FROM wings").fetchone()[0] == 0
    finally:
        conn.close()


def test_write_memory_snapshot_never_copies_full_db(tmp_path):
    # The report-builder helper must call the scoped export, never a whole-DB copy.
    from aria.agents.narrative.report_builder import ReportBuilderMixin

    mem = ARIAMemory(db_path=str(tmp_path / "lab.db"))
    _seed(mem, "expA")
    _seed(mem, "expB")

    class _Host(ReportBuilderMixin):
        def __init__(self, memory):
            self.memory = memory

    report_dir = tmp_path / "report"
    report_dir.mkdir()
    _Host(mem)._write_memory_snapshot(report_dir, "expA")

    snap = report_dir / "memory_snapshot.sqlite"
    assert snap.exists()
    conn = sqlite3.connect(str(snap))
    try:
        wings = {r[0] for r in conn.execute("SELECT id FROM wings")}
        assert wings == {"expA"}  # expB absent -> no full-DB copy
    finally:
        conn.close()


def test_reproducible_report_dir_is_wiped_clean_on_rerun(tmp_path):
    from aria.agents.narrative.report_builder import ReportBuilderMixin

    class _Host(ReportBuilderMixin):
        def __init__(self, reports_dir):
            self.reports_dir = Path(reports_dir)

        def _build_slug(self, intent, exp_ctx):
            return "slug"

    host = _Host(tmp_path)
    exp_ctx = {"reproducible_mode": True,
               "input_files": [{"sha256": "abc123def4567890"}]}

    d1 = host._build_report_dir("exp1234", {}, exp_ctx)
    stale = d1 / "figures" / "stale_from_prior_run.png"
    stale.write_bytes(b"old")
    assert stale.exists()

    d2 = host._build_report_dir("exp1234", {}, exp_ctx)  # deterministic same name
    assert d2 == d1                      # same deterministic dir
    assert not stale.exists()            # ... but wiped clean, no residue
    assert (d2 / "figures").is_dir() and (d2 / "tables").is_dir()
