"""Preprint-readiness audit A2: transactional, experiment-scoped capsules.

The public export must never copy the live lab database, leak another experiment,
publish a torn SQLite/ZIP, or serialize private absolute paths.  These guards cover
the full memory snapshot -> RO-Crate -> capsule -> verification path.
"""
from __future__ import annotations

import json
import sqlite3
import zipfile
from pathlib import Path

import pytest

from aria.agents.narrative.ledger_export import (
    build_ro_crate,
    verify_reproducible_capsule,
    write_reproducible_capsule,
)
from aria.memory.memory import ARIAMemory


def _seed(mem: ARIAMemory, experiment_id: str) -> None:
    mem.create_wing(experiment_id, f"{experiment_id}-name", "human", "hg38")
    hall_id = f"{experiment_id}-hall"
    room_id = f"{experiment_id}-room"
    mem.create_hall(hall_id, experiment_id, "bulk_RNA")
    mem.create_room(room_id, hall_id, "differential_expression")
    mem.store_finding(
        f"{experiment_id}-finding", room_id, "measured result", "high"
    )
    mem.store_decision(
        f"{experiment_id}-decision",
        experiment_id,
        checkpoint="1",
        question="Proceed?",
        decision="yes",
    )


def _methodology(private_input: Path) -> dict:
    return {
        "provenance": {
            "version": "4.7.0",
            "git_commit": "abc123",
            "git_dirty": False,
            "workflow_hash": "wf-abc123",
            "workspace": str(private_input.parent),
        },
        "inputs": [
            {
                "path": str(private_input),
                "sha256": "deadbeef",
                "bytes": 10,
                "modality": "bulk_RNA",
            }
        ],
        "claims": [],
        "run_ledger": {"entries": [], "n_divergences": 0},
    }


def _snapshot_from_zip(capsule: Path, member: str, dest: Path) -> Path:
    with zipfile.ZipFile(capsule) as archive:
        dest.write_bytes(archive.read(member))
    return dest


def test_scoped_snapshot_publish_is_atomic_on_export_failure(tmp_path, monkeypatch):
    mem = ARIAMemory(str(tmp_path / "lab.sqlite"))
    _seed(mem, "expA")
    dest = tmp_path / "snapshot.sqlite"
    dest.write_bytes(b"previous-good-snapshot")

    def _fail(*_args, **_kwargs):
        raise sqlite3.OperationalError("injected copy failure")

    monkeypatch.setattr(mem, "_copy_scoped_rows", _fail)

    with pytest.raises(sqlite3.OperationalError, match="injected copy failure"):
        mem.export_experiment_snapshot("expA", dest)

    assert dest.read_bytes() == b"previous-good-snapshot"
    assert not list(tmp_path.glob(".snapshot.sqlite.*.tmp"))


def test_ro_crate_replaces_private_absolute_input_paths():
    private = Path("/private/lab/patient-001/counts.tsv")

    crate = build_ro_crate(_methodology(private))
    blob = json.dumps(crate, sort_keys=True)

    assert str(private) not in blob
    assert str(private.parent) not in blob
    assert "input://sha256/deadbeef" in blob


def test_capsule_is_scoped_portable_and_sqlite_verified_end_to_end(tmp_path):
    mem = ARIAMemory(str(tmp_path / "lab.sqlite"))
    # Disable automatic WAL checkpointing so committed rows remain in the WAL
    # while the scoped export is read from the live connection.
    with mem._lock:
        mem._conn.execute("PRAGMA wal_autocheckpoint=0")
    _seed(mem, "expA")
    _seed(mem, "expB")
    assert (tmp_path / "lab.sqlite-wal").exists()

    report_dir = tmp_path / "report-expA"
    report_dir.mkdir()
    private_input = Path("/private/lab/patient-001/counts.tsv")
    (report_dir / "methodology.json").write_text(
        json.dumps(_methodology(private_input)), encoding="utf-8"
    )
    (report_dir / "report.html").write_text(
        f"<html><code>{private_input}</code></html>", encoding="utf-8"
    )
    mem.export_experiment_snapshot("expA", report_dir / "memory_snapshot.sqlite")

    capsule = write_reproducible_capsule(
        report_dir, repo_root=tmp_path / "empty-repo"
    )
    verified = verify_reproducible_capsule(
        capsule, repo_root=tmp_path / "empty-repo"
    )

    assert verified["status"] == "warning"  # portable input needs relocation
    assert verified["memory_snapshot"] == {
        "present": True,
        "path": "report-expA/memory_snapshot.sqlite",
        "integrity_check": "ok",
        "experiment_id": "expA",
        "valid": True,
    }
    assert verified["inputs"]["relocation_required"] == [
        "input://sha256/deadbeef"
    ]

    with zipfile.ZipFile(capsule) as archive:
        manifest = json.loads(archive.read("capsule_manifest.json"))
        methodology_blob = archive.read("report-expA/methodology.json").decode()
        crate_blob = archive.read("report-expA/ro-crate-metadata.json").decode()
        names = archive.namelist()
        public_text = "\n".join(
            archive.read(name).decode("utf-8")
            for name in names
            if Path(name).suffix in {".json", ".html", ".md", ".txt"}
        )

    assert manifest["memory_snapshot"]["experiment_id"] == "expA"
    assert manifest["memory_snapshot"]["integrity_check"] == "ok"
    assert str(private_input) not in methodology_blob
    assert str(private_input) not in crate_blob
    assert str(private_input) not in public_text
    assert str(private_input.parent) not in public_text
    assert "memory_snapshot.sqlite" in crate_blob
    assert not any(name.startswith("/") for name in names)

    snap_name = manifest["memory_snapshot"]["path"]
    extracted = _snapshot_from_zip(capsule, snap_name, tmp_path / "extracted.sqlite")
    with sqlite3.connect(extracted) as conn:
        assert conn.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
        assert {row[0] for row in conn.execute("SELECT id FROM wings")} == {"expA"}


def test_invalid_memory_snapshot_never_replaces_existing_capsule(tmp_path):
    report_dir = tmp_path / "report"
    report_dir.mkdir()
    (report_dir / "methodology.json").write_text(
        json.dumps(_methodology(Path("/private/input.tsv"))), encoding="utf-8"
    )
    (report_dir / "report.html").write_text("<html></html>", encoding="utf-8")
    (report_dir / "memory_snapshot.sqlite").write_bytes(b"not sqlite")
    out = tmp_path / "report_capsule.zip"
    out.write_bytes(b"previous-good-capsule")

    with pytest.raises(ValueError, match="memory snapshot"):
        write_reproducible_capsule(
            report_dir, out, repo_root=tmp_path / "empty-repo"
        )

    assert out.read_bytes() == b"previous-good-capsule"
    assert not list(tmp_path.glob(".report_capsule.zip.*.tmp"))
