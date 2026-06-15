"""U0 read-model (ExperimentView) tests.

The interactive TUI (`aria.tui`) and the headless runner (`aria.headless`) each
independently polled the global MessageBus and re-derived the same run state. A
Textual cockpit would have been a third copy. `aria.runtime.experiment_view`
turns the bus state (plus an optional ExperimentSession for the run-ledger) into
one immutable, UI-agnostic snapshot. These tests pin that contract:

- status-text normalization across the ``message`` vs ``status`` payload keys
  (agents publish ``message``; the legacy consumers read ``status``);
- phase transitions audit -> design -> dispatch -> done;
- confidence grouping of findings;
- pending-checkpoint surfacing (resolved excluded);
- report-path detection from the bus (text + payload key);
- run-ledger mapping from an ExperimentSession.
"""

from __future__ import annotations

import sys
import types
import uuid
from datetime import datetime, timedelta

# Keep imports cheap/offline like the other agent tests.
litellm_stub = types.ModuleType("litellm")
litellm_stub.completion = lambda *args, **kwargs: None
sys.modules.setdefault("litellm", litellm_stub)

from aria.bus.message_bus import bus, Message, MessageType, Confidence
from aria.runtime.experiment_session import ExperimentSession
from aria.runtime.experiment_view import (
    build_snapshot,
    status_text,
    ExperimentSnapshot,
    CheckpointView,
    CHECKPOINT_TITLES,
)


def _eid() -> str:
    return "view-" + uuid.uuid4().hex[:8]


def _status(eid, sender, text, progress, *, ts=None):
    bus.publish(Message(
        sender=sender, receiver="orchestrator", type=MessageType.STATUS,
        payload={"message": text, "progress": progress},
        experiment_id=eid, timestamp=ts or datetime.now(),
    ))


def _finding(eid, sender, summary, confidence):
    bus.publish(Message(
        sender=sender, receiver="all", type=MessageType.FINDING,
        confidence=confidence, payload={"summary": summary},
        experiment_id=eid,
    ))


def _escalation(eid, checkpoint, question, options):
    m = Message(
        sender="orchestrator", receiver="all", type=MessageType.ESCALATION,
        payload={"checkpoint": checkpoint, "question": question,
                 "options": options},
        checkpoint=checkpoint, experiment_id=eid,
    )
    bus.publish(m)
    return m


# ── status-text normalization ────────────────────────────────────────────────

def test_status_text_reads_either_key():
    assert status_text({"message": "hi"}) == "hi"
    assert status_text({"status": "yo"}) == "yo"
    # When both exist, the explicit ``status`` wins (legacy consumer key).
    assert status_text({"status": "yo", "message": "hi"}) == "yo"
    assert status_text({}) == ""
    assert status_text(None) == ""


# ── phase transitions ────────────────────────────────────────────────────────

def test_audit_phase_before_dispatch():
    eid = _eid()
    _status(eid, "orchestrator", "ARIA starting analysis...", 0.0)
    _escalation(eid, 1, "Confirm modality", ["Continue", "Correct"])
    snap = build_snapshot(eid)
    assert isinstance(snap, ExperimentSnapshot)
    assert snap.phase == "audit"
    assert snap.done is False
    assert isinstance(snap.pending_checkpoint, CheckpointView)
    assert snap.pending_checkpoint.number == 1
    assert snap.pending_checkpoint.title == CHECKPOINT_TITLES[1]


def test_design_phase_on_design_checkpoint():
    eid = _eid()
    _status(eid, "orchestrator", "Dispatching agents...", 0.1)
    _escalation(eid, 2.1, "Confirm experimental groups", ["Continue"])
    snap = build_snapshot(eid)
    assert snap.phase == "design"
    assert snap.pending_checkpoint.number == 2.1


def test_dispatch_phase_and_findings_grouped():
    eid = _eid()
    _status(eid, "scrna_agent", "Running clustering...", 0.3)
    _finding(eid, "scrna_agent", "8 clusters resolved", Confidence.HIGH)
    _finding(eid, "scrna_agent", "weak batch signal", Confidence.LOW)
    _finding(eid, "scrna_agent", "ambiguous label", Confidence.INSUFFICIENT)
    snap = build_snapshot(eid)
    assert snap.phase == "dispatch"
    assert abs(snap.progress - 0.3) < 1e-9
    assert snap.last_status is not None
    assert snap.last_status.text == "Running clustering..."
    assert len(snap.findings_by_confidence["HIGH"]) == 1
    assert len(snap.findings_by_confidence["LOW"]) == 1
    assert len(snap.findings_by_confidence["INSUFFICIENT"]) == 1
    assert snap.findings_by_confidence["HIGH"][0].summary == "8 clusters resolved"


def test_done_and_report_path_from_status_text():
    eid = _eid()
    _status(eid, "scrna_agent", "Running clustering...", 0.3)
    _status(eid, "narrative_agent",
            "Report saved: /tmp/aria/run/report.html", 1.0)
    snap = build_snapshot(eid)
    assert snap.done is True
    assert snap.phase == "done"
    assert snap.report_path == "/tmp/aria/run/report.html"


def test_report_path_from_payload_key():
    eid = _eid()
    bus.publish(Message(
        sender="narrative_agent", receiver="orchestrator",
        type=MessageType.STATUS,
        payload={"message": "done", "progress": 1.0,
                 "report_path": "/p/report.html"},
        experiment_id=eid,
    ))
    snap = build_snapshot(eid)
    assert snap.done is True
    assert snap.report_path == "/p/report.html"


def test_resolved_checkpoint_not_pending():
    eid = _eid()
    m = _escalation(eid, 3, "QC params", ["Accept", "Adjust"])
    bus.resolve_checkpoint(m.id, {"decision": "Accept"})
    snap = build_snapshot(eid)
    assert snap.pending_checkpoint is None


def test_elapsed_and_silent_seconds():
    eid = _eid()
    t0 = datetime.now()
    _status(eid, "scrna_agent", "working", 0.2, ts=t0 - timedelta(seconds=30))
    now = t0
    snap = build_snapshot(eid, start_time=t0 - timedelta(seconds=120), now=now)
    assert 110 <= snap.elapsed_s <= 130
    assert 25 <= snap.silent_s <= 35


# ── run-ledger mapping from a session ────────────────────────────────────────

def test_ledger_maps_planned_but_not_run():
    eid = _eid()
    _status(eid, "scrna_agent", "working", 0.4)
    session = ExperimentSession(experiment_id=eid)
    session.exp_context = {"design_intelligence": {
        "recommended": ["Donor-level pseudobulk DESeq2 between conditions."],
        "optional": [],
    }}
    session.agent_results = {"scrna_agent": {"findings": {
        "qc": {"status": "success", "n_cells": 100},
        # pseudobulk_de absent -> planned-but-not-run divergence
    }}}
    snap = build_snapshot(eid, session=session)
    nodes = {n.analysis: n for n in snap.ledger if n.modality == "scRNA"}
    assert "pseudobulk_de" in nodes
    assert nodes["pseudobulk_de"].planned is True
    assert nodes["pseudobulk_de"].status == "not_run"
    assert nodes["pseudobulk_de"].divergence is True
    assert nodes["qc"].status == "ran"


def test_ledger_empty_without_session():
    eid = _eid()
    _status(eid, "scrna_agent", "working", 0.4)
    snap = build_snapshot(eid)
    assert snap.ledger == []


def test_readiness_maps_from_session():
    eid = _eid()
    _status(eid, "scrna_agent", "working", 0.4)
    session = ExperimentSession(experiment_id=eid)
    session.exp_context = {"readiness_cards": {
        "scRNA": {"modality": "scRNA", "validation_level": "validated",
                  "status": "green", "dispatch_policy": "allowed",
                  "reason": "", "findings": []},
        "scATAC": {"modality": "scATAC", "validation_level": "alpha",
                   "status": "yellow", "dispatch_policy": "requires_ack",
                   "reason": "alpha lane",
                   "findings": [{"message": "requires ack", "severity": "warning"}]},
    }}
    snap = build_snapshot(eid, session=session)
    cards = {c.modality: c for c in snap.readiness}
    assert cards["scATAC"].dispatch_policy == "requires_ack"
    assert cards["scATAC"].validation_level == "alpha"
    assert cards["scATAC"].findings == ["requires ack"]
    assert cards["scRNA"].status == "green"


def test_readiness_empty_without_session():
    eid = _eid()
    _status(eid, "scrna_agent", "working", 0.4)
    assert build_snapshot(eid).readiness == []


def test_resources_map_local_stores_and_airgap(tmp_path, monkeypatch):
    eid = _eid()
    _status(eid, "scrna_agent", "working", 0.4)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))

    gmt_base = tmp_path / "genesets"
    gmt_lib = gmt_base / "GO_TEST"
    gmt_lib.mkdir(parents=True)
    (gmt_lib / "GO_TEST.gmt").write_text(
        "term\tdesc\tGENE1\tGENE2\n", encoding="utf-8")
    monkeypatch.setenv("ARIA_GMT_DIR", str(gmt_base))

    motif_base = tmp_path / "motifs"
    motif_lib = motif_base / "JASPAR_TEST"
    motif_lib.mkdir(parents=True)
    (motif_lib / "JASPAR_TEST.meme").write_text(
        "MEME version 4\n\nMOTIF X\n", encoding="utf-8")
    monkeypatch.setenv("ARIA_MOTIF_DIR", str(motif_base))

    genome_base = tmp_path / "genomes"
    (genome_base / "hg38").mkdir(parents=True)
    (genome_base / "hg38" / "hg38.fa").write_text(">chr1\nACGT\n",
                                                    encoding="utf-8")
    monkeypatch.setenv("ARIA_GENOME_DIR", str(genome_base))
    monkeypatch.delenv("ARIA_GENOME_FASTA", raising=False)
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")

    session = ExperimentSession(experiment_id=eid)
    session.exp_context = {
        "modalities": {"scRNA": ["sample.h5ad"], "scATAC": ["sample.h5mu"]},
        "organism": "Homo sapiens",
        "genome": "hg38",
        "motif_collection": "JASPAR_TEST",
        "celltypist_tissue_hint": "brain",
        "air_gapped": True,
    }
    session.agent_results = {"setup_agent": {"status": "done"}}

    snap = build_snapshot(eid, session=session)
    resources = {(r.category, r.name): r for r in snap.resources}
    assert resources[("env", "aria-rna-env")].status == "ready"
    assert resources[("env", "aria-chromatin-env")].status == "ready"

    by_category = {r.category: r for r in snap.resources}
    assert by_category["geneset"].status == "ready"
    assert by_category["motif"].status == "ready"
    assert by_category["genome"].status == "ready"
    assert by_category["privacy"].status == "blocked"

    celltypist = [r for r in snap.resources if r.category == "celltypist"]
    assert celltypist
    assert celltypist[0].name == "CellTypist: Adult_Human_PrefrontalCortex.pkl"
    assert celltypist[0].status == "missing"


def test_artifacts_read_report_bundle_from_disk(tmp_path):
    import json
    from aria.runtime.experiment_view import _artifact_views

    report_dir = tmp_path / "report_x"
    (report_dir / "figures").mkdir(parents=True)
    (report_dir / "tables").mkdir(parents=True)
    report_file = report_dir / "report.html"
    report_file.write_text("<html></html>", encoding="utf-8")
    (report_dir / "figures" / "volcano.png").write_bytes(b"\x89PNG")
    (report_dir / "tables" / "de.tsv").write_text("gene\tlfc\n", encoding="utf-8")
    methodology = {
        "claims": [
            {"claim_id": "bulk.de.exec", "text": "DEGs detected.",
             "tier": "descriptive", "confidence": "high",
             "ledger_status": "ran"},
            {"claim_id": "bulk.de.path", "text": "Pathway X is causal.",
             "tier": "associative", "confidence": "low",
             "ledger_status": "not_run"},
        ],
        "run_ledger": {
            "claim_linkage": {
                "violations": [{"claim_id": "bulk.de.path",
                                "ledger_status": "not_run"}],
            }
        },
    }
    (report_dir / "methodology.json").write_text(
        json.dumps(methodology), encoding="utf-8")

    arts = {(a.category, a.name): a for a in _artifact_views(str(report_file))}
    assert arts[("report", "report.html")].status == "present"
    assert arts[("methodology", "methodology.json")].status == "present"
    assert arts[("figure", "volcano.png")].status == "present"
    assert arts[("table", "de.tsv")].status == "present"
    assert arts[("claim", "bulk.de.exec")].status == "ok"
    # The not-run associative claim is flagged via the existing claim_linkage.
    assert arts[("claim", "bulk.de.path")].status == "violation"


def test_artifacts_empty_until_report_path_known():
    from aria.runtime.experiment_view import _artifact_views
    assert _artifact_views(None) == []
    assert _artifact_views("/nonexistent/dir/report.html") == []
