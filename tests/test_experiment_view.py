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
