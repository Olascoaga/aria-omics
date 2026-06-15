"""U1+U5 pure-renderer tests (no Textual needed).

The cockpit's Rich renderers (`aria.ui.render`) are pure functions over an
ExperimentSnapshot, so they are unit-tested here in the standard env. The
Textual shell itself is covered (skip-if-missing) in test_cockpit_app.py.
"""

from __future__ import annotations

from rich.console import Console

from aria.runtime.experiment_view import (
    ExperimentSnapshot, FindingView, CheckpointView, LedgerNodeView,
    ModalityCardView, ProgressEvent,
)
from aria.ui import render


def _capture(renderable) -> str:
    console = Console(width=120, record=True)
    with console.capture() as cap:
        console.print(renderable)
    return cap.get()


def _snap(**over) -> ExperimentSnapshot:
    base = dict(
        experiment_id="exp123",
        phase="dispatch",
        progress=0.3,
        last_status=ProgressEvent(
            ts=None, sender="scrna_agent", text="Running clustering...",
            progress=0.3),
        findings_by_confidence={"HIGH": [], "MEDIUM": [], "LOW": [],
                                "INSUFFICIENT": []},
        pending_checkpoint=None,
        ledger=[],
        report_path=None,
        done=False,
        elapsed_s=125.0,
        silent_s=5.0,
    )
    base.update(over)
    return ExperimentSnapshot(**base)


def test_run_header_shows_identity_and_airgapped():
    out = _capture(render.render_run_header(
        _snap(), version="4.6.1", data_dir="/data/run",
        modalities=["scRNA"], organism="Homo sapiens", air_gapped=True))
    assert "exp123" in out
    assert "4.6.1" in out
    assert "scRNA" in out
    assert "Homo sapiens" in out
    assert "yes" in out          # air-gapped on
    assert "2m05s" in out        # elapsed 125s


def test_timeline_marks_current_phase():
    out = _capture(render.render_timeline(_snap(phase="dispatch", progress=0.3)))
    assert "dispatch" in out
    assert "30%" in out
    assert "Running clustering..." in out


def test_findings_counts_and_stream():
    snap = _snap(findings_by_confidence={
        "HIGH": [FindingView(ts=1, sender="a", confidence="HIGH",
                             summary="8 clusters resolved")],
        "MEDIUM": [],
        "LOW": [FindingView(ts=2, sender="a", confidence="LOW",
                            summary="weak batch signal")],
        "INSUFFICIENT": [],
    })
    out = _capture(render.render_findings(snap))
    assert "HIGH:1" in out
    assert "LOW:1" in out
    assert "8 clusters resolved" in out
    assert "weak batch signal" in out


def test_checkpoint_lists_numbered_options():
    cp = CheckpointView(message_id="m1", number=3,
                        title="Quality Control / Parameter Decision",
                        question="Accept inferred QC thresholds?",
                        options=["Accept", "Adjust", "Cancel"])
    out = _capture(render.render_checkpoint(_snap(pending_checkpoint=cp)))
    assert "CP3" in out
    assert "Accept inferred QC thresholds?" in out
    assert "[1]" in out and "Accept" in out
    assert "[3]" in out and "Cancel" in out


def test_checkpoint_idle_and_complete():
    assert "No pending checkpoint" in _capture(
        render.render_checkpoint(_snap()))
    assert "Run complete" in _capture(
        render.render_checkpoint(_snap(done=True)))


def test_ledger_highlights_divergence():
    ledger = [
        LedgerNodeView(modality="scRNA", analysis="qc", label="Quality control",
                       status="ran", planned=True, divergence=False,
                       reason=None, node_id="ledger://scRNA/qc"),
        LedgerNodeView(modality="scRNA", analysis="pseudobulk_de",
                       label="Pseudobulk differential expression",
                       status="not_run", planned=True, divergence=True,
                       reason="planned but not run",
                       node_id="ledger://scRNA/pseudobulk_de"),
    ]
    out = _capture(render.render_ledger(_snap(ledger=ledger)))
    assert "Quality control" in out
    assert "Pseudobulk differential expression" in out
    assert "not_run" in out
    assert "1 divergence" in out


def test_ledger_empty_state():
    out = _capture(render.render_ledger(_snap(ledger=[])))
    assert "no ledger yet" in out
    assert "0 divergence" in out


def test_readiness_cards_render():
    cards = [
        ModalityCardView(modality="scRNA", validation_level="validated",
                         status="green", dispatch_policy="allowed", reason=""),
        ModalityCardView(modality="scATAC", validation_level="alpha",
                         status="yellow", dispatch_policy="requires_ack",
                         reason="alpha lane",
                         findings=["scATAC requires explicit acknowledgement."]),
        ModalityCardView(modality="bulk_ATAC", validation_level="scaffold",
                         status="red", dispatch_policy="blocked",
                         reason="not validated for dispatch"),
    ]
    out = _capture(render.render_readiness(_snap(readiness=cards)))
    assert "scRNA" in out and "allowed" in out
    assert "scATAC" in out and "requires_ack" in out
    assert "alpha" in out
    assert "bulk_ATAC" in out and "blocked" in out
    assert "requires explicit acknowledgement" in out


def test_readiness_empty_state():
    out = _capture(render.render_readiness(_snap(readiness=[])))
    assert "No readiness cards yet" in out
