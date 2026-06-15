"""U1+U5 pure-renderer tests (no Textual needed).

The cockpit's Rich renderers (`aria.ui.render`) are pure functions over an
ExperimentSnapshot, so they are unit-tested here in the standard env. The
Textual shell itself is covered (skip-if-missing) in test_cockpit_app.py.
"""

from __future__ import annotations

from rich.console import Console

from aria.runtime.experiment_view import (
    ArtifactView, ExperimentHistoryView, ExperimentSnapshot, FindingView,
    CheckpointView, LedgerNodeView, ModalityCardView, ProgressEvent,
    ResourceView,
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


def test_resources_render_local_state():
    resources = [
        ResourceView(category="env", name="aria-rna-env", status="ready",
                     detail="SetupAgent completed environment checks."),
        ResourceView(category="geneset", name="Local GMT libraries",
                     status="missing", detail="No local GMT libraries staged.",
                     path="/tmp/genesets",
                     action="Run scripts/fetch_genesets.py explicitly."),
        ResourceView(category="privacy", name="Network egress",
                     status="blocked",
                     detail="Air-gapped mode is active; network fetches are blocked.",
                     action="Use local staged resources/caches."),
    ]
    out = _capture(render.render_resources(_snap(resources=resources)))
    assert "aria-rna-env" in out
    assert "Local GMT libraries" in out
    assert "missing" in out
    assert "Network egress" in out
    assert "blocked" in out
    assert "/tmp/genesets" in out


def test_resources_empty_state():
    out = _capture(render.render_resources(_snap(resources=[])))
    assert "No resource snapshot yet" in out


def test_artifacts_render_bundle_and_claims():
    artifacts = [
        ArtifactView(category="report", name="report.html", status="present",
                     detail="Narrative HTML report.", path="/tmp/r/report.html"),
        ArtifactView(category="methodology", name="methodology.json",
                     status="present", detail="Provenance, claims, run-ledger.",
                     path="/tmp/r/methodology.json"),
        ArtifactView(category="figure", name="volcano.png", status="present",
                     detail="12,345 bytes", path="/tmp/r/figures/volcano.png"),
        ArtifactView(category="claim", name="bulk.de.executive_summary",
                     status="ok", detail="tier=descriptive · HIGH · ledger=ran"),
        ArtifactView(category="claim", name="bulk.de.pathway",
                     status="violation",
                     detail="tier=associative · LOW · ledger=not_run"),
    ]
    out = _capture(render.render_artifacts(_snap(artifacts=artifacts)))
    assert "report.html" in out
    assert "methodology.json" in out
    assert "volcano.png" in out
    assert "violation" in out
    assert "claim violation(s)" in out


def test_artifacts_empty_state():
    out = _capture(render.render_artifacts(_snap(artifacts=[])))
    assert "No artifacts yet" in out


def test_history_renders_prior_runs_and_resume_points():
    history = [
        ExperimentHistoryView(
            experiment_id="abcd1234efgh", name="H9 RNA timecourse",
            organism="Homo sapiens", genome="hg38",
            updated_at="2026-06-10T09:30:00", summary="DEGs found.",
            n_decisions=4, last_decision="CP2.6: confirm",
            modalities=["bulk_RNA"], report_path="/r/report.html",
            has_report=True),
        ExperimentHistoryView(
            experiment_id=" zz99", name="ATAC pilot",
            organism="Mus musculus", genome="mm10",
            updated_at="2026-06-09T12:00:00", summary="",
            n_decisions=1, last_decision=None,
            modalities=["scATAC"], report_path=None, has_report=False),
    ]
    out = _capture(render.render_history(history))
    assert "H9 RNA timecourse" in out
    assert "bulk_RNA" in out
    assert "report" in out
    assert "no report" in out
    assert "1 with report" in out


def test_history_empty_state():
    out = _capture(render.render_history([]))
    assert "No prior experiments yet" in out


def test_findings_disclose_truncated_older():
    many = [
        FindingView(ts=i, sender="a", confidence="MEDIUM",
                    summary=f"finding {i}")
        for i in range(20)
    ]
    snap = _snap(findings_by_confidence={
        "HIGH": [], "MEDIUM": many, "LOW": [], "INSUFFICIENT": [],
    })
    out = _capture(render.render_findings(snap, limit=12))
    # Older findings are disclosed, not silently dropped.
    assert "+8 older" in out
    assert "20 total" in out
    assert "finding 19" in out      # most recent shown
    assert "finding 0" not in out   # oldest is off-screen


def test_checkpoint_complete_shows_report_path():
    out = _capture(render.render_checkpoint(
        _snap(done=True, report_path="/home/x/.aria/reports/run/report.html")))
    assert "Run complete" in out
    assert "/home/x/.aria/reports/run/report.html" in out


def test_status_banner_shows_message_and_version():
    out = _capture(render.render_status_banner("Starting analysis…", "4.6.1"))
    assert "Starting analysis" in out
    assert "4.6.1" in out
    err = _capture(render.render_status_banner("Could not start", "4.6.1",
                                               error=True))
    assert "Could not start" in err


def test_mode_bar_highlights_active_view():
    out = _capture(render.render_mode_bar("ledger"))
    for mode in ("findings", "ledger", "readiness", "resources", "artifacts"):
        assert mode in out
    # Toggle keys are surfaced for discoverability.
    assert "[l] ledger" in out
    assert "[a] artifacts" in out
