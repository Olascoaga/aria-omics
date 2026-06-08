"""v4.6 scATAC step 7 — focused chromatin report/ledger smoke.

Integration-level checks above the individual scripts/narrator:

- the Chromatin section actually renders in the report findings HTML (exercises
  the report_builder section routing fixed in step 6);
- the v4.6 matrix-flow finding keys reconcile in the run ledger as "ran"
  (qc / dimensionality_reduction / differential_accessibility / motif_enrichment);
- a thin chromatin report (only QC ran) is surfaced as a divergence;
- the agent's stored finding keys stay aligned with the ledger's finding_keys
  (a guard against future agent/ledger drift).
"""

import pytest

from tests.test_chromatin_narrator_v46 import _findings
from aria.agents.narrative.narrators.chromatin import ChromatinNarrator


def _chromatin_blocks(findings):
    return ChromatinNarrator().collect(
        "chromatin_agent", {"status": "done", "findings": findings}, {})


# ── Report findings section renders the Chromatin section ──────────────────────

def test_chromatin_section_renders_in_findings_html():
    from aria.agents.narrative_agent import NarrativeAgent
    agent = NarrativeAgent.__new__(NarrativeAgent)

    findings = _findings()
    agent_results = {"chromatin_agent": {"status": "done", "findings": findings}}
    blocks = _chromatin_blocks(findings)

    html = agent._build_findings_section(
        {}, agent_results=agent_results, narrative_blocks=blocks,
        report_dir=None)

    # The section label and the real cluster / DA / motif content are present.
    assert "Chromatin" in html
    assert "8 clusters" in html
    assert "differentially accessible" in html
    assert "JASPAR2024_CORE_vertebrates" in html
    # honest no-result path is visible, not silently dropped
    assert "not run" in html.lower()


def test_chromatin_section_absent_without_chromatin_results():
    # Regression: the chromatin section must not appear for an RNA-only run.
    from aria.agents.narrative_agent import NarrativeAgent
    agent = NarrativeAgent.__new__(NarrativeAgent)
    html = agent._build_findings_section({}, agent_results={},
                                         narrative_blocks=[], report_dir=None)
    assert "Chromatin" not in html


# ── Run-ledger reconciliation for the v4.6 matrix flow ────────────────────────

def test_matrix_findings_reconcile_as_ran():
    from aria.agents.narrative.run_ledger import build_run_ledger
    exp_ctx = {"design_intelligence": {
        "recommended": ["LSI dimensionality reduction over the peak matrix.",
                        "Differential accessibility per cluster.",
                        "TF motif enrichment over accessible peaks."],
        "optional": []}}
    agent_results = {"chromatin_agent": {"status": "done",
                                         "findings": _findings()}}
    ledger = build_run_ledger(exp_ctx, agent_results)
    by_key = {e["analysis"]: e for e in ledger["entries"]}
    for analysis in ("qc", "dimensionality_reduction",
                     "differential_accessibility", "motif_enrichment"):
        assert by_key[analysis]["status"] == "ran", analysis
        assert by_key[analysis]["divergence"] is False, analysis


def test_thin_chromatin_report_is_a_divergence():
    from aria.agents.narrative.run_ledger import build_run_ledger
    exp_ctx = {"design_intelligence": {
        "recommended": ["LSI dimensionality reduction over the peak matrix.",
                        "Differential accessibility per cluster."],
        "optional": []}}
    # Only QC ran (the v4.6 stack stalled after QC).
    agent_results = {"chromatin_agent": {"status": "done", "findings": {
        "qc": {"status": "success", "n_cells": 3143, "n_peaks": 60990}}}}
    ledger = build_run_ledger(exp_ctx, agent_results)
    by_key = {e["analysis"]: e for e in ledger["entries"]}
    assert by_key["qc"]["status"] == "ran"
    assert by_key["dimensionality_reduction"]["divergence"] is True
    assert by_key["differential_accessibility"]["divergence"] is True
    assert ledger["n_divergences"] >= 2


# ── Agent finding keys stay aligned with the ledger finding_keys ──────────────

def test_agent_finding_keys_match_ledger_finding_keys():
    """The matrix flow stores findings under qc/lsi/differential_accessibility/
    motifs; the ledger must recognise those exact keys or a real run would show
    spurious divergences."""
    from aria.agents.narrative.run_ledger import _CHROMATIN_ANALYSES
    spec = {s["key"]: s["finding_keys"] for s in _CHROMATIN_ANALYSES}
    assert "lsi" in spec["dimensionality_reduction"]
    assert "differential_accessibility" in spec["differential_accessibility"]
    assert "motifs" in spec["motif_enrichment"]
    assert "qc" in spec["qc"]
