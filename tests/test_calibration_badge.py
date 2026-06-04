"""W-CALIB: the numerical-calibration report badge.

The badge is a pure, deterministic renderer over the calibration manifest
produced by ``aria.benchmarks.run_calibration_suite``. It must:

1. honestly say "not measured" with NO metric when no manifest is attached (the
   normal report path — calibration is a build property, not run per report);
2. render the measured recall / empirical-FDR / null false-positive-rate numbers
   and a pass/fail badge when a real manifest IS attached.

The helper itself is pure, but importing it pulls the `aria.agents` package
(litellm at import time), so the test gates on litellm like the other
report-render tests; it runs in the light / PR lane where litellm is present.
"""

import pytest

pytest.importorskip("litellm")

from aria.agents.narrative.report_sections import _build_calibration_badge


def test_badge_is_honest_when_no_manifest():
    html = _build_calibration_badge(None)
    assert "not measured in this run" in html
    # Must NOT assert any calibration metric / pass status when nothing was run.
    assert "PASS" not in html
    assert "W-CALIB" in html


def test_badge_handles_non_measured_manifest():
    # A dict without measured=True is treated as not measured (no fabrication).
    html = _build_calibration_badge({"status": "pass"})
    assert "not measured in this run" in html
    assert "PASS" not in html


def test_badge_renders_measured_metrics_and_pass_status():
    manifest = {
        "status": "pass",
        "measured": True,
        "seed": 11,
        "summary": {
            "bulk_recall": 0.95,
            "bulk_empirical_fdr": 0.02,
            "bulk_null_fpr": 0.001,
            "pseudobulk_recall": 0.93,
            "pseudobulk_empirical_fdr": 0.03,
            "pseudobulk_null_fpr": 0.0,
        },
    }
    html = _build_calibration_badge(manifest)
    assert "PASS" in html
    assert "badge high" in html
    assert "0.950" in html          # bulk recall rendered
    assert "0.001" in html          # bulk null FPR rendered
    assert "Seed" in html and "11" in html


def test_badge_marks_failure():
    manifest = {
        "status": "fail",
        "measured": True,
        "seed": 11,
        "summary": {
            "bulk_recall": 0.4,
            "bulk_empirical_fdr": 0.3,
            "bulk_null_fpr": 0.2,
            "pseudobulk_recall": 0.5,
            "pseudobulk_empirical_fdr": 0.1,
            "pseudobulk_null_fpr": 0.15,
        },
    }
    html = _build_calibration_badge(manifest)
    assert "FAIL" in html
    assert "badge low" in html
