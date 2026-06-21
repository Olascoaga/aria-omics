"""S8 (pre-integration audit): scientific calibration gate.

The recovery benchmarks compute recall + empirical FDR on synthetic ground truth.
Their bounds used to be hardcoded per test; this fence makes them a single
versioned baseline (calibration_baseline.json) checked by check_calibration, so a
refactor that erodes recall or inflates empirical FDR fails the build.

Light arm runs anywhere (parses the baseline + exercises the gate logic). The heavy
arm runs the real bulk ATAC DA recovery through the gate and needs pydeseq2 (heavy
CI lane / aria-rna-env).
"""

from __future__ import annotations

import pytest

from aria.benchmarks.calibration import (
    baseline_bounds,
    check_calibration,
    load_baseline,
)

# The recovery benchmarks that MUST have a versioned baseline (no silent gaps).
_REQUIRED_BENCHMARKS = {
    "synthetic_de_pseudobulk",
    "scatac_pseudobulk_da",
    "bulk_atac_da",
}


# ── light arm: baseline structure + gate logic (any env) ──────────────────

def test_baseline_covers_required_benchmarks():
    bm = load_baseline()["benchmarks"]
    missing = _REQUIRED_BENCHMARKS - set(bm)
    assert not missing, f"calibration baseline missing benchmarks: {missing}"
    for name, entry in bm.items():
        assert "min_recall" in entry and "max_empirical_fdr" in entry, name
        assert 0.0 <= entry["min_recall"] <= 1.0, name
        assert 0.0 <= entry["max_empirical_fdr"] <= 1.0, name
        assert entry.get("provenance"), f"{name} baseline lacks provenance"


def test_gate_passes_at_or_above_baseline():
    m = check_calibration(
        "bulk_atac_da",
        {"recall": 0.667, "empirical_fdr": 0.0, "n_false_positive": 0, "n_called": 200},
    )
    assert m["status"] == "pass", m


def test_gate_fails_on_recall_regression():
    m = check_calibration(
        "bulk_atac_da",
        {"recall": 0.30, "empirical_fdr": 0.0, "n_false_positive": 0, "n_called": 200},
    )
    assert m["status"] == "fail"
    assert any("recall" in v for v in m["violations"])


def test_gate_fails_on_fdr_inflation():
    m = check_calibration(
        "synthetic_de_pseudobulk", {"recall": 0.95, "empirical_fdr": 0.4})
    assert m["status"] == "fail"
    assert any("empirical_fdr" in v for v in m["violations"])


def test_gate_fails_on_missing_metric():
    # A gate must not pass on absent evidence.
    m = check_calibration("scatac_pseudobulk_da", {"empirical_fdr": 0.0})
    assert m["status"] == "fail"
    assert any("missing metric 'recall'" in v for v in m["violations"])


def test_unknown_benchmark_raises():
    with pytest.raises(KeyError):
        baseline_bounds("does_not_exist")


# ── heavy arm: real bulk ATAC DA recovery through the gate ─────────────────

def test_bulk_atac_da_meets_calibration_baseline(tmp_path):
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.synthetic_atac_da import run_bulk_atac_da_benchmark

    bounds = baseline_bounds("bulk_atac_da")
    m = run_bulk_atac_da_benchmark(
        seed=11,
        work_dir=str(tmp_path),
        min_recall=bounds["min_recall"],
        max_empirical_fdr=bounds["max_empirical_fdr"],
    )
    gate = check_calibration("bulk_atac_da", m)
    assert gate["status"] == "pass", {"gate": gate, "metrics": m}
