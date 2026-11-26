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
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.synthetic_atac_da import run_bulk_atac_da_benchmark
    from aria.scripts.chromatin_bulk_diffacc import chromatin_bulk_diffacc

    bounds = baseline_bounds("bulk_atac_da")
    m = run_bulk_atac_da_benchmark(
        seed=11,
        work_dir=str(tmp_path),
        min_recall=bounds["min_recall"],
        max_empirical_fdr=bounds["max_empirical_fdr"],
    )
    gate = check_calibration("bulk_atac_da", m)
    assert gate["status"] == "pass", {"gate": gate, "metrics": m}

    # Signed LFC direction oracle (kills log2FoldChange negation and (test, ref) swap).
    rng = np.random.default_rng(11)
    samples = ["ctrl_1", "ctrl_2", "ctrl_3", "ctrl_4", "trt_1", "trt_2", "trt_3", "trt_4"]
    n_peaks = 120
    peaks = [f"chr1:{i * 1000 + 100}-{i * 1000 + 600}" for i in range(n_peaks)]
    base = rng.uniform(60.0, 180.0, size=n_peaks)
    mat = np.zeros((n_peaks, len(samples)), dtype=int)
    for j, s in enumerate(samples):
        mu = base.copy()
        if s.startswith("trt"):
            mu[:15] *= 4.0     # planted up in trt vs ctrl
            mu[15:30] *= 0.25  # planted down in trt vs ctrl
        mat[:, j] = rng.negative_binomial(n=25, p=25.0 / (25.0 + mu))
    counts_p = tmp_path / "dir_counts.tsv"
    pd.DataFrame(mat, index=peaks, columns=samples).to_csv(
        counts_p, sep="\t", index_label="peak_id"
    )
    meta_p = tmp_path / "dir_meta.tsv"
    pd.DataFrame({
        "sample_id": samples,
        "condition": ["ctrl"] * 4 + ["trt"] * 4,
        "replicate": ["r1", "r2", "r3", "r4"] * 2,
    }).to_csv(meta_p, sep="\t", index=False)

    res = chromatin_bulk_diffacc({
        "data_type": "bulk_ATAC",
        "counts_matrix_path": str(counts_p),
        "sample_metadata_path": str(meta_p),
        "comparisons": [["trt", "ctrl"]],
        "padj_threshold": 0.05,
        "lfc_threshold": 0.5,
        "top_n": 100,
        "n_cpus": 1,
        "output_dir": str(tmp_path / "dir_out"),
    })
    assert res["status"] == "success" and res["ran"] is True, res
    comp = res["comparisons"][0]
    assert comp["status"] == "success", comp
    by_peak = {p["peak"]: p for p in comp["top_peaks"]}
    up_hits = [by_peak[pk] for pk in peaks[:15] if pk in by_peak]
    dn_hits = [by_peak[pk] for pk in peaks[15:30] if pk in by_peak]
    assert len(up_hits) >= 12 and all(p["log2fc"] > 0.5 for p in up_hits)
    assert len(dn_hits) >= 12 and all(p["log2fc"] < -0.5 for p in dn_hits)
