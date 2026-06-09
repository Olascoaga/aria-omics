"""Guards for the A1 SEQC/MAQC reference-data lane (aria/benchmarks/reference_seqc.py).

Light unit tests (numpy/pandas) run everywhere; the scorer end-to-end test is
pydeseq2-gated (heavy lane); the real-bundle test runs only when
``ARIA_SEQC_MAQC_BUNDLE`` points at a staged reference bundle.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest


def test_auc_separates_classes_and_handles_ties():
    from aria.benchmarks.reference_seqc import _auc

    assert _auc([3.0, 4.0, 5.0], [0.0, 1.0, 2.0]) == 1.0   # perfect
    assert _auc([0.0, 1.0, 2.0], [3.0, 4.0, 5.0]) == 0.0   # reversed
    assert _auc([1.0, 1.0], [1.0, 1.0]) == 0.5             # all ties
    assert _auc([1.0], []) is None                          # empty class


def test_load_bundle_returns_none_when_absent(tmp_path):
    from aria.benchmarks.reference_seqc import load_seqc_maqc_bundle

    assert load_seqc_maqc_bundle(tmp_path / "missing") is None


def _write_bundle(bundle_dir: Path, counts, samples, taqman):
    import pandas as pd

    bundle_dir.mkdir(parents=True, exist_ok=True)
    c = counts.copy()
    c.insert(0, "gene", c.index)
    c.to_csv(bundle_dir / "counts.tsv", sep="\t", index=False)
    pd.DataFrame(samples).to_csv(bundle_dir / "samples.tsv", sep="\t", index=False)
    pd.DataFrame(taqman).to_csv(bundle_dir / "taqman.tsv", sep="\t", index=False)


def test_load_bundle_parses_staged_files(tmp_path):
    import pandas as pd
    from aria.benchmarks.reference_seqc import load_seqc_maqc_bundle

    counts = pd.DataFrame(
        {"A_1": [10, 20], "B_1": [5, 40]}, index=["GENE_1", "GENE_2"]
    )
    _write_bundle(
        tmp_path,
        counts,
        {"sample": ["A_1", "B_1"], "group": ["A", "B"]},
        {"gene": ["GENE_1", "GENE_2"], "log2_ab": [1.0, -3.0]},
    )
    bundle = load_seqc_maqc_bundle(tmp_path)
    assert bundle is not None
    assert list(bundle["counts"].index) == ["GENE_1", "GENE_2"]
    assert bundle["taqman_log2_ab"]["GENE_2"] == -3.0
    assert set(bundle["samples"]["group"]) == {"A", "B"}


def test_runner_skips_honestly_without_bundle(tmp_path):
    from aria.benchmarks.reference_seqc import run_seqc_maqc_a1_benchmark

    out = run_seqc_maqc_a1_benchmark(tmp_path / "absent")
    assert out["status"] == "skipped"
    assert "fabricat" not in str(out).lower() or "Nothing is" in str(out)
    assert out["expected_bundle"]["files"] == ["counts.tsv", "samples.tsv", "taqman.tsv"]


def test_score_seqc_maqc_recovers_reference_truth(tmp_path):
    """End-to-end: build a bundle from the simulated bulk dataset (treating its
    known true log2FC as the TaqMan truth) and confirm ARIA's real DE path scores
    high LFC concordance and AUC against it. Maps COND_B->A, COND_A->B so ARIA's
    log2(A/B) aligns with true_log2fc = log2(COND_B/COND_A)."""
    pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.reference_seqc import (
        load_seqc_maqc_bundle,
        score_seqc_maqc_a1,
    )
    from aria.benchmarks.synthetic_de import simulate_bulk_dataset

    ds = simulate_bulk_dataset(n_genes=400, n_de=60, replicates_per_condition=5, seed=7)
    cond = ds.metadata["condition"]
    group_of = {"COND_B": "A", "COND_A": "B"}
    samples = {
        "sample": list(ds.counts.columns),
        "group": [group_of[cond[s]] for s in ds.counts.columns],
    }
    taqman = {
        "gene": list(ds.true_log2fc),
        "log2_ab": [ds.true_log2fc[g] for g in ds.true_log2fc],
    }
    _write_bundle(tmp_path, ds.counts, samples, taqman)

    bundle = load_seqc_maqc_bundle(tmp_path)
    manifest = score_seqc_maqc_a1(bundle, numerator="A", denominator="B")

    conc = manifest["axes"]["lfc_concordance"]
    det = manifest["axes"]["taqman_de_detection"]
    assert manifest["status"] == "pass", manifest
    # SEQC-standard log-ratio concordance (Pearson) is the gated metric.
    assert conc["pearson"] >= 0.7
    # Rank concordance on real-signal genes is robust to the null-tie mass.
    assert conc["spearman_signal"] >= 0.5
    # TaqMan-DE detection AUC.
    assert det["auc"] >= 0.8
    # No titration mixtures in this bundle -> honestly not computed.
    assert manifest["axes"]["titration_monotonicity"]["status"] == "not_computed"


@pytest.mark.skipif(
    not os.environ.get("ARIA_SEQC_MAQC_BUNDLE"),
    reason="set ARIA_SEQC_MAQC_BUNDLE to a staged real reference bundle",
)
def test_a1_seqc_maqc_real_bundle():
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.reference_seqc import run_seqc_maqc_a1_benchmark

    manifest = run_seqc_maqc_a1_benchmark(os.environ["ARIA_SEQC_MAQC_BUNDLE"])
    assert manifest["status"] in ("pass", "fail"), manifest
    assert manifest["axes"]["lfc_concordance"]["n_genes_scored"] > 100
    assert manifest["axes"]["taqman_de_detection"]["auc"] is not None
