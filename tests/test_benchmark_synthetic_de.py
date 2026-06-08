"""X6 numerical-accuracy benchmark: synthetic DE ground-truth recovery.

The simulator tests run anywhere (numpy/pandas/anndata). The recovery test
needs pydeseq2 and is skipped where it is absent (e.g. aria-env), like the
other pseudobulk-gated tests.
"""

import pytest


def test_simulator_is_deterministic_and_labels_truth():
    np = pytest.importorskip("numpy")
    pytest.importorskip("anndata")
    from aria.benchmarks.synthetic_de import simulate_pseudobulk_dataset

    ds1 = simulate_pseudobulk_dataset(n_genes=300, n_de=30, cells_per_donor=20, seed=1)
    ds2 = simulate_pseudobulk_dataset(n_genes=300, n_de=30, cells_per_donor=20, seed=1)

    # Deterministic for a fixed seed.
    assert np.array_equal(ds1.adata.X, ds2.adata.X)
    # Truth bookkeeping is consistent and disjoint.
    assert ds1.n_de == 30
    assert len(ds1.null_genes) == 300 - 30
    assert set(ds1.de_genes).isdisjoint(set(ds1.null_genes))
    # Design columns exist with the expected structure (unpaired, 6+6 donors).
    assert set(ds1.adata.obs["condition"].unique()) == {"COND_A", "COND_B"}
    assert ds1.adata.obs["donor"].nunique() == 12


def test_true_de_genes_have_separated_pseudobulk_means():
    np = pytest.importorskip("numpy")
    pytest.importorskip("anndata")
    import pandas as pd
    from aria.benchmarks.synthetic_de import simulate_pseudobulk_dataset

    ds = simulate_pseudobulk_dataset(n_genes=400, n_de=40, cells_per_donor=40, seed=3)
    X = np.asarray(ds.adata.X)
    cond = ds.adata.obs["condition"].values
    mean_a = X[cond == "COND_A"].mean(axis=0)
    mean_b = X[cond == "COND_B"].mean(axis=0)
    genes = list(ds.adata.var_names)
    de_idx = [genes.index(g) for g in ds.de_genes]
    null_idx = [genes.index(g) for g in ds.null_genes]

    # On average, true-DE genes shift far more than null genes between groups.
    eps = 1.0
    de_lfc = np.abs(np.log2((mean_b[de_idx] + eps) / (mean_a[de_idx] + eps)))
    null_lfc = np.abs(np.log2((mean_b[null_idx] + eps) / (mean_a[null_idx] + eps)))
    assert np.median(de_lfc) > np.median(null_lfc) * 2


def test_pseudobulk_simulator_tracks_true_log2fc():
    pytest.importorskip("anndata")
    from aria.benchmarks.synthetic_de import simulate_pseudobulk_dataset

    ds = simulate_pseudobulk_dataset(n_genes=180, n_de=20, cells_per_donor=12, seed=4)

    assert set(ds.true_log2fc) == set(ds.adata.var_names)
    assert all(ds.true_log2fc[g] != 0 for g in ds.de_genes)
    assert all(ds.true_log2fc[g] == 0 for g in ds.null_genes)


def test_pseudoreplication_null_has_no_true_condition_effect():
    np = pytest.importorskip("numpy")
    pytest.importorskip("anndata")
    from aria.benchmarks.synthetic_de import (
        _naive_cell_level_null_calls,
        simulate_pseudoreplication_null_dataset,
    )

    ds = simulate_pseudoreplication_null_dataset(
        n_genes=120, donors_per_condition=3, cells_per_donor=35,
        donor_gene_sd=1.0, seed=9,
    )

    assert ds.n_de == 0
    assert len(ds.null_genes) == 120
    assert set(ds.adata.obs["condition"].unique()) == {"COND_A", "COND_B"}
    assert ds.adata.obs["donor"].nunique() == 6
    assert np.allclose(list(ds.true_log2fc.values()), 0.0)

    naive = _naive_cell_level_null_calls(ds, nominal_alpha=0.05)
    assert naive["status"] == "success"
    assert naive["n_tested"] == 120
    # This is intentionally the anti-pattern: cells are not independent donors.
    assert naive["false_positive_rate"] >= 0.2, naive


def test_pseudobulk_de_recovers_ground_truth():
    """ARIA's real pseudobulk DE must recover the planted truth within
    tolerance (recall) while controlling false positives (empirical FDR)."""
    pytest.importorskip("anndata")
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.synthetic_de import run_pseudobulk_de_benchmark

    res = run_pseudobulk_de_benchmark(
        n_genes=1200, n_de=120, donors_per_condition=6, cells_per_donor=80, seed=11,
        min_recall=0.5, max_empirical_fdr=0.2,
    )
    assert res.status == "pass", res.as_dict()
    assert res.recall >= 0.5, res.as_dict()
    assert res.empirical_fdr <= 0.2, res.as_dict()


# ── W-CALIB: the same ground-truth gate for the BULK path ────────────────────

def test_bulk_simulator_is_deterministic_and_labels_truth():
    np = pytest.importorskip("numpy")
    pytest.importorskip("pandas")
    from aria.benchmarks.synthetic_de import simulate_bulk_dataset

    ds1 = simulate_bulk_dataset(n_genes=300, n_de=30, replicates_per_condition=5, seed=1)
    ds2 = simulate_bulk_dataset(n_genes=300, n_de=30, replicates_per_condition=5, seed=1)

    assert np.array_equal(ds1.counts.values, ds2.counts.values)   # deterministic
    assert ds1.n_de == 30
    assert len(ds1.null_genes) == 300 - 30
    assert set(ds1.de_genes).isdisjoint(set(ds1.null_genes))
    # genes x samples, 5 + 5 samples, integer counts.
    assert ds1.counts.shape == (300, 10)
    assert set(ds1.metadata["condition"].unique()) == {"COND_A", "COND_B"}
    assert str(ds1.counts.values.dtype).startswith("int")
    assert set(ds1.true_log2fc) == set(ds1.counts.index)
    assert all(ds1.true_log2fc[g] != 0 for g in ds1.de_genes)
    assert all(ds1.true_log2fc[g] == 0 for g in ds1.null_genes)


def test_bulk_a1_scorer_reports_frozen_axes():
    pd = pytest.importorskip("pandas")
    from aria.benchmarks.synthetic_de import (
        NegativeControlResult,
        score_bulk_de_a1,
        simulate_bulk_dataset,
    )

    ds = simulate_bulk_dataset(n_genes=60, n_de=10, replicates_per_condition=4, seed=5)
    truth = pd.Series(ds.true_log2fc, dtype=float)
    called = truth.abs().sort_values(ascending=False).head(8).index.tolist()
    results = pd.DataFrame({
        "log2FoldChange": truth,
        "padj": [0.001 if gene in called else 0.9 for gene in ds.counts.index],
    }, index=ds.counts.index)
    de_result = {
        "status": "success",
        "results": results,
        "sig_genes": called,
    }
    neg = NegativeControlResult(
        status="pass",
        false_positive_rate=0.0,
        max_false_positive_rate=0.0,
        nominal_alpha=0.05,
        n_permutations=3,
        n_tested=60,
        calls_per_permutation=[0, 0, 0],
        tolerances={},
    )

    manifest = score_bulk_de_a1(
        ds,
        de_result,
        neg,
        min_lfc_spearman=0.9,
        min_top_k_jaccard=0.2,
    )

    assert manifest["status"] == "pass"
    assert set(manifest["axes"]) == {
        "fdr_calibration",
        "lfc_concordance",
        "ranking_concordance",
        "significant_call_concordance",
    }
    assert manifest["axes"]["significant_call_concordance"]["recall"] == 0.8
    assert manifest["axes"]["fdr_calibration"]["false_positive_rate"] == 0.0
    assert manifest["axes"]["lfc_concordance"]["spearman_true_de"] >= 0.9


def test_bulk_de_recovers_ground_truth():
    """ARIA's real bulk DE (_run_deseq2, with apeGLM) must recover the planted
    truth within tolerance while controlling empirical FDR — the bulk analog of
    the pseudobulk gate. Protects the deferred P1-1 (b)/(c) changes."""
    pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.synthetic_de import run_bulk_de_benchmark

    res = run_bulk_de_benchmark(
        n_genes=1000, n_de=120, replicates_per_condition=6, seed=11,
        min_recall=0.5, max_empirical_fdr=0.2,
    )
    assert res.status == "pass", res.as_dict()
    assert res.recall >= 0.5, res.as_dict()
    assert res.empirical_fdr <= 0.2, res.as_dict()


# ── W-CALIB: label-permutation negative controls (empirical type-I) ───────────

def test_negative_control_result_shape():
    """The negative-control result serializes the calibration fields a report
    badge / methodology manifest consumes (runs without pydeseq2)."""
    from aria.benchmarks.synthetic_de import NegativeControlResult

    res = NegativeControlResult(
        status="pass", false_positive_rate=0.001, max_false_positive_rate=0.004,
        nominal_alpha=0.05, n_permutations=5, n_tested=1000,
        calls_per_permutation=[0, 1, 0, 2, 1],
        tolerances={"nominal_alpha": 0.05, "max_false_positive_rate": 0.05},
        messages=["ok"],
    )
    d = res.as_dict()
    assert d["status"] == "pass"
    assert d["false_positive_rate"] == 0.001
    assert d["nominal_alpha"] == 0.05
    assert d["n_permutations"] == 5
    assert d["calls_per_permutation"] == [0, 1, 0, 2, 1]


def test_bulk_de_negative_control_stays_quiet_under_permuted_null():
    """Permuting bulk condition labels destroys all signal; the real bulk DE
    path must NOT over-call — the mean false-positive rate must stay at/below
    the nominal alpha (empirical type-I control, the W-CALIB statement)."""
    pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.synthetic_de import run_bulk_de_negative_control

    res = run_bulk_de_negative_control(
        n_genes=800, n_de=80, replicates_per_condition=6,
        n_permutations=4, seed=11, nominal_alpha=0.05,
        max_false_positive_rate=0.05,
    )
    assert res.status == "pass", res.as_dict()
    assert res.false_positive_rate <= 0.05, res.as_dict()
    assert res.n_permutations >= 1, res.as_dict()


def test_pseudobulk_de_negative_control_stays_quiet_under_permuted_null():
    """Permuting the donor→condition map destroys all signal; the real
    pseudobulk DE path must not over-call under the null."""
    pytest.importorskip("anndata")
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.synthetic_de import run_pseudobulk_de_negative_control

    res = run_pseudobulk_de_negative_control(
        n_genes=800, n_de=80, donors_per_condition=6, cells_per_donor=60,
        n_permutations=3, seed=11, nominal_alpha=0.05,
        max_false_positive_rate=0.05,
    )
    assert res.status == "pass", res.as_dict()
    assert res.false_positive_rate <= 0.05, res.as_dict()


def test_calibration_suite_assembles_recovery_and_negative_control():
    """run_calibration_suite combines recall + empirical FDR + negative-control
    FP-rate for bulk and pseudobulk into one honest manifest the report badge
    consumes; a clean build is status=pass."""
    pytest.importorskip("anndata")
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.synthetic_de import run_calibration_suite

    # Use the full, genuinely-powered gate config (quick=False) for the pass
    # assertion; the quick path is a looser doctor smoke.
    manifest = run_calibration_suite(seed=11, quick=False)
    assert manifest["measured"] is True
    assert manifest["status"] in {"pass", "fail", "error"}
    assert set(manifest["paths"]) == {"bulk", "pseudobulk"}
    for path in ("bulk", "pseudobulk"):
        assert "recovery" in manifest["paths"][path]
        assert "negative_control" in manifest["paths"][path]
    summary = manifest["summary"]
    for key in (
        "bulk_recall", "bulk_empirical_fdr", "bulk_null_fpr",
        "pseudobulk_recall", "pseudobulk_empirical_fdr", "pseudobulk_null_fpr",
    ):
        assert key in summary
    # A clean build should pass calibration end to end: it recovers true effects
    # AND stays quiet (null false-positive rate ~ 0) under the permuted null.
    assert manifest["status"] == "pass", manifest
    assert summary["bulk_null_fpr"] <= 0.05, manifest
    assert summary["pseudobulk_null_fpr"] <= 0.05, manifest
