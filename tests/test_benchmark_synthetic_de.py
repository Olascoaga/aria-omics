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
