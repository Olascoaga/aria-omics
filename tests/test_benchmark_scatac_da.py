"""Guards for the scATAC DA calibration (synthetic_atac_da.py)."""

from __future__ import annotations

import pytest


def test_simulator_has_replicates_and_truth():
    pytest.importorskip("anndata")
    from aria.benchmarks.synthetic_atac_da import simulate_atac_da_dataset

    ds = simulate_atac_da_dataset(n_peaks=200, n_da=20, n_cells_per_condition=60,
                                  n_replicates_per_condition=3, seed=1)
    assert "replicate" in ds.adata.obs.columns
    assert ds.adata.obs["replicate"].nunique() == 6           # 3 donors x 2 conditions
    assert set(ds.adata.obs["condition"].unique()) == {"COND_A", "COND_B"}
    assert ds.n_da == 20
    assert set(ds.da_peaks).isdisjoint(set(ds.null_peaks))


def test_no_fabrication_without_da_fn():
    pytest.importorskip("anndata")
    from aria.benchmarks.synthetic_atac_da import run_atac_da_benchmark

    with pytest.raises(NotImplementedError):
        run_atac_da_benchmark(da_fn=None, n_peaks=100, n_da=10,
                              n_cells_per_condition=40)


def test_aria_pseudobulk_da_recovers_planted_peaks():
    """The real chromatin pseudobulk DA path (shared DESeq2 core) must recover
    planted true-DA peaks with controlled empirical FDR and signed log2FC
    concordance (D3-003)."""
    np = pytest.importorskip("numpy")
    pytest.importorskip("anndata")
    pytest.importorskip("scipy")
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.synthetic_atac_da import (
        run_atac_pseudobulk_da_benchmark,
        simulate_atac_da_dataset,
    )
    from aria.scripts.chromatin_diffacc import _pseudobulk_da

    m = run_atac_pseudobulk_da_benchmark(seed=11)
    assert m["status"] == "pass", m
    assert m["recall"] >= 0.5, m
    assert m["empirical_fdr"] <= 0.2, m
    assert m["n_false_positive"] <= 0.05 * max(m["n_called"], 1)

    # Signed LFC direction oracle (kills chromatin_diffacc log2FoldChange negation).
    ds = simulate_atac_da_dataset(
        n_peaks=300, n_da=40, n_cells_per_condition=320,
        n_replicates_per_condition=4, base_rate=0.15, seed=11,
    )
    pb = _pseudobulk_da(
        ds.adata, "condition", "replicate",
        comparisons=[{"test": "COND_B", "reference": "COND_A"}],
        padj_max=0.05, lfc_min=0.5, top_n=200, min_cells=5, min_reps=3,
        allow_mock=False, warnings=[], n_cpus=1,
    )
    assert pb["ran"] is True, pb
    comp = pb["comparisons"][0]
    assert comp["status"] == "success", comp
    tp_recs = [p for p in comp["top_peaks"] if p["peak"] in ds.da_peaks]
    assert len(tp_recs) >= 15, f"Expected >=15 true DA hits, got {len(tp_recs)}"
    sign_hits = [
        ("up" if float(p["log2fc"]) > 0 else "down") == ds.da_peaks[p["peak"]]
        for p in tp_recs
    ]
    assert sum(sign_hits) / len(sign_hits) >= 0.99, (
        f"scATAC pseudobulk DA sign concordance failed: {sum(sign_hits)}/{len(sign_hits)}"
    )
