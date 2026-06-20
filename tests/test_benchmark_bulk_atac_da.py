"""F12 (preprint audit 2026-06-19): live recall/FDR gate for the BULK ATAC DA wrapper.

The synthetic DA gate covered the scATAC pseudobulk path + the shared `_run_deseq2`
core, but NOT `chromatin_bulk_diffacc` (technical-sample -> biological-replicate
aggregation + the convergence gate). These guards close that gap: a light structural
guard (no pydeseq2) proves the simulator emits the wrapper's TSV contract and the
wrapper's own loaders accept it; a heavy guard runs the real wrapper end-to-end and
asserts recovery on planted peaks with controlled empirical FDR.
"""
from __future__ import annotations

import pytest


def test_emit_bulk_atac_tsvs_builds_valid_wrapper_contract(tmp_path):
    """The simulated dataset emits a peak x technical-sample matrix + sample metadata
    that chromatin_bulk_diffacc's own loaders accept (technical samples per donor)."""
    pytest.importorskip("anndata")
    pytest.importorskip("pandas")
    from aria.benchmarks.synthetic_atac_da import (
        _emit_bulk_atac_tsvs, simulate_atac_da_dataset,
    )
    from aria.scripts.chromatin_bulk_diffacc import (
        _load_counts_matrix, _load_sample_metadata,
    )

    ds = simulate_atac_da_dataset(n_peaks=120, n_da=12, n_cells_per_condition=60,
                                  n_replicates_per_condition=3, seed=1)
    counts_path, samples_path = _emit_bulk_atac_tsvs(
        ds, tmp_path, n_technical_per_replicate=2)

    counts = _load_counts_matrix(counts_path)        # raises if contract is wrong
    meta = _load_sample_metadata(samples_path)
    # 3 donors x 2 conditions x 2 technical samples = 12 technical samples
    assert counts.shape == (120, 12)
    assert set(meta.columns) >= {"condition", "replicate"}
    assert meta["replicate"].nunique() == 6          # 3 donors x 2 conditions
    # each biological replicate (donor) contributes exactly 2 technical samples
    assert (meta.groupby("replicate").size() == 2).all()
    assert (counts.values >= 0).all()


def test_bulk_atac_da_recovers_planted_peaks(tmp_path):
    """The REAL bulk ATAC wrapper (aggregation + replicate-gated DESeq2 + convergence
    gate) must recover planted true-DA peaks with controlled empirical FDR."""
    pytest.importorskip("anndata")
    pytest.importorskip("scipy")
    pytest.importorskip("pydeseq2")
    from aria.benchmarks.synthetic_atac_da import run_bulk_atac_da_benchmark

    m = run_bulk_atac_da_benchmark(seed=11, work_dir=str(tmp_path))
    assert m["status"] == "pass", m
    assert m["recall"] >= 0.5, m
    assert m["empirical_fdr"] <= 0.2, m
    assert m["n_false_positive"] <= 0.05 * max(m["n_called"], 1), m
    assert m["benchmark"] == "F12_bulk_atac_da"
