"""P3-2: typed scATAC differential-accessibility (DA) calibration hook.

The hook must (1) simulate a scATAC accessibility matrix with a KNOWN set of
differentially-accessible peaks, (2) score recovery (recall + empirical FDR)
with the same contract as the DE benchmarks, and (3) be an honest SCAFFOLD — it
must refuse to fabricate calibration numbers when no validated DA backend is
wired, so v4.6 is born with a real numerical regression the moment it injects a
real DA caller. Dependency-light (numpy/pandas/anndata only); no pydeseq2.
"""

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("anndata")

from aria.benchmarks import (
    SyntheticATACDADataset,
    simulate_atac_da_dataset,
    run_atac_da_benchmark,
)
from aria.benchmarks.synthetic_de import DEBenchmarkResult


def test_simulate_atac_da_dataset_has_known_truth():
    ds = simulate_atac_da_dataset(n_peaks=400, n_da=40, n_cells_per_condition=80,
                                  seed=7)
    assert isinstance(ds, SyntheticATACDADataset)
    assert ds.n_da == 40
    assert set(ds.da_peaks).isdisjoint(set(ds.null_peaks))
    assert len(ds.da_peaks) + len(ds.null_peaks) == 400
    # cells x peaks, two conditions, non-negative integer accessibility counts.
    # T10: the simulator must honor n_cells_per_condition EXACTLY even when it is
    # not divisible by the donor count (80 / 3 donors -> 26+27+27 = 80 per cond).
    assert ds.adata.shape == (160, 400)
    assert set(ds.adata.obs["condition"].unique()) == {"COND_A", "COND_B"}
    assert (ds.adata.obs["condition"] == "COND_A").sum() == 80
    assert (ds.adata.obs["condition"] == "COND_B").sum() == 80
    assert ds.adata.X.min() >= 0
    # deterministic for a fixed seed
    ds2 = simulate_atac_da_dataset(n_peaks=400, n_da=40, n_cells_per_condition=80,
                                   seed=7)
    assert ds.da_peaks == ds2.da_peaks


def test_run_atac_da_benchmark_refuses_to_fabricate_without_backend():
    ds = simulate_atac_da_dataset(n_peaks=200, n_da=20, n_cells_per_condition=40,
                                  seed=3)
    with pytest.raises(NotImplementedError):
        run_atac_da_benchmark(ds)  # no da_fn -> scaffold, must not fabricate


def test_run_atac_da_benchmark_scores_an_injected_da_caller():
    ds = simulate_atac_da_dataset(n_peaks=300, n_da=30, n_cells_per_condition=60,
                                  seed=5)

    # A perfect oracle recovers exactly the truth -> recall 1.0 / FDR 0.0.
    def oracle(dataset):
        return list(dataset.da_peaks)

    res = run_atac_da_benchmark(ds, da_fn=oracle)
    assert isinstance(res, DEBenchmarkResult)
    assert res.status == "pass"
    assert res.recall == 1.0
    assert res.empirical_fdr == 0.0
    assert res.n_true_de == 30

    # A caller that finds nothing fails the recall tolerance (still no crash).
    res_null = run_atac_da_benchmark(ds, da_fn=lambda d: [])
    assert res_null.status == "fail"
    assert res_null.recall == 0.0
