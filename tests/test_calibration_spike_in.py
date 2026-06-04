"""W-CALIB spike-ins (residual closure): a dose-response effect-size calibration.

Recovery + negative controls answer "does ARIA find true effects?" and "does it
stay quiet under the null?". Spike-ins add the third calibration question — "are
effect sizes recovered across a LADDER of known fold-changes, and are level-0
(null) spike-ins kept below alpha?" — an ERCC-style ground truth. The simulator
is dependency-light; running it on the REAL bulk DE path needs pydeseq2.
"""

import importlib

import pytest

from aria.benchmarks.synthetic_de import (
    simulate_spike_in_bulk_dataset,
    SpikeInDataset,
)

pydeseq2 = importlib.util.find_spec("pydeseq2")
requires_pydeseq2 = pytest.mark.skipif(pydeseq2 is None, reason="needs pydeseq2")


# ── simulator (runs anywhere) ────────────────────────────────────────────────

def test_spike_in_dataset_is_deterministic_and_well_formed():
    levels = (0.0, 1.0, 2.0)
    ds = simulate_spike_in_bulk_dataset(
        levels=levels, genes_per_level=8, n_background=200,
        replicates_per_condition=4, seed=17)
    assert isinstance(ds, SpikeInDataset)
    # genes_per_level spike genes at each level; level 0 has |log2fc| == 0.
    for lvl in levels:
        at = [g for g, v in ds.spike_level.items() if v == lvl]
        assert len(at) == 8
    assert all(v == 0.0 for g, v in ds.spike_true_log2fc.items()
               if ds.spike_level[g] == 0.0)
    # Determinism.
    ds2 = simulate_spike_in_bulk_dataset(
        levels=levels, genes_per_level=8, n_background=200,
        replicates_per_condition=4, seed=17)
    assert (ds.counts.values == ds2.counts.values).all()


# ── real DE path (pydeseq2) ──────────────────────────────────────────────────

@requires_pydeseq2
def test_spike_in_recovery_is_a_dose_response_and_nulls_stay_quiet():
    from aria.benchmarks.synthetic_de import run_bulk_de_spike_in, SpikeInResult
    res = run_bulk_de_spike_in(
        levels=(0.0, 1.0, 2.0, 3.0), genes_per_level=12, n_background=800,
        replicates_per_condition=6, seed=17)
    assert isinstance(res, SpikeInResult)
    assert res.status in {"pass", "fail"}
    rate = res.detection_rate_by_level
    # Level-0 spike-ins are true nulls -> detection at or below the FP tolerance.
    assert res.null_spike_fpr <= res.tolerances["max_null_fpr"]
    # The strongest level is recovered much better than the null level.
    assert rate[str(3.0)] >= rate[str(0.0)]
    assert rate[str(3.0)] >= res.tolerances["min_top_detection"]
    # Effect sizes are recovered within tolerance (apeGLM shrinkage accepted).
    assert res.lfc_mae <= res.tolerances["max_lfc_mae"]
    assert res.status == "pass"


@requires_pydeseq2
def test_calibration_suite_includes_spike_in_block():
    from aria.benchmarks.synthetic_de import run_calibration_suite
    suite = run_calibration_suite(quick=True, seed=17)
    assert "spike_in" in suite["paths"]["bulk"]
    assert "bulk_spike_null_fpr" in suite["summary"]
    assert "bulk_spike_top_detection" in suite["summary"]
