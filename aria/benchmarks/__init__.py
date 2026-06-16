"""Numerical-accuracy benchmarks for ARIA (audit item X6).

These benchmarks generate synthetic data with a KNOWN ground truth and check
that ARIA's real analysis code recovers it within tolerances — so a dependency
bump (pydeseq2, numpy, scanpy) that silently degrades biological accuracy is
caught, not just flow regressions.
"""

from aria.benchmarks.synthetic_de import (
    SyntheticDEDataset,
    SyntheticBulkDEDataset,
    DEBenchmarkResult,
    NegativeControlResult,
    simulate_pseudobulk_dataset,
    simulate_bulk_dataset,
    run_pseudobulk_de_benchmark,
    run_bulk_de_benchmark,
    run_bulk_de_negative_control,
    run_pseudobulk_de_negative_control,
    SpikeInDataset,
    SpikeInResult,
    simulate_spike_in_bulk_dataset,
    run_bulk_de_spike_in,
    run_calibration_suite,
)
from aria.benchmarks.synthetic_atac_da import (
    SyntheticATACDADataset,
    ATACDACaller,
    simulate_atac_da_dataset,
    run_atac_da_benchmark,
)
from aria.benchmarks.concordance_atac import (
    ConcordanceResult,
    adjusted_rand_index,
    normalized_mutual_info,
    peak_set_overlap,
    motif_concordance,
    seed_stability,
    score_atac_concordance,
)

__all__ = [
    "SyntheticDEDataset",
    "SyntheticBulkDEDataset",
    "DEBenchmarkResult",
    "simulate_pseudobulk_dataset",
    "simulate_bulk_dataset",
    "run_pseudobulk_de_benchmark",
    "run_bulk_de_benchmark",
    # W-CALIB: label-permutation negative controls + calibration manifest.
    "NegativeControlResult",
    "run_bulk_de_negative_control",
    "run_pseudobulk_de_negative_control",
    # W-CALIB: spike-in dose-response effect-size calibration ladder.
    "SpikeInDataset",
    "SpikeInResult",
    "simulate_spike_in_bulk_dataset",
    "run_bulk_de_spike_in",
    "run_calibration_suite",
    # P3-2: scATAC differential-accessibility calibration hook (v4.6 scaffold).
    "SyntheticATACDADataset",
    "ATACDACaller",
    "simulate_atac_da_dataset",
    "run_atac_da_benchmark",
    # P3a: scATAC external-concordance scoring (ARIA vs SnapATAC2/ArchR/Signac).
    "ConcordanceResult",
    "adjusted_rand_index",
    "normalized_mutual_info",
    "peak_set_overlap",
    "motif_concordance",
    "seed_stability",
    "score_atac_concordance",
]
