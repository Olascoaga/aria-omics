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
    simulate_pseudobulk_dataset,
    simulate_bulk_dataset,
    run_pseudobulk_de_benchmark,
    run_bulk_de_benchmark,
)
from aria.benchmarks.synthetic_atac_da import (
    SyntheticATACDADataset,
    ATACDACaller,
    simulate_atac_da_dataset,
    run_atac_da_benchmark,
)

__all__ = [
    "SyntheticDEDataset",
    "SyntheticBulkDEDataset",
    "DEBenchmarkResult",
    "simulate_pseudobulk_dataset",
    "simulate_bulk_dataset",
    "run_pseudobulk_de_benchmark",
    "run_bulk_de_benchmark",
    # P3-2: scATAC differential-accessibility calibration hook (v4.6 scaffold).
    "SyntheticATACDADataset",
    "ATACDACaller",
    "simulate_atac_da_dataset",
    "run_atac_da_benchmark",
]
