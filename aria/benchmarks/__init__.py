"""Numerical-accuracy benchmarks for ARIA (audit item X6).

These benchmarks generate synthetic data with a KNOWN ground truth and check
that ARIA's real analysis code recovers it within tolerances — so a dependency
bump (pydeseq2, numpy, scanpy) that silently degrades biological accuracy is
caught, not just flow regressions.
"""

from aria.benchmarks.synthetic_de import (
    SyntheticDEDataset,
    DEBenchmarkResult,
    simulate_pseudobulk_dataset,
    run_pseudobulk_de_benchmark,
)

__all__ = [
    "SyntheticDEDataset",
    "DEBenchmarkResult",
    "simulate_pseudobulk_dataset",
    "run_pseudobulk_de_benchmark",
]
