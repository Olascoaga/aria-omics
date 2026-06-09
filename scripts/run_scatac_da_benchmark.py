#!/usr/bin/env python3
"""Run ARIA scATAC differential-accessibility calibration (X6/W-CALIB for chromatin).

Wires ARIA's real chromatin pseudobulk DA (_pseudobulk_da -> shared DESeq2 core)
to the synthetic ground-truth simulator and scores recall on true-DA peaks and
empirical FDR on null peaks. Validates the pseudobulk DA lane that single-sample
real inputs (e.g. HC11) cannot exercise. Needs pydeseq2 (aria-rna-env).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.synthetic_atac_da import run_atac_pseudobulk_da_benchmark  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="docs/benchmark_results/scatac_da")
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args(argv)
    manifest = run_atac_pseudobulk_da_benchmark(seed=args.seed, output_dir=args.output_dir)
    print(json.dumps({k: manifest.get(k) for k in
                      ("status", "recall", "empirical_fdr", "n_true_da",
                       "n_called", "n_false_positive", "artifacts")},
                     indent=2, sort_keys=True))
    return 0 if manifest.get("status") == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
