#!/usr/bin/env python3
"""Run ARIA Benchmark A1 ERCC dose-response against a real SEQC bundle.

Scores ARIA's recovery of the ERCC spike-in design: fold-change recovery
(measured log2(A/B) vs known Mix1/Mix2 ratio, per subgroup) and dynamic-range
linearity (measured CPM vs known input concentration). Needs a bundle with
ercc_counts.tsv + ercc_truth.tsv (see scripts/fetch_seqc_maqc_reference.py).
Data-gated, honest-skip when absent.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.reference_seqc import run_ercc_dose_response  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bundle",
        default=os.environ.get(
            "ARIA_SEQC_MAQC_BUNDLE",
            os.path.expanduser("~/.aria/benchmarks/seqc_maqc_BGI"),
        ),
    )
    parser.add_argument("--output-dir", default="docs/benchmark_results/a1_seqc_maqc")
    args = parser.parse_args(argv)

    manifest = run_ercc_dose_response(args.bundle, output_dir=args.output_dir)
    print(json.dumps({
        "status": manifest.get("status"),
        "axes": manifest.get("axes"),
        "reason": manifest.get("reason"),
        "artifacts": manifest.get("artifacts"),
    }, indent=2, sort_keys=True))
    return 0 if manifest.get("status") in ("pass", "skipped") else 1


if __name__ == "__main__":
    raise SystemExit(main())
