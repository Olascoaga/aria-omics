#!/usr/bin/env python3
"""Run ARIA Benchmark B1: DesignAgent governance on the adversarial corpus.

Pure-Python: drives ARIA's real readiness audit + design-matrix validator over a
labelled adversarial design corpus and reports correct inference/refusal rates
and the headline unsafe-execution rate (target ~0).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.governance_b1 import run_b1_benchmark  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="docs/benchmark_results/b1_design")
    args = parser.parse_args(argv)

    manifest = run_b1_benchmark(output_dir=args.output_dir)
    print(json.dumps({
        "status": manifest.get("status"),
        "summary": manifest.get("summary"),
        "confusion_matrix": manifest.get("confusion_matrix"),
        "artifacts": manifest.get("artifacts"),
    }, indent=2, sort_keys=True))
    return 0 if manifest.get("status") == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
