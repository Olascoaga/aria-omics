#!/usr/bin/env python3
"""Run ARIA Benchmark B2 (preliminary): claim/narrative governance ablation.

Pure-Python: scores ARIA's claim governance vs an ungoverned arm over a
hand-labelled claim corpus (false-narrative + causal-overreach rates).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.governance_b2 import run_b2_benchmark  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="docs/benchmark_results/b2_claim")
    args = parser.parse_args(argv)
    manifest = run_b2_benchmark(output_dir=args.output_dir)
    print(json.dumps({"status": manifest.get("status"),
                      "summary": manifest.get("summary"),
                      "artifacts": manifest.get("artifacts")}, indent=2, sort_keys=True))
    return 0 if manifest.get("status") == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
