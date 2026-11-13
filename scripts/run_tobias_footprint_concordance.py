#!/usr/bin/env python3
"""ADR-078 TOBIAS differential-footprint concordance benchmark runner.

Runs ``aria.benchmarks.tobias_footprint_concordance.run_benchmark`` (external R
``stats`` in ``aria-bench-env`` + live TOBIAS 0.17.3 in ``aria-tobias-env``) and
writes the path-clean JSON receipt. Exits non-zero when any pre-registered gate
fails, so a regression cannot be committed silently.

Usage:
  python scripts/run_tobias_footprint_concordance.py \
      --output docs/benchmark_results/bulk_atac_footprint/adr078_tobias_replicate_concordance.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.tobias_footprint_concordance import run_benchmark  # noqa: E402

DEFAULT_OUTPUT = ("docs/benchmark_results/bulk_atac_footprint/"
                  "adr078_tobias_replicate_concordance.json")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                        help="receipt JSON path (repo-relative or absolute)")
    parser.add_argument("--workdir", default=None,
                        help="keep intermediate files here (default: temp dir)")
    parser.add_argument("--bench-env", default="aria-bench-env")
    parser.add_argument("--tobias-env", default="aria-tobias-env")
    parser.add_argument("--skip-live-tobias", action="store_true",
                        help="skip suite 5 (the receipt will then not pass)")
    args = parser.parse_args(argv)

    receipt = run_benchmark(args.workdir, bench_env=args.bench_env,
                            tobias_env=args.tobias_env,
                            run_live_tobias=not args.skip_live_tobias)
    output = Path(args.output)
    if not output.is_absolute():
        output = ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    for name, suite in receipt["suites"].items():
        print(f"{name:22s} pass={suite.get('pass')} gates={suite.get('gates')}")
    print(f"paths_clean={receipt['paths_clean']} "
          f"all_gates_pass={receipt['all_gates_pass']}")
    try:
        shown = output.relative_to(ROOT)
    except ValueError:
        shown = output.name
    print(f"receipt -> {shown}")
    return 0 if receipt["all_gates_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
