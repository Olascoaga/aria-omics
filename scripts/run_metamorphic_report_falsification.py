#!/usr/bin/env python3
"""CLI runner for the 600-case deterministic Claim-Boundary & Metamorphic Report Falsification Suite."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from aria.benchmarks.metamorphic_report_falsification import (
    run_metamorphic_falsification_suite,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="docs/benchmark_results/preprint_v1/claim_4/metamorphic_600",
        help="Directory to write the JSON benchmark receipt.",
    )
    parser.add_argument(
        "--manifest-name",
        default="metamorphic_600_falsification.json",
        help="Filename for the JSON benchmark receipt.",
    )
    parser.add_argument("--seeds", type=int, default=10)
    args = parser.parse_args()

    result = run_metamorphic_falsification_suite(n_seeds=args.seeds)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / args.manifest_name
    tmp_path = out_path.with_suffix(".json.tmp")
    with open(tmp_path, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)
        fh.write("\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp_path, out_path)

    if not result.get("all_gates_pass"):
        print(f"FAIL: {len(result.get('failures', []))} metamorphic cases failed", file=sys.stderr)
        return 1
    print(
        f"PASS: {result['total_cases']} cases "
        f"({result['n_invalid_caught']}/{result['n_invalid_mutations']} invalid blocked, "
        f"{result['n_valid_published']}/{result['n_valid_controls']} valid retained)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
