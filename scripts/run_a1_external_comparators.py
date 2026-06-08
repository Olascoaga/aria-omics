#!/usr/bin/env python3
"""Run ARIA Benchmark A1 external R comparators through JSON IPC."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.utils.environment_manager import EnvironmentManager  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="docs/benchmark_results/a1_external")
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--timeout", type=int, default=14400)
    parser.add_argument("--manifest-name", default="a1_external_comparators_v4.5.5.json")
    args = parser.parse_args(argv)

    mgr = EnvironmentManager()
    result = mgr.run_in_stack(
        stack="benchmark",
        script_path="aria/scripts/benchmark_a1_external_comparators.py",
        params={
            "output_dir": args.output_dir,
            "seed": args.seed,
            "quick": args.quick,
            "manifest_name": args.manifest_name,
            "timeout": args.timeout,
        },
        timeout=args.timeout,
    )
    print(json.dumps({
        "status": result.get("status"),
        "external_comparator_status": result.get("external_comparator_status"),
        "benchmark": result.get("benchmark"),
        "methods": result.get("method_status"),
        "artifacts": result.get("artifacts"),
        "error_type": result.get("error_type"),
        "details": result.get("details"),
    }, indent=2, sort_keys=True))
    return 0 if result.get("status") == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())

