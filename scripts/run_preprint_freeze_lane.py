#!/usr/bin/env python3
"""Execute one registered preprint-v1 lane and write its hashed receipt."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.preprint_freeze import execute_lane  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lane_id")
    parser.add_argument(
        "--output-root", default="docs/benchmark_results/preprint_v1"
    )
    parser.add_argument("--timeout", type=int, default=14_400)
    args = parser.parse_args(argv)
    receipt = execute_lane(
        ROOT, args.output_root, args.lane_id, timeout=args.timeout
    )
    print(json.dumps({
        "lane_id": receipt["lane_id"],
        "status": receipt["status"],
        "git_commit": receipt["git_commit"],
        "test_summary": receipt["test_summary"],
        "artifacts": receipt["artifacts"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
