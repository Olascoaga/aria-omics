#!/usr/bin/env python3
"""Run ARIA Benchmark A2: preliminary donor-aware pseudobulk validation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.synthetic_de import run_pseudobulk_a2_benchmark  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="docs/benchmark_results",
        help="Directory for the A2 manifest JSON and Fig. 2 SVG.",
    )
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--quick", action="store_true", help="Use the smaller smoke config.")
    parser.add_argument("--manifest-name", default="a2_pseudobulk_v4.5.5.json")
    parser.add_argument("--figure-name", default="fig2_a2_pseudobulk_v4.5.5.svg")
    parser.add_argument("--artifact-version", default="v4.5.5")
    args = parser.parse_args(argv)

    manifest = run_pseudobulk_a2_benchmark(
        seed=args.seed,
        quick=args.quick,
        output_dir=args.output_dir,
        manifest_name=args.manifest_name,
        figure_name=args.figure_name,
        artifact_version=args.artifact_version,
    )
    print(json.dumps({
        "status": manifest.get("status"),
        "artifact_version": manifest.get("artifact_version"),
        "artifacts": manifest.get("artifacts"),
        "axes": manifest.get("axes"),
    }, indent=2, sort_keys=True))
    return 0 if manifest.get("status") == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
