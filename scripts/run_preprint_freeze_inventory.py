#!/usr/bin/env python3
"""Build the fail-closed Claim 1-7 inventory for the preprint-v1 freeze."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.preprint_freeze import build_inventory, write_inventory  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        default="docs/benchmark_results/preprint_v1",
    )
    parser.add_argument("--inventory-name", default="inventory.json")
    args = parser.parse_args(argv)
    output_root = Path(args.output_root)
    payload = build_inventory(ROOT, output_root)
    path = write_inventory(payload, output_root / args.inventory_name)
    print(json.dumps({
        "inventory": path.as_posix(),
        "freeze_ready": payload["freeze_gate"]["ready"],
        "n_required_lanes": payload["freeze_gate"]["n_required_lanes"],
        "n_verified_lanes": payload["freeze_gate"]["n_verified_lanes"],
        "n_blockers": len(payload["freeze_gate"]["blockers"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
