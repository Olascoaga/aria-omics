#!/usr/bin/env python3
"""Run ARIA Benchmark A2 external reference: ARIA vs muscat on the Kang dataset.

Two phases:
1. Export (aria-bench-env): muscat aggregates the real Kang data to per-cluster
   pseudobulk and runs its DE, exporting the matrices + reference tables.
2. Score (this process, needs pydeseq2 -> aria-rna-env): ARIA's pseudobulk DE
   runs on the same matrices and is compared to muscat per cluster.

Data-gated: if the export is absent and the R step cannot run, it skips honestly.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.reference_kang import run_kang_muscat_benchmark  # noqa: E402

R_SCRIPT = ROOT / "aria" / "scripts" / "benchmark_a2_external_muscat.R"


def _export(export_dir: Path, bench_env: str) -> bool:
    cmd = [
        "conda", "run", "--name", bench_env, "--no-capture-output",
        "Rscript", str(R_SCRIPT), str(export_dir),
    ]
    print(f"[a2] {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0 or "A2_MUSCAT_EXPORT_DONE" not in proc.stdout:
        sys.stderr.write(proc.stderr[-4000:])
        print("[a2] muscat export FAILED — is `muscat` installed in the env?")
        return False
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--export-dir",
        default=os.path.expanduser("~/.aria/benchmarks/kang_muscat"),
    )
    parser.add_argument("--bench-env", default="aria-bench-env")
    parser.add_argument("--output-dir", default="docs/benchmark_results/a2_kang_muscat")
    parser.add_argument("--refresh", action="store_true",
                        help="Re-run the muscat export even if it exists.")
    args = parser.parse_args(argv)

    export_dir = Path(os.path.expanduser(args.export_dir))
    if args.refresh or not (export_dir / "clusters.json").exists():
        if not _export(export_dir, args.bench_env):
            print(json.dumps({"status": "skipped",
                              "reason": "muscat export unavailable"}, indent=2))
            return 0

    manifest = run_kang_muscat_benchmark(export_dir, output_dir=args.output_dir)
    print(json.dumps({
        "status": manifest.get("status"),
        "summary": manifest.get("summary"),
        "dataset": manifest.get("dataset"),
        "artifacts": manifest.get("artifacts"),
        "reason": manifest.get("reason"),
    }, indent=2, sort_keys=True))
    return 0 if manifest.get("status") in ("pass", "skipped") else 1


if __name__ == "__main__":
    raise SystemExit(main())
