#!/usr/bin/env python3
"""Re-stamp provenance onto the V47 bulk ATAC benchmark artifacts (F6).

Thin CLI over ``aria.benchmarks.bulk_atac_v47_provenance.stamp_artifacts``. The
scientific body of each artifact is preserved verbatim; only a fresh provenance block
(current repo state + source commit + tool versions + run date) is (re)written. This
is a re-emission stamper, not a re-computation of the DA/motif pipeline.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.bulk_atac_v47_provenance import stamp_artifacts  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bench-dir",
        default=None,
        help="Directory holding the V47 artifacts (default: docs/benchmark_results/bulk_atac).",
    )
    args = parser.parse_args(argv)
    manifest = stamp_artifacts(bench_dir=args.bench_dir, repo_root=ROOT)
    print(json.dumps(manifest, indent=2))
    return 0 if manifest["status"] == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
