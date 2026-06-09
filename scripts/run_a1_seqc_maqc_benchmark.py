#!/usr/bin/env python3
"""Run ARIA Benchmark A1 reference lane: bulk DE vs MAQC/SEQC TaqMan truth.

Data-gated: stage a reference bundle (see scripts/fetch_seqc_maqc_reference.py)
under ``ARIA_SEQC_MAQC_BUNDLE`` or pass ``--bundle``. Skips honestly when absent.
Requires pydeseq2 (run inside ``aria-rna-env``).
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.reference_seqc import run_seqc_maqc_a1_benchmark  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bundle",
        default=os.environ.get("ARIA_SEQC_MAQC_BUNDLE", ""),
        help="Reference bundle dir (counts.tsv, samples.tsv, taqman.tsv).",
    )
    parser.add_argument("--output-dir", default="docs/benchmark_results/a1_seqc_maqc")
    parser.add_argument("--numerator", default="A")
    parser.add_argument("--denominator", default="B")
    args = parser.parse_args(argv)

    if not args.bundle:
        print(json.dumps({
            "status": "skipped",
            "reason": "no bundle given (set ARIA_SEQC_MAQC_BUNDLE or --bundle)",
            "bootstrap": "scripts/fetch_seqc_maqc_reference.py",
        }, indent=2))
        return 0

    manifest = run_seqc_maqc_a1_benchmark(
        args.bundle,
        output_dir=args.output_dir,
        numerator=args.numerator,
        denominator=args.denominator,
    )
    print(json.dumps({
        "status": manifest.get("status"),
        "axes": manifest.get("axes"),
        "samples": manifest.get("samples"),
        "artifacts": manifest.get("artifacts"),
        "reason": manifest.get("reason"),
    }, indent=2, sort_keys=True))
    # Honest skip is not a failure; only fail/error are non-zero.
    return 0 if manifest.get("status") in ("pass", "skipped") else 1


if __name__ == "__main__":
    raise SystemExit(main())
