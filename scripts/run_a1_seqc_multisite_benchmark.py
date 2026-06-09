#!/usr/bin/env python3
"""Run ARIA Benchmark A1 cross-site SEQC reproducibility.

Runs ARIA's bulk DE (A vs B) at each SEQC sequencing site and reports the
pairwise log2FC concordance between sites (the SEQC reproducibility metric) plus
each site's TaqMan concordance. Stage per-site bundles first, e.g.:

    for S in BGI CNL MAY AGR NVS; do
      python scripts/fetch_seqc_maqc_reference.py \\
        --out ~/.aria/benchmarks/seqc_maqc_$S --count-table ILM_refseq_gene_$S
    done

Data-gated and honest-skip: sites with no bundle are skipped. Requires pydeseq2
(run in aria-rna-env).
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.reference_seqc import run_seqc_maqc_multisite  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bundle-root",
        default=os.path.expanduser("~/.aria/benchmarks"),
        help="Root holding per-site bundles named seqc_maqc_<SITE>.",
    )
    parser.add_argument("--sites", default="BGI,CNL,MAY,AGR,NVS")
    parser.add_argument("--output-dir", default="docs/benchmark_results/a1_seqc_maqc")
    args = parser.parse_args(argv)

    root = Path(os.path.expanduser(args.bundle_root))
    site_bundles = {
        s.strip(): str(root / f"seqc_maqc_{s.strip()}")
        for s in args.sites.split(",") if s.strip()
    }

    manifest = run_seqc_maqc_multisite(site_bundles, output_dir=args.output_dir)
    print(json.dumps({
        "status": manifest.get("status"),
        "cross_site": manifest.get("cross_site"),
        "per_site": {k: {"status": v.get("status"),
                         "taqman": v.get("taqman")}
                     for k, v in manifest.get("per_site", {}).items()},
        "artifacts": manifest.get("artifacts"),
    }, indent=2, sort_keys=True))
    return 0 if manifest.get("status") in ("pass", "skipped") else 1


if __name__ == "__main__":
    raise SystemExit(main())
