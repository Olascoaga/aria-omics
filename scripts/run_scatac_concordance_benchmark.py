#!/usr/bin/env python3
"""Run the scATAC P3a external-concordance benchmark (ARIA vs SnapATAC2/...).

Scores ARIA's scATAC outputs against a reference tool's outputs (cluster ARI/NMI,
DA-peak genomic overlap, motif top-k/RBO, seed stability). Both ARIA and the
external tool's outputs are plain JSON with ``clusters`` / ``da_peaks`` /
``motifs`` (clusters must be in ALIGNED cell order across the two).

P3a is the buildable-now scoring scaffold: it CONSUMES the comparator outputs.
Producing the SnapATAC2/ArchR/Signac outputs in aria-bench-env is the infra-gated
P3b step (installing those tools + public datasets is Samael's call). With no
external outputs the run reports an honest ``not_run`` instead of a fake score.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.scripts.benchmark_scatac_concordance import benchmark_scatac_concordance  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aria-outputs-json", required=True)
    parser.add_argument("--external-outputs-json", default=None,
                        help="reference tool outputs; omit/absent -> honest not_run")
    parser.add_argument("--tool", default="SnapATAC2")
    parser.add_argument("--output-dir", default="docs/benchmark_results/scatac_concordance")
    parser.add_argument("--top-k-motifs", type=int, default=None)
    parser.add_argument("--slack", type=int, default=0)
    args = parser.parse_args(argv)

    manifest = benchmark_scatac_concordance({
        "aria_outputs_json": args.aria_outputs_json,
        "external_outputs_json": args.external_outputs_json,
        "tool": args.tool,
        "output_dir": args.output_dir,
        "top_k_motifs": args.top_k_motifs,
        "slack": args.slack,
    })
    print(json.dumps(
        {k: manifest.get(k) for k in
         ("status", "tool", "cluster_concordance", "da_peak_concordance",
          "motif_concordance", "seed_stability", "reason", "artifacts")},
        indent=2, sort_keys=True))
    return 0 if manifest.get("status") == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
