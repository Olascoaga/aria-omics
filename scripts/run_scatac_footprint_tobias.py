#!/usr/bin/env python3
"""scATAC P4.3 — Tn5-bias-corrected footprinting + differential TF binding (TOBIAS).

Orchestrates the dedicated TOBIAS driver
(``aria/scripts/chromatin_footprint_tobias.py``) in ``aria-tobias-env`` over an ATAC
fragments file split into pseudobulk-per-cell-type BAMs, then wraps the differential
TF-binding result in a provenance-stamped manifest (same pattern as the peak2gene and
SnapATAC2 lanes). Differential TF binding is associative (a footprint-signal difference
between two cell-type groups), never causal regulation.

No fabrication (ADR-002 / W2.2): a non-success driver result (TOBIAS / genome / motifs
absent) is recorded honestly, never an invented footprint. Data is the caller's local
fragments + genome + motif collection; nothing is downloaded.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fragments-file", required=True)
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--peaks-bed", required=True)
    p.add_argument("--motif-meme", required=True)
    p.add_argument("--barcode-groups", required=True, help="barcode<TAB>group TSV")
    p.add_argument("--group-a", required=True)
    p.add_argument("--group-b", required=True)
    p.add_argument("--work-dir", default="/tmp/p43_tobias")
    p.add_argument("--output-dir", default="docs/benchmark_results/scatac_footprint")
    p.add_argument("--manifest-name", default="p4_footprint_tobias_bindetect.json")
    p.add_argument("--tobias-env", default="aria-tobias-env")
    p.add_argument("--dataset", default=None)
    p.add_argument("--skip-run", action="store_true",
                   help="reuse an existing driver result.json in --work-dir")
    args = p.parse_args(argv)

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    result_json = work / "result.json"

    if not args.skip_run:
        subprocess.run([
            "conda", "run", "--no-capture-output", "-n", args.tobias_env,
            "python", str(ROOT / "aria" / "scripts" / "chromatin_footprint_tobias.py"),
            "--fragments-file", args.fragments_file, "--genome-fasta", args.genome_fasta,
            "--peaks-bed", args.peaks_bed, "--motif-meme", args.motif_meme,
            "--barcode-groups", args.barcode_groups,
            "--group-a", args.group_a, "--group-b", args.group_b,
            "--output-dir", str(work), "--output-json", str(result_json),
        ], check=True)

    driver = (json.loads(result_json.read_text()) if result_json.is_file()
              else {"ran": False, "reason": "no driver result.json"})

    # Keep paths/large blobs out of the committed manifest; keep the science summary.
    group_bams = driver.get("group_bams") or {}
    n_frag = {g: (group_bams.get(g) or {}).get("n_fragments") for g in group_bams}

    from aria.version import __version__, collect_version_metadata
    manifest = {
        "benchmark": "P4.3_footprint_tobias_bindetect",
        "scope": "differential_tf_binding_between_cell_types_via_tobias",
        "aria_version": __version__,
        "provenance": collect_version_metadata(),
        "dataset": args.dataset or Path(args.fragments_file).stem,
        "design": {
            "method": "TOBIAS ATACorrect -> ScoreBigwig -> BINDetect "
                      "(Tn5-bias-corrected footprinting, pseudobulk per cell type)",
            "group_a": args.group_a, "group_b": args.group_b,
            "contrast": f"{args.group_a}_vs_{args.group_b}",
            "n_fragments_per_group": n_frag,
            "motifs": "JASPAR2024 CORE vertebrates (MEME)",
            "feature_space": "called peaks (ATACorrect restricted to peaks)",
        },
        "ran": driver.get("ran", False),
        "reason": driver.get("reason"),
        "differential_summary": driver.get("differential_summary"),
        "caveats": [
            "Differential TF binding is an associative footprint-signal difference "
            "between cell-type groups, not causal regulation.",
            "Footprints are Tn5-bias-corrected (ATACorrect); single-cell footprinting "
            "is noisy so groups are pseudobulk.",
            "Cell-type groups are marker-derived from the paired RNA (benchmark prep); "
            "the driver itself is dataset-agnostic (consumes any barcode->group TSV).",
        ],
    }
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / args.manifest_name).write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    s = driver.get("differential_summary") or {}
    print(f"ran={manifest['ran']} n_sig={s.get('n_significant')} "
          f"-> wrote {out_dir / args.manifest_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
