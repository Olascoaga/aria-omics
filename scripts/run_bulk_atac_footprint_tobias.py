#!/usr/bin/env python
"""B4 (bulk ATAC pre-print): condition-level differential TF footprinting orchestrator.

Runs the TOBIAS backend (``aria/scripts/chromatin_footprint_tobias.py``, bulk mode) in
``aria-tobias-env`` over explicit biological-replicate BAMs, then writes a science-only
manifest (no large paths). Replicate mean scores are tested with Welch tests, corrected
with BH across motifs, and accompanied by label-permutation diagnostics. A legacy
condition-only BAM map remains descriptive.

Example:
    python scripts/run_bulk_atac_footprint_tobias.py \
        --replicate-bams replicate_bams.json \
        --genome-fasta ~/.aria/genomes/hg38/genome.fa \
        --peaks-bed consensus_peaks.bed \
        --motif-meme JASPAR2024_CORE_vertebrates.meme \
        --group-a K562 --group-b GM12878

where replicate_bams.json = {
  "K562": {"rep1": ["k562_r1.bam"], "rep2": ["k562_r2.bam"]},
  "GM12878": {"rep1": ["gm_r1.bam"], "rep2": ["gm_r2.bam"]}
}.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--condition-bams",
                   help="legacy descriptive JSON {condition: [bam, ...]}")
    p.add_argument("--replicate-bams",
                   help="inferential JSON {condition: {biological_replicate: [bam]}}")
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--peaks-bed", required=True)
    p.add_argument("--motif-meme", required=True)
    p.add_argument("--group-a", required=True)
    p.add_argument("--group-b", required=True)
    p.add_argument("--work-dir", default="/tmp/b4_bulk_tobias")
    p.add_argument("--output-dir",
                   default="docs/benchmark_results/bulk_atac_footprint")
    p.add_argument("--manifest-name",
                   default="b4_bulk_footprint_tobias_bindetect.json")
    p.add_argument("--tobias-env", default="aria-tobias-env")
    p.add_argument("--dataset", default=None)
    p.add_argument("--region-scope", default="genome_wide",
                   help="declared genomic scope of the supplied peak universe")
    p.add_argument("--cores", type=int, default=8)
    p.add_argument("--min-replicates-per-condition", type=int, default=3)
    p.add_argument("--footprint-fdr", type=float, default=0.05)
    p.add_argument("--max-label-permutations", type=int, default=100)
    p.add_argument("--skip-aggregate-plots", action="store_true")
    p.add_argument("--skip-run", action="store_true",
                   help="reuse an existing driver result.json in --work-dir")
    args = p.parse_args(argv)
    if not args.condition_bams and not args.replicate_bams:
        p.error("--replicate-bams or --condition-bams is required")

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    result_json = work / "result.json"

    if not args.skip_run:
        child_env = dict(os.environ)
        child_env["PYTHONPATH"] = os.pathsep.join(
            part for part in (str(ROOT), child_env.get("PYTHONPATH")) if part)
        subprocess.run([
            "conda", "run", "--no-capture-output", "-n", args.tobias_env,
            "python", str(ROOT / "aria" / "scripts" / "chromatin_footprint_tobias.py"),
            "--mode", "bulk",
            "--genome-fasta", args.genome_fasta, "--peaks-bed", args.peaks_bed,
            "--motif-meme", args.motif_meme,
            "--group-a", args.group_a, "--group-b", args.group_b,
            "--output-dir", str(work), "--output-json", str(result_json),
            "--cores", str(args.cores),
            "--min-replicates-per-condition", str(args.min_replicates_per_condition),
            "--footprint-fdr", str(args.footprint_fdr),
            "--max-label-permutations", str(args.max_label_permutations),
            *(["--skip-aggregate-plots"] if args.skip_aggregate_plots else []),
            *(["--condition-bams", args.condition_bams] if args.condition_bams else []),
            *(["--replicate-bams", args.replicate_bams] if args.replicate_bams else []),
        ], check=True, env=child_env)

    driver = (json.loads(result_json.read_text()) if result_json.is_file()
              else {"ran": False, "reason": "no driver result.json"})

    bam_info = driver.get("replicate_bams") or driver.get("group_bams") or {}
    n_frag = {g: (bam_info.get(g) or {}).get("n_fragments") for g in bam_info}
    summary = json.loads(json.dumps(driver.get("differential_summary") or {}))
    results_table = summary.pop("results_table", None)
    summary["full_results_table_emitted"] = bool(results_table)
    inference = summary.get("inference") or {}

    from aria.version import __version__, collect_version_metadata
    provenance = collect_version_metadata()
    if isinstance(provenance.get("environment"), dict):
        provenance["environment"]["conda_prefix"] = None
    manifest = {
        "benchmark": "B4_bulk_footprint_tobias_bindetect",
        "scope": "differential_tf_binding_between_conditions_via_tobias",
        "aria_version": __version__,
        "provenance": provenance,
        "dataset": args.dataset or Path(args.replicate_bams or args.condition_bams).stem,
        "design": {
            "method": driver.get("method"),
            "group_a": args.group_a, "group_b": args.group_b,
            "contrast": f"{args.group_a}_vs_{args.group_b}",
            "n_reads_per_analysis_bam": n_frag,
            "inference": inference,
            "motifs": "JASPAR2024 CORE vertebrates (MEME)",
            "feature_space": "called peaks (ATACorrect restricted to peaks)",
            "region_scope": args.region_scope,
        },
        "validation_scope": {
            "replicate_inference": True,
            "genome_wide_sensitivity": args.region_scope == "genome_wide",
            "note": ("A restricted region scope validates replicate inference, null "
                     "controls, and real-data plumbing; it does not estimate "
                     "genome-wide sensitivity." if args.region_scope != "genome_wide"
                     else "Genome-wide supplied peak universe."),
        },
        "ran": driver.get("ran", False),
        "reason": driver.get("reason"),
        "differential_summary": summary,
        "caveats": [
            "Differential TF binding is an associative footprint-signal difference "
            "between conditions, not causal regulation.",
            "Inferential calls require explicit biological-replicate identity; legacy "
            "per-condition pools remain descriptive.",
        ],
    }
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / args.manifest_name).write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    count_label = (f"n_sig={summary.get('n_significant')}" if
                   inference.get("status") == "success" else
                   f"n_ranked={summary.get('n_ranked_candidates')}")
    print(f"ran={manifest['ran']} {count_label} "
          f"n_tested={summary.get('n_motifs_tested')} reason={driver.get('reason')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
