#!/usr/bin/env python
"""B4 (bulk ATAC pre-print): condition-level differential TF footprinting orchestrator.

Runs the TOBIAS backend (``aria/scripts/chromatin_footprint_tobias.py``, bulk mode) in
``aria-tobias-env`` over per-condition merged replicate BAMs, then writes a science-only
manifest (no large paths). Like the scATAC P4.3 driver, TOBIAS lives in a dedicated env,
so this is a dedicated CLI rather than live agent dispatch; a missing env / asset is an
honest ``ran: false``.

Example:
    python scripts/run_bulk_atac_footprint_tobias.py \
        --condition-bams cond_bams.json \
        --genome-fasta ~/.aria/genomes/hg38/genome.fa \
        --peaks-bed consensus_peaks.bed \
        --motif-meme JASPAR2024_CORE_vertebrates.meme \
        --group-a K562 --group-b GM12878

where cond_bams.json = {"K562": ["k562_r1.bam", "k562_r2.bam"],
                        "GM12878": ["gm_r1.bam", "gm_r2.bam"]}.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--condition-bams", required=True,
                   help="JSON {condition: [bam, ...]}")
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
    p.add_argument("--cores", type=int, default=8)
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
            "--mode", "bulk",
            "--condition-bams", args.condition_bams,
            "--genome-fasta", args.genome_fasta, "--peaks-bed", args.peaks_bed,
            "--motif-meme", args.motif_meme,
            "--group-a", args.group_a, "--group-b", args.group_b,
            "--output-dir", str(work), "--output-json", str(result_json),
            "--cores", str(args.cores),
        ], check=True)

    driver = (json.loads(result_json.read_text()) if result_json.is_file()
              else {"ran": False, "reason": "no driver result.json"})

    group_bams = driver.get("group_bams") or {}
    n_frag = {g: (group_bams.get(g) or {}).get("n_fragments") for g in group_bams}

    from aria.version import __version__, collect_version_metadata
    manifest = {
        "benchmark": "B4_bulk_footprint_tobias_bindetect",
        "scope": "differential_tf_binding_between_conditions_via_tobias",
        "aria_version": __version__,
        "provenance": collect_version_metadata(),
        "dataset": args.dataset or Path(args.condition_bams).stem,
        "design": {
            "method": "TOBIAS ATACorrect -> ScoreBigwig -> BINDetect "
                      "(Tn5-bias-corrected footprinting, per-condition merged "
                      "replicate BAMs)",
            "group_a": args.group_a, "group_b": args.group_b,
            "contrast": f"{args.group_a}_vs_{args.group_b}",
            "n_reads_per_condition": n_frag,
            "motifs": "JASPAR2024 CORE vertebrates (MEME)",
            "feature_space": "called peaks (ATACorrect restricted to peaks)",
        },
        "ran": driver.get("ran", False),
        "reason": driver.get("reason"),
        "differential_summary": driver.get("differential_summary"),
        "caveats": [
            "Differential TF binding is an associative footprint-signal difference "
            "between conditions, not causal regulation.",
            "Footprints are Tn5-bias-corrected (ATACorrect) over per-condition merged "
            "replicate BAMs.",
        ],
    }
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / args.manifest_name).write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    s = driver.get("differential_summary") or {}
    print(f"ran={manifest['ran']} n_sig={s.get('n_significant')} "
          f"n_tested={s.get('n_motifs_tested')} reason={driver.get('reason')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
