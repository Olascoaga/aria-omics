#!/usr/bin/env python3
"""scATAC P4.3 — Tn5-bias-corrected footprinting + differential TF binding (TOBIAS).

Orchestrates the dedicated TOBIAS driver
(``aria/scripts/chromatin_footprint_tobias.py``) in ``aria-tobias-env`` over an ATAC
fragments file split using an explicit barcode/group/replicate design, then wraps the
differential TF-binding result in a provenance-stamped manifest. A two-column barcode
map remains supported but descriptive; replicate/donor inference requires the third
column, BH across motifs, and label-permutation diagnostics.

No fabrication (ADR-002 / W2.2): a non-success driver result (TOBIAS / genome / motifs
absent) is recorded honestly, never an invented footprint. Data is the caller's local
fragments + genome + motif collection; nothing is downloaded.
"""
from __future__ import annotations

import argparse
import json
import os
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
    p.add_argument("--barcode-groups", required=True,
                   help="barcode<TAB>group[<TAB>replicate/donor] TSV")
    p.add_argument("--group-a", required=True)
    p.add_argument("--group-b", required=True)
    p.add_argument("--work-dir", default="/tmp/p43_tobias")
    p.add_argument("--output-dir", default="docs/benchmark_results/scatac_footprint")
    p.add_argument("--manifest-name", default="p4_footprint_tobias_bindetect.json")
    p.add_argument("--tobias-env", default="aria-tobias-env")
    p.add_argument("--dataset", default=None)
    p.add_argument("--rna-group-means", default=None,
                   help="optional JSON {gene: {group_a: mean, group_b: mean}} from the "
                        "paired RNA -> footprint<->RNA cross-evidence in the manifest")
    p.add_argument("--skip-run", action="store_true",
                   help="reuse an existing driver result.json in --work-dir")
    p.add_argument("--min-replicates-per-condition", type=int, default=3)
    p.add_argument("--footprint-fdr", type=float, default=0.05)
    p.add_argument("--max-label-permutations", type=int, default=100)
    p.add_argument("--skip-aggregate-plots", action="store_true")
    args = p.parse_args(argv)

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
            "--fragments-file", args.fragments_file, "--genome-fasta", args.genome_fasta,
            "--peaks-bed", args.peaks_bed, "--motif-meme", args.motif_meme,
            "--barcode-groups", args.barcode_groups,
            "--group-a", args.group_a, "--group-b", args.group_b,
            "--output-dir", str(work), "--output-json", str(result_json),
            "--min-replicates-per-condition", str(args.min_replicates_per_condition),
            "--footprint-fdr", str(args.footprint_fdr),
            "--max-label-permutations", str(args.max_label_permutations),
            *(["--skip-aggregate-plots"] if args.skip_aggregate_plots else []),
        ], check=True, env=child_env)

    driver = (json.loads(result_json.read_text()) if result_json.is_file()
              else {"ran": False, "reason": "no driver result.json"})

    # Keep paths/large blobs out of the committed manifest; keep the science summary.
    bam_info = driver.get("replicate_bams") or driver.get("group_bams") or {}
    n_frag = {g: (bam_info.get(g) or {}).get("n_fragments") for g in bam_info}

    # Optional cross-modal footprint<->RNA concordance (the governance differentiator).
    cross_evidence = None
    summary = json.loads(json.dumps(driver.get("differential_summary") or {}))
    results_table = summary.pop("results_table", None)
    summary["full_results_table_emitted"] = bool(results_table)
    inference = summary.get("inference") or {}
    if args.rna_group_means and summary.get("parsed"):
        from aria.agents.narrative.synthesis.footprint_rna import footprint_rna_concordance
        rna_means = json.loads(Path(args.rna_group_means).read_text())
        cross_evidence = footprint_rna_concordance(
            summary, rna_means, args.group_a, args.group_b)

    from aria.version import __version__, collect_version_metadata
    provenance = collect_version_metadata()
    if isinstance(provenance.get("environment"), dict):
        provenance["environment"]["conda_prefix"] = None
    manifest = {
        "benchmark": "P4.3_footprint_tobias_bindetect",
        "scope": "differential_tf_binding_between_cell_types_via_tobias",
        "aria_version": __version__,
        "provenance": provenance,
        "dataset": args.dataset or Path(args.fragments_file).stem,
        "design": {
            "method": driver.get("method"),
            "group_a": args.group_a, "group_b": args.group_b,
            "contrast": f"{args.group_a}_vs_{args.group_b}",
            "n_fragments_per_analysis_bam": n_frag,
            "inference": inference,
            "motifs": "JASPAR2024 CORE vertebrates (MEME)",
            "feature_space": "called peaks (ATACorrect restricted to peaks)",
        },
        "ran": driver.get("ran", False),
        "reason": driver.get("reason"),
        "differential_summary": summary,
        "rna_cross_evidence": cross_evidence,
        "caveats": [
            "Differential TF binding is an associative footprint-signal difference "
            "between cell-type groups, not causal regulation.",
            "Inferential calls require explicit biological-replicate/donor identity; "
            "two-column barcode groups remain descriptive pseudobulks.",
            "Cell-type groups are marker-derived from the paired RNA (benchmark prep); "
            "the driver itself is dataset-agnostic (consumes any barcode->group TSV).",
        ],
    }
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / args.manifest_name).write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    count_label = (f"n_sig={summary.get('n_significant')}" if
                   inference.get("status") == "success" else
                   f"n_ranked={summary.get('n_ranked_candidates')}")
    print(f"ran={manifest['ran']} {count_label} -> wrote "
          f"{out_dir / args.manifest_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
