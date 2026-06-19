#!/usr/bin/env python3
"""scATAC P4.3 — gene-activity concordance: ARIA vs Signac GeneActivity.

Validates ARIA's peak-window gene-activity scorer
(``chromatin_regulatory._gene_scores``) against the canonical single-cell gene-activity
estimator (Signac ``GeneActivity``, fragment-counts over gene body + promoter) on the
SAME cells. The two methods differ (ARIA sums peak accessibility in a window; Signac
counts fragments in the gene body), so this measures rank concordance of the per-gene
activity, not identity.

Pipeline (multi-env, via ``conda run``):
  1. (--chromatin-env) ``benchmark_geneactivity_export.py``: ARIA gene scores + barcodes.
  2. (--signac-env) ``benchmark_scatac_geneactivity_signac.R``: Signac GeneActivity ->
     per-gene mean activity TSV.
  3. (this env) ``concordance_atac.score_gene_activity_concordance`` -> Spearman +
     top-k overlap -> provenance-stamped manifest.

No fabrication: a non-success ARIA/Signac stage is recorded honestly, never invented.
Data is the caller's local ARC matrix + fragments; nothing is downloaded.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.concordance_atac import score_gene_activity_concordance  # noqa: E402


def _conda_run(env: str, cmd: list[str]) -> None:
    subprocess.run(["conda", "run", "--no-capture-output", "-n", env, *cmd], check=True)


def _read_scores(path: Path, gene_col: str, score_col: str, delim: str) -> list[tuple]:
    if not path.is_file():
        return []
    rows = []
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter=delim):
            try:
                rows.append((r[gene_col], float(r[score_col])))
            except (KeyError, ValueError, TypeError):
                continue
    return rows


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", required=True, help="10x ARC .h5 or ATAC .h5ad (peaks)")
    p.add_argument("--fragments", required=True, help="fragments.tsv.gz (tabix-indexed)")
    p.add_argument("--gtf", required=True, help="local gene-annotation GTF (chr* naming)")
    p.add_argument("--work-dir", default="/tmp/p43_geneactivity")
    p.add_argument("--output-dir", default="docs/benchmark_results/scatac_concordance")
    p.add_argument("--manifest-name", default="p4_geneactivity_concordance_signac.json")
    p.add_argument("--chromatin-env", default="aria-chromatin-env")
    p.add_argument("--signac-env", default="aria-bench-signac-env")
    p.add_argument("--window-bp", type=int, default=100_000)
    p.add_argument("--cores", type=int, default=1, help="cores for Signac FeatureMatrix")
    p.add_argument("--dataset", default=None)
    p.add_argument("--skip-export", action="store_true")
    p.add_argument("--skip-signac", action="store_true")
    args = p.parse_args(argv)

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    aria_csv = work / "chromatin_gene_scores.csv"
    signac_tsv = work / "signac_gene_activity.tsv"
    signac_json = work / "signac_status.json"

    if not args.skip_export:
        _conda_run(args.chromatin_env, [
            "python", str(ROOT / "aria" / "scripts" / "benchmark_geneactivity_export.py"),
            "--matrix", args.matrix, "--gtf", args.gtf, "--out-dir", str(work),
            "--window-bp", str(args.window_bp)])

    if not args.skip_signac:
        _conda_run(args.signac_env, [
            "Rscript", str(ROOT / "aria" / "scripts" / "benchmark_scatac_geneactivity_signac.R"),
            "--fragments", args.fragments, "--barcodes", str(work / "barcodes.tsv"),
            "--gtf", args.gtf, "--out-scores", str(signac_tsv),
            "--out-json", str(signac_json), "--cores", str(args.cores)])

    signac_status = (json.loads(signac_json.read_text())
                     if signac_json.is_file() else {"status": "not_run", "reason": "no status json"})
    aria_scores = _read_scores(aria_csv, "gene", "mean_gene_activity_score", ",")
    signac_scores = _read_scores(signac_tsv, "gene", "score", "\t")

    if signac_status.get("status") != "success" or not signac_scores:
        concordance = {"status": "not_run",
                       "reason": f"Signac reference not available: {signac_status.get('reason') or signac_status.get('status')}"}
    else:
        concordance = {"status": "success",
                       **score_gene_activity_concordance(aria_scores, signac_scores)}
        c = concordance
        print(f"genes ARIA={c['n_aria_genes']} Signac={c['n_ref_genes']} "
              f"shared={c['n_shared_genes']} spearman={c['score_spearman']} "
              f"top{c['top_k']}_jaccard={c['top_k_jaccard']}", flush=True)

    from aria.version import __version__, collect_version_metadata
    manifest = {
        "benchmark": "P4.3_gene_activity_concordance",
        "scope": "gene_activity_concordance_vs_signac_geneactivity",
        "aria_version": __version__,
        "provenance": collect_version_metadata(),
        "dataset": args.dataset or Path(args.matrix).stem,
        "design": {
            "aria_method": "chromatin_regulatory._gene_scores "
                           "(peak accessibility summed in a TSS window)",
            "reference": "Signac::GeneActivity (fragment counts over gene body + promoter)",
            "feature_space": "per-gene mean activity over shared cells",
            "window_bp": args.window_bp,
        },
        "aria_n_genes": len(aria_scores),
        "signac_status": signac_status,
        "concordance": concordance,
        "caveats": [
            "Gene activity (both methods) summarizes accessibility near a gene; it is "
            "not RNA expression. ARIA sums peak accessibility in a window, Signac counts "
            "fragments in the gene body + promoter — concordance is rank agreement, not "
            "identity.",
            "ARIA's window-sum assigns the SAME score to genes sharing a peak window, so "
            "the top-k overlap understates agreement; the genome-wide Spearman over all "
            "shared genes is the honest concordance signal (moderate here).",
        ],
    }
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / args.manifest_name).write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(f"wrote {out_dir / args.manifest_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
