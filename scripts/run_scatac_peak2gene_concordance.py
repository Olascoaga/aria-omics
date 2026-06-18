#!/usr/bin/env python3
"""scATAC P4.2 — peak-to-gene concordance: ARIA vs Signac LinkPeaks.

Validates ARIA's paired-cell peak-to-gene linker
(``chromatin_regulatory._peak_to_gene``) against the canonical single-cell linker
(Signac ``LinkPeaks``) on the SAME paired RNA+ATAC matrix. LinkPeaks works on the
peak-count + RNA assays (NO fragments), so it is runnable on the fragment-less
GSE278576 paired data — unlike ArchR/Signac GeneActivity, which needs fragments and
is therefore deferred (the P4.3 / fragment track).

Pipeline (multi-env, orchestrated via ``conda run``):
  1. (--chromatin-env) ``aria/scripts/benchmark_peak2gene_export.py``: run ARIA
     peak-to-gene on the paired cells + export the identical matrices for Signac.
  2. (--signac-env) ``aria/scripts/benchmark_scatac_peak2gene_signac.R``: Signac
     LinkPeaks on the same matrices -> reference links TSV.
  3. (this env) ``concordance_atac.score_peak2gene_concordance`` -> link-set overlap,
     sign agreement, score Spearman -> provenance-stamped manifest.

No fabrication: a non-success ARIA/Signac stage is recorded honestly, never invented.
Data is the caller's local ``*_paired.h5mu``; nothing is downloaded.
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

from aria.benchmarks.concordance_atac import score_peak2gene_concordance  # noqa: E402


def _conda_run(env: str, cmd: list[str]) -> None:
    subprocess.run(["conda", "run", "--no-capture-output", "-n", env, *cmd], check=True)


def _read_aria_links(path: Path) -> list[tuple]:
    if not path.is_file():
        return []
    rows = []
    with path.open() as fh:
        for r in csv.DictReader(fh):
            rows.append((r["gene"], r["peak"], r["correlation"]))
    return rows


def _read_signac_links(path: Path) -> list[tuple]:
    if not path.is_file():
        return []
    rows = []
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append((r["gene"], r["peak"], r["score"]))
    return rows


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5mu", required=True, help="paired RNA+ATAC .h5mu")
    p.add_argument("--gtf", required=True, help="local gene-annotation GTF (matches peak naming)")
    p.add_argument("--work-dir", default="/tmp/p42_peak2gene")
    p.add_argument("--output-dir", default="docs/benchmark_results/scatac_concordance")
    p.add_argument("--manifest-name", default="p4_peak2gene_concordance_signac.json")
    p.add_argument("--chromatin-env", default="aria-chromatin-env")
    p.add_argument("--signac-env", default="aria-bench-signac-env")
    p.add_argument("--distance-bp", type=int, default=500_000)
    p.add_argument("--min-corr", type=float, default=0.3)
    p.add_argument("--dataset", default=None, help="dataset label for the manifest")
    p.add_argument("--skip-export", action="store_true", help="reuse an existing work dir")
    p.add_argument("--skip-signac", action="store_true", help="reuse an existing Signac links TSV")
    args = p.parse_args(argv)

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    aria_csv = work / "chromatin_peak_to_gene.csv"
    signac_tsv = work / "signac_links.tsv"
    signac_json = work / "signac_status.json"

    # 1) export + ARIA peak-to-gene (chromatin stack)
    if not args.skip_export:
        _conda_run(args.chromatin_env, [
            "python", str(ROOT / "aria" / "scripts" / "benchmark_peak2gene_export.py"),
            "--h5mu", args.h5mu, "--gtf", args.gtf, "--out-dir", str(work),
            "--distance-bp", str(args.distance_bp), "--min-corr", str(args.min_corr)])

    # 2) Signac LinkPeaks (signac stack)
    if not args.skip_signac:
        _conda_run(args.signac_env, [
            "Rscript", str(ROOT / "aria" / "scripts" / "benchmark_scatac_peak2gene_signac.R"),
            "--in-dir", str(work), "--gtf", args.gtf,
            "--out-links", str(signac_tsv), "--out-json", str(signac_json),
            "--distance", str(args.distance_bp)])

    signac_status = (json.loads(signac_json.read_text())
                     if signac_json.is_file() else {"status": "not_run", "reason": "no status json"})

    aria_links = _read_aria_links(aria_csv)
    signac_links = _read_signac_links(signac_tsv)

    if signac_status.get("status") != "success" or not signac_links:
        concordance = {"status": "not_run",
                       "reason": f"Signac reference not available: {signac_status.get('reason') or signac_status.get('status')}"}
    else:
        concordance = {"status": "success",
                       **score_peak2gene_concordance(aria_links, signac_links)}
        c = concordance
        print(f"links ARIA={c['n_aria_links']} Signac={c['n_ref_links']} "
              f"shared={c['n_shared_links']} jaccard={c['link_jaccard']:.3f} "
              f"recall(ARIA in Signac)={c['recall_aria_in_ref']:.3f} "
              f"sign={c['sign_agreement']} spearman={c['score_spearman']}", flush=True)

    from aria.version import __version__, collect_version_metadata
    manifest = {
        "benchmark": "P4.2_peak_to_gene_concordance",
        "scope": "peak_to_gene_link_concordance_vs_signac_linkpeaks",
        "aria_version": __version__,
        "provenance": collect_version_metadata(),
        "dataset": args.dataset or Path(args.h5mu).stem,
        "design": {
            "aria_method": "chromatin_regulatory._peak_to_gene "
                           "(paired-cell peak-gene Pearson within distance, min_corr)",
            "reference": "Signac::LinkPeaks (single-cell peak-gene linker; no fragments)",
            "feature_space": "identical paired peak + RNA matrices (shared barcodes)",
            "distance_bp": args.distance_bp, "min_corr": args.min_corr,
        },
        "aria_n_links": len(aria_links),
        "signac_status": signac_status,
        "concordance": concordance,
        "caveats": [
            "Peak-to-gene links are associative correlations, not regulatory proof.",
            "ArchR/Signac GeneActivity (gene scores) needs fragments and is deferred "
            "to the fragment track; only peak-to-gene is validated here.",
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
