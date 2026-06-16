#!/usr/bin/env python3
"""scATAC multi-sample DA concordance — ARIA pseudobulk DA vs edgeR/limma/DESeq2.

Validates ARIA's pseudobulk CONDITION-DA lane on REAL biological replicates (the
gap a single sample cannot close; ADR-048). Per-donor scATAC peak matrices do not
share a feature space and (absent fragments) cannot be re-quantified, so the shared
space is a genomic-overlap CONSENSUS of the donors' peaks, with each donor's
pseudobulk re-mapped by peak-midpoint containment (approximation, labelled as such).

Pipeline (multi-env, orchestrated via ``conda run``):
  1. (this env) read each donor's ATAC matrix (h5py), sum over cells -> pseudobulk;
     build the consensus pseudobulk matrix; filter; write counts/metadata TSVs.
  2. (--rna-env) ARIA DESeq2 via ``scripts/aria_pseudobulk_da_from_tsv.py``.
  3. (--bench-env) ``aria/scripts/benchmark_a1_external_comparators.R`` on the SAME
     matrix (R-DESeq2 / edgeR-QLF / limma-voom).
  4. (this env) score each comparator vs ARIA with
     ``concordance_atac.score_da_lfc_concordance`` (LFC Spearman + sign agreement +
     significant-set overlap) and write a provenance-stamped manifest.

Data is the caller's local scATAC (e.g. per-donor ``*_paired.h5mu``); nothing is
downloaded. No fabrication: a failed ARIA/R step is recorded, never invented.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.consensus_pseudobulk import build_consensus_pseudobulk  # noqa: E402
from aria.benchmarks.concordance_atac import score_da_lfc_concordance  # noqa: E402


def _index_key(group) -> str:
    vi = group["var"].attrs.get("_index", "_index")
    return vi if isinstance(vi, str) else vi.decode()


def _donor_atac_pseudobulk(h5mu_path: Path, atac_mod: str) -> tuple[list[str], np.ndarray]:
    """Sum a donor's ATAC peak matrix over cells via direct h5py CSR read."""
    import h5py
    import scipy.sparse as sp
    with h5py.File(h5mu_path, "r") as f:
        g = f[f"mod/{atac_mod}"]
        names = [v.decode() if isinstance(v, bytes) else str(v)
                 for v in g[f"var/{_index_key(g)}"][:]]
        n_cells = g["X/indptr"].shape[0] - 1
        X = sp.csr_matrix((g["X/data"][:], g["X/indices"][:], g["X/indptr"][:]),
                          shape=(n_cells, len(names)))
    return names, np.asarray(X.sum(axis=0)).ravel()


def _conda_run(env: str, cmd: list[str]) -> None:
    subprocess.run(["conda", "run", "-n", env, *cmd], check=True)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", required=True, help="dir holding the per-donor h5mu files")
    p.add_argument("--h5mu-template", default="{donor}_paired.h5mu",
                   help="filename template within --data-dir; {donor} is substituted")
    p.add_argument("--young", required=True, help="comma-separated reference-group donors (COND_A)")
    p.add_argument("--old", required=True, help="comma-separated test-group donors (COND_B)")
    p.add_argument("--atac-mod", default="atac", help="ATAC modality key inside the h5mu")
    p.add_argument("--output-dir", default="docs/benchmark_results/scatac_concordance")
    p.add_argument("--manifest-name", default="p3_multisample_aging_da_concordance.json")
    p.add_argument("--work-dir", default="/tmp/p3_multisample")
    p.add_argument("--rna-env", default="aria-rna-env")
    p.add_argument("--bench-env", default="aria-bench-env")
    p.add_argument("--min-detect", type=int, default=3, help="keep peaks detected in >= this many donors")
    p.add_argument("--min-total", type=int, default=10, help="keep peaks with total count >= this")
    p.add_argument("--padj", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=0.5)
    args = p.parse_args(argv)

    import pandas as pd

    young = [d.strip() for d in args.young.split(",") if d.strip()]
    old = [d.strip() for d in args.old.split(",") if d.strip()]
    donors = young + old
    data_dir = Path(args.data_dir)
    work = Path(args.work_dir); work.mkdir(parents=True, exist_ok=True)

    # 1) per-donor pseudobulk -> consensus matrix
    donor_pb: dict[str, tuple[list[str], np.ndarray]] = {}
    for d in donors:
        path = data_dir / args.h5mu_template.format(donor=d)
        names, vec = _donor_atac_pseudobulk(path, args.atac_mod)
        donor_pb[d] = (names, vec)
        print(f"{d}: {len(names)} peaks, depth {vec.sum():.0f}", flush=True)
    mat, consensus_ids, order = build_consensus_pseudobulk(donor_pb)

    counts = pd.DataFrame(np.rint(mat).astype(int), index=consensus_ids, columns=order)
    keep = ((counts > 0).sum(axis=1) >= args.min_detect) & (counts.sum(axis=1) >= args.min_total)
    counts = counts[keep]
    print(f"consensus peaks {len(consensus_ids)} -> {len(counts)} after filter", flush=True)

    cond = {d: "COND_A" for d in young} | {d: "COND_B" for d in old}
    counts_tsv = work / "counts.tsv"; meta_tsv = work / "metadata.tsv"
    cc = counts.copy(); cc.insert(0, "gene", cc.index); cc.to_csv(counts_tsv, sep="\t", index=False)
    pd.DataFrame({"sample": order, "condition": [cond[d] for d in order]}).to_csv(
        meta_tsv, sep="\t", index=False)

    # 2) ARIA pseudobulk DA (RNA stack) on the consensus matrix
    aria_json = work / "aria_da.json"
    _conda_run(args.rna_env, [
        "python", str(ROOT / "scripts" / "aria_pseudobulk_da_from_tsv.py"),
        "--counts", str(counts_tsv), "--metadata", str(meta_tsv),
        "--output-json", str(aria_json), "--padj", str(args.padj), "--lfc", str(args.lfc)])
    aria = json.loads(aria_json.read_text())
    if aria.get("status") != "success":
        print("ARIA DA did not succeed:", aria.get("error_type"), aria.get("details"))
        return 1

    # 3) R comparators (benchmark stack) on the SAME matrix
    r_json = work / "r_comparators.json"
    _conda_run(args.bench_env, [
        "Rscript", str(ROOT / "aria" / "scripts" / "benchmark_a1_external_comparators.R"),
        "--counts", str(counts_tsv), "--metadata", str(meta_tsv),
        "--output-dir", str(work / "r_out"), "--output-json", str(r_json),
        "--padj", str(args.padj), "--lfc", str(args.lfc)])
    rj = json.loads(r_json.read_text())

    # 4) score each comparator vs ARIA
    from aria.version import __version__, collect_version_metadata
    methods: dict[str, dict] = {}
    for name, m in (rj.get("methods") or {}).items():
        if m.get("status") != "success":
            methods[name] = {"status": m.get("status"), "error": m.get("details")}
            continue
        rdf = pd.read_csv(m["results_tsv"], sep="\t")
        cols = {c.lower(): c for c in rdf.columns}
        gcol = cols.get("gene") or rdf.columns[0]
        lcol = cols.get("log2foldchange") or cols.get("logfc")
        pcol = cols.get("padj") or cols.get("fdr") or cols.get("adj.p.val")
        rdf[gcol] = rdf[gcol].astype(str)
        r_lfc = dict(zip(rdf[gcol], rdf[lcol]))
        r_sig = list(rdf[(rdf[pcol] < args.padj) & (rdf[lcol].abs() >= args.lfc)][gcol])
        methods[name] = score_da_lfc_concordance(
            aria["lfc_by_peak"], r_lfc, aria_sig=aria["sig_peaks"], ref_sig=r_sig)
        s = methods[name]
        print(f"{name}: spearman={s['lfc_spearman']:.3f} sign={s['lfc_sign_agreement']:.3f} "
              f"recall(ARIA in {name})={s.get('recall_aria_in_ref')}", flush=True)

    manifest = {
        "benchmark": "P3_multisample_DA_concordance",
        "scope": "external_DE_framework_concordance_on_multireplicate_pseudobulk",
        "aria_version": __version__,
        "provenance": collect_version_metadata(),
        "design": {
            "reference_group_COND_A": young, "test_group_COND_B": old,
            "n_per_group": {"COND_A": len(young), "COND_B": len(old)},
            "feature_space": "genomic-overlap consensus peaks; pseudobulk re-mapped "
                             "by peak-midpoint containment (no fragments)",
            "consensus_peaks": len(consensus_ids), "tested_after_filter": int(len(counts)),
            "thresholds": {"padj": args.padj, "lfc": args.lfc},
            "aria_da_core": "rna_bulk_de._run_deseq2 (pydeseq2)",
            "comparators": "benchmark_a1_external_comparators.R",
        },
        "aria_n_tested": aria.get("n_tested"), "aria_n_sig": aria.get("n_sig"),
        "methods": methods,
        "caveats": [
            "Validates the pseudobulk CONDITION-DA path (donor-level, biological "
            "replication, min_replicates_per_condition>=3). NOT the per-cluster "
            "wilcoxon marker path (single-sample fragile).",
            "Consensus is a genomic-overlap merge with midpoint-containment count "
            "re-mapping (approximation vs fragment-level re-quantification).",
        ],
    }
    out_dir = Path(args.output_dir); out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / args.manifest_name).write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(f"wrote {out_dir / args.manifest_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
