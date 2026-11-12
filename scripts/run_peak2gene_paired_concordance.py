#!/usr/bin/env python3
"""V48 Slice I5 — Paired Peak-to-Gene external concordance & determinism runner.

Executes end-to-end validation on real paired 10x Multiome data (`GSE278576` `hc11`):
  1. `aria/scripts/integration_paired_qc.py` (`aria-chromatin-env`) for joint
     RNA+ATAC QC and optional bijective reference `/obs` label join.
  2. Matrix export (`atac_counts.mtx`, `rna_counts.mtx`, `peaks.tsv`, `genes.tsv`,
     `barcodes.tsv`) of the post-QC standard-chromosome (`chr1..chr22,chrX,chrY`)
     paired universe so R `Signac::LinkPeaks` (`aria-bench-signac-env`) runs on the
     IDENTICAL cells, genes, and peaks.
  3. External reference execution via `aria/scripts/benchmark_scatac_peak2gene_signac.R`
     (`Signac::LinkPeaks`, `BSgenome.Hsapiens.UCSC.hg38` GC stats, `distance=500000`,
     `score_cutoff=0.05`, `pvalue_cutoff=0.05`).
  4. ARIA `aria/scripts/integration_peak2gene.py` (`aria-integration-env`) runs:
     - Signac-matched cutoffs (`min_corr=0.05`, `max_padj=None`, `max_pvalue_signac=0.05`)
       for direct apples-to-apples external concordance (`score_peak2gene_concordance`
       + Pearson/Spearman of `r` and `zscore`).
     - Production V48 `500kb` BH-FDR cutoffs (`min_corr=0.30`, `max_padj=0.05`) at
       `seed=0` x 2 (bit-exact determinism) and `seed=1` (background-sampling stability),
       plus directional recall inside `Signac::LinkPeaks`.
     - Heterogeneous-tissue `250kb` BH-FDR window (`distance_kb=250`, `min_corr=0.30`,
       `max_padj=0.05`).
     - Optional Ensembl (`1..22,X,Y`) vs UCSC (`chr1..chr22,chrX,chrY`) GTF contig
       canonicalization check when `--ensembl-gtf` is provided.

Emits a sanitized public JSON artifact with zero machine-absolute paths
(guarded by `tests/test_public_artifacts_no_absolute_paths.py`).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.concordance_atac import (  # noqa: E402
    _link_key,
    _spearman,
    score_peak2gene_concordance,
)
from aria.utils.peak2gene_core import is_standard_contig, parse_peak_coord  # noqa: E402


def _resolve_conda_bin() -> str:
    found = shutil.which("conda")
    if found:
        return found
    cand = Path(sys.executable).resolve().parents[2] / "bin" / "conda"
    if cand.is_file():
        return str(cand)
    fallback = Path.home() / "anaconda3" / "bin" / "conda"
    if fallback.is_file():
        return str(fallback)
    return "conda"


def _run_ipc_script(
    env: str,
    script_rel: str,
    params: dict[str, Any],
    work_dir: Path,
    tag: str,
) -> dict[str, Any]:
    in_json = work_dir / f"{tag}_in.json"
    out_json = work_dir / f"{tag}_out.json"
    in_json.write_text(json.dumps(params, indent=2), encoding="utf-8")
    proc = subprocess.run(
        [
            _resolve_conda_bin(),
            "run",
            "--no-capture-output",
            "-n",
            env,
            "python",
            str(ROOT / script_rel),
            str(in_json),
            str(out_json),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not out_json.is_file():
        raise RuntimeError(
            f"{script_rel} in {env} exited with {proc.returncode}:\n"
            f"STDOUT:\n{proc.stdout[-2000:]}\nSTDERR:\n{proc.stderr[-2000:]}"
        )
    return json.loads(out_json.read_text(encoding="utf-8"))


def _export_signac_mtx_inline(
    env: str,
    h5mu_path: Path,
    export_dir: Path,
) -> dict[str, int]:
    """Export post-QC `.h5mu` standard-chrom matrices for `benchmark_scatac_peak2gene_signac.R`."""
    export_dir.mkdir(parents=True, exist_ok=True)
    helper = export_dir / "_export_mtx.py"
    meta_json = export_dir / "_export_meta.json"
    helper.write_text(
        """import json, sys
from pathlib import Path
import mudata as md
import scipy.io as sio
import scipy.sparse as sp

sys.path.insert(0, sys.argv[1])
from aria.utils.mudata_io import read_h5mu_paired
from aria.utils.peak2gene_core import is_standard_contig, parse_peak_coord

h5mu_p = sys.argv[2]
out = Path(sys.argv[3])
meta_p = Path(sys.argv[4])
out.mkdir(parents=True, exist_ok=True)

paired = read_h5mu_paired(h5mu_p, require_raw_counts=True)
assert paired["status"] == "success", paired
barcodes = list(paired["shared_barcodes"])
rna = paired["rna_adata"]
atac = paired["atac_adata"]
r_map = {str(b): i for i, b in enumerate(rna.obs_names)}
a_map = {str(b): i for i, b in enumerate(atac.obs_names)}
r_mat = sp.csc_matrix(paired["rna_counts_matrix"][[r_map[b] for b in barcodes], :])
a_mat_all = sp.csc_matrix(paired["atac_counts_matrix"][[a_map[b] for b in barcodes], :])

keep_idx = []
peak_rows = []
for i, name in enumerate(map(str, atac.var_names)):
    c = parse_peak_coord(name)
    if c is None or not is_standard_contig(c[0]):
        continue
    keep_idx.append(i)
    peak_rows.append((name, c[0], c[1], c[2]))

a_mat = a_mat_all[:, keep_idx]
sio.mmwrite(str(out / "atac_counts.mtx"), a_mat.T)
sio.mmwrite(str(out / "rna_counts.mtx"), r_mat.T)
(out / "peaks.tsv").write_text("\\n".join(f"{c}\\t{s}\\t{e}\\t{n}" for (n, c, s, e) in peak_rows) + "\\n", encoding="utf-8")
(out / "genes.tsv").write_text("\\n".join(map(str, rna.var_names)) + "\\n", encoding="utf-8")
(out / "barcodes.tsv").write_text("\\n".join(barcodes) + "\\n", encoding="utf-8")
meta_p.write_text(json.dumps({"n_cells": len(barcodes), "n_genes": int(rna.n_vars), "n_peaks": len(peak_rows)}), encoding="utf-8")
""",
        encoding="utf-8",
    )
    proc = subprocess.run(
        [
            _resolve_conda_bin(),
            "run",
            "--no-capture-output",
            "-n",
            env,
            "python",
            str(helper),
            str(ROOT),
            str(h5mu_path),
            str(export_dir),
            str(meta_json),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not meta_json.is_file():
        raise RuntimeError(
            f"Export helper failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}"
        )
    return json.loads(meta_json.read_text(encoding="utf-8"))


def _read_aria_csv_detailed(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows.append(
                {
                    "gene": r["gene"],
                    "peak": r["peak"],
                    "correlation": float(r["correlation"]),
                    "zscore": float(r["zscore"]) if r.get("zscore") not in (None, "") else None,
                    "p_value": float(r["p_value"]) if r.get("p_value") not in (None, "") else None,
                    "padj": float(r["padj"]) if r.get("padj") not in (None, "") else None,
                    "p_value_signac": (
                        float(r["p_value_signac"])
                        if r.get("p_value_signac") not in (None, "")
                        else None
                    ),
                    "distance_bp": int(r["distance_bp"]),
                    "direction": r["direction"],
                    "region_class": r["region_class"],
                }
            )
    return rows


def _read_signac_tsv_detailed(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(
                {
                    "gene": r["gene"],
                    "peak": r["peak"],
                    "score": float(r["score"]),
                    "zscore": float(r["zscore"]),
                    "pvalue": float(r["pvalue"]),
                }
            )
    return rows


def _score_detailed_concordance(
    aria_rows: list[dict[str, Any]],
    signac_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    aria_tuples = [(r["gene"], r["peak"], r["correlation"]) for r in aria_rows]
    ref_tuples = [(r["gene"], r["peak"], r["score"]) for r in signac_rows]
    base = score_peak2gene_concordance(aria_tuples, ref_tuples)

    a_map = {_link_key(r["gene"], r["peak"]): r for r in aria_rows if _link_key(r["gene"], r["peak"])}
    b_map = {_link_key(r["gene"], r["peak"]): r for r in signac_rows if _link_key(r["gene"], r["peak"])}
    shared_keys = sorted(set(a_map) & set(b_map))

    if len(shared_keys) >= 2:
        a_r = np.asarray([a_map[k]["correlation"] for k in shared_keys], dtype=np.float64)
        b_r = np.asarray([b_map[k]["score"] for k in shared_keys], dtype=np.float64)
        base["score_pearson"] = (
            round(float(np.corrcoef(a_r, b_r)[0, 1]), 6)
            if np.std(a_r) > 0 and np.std(b_r) > 0
            else 1.0
        )
        base["score_mae"] = round(float(np.mean(np.abs(a_r - b_r))), 6)
        base["score_spearman"] = round(float(base["score_spearman"]), 6)
        base["sign_agreement"] = round(float(base["sign_agreement"]), 6)

        z_pairs = [
            (a_map[k]["zscore"], b_map[k]["zscore"])
            for k in shared_keys
            if a_map[k].get("zscore") is not None
        ]
        if len(z_pairs) >= 2:
            az = np.asarray([p[0] for p in z_pairs], dtype=np.float64)
            bz = np.asarray([p[1] for p in z_pairs], dtype=np.float64)
            base["zscore_pearson"] = round(float(np.corrcoef(az, bz)[0, 1]), 6)
            base["zscore_spearman"] = round(float(_spearman(az, bz)), 6)
            base["zscore_sign_agreement"] = round(float(np.mean(np.sign(az) == np.sign(bz))), 6)
    for k in ("link_jaccard", "recall_aria_in_ref", "recall_ref_in_aria"):
        if base.get(k) is not None:
            base[k] = round(float(base[k]), 6)
    return base


def _sanitize_public_obj(obj: Any) -> Any:
    """Ensure zero machine-absolute paths leak into committed benchmark JSON."""
    if isinstance(obj, dict):
        out: dict[str, Any] = {}
        for k, v in obj.items():
            if k in ("conda_prefix", "output_path", "h5mu_output_path", "input_path"):
                continue
            out[str(k)] = _sanitize_public_obj(v)
        return out
    if isinstance(obj, list):
        return [_sanitize_public_obj(x) for x in obj]
    if isinstance(obj, str):
        if obj.startswith(("/home/", "/tmp/", "/var/", "/opt/", "/mnt/", "/Users/")):
            return Path(obj).name
        return obj
    return obj


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--h5mu",
        default="/home/medusa/Samael/Oligos/data/interim/paired/hc11_paired.h5mu",
    )
    ap.add_argument(
        "--reference-h5ad",
        default="/home/medusa/Samael/Oligos/data/raw/reference/GSE278576_hippocampus_RNA.h5ad",
    )
    ap.add_argument("--sample-id", default="hc11")
    ap.add_argument(
        "--gtf",
        default="/home/medusa/Samael/H9-RNA/reference/gencode.v49.primary_assembly.annotation.gtf",
    )
    ap.add_argument(
        "--ensembl-gtf",
        default="/home/medusa/.aria/genomes/hg38/annotation.gtf",
    )
    ap.add_argument(
        "--genome-fasta",
        default="/home/medusa/.aria/genomes/hg38/GRCh38.primary_assembly.genome.fa",
    )
    ap.add_argument("--chromatin-env", default="aria-chromatin-env")
    ap.add_argument("--integration-env", default="aria-integration-env")
    ap.add_argument("--signac-env", default="aria-bench-signac-env")
    ap.add_argument(
        "--out-json",
        default="docs/benchmark_results/multiome_concordance/v48_i5_peak2gene_concordance_hc11.json",
    )
    args = ap.parse_args(argv)

    with tempfile.TemporaryDirectory(prefix="v48_i5_p2g_bench_") as td:
        work = Path(td)
        print("[1/6] Running Slice I1 integration_paired_qc...", flush=True)
        qc_params: dict[str, Any] = {
            "h5mu_path": str(args.h5mu),
            "output_dir": str(work / "qc"),
            "sample_id": args.sample_id,
            "organism": "Homo sapiens",
            "genome": "hg38",
        }
        if args.reference_h5ad and Path(args.reference_h5ad).is_file():
            qc_params.update(
                {
                    "reference_h5ad": str(args.reference_h5ad),
                    "reference_sample_col": "orig.ident",
                    "reference_label_col": "subclass",
                    "target_label_col": "subclass",
                }
            )
        qc_res = _run_ipc_script(
            args.chromatin_env,
            "aria/scripts/integration_paired_qc.py",
            qc_params,
            work,
            "i1_qc",
        )
        if qc_res.get("status") != "success":
            raise RuntimeError(f"Paired QC failed: {qc_res}")
        post_qc_h5mu = Path(qc_res["output_path"])
        print(
            f"      Post-QC cells: {qc_res['n_cells_after']}/{qc_res['n_cells_before']}",
            flush=True,
        )

        print("[2/6] Exporting identical post-QC matrices & running R Signac::LinkPeaks...", flush=True)
        export_dir = work / "signac_export"
        export_meta = _export_signac_mtx_inline(args.chromatin_env, post_qc_h5mu, export_dir)
        signac_tsv = work / "signac_links.tsv"
        signac_status_json = work / "signac_status.json"
        proc_r = subprocess.run(
            [
                _resolve_conda_bin(),
                "run",
                "--no-capture-output",
                "-n",
                args.signac_env,
                "Rscript",
                str(ROOT / "aria" / "scripts" / "benchmark_scatac_peak2gene_signac.R"),
                "--in-dir",
                str(export_dir),
                "--gtf",
                str(args.gtf),
                "--out-links",
                str(signac_tsv),
                "--out-json",
                str(signac_status_json),
                "--distance",
                "500000",
                "--score-cutoff",
                "0.05",
                "--pvalue-cutoff",
                "0.05",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc_r.returncode != 0 or not signac_status_json.is_file():
            raise RuntimeError(
                f"Signac R runner failed ({proc_r.returncode}):\n{proc_r.stdout}\n{proc_r.stderr}"
            )
        signac_status = json.loads(signac_status_json.read_text(encoding="utf-8"))
        if signac_status.get("status") != "success":
            raise RuntimeError(f"Signac LinkPeaks returned non-success: {signac_status}")
        signac_rows = _read_signac_tsv_detailed(signac_tsv)
        signac_rows_r03 = [r for r in signac_rows if abs(r["score"]) >= 0.30]
        print(
            f"      Signac LinkPeaks (|r|>=0.05, p<=0.05): {len(signac_rows)} links "
            f"(|r|>=0.30 subset: {len(signac_rows_r03)} links)",
            flush=True,
        )

        print("[3/6] Running ARIA integration_peak2gene (Signac-matched cutoffs)...", flush=True)
        aria_sig_cut = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_peak2gene.py",
            {
                "h5mu_path": str(post_qc_h5mu),
                "gtf_path": str(args.gtf),
                "genome_fasta": str(args.genome_fasta),
                "distance_kb": 500,
                "min_corr": 0.05,
                "max_padj": None,
                "max_pvalue_signac": 0.05,
                "standard_chroms_only": True,
                "seed": 0,
                "output_dir": str(work / "aria_signac_cut"),
            },
            work,
            "aria_signac_cut",
        )
        aria_sig_rows = _read_aria_csv_detailed(Path(aria_sig_cut["output_path"]))
        conc_matched = _score_detailed_concordance(aria_sig_rows, signac_rows)
        print(
            f"      Matched-cutoff concordance: ARIA={conc_matched['n_aria_links']}, "
            f"Signac={conc_matched['n_ref_links']}, Jaccard={conc_matched['link_jaccard']}, "
            f"r_Pearson={conc_matched.get('score_pearson')}, "
            f"z_Spearman={conc_matched.get('zscore_spearman')}, "
            f"sign={conc_matched.get('sign_agreement')}",
            flush=True,
        )

        print("[4/6] Running ARIA production 500kb BH-FDR (seed=0 x2 + seed=1)...", flush=True)
        base_500_params = {
            "h5mu_path": str(post_qc_h5mu),
            "gtf_path": str(args.gtf),
            "genome_fasta": str(args.genome_fasta),
            "distance_kb": 500,
            "min_corr": 0.30,
            "max_padj": 0.05,
            "standard_chroms_only": True,
        }
        aria_500_s0_a = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_peak2gene.py",
            {**base_500_params, "seed": 0, "output_dir": str(work / "aria_500_s0_a")},
            work,
            "aria_500_s0_a",
        )
        aria_500_s0_b = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_peak2gene.py",
            {**base_500_params, "seed": 0, "output_dir": str(work / "aria_500_s0_b")},
            work,
            "aria_500_s0_b",
        )
        aria_500_s1 = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_peak2gene.py",
            {**base_500_params, "seed": 1, "output_dir": str(work / "aria_500_s1")},
            work,
            "aria_500_s1",
        )
        rows_s0_a = _read_aria_csv_detailed(Path(aria_500_s0_a["output_path"]))
        rows_s0_b = _read_aria_csv_detailed(Path(aria_500_s0_b["output_path"]))
        rows_s1 = _read_aria_csv_detailed(Path(aria_500_s1["output_path"]))

        exact_same_pairs = [(r["gene"], r["peak"]) for r in rows_s0_a] == [
            (r["gene"], r["peak"]) for r in rows_s0_b
        ]
        max_delta_r = (
            float(
                np.max(
                    np.abs(
                        np.asarray([r["correlation"] for r in rows_s0_a])
                        - np.asarray([r["correlation"] for r in rows_s0_b])
                    )
                )
            )
            if rows_s0_a and exact_same_pairs
            else 0.0
        )
        max_delta_z = (
            float(
                np.max(
                    np.abs(
                        np.asarray([r["zscore"] or 0.0 for r in rows_s0_a])
                        - np.asarray([r["zscore"] or 0.0 for r in rows_s0_b])
                    )
                )
            )
            if rows_s0_a and exact_same_pairs
            else 0.0
        )
        max_delta_padj = (
            float(
                np.max(
                    np.abs(
                        np.asarray([r["padj"] or 0.0 for r in rows_s0_a])
                        - np.asarray([r["padj"] or 0.0 for r in rows_s0_b])
                    )
                )
            )
            if rows_s0_a and exact_same_pairs
            else 0.0
        )

        cross_seed = score_peak2gene_concordance(
            [(r["gene"], r["peak"], r["correlation"]) for r in rows_s0_a],
            [(r["gene"], r["peak"], r["correlation"]) for r in rows_s1],
        )
        conc_fdr_vs_signac_all = _score_detailed_concordance(rows_s0_a, signac_rows)
        conc_fdr_vs_signac_r03 = _score_detailed_concordance(rows_s0_a, signac_rows_r03)
        print(
            f"      ARIA 500kb BH-FDR: {len(rows_s0_a)} links "
            f"(pos={aria_500_s0_a['n_positive_correlations']}, "
            f"neg={aria_500_s0_a['n_negative_correlations']}), "
            f"recall_in_Signac={conc_fdr_vs_signac_all['recall_aria_in_ref']}, "
            f"cross_seed_jaccard={cross_seed['link_jaccard']:.4f}",
            flush=True,
        )

        print("[5/6] Running ARIA 250kb heterogeneous-tissue window...", flush=True)
        aria_250 = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_peak2gene.py",
            {
                **base_500_params,
                "distance_kb": 250,
                "seed": 0,
                "output_dir": str(work / "aria_250_s0"),
            },
            work,
            "aria_250_s0",
        )
        rows_250 = _read_aria_csv_detailed(Path(aria_250["output_path"]))
        set_500 = {_link_key(r["gene"], r["peak"]) for r in rows_s0_a}
        set_250 = {_link_key(r["gene"], r["peak"]) for r in rows_250}

        ensembl_check: dict[str, Any] = {"ran": False}
        if args.ensembl_gtf and Path(args.ensembl_gtf).is_file():
            print("[6/6] Verifying Ensembl ('1..22,X,Y') vs UCSC ('chr*') GTF canonicalization...", flush=True)
            aria_ens = _run_ipc_script(
                args.integration_env,
                "aria/scripts/integration_peak2gene.py",
                {
                    "h5mu_path": str(post_qc_h5mu),
                    "gtf_path": str(args.ensembl_gtf),
                    "genome_fasta": str(args.genome_fasta),
                    "distance_kb": 500,
                    "min_corr": 0.30,
                    "max_padj": 0.05,
                    "standard_chroms_only": True,
                    "seed": 0,
                    "output_dir": str(work / "aria_ensembl"),
                },
                work,
                "aria_ensembl",
            )
            rows_ens = _read_aria_csv_detailed(Path(aria_ens["output_path"]))
            # Compare on shared gene symbols between GENCODE v49 and Ensembl GTFs
            ens_genes = {r["gene"].upper() for r in rows_ens}
            genc_genes = {r["gene"].upper() for r in rows_s0_a}
            shared_sym = ens_genes & genc_genes
            conc_gtf = score_peak2gene_concordance(
                [(r["gene"], r["peak"], r["correlation"]) for r in rows_s0_a if r["gene"].upper() in shared_sym],
                [(r["gene"], r["peak"], r["correlation"]) for r in rows_ens if r["gene"].upper() in shared_sym],
            )
            ensembl_check = {
                "ran": True,
                "status": aria_ens["status"],
                "n_genes_tested_ensembl": aria_ens["n_genes_tested"],
                "n_links_ensembl": len(rows_ens),
                "shared_linked_gene_concordance": {
                    k: (round(float(v), 6) if isinstance(v, float) else v)
                    for k, v in conc_gtf.items()
                },
            }

        from aria.version import __version__, collect_version_metadata

        manifest = _sanitize_public_obj(
            {
                "benchmark": "v48_i5_paired_peak2gene_concordance",
                "slice": "V48_I5_peak2gene_paired",
                "aria_version": __version__,
                "provenance": collect_version_metadata(),
                "dataset": {
                    "accession": "GSE278576",
                    "sample_id": args.sample_id,
                    "tissue": "human_hippocampus_multiome_10x",
                    "n_cells_raw": qc_res["n_cells_before"],
                    "n_cells_post_qc": qc_res["n_cells_after"],
                    "n_genes_total": export_meta["n_genes"],
                    "n_standard_peaks": export_meta["n_peaks"],
                },
                "design": {
                    "aria_script": "aria/scripts/integration_peak2gene.py",
                    "shared_core": "aria/utils/peak2gene_core.py",
                    "reconciled_scatac_layer": "aria/scripts/chromatin_regulatory.py::_peak_to_gene",
                    "external_reference": "Signac::LinkPeaks (Seurat 5.3.0 / BSgenome.Hsapiens.UCSC.hg38)",
                    "rna_normalization": aria_500_s0_a["rna_normalization"],
                    "atac_normalization": aria_500_s0_a["atac_normalization"],
                    "min_cells_detected": aria_500_s0_a["min_cells_detected"],
                    "background_method": aria_500_s0_a["background_method"],
                    "background_match": aria_500_s0_a["background_match"],
                    "gc_matching": aria_500_s0_a["gc_matching"],
                },
                "seed_determinism": {
                    "seed": 0,
                    "exact_link_order_identical": bool(exact_same_pairs),
                    "max_abs_delta_correlation": max_delta_r,
                    "max_abs_delta_zscore": max_delta_z,
                    "max_abs_delta_padj": max_delta_padj,
                    "cross_seed_0_vs_1": {
                        k: (round(float(v), 6) if isinstance(v, float) else v)
                        for k, v in cross_seed.items()
                    },
                },
                "external_concordance_vs_signac": {
                    "signac_reference_summary": signac_status,
                    "matched_signac_cutoffs_r005_p005": conc_matched,
                    "aria_production_fdr005_r030_in_signac_all": conc_fdr_vs_signac_all,
                    "aria_production_fdr005_r030_vs_signac_r030": conc_fdr_vs_signac_r03,
                },
                "production_windows": {
                    "window_500kb": {
                        "distance_bp": 500_000,
                        "min_corr": 0.30,
                        "max_padj": 0.05,
                        "n_genes_tested": aria_500_s0_a["n_genes_tested"],
                        "n_peaks_tested": aria_500_s0_a["n_peaks_tested"],
                        "n_tested": aria_500_s0_a["n_tested"],
                        "n_fdr_significant_all_tested": aria_500_s0_a["n_fdr_significant"],
                        "n_links": aria_500_s0_a["n_links"],
                        "n_positive_correlations": aria_500_s0_a["n_positive_correlations"],
                        "n_negative_correlations": aria_500_s0_a["n_negative_correlations"],
                        "n_promoter_proximal": aria_500_s0_a["n_promoter_proximal"],
                        "n_distal_cis": aria_500_s0_a["n_distal_cis"],
                        "top_links": aria_500_s0_a["top_links"][:10],
                    },
                    "window_250kb": {
                        "distance_bp": 250_000,
                        "min_corr": 0.30,
                        "max_padj": 0.05,
                        "n_tested": aria_250["n_tested"],
                        "n_fdr_significant_all_tested": aria_250["n_fdr_significant"],
                        "n_links": aria_250["n_links"],
                        "n_positive_correlations": aria_250["n_positive_correlations"],
                        "n_negative_correlations": aria_250["n_negative_correlations"],
                        "n_promoter_proximal": aria_250["n_promoter_proximal"],
                        "n_distal_cis": aria_250["n_distal_cis"],
                        "fraction_of_250kb_links_in_500kb": (
                            round(float(len(set_250 & set_500)) / float(len(set_250)), 6)
                            if set_250
                            else 1.0
                        ),
                    },
                },
                "ensembl_vs_ucsc_gtf_canonicalization": ensembl_check,
                "governance": {
                    "associative_only": True,
                    "ctcf_validated": False,
                    "hic_corroborated": False,
                    "caveat": aria_500_s0_a["caveat"],
                },
            }
        )

        out_path = ROOT / args.out_json
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(f"Wrote sanitized concordance manifest: {args.out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
