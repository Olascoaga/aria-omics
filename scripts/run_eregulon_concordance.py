#!/usr/bin/env python3
"""V48 Slice I6 — Paired Candidate eRegulons external concordance & determinism runner.

Executes end-to-end real-data validation on paired 10x Multiome (`GSE278576` `hc11`):
  1. `aria/scripts/integration_paired_qc.py` (`aria-chromatin-env`) for joint
     RNA+ATAC QC and bijective `/obs` cell-type label join.
  2. `aria/scripts/integration_eregulons.py` (`aria-integration-env`) at:
     - `motif_pvalue=5e-5`, `seed=0` x 2 (bit-exact determinism check)
     - `motif_pvalue=5e-5`, `seed=1` (null-draw stability check)
     - `motif_pvalue=1e-4`, `seed=0` (PWM threshold sensitivity check)
  3. Exports the identical post-QC annotated expressed gene universe, standard-chrom
     accessible peak universe (`chr:start-end`), GTF gene annotations, `.fai`
     chromsizes, and fixed-500bp JASPAR2024 `5e-5` motif cistromes so external
     `scenicplus==1.0.0` (`scenicplus_env`) runs offline on the exact same inputs.
  4. Runs `aria/scripts/benchmark_eregulons_scenicplus_ref.py` across `seed=0` and
     `seed=1` to quantify both SCENIC+'s multi-seed GBM+GSEA stochasticity ceiling
     and ARIA vs SCENIC+ 4-tier concordance (TFs, modules, target edges, cis
     triplets, and per-cell AUCell Spearman vs random-gene-set baseline).

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
import time
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.benchmarks.concordance_eregulons import (  # noqa: E402
    canonical_sign_class,
    score_eregulon_concordance,
    score_eregulon_determinism,
    spearman_corr,
)

_FORBIDDEN_PATH_PREFIXES = (
    "/home/",
    "/Users/",
    "/root/",
    "/tmp/",
    "/mnt/",
    "/private/",
    "/data/",
    "/var/",
    "/opt/",
)


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
            f"STDOUT:\n{proc.stdout[-2500:]}\nSTDERR:\n{proc.stderr[-2500:]}"
        )
    return json.loads(out_json.read_text(encoding="utf-8"))


def _read_csv_rows(path: Path, delimiter: str = ",") -> list[dict[str, Any]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def _export_scenicplus_bundle(
    integration_env: str,
    h5mu_path: Path,
    gtf_path: Path,
    fasta_path: Path,
    motif_collection: str,
    motif_pvalue: float,
    min_cells_detected: int,
    export_dir: Path,
) -> dict[str, Any]:
    """Export post-QC matrices, GTF PyRanges table, chromsizes & JASPAR cistromes."""
    export_dir.mkdir(parents=True, exist_ok=True)
    helper = export_dir / "_export_scplus.py"
    meta_json = export_dir / "_export_meta.json"
    helper.write_text(
        """import csv, json, sys
from pathlib import Path
import numpy as np
import scipy.io as sio
import scipy.sparse as sp

sys.path.insert(0, sys.argv[1])
from aria.utils.eregulon_core import (
    _read_fai_index,
    compute_aucell_scores,
    fetch_peak_center_sequences_from_indexed_fasta,
    parse_meme_motifs,
    resolve_motifs_to_expressed_tfs,
    scan_motifs_in_sequences,
)
from aria.utils.motifs import load_local_motif_collection
from aria.utils.mudata_io import read_h5mu_paired
from aria.utils.peak2gene_core import (
    canonical_contig,
    is_standard_contig,
    parse_gtf_tss_detailed,
    parse_peak_coord,
)

h5mu_p = sys.argv[2]
gtf_p = sys.argv[3]
fa_p = sys.argv[4]
motif_coll = sys.argv[5]
motif_pval = float(sys.argv[6])
min_cells = int(sys.argv[7])
out = Path(sys.argv[8])
meta_p = Path(sys.argv[9])

paired = read_h5mu_paired(h5mu_p, require_raw_counts=True)
assert paired["status"] == "success", paired
barcodes = list(paired["shared_barcodes"])
rna = paired["rna_adata"]
atac = paired["atac_adata"]
r_map = {str(b): i for i, b in enumerate(rna.obs_names)}
a_map = {str(b): i for i, b in enumerate(atac.obs_names)}
r_all = sp.csc_matrix(paired["rna_counts_matrix"][[r_map[b] for b in barcodes], :])
a_all = sp.csc_matrix(paired["atac_counts_matrix"][[a_map[b] for b in barcodes], :])

all_genes = [str(g) for g in rna.var_names]
gene_nnz = np.diff(r_all.indptr)
gene_records, _ = parse_gtf_tss_detailed(gtf_p, wanted=set(all_genes))

# Annotated standard-contig expressed genes (>= min_cells, non-zero variance)
keep_g_idx = []
keep_genes = []
annot_rows = []
for gi, g in enumerate(all_genes):
    if int(gene_nnz[gi]) < min_cells:
        continue
    rec = gene_records.get(g)
    if rec is None:
        continue
    c = rec.get("canon_chrom") or canonical_contig(rec["chrom"])
    if not is_standard_contig(c):
        continue
    col = r_all[:, gi].toarray().ravel()
    if float(np.std(col)) <= 1e-12:
        continue
    keep_g_idx.append(gi)
    keep_genes.append(g)
    width = max(1, int(rec["end"]) - int(rec["start"]))
    annot_rows.append({
        "Chromosome": c,
        "Start": int(rec["start"]),
        "End": int(rec["end"]),
        "Strand": str(rec["strand"]),
        "Gene": g,
        "Transcription_Start_Site": int(rec["tss"]),
        "Transcript_type": str(rec.get("gene_biotype") or "protein_coding"),
        "Gene_size_weight": 1.0,
    })

r_sub = r_all[:, keep_g_idx].tocsc()

# Standard-contig accessible peaks (>= min_cells) with valid 500bp FASTA sequences
all_peaks = [str(p) for p in atac.var_names]
peak_nnz = np.diff(a_all.indptr)
cand_p_idx = []
cand_p_raw = []
cand_p_canon = []
cand_coords = []
for pi, pname in enumerate(all_peaks):
    if int(peak_nnz[pi]) < min_cells:
        continue
    parsed = parse_peak_coord(pname)
    if parsed is None or not is_standard_contig(parsed[0]):
        continue
    c = canonical_contig(parsed[0]) or parsed[0]
    s, e = int(parsed[1]), int(parsed[2])
    cand_p_idx.append(pi)
    cand_p_raw.append(pname)
    cand_p_canon.append(f"{c}:{s}-{e}")
    cand_coords.append((c, s, e))

seq_bundle = fetch_peak_center_sequences_from_indexed_fasta(
    cand_p_raw, fa_p, scan_window_bp=500
)
valid_mask = seq_bundle["valid_mask"]
valid_pos = np.flatnonzero(valid_mask)

# Deduplicate canonical `chr:start-end` if any
seen_reg = set()
final_pos = []
for vp in valid_pos:
    reg_id = cand_p_canon[int(vp)]
    if reg_id not in seen_reg:
        seen_reg.add(reg_id)
        final_pos.append(int(vp))

keep_p_idx = [cand_p_idx[vp] for vp in final_pos]
keep_regions = [cand_p_canon[vp] for vp in final_pos]
keep_coords = [cand_coords[vp] for vp in final_pos]
keep_seqs = [str(seq_bundle["sequences"][vp]) for vp in final_pos]
a_sub = a_all[:, keep_p_idx].tocsc()

# Resolve motifs on the shared expressed gene universe & scan fixed 500bp windows
loaded = load_local_motif_collection(motif_coll)
assert loaded is not None, f"Missing motif collection: {motif_coll}"
meme_path, _ = loaded
motifs = parse_meme_motifs(meme_path)
sub_gene_nnz = np.diff(r_sub.indptr)
resolved_tfs, unmapped = resolve_motifs_to_expressed_tfs(
    motifs, keep_genes, sub_gene_nnz, min_cells_detected=min_cells
)
unique_m_idx = sorted({r.motif_idx for r in resolved_tfs})
m_to_row = {m_i: r_i for r_i, m_i in enumerate(unique_m_idx)}
hits_bool, scanner_backend = scan_motifs_in_sequences(
    meme_path=meme_path,
    motifs=motifs,
    sequences=keep_seqs,
    motif_Indices=unique_m_idx,
    pvalue=motif_pval,
    rc=True,
    backend="auto",
)

# Aggregate union cistrome per testable TF (skipping any symbol containing `_`)
tf_to_peaks = {}
for r in resolved_tfs:
    if "_" in r.tf_gene:
        continue
    row = hits_bool[m_to_row[r.motif_idx]]
    hit_idxs = np.flatnonzero(row)
    if hit_idxs.size == 0:
        continue
    tf_to_peaks.setdefault(r.tf_gene, set()).update(hit_idxs.tolist())

testable_tfs = sorted(tf_to_peaks.keys())

# Write exported matrices & tables
(out / "barcodes.tsv").write_text("\\n".join(barcodes) + "\\n", encoding="utf-8")
(out / "genes.tsv").write_text("\\n".join(keep_genes) + "\\n", encoding="utf-8")
(out / "regions.tsv").write_text("\\n".join(keep_regions) + "\\n", encoding="utf-8")
(out / "tfs.txt").write_text("\\n".join(testable_tfs) + "\\n", encoding="utf-8")
sio.mmwrite(str(out / "rna_counts.mtx"), r_sub)
sio.mmwrite(str(out / "atac_counts.mtx"), a_sub)

with open(out / "annot.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(
        fh,
        fieldnames=[
            "Chromosome", "Start", "End", "Strand", "Gene",
            "Transcription_Start_Site", "Transcript_type", "Gene_size_weight",
        ],
        delimiter="\\t",
    )
    w.writeheader()
    for row in annot_rows:
        w.writerow(row)

fai = _read_fai_index(Path(f"{fa_p}.fai"))
chrom_rows = []
seen_c = set()
for raw_c, (clen, _, _, _) in fai.items():
    cc = canonical_contig(raw_c) or raw_c
    if is_standard_contig(cc) and cc not in seen_c:
        seen_c.add(cc)
        chrom_rows.append({"Chromosome": cc, "Start": 0, "End": int(clen)})
with open(out / "chromsizes.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["Chromosome", "Start", "End"], delimiter="\\t")
    w.writeheader()
    for row in chrom_rows:
        w.writerow(row)

n_cistrome_rows = 0
with open(out / "cistromes.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(
        fh, fieldnames=["tf", "Chromosome", "Start", "End"], delimiter="\\t"
    )
    w.writeheader()
    for tf in testable_tfs:
        for p_i in sorted(tf_to_peaks[tf]):
            c, s, e = keep_coords[p_i]
            w.writerow({"tf": tf, "Chromosome": c, "Start": s, "End": e})
            n_cistrome_rows += 1

meta_p.write_text(json.dumps({
    "n_cells": len(barcodes),
    "n_genes": len(keep_genes),
    "n_regions": len(keep_regions),
    "n_testable_tfs": len(testable_tfs),
    "testable_tfs": testable_tfs,
    "genes": keep_genes,
    "n_cistrome_rows": n_cistrome_rows,
    "scanner_backend": scanner_backend,
}), encoding="utf-8")
""",
        encoding="utf-8",
    )
    proc = subprocess.run(
        [
            _resolve_conda_bin(),
            "run",
            "--no-capture-output",
            "-n",
            integration_env,
            "python",
            str(helper),
            str(ROOT),
            str(h5mu_path),
            str(gtf_path),
            str(fasta_path),
            str(motif_collection),
            str(motif_pvalue),
            str(min_cells_detected),
            str(export_dir),
            str(meta_json),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not meta_json.is_file():
        raise RuntimeError(
            f"SCENIC+ export helper failed ({proc.returncode}):\n"
            f"STDOUT:\n{proc.stdout[-2500:]}\nSTDERR:\n{proc.stderr[-2500:]}"
        )
    return json.loads(meta_json.read_text(encoding="utf-8"))


def _compute_aria_aucell_and_random_baseline(
    integration_env: str,
    export_dir: Path,
    aria_eregulons: list[dict[str, Any]],
    ref_aucell_npz: Path,
    work_dir: Path,
    seed: int = 0,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], list[float]]:
    """Compute ARIA RNA AUCell on identical `rna_counts.mtx` + random size-matched baseline."""
    helper = work_dir / "_aucell_eval.py"
    in_sets = work_dir / "_aria_sets.json"
    out_npz = work_dir / "_aria_aucell.npz"
    aria_sets = {
        f"{r['tf']}__{canonical_sign_class(r.get('sign_class'))}": [
            g.strip() for g in str(r.get("target_genes") or "").split(";") if g.strip()
        ]
        if isinstance(r.get("target_genes"), str)
        else list(r.get("target_genes") or [])
        for r in aria_eregulons
    }
    in_sets.write_text(json.dumps(aria_sets), encoding="utf-8")
    helper.write_text(
        """import json, sys
from pathlib import Path
import numpy as np
import scipy.io as sio
import scipy.sparse as sp

sys.path.insert(0, sys.argv[1])
from aria.utils.eregulon_core import compute_aucell_scores
from aria.benchmarks.concordance_eregulons import spearman_corr

exp_dir = Path(sys.argv[2])
aria_sets = json.loads(Path(sys.argv[3]).read_text(encoding="utf-8"))
ref_npz = np.load(sys.argv[4], allow_pickle=True)
out_p = sys.argv[5]
seed = int(sys.argv[6])

genes = [x.strip() for x in (exp_dir / "genes.tsv").read_text(encoding="utf-8").splitlines() if x.strip()]
rna_csc = sp.csc_matrix(sio.mmread(str(exp_dir / "rna_counts.mtx")))

a_keys = sorted(aria_sets.keys())
a_mat = compute_aucell_scores(rna_csc, genes, [aria_sets[k] for k in a_keys], top_frac=0.05)

r_keys = [str(k) for k in ref_npz["module_keys"].tolist()]
r_mat = np.asarray(ref_npz["aucell"], dtype=np.float32)
r_map = {k: r_mat[:, i] for i, k in enumerate(r_keys)}
a_map = {k: a_mat[:, i] for i, k in enumerate(a_keys)}

# Size-matched random gene-set AUCell Spearman baseline over shared modules
rng = np.random.default_rng(seed)
shared = sorted(set(a_map.keys()) & set(r_map.keys()))
rand_sets_a = []
rand_sets_b = []
for k in shared:
    sz_a = max(5, len(aria_sets.get(k, [])))
    for _ in range(5):
        pick_a = [genes[int(i)] for i in rng.choice(len(genes), size=min(sz_a, len(genes)), replace=False)]
        pick_b = [genes[int(i)] for i in rng.choice(len(genes), size=min(sz_a, len(genes)), replace=False)]
        rand_sets_a.append(pick_a)
        rand_sets_b.append(pick_b)

rand_corrs = []
if rand_sets_a:
    ra = compute_aucell_scores(rna_csc, genes, rand_sets_a, top_frac=0.05)
    rb = compute_aucell_scores(rna_csc, genes, rand_sets_b, top_frac=0.05)
    for j in range(ra.shape[1]):
        rho = spearman_corr(ra[:, j], rb[:, j])
        if rho is not None:
            rand_corrs.append(float(rho))

np.savez_compressed(
    out_p,
    aria_keys=np.asarray(a_keys, dtype=object),
    aria_aucell=a_mat,
    rand_corrs=np.asarray(rand_corrs, dtype=np.float64),
)
""",
        encoding="utf-8",
    )
    proc = subprocess.run(
        [
            _resolve_conda_bin(),
            "run",
            "--no-capture-output",
            "-n",
            integration_env,
            "python",
            str(helper),
            str(ROOT),
            str(export_dir),
            str(in_sets),
            str(ref_aucell_npz),
            str(out_npz),
            str(seed),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not out_npz.is_file():
        raise RuntimeError(f"AUCell baseline helper failed: {proc.stderr}")

    aria_z = np.load(str(out_npz), allow_pickle=True)
    ref_z = np.load(str(ref_aucell_npz), allow_pickle=True)
    a_keys = [str(k) for k in aria_z["aria_keys"].tolist()]
    a_mat = np.asarray(aria_z["aria_aucell"], dtype=np.float32)
    r_keys = [str(k) for k in ref_z["module_keys"].tolist()]
    r_mat = np.asarray(ref_z["aucell"], dtype=np.float32)
    rand_corrs = [float(x) for x in aria_z["rand_corrs"].tolist()]
    return (
        {k: a_mat[:, i] for i, k in enumerate(a_keys)},
        {k: r_mat[:, i] for i, k in enumerate(r_keys)},
        rand_corrs,
    )


def _run_scenicplus_seed(
    scenicplus_env: str,
    export_dir: Path,
    work_dir: Path,
    seed: int,
    ray_cpus: int,
) -> tuple[dict[str, Any], list[dict[str, Any]], Path]:
    out_trip = work_dir / f"scplus_seed{seed}_triplets.tsv"
    out_npz = work_dir / f"scplus_seed{seed}_aucell.npz"
    out_json = work_dir / f"scplus_seed{seed}_summary.json"
    proc = subprocess.run(
        [
            _resolve_conda_bin(),
            "run",
            "--no-capture-output",
            "-n",
            scenicplus_env,
            "python",
            str(ROOT / "aria" / "scripts" / "benchmark_eregulons_scenicplus_ref.py"),
            "--in-dir",
            str(export_dir),
            "--seed",
            str(seed),
            "--ray-cpus",
            str(ray_cpus),
            "--distance-bp",
            "500000",
            "--min-target-genes",
            "5",
            "--out-triplets",
            str(out_trip),
            "--out-aucell-npz",
            str(out_npz),
            "--out-json",
            str(out_json),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not out_json.is_file():
        raise RuntimeError(
            f"SCENIC+ seed={seed} failed ({proc.returncode}):\n"
            f"STDOUT:\n{proc.stdout[-2500:]}\nSTDERR:\n{proc.stderr[-2500:]}"
        )
    summary = json.loads(out_json.read_text(encoding="utf-8"))
    rows = _read_csv_rows(out_trip, delimiter="\t")
    return summary, rows, out_npz


def _sanitize_public_obj(obj: Any) -> Any:
    """Strip any machine-absolute filesystem paths from public benchmark JSON."""
    if isinstance(obj, dict):
        out: dict[str, Any] = {}
        for k, v in obj.items():
            if k in (
                "conda_prefix",
                "output_path",
                "triplets_path",
                "h5mu_output_path",
                "input_path",
            ):
                continue
            out[str(k)] = _sanitize_public_obj(v)
        return out
    if isinstance(obj, list):
        return [_sanitize_public_obj(x) for x in obj]
    if isinstance(obj, str):
        if obj.startswith(_FORBIDDEN_PATH_PREFIXES):
            return Path(obj).name
        return obj
    return obj


def main(argv: list[str] | None = None) -> int:
    home = Path.home()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--h5mu",
        default=str(home / "Samael/Oligos/data/interim/paired/hc11_paired.h5mu"),
    )
    ap.add_argument(
        "--reference-h5ad",
        default=str(
            home / "Samael/Oligos/data/raw/reference/GSE278576_hippocampus_RNA.h5ad"
        ),
    )
    ap.add_argument("--sample-id", default="hc11")
    ap.add_argument(
        "--gtf",
        default=str(
            home / "Samael/H9-RNA/reference/gencode.v49.primary_assembly.annotation.gtf"
        ),
    )
    ap.add_argument(
        "--genome-fasta",
        default=str(home / ".aria/genomes/hg38/GRCh38.primary_assembly.genome.fa"),
    )
    ap.add_argument("--motif-collection", default="JASPAR2024_CORE_vertebrates")
    ap.add_argument("--chromatin-env", default="aria-chromatin-env")
    ap.add_argument("--integration-env", default="aria-integration-env")
    ap.add_argument("--scenicplus-env", default="scenicplus_env")
    ap.add_argument("--ray-cpus", type=int, default=24)
    ap.add_argument(
        "--out-json",
        default="docs/benchmark_results/multiome_concordance/v48_i6_eregulon_concordance_hc11.json",
    )
    args = ap.parse_args(argv)

    with tempfile.TemporaryDirectory(prefix="v48_i6_ereg_bench_") as td:
        work = Path(td)
        print("[1/5] Running Slice I1 integration_paired_qc...", flush=True)
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

        print(
            "[2/5] Running ARIA integration_eregulons (5e-5 seed=0 x2, seed=1, and 1e-4)...",
            flush=True,
        )
        base_params = {
            "h5mu_path": str(post_qc_h5mu),
            "gtf_path": str(args.gtf),
            "genome_fasta": str(args.genome_fasta),
            "motif_collection": str(args.motif_collection),
            "distance_kb": 500,
            "p2g_min_corr": 0.05,
            "p2g_max_padj": 0.05,
            "min_tf_gene_corr": 0.05,
            "min_target_genes": 5,
            "max_regulon_padj": 0.05,
            "standard_chroms_only": True,
        }

        t0 = time.time()
        aria_s0_a = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_eregulons.py",
            {
                **base_params,
                "motif_pvalue": 5e-5,
                "seed": 0,
                "output_dir": str(work / "aria_5e5_s0_a"),
            },
            work,
            "aria_5e5_s0_a",
        )
        aria_s0_runtime = round(time.time() - t0, 2)
        aria_s0_b = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_eregulons.py",
            {
                **base_params,
                "motif_pvalue": 5e-5,
                "seed": 0,
                "output_dir": str(work / "aria_5e5_s0_b"),
            },
            work,
            "aria_5e5_s0_b",
        )
        aria_s1 = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_eregulons.py",
            {
                **base_params,
                "motif_pvalue": 5e-5,
                "seed": 1,
                "output_dir": str(work / "aria_5e5_s1"),
            },
            work,
            "aria_5e5_s1",
        )
        aria_1e4 = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_eregulons.py",
            {
                **base_params,
                "motif_pvalue": 1e-4,
                "seed": 0,
                "output_dir": str(work / "aria_1e4_s0"),
            },
            work,
            "aria_1e4_s0",
        )

        ereg_s0_a = _read_csv_rows(Path(aria_s0_a["output_path"]))
        trip_s0_a = _read_csv_rows(Path(aria_s0_a["triplets_path"]))
        ereg_s0_b = _read_csv_rows(Path(aria_s0_b["output_path"]))
        trip_s0_b = _read_csv_rows(Path(aria_s0_b["triplets_path"]))
        ereg_s1 = _read_csv_rows(Path(aria_s1["output_path"]))
        trip_s1 = _read_csv_rows(Path(aria_s1["triplets_path"]))
        ereg_1e4 = _read_csv_rows(Path(aria_1e4["output_path"]))
        trip_1e4 = _read_csv_rows(Path(aria_1e4["triplets_path"]))

        det_same_seed = score_eregulon_determinism(
            ereg_s0_a, trip_s0_a, ereg_s0_b, trip_s0_b
        )
        det_cross_seed = score_eregulon_determinism(
            ereg_s0_a, trip_s0_a, ereg_s1, trip_s1
        )
        sens_5e5_vs_1e4 = score_eregulon_determinism(
            ereg_s0_a, trip_s0_a, ereg_1e4, trip_1e4
        )
        print(
            f"      ARIA 5e-5: {aria_s0_a['n_eregulons']} eRegulons, "
            f"{aria_s0_a['n_triplets']} triplets in {aria_s0_runtime}s "
            f"(bit_exact={det_same_seed['bit_exact']}, "
            f"seed0_vs_1_mod_jacc={det_cross_seed['module_jaccard']})",
            flush=True,
        )

        print(
            "[3/5] Exporting shared post-QC universe & JASPAR2024 5e-5 cistromes for SCENIC+...",
            flush=True,
        )
        export_dir = work / "scplus_export"
        export_meta = _export_scenicplus_bundle(
            args.integration_env,
            post_qc_h5mu,
            Path(args.gtf),
            Path(args.genome_fasta),
            str(args.motif_collection),
            5e-5,
            10,
            export_dir,
        )
        print(
            f"      Exported {export_meta['n_cells']} cells, {export_meta['n_genes']} genes, "
            f"{export_meta['n_regions']} regions, {export_meta['n_testable_tfs']} TFs",
            flush=True,
        )

        print(
            f"[4/5] Running external SCENIC+ 1.0.0 (seed=0 & seed=1, ray_cpus={args.ray_cpus})...",
            flush=True,
        )
        scplus_s0_sum, scplus_s0_trips, scplus_s0_npz = _run_scenicplus_seed(
            args.scenicplus_env, export_dir, work, seed=0, ray_cpus=args.ray_cpus
        )
        print(
            f"      SCENIC+ seed=0: {scplus_s0_sum['n_eregulons_post_selfloop_filter']} eRegulons, "
            f"{scplus_s0_sum['n_triplets']} triplets ({scplus_s0_sum['timings_sec']['total']}s)",
            flush=True,
        )
        scplus_s1_sum, scplus_s1_trips, scplus_s1_npz = _run_scenicplus_seed(
            args.scenicplus_env, export_dir, work, seed=1, ray_cpus=args.ray_cpus
        )
        print(
            f"      SCENIC+ seed=1: {scplus_s1_sum['n_eregulons_post_selfloop_filter']} eRegulons, "
            f"{scplus_s1_sum['n_triplets']} triplets ({scplus_s1_sum['timings_sec']['total']}s)",
            flush=True,
        )

        print("[5/5] Scoring 4-tier concordance & SCENIC+ multi-seed ceiling...", flush=True)
        shared_tfs = export_meta["testable_tfs"]
        shared_genes = export_meta["genes"]

        aria_auc_map, ref_s0_auc_map, rand_corrs = (
            _compute_aria_aucell_and_random_baseline(
                args.integration_env,
                export_dir,
                ereg_s0_a,
                scplus_s0_npz,
                work,
                seed=0,
            )
        )
        ref_s1_z = np.load(str(scplus_s1_npz), allow_pickle=True)
        ref_s1_keys = [str(k) for k in ref_s1_z["module_keys"].tolist()]
        ref_s1_mat = np.asarray(ref_s1_z["aucell"], dtype=np.float32)
        ref_s1_auc_map = {k: ref_s1_mat[:, i] for i, k in enumerate(ref_s1_keys)}

        conc_vs_s0 = score_eregulon_concordance(
            trip_s0_a,
            scplus_s0_trips,
            shared_tf_universe=shared_tfs,
            shared_gene_universe=shared_genes,
            aria_cell_aucell=aria_auc_map,
            ref_cell_aucell=ref_s0_auc_map,
            random_aucell_baseline_spearmans=rand_corrs,
        )
        conc_vs_s1 = score_eregulon_concordance(
            trip_s0_a,
            scplus_s1_trips,
            shared_tf_universe=shared_tfs,
            shared_gene_universe=shared_genes,
            aria_cell_aucell=aria_auc_map,
            ref_cell_aucell=ref_s1_auc_map,
            random_aucell_baseline_spearmans=rand_corrs,
        )
        scplus_ceiling = score_eregulon_concordance(
            scplus_s0_trips,
            scplus_s1_trips,
            shared_tf_universe=shared_tfs,
            shared_gene_universe=shared_genes,
            aria_cell_aucell=ref_s0_auc_map,
            ref_cell_aucell=ref_s1_auc_map,
            random_aucell_baseline_spearmans=rand_corrs,
        )

        from aria.version import __version__, collect_version_metadata

        manifest = _sanitize_public_obj(
            {
                "benchmark": "v48_i6_paired_eregulon_concordance",
                "slice": "V48_I6_candidate_eregulons",
                "aria_version": __version__,
                "provenance": collect_version_metadata(),
                "dataset": {
                    "accession": "GSE278576",
                    "sample_id": args.sample_id,
                    "tissue": "human_hippocampus_multiome_10x",
                    "n_cells_raw": qc_res["n_cells_before"],
                    "n_cells_post_qc": qc_res["n_cells_after"],
                    "n_shared_annotated_expressed_genes": export_meta["n_genes"],
                    "n_shared_standard_peaks": export_meta["n_regions"],
                    "n_shared_testable_tfs": export_meta["n_testable_tfs"],
                },
                "design": {
                    "aria_script": "aria/scripts/integration_eregulons.py",
                    "shared_core": "aria/utils/eregulon_core.py",
                    "reused_p2g_core": "aria/utils/peak2gene_core.py",
                    "external_reference": (
                        f"scenicplus=={scplus_s0_sum['scenicplus_version']} "
                        "(GBM r2g + arboreto GBM tf2g + GSEA eRegulons)"
                    ),
                    "motif_collection": aria_s0_a["motif_source"],
                    "scanner_backend": aria_s0_a["scanner_backend"],
                    "scan_window_bp": aria_s0_a["scan_window_bp"],
                    "motif_pvalue_primary": 5e-5,
                    "motif_pvalue_sensitivity": 1e-4,
                    "distance_bp": 500_000,
                    "p2g_min_corr": 0.05,
                    "p2g_max_padj": 0.05,
                    "min_tf_gene_corr": 0.05,
                    "min_target_genes": 5,
                    "max_regulon_padj": 0.05,
                    "regulon_padj_scope": "module_level",
                },
                "aria_summary": {
                    "runtime_sec": aria_s0_runtime,
                    "p2g_summary": aria_s0_a["p2g_summary"],
                    "n_motif_tf_mappings": aria_s0_a["n_motif_tf_mappings"],
                    "n_unmapped_motifs": aria_s0_a["n_unmapped_motifs"],
                    "unmapped_motifs_sample": aria_s0_a["unmapped_motifs_sample"],
                    "n_candidate_modules_tested": aria_s0_a["n_candidate_modules_tested"],
                    "n_eregulons": aria_s0_a["n_eregulons"],
                    "n_triplets": aria_s0_a["n_triplets"],
                    "n_unique_tfs": aria_s0_a["n_unique_tfs"],
                    "n_unique_enhancers": aria_s0_a["n_unique_enhancers"],
                    "n_unique_target_genes": aria_s0_a["n_unique_target_genes"],
                    "n_pos_tf_pos_p2g": aria_s0_a["n_pos_tf_pos_p2g"],
                    "n_neg_tf_pos_p2g": aria_s0_a["n_neg_tf_pos_p2g"],
                    "n_paralog_ambiguous": aria_s0_a["n_paralog_ambiguous"],
                    "n_depth_confounded": aria_s0_a["n_depth_confounded"],
                    "top_eregulons": aria_s0_a["top_eregulons"][:15],
                },
                "determinism_and_sensitivity": {
                    "aria_seed0_repeat_determinism": det_same_seed,
                    "aria_seed0_vs_seed1_null_stability": det_cross_seed,
                    "aria_motif_pvalue_5e5_vs_1e4_sensitivity": sens_5e5_vs_1e4,
                },
                "scenicplus_reference": {
                    "seed_0": scplus_s0_sum,
                    "seed_1": scplus_s1_sum,
                    "multi_seed_stochastic_ceiling_seed0_vs_seed1": scplus_ceiling,
                },
                "external_concordance_aria_vs_scenicplus": {
                    "vs_seed_0": conc_vs_s0,
                    "vs_seed_1": conc_vs_s1,
                },
                "governance": {
                    "associative_only": True,
                    "perturbation_validated": False,
                    "chromatin_conformation_validated": False,
                    "caveat": aria_s0_a["caveat"],
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
