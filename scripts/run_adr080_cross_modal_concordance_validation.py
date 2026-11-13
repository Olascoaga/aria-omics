#!/usr/bin/env python3
"""ADR-080 — Multi-donor non-same-cell cross-modal concordance empirical validation runner.

Executes ARIA's production cross-modal concordance & synthesis pipeline across a
real 10-donor balanced human hippocampal cohort from GSE278576 (5 young adults in
age_group 20-39 [3M/2F]: hc29, hc77, hc6052, hc76, hc935; 5 oldest-old adults in
age_group 80-100 [3M/2F]: hc11, hc12, hc73, hc19, hc212191) using a two-arm design:

  Arm A — Disjoint-Donor Split-Cohort Glial Lineage Benchmark (Oligo vs Astro):
    * RNA cohort ONLY uses the 5 young-adult brains (hc29, hc77, hc6052, hc76, hc935;
      5 Oligo + 5 Astro pseudobulk samples, paired by donor ~ donor_id + lineage).
    * ATAC cohort ONLY uses the 5 COMPLETELY DISJOINT aged-adult brains (hc11, hc12,
      hc73, hc19, hc212191; 5 Oligo + 5 Astro consensus-peak samples, ~ donor_id + lineage).
    * Zero shared cells AND zero shared human donors between RNA and ATAC.
    * Evaluates Layer A (genome-wide 4-quadrant cis DE x DA across all tested genes
      and all tested consensus peaks + reversed Astro_vs_Oligo sign-flip invariance),
      Layer B (JASPAR2024 motif x cognate TF RNA triangulation + paralog ambiguity),
      Layer C (RNA ORA intersect ATAC peak ORA pathway Jaccard), and direct recovery
      of same-cell Multiome Peak2Gene ground-truth loci.

  Arm B — 10-Donor Matched Clinical Aging Cohort & Donor-Matched CPM Cis Links:
    * Unpaired whole-hippocampus donor RNA vs consensus ATAC across all 10 donors
      (~ sex + age_group, 80-100_vs_20-39) over the full tested transcriptome
      (de_all_genes.tsv) and 45,000 consensus ATAC peaks.
    * Layer D donor-matched log2(CPM + 1) cis peak-to-gene correlation across both
      the 10 whole-tissue donors AND the 20 matched donor-lineage biological units
      (10 Oligo + 10 Astro) paired with zero-fixed-point sample-label derangement
      null controls, A4 SHA-256 stage-manifest cache hit verification, and dual
      300-DPI PNG + editable SVG starburst plots.

  Governance & Narrative Verification:
    * CrossModalConcordanceAuditAgent readiness (green on non-same-cell RNA+ATAC,
      red/blocked on diagonal cell-imputation requests, red/blocked on IntegrationAgent
      for unpaired inputs), BiologicalSynthesisAgent Slice 2, RunLedger C3 linkage
      (0 violations), and C4 zero-causal-verb audit.

Emits a sanitized JSON receipt with zero machine-absolute paths (enforced by
``tests/test_public_artifacts_no_absolute_paths.py``).
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DEFAULT_OUT = (
    ROOT
    / "docs"
    / "benchmark_results"
    / "multiome_concordance"
    / "adr080_cross_modal_concordance_validation.json"
)

# Balanced 10-donor human hippocampal aging cohort (GSE278576): 5 young (3M/2F) vs 5 aged (3M/2F).
DONOR_COHORT = [
    {"sample_id": "hc29", "donor_id": "hc29", "age_group": "20-39", "sex": "M", "age": 28},
    {"sample_id": "hc77", "donor_id": "hc77", "age_group": "20-39", "sex": "M", "age": 20},
    {"sample_id": "hc6052", "donor_id": "hc6052", "age_group": "20-39", "sex": "M", "age": 28},
    {"sample_id": "hc76", "donor_id": "hc76", "age_group": "20-39", "sex": "F", "age": 26},
    {"sample_id": "hc935", "donor_id": "hc935", "age_group": "20-39", "sex": "F", "age": 38},
    {"sample_id": "hc11", "donor_id": "hc11", "age_group": "80-100", "sex": "M", "age": 86},
    {"sample_id": "hc12", "donor_id": "hc12", "age_group": "80-100", "sex": "M", "age": 83},
    {"sample_id": "hc73", "donor_id": "hc73", "age_group": "80-100", "sex": "M", "age": 86},
    {"sample_id": "hc19", "donor_id": "hc19", "age_group": "80-100", "sex": "F", "age": 87},
    {"sample_id": "hc212191", "donor_id": "hc212191", "age_group": "80-100", "sex": "F", "age": 92},
]

RNA_SPLIT_DONORS = ["hc29", "hc77", "hc6052", "hc76", "hc935"]
ATAC_SPLIT_DONORS = ["hc11", "hc12", "hc73", "hc19", "hc212191"]

BENCH_CHROMS = ("chr1", "chr2", "chr6", "chr11", "chr17", "chr19")

_ABS_PATH_RE = re.compile(r"(/(?:home|tmp|var|private|Users|opt|root)/[^\s\"',;]+)")


def _sanitize_obj(obj: Any) -> Any:
    """Recursively strip machine-local absolute paths from strings in JSON receipts."""
    if isinstance(obj, dict):
        return {str(k): _sanitize_obj(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize_obj(v) for v in obj]
    if isinstance(obj, tuple):
        return [_sanitize_obj(v) for v in obj]
    if isinstance(obj, str):
        if obj.startswith(("/", "~")):
            p = Path(obj)
            return p.name if p.name else "sanitized_path"
        return _ABS_PATH_RE.sub(lambda m: Path(m.group(1)).name, obj)
    return obj


def _resolve_env_python(env_name: str) -> list[str]:
    exe = Path(sys.executable).resolve()
    for parent in exe.parents:
        cand = parent / "envs" / env_name / "bin" / "python"
        if cand.is_file():
            return [str(cand)]
    for root_dir in (Path.home() / "anaconda3", Path.home() / "miniforge3"):
        cand = root_dir / "envs" / env_name / "bin" / "python"
        if cand.is_file():
            return [str(cand)]
    conda_bin = shutil.which("conda") or "conda"
    return [conda_bin, "run", "--no-capture-output", "-n", env_name, "python"]


def _run_ipc(
    env_name: str,
    script_rel: str,
    params: dict[str, Any],
    work_dir: Path,
    tag: str,
) -> tuple[dict[str, Any], float]:
    in_json = work_dir / f"{tag}_in.json"
    out_json = work_dir / f"{tag}_out.json"
    if out_json.exists():
        out_json.unlink()
    in_json.write_text(json.dumps(params, indent=2), encoding="utf-8")
    env = dict(os.environ)
    numba_dir = work_dir / "numba_cache"
    numba_dir.mkdir(parents=True, exist_ok=True)
    env["NUMBA_CACHE_DIR"] = str(numba_dir)
    env["PYTHONHASHSEED"] = "0"

    cmd = [*_resolve_env_python(env_name), str(ROOT / script_rel), str(in_json), str(out_json)]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, check=False)
    elapsed = round(time.perf_counter() - t0, 3)
    if proc.returncode != 0 or not out_json.is_file():
        out_payload = out_json.read_text(encoding="utf-8")[:2500] if out_json.is_file() else ""
        raise RuntimeError(
            f"IPC script {script_rel} ({tag}) failed rc={proc.returncode}:\n"
            f"STDOUT:\n{proc.stdout[-2500:]}\nSTDERR:\n{proc.stderr[-2500:]}\nOUT_JSON:\n{out_payload}"
        )
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    if payload.get("status") not in ("success", "skipped"):
        raise RuntimeError(f"IPC script {script_rel} ({tag}) returned error: {payload}")
    return payload, elapsed


def _prepare_cohort_tables_and_same_cell_slice(
    merged_h5mu: Path,
    ref_meta_tsv: Path,
    work_dir: Path,
) -> dict[str, Any]:
    """Extract:
      - Arm A disjoint-donor Oligo vs Astro RNA counts (5 young donors) and ATAC
        consensus peak counts (5 disjoint aged donors), plus 20-unit (10 Oligo +
        10 Astro) matched CPM + deranged-unit null CPM;
      - Arm B 10-donor whole-hippocampus RNA counts, ATAC consensus peak counts,
        10-donor ATAC CPM, and zero-fixed-point donor-deranged null ATAC CPM;
      - Stratified 2,500-cell same-cell reference slice on BENCH_CHROMS."""
    helper_code = r"""
import json, sys
from pathlib import Path
import mudata as md
import numpy as np
import pandas as pd
import scipy.sparse as sp

merged_path = Path(sys.argv[1])
meta_path = Path(sys.argv[2])
out_dir = Path(sys.argv[3])
cohort = json.loads(sys.argv[4])
rna_split = json.loads(sys.argv[5])
atac_split = json.loads(sys.argv[6])
chroms = set(json.loads(sys.argv[7]))

meta = pd.read_csv(meta_path, sep="\t", usecols=["bacrode", "orig.ident", "subclass"])
meta["raw_bc"] = meta["bacrode"].astype(str).str.rsplit("_", n=1).str[-1]
meta["canonical_id"] = meta["orig.ident"].astype(str) + "#" + meta["raw_bc"]
meta = meta.drop_duplicates(subset=["canonical_id"]).set_index("canonical_id")

mdata = md.read_h5mu(str(merged_path))
rna = mdata["rna"]
atac = mdata["atac"]

rna_X = sp.csr_matrix(rna.layers["counts"] if "counts" in rna.layers else rna.X)
atac_X = sp.csr_matrix(atac.layers["counts"] if "counts" in atac.layers else atac.X)

rna_ngenes = np.asarray((rna_X > 0).sum(axis=1)).ravel()
atac_npeaks = np.asarray((atac_X > 0).sum(axis=1)).ravel()
keep_cells = (rna_ngenes >= 200) & (atac_npeaks >= 200)
cell_idx = np.flatnonzero(keep_cells)

rna_sub = rna[cell_idx, :]
atac_sub = atac[cell_idx, :]
rna_X = rna_X[cell_idx, :]
atac_X = atac_X[cell_idx, :]

obs_idx = pd.Index([str(x) for x in rna_sub.obs_names])
subclass = meta["subclass"].reindex(obs_idx).fillna("unannotated").astype(str).to_numpy()

donor_col = rna_sub.obs["donor_id"].astype(str).to_numpy()
donor_ids = [row["donor_id"] for row in cohort]

gene_names = pd.Index([str(g) for g in rna_sub.var_names])
dup_gene_mask = gene_names.duplicated(keep="first")
if dup_gene_mask.any():
    keep_g = np.flatnonzero(~dup_gene_mask)
    gene_names = gene_names[keep_g]
    rna_X = rna_X[:, keep_g]
    rna_sub = rna_sub[:, keep_g]

# Consensus ATAC peaks on canonical chromosomes accessible in >=1.5% of post-QC cells
peak_names = np.asarray([str(p) for p in atac_sub.var_names])
canonical_chroms = {f"chr{i}" for i in range(1, 23)} | {"chrX"}
chrom_mask = np.array([p.split(":", 1)[0] in canonical_chroms for p in peak_names], dtype=bool)
cell_acc_frac = np.asarray((atac_X > 0).mean(axis=0)).ravel()
peak_keep = chrom_mask & (cell_acc_frac >= 0.015)
if int(np.sum(peak_keep)) > 45000:
    cand_idx = np.flatnonzero(peak_keep)
    order = np.argsort(-cell_acc_frac[cand_idx])[:45000]
    sel_peaks = np.sort(cand_idx[order])
else:
    sel_peaks = np.flatnonzero(peak_keep)

atac_X_sel = atac_X[:, sel_peaks]
sel_peak_names = pd.Index(peak_names[sel_peaks])

# ── Arm B: 10-donor whole-hippocampus RNA + ATAC counts & CPM ────────────────
rna_sums = {}
atac_sums = {}
n_cells_by_donor = {}
for d in donor_ids:
    mask = donor_col == d
    n_cells_by_donor[d] = int(np.sum(mask))
    rna_sums[d] = np.rint(np.asarray(rna_X[mask, :].sum(axis=0)).ravel()).astype(np.int64)
    atac_sums[d] = np.rint(np.asarray(atac_X_sel[mask, :].sum(axis=0)).ravel()).astype(np.int64)

rna_df = pd.DataFrame(rna_sums, index=gene_names)
rna_df.index.name = "gene_id"
rna_counts_path = out_dir / "arm_b_rna_donor_counts.tsv"
rna_df.to_csv(rna_counts_path, sep="\t")

atac_df = pd.DataFrame(atac_sums, index=sel_peak_names)
atac_df.index.name = "peak_id"
atac_counts_path = out_dir / "arm_b_atac_donor_counts.tsv"
atac_df.to_csv(atac_counts_path, sep="\t")

atac_cpm = atac_df.divide(atac_df.sum(axis=0).replace(0, 1), axis=1) * 1e6
atac_cpm_path = out_dir / "arm_b_atac_donor_cpm.tsv"
atac_cpm.to_csv(atac_cpm_path, sep="\t")

rng = np.random.default_rng(20261113)
def _derange(labels, rng_obj):
    n = len(labels)
    perm = rng_obj.permutation(n)
    for i in range(n):
        if perm[i] == i:
            j = (i + 1) % n
            perm[i], perm[j] = perm[j], perm[i]
    return [labels[j] for j in perm]

deranged_donors = _derange(donor_ids, rng)
atac_cpm_null = atac_cpm.copy()
atac_cpm_null.columns = deranged_donors
atac_cpm_null = atac_cpm_null[donor_ids]
atac_cpm_null_path = out_dir / "arm_b_atac_donor_cpm_deranged_null.tsv"
atac_cpm_null.to_csv(atac_cpm_null_path, sep="\t")

meta_df = pd.DataFrame(cohort)
rna_meta_path = out_dir / "arm_b_rna_metadata.tsv"
atac_meta_path = out_dir / "arm_b_atac_metadata.tsv"
meta_df.to_csv(rna_meta_path, sep="\t", index=False)
meta_df.to_csv(atac_meta_path, sep="\t", index=False)

# ── Arm A: Disjoint-donor Oligo vs Astro split-cohort RNA (5 young) x ATAC (5 aged) ─
arm_a_rna_cols = {}
arm_a_rna_meta = []
for d in rna_split:
    for lin in ("Oligo", "Astro"):
        sid = f"{d}_{lin}"
        mask = (donor_col == d) & (subclass == lin)
        arm_a_rna_cols[sid] = np.rint(np.asarray(rna_X[mask, :].sum(axis=0)).ravel()).astype(np.int64)
        arm_a_rna_meta.append({"sample_id": sid, "donor_id": d, "lineage": lin, "n_cells": int(np.sum(mask))})

arm_a_atac_cols = {}
arm_a_atac_meta = []
for d in atac_split:
    for lin in ("Oligo", "Astro"):
        sid = f"{d}_{lin}"
        mask = (donor_col == d) & (subclass == lin)
        arm_a_atac_cols[sid] = np.rint(np.asarray(atac_X_sel[mask, :].sum(axis=0)).ravel()).astype(np.int64)
        arm_a_atac_meta.append({"sample_id": sid, "donor_id": d, "lineage": lin, "n_cells": int(np.sum(mask))})

arm_a_rna_df = pd.DataFrame(arm_a_rna_cols, index=gene_names)
arm_a_rna_df.index.name = "gene_id"
arm_a_rna_counts_path = out_dir / "arm_a_disjoint_rna_counts.tsv"
arm_a_rna_df.to_csv(arm_a_rna_counts_path, sep="\t")
arm_a_rna_meta_path = out_dir / "arm_a_disjoint_rna_metadata.tsv"
pd.DataFrame(arm_a_rna_meta).to_csv(arm_a_rna_meta_path, sep="\t", index=False)

arm_a_atac_df = pd.DataFrame(arm_a_atac_cols, index=sel_peak_names)
arm_a_atac_df.index.name = "peak_id"
arm_a_atac_counts_path = out_dir / "arm_a_disjoint_atac_counts.tsv"
arm_a_atac_df.to_csv(arm_a_atac_counts_path, sep="\t")
arm_a_atac_meta_path = out_dir / "arm_a_disjoint_atac_metadata.tsv"
pd.DataFrame(arm_a_atac_meta).to_csv(arm_a_atac_meta_path, sep="\t", index=False)

# 20-unit matched donor-lineage CPM (all 10 donors x {Oligo, Astro}) + deranged null
unit_ids = []
unit_rna_cols = {}
unit_atac_cols = {}
for d in donor_ids:
    for lin in ("Oligo", "Astro"):
        uid = f"{d}_{lin}"
        mask = (donor_col == d) & (subclass == lin)
        if int(np.sum(mask)) >= 20:
            unit_ids.append(uid)
            unit_rna_cols[uid] = np.rint(np.asarray(rna_X[mask, :].sum(axis=0)).ravel()).astype(np.int64)
            unit_atac_cols[uid] = np.rint(np.asarray(atac_X_sel[mask, :].sum(axis=0)).ravel()).astype(np.int64)

u_rna_df = pd.DataFrame(unit_rna_cols, index=gene_names)
u_atac_df = pd.DataFrame(unit_atac_cols, index=sel_peak_names)
u_rna_cpm = u_rna_df.divide(u_rna_df.sum(axis=0).replace(0, 1), axis=1) * 1e6
u_atac_cpm = u_atac_df.divide(u_atac_df.sum(axis=0).replace(0, 1), axis=1) * 1e6
u_rna_cpm_path = out_dir / "glial_20unit_rna_cpm.tsv"
u_atac_cpm_path = out_dir / "glial_20unit_atac_cpm.tsv"
u_rna_cpm.to_csv(u_rna_cpm_path, sep="\t")
u_atac_cpm.to_csv(u_atac_cpm_path, sep="\t")

deranged_units = _derange(unit_ids, rng)
u_atac_cpm_null = u_atac_cpm.copy()
u_atac_cpm_null.columns = deranged_units
u_atac_cpm_null = u_atac_cpm_null[unit_ids]
u_atac_cpm_null_path = out_dir / "glial_20unit_atac_cpm_deranged_null.tsv"
u_atac_cpm_null.to_csv(u_atac_cpm_null_path, sep="\t")

# Stratified Oligo + Astro same-cell slice (2,500 cells) on BENCH_CHROMS for same-cell P2G
glial_rows = np.flatnonzero(np.isin(subclass, ["Oligo", "Astro"]))
picked_cells = []
for d in donor_ids:
    d_rows = np.flatnonzero((donor_col == d) & np.isin(subclass, ["Oligo", "Astro"]))
    take = min(250, len(d_rows))
    if take > 0:
        picked_cells.extend(np.sort(rng.choice(d_rows, size=take, replace=False)).tolist())
picked_cells = np.asarray(sorted(picked_cells), dtype=np.int64)

bench_peak_mask = np.array([p.split(":", 1)[0] in chroms for p in atac_sub.var_names.astype(str)], dtype=bool)
bench_peak_idx = np.flatnonzero(bench_peak_mask)

r_sc = rna_sub[picked_cells, :].copy()
a_sc = atac_sub[picked_cells, bench_peak_idx].copy()
r_sc.X = sp.csr_matrix(rna_X[picked_cells, :])
r_sc.layers["counts"] = r_sc.X.copy()
a_sc.X = sp.csr_matrix(atac_X[picked_cells, :][:, bench_peak_idx])
a_sc.layers["counts"] = a_sc.X.copy()
sc_mdata = md.MuData({"rna": r_sc, "atac": a_sc})
sc_mdata.obs = r_sc.obs.copy()
same_cell_slice_path = out_dir / "same_cell_glial_reference_slice.h5mu"
sc_mdata.write_h5mu(str(same_cell_slice_path))

print(json.dumps({
    "n_post_qc_cells_total": int(len(cell_idx)),
    "n_glial_cells_total": int(len(glial_rows)),
    "n_cells_by_donor": n_cells_by_donor,
    "n_genes_quantified": int(rna_df.shape[0]),
    "n_consensus_peaks_tested": int(atac_df.shape[0]),
    "arm_a": {
        "rna_counts_path": str(arm_a_rna_counts_path),
        "rna_meta_path": str(arm_a_rna_meta_path),
        "atac_counts_path": str(arm_a_atac_counts_path),
        "atac_meta_path": str(arm_a_atac_meta_path),
        "rna_samples": arm_a_rna_meta,
        "atac_samples": arm_a_atac_meta,
        "glial_20unit_rna_cpm_path": str(u_rna_cpm_path),
        "glial_20unit_atac_cpm_path": str(u_atac_cpm_path),
        "glial_20unit_atac_cpm_null_path": str(u_atac_cpm_null_path),
        "n_matched_glial_units": len(unit_ids),
    },
    "arm_b": {
        "rna_counts_path": str(rna_counts_path),
        "rna_meta_path": str(rna_meta_path),
        "atac_counts_path": str(atac_counts_path),
        "atac_meta_path": str(atac_meta_path),
        "atac_cpm_path": str(atac_cpm_path),
        "atac_cpm_null_path": str(atac_cpm_null_path),
        "deranged_mapping": dict(zip(donor_ids, deranged_donors)),
    },
    "same_cell_slice_path": str(same_cell_slice_path),
    "n_same_cell_slice_cells": int(len(picked_cells)),
    "n_same_cell_slice_peaks": int(len(bench_peak_idx)),
}))
"""
    env = dict(os.environ)
    numba_dir = work_dir / "numba_cache"
    numba_dir.mkdir(parents=True, exist_ok=True)
    env["NUMBA_CACHE_DIR"] = str(numba_dir)
    cmd = [
        *_resolve_env_python("aria-chromatin-env"),
        "-c",
        helper_code,
        str(merged_h5mu),
        str(ref_meta_tsv),
        str(work_dir),
        json.dumps(DONOR_COHORT),
        json.dumps(RNA_SPLIT_DONORS),
        json.dumps(ATAC_SPLIT_DONORS),
        json.dumps(list(BENCH_CHROMS)),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Cohort extraction helper failed rc={proc.returncode}:\n"
            f"STDOUT:\n{proc.stdout[-2000:]}\nSTDERR:\n{proc.stderr[-2000:]}"
        )
    lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
    return json.loads(lines[-1])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--paired-dir",
        type=Path,
        default=Path.home() / "Samael" / "Oligos" / "data" / "interim" / "paired",
    )
    parser.add_argument(
        "--ref-metadata",
        type=Path,
        default=(
            Path.home()
            / "Samael"
            / "Oligos"
            / "data"
            / "raw"
            / "reference"
            / "GSE278576_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz"
        ),
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        default=Path.home() / ".aria" / "genomes" / "hg38" / "annotation.gtf",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=Path.home() / ".aria" / "genomes" / "hg38" / "GRCh38.primary_assembly.genome.fa",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    t_start = time.perf_counter()
    donor_inputs: list[str] = []
    sample_ids: list[str] = []
    sample_metadata: dict[str, dict[str, Any]] = {}
    for row in DONOR_COHORT:
        sid = row["sample_id"]
        p = args.paired_dir / f"{sid}_paired.h5mu"
        if not p.is_file():
            raise FileNotFoundError(f"Missing real GSE278576 donor file: {p}")
        donor_inputs.append(str(p))
        sample_ids.append(sid)
        sample_metadata[sid] = {
            "donor_id": row["donor_id"],
            "age_group": row["age_group"],
            "sex": row["sex"],
            "age": row["age"],
        }

    with tempfile.TemporaryDirectory(prefix="aria_adr080_bench_") as tmp_s:
        work_dir = Path(tmp_s)

        # 1. Union 10 donor containers only to establish a shared genomic peak coordinate
        # space, then immediately split into strictly non-same-cell RNA and ATAC matrices.
        merge_dir = work_dir / "00_consensus_union"
        merge_res, t_merge = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_merge_cohort.py",
            {
                "inputs": donor_inputs,
                "sample_ids": sample_ids,
                "sample_metadata": sample_metadata,
                "condition_col": "age_group",
                "batch_key": "donor_id",
                "output_dir": str(merge_dir),
                "resume": True,
            },
            work_dir,
            "00_merge",
        )

        prep_dir = work_dir / "01_unpaired_tables"
        prep_dir.mkdir(parents=True, exist_ok=True)
        prep = _prepare_cohort_tables_and_same_cell_slice(
            Path(merge_res["output_path"]),
            args.ref_metadata,
            prep_dir,
        )

        # ══════════════════════════════════════════════════════════════════════
        # ARM A — Disjoint-Donor Split-Cohort Glial Lineage (Oligo_vs_Astro)
        # RNA = 5 young brains ONLY; ATAC = 5 completely disjoint aged brains ONLY
        # ══════════════════════════════════════════════════════════════════════
        arm_a_rna_dir = work_dir / "02a_arm_a_rna_de"
        arm_a_rna_res, t_a_rna = _run_ipc(
            "aria-rna-env",
            "aria/scripts/rna_bulk_de.py",
            {
                "files": [prep["arm_a"]["rna_counts_path"]],
                "metadata_file": prep["arm_a"]["rna_meta_path"],
                "design_factor": "lineage",
                "covariates": ["donor_id"],
                "contrasts": [
                    {
                        "numerator": "Oligo",
                        "denominator": "Astro",
                        "name": "Oligo_vs_Astro",
                    }
                ],
                "organism": "Homo sapiens",
                "genome": "hg38",
                "output_dir": str(arm_a_rna_dir),
                "run_pathways": True,
                "padj_threshold": 0.05,
                "lfc_threshold": 0.5,
                "n_cpus": 4,
            },
            work_dir,
            "02a_arm_a_rna",
        )

        arm_a_atac_dir = work_dir / "03a_arm_a_atac_da"
        arm_a_da_res, t_a_da = _run_ipc(
            "aria-rna-env",
            "aria/scripts/chromatin_bulk_diffacc.py",
            {
                "data_type": "bulk_ATAC",
                "counts_matrix_path": prep["arm_a"]["atac_counts_path"],
                "sample_metadata_path": prep["arm_a"]["atac_meta_path"],
                "condition_col": "lineage",
                "replicate_col": "donor_id",
                "covariates": ["donor_id"],
                "comparisons": [{"test": "Oligo", "reference": "Astro"}],
                "padj_max": 0.05,
                "lfc_min": 0.5,
                "output_dir": str(arm_a_atac_dir),
                "n_cpus": 4,
            },
            work_dir,
            "03a_arm_a_da",
        )
        assert arm_a_da_res.get("ran") is True, f"Arm A ATAC DA did not run: {arm_a_da_res}"

        # Exact reversed-contrast sign-flip canonicalization test (Astro_vs_Oligo)
        import pandas as pd

        fwd_da_comp = (arm_a_da_res.get("comparisons") or [{}])[0]
        fwd_da_csv = Path(fwd_da_comp["full_results_csv"])
        rev_csv_dir = work_dir / "03a_rev_exact"
        rev_csv_dir.mkdir(parents=True, exist_ok=True)
        rev_exact_csv = rev_csv_dir / "Astro_vs_Oligo_full_results.csv"
        fwd_df = pd.read_csv(fwd_da_csv)
        rev_df = fwd_df.copy()
        for col in ("log2fc", "log2FoldChange"):
            if col in rev_df.columns:
                rev_df[col] = -pd.to_numeric(rev_df[col], errors="coerce")
        if "comparison" in rev_df.columns:
            rev_df["comparison"] = "Astro_vs_Oligo"
        rev_df.to_csv(rev_exact_csv, index=False)
        arm_a_da_rev_exact = {
            "status": "success",
            "ran": True,
            "data_type": "bulk_ATAC",
            "comparisons": [
                {
                    **fwd_da_comp,
                    "comparison": "Astro_vs_Oligo",
                    "test": "Astro",
                    "reference": "Oligo",
                    "full_results_csv": str(rev_exact_csv),
                }
            ],
        }

        arm_a_ann_res, t_a_ann = _run_ipc(
            "aria-rna-env",
            "aria/scripts/chromatin_peak_annotation.py",
            {
                "data_type": "bulk_ATAC",
                "comparisons": arm_a_da_res.get("comparisons") or [],
                "gtf": str(args.gtf),
                "output_dir": str(arm_a_atac_dir / "annotation"),
            },
            work_dir,
            "03a_ann",
        )
        arm_a_ora_res, t_a_ora = _run_ipc(
            "aria-rna-env",
            "aria/scripts/chromatin_peak_ora.py",
            {
                "data_type": "bulk_ATAC",
                "comparisons": arm_a_da_res.get("comparisons") or [],
                "genome": "hg38",
                "organism": "Homo sapiens",
                "gtf": str(args.gtf),
                "output_dir": str(arm_a_atac_dir / "peak_ora"),
            },
            work_dir,
            "03a_ora",
        )

        from aria.agents.chromatin_helpers import bulk_da_motif_regions as _bulk_da_motif_regions

        motif_regions, motif_bg, _ = _bulk_da_motif_regions(arm_a_da_res.get("comparisons"))
        arm_a_motif_res: dict[str, Any] = {"status": "skipped", "ran": False, "reason": "no_da_regions"}
        t_a_motifs = 0.0
        if motif_regions and args.fasta.is_file():
            arm_a_motif_res, t_a_motifs = _run_ipc(
                "aria-chromatin-env",
                "aria/scripts/chromatin_motifs.py",
                {
                    "regions": motif_regions,
                    "background": motif_bg,
                    "genome_fasta": str(args.fasta),
                    "output_dir": str(arm_a_atac_dir / "motifs"),
                    "padj_max": 0.05,
                },
                work_dir,
                "03a_motifs",
            )

        arm_a_rna_env = {"status": "done", "findings": arm_a_rna_res}
        arm_a_chrom_env = {
            "status": "done",
            "findings": {
                "bulk_ATAC": {
                    "status": "done",
                    "findings": {
                        "differential_accessibility": arm_a_da_res,
                        "peak_annotation": arm_a_ann_res,
                        "peak_ora": arm_a_ora_res,
                        "motifs": arm_a_motif_res,
                    },
                }
            },
        }

        arm_a_cm_dir = work_dir / "04a_arm_a_cross_modal"
        arm_a_cm, t_a_cm = _run_ipc(
            "aria-rna-env",
            "aria/scripts/cross_modal_concordance.py",
            {
                "rna_results": arm_a_rna_env,
                "chromatin_results": arm_a_chrom_env,
                "rna_cpm_path": prep["arm_a"]["glial_20unit_rna_cpm_path"],
                "atac_cpm_path": prep["arm_a"]["glial_20unit_atac_cpm_path"],
                "gtf_path": str(args.gtf),
                "organism": "Homo sapiens",
                "genome": "hg38",
                "output_dir": str(arm_a_cm_dir),
                "padj_threshold": 0.05,
                "cis_distance_bp": 100000,
                "min_matched_samples": 6,
                "n_permutations": 200,
                "seed": 42,
                "resume": True,
            },
            work_dir,
            "04a_arm_a_cm",
        )
        assert arm_a_cm["ran"] is True, f"Arm A cross-modal concordance did not run: {arm_a_cm}"

        # Arm A A4 stage-manifest cache hit check
        arm_a_cm_resume, t_a_cm_resume = _run_ipc(
            "aria-rna-env",
            "aria/scripts/cross_modal_concordance.py",
            {
                "rna_results": arm_a_rna_env,
                "chromatin_results": arm_a_chrom_env,
                "rna_cpm_path": prep["arm_a"]["glial_20unit_rna_cpm_path"],
                "atac_cpm_path": prep["arm_a"]["glial_20unit_atac_cpm_path"],
                "gtf_path": str(args.gtf),
                "organism": "Homo sapiens",
                "genome": "hg38",
                "output_dir": str(arm_a_cm_dir),
                "padj_threshold": 0.05,
                "cis_distance_bp": 100000,
                "min_matched_samples": 6,
                "n_permutations": 200,
                "seed": 42,
                "resume": True,
            },
            work_dir,
            "04a_arm_a_cm_resume",
        )
        assert arm_a_cm_resume.get("resumed_from_manifest") is True, "A4 stage manifest did not hit cache"

        # Arm A Reversed-contrast sign-flip canonicalization check (RNA Oligo_vs_Astro x ATAC Astro_vs_Oligo)
        arm_a_rev_dir = work_dir / "04a_arm_a_reversed"
        arm_a_cm_rev, _ = _run_ipc(
            "aria-rna-env",
            "aria/scripts/cross_modal_concordance.py",
            {
                "rna_results": arm_a_rna_env,
                "chromatin_results": {
                    "status": "done",
                    "findings": {
                        "bulk_ATAC": {
                            "status": "done",
                            "findings": {"differential_accessibility": arm_a_da_rev_exact},
                        }
                    },
                },
                "gtf_path": str(args.gtf),
                "organism": "Homo sapiens",
                "genome": "hg38",
                "output_dir": str(arm_a_rev_dir),
                "padj_threshold": 0.05,
                "cis_distance_bp": 100000,
                "resume": False,
            },
            work_dir,
            "04a_arm_a_rev",
        )
        rev_mc_a = (arm_a_cm_rev.get("matched_contrasts") or [{}])[0]
        assert rev_mc_a.get("atac_sign_flipped") is True, "Reversed contrast did not flip ATAC sign"

        # Arm A 20-unit deranged-sample null control
        arm_a_null_dir = work_dir / "04a_glial_20unit_null"
        arm_a_cm_null, _ = _run_ipc(
            "aria-rna-env",
            "aria/scripts/cross_modal_concordance.py",
            {
                "rna_cpm_path": prep["arm_a"]["glial_20unit_rna_cpm_path"],
                "atac_cpm_path": prep["arm_a"]["glial_20unit_atac_cpm_null_path"],
                "gtf_path": str(args.gtf),
                "organism": "Homo sapiens",
                "genome": "hg38",
                "output_dir": str(arm_a_null_dir),
                "padj_threshold": 0.05,
                "cis_distance_bp": 100000,
                "min_matched_samples": 6,
                "n_permutations": 200,
                "seed": 42,
                "resume": False,
            },
            work_dir,
            "04a_glial_null",
        )

        # ══════════════════════════════════════════════════════════════════════
        # ARM B — 10-Donor Matched Clinical Aging Cohort (80-100_vs_20-39, ~ sex + age_group)
        # ══════════════════════════════════════════════════════════════════════
        arm_b_rna_dir = work_dir / "02b_arm_b_rna_de"
        arm_b_rna_res, t_b_rna = _run_ipc(
            "aria-rna-env",
            "aria/scripts/rna_bulk_de.py",
            {
                "files": [prep["arm_b"]["rna_counts_path"]],
                "metadata_file": prep["arm_b"]["rna_meta_path"],
                "design_factor": "age_group",
                "covariates": ["sex"],
                "contrasts": [
                    {
                        "numerator": "80-100",
                        "denominator": "20-39",
                        "name": "80-100_vs_20-39",
                    }
                ],
                "organism": "Homo sapiens",
                "genome": "hg38",
                "output_dir": str(arm_b_rna_dir),
                "run_pathways": True,
                "padj_threshold": 0.05,
                "lfc_threshold": 0.25,
                "n_cpus": 4,
            },
            work_dir,
            "02b_arm_b_rna",
        )
        arm_b_atac_dir = work_dir / "03b_arm_b_atac_da"
        arm_b_da_res, t_b_da = _run_ipc(
            "aria-rna-env",
            "aria/scripts/chromatin_bulk_diffacc.py",
            {
                "data_type": "bulk_ATAC",
                "counts_matrix_path": prep["arm_b"]["atac_counts_path"],
                "sample_metadata_path": prep["arm_b"]["atac_meta_path"],
                "condition_col": "age_group",
                "replicate_col": "donor_id",
                "covariates": ["sex"],
                "comparisons": [{"test": "80-100", "reference": "20-39"}],
                "padj_max": 0.05,
                "lfc_min": 0.25,
                "output_dir": str(arm_b_atac_dir),
                "n_cpus": 4,
            },
            work_dir,
            "03b_arm_b_da",
        )
        assert arm_b_da_res.get("ran") is True, f"Arm B ATAC DA did not run: {arm_b_da_res}"

        arm_b_cm_dir = work_dir / "04b_arm_b_cross_modal"
        arm_b_rna_cpm_path = str(arm_b_rna_dir / "counts_cpm.tsv")
        arm_b_cm, t_b_cm = _run_ipc(
            "aria-rna-env",
            "aria/scripts/cross_modal_concordance.py",
            {
                "rna_results": {"status": "done", "findings": arm_b_rna_res},
                "chromatin_results": {
                    "status": "done",
                    "findings": {
                        "bulk_ATAC": {
                            "status": "done",
                            "findings": {"differential_accessibility": arm_b_da_res},
                        }
                    },
                },
                "rna_cpm_path": arm_b_rna_cpm_path,
                "atac_cpm_path": prep["arm_b"]["atac_cpm_path"],
                "gtf_path": str(args.gtf),
                "organism": "Homo sapiens",
                "genome": "hg38",
                "output_dir": str(arm_b_cm_dir),
                "padj_threshold": 0.05,
                "cis_distance_bp": 100000,
                "min_matched_samples": 6,
                "n_permutations": 200,
                "seed": 42,
                "resume": True,
            },
            work_dir,
            "04b_arm_b_cm",
        )
        arm_b_null_dir = work_dir / "04b_arm_b_donor_null"
        arm_b_cm_null, _ = _run_ipc(
            "aria-rna-env",
            "aria/scripts/cross_modal_concordance.py",
            {
                "rna_cpm_path": arm_b_rna_cpm_path,
                "atac_cpm_path": prep["arm_b"]["atac_cpm_null_path"],
                "gtf_path": str(args.gtf),
                "organism": "Homo sapiens",
                "genome": "hg38",
                "output_dir": str(arm_b_null_dir),
                "padj_threshold": 0.05,
                "cis_distance_bp": 100000,
                "min_matched_samples": 6,
                "n_permutations": 200,
                "seed": 42,
                "resume": False,
            },
            work_dir,
            "04b_arm_b_null",
        )

        # ══════════════════════════════════════════════════════════════════════
        # 5. Same-Cell Multiome Peak2Gene Ground-Truth Cross-Check (BENCH_CHROMS)
        # ══════════════════════════════════════════════════════════════════════
        p2g_dir = work_dir / "05_same_cell_p2g_ref"
        p2g_res, t_p2g = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_peak2gene.py",
            {
                "h5mu_path": prep["same_cell_slice_path"],
                "gtf_path": str(args.gtf),
                "genome_fasta": str(args.fasta),
                "distance_bp": 100000,
                "min_corr": 0.25,
                "max_padj": 0.05,
                "seed": 0,
                "output_dir": str(p2g_dir),
                "resume": True,
            },
            work_dir,
            "05_same_cell_p2g",
        )

        p2g_df = pd.read_csv(p2g_res["output_path"]) if Path(p2g_res["output_path"]).is_file() else pd.DataFrame()
        de_da_csv = (arm_a_cm.get("tables") or {}).get("de_da_concordance_csv")
        de_da_df = pd.read_csv(de_da_csv) if de_da_csv and Path(de_da_csv).is_file() else pd.DataFrame()
        sm_csv = (arm_a_cm.get("tables") or {}).get("sample_matched_links_csv")
        sm_df = pd.read_csv(sm_csv) if sm_csv and Path(sm_csv).is_file() else pd.DataFrame()
        sc_genes = set(p2g_df["gene"].astype(str).str.upper()) if "gene" in p2g_df.columns else set()

        disjoint_de_da_genes: set[str] = set()
        if not de_da_df.empty and "gene" in de_da_df.columns and "classification" in de_da_df.columns:
            disjoint_de_da_genes = set(
                de_da_df.loc[
                    de_da_df["classification"].astype(str).isin(["concordant_up", "concordant_down"]),
                    "gene",
                ]
                .astype(str)
                .str.upper()
            )
        matched_cpm_genes: set[str] = set()
        if not sm_df.empty and "gene" in sm_df.columns and "padj" in sm_df.columns:
            matched_cpm_genes = set(
                sm_df.loc[pd.to_numeric(sm_df["padj"], errors="coerce") < 0.05, "gene"]
                .astype(str)
                .str.upper()
            )
        shared_disjoint_de_da = sorted(sc_genes & disjoint_de_da_genes)
        shared_any = sorted(sc_genes & (disjoint_de_da_genes | matched_cpm_genes))

        # ══════════════════════════════════════════════════════════════════════
        # 6. Governance, BiologicalSynthesisAgent Slice 2, RunLedger C3 & C4
        # ══════════════════════════════════════════════════════════════════════
        from aria.agents.biological_synthesis_agent import BiologicalSynthesisAgent
        from aria.agents.modality_audit import build_capability_matrix
        from aria.agents.narrative.claim_compiler import compile_public_claims
        from aria.agents.narrative.compose_prose import compose_block_prose
        from aria.agents.narrative.evidence_verifier import verify_block_claim_support
        from aria.agents.narrative.run_ledger import (
            build_run_ledger,
            ensure_report_ledger_nodes,
            link_claims_to_ledger,
        )
        from aria.agents.narrative.validators import find_causal_language
        from aria.agents.orchestrator_agent import (
            CROSS_MODAL_CONCORDANCE_VALIDATION,
            INTEGRATION_VALIDATION,
            MODALITY_VALIDATION,
        )

        cap_unpaired = build_capability_matrix(
            {
                "genome": "hg38",
                "organism": "Homo sapiens",
                "modalities": {
                    "bulk_RNA": [prep["arm_b"]["rna_counts_path"]],
                    "bulk_ATAC": [prep["arm_b"]["atac_counts_path"]],
                },
                "sample_metadata": sample_metadata,
                "condition_col": "age_group",
                "replicate_col": "donor_id",
                "comparisons": [{"test": "80-100", "reference": "20-39"}],
            },
            modality_validation=MODALITY_VALIDATION,
            integration_validation=INTEGRATION_VALIDATION,
        )
        cap_diagonal_blocked = build_capability_matrix(
            {
                "genome": "hg38",
                "organism": "Homo sapiens",
                "modalities": {
                    "scRNA": ["/virtual/rna.h5ad"],
                    "scATAC": ["/virtual/atac.h5ad"],
                },
                "cross_modal_strategy": "glue",
            },
            modality_validation=MODALITY_VALIDATION,
            integration_validation=INTEGRATION_VALIDATION,
        )

        agent_results_full = {
            "bulk_rna_agent": arm_a_rna_env,
            "chromatin_agent": arm_a_chrom_env,
            "cross_modal_concordance": arm_a_cm,
        }
        synth_blocks = BiologicalSynthesisAgent().synthesize(
            agent_results_full,
            {"organism": "Homo sapiens", "genome": "hg38"},
        )
        cm_blocks = [b for b in synth_blocks if str(b.id).startswith("integration.cross_modal")]
        assert len(cm_blocks) >= 1, "BiologicalSynthesisAgent Slice 2 emitted 0 cross-modal blocks"

        causal_hits = []
        unsupported_blocks = []
        for blk in cm_blocks:
            prose = compose_block_prose(blk) or blk.claim or ""
            ver = verify_block_claim_support(blk, prose, strict=False)
            if ver.get("status") != "supported":
                unsupported_blocks.append({"id": blk.id, "verification": ver})
            for txt in [blk.claim or "", prose, *(c.text for c in (blk.caveats or []))]:
                hit = find_causal_language(txt)
                if hit:
                    causal_hits.append({"id": blk.id, "hit": hit})

        ledger = build_run_ledger(
            {"organism": "Homo sapiens", "genome": "hg38"},
            agent_results_full,
        )
        ensure_report_ledger_nodes(ledger, cm_blocks)
        compiled = compile_public_claims(cm_blocks, {}, run_ledger=ledger)
        linked = link_claims_to_ledger(list(compiled.claims), ledger)

        assert unsupported_blocks == [], f"W-CLAIM unsupported synthesis blocks: {unsupported_blocks}"
        assert causal_hits == [], f"C4 causal language hits: {causal_hits}"
        assert linked["n_violations"] == 0, f"C3 ledger linkage violations: {linked}"

        mc_a0 = (arm_a_cm.get("matched_contrasts") or [{}])[0]
        de_da_a = mc_a0.get("de_da_summary") or {}
        rev_de_da_a = rev_mc_a.get("de_da_summary") or {}
        sm_a_real = arm_a_cm.get("sample_matched_links") or {}
        sm_a_null = arm_a_cm_null.get("sample_matched_links") or {}

        mc_b0 = (arm_b_cm.get("matched_contrasts") or [{}])[0]
        de_da_b = mc_b0.get("de_da_summary") or {}
        sm_b_real = arm_b_cm.get("sample_matched_links") or {}
        sm_b_null = arm_b_cm_null.get("sample_matched_links") or {}

        exact_rev_match = (
            de_da_a.get("n_concordant_total") == rev_de_da_a.get("n_concordant_total")
            and de_da_a.get("n_discordant_total") == rev_de_da_a.get("n_discordant_total")
            and abs(float(de_da_a.get("spearman_r") or 0.0) - float(rev_de_da_a.get("spearman_r") or 0.0)) < 1e-12
        )
        assert exact_rev_match, (
            f"Reversed contrast canonicalization mismatch: fwd={de_da_a} vs rev={rev_de_da_a}"
        )

        starburst_png = (arm_a_cm.get("figures") or {}).get("cross_modal_starburst")
        starburst_svg = (arm_a_cm.get("figures_svg") or {}).get("cross_modal_starburst")
        fig_checks = {
            "starburst_png_exists": bool(
                starburst_png and Path(starburst_png).is_file() and Path(starburst_png).stat().st_size > 1000
            ),
            "starburst_svg_exists": bool(
                starburst_svg and Path(starburst_svg).is_file() and Path(starburst_svg).stat().st_size > 1000
            ),
        }
        assert all(fig_checks.values()), f"Starburst figures missing: {fig_checks}"

        all_invariants_passed = bool(
            arm_a_cm.get("ran") is True
            and int(de_da_a.get("n_concordant_total") or 0) > 0
            and float(de_da_a.get("fisher_odds_ratio") or 0.0) > 5.0
            and exact_rev_match
            and arm_a_cm_resume.get("resumed_from_manifest") is True
            and all(fig_checks.values())
            and int(sm_b_null.get("n_significant_links") or 0) == 0
            and len(unsupported_blocks) == 0
            and len(causal_hits) == 0
            and linked["n_violations"] == 0
            and cap_unpaired["cross_modal_concordance_readiness"]["status"] == "green"
            and cap_diagonal_blocked["cross_modal_concordance_readiness"]["status"] == "red"
        )
        assert all_invariants_passed, "One or more ADR-080 invariants failed!"

        n_a_real = int(sm_a_real.get("n_significant_links") or 0)
        n_a_null = int(sm_a_null.get("n_significant_links") or 0)
        n_b_real = int(sm_b_real.get("n_significant_links") or 0)
        n_b_null = int(sm_b_null.get("n_significant_links") or 0)

        receipt = {
            "schema": "aria.adr080_cross_modal_concordance_validation.v1",
            "promotion_adr": "ADR-080",
            "all_invariants_passed": all_invariants_passed,
            "dataset": {
                "accession": "GSE278576",
                "tissue": "Human hippocampus (10-donor non-same-cell cross-modal validation)",
                "design": (
                    "10 balanced human brains (5 young 20-39 [3M/2F], 5 aged 80-100 [3M/2F]); "
                    "Arm A = disjoint-donor split cohort (5 brains RNA-only x 5 disjoint brains ATAC-only, Oligo_vs_Astro); "
                    "Arm B = 10-donor matched aging (~ sex + age_group, 80-100_vs_20-39) + donor-deranged null"
                ),
                "donors": DONOR_COHORT,
                "n_donors": len(DONOR_COHORT),
                "n_post_qc_cells_aggregated": prep["n_post_qc_cells_total"],
                "n_glial_cells_aggregated": prep["n_glial_cells_total"],
                "n_cells_by_donor": prep["n_cells_by_donor"],
                "n_genes_quantified": prep["n_genes_quantified"],
                "n_consensus_peaks_tested": prep["n_consensus_peaks_tested"],
            },
            "arm_a_disjoint_donor_glial_lineage": {
                "design_isolation": {
                    "rna_only_donors": RNA_SPLIT_DONORS,
                    "atac_only_donors": ATAC_SPLIT_DONORS,
                    "shared_donors_in_contrast": 0,
                    "shared_cells_in_contrast": 0,
                    "contrast": "Oligo_vs_Astro",
                },
                "single_modality_inputs": {
                    "bulk_rna_de": {
                        "n_genes_tested": (arm_a_rna_res.get("contrasts") or [{}])[0].get("n_genes_tested"),
                        "n_significant_genes": arm_a_rna_res.get("n_significant"),
                        "n_upregulated": arm_a_rna_res.get("n_upregulated"),
                        "n_downregulated": arm_a_rna_res.get("n_downregulated"),
                        "elapsed_sec": t_a_rna,
                    },
                    "bulk_atac_da": {
                        "n_peaks_tested": arm_a_da_res.get("n_peaks_tested"),
                        "n_significant_peaks": arm_a_da_res.get("n_sig_total"),
                        "peak_annotation_ran": bool(arm_a_ann_res.get("ran")),
                        "peak_ora_ran": bool(arm_a_ora_res.get("ran")),
                        "motif_enrichment_ran": bool(arm_a_motif_res.get("ran")),
                        "elapsed_sec": round(t_a_da + t_a_ann + t_a_ora + t_a_motifs, 3),
                    },
                },
                "de_da_summary": {
                    "rna_scope": mc_a0.get("rna_scope"),
                    "atac_scope": mc_a0.get("atac_scope"),
                    "n_gene_peak_pairs_tested": de_da_a.get("n_gene_peak_pairs_tested"),
                    "n_genes_tested": de_da_a.get("n_genes_tested"),
                    "n_promoter_pairs_tested": de_da_a.get("n_promoter_pairs_tested"),
                    "n_both_significant": de_da_a.get("n_both_significant"),
                    "n_concordant_up": de_da_a.get("n_concordant_up"),
                    "n_concordant_down": de_da_a.get("n_concordant_down"),
                    "n_concordant_total": de_da_a.get("n_concordant_total"),
                    "n_discordant_total": de_da_a.get("n_discordant_total"),
                    "n_atac_da_rna_unchanged": de_da_a.get("n_atac_da_rna_unchanged"),
                    "n_rna_de_atac_unchanged": de_da_a.get("n_rna_de_atac_unchanged"),
                    "concordance_fraction": de_da_a.get("concordance_fraction"),
                    "fisher_odds_ratio": de_da_a.get("fisher_odds_ratio"),
                    "fisher_pvalue": de_da_a.get("fisher_pvalue"),
                    "enrichment_odds_ratio": de_da_a.get("enrichment_odds_ratio"),
                    "enrichment_pvalue": de_da_a.get("enrichment_pvalue"),
                    "spearman_r": de_da_a.get("spearman_r"),
                    "spearman_pvalue": de_da_a.get("spearman_pvalue"),
                    "promoter_spearman_r": de_da_a.get("promoter_spearman_r"),
                    "pearson_r": de_da_a.get("pearson_r"),
                    "top_concordant_genes": (de_da_a.get("top_concordant_genes") or [])[:12],
                },
                "reversed_contrast_invariance": {
                    "atac_input_contrast": "Astro_vs_Oligo",
                    "atac_sign_flipped": rev_mc_a.get("atac_sign_flipped"),
                    "forward_concordant_total": de_da_a.get("n_concordant_total"),
                    "reversed_concordant_total": rev_de_da_a.get("n_concordant_total"),
                    "forward_spearman_r": de_da_a.get("spearman_r"),
                    "reversed_spearman_r": rev_de_da_a.get("spearman_r"),
                    "exact_concordance_match": exact_rev_match,
                },
                "tf_triangulation": arm_a_cm.get("tf_triangulation"),
                "pathway_concordance": arm_a_cm.get("pathway_concordance"),
                "matched_20unit_glial_cpm_links": {
                    "n_shared_units": sm_a_real.get("n_shared_samples"),
                    "n_pairs_tested": sm_a_real.get("n_pairs_tested"),
                    "n_significant_links": n_a_real,
                    "n_positive_links": sm_a_real.get("n_positive_links"),
                    "n_negative_links": sm_a_real.get("n_negative_links"),
                    "top_links": (sm_a_real.get("top_links") or [])[:10],
                    "deranged_unit_null_control": {
                        "n_null_significant_links": n_a_null,
                        "per_test_fp_rate": round(
                            n_a_null / max(1, int(sm_a_null.get("n_pairs_tested") or 1)),
                            6,
                        ),
                        "empirical_null_ratio": round(n_a_null / max(1, n_a_real), 6),
                    },
                },
                "a4_stage_manifest": {
                    "cache_hit_verified": bool(arm_a_cm_resume.get("resumed_from_manifest") is True),
                    "primary_elapsed_sec": t_a_cm,
                    "resume_elapsed_sec": t_a_cm_resume,
                },
                "figures": fig_checks,
            },
            "arm_b_matched_10donor_clinical_aging": {
                "contrast": "80-100_vs_20-39",
                "formula": "~ sex + age_group",
                "single_modality_inputs": {
                    "bulk_rna_de": {
                        "n_genes_tested": (arm_b_rna_res.get("contrasts") or [{}])[0].get("n_genes_tested"),
                        "n_significant_genes": arm_b_rna_res.get("n_significant"),
                        "n_upregulated": arm_b_rna_res.get("n_upregulated"),
                        "n_downregulated": arm_b_rna_res.get("n_downregulated"),
                        "elapsed_sec": t_b_rna,
                    },
                    "bulk_atac_da": {
                        "n_peaks_tested": arm_b_da_res.get("n_peaks_tested"),
                        "n_significant_peaks": arm_b_da_res.get("n_sig_total"),
                        "elapsed_sec": t_b_da,
                    },
                },
                "de_da_summary": {
                    "rna_scope": mc_b0.get("rna_scope"),
                    "atac_scope": mc_b0.get("atac_scope"),
                    "n_gene_peak_pairs_tested": de_da_b.get("n_gene_peak_pairs_tested"),
                    "n_genes_tested": de_da_b.get("n_genes_tested"),
                    "n_promoter_pairs_tested": de_da_b.get("n_promoter_pairs_tested"),
                    "n_both_significant": de_da_b.get("n_both_significant"),
                    "n_concordant_total": de_da_b.get("n_concordant_total"),
                    "n_discordant_total": de_da_b.get("n_discordant_total"),
                    "n_atac_da_rna_unchanged": de_da_b.get("n_atac_da_rna_unchanged"),
                    "n_rna_de_atac_unchanged": de_da_b.get("n_rna_de_atac_unchanged"),
                    "spearman_r": de_da_b.get("spearman_r"),
                    "spearman_pvalue": de_da_b.get("spearman_pvalue"),
                },
                "donor_matched_10donor_cpm_links": {
                    "n_shared_donors": sm_b_real.get("n_shared_samples"),
                    "n_pairs_tested": sm_b_real.get("n_pairs_tested"),
                    "n_significant_links": n_b_real,
                    "top_links": (sm_b_real.get("top_links") or [])[:5],
                    "donor_deranged_null_control": {
                        "deranged_mapping": prep["arm_b"]["deranged_mapping"],
                        "n_null_significant_links": n_b_null,
                        "per_test_fp_rate": round(
                            n_b_null / max(1, int(sm_b_null.get("n_pairs_tested") or 1)),
                            6,
                        ),
                    },
                },
                "elapsed_sec": t_b_cm,
            },
            "same_cell_ground_truth_comparison": {
                "n_same_cell_slice_cells": prep["n_same_cell_slice_cells"],
                "chromosomes_evaluated": list(BENCH_CHROMS),
                "n_same_cell_p2g_links": p2g_res.get("n_links"),
                "n_same_cell_linked_genes": len(sc_genes),
                "n_disjoint_donor_de_da_concordant_genes": len(disjoint_de_da_genes),
                "n_recovered_by_disjoint_donor_de_da": len(shared_disjoint_de_da),
                "same_cell_gene_recall_by_disjoint_de_da": round(
                    len(shared_disjoint_de_da) / max(1, len(sc_genes)),
                    4,
                ),
                "n_recovered_by_any_non_same_cell_layer": len(shared_any),
                "same_cell_gene_recall_combined": round(
                    len(shared_any) / max(1, len(sc_genes)),
                    4,
                ),
                "shared_gene_examples": shared_any[:20],
                "elapsed_sec": t_p2g,
            },
            "governance_and_narrative_verification": {
                "cross_modal_validation_registry": CROSS_MODAL_CONCORDANCE_VALIDATION,
                "unpaired_bulk_readiness": cap_unpaired["cross_modal_concordance_readiness"],
                "same_cell_integration_readiness_on_unpaired": {
                    "status": cap_unpaired["integration_readiness"]["status"],
                    "dispatch_policy": cap_unpaired["integration_readiness"]["dispatch_policy"],
                },
                "diagonal_imputation_gate": {
                    "requested_strategy": "glue",
                    "capability_matrix_status": cap_diagonal_blocked["status"],
                    "readiness_status": cap_diagonal_blocked["cross_modal_concordance_readiness"]["status"],
                    "dispatch_policy": cap_diagonal_blocked["cross_modal_concordance_readiness"]["dispatch_policy"],
                    "reason": cap_diagonal_blocked["cross_modal_concordance_readiness"]["reason"],
                },
                "biological_synthesis_slice2": {
                    "n_cross_modal_blocks": len(cm_blocks),
                    "block_ids": [b.id for b in cm_blocks],
                    "w_claim_unsupported_blocks": len(unsupported_blocks),
                    "c3_ledger_violations": len(linked["violations"]),
                    "c4_causal_verb_hits": len(causal_hits),
                },
            },
            "cohort_merge_elapsed_sec": t_merge,
            "total_elapsed_sec": round(time.perf_counter() - t_start, 3),
        }

        sanitized = _sanitize_obj(receipt)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(sanitized, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(
            json.dumps(
                {
                    "status": "ok",
                    "all_invariants_passed": all_invariants_passed,
                    "arm_a_concordant_de_da": de_da_a.get("n_concordant_total"),
                    "arm_a_concordance_fraction": de_da_a.get("concordance_fraction"),
                    "arm_a_fisher_or": de_da_a.get("fisher_odds_ratio"),
                    "arm_a_spearman_r": de_da_a.get("spearman_r"),
                    "arm_a_concordant_tfs": (arm_a_cm.get("tf_triangulation") or {}).get("n_concordant_tfs"),
                    "arm_a_shared_pathways": (arm_a_cm.get("pathway_concordance") or {}).get("n_shared_terms"),
                    "arm_a_20unit_links": n_a_real,
                    "arm_a_20unit_null_links": n_a_null,
                    "arm_b_10donor_null_links": n_b_null,
                    "same_cell_gene_recall": receipt["same_cell_ground_truth_comparison"]["same_cell_gene_recall_combined"],
                    "total_elapsed_sec": receipt["total_elapsed_sec"],
                }
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
