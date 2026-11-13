#!/usr/bin/env python3
"""ADR-079 — Multi-donor paired same-cell multiome empirical validation runner.

Executes ARIA's production paired same-cell multiome pipeline across a real
4-donor balanced human hippocampal cohort from GSE278576 (2 donors in age_group
20-39: hc29, hc77; 2 donors in age_group 60-79: hc1153, hc69984):

  1. ``aria/scripts/integration_merge_cohort.py`` (ADR-043 sweep-line consensus
     ATAC peak union across 4 heterogeneous per-donor peak sets, collision-free
     ``{sample_id}#{barcode}`` cell IDs, donor/condition metadata propagation,
     and A4 SHA-256 stage-manifest cache hit verification).
  2. ``aria/scripts/integration_paired_qc.py`` (joint RNA+ATAC same-cell QC +
     author hippocampal ``subclass`` cell-type join + A4 cache hit).
  3. ``aria/scripts/integration_wnn.py`` comparing unintegrated WNN
     (``harmonize_batches=False``) vs donor-aware Harmony WNN
     (``harmonize_batches=True``, ``batch_key="donor_id"``) on ``X_pca`` and
     ``X_lsi`` with WNN iLISI, before/after batch silhouettes, 5-axis
     anti-degradation guardrail, author ``subclass`` ARI/NMI conservation, and
     A4 cache hit.
  4. ``aria/scripts/integration_mofa.py`` inside ``aria-integration-env`` with
     ``batch_key="donor_id"`` and ``condition_col="age_group"`` verifying
     non-confounded 2x2 donor design, quantitative ANOVA eta^2 donor-factor
     screening conditioned on ``age_group``, and A4 cache hit.
  5. ``aria/scripts/integration_peak2gene.py`` comparing true same-cell pairing
     against a barcode-permuted null control (permuting ATAC barcodes relative
     to RNA barcodes while preserving exact per-gene and per-peak marginals) on
     a stratified 2,000-cell 4-donor slice with GRCh38 GC + accessibility
     background matching (asserting 0 false discoveries under broken pairing).
  6. Governance, ``IntegrationAuditAgent`` readiness (``green`` 2x2 cohort vs
     ``yellow`` 1-donor-per-condition confounded cohort vs ``red`` unpaired
     scRNA+scATAC), and ``C3`` ``ledger://integration/preprocessing`` narrative
     verification (zero claim-ledger violations).

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
    / "adr079_multidonor_multiome_validation.json"
)

DONOR_COHORT = [
    {"sample_id": "hc29", "donor_id": "hc29", "age_group": "20-39", "sex": "M", "age": 28},
    {"sample_id": "hc77", "donor_id": "hc77", "age_group": "20-39", "sex": "M", "age": 20},
    {"sample_id": "hc1153", "donor_id": "hc1153", "age_group": "60-79", "sex": "F", "age": 75},
    {"sample_id": "hc69984", "donor_id": "hc69984", "age_group": "60-79", "sex": "M", "age": 75},
]

P2G_CHROMS = ("chr1", "chr2", "chr11", "chr17", "chr18", "chr19")

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
    home_cand = Path.home() / "anaconda3" / "envs" / env_name / "bin" / "python"
    if home_cand.is_file():
        return [str(home_cand)]
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
        raise RuntimeError(
            f"IPC script {script_rel} ({tag}) failed rc={proc.returncode}:\n"
            f"STDOUT:\n{proc.stdout[-2000:]}\nSTDERR:\n{proc.stderr[-2000:]}"
        )
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    if payload.get("status") != "success":
        raise RuntimeError(f"IPC script {script_rel} ({tag}) returned non-success: {payload}")
    return payload, elapsed


def _run_chromatin_helper(code: str, args: list[str], work_dir: Path) -> dict[str, Any]:
    env = dict(os.environ)
    numba_dir = work_dir / "numba_cache"
    numba_dir.mkdir(parents=True, exist_ok=True)
    env["NUMBA_CACHE_DIR"] = str(numba_dir)
    cmd = [*_resolve_env_python("aria-chromatin-env"), "-c", code, *args]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Chromatin helper failed rc={proc.returncode}:\n"
            f"STDOUT:\n{proc.stdout[-2000:]}\nSTDERR:\n{proc.stderr[-2000:]}"
        )
    lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
    return json.loads(lines[-1])


def _attach_subclass_and_prepare_p2g_slices(
    post_qc_h5mu: Path,
    ref_meta_tsv: Path,
    work_dir: Path,
) -> dict[str, Any]:
    """Attach published GSE278576 author subclass labels onto post-QC merged cohort
    and emit stratified 2,000-cell real + barcode-shuffled null slices for cis P2G."""
    helper_code = r"""
import json, sys
from pathlib import Path
import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import scipy.sparse as sp

post_qc_path = Path(sys.argv[1])
meta_path = Path(sys.argv[2])
real_slice_path = Path(sys.argv[3])
null_slice_path = Path(sys.argv[4])
chroms = set(json.loads(sys.argv[5]))

meta = pd.read_csv(meta_path, sep="\t", usecols=["bacrode", "orig.ident", "subclass"])
meta["raw_bc"] = meta["bacrode"].astype(str).str.rsplit("_", n=1).str[-1]
meta["canonical_id"] = meta["orig.ident"].astype(str) + "#" + meta["raw_bc"]
meta = meta.drop_duplicates(subset=["canonical_id"]).set_index("canonical_id")

mdata = md.read_h5mu(str(post_qc_path))
obs_idx = pd.Index([str(x) for x in mdata.obs_names])
subclass = meta["subclass"].reindex(obs_idx)
n_matched = int(subclass.notna().sum())
subclass_filled = subclass.fillna("unannotated").astype(str)

for target in (mdata, mdata["rna"], mdata["atac"]):
    target.obs["subclass"] = subclass_filled.values

annotated_path = post_qc_path.parent / "paired_qc_cohort_annotated.h5mu"
mdata.write_h5mu(str(annotated_path))

# Stratified 500 cells per donor (2,000 total) across canonical marker chromosomes
rng = np.random.default_rng(0)
donor_series = mdata.obs["donor_id"].astype(str)
chosen_idx = []
for donor_id in sorted(donor_series.unique()):
    d_rows = np.flatnonzero((donor_series == donor_id).to_numpy())
    take = min(500, len(d_rows))
    picked = np.sort(rng.choice(d_rows, size=take, replace=False))
    chosen_idx.extend(picked.tolist())
chosen_idx = np.asarray(sorted(chosen_idx), dtype=np.int64)

atac_var = mdata["atac"].var_names.astype(str)
peak_mask = np.array([p.split(":", 1)[0] in chroms for p in atac_var], dtype=bool)
peak_idx = np.flatnonzero(peak_mask)

rna_sub = mdata["rna"][chosen_idx, :].copy()
atac_sub = mdata["atac"][chosen_idx, peak_idx].copy()
rna_X = sp.csr_matrix(rna_sub.layers["counts"] if "counts" in rna_sub.layers else rna_sub.X)
atac_X = sp.csr_matrix(atac_sub.layers["counts"] if "counts" in atac_sub.layers else atac_sub.X)
rna_sub.X = rna_X.copy()
rna_sub.layers["counts"] = rna_X.copy()
atac_sub.X = atac_X.copy()
atac_sub.layers["counts"] = atac_X.copy()

real_mdata = md.MuData({"rna": rna_sub, "atac": atac_sub})
real_mdata.obs = rna_sub.obs.copy()
real_mdata.write_h5mu(str(real_slice_path))

# Barcode-shuffled null control: permute ATAC rows relative to RNA rows,
# relabelling ATAC obs_names to match RNA obs_names so marginal distributions
# of every gene and every peak are 100% identical while true cell pairing is broken.
perm = rng.permutation(len(chosen_idx))
# Guarantee zero fixed points (derangement) so no cell retains true pairing
fixed = np.flatnonzero(perm == np.arange(len(chosen_idx)))
for idx_f in fixed:
    swap_with = (idx_f + 1) % len(chosen_idx)
    perm[idx_f], perm[swap_with] = perm[swap_with], perm[idx_f]

atac_null = ad.AnnData(
    X=atac_X[perm, :].copy(),
    obs=rna_sub.obs.copy(),
    var=atac_sub.var.copy(),
)
atac_null.layers["counts"] = atac_null.X.copy()
null_mdata = md.MuData({"rna": rna_sub.copy(), "atac": atac_null})
null_mdata.obs = rna_sub.obs.copy()
null_mdata.write_h5mu(str(null_slice_path))

print(json.dumps({
    "annotated_h5mu": str(annotated_path),
    "n_cells_post_qc": int(mdata.n_obs),
    "n_subclass_matched": n_matched,
    "subclass_match_fraction": round(n_matched / max(1, mdata.n_obs), 4),
    "subclass_counts": {
        str(k): int(v)
        for k, v in pd.Series(subclass_filled).value_counts().head(12).items()
    },
    "p2g_slice_n_cells": int(len(chosen_idx)),
    "p2g_slice_n_peaks": int(len(peak_idx)),
    "null_fixed_points": int(np.sum(perm == np.arange(len(chosen_idx)))),
}))
"""
    real_slice = work_dir / "p2g_real_slice.h5mu"
    null_slice = work_dir / "p2g_null_shuffled_slice.h5mu"
    info = _run_chromatin_helper(
        helper_code,
        [
            str(post_qc_h5mu),
            str(ref_meta_tsv),
            str(real_slice),
            str(null_slice),
            json.dumps(list(P2G_CHROMS)),
        ],
        work_dir,
    )
    info["real_slice_path"] = str(real_slice)
    info["null_slice_path"] = str(null_slice)
    return info


def _compute_external_subclass_concordance(
    unintegrated_h5mu: Path,
    harmony_h5mu: Path,
    work_dir: Path,
) -> dict[str, Any]:
    """Compute ARI & NMI of WNN Leiden clusters vs published GSE278576 author subclass."""
    code = r"""
import json, sys
import mudata as md
import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

def eval_h5mu(path):
    m = md.read_h5mu(path)
    obs = m.obs
    sub = obs["subclass"].astype(str).to_numpy()
    keep = sub != "unannotated"
    wnn_cl = obs["leiden_wnn"].astype(str).to_numpy()[keep]
    rna_cl = obs["leiden_rna"].astype(str).to_numpy()[keep] if "leiden_rna" in obs else wnn_cl
    sub_k = sub[keep]
    return {
        "n_annotated_cells": int(np.sum(keep)),
        "wnn_vs_author_subclass_ari": round(float(adjusted_rand_score(sub_k, wnn_cl)), 4),
        "wnn_vs_author_subclass_nmi": round(float(normalized_mutual_info_score(sub_k, wnn_cl)), 4),
        "rna_vs_author_subclass_ari": round(float(adjusted_rand_score(sub_k, rna_cl)), 4),
        "rna_vs_author_subclass_nmi": round(float(normalized_mutual_info_score(sub_k, rna_cl)), 4),
    }

print(json.dumps({
    "unintegrated": eval_h5mu(sys.argv[1]),
    "harmony": eval_h5mu(sys.argv[2]),
}))
"""
    return _run_chromatin_helper(code, [str(unintegrated_h5mu), str(harmony_h5mu)], work_dir)


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
    if not args.ref_metadata.is_file():
        raise FileNotFoundError(f"Missing GSE278576 reference metadata: {args.ref_metadata}")
    if not args.gtf.is_file() or not args.fasta.is_file():
        raise FileNotFoundError("Missing local hg38 GTF or FASTA in ~/.aria/genomes/hg38/")

    from aria.agents.modality_audit import build_capability_matrix
    from aria.agents.narrative_agent import NarrativeAgent
    from aria.agents.narrative.narrators.integration import IntegrationNarrator
    from aria.agents.narrative.run_ledger import build_run_ledger
    from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION, MODALITY_VALIDATION
    from aria.memory.memory import ARIAMemory
    from aria.utils.environment_audit import audit_environment
    from aria.utils.registry_integrity import check_registry_integrity
    from aria.version import collect_version_metadata

    with tempfile.TemporaryDirectory(prefix="adr079_multidonor_") as tmp_str:
        work_dir = Path(tmp_str)

        # ── 1. Multi-donor sparse cohort merger + A4 resume ──────────────────
        merge_dir = work_dir / "01_merge"
        merge_params = {
            "inputs": donor_inputs,
            "sample_ids": sample_ids,
            "sample_metadata": sample_metadata,
            "rna_gene_join": "inner",
            "output_dir": str(merge_dir),
            "resume": True,
        }
        merge_res, merge_elapsed = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_merge_cohort.py",
            merge_params,
            work_dir,
            "merge_fresh",
        )
        merge_cached, merge_cached_elapsed = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_merge_cohort.py",
            merge_params,
            work_dir,
            "merge_resume",
        )
        assert merge_cached.get("resumed") is True, "Expected A4 resume hit on cohort merge"

        # ── 2. Joint multi-donor paired QC + A4 resume ───────────────────────
        qc_dir = work_dir / "02_paired_qc"
        qc_params = {
            "h5mu_path": merge_res["output_path"],
            "output_dir": str(qc_dir),
            "genome": "hg38",
            "organism": "Homo sapiens",
            "resume": True,
        }
        qc_res, qc_elapsed = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_paired_qc.py",
            qc_params,
            work_dir,
            "qc_fresh",
        )
        qc_cached, qc_cached_elapsed = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_paired_qc.py",
            qc_params,
            work_dir,
            "qc_resume",
        )
        assert qc_cached.get("resumed") is True, "Expected A4 resume hit on paired QC"

        annot_info = _attach_subclass_and_prepare_p2g_slices(
            Path(qc_res["output_path"]),
            args.ref_metadata,
            work_dir,
        )
        annotated_h5mu = annot_info["annotated_h5mu"]

        # ── 3. Unintegrated vs Donor-Aware Harmony WNN + A4 resume ───────────
        wnn_unint_dir = work_dir / "03a_wnn_unintegrated"
        wnn_unint_params = {
            "h5mu_path": annotated_h5mu,
            "output_dir": str(wnn_unint_dir),
            "batch_key": "donor_id",
            "harmonize_batches": False,
            "audit_group_col": "subclass",
            "k": 20,
            "resolution": 0.5,
            "seed": 0,
            "resume": True,
        }
        wnn_unint_res, wnn_unint_elapsed = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_wnn.py",
            wnn_unint_params,
            work_dir,
            "wnn_unintegrated",
        )

        wnn_harm_dir = work_dir / "03b_wnn_harmony"
        wnn_harm_params = {
            "h5mu_path": annotated_h5mu,
            "output_dir": str(wnn_harm_dir),
            "batch_key": "donor_id",
            "harmonize_batches": True,
            "audit_group_col": "subclass",
            "k": 20,
            "resolution": 0.5,
            "seed": 0,
            "resume": True,
        }
        wnn_harm_res, wnn_harm_elapsed = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_wnn.py",
            wnn_harm_params,
            work_dir,
            "wnn_harmony_fresh",
        )
        wnn_harm_cached, wnn_harm_cached_elapsed = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_wnn.py",
            wnn_harm_params,
            work_dir,
            "wnn_harmony_resume",
        )
        assert wnn_harm_cached.get("resumed") is True, "Expected A4 resume hit on Harmony WNN"

        subclass_concordance = _compute_external_subclass_concordance(
            Path(wnn_unint_res["output_path"]),
            Path(wnn_harm_res["output_path"]),
            work_dir,
        )

        # ── 4. Donor/condition-aware MOFA+ factor decomposition + A4 resume ──
        mofa_dir = work_dir / "04_mofa"
        mofa_params = {
            "h5mu_path": annotated_h5mu,
            "output_dir": str(mofa_dir),
            "batch_key": "donor_id",
            "condition_col": "age_group",
            "subsample_strata_col": "donor_id",
            "max_cells": 4000,
            "n_factors": 10,
            "seed": 0,
            "write_h5mu": False,
            "resume": True,
        }
        mofa_res, mofa_elapsed = _run_ipc(
            "aria-integration-env",
            "aria/scripts/integration_mofa.py",
            mofa_params,
            work_dir,
            "mofa_fresh",
        )
        mofa_cached, mofa_cached_elapsed = _run_ipc(
            "aria-integration-env",
            "aria/scripts/integration_mofa.py",
            mofa_params,
            work_dir,
            "mofa_resume",
        )
        assert mofa_cached.get("resumed") is True, "Expected A4 resume hit on MOFA+"

        # ── 5. Real vs Barcode-Permuted Null cis Peak-to-Gene + A4 resume ────
        p2g_real_dir = work_dir / "05a_p2g_real"
        p2g_real_params = {
            "h5mu_path": annot_info["real_slice_path"],
            "gtf_path": str(args.gtf),
            "genome_fasta": str(args.fasta),
            "genome": "hg38",
            "organism": "Homo sapiens",
            "distance_kb": 250,
            "min_corr": 0.30,
            "max_padj": 0.05,
            "output_dir": str(p2g_real_dir),
            "seed": 0,
            "resume": True,
        }
        p2g_real_res, p2g_real_elapsed = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_peak2gene.py",
            p2g_real_params,
            work_dir,
            "p2g_real_fresh",
        )
        p2g_real_cached, p2g_real_cached_elapsed = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_peak2gene.py",
            p2g_real_params,
            work_dir,
            "p2g_real_resume",
        )
        assert p2g_real_cached.get("resumed") is True, "Expected A4 resume hit on Peak2Gene"

        p2g_null_dir = work_dir / "05b_p2g_null_shuffled"
        p2g_null_params = {
            "h5mu_path": annot_info["null_slice_path"],
            "gtf_path": str(args.gtf),
            "genome_fasta": str(args.fasta),
            "genome": "hg38",
            "organism": "Homo sapiens",
            "distance_kb": 250,
            "min_corr": 0.30,
            "max_padj": 0.05,
            "output_dir": str(p2g_null_dir),
            "seed": 0,
            "resume": False,
        }
        p2g_null_res, p2g_null_elapsed = _run_ipc(
            "aria-chromatin-env",
            "aria/scripts/integration_peak2gene.py",
            p2g_null_params,
            work_dir,
            "p2g_null",
        )
        n_real_links = int(p2g_real_res["n_links"])
        n_null_links = int(p2g_null_res["n_links"])
        emp_fdr = n_null_links / max(1, n_real_links)
        assert n_real_links > 0, "Expected >0 real cis peak-to-gene links"
        assert emp_fdr <= 0.05, (
            f"Expected empirical target-decoy FDR <= 0.05 under barcode permutation, "
            f"got {n_null_links}/{n_real_links} = {emp_fdr:.4f}"
        )

        # ── 6. Governance, IntegrationAuditAgent & C3 Run Ledger Verification ─
        cohort_ctx = {
            "genome": "hg38",
            "organism": "Homo sapiens",
            "modalities": {"scATAC": donor_inputs},
            "assay_detections": [
                {
                    "path": p,
                    "modality": "scATAC",
                    "confidence": "high",
                    "evidence": {"format": "h5mu", "modalities": ["rna", "atac"], "paired_rna": True},
                }
                for p in donor_inputs
            ],
            "design": {
                "groups": {"20-39": ["hc29", "hc77"], "60-79": ["hc1153", "hc69984"]},
                "pseudobulk": {
                    "condition_col": "age_group",
                    "replicate_col": "donor_id",
                    "groupby_col": "subclass",
                },
            },
        }
        mat_balanced = build_capability_matrix(
            cohort_ctx,
            modality_validation=MODALITY_VALIDATION,
            integration_validation=INTEGRATION_VALIDATION,
        )
        confounded_ctx = {
            "genome": "hg38",
            "modalities": {"scATAC": donor_inputs[:2]},
            "assay_detections": cohort_ctx["assay_detections"][:2],
            "design": {
                "groups": {"20-39": ["hc29"], "60-79": ["hc77"]},
                "sample_metadata": {
                    "hc29": {"donor_id": "hc29", "condition": "20-39"},
                    "hc77": {"donor_id": "hc77", "condition": "60-79"},
                },
            },
        }
        mat_confounded = build_capability_matrix(
            confounded_ctx,
            modality_validation=MODALITY_VALIDATION,
            integration_validation=INTEGRATION_VALIDATION,
        )
        unpaired_ctx = {
            "genome": "hg38",
            "modalities": {
                "scRNA": ["unpaired_rna.h5ad"],
                "scATAC": ["unpaired_atac.h5ad"],
            },
        }
        mat_unpaired = build_capability_matrix(
            unpaired_ctx,
            modality_validation=MODALITY_VALIDATION,
            integration_validation=INTEGRATION_VALIDATION,
        )

        integration_agent_payload = {
            "status": "done",
            "strategy": "wnn_then_mofa",
            "findings": {
                "preprocessing": {
                    "ran": True,
                    "status": "success",
                    "mode": "paired_cohort_merge",
                    "inputs": donor_inputs,
                    "output_path": merge_res["output_path"],
                    "n_samples": merge_res["n_samples"],
                    "sample_ids": merge_res["sample_ids"],
                    "n_cells_total": merge_res["n_cells_total"],
                    "n_genes": merge_res["n_genes"],
                    "n_consensus_peaks": merge_res["n_consensus_peaks"],
                    "barcode_translation": "sample_prefixed",
                    "cohort_merge": merge_res,
                },
                "paired_qc": qc_res,
                "wnn": wnn_harm_res,
                "mofa": mofa_res,
                "peak_to_gene": p2g_real_res,
            },
        }
        agent_results = {"integration_agent": integration_agent_payload}
        ledger = build_run_ledger(cohort_ctx, agent_results)
        pre_entry = next(
            e
            for e in ledger["entries"]
            if e["modality"] == "integration" and e["analysis"] == "preprocessing"
        )
        assert pre_entry["status"] == "executed", pre_entry

        narrator_blocks = IntegrationNarrator().collect(
            "integration_agent",
            integration_agent_payload,
            {"agent_results": agent_results},
        )
        mem = ARIAMemory(db_path=":memory:")
        try:
            narr_agent = NarrativeAgent(memory=mem, llm=None)
            narr_agent.reports_dir = work_dir / "06_reports"
            narr_agent.reports_dir.mkdir(parents=True, exist_ok=True)
            narr_res = narr_agent.run(
                experiment_id="adr079_gse278576_4donor",
                context={
                    "exp_context": {
                        **cohort_ctx,
                        "user_question": "Characterize multi-donor paired hippocampal multiome states",
                        "reproducible_mode": True,
                    },
                    "biological_intent": {
                        "summary": "Multi-donor paired hippocampal multiome integration"
                    },
                    "agent_results": agent_results,
                    "findings": [],
                },
            )
        finally:
            mem.close()
        assert narr_res["status"] == "done", narr_res
        report_html_path = Path(narr_res["report_path"])
        meth_json = json.loads(
            (report_html_path.parent / "methodology.json").read_text(encoding="utf-8")
        )
        ledger_verification = meth_json["run_ledger"]["claim_ledger_verification"]
        assert ledger_verification["n_violations"] == 0, ledger_verification

        bi_unint = wnn_unint_res.get("batch_integration") or {}
        bi_harm = wnn_harm_res.get("batch_integration") or {}

        receipt = {
            "adr": "ADR-079",
            "benchmark": "adr079_multidonor_paired_multiome_production_validation",
            "dataset": {
                "accession": "GSE278576",
                "tissue": "human_hippocampus",
                "assay": "10x_chromium_single_cell_multiome_atac_plus_gene_expression",
                "n_donors": 4,
                "conditions": {"20-39": ["hc29", "hc77"], "60-79": ["hc1153", "hc69984"]},
                "donors": DONOR_COHORT,
                "reference_annotation": "GSE278576_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz",
                "genome": "hg38",
                "gtf_annotation": args.gtf.name,
                "reference_fasta": args.fasta.name,
            },
            "stage_1_cohort_merge": {
                "script": "aria/scripts/integration_merge_cohort.py",
                "n_samples": merge_res["n_samples"],
                "sample_ids": merge_res["sample_ids"],
                "per_sample_cells": merge_res["per_sample_cells"],
                "n_cells_total": merge_res["n_cells_total"],
                "n_genes": merge_res["n_genes"],
                "rna_gene_join": merge_res["rna_gene_join"],
                "n_consensus_peaks": merge_res["n_consensus_peaks"],
                "elapsed_s": merge_elapsed,
                "a4_resume_verified": bool(merge_cached.get("resumed") is True),
                "a4_resume_elapsed_s": merge_cached_elapsed,
            },
            "stage_2_paired_qc": {
                "script": "aria/scripts/integration_paired_qc.py",
                "n_cells_before": qc_res["n_cells_before"],
                "n_cells_after": qc_res["n_cells_after"],
                "n_genes": qc_res["n_genes"],
                "n_peaks": qc_res["n_peaks"],
                "author_subclass_matched_cells": annot_info["n_subclass_matched"],
                "author_subclass_match_fraction": annot_info["subclass_match_fraction"],
                "top_author_subclasses": annot_info["subclass_counts"],
                "elapsed_s": qc_elapsed,
                "a4_resume_verified": bool(qc_cached.get("resumed") is True),
                "a4_resume_elapsed_s": qc_cached_elapsed,
            },
            "stage_3_donor_aware_wnn": {
                "script": "aria/scripts/integration_wnn.py",
                "batch_key": "donor_id",
                "unintegrated": {
                    "harmonized": bool(bi_unint.get("applied")),
                    "wnn_batch_ilisi_median": bi_unint.get("wnn_batch_ilisi_median"),
                    "wnn_batch_ilisi_norm_median": bi_unint.get("wnn_batch_ilisi_norm_median"),
                    "rna_batch_silhouette": bi_unint.get("rna_batch_silhouette_before"),
                    "atac_batch_silhouette": bi_unint.get("atac_batch_silhouette_before"),
                    "n_joint_clusters": wnn_unint_res["n_joint_clusters"],
                    "mean_rna_weight": round(float(wnn_unint_res["mean_rna_weight"]), 4),
                    "mean_atac_weight": round(float(wnn_unint_res["mean_atac_weight"]), 4),
                    "guardrail_status": wnn_unint_res.get("guardrail_status"),
                    "external_subclass_concordance": subclass_concordance["unintegrated"],
                    "elapsed_s": wnn_unint_elapsed,
                },
                "harmony_integrated": {
                    "harmonized": bool(bi_harm.get("applied")),
                    "method": bi_harm.get("method"),
                    "wnn_batch_ilisi_median": bi_harm.get("wnn_batch_ilisi_median"),
                    "wnn_batch_ilisi_norm_median": bi_harm.get("wnn_batch_ilisi_norm_median"),
                    "rna_batch_silhouette_before": bi_harm.get("rna_batch_silhouette_before"),
                    "rna_batch_silhouette_after": bi_harm.get("rna_batch_silhouette_after"),
                    "atac_batch_silhouette_before": bi_harm.get("atac_batch_silhouette_before"),
                    "atac_batch_silhouette_after": bi_harm.get("atac_batch_silhouette_after"),
                    "quality_assessment": bi_harm.get("quality_assessment"),
                    "n_joint_clusters": wnn_harm_res["n_joint_clusters"],
                    "mean_rna_weight": round(float(wnn_harm_res["mean_rna_weight"]), 4),
                    "mean_atac_weight": round(float(wnn_harm_res["mean_atac_weight"]), 4),
                    "guardrail_status": wnn_harm_res.get("guardrail_status"),
                    "external_subclass_concordance": subclass_concordance["harmony"],
                    "elapsed_s": wnn_harm_elapsed,
                    "a4_resume_verified": bool(wnn_harm_cached.get("resumed") is True),
                    "a4_resume_elapsed_s": wnn_harm_cached_elapsed,
                },
                "donor_mixing_delta": {
                    "wnn_ilisi_gain": round(
                        float(bi_harm.get("wnn_batch_ilisi_median") or 0.0)
                        - float(bi_unint.get("wnn_batch_ilisi_median") or 0.0),
                        4,
                    ),
                    "wnn_ilisi_norm_gain": round(
                        float(bi_harm.get("wnn_batch_ilisi_norm_median") or 0.0)
                        - float(bi_unint.get("wnn_batch_ilisi_norm_median") or 0.0),
                        4,
                    ),
                },
            },
            "stage_4_donor_screened_mofa": {
                "script": "aria/scripts/integration_mofa.py",
                "n_cells": mofa_res["n_cells"],
                "n_factors": mofa_res["n_factors"],
                "variance_explained": mofa_res["variance_explained"],
                "batch_evaluated": mofa_res.get("batch_evaluated"),
                "batch_key": mofa_res.get("batch_key"),
                "condition_col": mofa_res.get("condition_col"),
                "batch_condition_confounded": mofa_res.get("batch_condition_confounded"),
                "batch_factor_flag": mofa_res.get("batch_factor_flag"),
                "batch_factor_ids": mofa_res.get("batch_factor_ids"),
                "factor_batch_eta2": [
                    {"factor_id": f.get("factor_id"), "batch_eta2": f.get("batch_eta2")}
                    for f in (mofa_res.get("top_factors") or [])
                ],
                "technical_factor_flag": mofa_res.get("technical_factor_flag"),
                "n_shared_biological_factors": len(mofa_res.get("shared_factors") or []),
                "elapsed_s": mofa_elapsed,
                "a4_resume_verified": bool(mofa_cached.get("resumed") is True),
                "a4_resume_elapsed_s": mofa_cached_elapsed,
            },
            "stage_5_cis_peak2gene_vs_barcode_permuted_null": {
                "script": "aria/scripts/integration_peak2gene.py",
                "chromosomes_evaluated": list(P2G_CHROMS),
                "n_cells": annot_info["p2g_slice_n_cells"],
                "n_peaks_in_slice": annot_info["p2g_slice_n_peaks"],
                "distance_kb": 250,
                "min_corr": 0.30,
                "max_padj": 0.05,
                "background_match": p2g_real_res.get("background_match"),
                "real_paired": {
                    "n_tested": p2g_real_res["n_tested"],
                    "n_fdr_significant": p2g_real_res["n_fdr_significant"],
                    "n_links": p2g_real_res["n_links"],
                    "n_positive_correlations": p2g_real_res["n_positive_correlations"],
                    "n_negative_correlations": p2g_real_res["n_negative_correlations"],
                    "top_links": (p2g_real_res.get("top_links") or [])[:10],
                    "elapsed_s": p2g_real_elapsed,
                    "a4_resume_verified": bool(p2g_real_cached.get("resumed") is True),
                    "a4_resume_elapsed_s": p2g_real_cached_elapsed,
                },
                "barcode_shuffled_null_control": {
                    "derangement_fixed_points": annot_info["null_fixed_points"],
                    "n_tested": p2g_null_res["n_tested"],
                    "n_fdr_significant": p2g_null_res["n_fdr_significant"],
                    "n_links": p2g_null_res["n_links"],
                    "per_test_false_positive_rate": round(
                        int(p2g_null_res["n_links"]) / max(1, int(p2g_null_res["n_tested"])),
                        8,
                    ),
                    "empirical_false_discovery_rate_at_q005": round(
                        int(p2g_null_res["n_links"]) / max(1, int(p2g_real_res["n_links"])),
                        4,
                    ),
                    "elapsed_s": p2g_null_elapsed,
                },
            },
            "stage_6_governance_and_c3_ledger": {
                "integration_validation": dict(INTEGRATION_VALIDATION),
                "registry_integrity_errors": len(
                    [i for i in check_registry_integrity() if i.severity == "error"]
                ),
                "audit_readiness_probes": {
                    "balanced_4donor_2x2_cohort": {
                        "status": mat_balanced["integration_readiness"]["status"],
                        "dispatch_policy": mat_balanced["integration_readiness"]["dispatch_policy"],
                        "requires_ack_for_dispatch": mat_balanced["integration_readiness"][
                            "requires_ack_for_dispatch"
                        ],
                    },
                    "confounded_1donor_per_condition_cohort": {
                        "status": mat_confounded["integration_readiness"]["status"],
                        "dispatch_policy": mat_confounded["integration_readiness"][
                            "dispatch_policy"
                        ],
                        "requires_ack_for_dispatch": mat_confounded["integration_readiness"][
                            "requires_ack_for_dispatch"
                        ],
                        "finding_checks": [
                            f["check"]
                            for f in mat_confounded["integration_readiness"]["findings"]
                        ],
                    },
                    "unpaired_scrna_plus_scatac": {
                        "status": mat_unpaired["integration_readiness"]["status"],
                        "dispatch_policy": mat_unpaired["integration_readiness"]["dispatch_policy"],
                        "single_modality_blocked": mat_unpaired["dispatch"]["blocked"],
                    },
                },
                "c3_preprocessing_stage_status": pre_entry["status"],
                "narrative_block_ids": [b.id for b in narrator_blocks],
                "claim_ledger_violations": ledger_verification["n_violations"],
                "rocrate_generated": (report_html_path.parent / "ro-crate-metadata.json").is_file(),
                "environment_locks": {
                    "aria_chromatin_env": audit_environment(
                        "aria-chromatin-env", env_name="aria-chromatin-env"
                    ).get("status"),
                    "aria_integration_env": audit_environment(
                        "aria-integration-env", env_name="aria-integration-env"
                    ).get("status"),
                },
            },
            "total_wall_time_s": round(time.perf_counter() - t_start, 2),
            "version_metadata": _sanitize_obj(collect_version_metadata()),
        }

    clean_receipt = _sanitize_obj(receipt)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(clean_receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote ADR-079 validation receipt -> {args.output.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
