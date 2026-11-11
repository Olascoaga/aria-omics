#!/usr/bin/env python3
"""V48 Slice I4 — End-to-end MOFA+ external concordance & determinism runner.

Executes on real paired 10x Multiome data:
  1. `aria/scripts/integration_paired_qc.py` (`aria-chromatin-env`) for joint
     RNA+ATAC QC and optional bijective reference label join.
  2. Repeated seeded executions of `aria/scripts/integration_mofa.py`
     (`aria-integration-env`):
     - `seed=0` x 2 to prove bit/float-exact determinism (`max |delta Z| == 0`).
     - `seed=1` to measure MOFA+'s intrinsic variational convergence ceiling.
     - `max_cells=1200` subsampled execution to validate deterministic
       subsampling + closed-form variational posterior projection onto all cells.
  3. Independent external reference execution via `muon.tl.mofa`
     (`aria/scripts/benchmark_mofa_muon_ref.py` in comparator environment, e.g.
     `oligos` `muon 0.1.9` + `mofapy2 0.7.4`).
  4. Full concordance scoring via `aria/benchmarks/concordance_mofa.py`
     (Hungarian 1-to-1 active-factor matching weighted by reference R^2,
     rotation-invariant canonical subspace correlations, sign-aligned loading
     concordance, and per-view R^2 agreement).

Emits a sanitized public JSON artifact with zero machine-absolute paths
(enforced by `tests/test_public_artifacts_no_absolute_paths.py`).
"""

from __future__ import annotations

import argparse
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

from aria.benchmarks.concordance_mofa import (
    score_mofa_concordance,
    score_mofa_seed_stability,
)


def _resolve_conda_bin() -> str:
    found = shutil.which("conda")
    if found:
        return found
    candidate = Path(sys.executable).resolve().parents[2] / "bin" / "conda"
    if candidate.is_file():
        return str(candidate)
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
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(out_json.read_text(encoding="utf-8"))


def _read_scores_csv(csv_path: Path) -> tuple[list[str], np.ndarray]:
    lines = csv_path.read_text(encoding="utf-8").strip().splitlines()
    barcodes: list[str] = []
    rows: list[list[float]] = []
    for line in lines[1:]:
        parts = line.split(",")
        barcodes.append(parts[0])
        rows.append([float(x) for x in parts[1:]])
    return barcodes, np.asarray(rows, dtype=np.float64)


def _extract_aria_loadings_from_h5mu(
    env: str,
    h5mu_path: Path,
    work_dir: Path,
    tag: str,
) -> dict[str, tuple[list[str], np.ndarray]]:
    out_npz = work_dir / f"{tag}_aria_loadings.npz"
    code = (
        "import mudata as md, numpy as np, sys; "
        "m = md.read_h5mu(sys.argv[1]); "
        "mods = list(m.mod.keys()); "
        "r = next(k for k in mods if 'rna' in k.lower()); "
        "a = next(k for k in mods if 'atac' in k.lower() or 'peak' in k.lower()); "
        "rm = m.mod[r]; am = m.mod[a]; "
        "ru = rm.var['mofa_feature_used'].values.astype(bool); "
        "au = am.var['mofa_feature_used'].values.astype(bool); "
        "np.savez_compressed(sys.argv[2], "
        "rna_mod=r, atac_mod=a, "
        "rna_features=np.asarray(rm.var_names[ru], dtype=object), "
        "rna_weights=np.asarray(rm.varm['LFs'][ru, :], dtype=np.float64), "
        "atac_features=np.asarray(am.var_names[au], dtype=object), "
        "atac_weights=np.asarray(am.varm['LFs'][au, :], dtype=np.float64))"
    )
    proc = subprocess.run(
        [
            _resolve_conda_bin(),
            "run",
            "-n",
            env,
            "python",
            "-c",
            code,
            str(h5mu_path),
            str(out_npz),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not out_npz.is_file():
        raise RuntimeError(f"Failed extracting ARIA LFs: {proc.stderr}")
    dat = np.load(out_npz, allow_pickle=True)
    rna_mod = str(dat["rna_mod"].item())
    atac_mod = str(dat["atac_mod"].item())
    return {
        rna_mod: (
            [str(x) for x in dat["rna_features"].tolist()],
            np.asarray(dat["rna_weights"], dtype=np.float64),
        ),
        atac_mod: (
            [str(x) for x in dat["atac_features"].tolist()],
            np.asarray(dat["atac_weights"], dtype=np.float64),
        ),
    }


def _sanitize_public_value(val: Any) -> Any:
    """Strip machine-absolute filesystem prefixes from public JSON artifacts."""
    if isinstance(val, dict):
        return {str(k): _sanitize_public_value(v) for k, v in val.items()}
    if isinstance(val, list):
        return [_sanitize_public_value(v) for v in val]
    if isinstance(val, str):
        p = Path(val)
        if val.startswith(
            ("/home/", "/Users/", "/root/", "/tmp/", "/mnt/", "/private/", "/data/")
        ):
            try:
                return str(p.resolve().relative_to(ROOT.resolve()))
            except Exception:
                return p.name
        return val
    return val


def _r2_per_view_from_aria(res: dict[str, Any]) -> dict[str, list[float]]:
    views = list((res.get("variance_explained") or {}).keys())
    out: dict[str, list[float]] = {v: [] for v in views}
    for f in res.get("top_factors", []):
        r2_map = f.get("r2_by_view") or {}
        for v in views:
            out[v].append(float(r2_map.get(v, 0.0)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5mu", required=True, help="Input paired .h5mu.")
    ap.add_argument("--reference-h5ad", default="")
    ap.add_argument("--audit-group-col", default="")
    ap.add_argument("--sample-id", default="")
    ap.add_argument("--reference-sample-col", default="")
    ap.add_argument("--chromatin-env", default="aria-chromatin-env")
    ap.add_argument("--integration-env", default="aria-integration-env")
    ap.add_argument("--reference-env", default="oligos")
    ap.add_argument("--n-factors", type=int, default=10)
    ap.add_argument("--n-top-genes", type=int, default=3000)
    ap.add_argument("--n-top-peaks", type=int, default=5000)
    ap.add_argument("--subsample-cells", type=int, default=1200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--skip-paired-qc", action="store_true")
    ap.add_argument("--work-dir", default="")
    ap.add_argument(
        "--output-json",
        default=str(
            ROOT
            / "docs"
            / "benchmark_results"
            / "multiome_concordance"
            / "v48_i4_mofa_concordance_hc11.json"
        ),
    )
    args = ap.parse_args()

    cleanup_tmp: tempfile.TemporaryDirectory[str] | None = None
    if args.work_dir:
        work_dir = Path(args.work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
    else:
        cleanup_tmp = tempfile.TemporaryDirectory(prefix="v48_i4_mofa_")
        work_dir = Path(cleanup_tmp.name)

    try:
        input_h5mu = Path(args.h5mu)
        qc_summary: dict[str, Any] | None = None
        if args.skip_paired_qc:
            mofa_input_h5mu = input_h5mu
        else:
            qc_dir = work_dir / "i1_qc"
            qc_params: dict[str, Any] = {
                "h5mu_path": str(input_h5mu),
                "output_dir": str(qc_dir),
            }
            if args.sample_id:
                qc_params["sample_id"] = str(args.sample_id)
            if args.reference_h5ad and args.audit_group_col:
                qc_params["reference_h5ad"] = str(args.reference_h5ad)
                qc_params["reference_label_col"] = str(args.audit_group_col)
                if args.reference_sample_col:
                    qc_params["reference_sample_col"] = str(args.reference_sample_col)
            qc_res = _run_ipc_script(
                args.chromatin_env,
                "aria/scripts/integration_paired_qc.py",
                qc_params,
                work_dir=work_dir,
                tag="i1_qc",
            )
            if qc_res.get("status") != "success":
                raise RuntimeError(f"I1 paired QC failed: {qc_res}")
            mofa_input_h5mu = Path(str(qc_res["output_path"]))
            lj = dict(qc_res.get("label_join") or {})
            if isinstance(lj.get("label_counts"), dict):
                raw_counts = lj.pop("label_counts")
                lj["n_unique_labels"] = len(raw_counts)
                lj["label_size_distribution"] = sorted(
                    int(v) for v in raw_counts.values()
                )
            qc_summary = {
                "n_cells_before": qc_res["n_cells_before"],
                "n_cells_after": qc_res["n_cells_after"],
                "n_genes": qc_res["n_genes"],
                "n_peaks": qc_res["n_peaks"],
                "label_join": lj,
            }

        base_mofa_params = {
            "h5mu_path": str(mofa_input_h5mu),
            "n_factors": int(args.n_factors),
            "n_top_genes": int(args.n_top_genes),
            "n_top_peaks": int(args.n_top_peaks),
            "max_cells": 100000,
            "scale_views": True,
            "write_h5mu": True,
        }

        # 1. ARIA seed=0 (run 1 & run 2 -> exact determinism)
        aria_r1 = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_mofa.py",
            {
                **base_mofa_params,
                "seed": int(args.seed),
                "output_dir": str(work_dir / "aria_seed0_r1"),
            },
            work_dir=work_dir,
            tag="aria_seed0_r1",
        )
        aria_r2 = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_mofa.py",
            {
                **base_mofa_params,
                "seed": int(args.seed),
                "output_dir": str(work_dir / "aria_seed0_r2"),
            },
            work_dir=work_dir,
            tag="aria_seed0_r2",
        )
        # 2. ARIA seed=1 (cross-seed variational convergence ceiling)
        aria_s1 = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_mofa.py",
            {
                **base_mofa_params,
                "seed": int(args.seed) + 1,
                "write_h5mu": False,
                "output_dir": str(work_dir / "aria_seed1"),
            },
            work_dir=work_dir,
            tag="aria_seed1",
        )
        # 3. ARIA subsampled run + posterior projection validation
        sub_params = {
            **base_mofa_params,
            "seed": int(args.seed),
            "max_cells": int(args.subsample_cells),
            "output_dir": str(work_dir / "aria_subsampled"),
        }
        if args.audit_group_col:
            sub_params["subsample_strata_col"] = str(args.audit_group_col)
        aria_sub = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_mofa.py",
            sub_params,
            work_dir=work_dir,
            tag="aria_subsampled",
        )
        for tag_name, r_obj in (
            ("aria_seed0_r1", aria_r1),
            ("aria_seed0_r2", aria_r2),
            ("aria_seed1", aria_s1),
            ("aria_subsampled", aria_sub),
        ):
            if r_obj.get("status") != "success":
                raise RuntimeError(f"{tag_name} failed: {r_obj}")

        # 4. External reference (`muon.tl.mofa` in comparator environment)
        muon_ref = _run_ipc_script(
            args.reference_env,
            "aria/scripts/benchmark_mofa_muon_ref.py",
            {
                "h5mu_path": str(mofa_input_h5mu),
                "output_dir": str(work_dir / "muon_ref"),
                "n_factors": int(args.n_factors),
                "n_top_genes": int(args.n_top_genes),
                "n_top_peaks": int(args.n_top_peaks),
                "scale_views": True,
                "seed": int(args.seed),
            },
            work_dir=work_dir,
            tag="muon_ref",
        )
        if muon_ref.get("status") != "success":
            raise RuntimeError(f"muon.tl.mofa reference failed: {muon_ref}")

        bcs_r1, z_r1 = _read_scores_csv(Path(aria_r1["factor_scores"]))
        _, z_r2 = _read_scores_csv(Path(aria_r2["factor_scores"]))
        _, z_s1 = _read_scores_csv(Path(aria_s1["factor_scores"]))
        bcs_sub, z_sub = _read_scores_csv(Path(aria_sub["factor_scores"]))
        bcs_ref, z_ref = _read_scores_csv(Path(muon_ref["scores_csv"]))

        aria_loadings = _extract_aria_loadings_from_h5mu(
            args.integration_env,
            Path(aria_r1["h5mu_output_path"]),
            work_dir=work_dir,
            tag="r1",
        )
        ref_npz = np.load(Path(muon_ref["loadings_npz"]), allow_pickle=True)
        ref_loadings = {
            str(muon_ref["rna_modality"]): (
                [str(x) for x in ref_npz["rna_features"].tolist()],
                np.asarray(ref_npz["rna_weights"], dtype=np.float64),
            ),
            str(muon_ref["atac_modality"]): (
                [str(x) for x in ref_npz["atac_features"].tolist()],
                np.asarray(ref_npz["atac_weights"], dtype=np.float64),
            ),
        }

        same_seed_determinism = score_mofa_seed_stability(
            [z_r1, z_r2],
            [aria_r1["variance_explained"], aria_r2["variance_explained"]],
        )
        cross_seed_ceiling = score_mofa_seed_stability(
            [z_r1, z_s1],
            [aria_r1["variance_explained"], aria_s1["variance_explained"]],
        )
        external_concordance = score_mofa_concordance(
            aria_barcodes=bcs_r1,
            aria_scores=z_r1,
            aria_r2_total=aria_r1["variance_explained"],
            aria_r2_per_view=_r2_per_view_from_aria(aria_r1),
            aria_loadings=aria_loadings,
            ref_barcodes=bcs_ref,
            ref_scores=z_ref,
            ref_r2_total=muon_ref["r2_total"],
            ref_r2_per_view=muon_ref["r2_per_view"],
            ref_loadings=ref_loadings,
            min_shared_features=500,
            top_k=50,
            reference_name="muon.tl.mofa",
        )
        subsampled_vs_full = score_mofa_concordance(
            aria_barcodes=bcs_sub,
            aria_scores=z_sub,
            aria_r2_total=aria_sub["variance_explained"],
            aria_r2_per_view=_r2_per_view_from_aria(aria_sub),
            ref_barcodes=bcs_r1,
            ref_scores=z_r1,
            ref_r2_total=aria_r1["variance_explained"],
            ref_r2_per_view=_r2_per_view_from_aria(aria_r1),
            reference_name="aria_full_cohort",
        )

        # Redact raw feature names in public summary while keeping statistical metrics.
        factors_public = []
        for f in aria_r1.get("top_factors", []):
            factors_public.append(
                {
                    "factor_id": f["factor_id"],
                    "active": f["active"],
                    "category": f["category"],
                    "active_views": f["active_views"],
                    "dominant_modality": f["dominant_modality"],
                    "sharing_balance": f["sharing_balance"],
                    "r2_by_view": f["r2_by_view"],
                    "max_depth_abs_r": f["max_depth_abs_r"],
                    "technical_depth_flag": f["technical_depth_flag"],
                    "cell_cycle_flag": f["cell_cycle_flag"],
                    "excluded_from_biological_interpretation": f[
                        "excluded_from_biological_interpretation"
                    ],
                    "exclusion_reasons": f["exclusion_reasons"],
                }
            )

        manifest = {
            "slice": "V48_I4_mofa_interpretable_factors",
            "dataset_id": input_h5mu.stem,
            "parameters": {
                "seed": int(args.seed),
                "n_factors_requested": int(args.n_factors),
                "n_top_genes": int(args.n_top_genes),
                "n_top_peaks": int(args.n_top_peaks),
                "scale_views": True,
                "subsample_validation_max_cells": int(args.subsample_cells),
                "integration_env": str(args.integration_env),
                "reference_comparator": "muon.tl.mofa",
            },
            "paired_qc_summary": qc_summary,
            "aria_mofa_summary": {
                "n_cells": aria_r1["n_cells"],
                "n_factors": aria_r1["n_factors"],
                "n_active_factors": aria_r1["n_active_factors"],
                "n_shared_factors": len(aria_r1["shared_factors"]),
                "modality_specific_counts": {
                    k: len(v)
                    for k, v in aria_r1["modality_specific_factors"].items()
                },
                "n_excluded_factors": len(aria_r1["excluded_factors"]),
                "technical_factor_flag": aria_r1["technical_factor_flag"],
                "technical_factor_id": aria_r1["technical_factor_id"],
                "cell_cycle_evaluated": aria_r1["cell_cycle_evaluated"],
                "cell_cycle_factor_flag": aria_r1["cell_cycle_factor_flag"],
                "variance_explained": aria_r1["variance_explained"],
                "variance_explained_additive": aria_r1[
                    "variance_explained_additive"
                ],
                "metrics_not_computed": aria_r1["metrics_not_computed"],
                "factors": factors_public,
            },
            "determinism_seed0": same_seed_determinism,
            "variational_cross_seed_ceiling": cross_seed_ceiling,
            "external_concordance_vs_muon_mofa": external_concordance,
            "subsampling_and_posterior_projection": {
                "subsampling_metadata": aria_sub["subsampling"],
                "all_cells_finite_no_nan": bool(np.all(np.isfinite(z_sub))),
                "concordance_vs_full_cohort": {
                    "ref_r2_weighted_abs_r": subsampled_vs_full["factor_matching"][
                        "ref_r2_weighted_abs_r"
                    ],
                    "mean_matched_abs_r": subsampled_vs_full["factor_matching"][
                        "mean_matched_abs_r"
                    ],
                    "subspace_mean_canonical_corr": subsampled_vs_full[
                        "subspace_concordance"
                    ]["mean_canonical_corr"],
                    "subspace_min_top3_canonical_corr": subsampled_vs_full[
                        "subspace_concordance"
                    ]["min_top3_canonical_corr"],
                },
            },
        }

        sanitized = _sanitize_public_value(manifest)
        out_path = Path(args.output_json)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            json.dumps(sanitized, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(
            json.dumps(
                {
                    "status": "success",
                    "output_json": str(out_path.name),
                    "n_cells": aria_r1["n_cells"],
                    "n_active_factors": aria_r1["n_active_factors"],
                    "n_shared_factors": len(aria_r1["shared_factors"]),
                    "technical_factor_flag": aria_r1["technical_factor_flag"],
                    "technical_factor_id": aria_r1["technical_factor_id"],
                    "seed0_max_abs_diff": same_seed_determinism[
                        "max_elementwise_abs_diff"
                    ],
                    "external_ref_r2_weighted_abs_r": external_concordance[
                        "factor_matching"
                    ]["ref_r2_weighted_abs_r"],
                    "external_subspace_mean_canonical_corr": external_concordance[
                        "subspace_concordance"
                    ]["mean_canonical_corr"],
                    "subsampled_projection_min_train_r": aria_sub["subsampling"][
                        "projection_min_train_abs_r"
                    ],
                },
                indent=2,
            )
        )
        return 0
    finally:
        if cleanup_tmp is not None:
            cleanup_tmp.cleanup()


if __name__ == "__main__":
    raise SystemExit(main())
