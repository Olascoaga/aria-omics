"""
ARIA WNN Anti-Degradation Diagnostics Script (IPC)
--------------------------------------------------
Evaluates a paired WNN `.h5mu` object produced by ``integration_wnn.py``
against the 5-axis "Lens, not replacement" guardrail (V48 Slice I3, ADR-066).

Executed inside ``aria-chromatin-env`` (``stack="chromatin"``) via
``EnvironmentManager.run_in_stack``.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from aria.scripts._base import run_script
from aria.utils.mudata_io import read_h5mu_paired
from aria.utils.wnn_diagnostics import (
    compute_wnn_anti_degradation_diagnostics,
    sanitize_for_hdf5_uns,
)


_DOUBLET_SCORE_CANDIDATES: tuple[str, ...] = (
    "doublet_score",
    "scrublet_score",
    "scdblfinder_score",
)


def _resolve_obs_column(
    col_name: str,
    sources: list[Any],
) -> list[str] | None:
    """Lookup an obs column across MuData global obs and modality obs tables."""
    for df in sources:
        if df is not None and col_name in df.columns:
            return [str(v) for v in df[col_name].tolist()]
    return None


def _resolve_numeric_obs(
    candidates: tuple[str, ...],
    sources: list[Any],
):
    """Return the first finite numeric obs array matching candidate names."""
    import numpy as np

    for df in sources:
        if df is None:
            continue
        for col in candidates:
            if col in df.columns:
                try:
                    arr = np.asarray(df[col].values, dtype=np.float64)
                except (TypeError, ValueError):
                    continue
                if arr.ndim == 1 and np.all(np.isfinite(arr)):
                    return arr, col
    return None, None


def integration_wnn_diagnostics(params: dict) -> dict[str, Any]:
    """Run the 5-axis WNN anti-degradation guardrail on a paired `.h5mu`."""
    h5mu_raw = params.get("h5mu_path") or params.get("data_path")
    output_dir = params.get("output_dir")
    output_json = params.get("output_json")
    update_h5mu = bool(params.get("update_h5mu", False))
    audit_group_col = params.get("audit_group_col")
    seed = int(params.get("seed", 0))
    k = int(params.get("k", 20))
    rare_fraction_threshold = float(params.get("rare_fraction_threshold", 0.05))
    min_rare_cells = int(params.get("min_rare_cells", 3))
    min_evaluable_rare_cells = int(params.get("min_evaluable_rare_cells", 5))
    n_folds = int(params.get("n_folds", 5))

    if not h5mu_raw:
        return {
            "status": "error",
            "error_type": "MissingH5MU",
            "details": "Parameter 'h5mu_path' is required.",
        }

    h5mu_path = Path(str(h5mu_raw))
    if not h5mu_path.exists():
        return {
            "status": "error",
            "error_type": "MissingH5MU",
            "details": f"Paired WNN .h5mu file does not exist: {h5mu_path}",
        }

    try:
        import numpy as np

        loaded = read_h5mu_paired(h5mu_path, require_raw_counts=False)
        if loaded.get("status") != "success":
            return {
                "status": "error",
                "error_type": str(loaded.get("error_type", "UnpairedMuData")),
                "details": str(
                    loaded.get("details", "Failed paired .h5mu validation.")
                ),
            }

        mdata = loaded["mdata"]
        rna = loaded["rna_adata"]
        atac = loaded["atac_adata"]
        if not loaded.get("exact_same_order", False) or int(
            loaded.get("n_rna_only", 0)
        ) > 0 or int(loaded.get("n_atac_only", 0)) > 0:
            return {
                "status": "error",
                "error_type": "UnalignedModalities",
                "details": (
                    "WNN .h5mu modalities must share identical barcodes in the "
                    "exact same cell order."
                ),
            }

        rna_emb = (
            rna.obsm["X_pca_l2"]
            if "X_pca_l2" in rna.obsm
            else (rna.obsm["X_pca"] if "X_pca" in rna.obsm else None)
        )
        atac_emb = (
            atac.obsm["X_lsi_l2"]
            if "X_lsi_l2" in atac.obsm
            else (atac.obsm["X_lsi"] if "X_lsi" in atac.obsm else None)
        )
        if rna_emb is None or atac_emb is None:
            return {
                "status": "error",
                "error_type": "MissingWNNArtifacts",
                "details": (
                    "WNN diagnostics require RNA ('X_pca_l2' or 'X_pca') and "
                    "ATAC ('X_lsi_l2' or 'X_lsi') latent embeddings."
                ),
            }

        if "wnn_distances" in mdata.obsp:
            wnn_graph = mdata.obsp["wnn_distances"]
            is_conn = False
        elif "wnn_connectivities" in mdata.obsp:
            wnn_graph = mdata.obsp["wnn_connectivities"]
            is_conn = True
        else:
            return {
                "status": "error",
                "error_type": "MissingWNNArtifacts",
                "details": (
                    "MuData is missing 'wnn_distances' and 'wnn_connectivities' "
                    "in mdata.obsp."
                ),
            }

        obs_sources = [mdata.obs, rna.obs, atac.obs]
        wnn_clusters = _resolve_obs_column("leiden_wnn", obs_sources)
        rna_clusters = _resolve_obs_column("leiden_rna", [rna.obs, mdata.obs])
        atac_clusters = _resolve_obs_column("leiden_atac", [atac.obs, mdata.obs])
        if wnn_clusters is None or rna_clusters is None or atac_clusters is None:
            return {
                "status": "error",
                "error_type": "MissingWNNArtifacts",
                "details": (
                    "MuData is missing one or more required cluster columns: "
                    "'leiden_wnn', 'leiden_rna', 'leiden_atac'."
                ),
            }

        rna_w = None
        atac_w = None
        for r_col, a_col in (
            ("rna_weight", "atac_weight"),
            ("rna:mod_weight", "atac:mod_weight"),
        ):
            if r_col in mdata.obs.columns and a_col in mdata.obs.columns:
                rna_w = np.asarray(mdata.obs[r_col].values, dtype=np.float64)
                atac_w = np.asarray(mdata.obs[a_col].values, dtype=np.float64)
                break

        audit_labels: list[str] | None = None
        if audit_group_col:
            audit_labels = _resolve_obs_column(str(audit_group_col), obs_sources)
            if audit_labels is None:
                return {
                    "status": "error",
                    "error_type": "MissingAuditColumn",
                    "details": (
                        f"Requested audit_group_col='{audit_group_col}' not "
                        f"found in MuData obs tables."
                    ),
                }

        dbl_scores, dbl_col = _resolve_numeric_obs(
            _DOUBLET_SCORE_CANDIDATES, obs_sources
        )

        uns_wnn = mdata.uns.get("aria_wnn", {})
        if "k" not in params and isinstance(uns_wnn, dict) and "k" in uns_wnn:
            k = int(uns_wnn["k"])

        diag = compute_wnn_anti_degradation_diagnostics(
            rna_embedding=np.asarray(rna_emb, dtype=np.float64),
            atac_embedding=np.asarray(atac_emb, dtype=np.float64),
            wnn_distances=wnn_graph,
            wnn_clusters=wnn_clusters,
            rna_clusters=rna_clusters,
            atac_clusters=atac_clusters,
            rna_weights=rna_w,
            atac_weights=atac_w,
            audit_labels=audit_labels,
            doublet_scores=dbl_scores,
            graph_is_connectivity=is_conn,
            k=k,
            seed=seed,
            rare_fraction_threshold=rare_fraction_threshold,
            min_rare_cells=min_rare_cells,
            min_evaluable_rare_cells=min_evaluable_rare_cells,
            n_folds=n_folds,
        )

        if update_h5mu:
            mdata.uns["aria_wnn_diagnostics"] = sanitize_for_hdf5_uns(diag)
            tmp_h5mu = h5mu_path.parent / f".{h5mu_path.stem}.diag.{os.getpid()}.h5mu"
            try:
                mdata.write_h5mu(str(tmp_h5mu))
                os.replace(str(tmp_h5mu), str(h5mu_path))
            finally:
                if tmp_h5mu.exists():
                    try:
                        tmp_h5mu.unlink()
                    except OSError:
                        pass

        written_report_path: str | None = None
        if output_json or output_dir:
            target_json = (
                Path(str(output_json))
                if output_json
                else Path(str(output_dir)) / "wnn_anti_degradation_diagnostics.json"
            )
            target_json.parent.mkdir(parents=True, exist_ok=True)
            tmp_json = target_json.parent / f".{target_json.name}.tmp.{os.getpid()}"
            try:
                tmp_json.write_text(
                    json.dumps(diag, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
                os.replace(str(tmp_json), str(target_json))
                written_report_path = str(target_json)
            finally:
                if tmp_json.exists():
                    try:
                        tmp_json.unlink()
                    except OSError:
                        pass

        return {
            "status": "success",
            "guardrail_status": diag["guardrail_status"],
            "n_cells": diag["n_cells"],
            "seed": diag["seed"],
            "issue_codes": diag["issue_codes"],
            "issues": diag["issues"],
            "doublet_score_source": dbl_col,
            "audit_group_col": str(audit_group_col) if audit_group_col else None,
            "diagnostics": diag,
            "report_path": written_report_path,
        }

    except ImportError as exc:
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": f"Chromatin/MuData dependencies required: {exc}",
        }
    except Exception as exc:
        return {
            "status": "error",
            "error_type": "WNNDiagnosticsFailed",
            "details": f"{type(exc).__name__}: {str(exc)[:500]}",
        }


if __name__ == "__main__":
    run_script(integration_wnn_diagnostics)
