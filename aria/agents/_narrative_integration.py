"""Publication-grade figure generation & table staging for Multi-Omics Integration (V48 I7).

Renders dual raster + vector (``.png`` + ``.svg``) figures inline using ``h5py``
and ``matplotlib`` (Agg backend, ``svg.fonttype='none'`` so axes/labels remain
editable text in vector editors while dense UMAP point clouds are rasterized
inside the SVG):

  1. ``wnn_joint_vs_modality_umap.{png,svg}`` — Joint WNN vs RNA-only vs ATAC-only
     UMAP colored consistently by cell label / joint cluster, with ARI/NMI and
     I3 5-axis guardrail badge.
  2. ``wnn_modality_weights_by_cluster.{png,svg}`` — Per-cluster RNA vs ATAC WNN
     modality weight distributions (with cell-level IQR when ``.h5mu`` is
     available) and 0.15/0.85 dominance guardrail thresholds.
  3. ``mofa_variance_explained.{png,svg}`` — Factor x modality R^2 (%) heatmap
     and bar decomposition explicitly hatching and labeling excluded technical
     sequencing-depth and cell-cycle confounder factors alongside non-additive
     total R^2 vs additive sum.
  4. ``regulatory_triangulation_summary.{png,svg}`` — Top Signac-faithful
     ``peak_to_gene`` cis links (distance vs Pearson r by region class) paired
     with top candidate ``eRegulons`` (``pos_tf_pos_p2g`` vs ``neg_tf_pos_p2g``,
     gene/peak counts, RNA/ATAC AUCell concordance, and paralog/depth badges).

No fabrication: figures render strictly from real pipeline outputs and skip
honestly when upstream data is absent. Zero hardcoded biological gene/TF panels.
"""

from __future__ import annotations

import csv
import json
import logging
import math
import shutil
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

log = logging.getLogger(__name__)

_PALETTE_20 = [
    "#0f766e", "#2563eb", "#d97706", "#dc2626", "#7c3aed",
    "#059669", "#db2777", "#4f46e5", "#ca8a04", "#0891b2",
    "#65a30d", "#9333ea", "#ea580c", "#0284c7", "#be123c",
    "#15803d", "#6d28d9", "#b45309", "#0e7490", "#475569",
]

_GUARDRAIL_COLORS = {
    "pass": "#059669",
    "warn": "#d97706",
    "degraded": "#dc2626",
    "skipped": "#64748b",
}


def _decode_str_array(raw: np.ndarray) -> list[str]:
    out: list[str] = []
    for item in raw.tolist():
        if isinstance(item, (bytes, bytearray, np.bytes_)):
            out.append(item.decode("utf-8", errors="replace"))
        elif item is None:
            out.append("NA")
        else:
            s = str(item)
            out.append("NA" if s in {"", "nan", "None", "<NA>"} else s)
    return out


def _read_h5_index(obs_group) -> list[str]:
    """Read barcode index from an AnnData/MuData HDF5 ``obs`` group safely."""
    idx_key = obs_group.attrs.get("_index", "_index")
    if isinstance(idx_key, (bytes, bytearray, np.bytes_)):
        idx_key = idx_key.decode("utf-8", errors="replace")
    idx_key = str(idx_key or "_index")
    if idx_key in obs_group:
        ds = obs_group[idx_key]
        return _decode_str_array(np.asarray(ds[()]))
    for fallback in ("_index", "index", "cell_id", "barcode"):
        if fallback in obs_group:
            return _decode_str_array(np.asarray(obs_group[fallback][()]))
    return []


def _decode_h5_obs_column(item) -> np.ndarray | list[str] | None:
    """Decode an HDF5 ``obs`` column (Dataset, categorical Group, or nullable).

    Critical invariant: categorical ``codes < 0`` (AnnData NA code ``-1``) MUST
    map to ``"NA"``, never ``categories[-1]``.
    """
    import h5py

    if isinstance(item, h5py.Group):
        if "codes" in item and "categories" in item:
            codes = np.asarray(item["codes"][()], dtype=np.int64).ravel()
            cats = _decode_str_array(np.asarray(item["categories"][()]))
            n_cats = len(cats)
            decoded: list[str] = []
            for c in codes.tolist():
                if c < 0 or c >= n_cats:
                    decoded.append("NA")
                else:
                    decoded.append(cats[c])
            return decoded
        if "values" in item and "mask" in item:
            vals = np.asarray(item["values"][()])
            mask = np.asarray(item["mask"][()], dtype=bool)
            if vals.dtype.kind in {"f", "i", "u"}:
                arr = vals.astype(np.float64)
                arr[mask] = np.nan
                return arr
            str_vals = _decode_str_array(vals)
            return ["NA" if bool(m) else v for v, m in zip(str_vals, mask.tolist())]
        return None

    if isinstance(item, h5py.Dataset):
        arr = np.asarray(item[()])
        if arr.dtype.kind in {"f"}:
            return arr.astype(np.float64)
        if arr.dtype.kind in {"i", "u"}:
            return arr
        return _decode_str_array(arr)
    return None


def _read_wnn_h5mu(h5mu_path: str | Path) -> dict[str, Any] | None:
    """Read aligned UMAP coordinates, cluster labels, and weights from ``.h5mu``."""
    path = Path(h5mu_path)
    if not path.exists():
        return None
    try:
        import h5py
    except ImportError:
        log.warning("h5py unavailable; skipping .h5mu embedding read")
        return None

    try:
        with h5py.File(path, "r") as f:
            if "obs" not in f or "obsm" not in f:
                return None
            obs_grp = f["obs"]
            barcodes = _read_h5_index(obs_grp)
            if not barcodes:
                return None
            n_cells = len(barcodes)

            obsm_grp = f["obsm"]
            wnn_key = (
                "X_wnn_umap"
                if "X_wnn_umap" in obsm_grp
                else ("X_umap" if "X_umap" in obsm_grp else None)
            )
            if wnn_key is None:
                return None
            wnn_umap = np.asarray(obsm_grp[wnn_key][()], dtype=np.float64)
            if wnn_umap.ndim != 2 or wnn_umap.shape[0] != n_cells or wnn_umap.shape[1] < 2:
                return None

            def _aligned_mod_umap(mod_name: str) -> np.ndarray | None:
                if "mod" not in f or mod_name not in f["mod"]:
                    return None
                mgrp = f["mod"][mod_name]
                if "obsm" not in mgrp or "X_umap" not in mgrp["obsm"]:
                    return None
                raw_umap = np.asarray(mgrp["obsm"]["X_umap"][()], dtype=np.float64)
                if raw_umap.ndim != 2 or raw_umap.shape[1] < 2:
                    return None
                mod_barcodes = _read_h5_index(mgrp["obs"]) if "obs" in mgrp else []
                if not mod_barcodes or len(mod_barcodes) != raw_umap.shape[0]:
                    return None
                pos_map = {bc: i for i, bc in enumerate(mod_barcodes)}
                aligned = np.full((n_cells, 2), np.nan, dtype=np.float64)
                for out_i, bc in enumerate(barcodes):
                    src_i = pos_map.get(bc)
                    if src_i is not None:
                        aligned[out_i] = raw_umap[src_i, :2]
                if np.all(np.isnan(aligned)):
                    return None
                return aligned

            rna_umap = _aligned_mod_umap("rna")
            atac_umap = _aligned_mod_umap("atac")

            obs_cols: dict[str, Any] = {}
            for col_name in obs_grp.keys():
                if str(col_name).startswith("_"):
                    continue
                decoded = _decode_h5_obs_column(obs_grp[col_name])
                if decoded is not None and len(decoded) == n_cells:
                    obs_cols[str(col_name)] = decoded

            return {
                "barcodes": barcodes,
                "wnn_umap": wnn_umap[:, :2],
                "rna_umap": rna_umap,
                "atac_umap": atac_umap,
                "obs": obs_cols,
            }
    except Exception as exc:
        log.warning("Failed reading WNN .h5mu %s: %s", path, exc)
        return None


def _save_dual(fig: plt.Figure, png_path: Path, dpi: int = 220) -> tuple[str, str] | None:
    """Save both ``.png`` and vector ``.svg`` alongside each other."""
    png_path = Path(png_path)
    svg_path = png_path.with_suffix(".svg")
    png_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with plt.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
            fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
            fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    finally:
        plt.close(fig)
    if (
        png_path.exists()
        and svg_path.exists()
        and png_path.stat().st_size > 0
        and svg_path.stat().st_size > 0
    ):
        return str(png_path), str(svg_path)
    return None


def _choose_label_column(obs: dict[str, Any]) -> tuple[str, list[str]]:
    """Select cell-state label column deterministically from ``obs``."""
    preferred = (
        "subclass", "cell_type", "celltype", "cell_label",
        "annotation", "cluster_label",
    )
    skip_prefixes = ("rna:", "atac:", "n_genes", "total_counts", "pct_", "leiden_")

    def _valid_categorical(vals: Any) -> list[str] | None:
        if isinstance(vals, np.ndarray) and vals.dtype.kind == "f":
            return None
        labels = [str(v) for v in (vals.tolist() if isinstance(vals, np.ndarray) else vals)]
        non_na = {x for x in labels if x not in {"NA", "nan", "None", ""}}
        if 2 <= len(non_na) <= 30:
            return labels
        return None

    for key in preferred:
        if key in obs:
            cand = _valid_categorical(obs[key])
            if cand is not None:
                return key, cand

    for key in sorted(obs.keys()):
        low = key.lower()
        if low.startswith(skip_prefixes) or low.endswith("_weight"):
            continue
        cand = _valid_categorical(obs[key])
        if cand is not None:
            return key, cand

    if "leiden_wnn" in obs:
        vals = obs["leiden_wnn"]
        labels = [str(v) for v in (vals.tolist() if isinstance(vals, np.ndarray) else vals)]
        return "leiden_wnn", labels

    n = len(next(iter(obs.values()))) if obs else 0
    return "cells", ["all"] * n


def render_wnn_umap_comparison(
    wnn_findings: dict[str, Any] | None,
    output_dir: Path,
) -> tuple[str, str] | None:
    """Render 3-panel Joint WNN vs RNA-only vs ATAC-only UMAP comparison."""
    if not isinstance(wnn_findings, dict):
        return None
    if str(wnn_findings.get("status") or "").lower() not in {"success", "ok", "done"}:
        return None
    h5mu_path = wnn_findings.get("output_path")
    if not h5mu_path:
        return None

    bundle = _read_wnn_h5mu(h5mu_path)
    if bundle is None:
        return None
    wnn_umap = bundle["wnn_umap"]
    rna_umap = bundle["rna_umap"]
    atac_umap = bundle["atac_umap"]
    if rna_umap is None or atac_umap is None:
        return None

    label_col, labels = _choose_label_column(bundle["obs"])
    unique_groups = sorted(
        {g for g in labels if g != "NA"},
        key=lambda x: (0, int(x)) if x.isdigit() else (1, x),
    )
    if "NA" in set(labels):
        unique_groups.append("NA")

    color_map = {
        grp: ("#cbd5e1" if grp == "NA" else _PALETTE_20[i % len(_PALETTE_20)])
        for i, grp in enumerate(unique_groups)
    }

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.9), constrained_layout=True)
    panels = [
        ("Joint WNN UMAP", wnn_umap),
        ("RNA-only UMAP (PCA L2)", rna_umap),
        ("ATAC-only UMAP (TF-IDF LSI L2)", atac_umap),
    ]
    n_cells = len(labels)
    pt_size = 9.0 if n_cells <= 3500 else (5.0 if n_cells <= 12000 else 2.5)

    labels_arr = np.asarray(labels, dtype=object)
    for ax, (title, coords) in zip(axes, panels):
        for grp in unique_groups:
            mask = (labels_arr == grp) & np.isfinite(coords[:, 0]) & np.isfinite(coords[:, 1])
            if not np.any(mask):
                continue
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                s=pt_size,
                c=color_map[grp],
                alpha=0.82 if grp != "NA" else 0.35,
                edgecolors="none",
                label=f"{grp} ({int(mask.sum())})",
                rasterized=True,
            )
        ax.set_title(title, fontsize=11, fontweight="bold", color="#0f172a")
        ax.set_xlabel("UMAP 1", fontsize=9, color="#475569")
        ax.set_ylabel("UMAP 2", fontsize=9, color="#475569")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color("#cbd5e1")

    handles, leg_labels = axes[0].get_legend_handles_labels()
    if handles:
        ncol = 1 if len(unique_groups) <= 12 else 2
        axes[2].legend(
            handles,
            leg_labels,
            title=f"Group ({label_col})",
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            frameon=True,
            fontsize=7.8,
            title_fontsize=8.5,
            markerscale=2.0,
        )

    ari_rna = wnn_findings.get("wnn_vs_rna_ari")
    ari_atac = wnn_findings.get("wnn_vs_atac_ari")
    nmi_rna = wnn_findings.get("wnn_vs_rna_nmi")
    nmi_atac = wnn_findings.get("wnn_vs_atac_nmi")
    guard = str(wnn_findings.get("guardrail_status") or "skipped")
    n_joint = wnn_findings.get("n_joint_clusters", "?")
    subtitle = (
        f"n={n_cells:,} shared cells | Joint clusters={n_joint} | "
        f"WNN-RNA ARI={ari_rna} (NMI={nmi_rna}) | "
        f"WNN-ATAC ARI={ari_atac} (NMI={nmi_atac}) | "
        f"Guardrail={guard.upper()} (independent UMAP axes; lens, not replacement)"
    )
    fig.suptitle(
        "Paired Multiome Manifold Comparison: Joint WNN vs Single-Modality Embeddings\n"
        + subtitle,
        fontsize=10.5,
        fontweight="bold",
        color="#0f172a",
    )
    return _save_dual(fig, Path(output_dir) / "wnn_joint_vs_modality_umap.png")


def render_wnn_modality_weights(
    wnn_findings: dict[str, Any] | None,
    output_dir: Path,
) -> tuple[str, str] | None:
    """Render per-cluster WNN RNA vs ATAC modality weight profile."""
    if not isinstance(wnn_findings, dict):
        return None
    if str(wnn_findings.get("status") or "").lower() not in {"success", "ok", "done"}:
        return None
    per_cl = wnn_findings.get("per_cluster_modality_weights")
    if not isinstance(per_cl, dict) or not per_cl:
        return None

    # Optional cell-level IQR from .h5mu when available.
    iqr_by_cluster: dict[str, tuple[float, float, float]] = {}
    h5mu_path = wnn_findings.get("output_path")
    if h5mu_path:
        bundle = _read_wnn_h5mu(h5mu_path)
        if bundle is not None:
            obs = bundle["obs"]
            cl_vals = obs.get("leiden_wnn")
            rw_vals = obs.get("rna_weight")
            if cl_vals is not None and isinstance(rw_vals, np.ndarray):
                cl_arr = np.asarray(cl_vals, dtype=object)
                for cl_key in per_cl.keys():
                    mask = (cl_arr == str(cl_key)) & np.isfinite(rw_vals)
                    if np.any(mask):
                        q25, q50, q75 = np.percentile(rw_vals[mask], [25, 50, 75])
                        iqr_by_cluster[str(cl_key)] = (float(q25), float(q50), float(q75))

    rows = []
    for cl_id, meta in per_cl.items():
        if not isinstance(meta, dict):
            continue
        rna_w = meta.get("mean_rna_weight")
        atac_w = meta.get("mean_atac_weight")
        if rna_w is None or atac_w is None:
            continue
        rows.append({
            "cluster": str(cl_id),
            "n_cells": int(meta.get("n_cells") or 0),
            "rna_w": float(rna_w),
            "atac_w": float(atac_w),
            "dominant": str(meta.get("dominant_modality") or "balanced"),
        })
    if not rows:
        return None

    rows.sort(key=lambda r: (r["rna_w"], r["cluster"]))
    n_rows = len(rows)
    fig_h = max(4.0, 0.42 * n_rows + 1.8)
    fig, ax = plt.subplots(figsize=(9.2, fig_h), constrained_layout=True)

    y_pos = np.arange(n_rows)
    rna_vals = np.array([r["rna_w"] for r in rows], dtype=np.float64)
    atac_vals = np.array([r["atac_w"] for r in rows], dtype=np.float64)

    ax.barh(
        y_pos,
        rna_vals,
        color="#2563eb",
        alpha=0.85,
        height=0.66,
        label="Mean RNA weight",
    )
    ax.barh(
        y_pos,
        atac_vals,
        left=rna_vals,
        color="#0d9488",
        alpha=0.85,
        height=0.66,
        label="Mean ATAC weight",
    )

    for idx, r in enumerate(rows):
        iqr = iqr_by_cluster.get(r["cluster"])
        if iqr is not None:
            q25, q50, q75 = iqr
            ax.plot([q25, q75], [idx, idx], color="#0f172a", lw=1.8, zorder=4)
            ax.scatter([q50], [idx], color="white", edgecolors="#0f172a", s=24, zorder=5)

    ax.axvline(0.50, color="#334155", ls="--", lw=1.1, label="Balanced (0.50)")
    ax.axvline(0.15, color="#dc2626", ls=":", lw=1.1, label="Dominance guardrail (0.15 / 0.85)")
    ax.axvline(0.85, color="#dc2626", ls=":", lw=1.1)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [f"Cluster {r['cluster']} (n={r['n_cells']})" for r in rows],
        fontsize=8.8,
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("WNN Modality Weight Fraction (RNA vs ATAC)", fontsize=9.5)
    guard = str(wnn_findings.get("guardrail_status") or "skipped")
    mean_r = wnn_findings.get("mean_rna_weight")
    mean_a = wnn_findings.get("mean_atac_weight")
    ax.set_title(
        f"WNN Per-Cluster Modality Weights (global mean RNA={mean_r}, ATAC={mean_a}; "
        f"guardrail={guard.upper()})",
        fontsize=10.5,
        fontweight="bold",
        color="#0f172a",
    )
    ax.legend(loc="lower right", fontsize=8.2, frameon=True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    return _save_dual(fig, Path(output_dir) / "wnn_modality_weights_by_cluster.png")


def render_mofa_variance_explained(
    mofa_findings: dict[str, Any] | None,
    output_dir: Path,
) -> tuple[str, str] | None:
    """Render MOFA+ variance explained by modality highlighting excluded confounders."""
    if not isinstance(mofa_findings, dict):
        return None
    if str(mofa_findings.get("status") or "").lower() not in {"success", "ok", "done"}:
        return None

    active = mofa_findings.get("top_factors") or []
    excluded = mofa_findings.get("excluded_factors") or []
    factor_map: dict[int, dict[str, Any]] = {}

    def _extract_fac_r2(row: dict[str, Any]) -> tuple[int | None, float, float]:
        raw_id = row.get("factor_id") if row.get("factor_id") is not None else row.get("factor")
        if raw_id is None:
            return None, 0.0, 0.0
        r2_view = row.get("r2_by_view") if isinstance(row.get("r2_by_view"), dict) else {}
        rna_v = row.get("rna_r2") if row.get("rna_r2") is not None else r2_view.get("rna")
        atac_v = row.get("atac_r2") if row.get("atac_r2") is not None else r2_view.get("atac")
        return int(raw_id), float(rna_v or 0.0) * 100.0, float(atac_v or 0.0) * 100.0

    for row in active:
        if not isinstance(row, dict):
            continue
        fid, r_pct, a_pct = _extract_fac_r2(row)
        if fid is None:
            continue
        is_ex_row = bool(row.get("excluded_from_biological_interpretation"))
        factor_map[fid] = {
            "factor": fid,
            "rna_pct": r_pct,
            "atac_pct": a_pct,
            "classification": (
                "excluded_confounder"
                if is_ex_row
                else str(row.get("category") or row.get("classification") or "active")
            ),
            "excluded": is_ex_row,
            "reasons": [str(r) for r in (row.get("exclusion_reasons") or [])],
        }
    for row in excluded:
        if not isinstance(row, dict):
            continue
        fid, r_pct, a_pct = _extract_fac_r2(row)
        if fid is None:
            continue
        factor_map[fid] = {
            "factor": fid,
            "rna_pct": r_pct,
            "atac_pct": a_pct,
            "classification": "excluded_confounder",
            "excluded": True,
            "reasons": [str(r) for r in (row.get("exclusion_reasons") or [])],
        }
    if not factor_map:
        return None

    ordered = [factor_map[k] for k in sorted(factor_map.keys())]
    n_fac = len(ordered)
    mat = np.array([[r["rna_pct"], r["atac_pct"]] for r in ordered], dtype=np.float64)

    fig_h = max(4.6, 0.42 * n_fac + 1.9)
    fig, (ax_hm, ax_bar) = plt.subplots(
        1, 2,
        figsize=(11.8, fig_h),
        gridspec_kw={"width_ratios": [1.0, 1.85]},
        constrained_layout=True,
    )

    im = ax_hm.imshow(mat, aspect="auto", cmap="YlGnBu", vmin=0.0)
    ax_hm.set_xticks([0, 1])
    ax_hm.set_xticklabels(["RNA R² (%)", "ATAC R² (%)"], fontsize=9, fontweight="bold")
    y_labels = []
    for r in ordered:
        tag = " [EXCLUDED]" if r["excluded"] else f" [{r['classification']}]"
        y_labels.append(f"Factor {r['factor']}{tag}")
    ax_hm.set_yticks(np.arange(n_fac))
    ax_hm.set_yticklabels(y_labels, fontsize=8.4)
    vmax = float(np.nanmax(mat)) if mat.size else 1.0
    for i in range(n_fac):
        for j in range(2):
            val = mat[i, j]
            txt_color = "white" if val > 0.55 * max(vmax, 1e-6) else "#0f172a"
            ax_hm.text(
                j, i, f"{val:.2f}%",
                ha="center", va="center",
                fontsize=7.8, color=txt_color,
            )
    ax_hm.set_title("Per-Factor Variance Explained", fontsize=10, fontweight="bold")
    fig.colorbar(im, ax=ax_hm, fraction=0.046, pad=0.04, label="Variance explained (%)")

    y_pos = np.arange(n_fac)
    height = 0.36
    for idx, r in enumerate(ordered):
        is_ex = r["excluded"]
        rna_color = "#94a3b8" if is_ex else "#2563eb"
        atac_color = "#cbd5e1" if is_ex else "#0d9488"
        hatch = "///" if is_ex else None
        edge = "#b91c1c" if is_ex else "none"
        ax_bar.barh(
            idx - height / 2,
            r["rna_pct"],
            height=height,
            color=rna_color,
            edgecolor=edge,
            hatch=hatch,
            linewidth=0.9 if is_ex else 0.0,
        )
        ax_bar.barh(
            idx + height / 2,
            r["atac_pct"],
            height=height,
            color=atac_color,
            edgecolor=edge,
            hatch=hatch,
            linewidth=0.9 if is_ex else 0.0,
        )
        if is_ex:
            reason_short = "depth/cell-cycle confounder"
            joined = " ".join(r["reasons"]).lower()
            if "depth" in joined and "cycle" in joined:
                reason_short = "excluded: depth + cell-cycle"
            elif "depth" in joined:
                reason_short = "excluded: depth |r|>=0.70"
            elif "cycle" in joined:
                reason_short = "excluded: cell-cycle"
            x_end = max(r["rna_pct"], r["atac_pct"])
            ax_bar.text(
                x_end + max(vmax * 0.02, 0.08),
                idx,
                reason_short,
                va="center",
                ha="left",
                fontsize=7.8,
                color="#b91c1c",
                fontweight="bold",
            )

    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels([f"F{r['factor']}" for r in ordered], fontsize=8.4)
    ax_bar.invert_yaxis()
    ax_hm.set_ylim(n_fac - 0.5, -0.5)
    ax_bar.set_xlim(0.0, max(vmax * 1.38, 1.0))
    ax_bar.set_xlabel("Variance Explained R² (%)", fontsize=9.5)
    ax_bar.set_title(
        "Active Biological Programs vs Excluded Confounder Factors",
        fontsize=10,
        fontweight="bold",
    )
    legend_handles = [
        Patch(facecolor="#2563eb", label="Active factor: RNA R² (%)"),
        Patch(facecolor="#0d9488", label="Active factor: ATAC R² (%)"),
        Patch(
            facecolor="#cbd5e1",
            edgecolor="#b91c1c",
            hatch="///",
            label="Excluded technical / cell-cycle confounder",
        ),
    ]
    ax_bar.legend(handles=legend_handles, loc="lower right", fontsize=8.0, frameon=True)
    for spine in ("top", "right"):
        ax_bar.spines[spine].set_visible(False)

    ve_non_add = mofa_findings.get("variance_explained") or {}
    ve_add = mofa_findings.get("variance_explained_additive") or {}
    rna_tot = float(ve_non_add.get("rna") or 0.0) * 100.0
    atac_tot = float(ve_non_add.get("atac") or 0.0) * 100.0
    rna_add = float(ve_add.get("rna") or 0.0) * 100.0
    atac_add = float(ve_add.get("atac") or 0.0) * 100.0
    fig.suptitle(
        "MOFA+ Joint Variance Decomposition "
        f"(non-additive total R²: RNA={rna_tot:.2f}%, ATAC={atac_tot:.2f}% | "
        f"additive sum: RNA={rna_add:.2f}%, ATAC={atac_add:.2f}%)",
        fontsize=10.5,
        fontweight="bold",
        color="#0f172a",
    )
    return _save_dual(fig, Path(output_dir) / "mofa_variance_explained.png")


def _load_p2g_rows(p2g_findings: dict[str, Any]) -> list[dict[str, Any]]:
    csv_path = p2g_findings.get("output_path")
    rows: list[dict[str, Any]] = []
    if csv_path and Path(csv_path).exists():
        try:
            with open(csv_path, "r", encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh)
                for i, row in enumerate(reader):
                    if i >= 600:
                        break
                    rows.append(dict(row))
        except Exception as exc:
            log.warning("Could not read peak_gene_links.csv (%s); using top_links", exc)
    if not rows:
        rows = [dict(r) for r in (p2g_findings.get("top_links") or []) if isinstance(r, dict)]
    return rows


def render_regulatory_triangulation(
    p2g_findings: dict[str, Any] | None,
    ereg_findings: dict[str, Any] | None,
    output_dir: Path,
) -> tuple[str, str] | None:
    """Render regulatory triangulation summary across peak2gene & candidate eRegulons."""
    p2g_ok = (
        isinstance(p2g_findings, dict)
        and str(p2g_findings.get("status") or "").lower() in {"success", "ok", "done"}
        and int(p2g_findings.get("n_links") or 0) > 0
    )
    ereg_ok = (
        isinstance(ereg_findings, dict)
        and str(ereg_findings.get("status") or "").lower() in {"success", "ok", "done"}
        and int(ereg_findings.get("n_eregulons") or 0) > 0
    )
    if not p2g_ok and not ereg_ok:
        return None

    if p2g_ok and ereg_ok:
        fig, (ax_p2g, ax_ereg) = plt.subplots(
            1, 2, figsize=(13.6, 5.4), constrained_layout=True
        )
    elif p2g_ok:
        fig, ax_p2g = plt.subplots(1, 1, figsize=(7.6, 5.0), constrained_layout=True)
        ax_ereg = None
    else:
        fig, ax_ereg = plt.subplots(1, 1, figsize=(7.8, 5.0), constrained_layout=True)
        ax_p2g = None

    if ax_p2g is not None and isinstance(p2g_findings, dict):
        p2g_rows = _load_p2g_rows(p2g_findings)
        xs_prom, ys_prom, ss_prom = [], [], []
        xs_dist, ys_dist, ss_dist = [], [], []
        ranked_annot: list[tuple[float, float, float, str]] = []
        for row in p2g_rows:
            try:
                raw_corr = (
                    row.get("correlation")
                    if row.get("correlation") is not None
                    else row.get("corr")
                )
                corr = float(raw_corr)
                dist_bp = float(row.get("distance_bp") or row.get("tss_distance_bp") or 0.0)
                padj = max(float(row.get("padj") or 0.05), 1e-300)
            except (TypeError, ValueError):
                continue
            if not (math.isfinite(corr) and math.isfinite(dist_bp)):
                continue
            dist_kb = dist_bp / 1000.0
            neg_log_p = min(-math.log10(padj), 25.0)
            sz = 18.0 + 4.5 * neg_log_p
            reg_cls = str(row.get("region_class") or "distal_cis")
            gene = str(row.get("gene") or "")
            if "promoter" in reg_cls:
                xs_prom.append(dist_kb)
                ys_prom.append(corr)
                ss_prom.append(sz)
            else:
                xs_dist.append(dist_kb)
                ys_dist.append(corr)
                ss_dist.append(sz)
            if gene:
                ranked_annot.append((abs(corr) * neg_log_p, dist_kb, corr, gene))

        if xs_dist:
            ax_p2g.scatter(
                xs_dist, ys_dist, s=ss_dist, c="#0d9488", alpha=0.72,
                edgecolors="white", linewidths=0.4,
                label=f"Distal cis links (n={len(xs_dist)})",
                rasterized=len(xs_dist) > 250,
            )
        if xs_prom:
            ax_p2g.scatter(
                xs_prom, ys_prom, s=ss_prom, c="#2563eb", alpha=0.80,
                edgecolors="white", linewidths=0.4,
                label=f"Promoter-proximal <=2kb (n={len(xs_prom)})",
                rasterized=len(xs_prom) > 250,
            )
        ax_p2g.axhline(0.0, color="#64748b", ls="--", lw=0.9)
        ranked_annot.sort(key=lambda t: t[0], reverse=True)
        seen_genes: set[str] = set()
        for _, x_kb, y_r, g_sym in ranked_annot:
            if g_sym in seen_genes:
                continue
            seen_genes.add(g_sym)
            ax_p2g.annotate(
                g_sym,
                (x_kb, y_r),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7.8,
                color="#0f172a",
                fontweight="bold",
            )
            if len(seen_genes) >= 6:
                break
        ax_p2g.set_xlabel("Distance to Canonical TSS (kb)", fontsize=9.2)
        ax_p2g.set_ylabel("Peak-Gene Pearson r (shared cells)", fontsize=9.2)
        ax_p2g.set_title(
            f"Candidate Cis Peak-to-Gene Links "
            f"(total BH-FDR={p2g_findings.get('n_links', len(p2g_rows))})",
            fontsize=10,
            fontweight="bold",
        )
        ax_p2g.legend(loc="best", fontsize=8.0, frameon=True)
        for spine in ("top", "right"):
            ax_p2g.spines[spine].set_visible(False)

    if ax_ereg is not None and isinstance(ereg_findings, dict):
        mods = [
            m for m in (ereg_findings.get("top_eregulons") or [])
            if isinstance(m, dict) and m.get("tf")
        ][:12]
        if mods:
            labels_e = []
            scores_e = []
            colors_e = []
            notes_e = []
            for m in mods:
                tf = str(m.get("tf"))
                sc_cls = str(m.get("sign_class") or "pos_tf_pos_p2g")
                raw_padj = (
                    m.get("regulon_padj")
                    if m.get("regulon_padj") is not None
                    else m.get("padj")
                )
                padj = max(float(raw_padj if raw_padj is not None else 0.05), 1e-50)
                neg_lp = min(-math.log10(padj), 30.0)
                ng = int(m.get("n_target_genes") or 0)
                npk = int(m.get("n_enhancer_peaks") or 0)
                au_r = (
                    m.get("rna_atac_activity_corr")
                    if m.get("rna_atac_activity_corr") is not None
                    else m.get("rna_atac_aucell_corr")
                )
                au_txt = f"{float(au_r):.2f}" if isinstance(au_r, (int, float)) else "NA"
                flags = []
                if m.get("motif_family_ambiguity"):
                    flags.append("paralog-motif")
                if m.get("depth_confounded"):
                    flags.append("depth-conf")
                flag_str = f" [{'/'.join(flags)}]" if flags else ""
                sign_tag = "+" if sc_cls == "pos_tf_pos_p2g" else "-"
                labels_e.append(f"{tf} ({sign_tag}){flag_str}")
                scores_e.append(neg_lp)
                colors_e.append(
                    "#0f766e" if sc_cls == "pos_tf_pos_p2g" else "#d97706"
                )
                notes_e.append(f"{ng}g / {npk}p | AUCell r={au_txt}")

            y_pos = np.arange(len(mods))
            ax_ereg.barh(y_pos, scores_e, color=colors_e, alpha=0.86, height=0.64)
            xmax = max(scores_e) if scores_e else 2.0
            for idx, (sc_val, note) in enumerate(zip(scores_e, notes_e)):
                ax_ereg.text(
                    sc_val + xmax * 0.02,
                    idx,
                    note,
                    va="center",
                    ha="left",
                    fontsize=7.6,
                    color="#334155",
                )
            ax_ereg.set_yticks(y_pos)
            ax_ereg.set_yticklabels(labels_e, fontsize=8.4)
            ax_ereg.invert_yaxis()
            ax_ereg.set_xlim(0.0, xmax * 1.48)
            ax_ereg.set_xlabel("-log10(BH-adjusted p-value)", fontsize=9.2)
            ax_ereg.set_title(
                f"Candidate eRegulon Modules "
                f"(n={ereg_findings.get('n_eregulons', len(mods))}; associative)",
                fontsize=10,
                fontweight="bold",
            )
            legend_handles = [
                Patch(facecolor="#0f766e", label="pos_tf_pos_p2g (concordant candidate)"),
                Patch(facecolor="#d97706", label="neg_tf_pos_p2g (inverse candidate)"),
            ]
            ax_ereg.legend(handles=legend_handles, loc="lower right", fontsize=7.8, frameon=True)
            for spine in ("top", "right"):
                ax_ereg.spines[spine].set_visible(False)

    fig.suptitle(
        "Cross-Modal Regulatory Triangulation (Observational Cis & Motif-Linked Candidates)",
        fontsize=10.5,
        fontweight="bold",
        color="#0f172a",
    )
    return _save_dual(fig, Path(output_dir) / "regulatory_triangulation_summary.png")


def _convert_mofa_loadings_json_to_csv(src_json: Path, dest_csv: Path) -> bool:
    """Convert ``mofa_loadings.json`` into a flat publication CSV table."""
    payload = json.loads(src_json.read_text(encoding="utf-8"))
    factors = payload.get("factors") if isinstance(payload, dict) else None
    if not isinstance(factors, dict):
        return False
    dest_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(dest_csv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["factor", "modality", "feature", "weight", "abs_weight", "direction"])
        for fac_name, views_map in factors.items():
            if not isinstance(views_map, dict):
                continue
            for mod_name, rows in views_map.items():
                for item in rows or []:
                    if not isinstance(item, dict):
                        continue
                    writer.writerow([
                        fac_name,
                        mod_name,
                        item.get("feature", ""),
                        item.get("weight", ""),
                        item.get("abs_weight", ""),
                        item.get("direction", ""),
                    ])
    return dest_csv.is_file() and dest_csv.stat().st_size > 0


def stage_integration_tables(
    findings: dict[str, Any],
    tables_dir: Path,
) -> dict[str, str]:
    """Copy integration CSV tables into ``report_dir/tables/`` and return relative paths."""
    if not isinstance(findings, dict):
        return {}
    tables_dir = Path(tables_dir)
    tables_dir.mkdir(parents=True, exist_ok=True)
    staged: dict[str, str] = findings.setdefault("staged_tables", {})

    candidates = [
        ("mofa_factor_loadings", (findings.get("mofa") or {}).get("factor_loadings_path"), "mofa_factor_loadings.csv"),
        ("peak_gene_links", (findings.get("peak_to_gene") or {}).get("output_path"), "peak_gene_links.csv"),
        ("eregulons", (findings.get("eregulons") or {}).get("output_path"), "eregulons.csv"),
        ("eregulon_triplets", (findings.get("eregulons") or {}).get("triplets_path"), "eregulon_triplets.csv"),
    ]
    for key, src_raw, dest_name in candidates:
        if not src_raw:
            continue
        src = Path(str(src_raw))
        if not src.exists() or not src.is_file():
            continue
        dest = tables_dir / dest_name
        try:
            if key == "mofa_factor_loadings" and src.suffix.lower() == ".json":
                if _convert_mofa_loadings_json_to_csv(src, dest):
                    staged[key] = f"tables/{dest_name}"
                continue
            if src.resolve() != dest.resolve():
                shutil.copy2(src, dest)
            staged[key] = f"tables/{dest_name}"
        except Exception as exc:
            log.warning("Failed staging integration table %s: %s", src, exc)
    return staged


def generate_figures(
    findings: dict[str, Any],
    output_dir: str | Path,
) -> dict[str, Any]:
    """Render all data-backed integration figures (dual PNG + SVG) into ``output_dir``."""
    if not isinstance(findings, dict):
        return {}
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    figs = findings.setdefault("figures", {})
    figs_svg = findings.setdefault("figures_svg", {})

    wnn = findings.get("wnn")
    mofa = findings.get("mofa")
    p2g = findings.get("peak_to_gene")
    ereg = findings.get("eregulons")

    render_jobs = [
        ("wnn_joint_vs_modality_umap", lambda: render_wnn_umap_comparison(wnn, out_dir)),
        ("wnn_modality_weights_by_cluster", lambda: render_wnn_modality_weights(wnn, out_dir)),
        ("mofa_variance_explained", lambda: render_mofa_variance_explained(mofa, out_dir)),
        ("regulatory_triangulation_summary", lambda: render_regulatory_triangulation(p2g, ereg, out_dir)),
    ]
    for key, fn in render_jobs:
        try:
            dual = fn()
            if dual is not None:
                png_p, svg_p = dual
                figs[key] = png_p
                figs_svg[key] = svg_p
        except Exception as exc:
            log.warning("Integration figure %s failed: %s", key, exc, exc_info=True)

    return findings
