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
     with top candidate ``eRegulons`` (all four TF x peak2gene sign classes:
     ``pos_tf_pos_p2g``, ``neg_tf_pos_p2g``, ``pos_tf_neg_p2g``, ``neg_tf_neg_p2g``;
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
from matplotlib.lines import Line2D
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


# Wave 2A / ADR-064: strictly non-causal TF x peak2gene sign classes -> (colour, short tag, legend label).
_SIGN_CLASS_STYLE = {
    "pos_tf_pos_p2g": ("#0f766e", "TF+/p2g+", "pos_tf_pos_p2g (concordant TF+ / p2g+)"),
    "neg_tf_pos_p2g": ("#d97706", "TF-/p2g+", "neg_tf_pos_p2g (inverse TF- / p2g+)"),
    "pos_tf_neg_p2g": ("#7c3aed", "TF+/p2g-", "pos_tf_neg_p2g (discordant TF+ / p2g-)"),
    "neg_tf_neg_p2g": ("#be123c", "TF-/p2g-", "neg_tf_neg_p2g (double-negative TF- / p2g-)"),
}
_SIGN_CLASS_OTHER = ("#64748b", "?", "other / unclassified sign")

_RNA_VIEW_HINTS = ("rna", "gex", "expr")
_ATAC_VIEW_HINTS = ("atac", "peak", "chrom", "acc")

# Columns that describe donors / batches / QC flags, never cell state.
_NON_STATE_LABEL_TOKENS = (
    "donor", "sample", "batch", "library", "lane", "replicate", "barcode",
    "pass_", "_qc", "subset", "doublet",
)


def _finite_float(value: Any) -> float | None:
    """``float(value)`` when finite, else ``None`` (no falsy-``0.0`` traps)."""
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _classify_view_kind(name: str, views_meta: Any = None) -> str:
    """Resolve a MOFA+ view name to ``"rna"`` / ``"atac"`` / ``"other"``.

    Prefers the recorded ``mofa["views"][name]["kind"]``; otherwise mirrors
    ``integration_mofa._classify_modality_kind`` (case-insensitive hints), so
    ``"scRNA"``/``"GEX"``/``"RNA"`` and ``"scATAC"``/``"Peaks"``/``"ATAC"`` work."""
    if isinstance(views_meta, dict):
        meta = views_meta.get(name)
        kind = str((meta or {}).get("kind") or "").lower() if isinstance(meta, dict) else ""
        if kind in {"rna", "atac"}:
            return kind
    low = str(name).strip().lower()
    if any(h in low for h in _RNA_VIEW_HINTS):
        return "rna"
    if any(h in low for h in _ATAC_VIEW_HINTS):
        return "atac"
    return "other"


def _view_value(mapping: Any, kind: str, views_meta: Any = None) -> Any:
    """Value of the first key in ``mapping`` whose view resolves to ``kind``."""
    if not isinstance(mapping, dict):
        return None
    if kind in mapping:  # canonical lowercase key wins when present
        return mapping[kind]
    for key, val in mapping.items():
        if _classify_view_kind(key, views_meta) == kind:
            return val
    return None


def _view_key(mofa: dict[str, Any], kind: str) -> str | None:
    """Original view name (e.g. ``"scRNA"``) resolving to ``kind``, if any."""
    views_meta = mofa.get("views") if isinstance(mofa.get("views"), dict) else None
    candidates: list[str] = list(views_meta or {})
    for src in (mofa.get("variance_explained"), mofa.get("variance_explained_additive")):
        if isinstance(src, dict):
            candidates.extend(k for k in src if k not in candidates)
    for key in candidates:
        if _classify_view_kind(key, views_meta) == kind:
            return str(key)
    return None


def _label_str(value: Any) -> str:
    """Stable string label: ``3``, ``np.int64(3)``, ``3.0`` and ``"3"`` -> ``"3"``."""
    if isinstance(value, (bytes, bytearray, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(float(value)):
            return "NA"
        if float(value).is_integer():
            return str(int(value))
    if isinstance(value, (np.integer,)):
        return str(int(value))
    return str(value)


def _neglog10_display(pvals: list[float], headroom: float = 1.15,
                      all_zero_floor: float = 10.0) -> tuple[list[float], float | None]:
    """Data-adaptive ``-log10(p)`` via canonical ``stats.adaptive_neglog10``."""
    from aria.utils.stats import adaptive_neglog10

    vals, cap, has_zero = adaptive_neglog10(
        pvals, headroom=headroom, all_zero_floor=all_zero_floor
    )
    return [float(v) for v in vals.tolist()], (cap if has_zero else None)


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

            from aria.utils.mudata_io import _ATAC_KEYS, _RNA_KEYS, _match_modality

            mod_names = list(f["mod"].keys()) if "mod" in f else []
            rna_mod_key = _match_modality(mod_names, _RNA_KEYS) or "rna"
            atac_mod_key = _match_modality(mod_names, _ATAC_KEYS) or "atac"
            rna_umap = _aligned_mod_umap(rna_mod_key)
            atac_umap = _aligned_mod_umap(atac_mod_key)

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

    def _is_boolean(vals: Any) -> bool:
        if isinstance(vals, np.ndarray) and vals.dtype.kind == "b":
            return True
        seen = {str(v).strip().lower() for v in
                (vals.tolist() if isinstance(vals, np.ndarray) else vals)}
        seen -= {"na", "nan", "none", ""}
        return bool(seen) and seen <= {"true", "false", "0", "1"}

    for key in preferred:
        if key in obs:
            cand = _valid_categorical(obs[key])
            if cand is not None:
                return key, cand

    # Wave 2A: the biological joint clustering beats any generic obs column —
    # otherwise multi-donor cohorts were coloured by donor_id / sample_id.
    if "leiden_wnn" in obs:
        cand = _valid_categorical(obs["leiden_wnn"])
        if cand is not None:
            return "leiden_wnn", [_label_str(v) for v in (
                obs["leiden_wnn"].tolist() if isinstance(obs["leiden_wnn"], np.ndarray)
                else obs["leiden_wnn"])]

    for key in sorted(obs.keys()):
        low = key.lower()
        if low.startswith(skip_prefixes) or low.endswith("_weight"):
            continue
        if low.startswith(("rna_", "atac_")) or any(t in low for t in _NON_STATE_LABEL_TOKENS):
            continue
        vals = obs[key]
        if isinstance(vals, np.ndarray) and vals.dtype.kind in {"i", "u"}:
            continue  # integer obs columns are counts, not cell states
        if _is_boolean(vals):
            continue
        cand = _valid_categorical(vals)
        if cand is not None:
            return key, cand

    if "leiden_wnn" in obs:
        vals = obs["leiden_wnn"]
        labels = [_label_str(v) for v in (vals.tolist() if isinstance(vals, np.ndarray) else vals)]
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
            ncol=ncol,
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
                # Wave 2A: int64 / float-coded leiden_wnn never equalled the str
                # per-cluster keys, silently dropping every IQR whisker.
                cl_list = cl_vals.tolist() if isinstance(cl_vals, np.ndarray) else list(cl_vals)
                cl_arr = np.asarray([_label_str(x) for x in cl_list], dtype=object)
                rw_vals = np.asarray(rw_vals, dtype=np.float64)
                for cl_key in per_cl.keys():
                    mask = (cl_arr == _label_str(cl_key)) & np.isfinite(rw_vals)
                    if np.any(mask):
                        q25, q50, q75 = np.percentile(rw_vals[mask], [25, 50, 75])
                        iqr_by_cluster[_label_str(cl_key)] = (float(q25), float(q50), float(q75))

    rows = []
    for cl_id, meta in per_cl.items():
        if not isinstance(meta, dict):
            continue
        rna_w = _finite_float(meta.get("mean_rna_weight"))
        atac_w = _finite_float(meta.get("mean_atac_weight"))
        if rna_w is None or atac_w is None:
            continue
        rows.append({
            "cluster": _label_str(cl_id),
            "n_cells": int(meta.get("n_cells") or 0),
            "rna_w": rna_w,
            "atac_w": atac_w,
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
    handles, leg_labels = ax.get_legend_handles_labels()
    if iqr_by_cluster:
        handles.append(Line2D(
            [0], [0], color="#0f172a", lw=1.8, marker="o", markersize=5,
            markerfacecolor="white", markeredgecolor="#0f172a",
        ))
        leg_labels.append("Median & IQR (cell-level RNA weight)")
    ax.legend(handles, leg_labels, loc="lower right", fontsize=8.2, frameon=True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    return _save_dual(fig, Path(output_dir) / "wnn_modality_weights_by_cluster.png")


def _exclusion_label(reasons: list[str]) -> str:
    """Short, accurate figure label from the real ``exclusion_reasons``.

    ``sequencing_depth_correlated:atac(|r|=0.812)`` -> ``depth |r|=0.81``;
    ``cell_cycle_correlated(|r|=0.61)`` -> ``cell-cycle |r|=0.61``;
    ``batch_or_donor_confounded:donor_id(eta2=0.40)`` -> ``batch/donor (donor_id)``.
    No hardcoded cutoffs and no batch factor mislabelled as depth/cell-cycle."""
    import re

    parts: list[str] = []
    for raw in reasons:
        low = str(raw).lower()
        r_match = re.search(r"\|r\|\s*=\s*([0-9.]+)", low)
        r_txt = f" |r|={float(r_match.group(1)):.2f}" if r_match else ""
        if "batch" in low or "donor" in low:
            col = re.search(r"confounded:([^(\s]+)", str(raw))
            tag = f"batch/donor ({col.group(1)})" if col else "batch/donor"
        elif "depth" in low:
            tag = f"depth{r_txt}"
        elif "cycle" in low:
            tag = f"cell-cycle{r_txt}"
        else:
            tag = str(raw)[:32]
        if tag not in parts:
            parts.append(tag)
    if not parts:
        return "excluded (reason not recorded)"
    return "excluded: " + " + ".join(parts)


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
    # Wave 2A/6D: resolve all recorded MOFA+ views dynamically (2 views or 3+ views,
    # e.g. RNA + ATAC + ADT/Protein) preserving original modality names.
    views_meta = mofa_findings.get("views") if isinstance(mofa_findings.get("views"), dict) else None
    ve_non_add = mofa_findings.get("variance_explained") or {}
    ve_add = mofa_findings.get("variance_explained_additive") or {}

    view_keys: list[str] = []

    def _add_view_key(k: Any) -> None:
        ks = str(k).strip()
        if ks and ks not in view_keys:
            view_keys.append(ks)

    if isinstance(views_meta, dict):
        for k in views_meta:
            _add_view_key(k)
    for src in (ve_non_add, ve_add):
        if isinstance(src, dict):
            for k in src:
                _add_view_key(k)
    for row in [*active, *excluded]:
        if not isinstance(row, dict):
            continue
        r2_map = row.get("r2_by_view")
        if isinstance(r2_map, dict):
            for k in r2_map:
                _add_view_key(k)
        for k in row:
            if (
                isinstance(k, str)
                and k.startswith("r2_")
                and k not in {"r2_by_view", "r2_raw_by_view"}
            ):
                _add_view_key(k[3:])
    existing_kinds = {_classify_view_kind(k, views_meta) for k in view_keys}
    for row in [*active, *excluded]:
        if not isinstance(row, dict):
            continue
        if "rna_r2" in row and "rna" not in existing_kinds:
            _add_view_key("rna")
            existing_kinds.add("rna")
        if "atac_r2" in row and "atac" not in existing_kinds:
            _add_view_key("atac")
            existing_kinds.add("atac")
    if not view_keys:
        view_keys = ["rna", "atac"]

    view_kinds = [_classify_view_kind(k, views_meta) for k in view_keys]
    kind_counts = {
        kind: sum(1 for vk in view_kinds if vk == kind)
        for kind in set(view_kinds)
    }

    def _view_r2_for_key(row: dict[str, Any], v_key: str, v_kind: str) -> float:
        r2_map = row.get("r2_by_view")
        direct = r2_map.get(v_key) if isinstance(r2_map, dict) else None
        if direct is None:
            direct = row.get(f"r2_{v_key}")
        if direct is None and v_kind in {"rna", "atac"} and kind_counts.get(v_kind, 0) == 1:
            direct = row.get(f"{v_kind}_r2")
            if direct is None:
                direct = _view_value(r2_map, v_kind, views_meta)
            if direct is None:
                flat = {
                    k[3:]: v
                    for k, v in row.items()
                    if isinstance(k, str)
                    and k.startswith("r2_")
                    and k not in {"r2_by_view", "r2_raw_by_view"}
                }
                direct = _view_value(flat, v_kind, views_meta)
        val = _finite_float(direct)
        return (val if val is not None else 0.0) * 100.0

    def _extract_fac_r2(row: dict[str, Any]) -> tuple[int | None, list[float]]:
        raw_id = row.get("factor_id") if row.get("factor_id") is not None else row.get("factor")
        if raw_id is None:
            return None, []
        try:
            fid = int(raw_id)
        except (TypeError, ValueError):
            return None, []
        return fid, [_view_r2_for_key(row, vk, vkind) for vk, vkind in zip(view_keys, view_kinds)]

    batch_ids_top: set[int] = set()
    for bkey in ("batch_factor_ids", "batch_confounded_factors"):
        for item in mofa_findings.get(bkey) or []:
            raw_b = item
            if isinstance(item, dict):
                raw_b = item.get("factor_id") if item.get("factor_id") is not None else item.get("factor")
            try:
                if raw_b is not None:
                    batch_ids_top.add(int(raw_b))
            except (TypeError, ValueError):
                pass

    def _flag_reasons(row: dict[str, Any], fid: int | None = None) -> list[str]:
        reasons = [str(r) for r in (row.get("exclusion_reasons") or [])]
        low_joined = " ".join(reasons).lower()
        if row.get("technical_depth_flag") and "depth" not in low_joined:
            reasons.append("sequencing_depth_correlated")
        if row.get("cell_cycle_flag") and "cycle" not in low_joined:
            reasons.append("cell_cycle_correlated")
        if (row.get("batch_flag") or (fid is not None and fid in batch_ids_top)) and (
            "batch" not in low_joined and "donor" not in low_joined
        ):
            ba = row.get("batch_association")
            bcol = ba.get("batch_key") if isinstance(ba, dict) else None
            reasons.append(
                f"batch_or_donor_confounded:{bcol}" if bcol else "batch_or_donor_confounded"
            )
        return reasons

    for row in active:
        if not isinstance(row, dict):
            continue
        fid, v_pcts = _extract_fac_r2(row)
        if fid is None:
            continue
        is_ex_row = bool(
            row.get("excluded_from_biological_interpretation")
            or row.get("batch_flag")
            or row.get("technical_depth_flag")
            or row.get("cell_cycle_flag")
            or fid in batch_ids_top
        )
        factor_map[fid] = {
            "factor": fid,
            "view_pcts": v_pcts,
            "classification": (
                "excluded_confounder"
                if is_ex_row
                else str(row.get("category") or row.get("classification") or "active")
            ),
            "excluded": is_ex_row,
            "reasons": _flag_reasons(row, fid),
        }
    for row in excluded:
        if not isinstance(row, dict):
            continue
        fid, v_pcts = _extract_fac_r2(row)
        if fid is None:
            continue
        factor_map[fid] = {
            "factor": fid,
            "view_pcts": v_pcts,
            "classification": "excluded_confounder",
            "excluded": True,
            "reasons": _flag_reasons(row, fid),
        }
    if not factor_map:
        return None

    ordered = [factor_map[k] for k in sorted(factor_map.keys())]
    n_fac = len(ordered)
    n_views = len(view_keys)
    mat = np.array([r["view_pcts"] for r in ordered], dtype=np.float64)
    mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0)
    finite_max = float(np.max(mat)) if mat.size else 0.0
    # Guard degenerate Normalize(0, 0) when every R² is zero.
    vmax = max(finite_max, 1e-6)

    view_labels: list[str] = []
    view_short_labels: list[str] = []
    view_colors: list[str] = []
    extra_palette = ["#7c3aed", "#d97706", "#db2777", "#059669", "#ea580c", "#0284c7"]
    muted_palette = ["#94a3b8", "#cbd5e1", "#a8a29e", "#d4d4d8"]
    extra_idx = 0
    for vk, vkind in zip(view_keys, view_kinds):
        if vkind == "rna" and kind_counts.get("rna", 0) == 1:
            lab = "RNA" if vk.lower() == "rna" else f"RNA [{vk}]"
            short_lab = "RNA"
            col = "#2563eb"
        elif vkind == "atac" and kind_counts.get("atac", 0) == 1:
            lab = "ATAC" if vk.lower() == "atac" else f"ATAC [{vk}]"
            short_lab = "ATAC"
            col = "#0d9488"
        else:
            lab = vk
            short_lab = vk
            col = extra_palette[extra_idx % len(extra_palette)]
            extra_idx += 1
        view_labels.append(lab)
        view_short_labels.append(short_lab)
        view_colors.append(col)

    fig_h = max(4.6, max(0.42, 0.18 * n_views) * n_fac + 1.9)
    fig_w = max(11.8, 9.6 + 1.1 * n_views)
    fig, (ax_hm, ax_bar) = plt.subplots(
        1, 2,
        figsize=(fig_w, fig_h),
        gridspec_kw={"width_ratios": [max(1.0, 0.5 * n_views), 1.85]},
        constrained_layout=True,
    )

    im = ax_hm.imshow(mat, aspect="auto", cmap="YlGnBu", vmin=0.0, vmax=vmax)
    ax_hm.set_xticks(np.arange(n_views))
    ax_hm.set_xticklabels(
        [f"{lab} R² (%)" for lab in view_labels],
        fontsize=9,
        fontweight="bold",
    )
    y_labels = []
    for r in ordered:
        tag = " [EXCLUDED]" if r["excluded"] else f" [{r['classification']}]"
        y_labels.append(f"Factor {r['factor']}{tag}")
    ax_hm.set_yticks(np.arange(n_fac))
    ax_hm.set_yticklabels(y_labels, fontsize=8.4)
    for i in range(n_fac):
        for j in range(n_views):
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
    group_h = 0.72
    height = group_h / max(n_views, 1)
    for idx, r in enumerate(ordered):
        is_ex = r["excluded"]
        hatch = "///" if is_ex else None
        edge = "#b91c1c" if is_ex else "none"
        for j in range(n_views):
            y_off = idx - group_h / 2.0 + (j + 0.5) * height
            bar_color = muted_palette[j % len(muted_palette)] if is_ex else view_colors[j]
            ax_bar.barh(
                y_off,
                r["view_pcts"][j],
                height=height,
                color=bar_color,
                edgecolor=edge,
                hatch=hatch,
                linewidth=0.9 if is_ex else 0.0,
            )
        if is_ex:
            reason_short = _exclusion_label(r["reasons"])
            x_end = max(r["view_pcts"]) if r["view_pcts"] else 0.0
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
        Patch(facecolor=col, label=f"Active factor: {lab} R² (%)")
        for col, lab in zip(view_colors, view_short_labels)
    ] + [
        Patch(
            facecolor="#cbd5e1",
            edgecolor="#b91c1c",
            hatch="///",
            label="Excluded confounder (depth / cell-cycle / batch-donor)",
        ),
    ]
    ax_bar.legend(handles=legend_handles, loc="lower right", fontsize=8.0, frameon=True)
    for spine in ("top", "right"):
        ax_bar.spines[spine].set_visible(False)

    def _tot_pct(mapping: Any, v_key: str, v_kind: str) -> float:
        raw = mapping.get(v_key) if isinstance(mapping, dict) else None
        if raw is None and v_kind in {"rna", "atac"} and kind_counts.get(v_kind, 0) == 1:
            raw = _view_value(mapping, v_kind, views_meta)
        val = _finite_float(raw)
        return (val if val is not None else 0.0) * 100.0

    non_add_str = ", ".join(
        f"{sl}={_tot_pct(ve_non_add, vk, vkind):.2f}%"
        for sl, vk, vkind in zip(view_short_labels, view_keys, view_kinds)
    )
    add_str = ", ".join(
        f"{sl}={_tot_pct(ve_add, vk, vkind):.2f}%"
        for sl, vk, vkind in zip(view_short_labels, view_keys, view_kinds)
    )
    fig.suptitle(
        "MOFA+ Joint Variance Decomposition "
        f"(non-additive total R²: {non_add_str} | "
        f"additive sum: {add_str})",
        fontsize=10.5,
        fontweight="bold",
        color="#0f172a",
    )
    return _save_dual(fig, Path(output_dir) / "mofa_variance_explained.png")


_P2G_ROW_CAP = 600
_EREG_ROW_CAP = 12


def _load_p2g_rows_with_total(p2g_findings: dict[str, Any]) -> tuple[list[dict[str, Any]], int]:
    """Rows to plot (capped at ``_P2G_ROW_CAP``) plus the TOTAL row count, so the
    figure can disclose ``top N of M shown`` instead of silently truncating."""
    csv_path = p2g_findings.get("output_path")
    rows: list[dict[str, Any]] = []
    total = 0
    if csv_path and Path(csv_path).exists():
        try:
            with open(csv_path, "r", encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    total += 1
                    if len(rows) < _P2G_ROW_CAP:
                        rows.append(dict(row))
        except Exception as exc:
            log.warning("Could not read peak_gene_links.csv (%s); using top_links", exc)
            rows, total = [], 0
    if not rows:
        rows = [dict(r) for r in (p2g_findings.get("top_links") or []) if isinstance(r, dict)]
        total = len(rows)
    return rows, total


def _load_p2g_rows(p2g_findings: dict[str, Any]) -> list[dict[str, Any]]:
    return _load_p2g_rows_with_total(p2g_findings)[0]


def _first_present(row: dict[str, Any], *keys: str) -> Any:
    for key in keys:
        val = row.get(key)
        if val is not None and not (isinstance(val, str) and not val.strip()):
            return val
    return None


def _parse_padj(raw: Any, default: float = 0.05) -> float | None:
    """Explicit padj parsing: ``0.0`` stays ``0.0`` (no ``padj or 0.05`` trap);
    missing -> ``default``; non-finite / outside [0, 1] -> ``None`` (unplottable)."""
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        return default
    try:
        val = float(raw)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(val) or val < 0.0 or val > 1.0:
        return None
    return val


def _collect_p2g_points(p2g_findings: dict[str, Any]) -> tuple[list[dict[str, Any]], int]:
    rows, total = _load_p2g_rows_with_total(p2g_findings)
    pts: list[dict[str, Any]] = []
    for row in rows:
        corr = _finite_float(_first_present(row, "correlation", "corr"))
        dist_bp = _finite_float(_first_present(row, "distance_bp", "tss_distance_bp"))
        if dist_bp is None and _first_present(row, "distance_bp", "tss_distance_bp") is None:
            dist_bp = 0.0
        padj = _parse_padj(row.get("padj"))
        if corr is None or dist_bp is None or padj is None:
            continue
        pts.append({
            "x_kb": dist_bp / 1000.0, "r": corr, "padj": padj,
            "promoter": "promoter" in str(row.get("region_class") or "distal_cis"),
            "gene": str(row.get("gene") or ""),
        })
    return pts, total


def _collect_eregulon_bars(ereg_findings: dict[str, Any]) -> tuple[list[dict[str, Any]], int]:
    mods_all = [
        m for m in (ereg_findings.get("top_eregulons") or [])
        if isinstance(m, dict) and m.get("tf")
    ]
    bars: list[dict[str, Any]] = []
    for m in mods_all:
        padj = _parse_padj(_first_present(m, "regulon_padj", "padj"))
        if padj is None:
            continue  # NaN / inf / out-of-range padj is not plottable
        bars.append({"m": m, "padj": padj})
        if len(bars) >= _EREG_ROW_CAP:
            break
    return bars, len(mods_all)


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
    # Wave 2A: collect plottable data BEFORE allocating panels, so an empty side
    # never leaves a blank [0, 1] axes box (and no data -> no figure at all).
    p2g_pts: list[dict[str, Any]] = []
    p2g_total = 0
    if p2g_ok and isinstance(p2g_findings, dict):
        p2g_pts, p2g_total = _collect_p2g_points(p2g_findings)
    ereg_bars: list[dict[str, Any]] = []
    ereg_total = 0
    if ereg_ok and isinstance(ereg_findings, dict):
        ereg_bars, ereg_total = _collect_eregulon_bars(ereg_findings)
    has_p2g, has_ereg = bool(p2g_pts), bool(ereg_bars)
    if not has_p2g and not has_ereg:
        return None

    if has_p2g and has_ereg:
        fig, (ax_p2g, ax_ereg) = plt.subplots(
            1, 2, figsize=(13.6, 5.4), constrained_layout=True
        )
    elif has_p2g:
        fig, ax_p2g = plt.subplots(1, 1, figsize=(7.6, 5.0), constrained_layout=True)
        ax_ereg = None
    else:
        fig, ax_ereg = plt.subplots(1, 1, figsize=(7.8, 5.0), constrained_layout=True)
        ax_p2g = None

    if ax_p2g is not None and isinstance(p2g_findings, dict):
        neg_lps, p_cap = _neglog10_display([p["padj"] for p in p2g_pts])
        lp_top = max(neg_lps) if neg_lps else 0.0
        s_min, s_max = 14.0, 140.0

        def _p2g_area(v: float) -> float:
            return s_min if lp_top <= 0 else s_min + (s_max - s_min) * min(max(v / lp_top, 0.0), 1.0)

        groups = {True: ([], [], []), False: ([], [], [])}
        ranked_annot: list[tuple[float, float, float, str]] = []
        for pt, lp in zip(p2g_pts, neg_lps):
            xs, ys, ss = groups[pt["promoter"]]
            xs.append(pt["x_kb"])
            ys.append(pt["r"])
            ss.append(_p2g_area(lp))
            if pt["gene"]:
                ranked_annot.append((abs(pt["r"]) * lp, pt["x_kb"], pt["r"], pt["gene"]))
        xs_dist, ys_dist, ss_dist = groups[False]
        xs_prom, ys_prom, ss_prom = groups[True]
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
        # Deterministic ranking (score desc, then gene) + light label staggering.
        ranked_annot.sort(key=lambda t: (-t[0], t[3]))
        xs_all = [p["x_kb"] for p in p2g_pts]
        ys_all = [p["r"] for p in p2g_pts]
        x_span = (max(xs_all) - min(xs_all)) or 1.0
        y_span = (max(ys_all) - min(ys_all)) or 1.0
        offsets = [(4, 4), (4, -11), (-6, 12), (-6, -18), (8, 20), (8, -26)]
        placed: list[tuple[float, float]] = []
        seen_genes: set[str] = set()
        for _, x_kb, y_r, g_sym in ranked_annot:
            if g_sym in seen_genes:
                continue
            seen_genes.add(g_sym)
            n_near = sum(
                1 for px, py in placed
                if abs(px - x_kb) < 0.06 * x_span and abs(py - y_r) < 0.08 * y_span
            )
            dx, dy = offsets[n_near % len(offsets)]
            placed.append((x_kb, y_r))
            ax_p2g.annotate(
                g_sym,
                (x_kb, y_r),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=7.8,
                color="#0f172a",
                fontweight="bold",
                ha="right" if dx < 0 else "left",
                arrowprops=(dict(arrowstyle="-", color="#94a3b8", lw=0.5)
                            if n_near else None),
            )
            if len(seen_genes) >= 6:
                break
        ax_p2g.set_xlabel("Distance to Canonical TSS (kb)", fontsize=9.2)
        ax_p2g.set_ylabel("Peak-Gene Pearson r (shared cells)", fontsize=9.2)
        n_links = int(p2g_findings.get("n_links") or 0)
        total = max(n_links, p2g_total)
        shown = len(p2g_pts)
        trunc = f"; top {shown} of {total} shown" if shown < total else ""
        bg_none = str(p2g_findings.get("background_method") or "").strip().lower() == "none"
        no_fdr = bg_none or ("max_padj" in p2g_findings and p2g_findings.get("max_padj") is None)
        total_tag = f"total={n_links or total}" if no_fdr else f"total BH-FDR={n_links or total}"
        ax_p2g.set_title(
            f"Candidate Cis Peak-to-Gene Links ({total_tag}{trunc})",
            fontsize=10,
            fontweight="bold",
        )
        # Single legend (link classes + bubble-size key for -log10(padj)) placed
        # by matplotlib away from the data, so no point is hidden behind it.
        cls_handles, cls_labels = ax_p2g.get_legend_handles_labels()
        if lp_top > 0:
            levels = sorted({round(v, 1) for v in (lp_top / 4.0, lp_top / 2.0, lp_top) if v > 0})
        else:
            levels = [0.0]
        size_handles = [
            Line2D([0], [0], marker="o", ls="", markerfacecolor="#cbd5e1",
                   markeredgecolor="#475569", markersize=math.sqrt(_p2g_area(lv)))
            for lv in levels
        ]
        size_labels = [f"-log10(padj) = {lv:g}" for lv in levels]
        if p_cap is not None:
            size_labels[-1] = f"-log10(padj) = {levels[-1]:g} (padj=0 capped)"
        ax_p2g.legend(cls_handles + size_handles, cls_labels + size_labels,
                      loc="best", fontsize=7.2, frameon=True, labelspacing=0.9)
        for spine in ("top", "right"):
            ax_p2g.spines[spine].set_visible(False)

    if ax_ereg is not None and isinstance(ereg_findings, dict):
        scores_e, e_cap = _neglog10_display([b["padj"] for b in ereg_bars])
        labels_e = []
        colors_e = []
        notes_e = []
        classes_seen: list[str] = []
        for bar, sc_val in zip(ereg_bars, scores_e):
            m = bar["m"]
            tf = str(m.get("tf"))
            sc_cls = str(m.get("sign_class") or "pos_tf_pos_p2g")
            color, tag, _ = _SIGN_CLASS_STYLE.get(sc_cls, _SIGN_CLASS_OTHER)
            key = sc_cls if sc_cls in _SIGN_CLASS_STYLE else "__other__"
            if key not in classes_seen:
                classes_seen.append(key)
            ng = int(_finite_float(m.get("n_target_genes")) or 0)
            npk = int(_finite_float(m.get("n_enhancer_peaks")) or 0)
            au_r = _finite_float(_first_present(m, "rna_atac_activity_corr", "rna_atac_aucell_corr"))
            au_txt = f"{au_r:.2f}" if au_r is not None else "NA"
            flags = []
            if m.get("motif_family_ambiguity"):
                flags.append("paralog-motif")
            if m.get("depth_confounded"):
                flags.append("depth-conf")
            flag_str = f" [{'/'.join(flags)}]" if flags else ""
            labels_e.append(f"{tf} ({tag}){flag_str}")
            colors_e.append(color)
            cap_txt = " | padj=0 (capped)" if bar["padj"] <= 0 else ""
            notes_e.append(f"{ng}g / {npk}p | AUCell r={au_txt}{cap_txt}")

        y_pos = np.arange(len(ereg_bars))
        ax_ereg.barh(y_pos, scores_e, color=colors_e, alpha=0.86, height=0.64)
        xmax = max(max(scores_e) if scores_e else 0.0, 1.0)
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
        ax_ereg.set_xlabel(
            "-log10(BH-adjusted p-value)" + (" [padj=0 capped]" if e_cap is not None else ""),
            fontsize=9.2,
        )
        n_ereg = int(ereg_findings.get("n_eregulons") or 0) or ereg_total
        shown = len(ereg_bars)
        trunc = f"; top {shown} of {n_ereg} shown" if shown < n_ereg else ""
        ax_ereg.set_title(
            f"Candidate eRegulon Modules "
            f"(n={n_ereg}{trunc}; associative)",
            fontsize=10,
            fontweight="bold",
        )
        legend_handles = []
        for key in classes_seen:
            color, _, label = _SIGN_CLASS_STYLE.get(key, _SIGN_CLASS_OTHER)
            legend_handles.append(Patch(facecolor=color, label=label))
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
