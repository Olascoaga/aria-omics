"""
ARIA scATAC cluster QC figures (scATAC P0 W0.2).

Reads a clustered scATAC peak `.h5ad` and renders, in the chromatin stack:

  - a MARKER-PEAK accessibility heatmap (clusters x top per-cluster marker peaks,
    mean accessibility, z-scored across clusters), and
  - a per-cluster QC depth VIOLIN (per-cell total peak counts by cluster).

Both need the AnnData matrix, so this runs as a subprocess in `aria-chromatin-env`
(via EnvironmentManager). No fabrication: a panel is produced only when its inputs
exist (marker peaks for the heatmap; >=1 cluster for the violin).

Subprocess interface (aria.scripts._base.run_script):
    python chromatin_figure_clusters.py <in.json> <out.json>

Input params:
    h5ad_path:     str  — clustered peak .h5ad (with `cluster_key` in obs)
    output_dir:    str  — directory to write PNGs into
    cluster_key:   str (optional, default "leiden")
    marker_peaks:  {cluster_id: [peak, ...]} (optional) — for the heatmap
    top_per_cluster: int (optional, default 5) — marker peaks/cluster in the heatmap

Output:
    {status, figures: {"marker_heatmap": path, "qc_depth_violin": path}}
"""

from __future__ import annotations

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from pathlib import Path
from aria.scripts._base import run_script


def _save_dual(fig, png_path) -> str:
    """P4.1: save a line-art figure as raster PNG (inlined in the report) AND a
    publication-grade vector SVG alongside it (`<stem>.svg`), matching the RNA side.
    Only the PNG path is returned/surfaced; the SVG sits next to it for manuscript use."""
    png_path = Path(png_path)
    fig.savefig(png_path, dpi=160, bbox_inches="tight", facecolor="white")
    fig.savefig(png_path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    return str(png_path)


def _marker_heatmap(adata, cluster_key, marker_peaks, top_per_cluster, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy import sparse

    # Ordered, de-duplicated peak list (top N per cluster).
    peaks: list = []
    for cl in sorted(marker_peaks):
        for p in list(marker_peaks[cl])[:top_per_cluster]:
            if p not in peaks and p in adata.var_names:
                peaks.append(p)
    if not peaks:
        return None
    clusters = sorted(adata.obs[cluster_key].astype(str).unique())
    col = {p: i for i, p in enumerate(adata.var_names)}
    idx = [col[p] for p in peaks]
    X = adata.X
    X = X.tocsc() if sparse.issparse(X) else X

    # Mean accessibility per cluster x peak.
    mat = np.zeros((len(clusters), len(peaks)), dtype=float)
    obs_cl = adata.obs[cluster_key].astype(str).values
    for ci, cl in enumerate(clusters):
        rows = np.where(obs_cl == cl)[0]
        if len(rows) == 0:
            continue
        sub = X[rows][:, idx] if sparse.issparse(X) else X[np.ix_(rows, idx)]
        mat[ci, :] = (np.asarray(sub.mean(axis=0)).ravel()
                      if sparse.issparse(sub) else sub.mean(axis=0))
    # Z-score per peak (column) so the cluster that "owns" each marker stands out.
    mu = mat.mean(axis=0, keepdims=True)
    sd = mat.std(axis=0, keepdims=True)
    z = (mat - mu) / np.where(sd > 0, sd, 1.0)

    fig, ax = plt.subplots(
        figsize=(max(5.0, 0.18 * len(peaks) + 2),
                 max(3.0, 0.4 * len(clusters) + 1)), dpi=160)
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(clusters, fontsize=7)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(peaks)} top marker peaks (cluster-ordered)", fontsize=8)
    ax.set_ylabel(cluster_key, fontsize=9)
    ax.set_title("Marker-peak accessibility (z-scored per peak)", fontsize=10,
                 fontweight="bold")
    plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="z(mean accessibility)")
    out = _save_dual(fig, out_path)
    plt.close(fig)
    return out


def _qc_depth_violin(adata, cluster_key, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import sparse

    X = adata.X
    depth = (np.asarray(X.sum(axis=1)).ravel() if sparse.issparse(X)
             else np.asarray(X).sum(axis=1))
    log_depth = np.log10(np.clip(depth, 1, None))
    clusters = sorted(adata.obs[cluster_key].astype(str).unique())
    obs_cl = adata.obs[cluster_key].astype(str).values
    # Pair each cluster with its data so tick labels can never drift from the
    # plotted positions when an empty cluster is dropped.
    pairs = [(cl, log_depth[obs_cl == cl]) for cl in clusters]
    pairs = [(cl, d[np.isfinite(d)]) for cl, d in pairs]
    pairs = [(cl, d) for cl, d in pairs if len(d) > 0]
    if not pairs:
        return None
    labels = [cl for cl, _ in pairs]
    positions = list(range(1, len(pairs) + 1))

    def _kde_ok(d) -> bool:
        # Gaussian KDE needs >= 2 points and non-zero variance (a singleton or a
        # constant-depth cluster raises LinAlgError / singular covariance).
        return len(d) >= 2 and float(np.ptp(d)) > 0.0

    fig, ax = plt.subplots(figsize=(max(4.0, 0.5 * len(labels) + 2), 4),
                           dpi=160)
    kde_pos = [p for p, (_, d) in zip(positions, pairs) if _kde_ok(d)]
    kde_dat = [d for _, d in pairs if _kde_ok(d)]
    if kde_dat:
        parts = ax.violinplot(kde_dat, positions=kde_pos, showmedians=True,
                              showextrema=False)
        for body in parts["bodies"]:
            body.set_alpha(0.7)
    degenerate = [(p, d) for p, (_, d) in zip(positions, pairs) if not _kde_ok(d)]
    if degenerate:
        # Honest fallback: raw points + a median tick (no fabricated density).
        for p, d in degenerate:
            ax.scatter(np.full(len(d), p), d, s=10, color="#1f77b4",
                       alpha=0.8, zorder=3, linewidths=0)
            ax.hlines(float(np.median(d)), p - 0.2, p + 0.2, color="#1f77b4",
                      lw=1.2, zorder=3)
        ax.text(0.99, 0.01, "points + median tick: n<2 or zero-variance cluster",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=6,
                color="#555555")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_xlim(0.4, len(positions) + 0.6)
    ax.set_xlabel(cluster_key, fontsize=9)
    ax.set_ylabel("log10 per-cell peak counts", fontsize=9)
    ax.set_title("Per-cluster accessibility depth", fontsize=10,
                 fontweight="bold")
    ax.tick_params(labelsize=7)
    out = _save_dual(fig, out_path)
    plt.close(fig)
    return out


def make_cluster_figures(params: dict) -> dict:
    from aria.utils.safe_h5ad import read_h5ad

    h5ad_path = params.get("h5ad_path")
    out_dir = Path(params["output_dir"])
    cluster_key = params.get("cluster_key", "leiden")
    marker_peaks = params.get("marker_peaks") or {}
    top_per_cluster = int(params.get("top_per_cluster", 5))

    if not h5ad_path or not Path(h5ad_path).exists():
        return {"status": "error", "error_type": "MissingInput",
                "details": f"h5ad_path not found: {h5ad_path}"}
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = read_h5ad(h5ad_path)
    if cluster_key not in adata.obs.columns:
        return {"status": "error", "error_type": "NoClusterKey",
                "details": f"{cluster_key!r} not in obs"}

    figures: dict = {}
    if marker_peaks:
        try:
            p = _marker_heatmap(adata, cluster_key, marker_peaks,
                                top_per_cluster, out_dir / "marker_heatmap.png")
            if p:
                figures["marker_heatmap"] = p
        except Exception as e:
            figures.setdefault("_errors", []).append(f"heatmap: {e}")
    try:
        p = _qc_depth_violin(adata, cluster_key, out_dir / "qc_depth_violin.png")
        if p:
            figures["qc_depth_violin"] = p
    except Exception as e:
        figures.setdefault("_errors", []).append(f"violin: {e}")

    return {"status": "success", "figures": figures}


if __name__ == "__main__":
    run_script(make_cluster_figures)
