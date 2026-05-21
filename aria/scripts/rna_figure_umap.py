"""
ARIA scRNA 2D embedding figure generator.

Reads an h5ad with a precomputed 2D embedding (UMAP preferred; t-SNE and
other `X_*` 2D embeddings accepted) and writes one PNG per requested `obs`
column to colour-by. If no 2D embedding exists, it attempts to compute UMAP.

Subprocess interface (aria.scripts._base.run_script):
    python rna_figure_umap.py <in.json> <out.json>

Input params:
    h5ad_path:       str  — path to AnnData
    color_by:        list[str]  — obs column names to color cells by
    output_dir:      str  — directory to write PNGs into
    embedding_key:   str (optional)  — explicit obsm key (default: auto)
    max_cells:       int (optional)  — subsample for plotting (default: 60_000)
    point_size:      float (optional) — matplotlib s= (default: auto from N)

Output:
    {status, embedding_key, embedding_label, n_cells_plotted,
     figures: {color_key: path, ...}}
"""

from __future__ import annotations

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from pathlib import Path
from aria.scripts._base import run_script


def _is_2d_embedding(adata, key: str) -> bool:
    try:
        arr = adata.obsm[key]
        return getattr(arr, "ndim", 0) == 2 and arr.shape[1] >= 2
    except Exception:
        return False


def _pick_embedding(adata, requested: str | None) -> str:
    """Find a 2D obsm embedding. Prefer explicit > UMAP > t-SNE > any X_*."""
    obsm_keys = list(adata.obsm.keys())
    if requested and requested in obsm_keys and _is_2d_embedding(adata, requested):
        return requested
    if "X_umap" in obsm_keys and _is_2d_embedding(adata, "X_umap"):
        return "X_umap"
    for k in obsm_keys:
        if k.lower().startswith("x_umap") and _is_2d_embedding(adata, k):
            return k
    if "X_tsne" in obsm_keys and _is_2d_embedding(adata, "X_tsne"):
        return "X_tsne"
    for k in obsm_keys:
        if k.lower().startswith("x_tsne") and _is_2d_embedding(adata, k):
            return k
    for k in obsm_keys:
        if k.lower().startswith("x_") and _is_2d_embedding(adata, k):
            return k
    raise KeyError(
        f"No 2D embedding found in obsm. Available keys: {obsm_keys}"
    )


def _embedding_label(key: str) -> str:
    key_l = key.lower()
    if "umap" in key_l:
        return "UMAP"
    if "tsne" in key_l or "t-sne" in key_l:
        return "t-SNE"
    return key.replace("X_", "").replace("_", " ").strip() or key


def _compute_umap(adata) -> str:
    import scanpy as sc
    if "X_pca" not in adata.obsm:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.pca(adata, n_comps=min(50, max(2, adata.n_vars - 1)))
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
    sc.tl.umap(adata, random_state=0)
    return "X_umap"


def _subsample_indices(n_total: int, n_max: int, seed: int = 0):
    """Stable random subsample for plotting only (does not modify adata)."""
    import numpy as np
    if n_total <= n_max:
        return None
    rng = np.random.default_rng(seed)
    return rng.choice(n_total, size=n_max, replace=False)


def _categorical_colors(n: int) -> list:
    """Generate up to N distinct colors from a categorical-friendly palette."""
    import matplotlib.pyplot as plt
    if n <= 10:
        return list(plt.get_cmap("tab10").colors)[:n]
    if n <= 20:
        return list(plt.get_cmap("tab20").colors)[:n]
    # Cycle through several palettes for >20 categories
    base = (list(plt.get_cmap("tab20").colors)
            + list(plt.get_cmap("tab20b").colors)
            + list(plt.get_cmap("tab20c").colors))
    while len(base) < n:
        base = base + base
    return base[:n]


def _plot_one(coords, values, key: str, out_path: Path,
              point_size: float, title: str, embedding_label: str) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    fig, ax = plt.subplots(figsize=(6.5, 6), dpi=160)

    # Decide categorical vs continuous
    values_dtype = getattr(values, "dtype", None)
    is_categorical = (
        isinstance(values_dtype, pd.CategoricalDtype)
        or pd.api.types.is_object_dtype(values)
        or pd.api.types.is_string_dtype(values)
        or pd.api.types.is_bool_dtype(values)
    )

    if is_categorical:
        cats = pd.Series(values).astype(str).fillna("NA")
        unique = sorted(cats.unique())
        colors = _categorical_colors(len(unique))
        for cat, color in zip(unique, colors):
            mask = (cats == cat).values
            ax.scatter(coords[mask, 0], coords[mask, 1],
                       s=point_size, c=[color], label=cat,
                       linewidths=0, alpha=0.85)
        # Compact legend outside plot if many categories
        ncol = 1 if len(unique) <= 14 else 2
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
                  frameon=False, fontsize=7, ncol=ncol,
                  handletextpad=0.3, columnspacing=0.6,
                  markerscale=2.0)
    else:
        vals = np.asarray(values, dtype=float)
        sc = ax.scatter(coords[:, 0], coords[:, 1], s=point_size,
                        c=vals, cmap="viridis",
                        linewidths=0, alpha=0.9)
        plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02,
                     label=key)

    ax.set_xlabel(f"{embedding_label}-1", fontsize=9)
    ax.set_ylabel(f"{embedding_label}-2", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.tick_params(labelsize=7)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    return str(out_path)


def make_umap_figures(params: dict) -> dict:
    import numpy as np
    from aria.utils.safe_h5ad import read_h5ad

    h5ad_path = Path(params["h5ad_path"])
    color_by  = list(params.get("color_by") or [])
    out_dir   = Path(params["output_dir"])
    requested = params.get("embedding_key")
    max_cells = int(params.get("max_cells", 60_000))
    point_size = params.get("point_size")

    if not color_by:
        return {"status": "error",
                "error_type": "MissingParam",
                "details": "color_by must be a non-empty list"}

    adata = read_h5ad(h5ad_path)
    try:
        emb_key = _pick_embedding(adata, requested)
        embedding_was_computed = False
    except KeyError:
        emb_key = _compute_umap(adata)
        embedding_was_computed = True
    label = _embedding_label(emb_key)
    coords  = np.asarray(adata.obsm[emb_key])
    if coords.ndim != 2 or coords.shape[1] < 2:
        return {"status": "error",
                "error_type": "InvalidEmbedding",
                "details": f"{emb_key} is not a 2D embedding"}
    coords = coords[:, :2]
    n_total = coords.shape[0]

    idx = _subsample_indices(n_total, max_cells)
    if idx is not None:
        coords_plot = coords[idx]
    else:
        coords_plot = coords

    n_plot = coords_plot.shape[0]
    if point_size is None:
        # Auto-size: smaller dots for denser plots
        point_size = max(0.6, min(6.0, 6.0 - np.log10(max(n_plot, 100))))

    figures: dict = {}
    missing: list = []
    file_prefix = "umap" if label == "UMAP" else label.lower().replace("-", "")
    for key in color_by:
        if key not in adata.obs.columns:
            missing.append(key)
            continue
        values = adata.obs[key].values
        if idx is not None:
            values = values[idx]
        safe = (
            str(key).replace("/", "_").replace(" ", "_").replace(":", "_")
        )
        out_path = out_dir / f"{file_prefix}_{safe}.png"
        title = f"{label} — {key}"
        figures[key] = _plot_one(
            coords_plot, values, key, out_path, point_size, title, label
        )

    return {
        "status":          "success",
        "embedding_key":   emb_key,
        "embedding_label": label,
        "embedding_was_computed": embedding_was_computed,
        "n_cells_total":   int(n_total),
        "n_cells_plotted": int(n_plot),
        "figures":         figures,
        "missing_keys":    missing,
    }


if __name__ == "__main__":
    run_script(make_umap_figures)
