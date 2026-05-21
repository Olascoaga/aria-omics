"""
ARIA scRNA PAGA + DPT figure generator.

Reads a trajectory.h5ad (from rna_trajectory.py) and writes:

  1. PAGA force-directed graph PNG — nodes coloured by group, sized by
     population, edge widths ∝ connectivity. Self-rendered with matplotlib
     (no scanpy.pl.paga dependency) to keep the figure publication-grade.

  2. (Optional) DPT-coloured UMAP, delegated to the same plotting helpers
     used by rna_figure_umap.py — but accepts the trajectory h5ad whose
     UMAP may be `X_umap` (just computed via tl.umap(init_pos='paga')) or
     an inherited `X_umap.rpca`.

Subprocess interface (aria.scripts._base.run_script):
    python rna_figure_paga.py <in.json> <out.json>

Input params:
    h5ad_path:        str  — path to trajectory.h5ad
    output_dir:       str  — directory to write PNGs into
    groupby:          str (optional) — obs column for PAGA nodes
                                       (default: read from uns['paga']['groups'])
    min_connectivity: float (optional) — drop edges weaker than this
                                         (default: 0 = keep all)

Output:
    {status, figures: {paga_graph: path, paga_log10_graph: path?}}
"""

from __future__ import annotations

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from pathlib import Path
from aria.scripts._base import run_script


def _categorical_colors(n: int) -> list:
    import matplotlib.pyplot as plt
    if n <= 10:
        return list(plt.get_cmap("tab10").colors)[:n]
    if n <= 20:
        return list(plt.get_cmap("tab20").colors)[:n]
    base = (list(plt.get_cmap("tab20").colors)
            + list(plt.get_cmap("tab20b").colors))
    while len(base) < n:
        base = base + base
    return base[:n]


def _spring_layout(conn, seed: int = 0, iterations: int = 100):
    """
    Force-directed layout for the PAGA group graph. Uses networkx if
    available (more stable), otherwise falls back to a simple Fruchterman-
    Reingold implementation on the dense matrix.
    """
    import numpy as np
    try:
        import networkx as nx
        G = nx.from_numpy_array(conn)
        pos = nx.spring_layout(G, seed=seed, iterations=iterations,
                               k=1.0 / max(1, conn.shape[0]) ** 0.5)
        return np.array([pos[i] for i in range(conn.shape[0])])
    except Exception:
        # Minimal FR fallback
        n = conn.shape[0]
        rng = np.random.default_rng(seed)
        pos = rng.uniform(-1, 1, size=(n, 2))
        k = 1.0 / (n ** 0.5)
        for _ in range(iterations):
            diff = pos[:, None, :] - pos[None, :, :]
            dist = np.linalg.norm(diff, axis=-1) + 1e-9
            # Repulsion ∝ 1/dist²
            rep = (k * k / dist)[..., None] * diff / dist[..., None]
            rep[np.eye(n, dtype=bool)] = 0
            # Attraction ∝ dist · weight
            attr = -(conn[..., None] * diff)
            disp = rep.sum(axis=1) + attr.sum(axis=1)
            d_norm = np.linalg.norm(disp, axis=1, keepdims=True) + 1e-9
            pos += disp / d_norm * 0.05
        return pos


def _draw_paga(conn, cats, sizes, output_path: Path,
               title: str = "PAGA — cluster graph",
               use_log: bool = False) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    pos = _spring_layout(conn)
    fig, ax = plt.subplots(figsize=(7, 6.5), dpi=160)
    colors = _categorical_colors(len(cats))

    # Edges first (so nodes overlay)
    if use_log:
        # log10(1 + 1000 * conn) compresses tiny values into visible widths
        edge_weights = np.log10(1.0 + 1000.0 * conn)
        w_max = edge_weights.max() if edge_weights.max() > 0 else 1.0
    else:
        edge_weights = conn
        w_max = edge_weights.max() if edge_weights.max() > 0 else 1.0

    n = conn.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            if conn[i, j] <= 0:
                continue
            lw = 0.5 + 4.5 * (edge_weights[i, j] / w_max)
            alpha = 0.30 + 0.65 * (edge_weights[i, j] / w_max)
            ax.plot(
                [pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]],
                color="#475569", linewidth=lw, alpha=alpha, zorder=1,
            )

    # Nodes — sized by sqrt(population) for readability
    s_min, s_max = 350.0, 2800.0
    sz_arr = np.array(sizes, dtype=float)
    if sz_arr.max() > sz_arr.min():
        sz_norm = (np.sqrt(sz_arr) - np.sqrt(sz_arr.min())) / (
            np.sqrt(sz_arr.max()) - np.sqrt(sz_arr.min())
        )
    else:
        sz_norm = np.full_like(sz_arr, 0.5)
    sz_plot = s_min + sz_norm * (s_max - s_min)

    ax.scatter(pos[:, 0], pos[:, 1], s=sz_plot, c=colors,
               edgecolors="#1e293b", linewidths=1.0, zorder=2)
    for i, cat in enumerate(cats):
        ax.annotate(cat, (pos[i, 0], pos[i, 1]),
                    ha="center", va="center", fontsize=9,
                    fontweight="bold", color="#0f172a", zorder=3)

    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    # Edge-width legend
    legend_lines = []
    if w_max > 0:
        from matplotlib.lines import Line2D
        anchors = [0.25, 0.5, 1.0]
        for a in anchors:
            lw = 0.5 + 4.5 * a
            wv = a * (conn.max() if not use_log else conn.max())
            legend_lines.append(
                Line2D([0], [0], color="#475569", linewidth=lw,
                       label=f"~{wv:.4f}" if conn.max() < 0.1
                             else f"~{wv:.2f}")
            )
        ax.legend(handles=legend_lines, loc="lower left",
                  frameon=False, fontsize=7,
                  title="edge weight" + (" (log-scaled)" if use_log else ""),
                  title_fontsize=7)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    return str(output_path)


def make_paga_figures(params: dict) -> dict:
    import numpy as np
    from aria.utils.safe_h5ad import read_h5ad

    h5ad_path = Path(params["h5ad_path"])
    out_dir   = Path(params["output_dir"])
    groupby   = params.get("groupby")
    min_conn  = float(params.get("min_connectivity", 0))

    adata = read_h5ad(h5ad_path)
    if "paga" not in adata.uns:
        return {"status": "error", "error_type": "NoPAGA",
                "details": "uns['paga'] not present in h5ad"}

    if groupby is None:
        groupby = adata.uns["paga"].get("groups", "leiden")
    if groupby not in adata.obs.columns:
        return {"status": "error", "error_type": "GroupNotFound",
                "details": f"obs column '{groupby}' missing"}

    conn = adata.uns["paga"]["connectivities"]
    if hasattr(conn, "toarray"):
        conn = conn.toarray()
    conn = np.asarray(conn, dtype=float)

    if min_conn > 0:
        conn = np.where(conn >= min_conn, conn, 0.0)

    cats = (list(adata.obs[groupby].cat.categories)
            if hasattr(adata.obs[groupby], "cat")
            else sorted(adata.obs[groupby].unique()))
    sizes = [int((adata.obs[groupby] == c).sum()) for c in cats]
    max_off = float((conn - np.diag(np.diag(conn))).max())

    out_dir.mkdir(parents=True, exist_ok=True)
    figures: dict = {}

    figures["paga_graph"] = _draw_paga(
        conn, cats, sizes, out_dir / "paga_graph.png",
        title=f"PAGA — cluster graph  (max edge = {max_off:.4f})",
        use_log=False,
    )
    if max_off > 0 and max_off < 0.05:
        # Adult / mature populations have absolute connectivities orders
        # of magnitude below the developmental case. A log-scaled twin
        # gives the reader a visible structure without misleading them
        # into thinking the absolute weights are large.
        figures["paga_log10_graph"] = _draw_paga(
            conn, cats, sizes, out_dir / "paga_graph_log10.png",
            title=("PAGA — cluster graph (log-scaled edges; "
                   "absolute connectivities are weak — see legend)"),
            use_log=True,
        )

    return {
        "status":            "success",
        "groupby":           groupby,
        "max_connectivity":  round(max_off, 5),
        "n_groups":          len(cats),
        "figures":           figures,
    }


if __name__ == "__main__":
    run_script(make_paga_figures)
