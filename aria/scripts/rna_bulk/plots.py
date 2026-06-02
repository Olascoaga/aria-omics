"""Bulk RNA-seq plotting helpers (P2-8 split from rna_bulk_de.py).

Behavior-preserving extraction; re-exported from `aria.scripts.rna_bulk_de`.
Matplotlib/numpy/pandas are imported lazily inside each function."""

from __future__ import annotations

from pathlib import Path


_P = dict(
    fig     = "white",
    ax      = "white",
    text    = "#1e293b",
    muted   = "#475569",
    dim     = "#94a3b8",
    border  = "#e2e8f0",
    title   = "#0f172a",
    up      = "#dc2626",   # upregulated  (red  — convention)
    down    = "#2563eb",   # downregulated (blue — convention)
    ns      = "#d1d5db",   # non-significant (light gray)
    ref     = "#94a3b8",   # threshold reference lines
    annot   = "#374151",   # gene label annotations
    palette = ["#1d4ed8", "#dc2626", "#059669", "#d97706",
               "#7c3aed", "#db2777", "#0891b2"],
)


def _plot_pca_mds(vst_variable, metadata, output_dir: str,
                    warnings: list) -> tuple:
    """
    Generate PCA + MDS plots from a VST-transformed, variable-filtered matrix.

    Both plots operate on the SAME data matrix → apples-to-apples comparison.
    PCA: linear orthogonal components capturing variance.
    MDS: preserves Euclidean distances (non-linear relationships).

    Returns (pca_path, mds_path, pca_coords, pca_variance_pct).
    """
    try:
        import numpy as np
        import pandas as pd
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sklearn.decomposition import PCA
        from sklearn.manifold import MDS
        from scipy.spatial.distance import pdist, squareform

        # Samples × genes for PCA/MDS (standard orientation)
        X = vst_variable.T.values

        # Identify the condition column
        condition_col = None
        for c in metadata.columns:
            if c not in ("sample", "batch", "replicate"):
                condition_col = c
                break
        conditions = (metadata[condition_col].astype(str)
                       if condition_col else
                       pd.Series(["all"] * X.shape[0], index=metadata.index))

        # ── PCA ──────────────────────────────────────────────────────────
        # NO StandardScaler — VST already stabilizes variance.
        # PCA on genes × samples with centered (not z-scored) values.
        pca     = PCA(n_components=min(5, X.shape[0] - 1))
        coords  = pca.fit_transform(X)
        var_pct = pca.explained_variance_ratio_ * 100

        # ── MDS (classical / metric) ─────────────────────────────────────
        # Euclidean distance on VST — same convention as edgeR's plotMDS.
        dist_matrix = squareform(pdist(X, metric="euclidean"))
        # Use future defaults explicitly to avoid FutureWarnings on sklearn 1.9+
        import warnings as _w
        with _w.catch_warnings():
            _w.simplefilter("ignore", FutureWarning)
            try:
                # sklearn 1.9+ signature
                mds = MDS(n_components=2, metric=True,
                           normalized_stress="auto", random_state=42,
                           n_init=4, dissimilarity="precomputed")
            except TypeError:
                mds = MDS(n_components=2, dissimilarity="precomputed",
                           random_state=42)
            mds_coords = mds.fit_transform(dist_matrix)

        # ── Plot both side by side in the dark theme ─────────────────────
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.patch.set_facecolor(_P["fig"])

        groups = sorted(conditions.unique())
        color_map = {g: _P["palette"][i % len(_P["palette"])]
                     for i, g in enumerate(groups)}

        for ax, coords_plot, title, xlbl, ylbl in [
            (axes[0], coords[:, :2], "PCA (VST, top variable PC genes)",
             f"PC1 ({var_pct[0]:.1f}%)", f"PC2 ({var_pct[1]:.1f}%)"),
            (axes[1], mds_coords,    "MDS (VST, Euclidean distance)",
             "MDS1", "MDS2"),
        ]:
            ax.set_facecolor(_P["ax"])
            for g in groups:
                mask = (conditions.values == g)
                ax.scatter(coords_plot[mask, 0], coords_plot[mask, 1],
                           s=120, c=color_map[g], label=g,
                           edgecolor="white", linewidth=0.8, alpha=0.9)
            # Sample labels
            for i, s in enumerate(metadata.index):
                ax.annotate(s, (coords_plot[i, 0], coords_plot[i, 1]),
                            fontsize=8, color=_P["annot"],
                            xytext=(5, 5), textcoords="offset points")
            ax.set_title(title, color=_P["title"], fontsize=11, fontweight="bold")
            ax.set_xlabel(xlbl, color=_P["muted"])
            ax.set_ylabel(ylbl, color=_P["muted"])
            ax.tick_params(colors=_P["muted"])
            ax.legend(facecolor=_P["fig"], edgecolor=_P["border"],
                      labelcolor=_P["text"], fontsize=9, loc="best")
            for spine in ax.spines.values():
                spine.set_edgecolor(_P["border"])

        plt.tight_layout()
        combined_path = str(Path(output_dir) / "pca_mds.svg")
        plt.savefig(combined_path, format="svg",
                     facecolor=_P["fig"])
        plt.close(fig)

        # Also save individual SVGs for manuscript use
        pca_path = _save_single_dr_plot(
            coords[:, :2], conditions, metadata.index, color_map,
            f"PCA (VST, top variable PC genes)",
            f"PC1 ({var_pct[0]:.1f}%)", f"PC2 ({var_pct[1]:.1f}%)",
            str(Path(output_dir) / "pca.svg"),
        )
        mds_path = _save_single_dr_plot(
            mds_coords, conditions, metadata.index, color_map,
            "MDS (VST, Euclidean distance)",
            "MDS1", "MDS2",
            str(Path(output_dir) / "mds.svg"),
        )

        return pca_path, mds_path, coords, [round(float(v), 3)
                                               for v in var_pct[:2]]

    except Exception as e:
        warnings.append(f"PCA/MDS plotting failed: {e}")
        return None, None, None, []

def _save_single_dr_plot(coords, conditions, sample_names, color_map,
                          title, xlbl, ylbl, output_path):
    """Helper — single-axis DR plot for manuscript figures."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 5))
        fig.patch.set_facecolor(_P["fig"])
        ax.set_facecolor(_P["ax"])
        for g in sorted(color_map.keys()):
            mask = (conditions.values == g)
            ax.scatter(coords[mask, 0], coords[mask, 1],
                       s=140, c=color_map[g], label=g,
                       edgecolor="white", linewidth=0.8, alpha=0.9)
        for i, s in enumerate(sample_names):
            ax.annotate(s, (coords[i, 0], coords[i, 1]),
                        fontsize=8, color=_P["annot"],
                        xytext=(5, 5), textcoords="offset points")
        ax.set_title(title, color=_P["title"], fontsize=12, fontweight="bold")
        ax.set_xlabel(xlbl, color=_P["muted"])
        ax.set_ylabel(ylbl, color=_P["muted"])
        ax.tick_params(colors=_P["muted"])
        ax.legend(facecolor=_P["fig"], edgecolor=_P["border"],
                  labelcolor=_P["text"], fontsize=9, loc="best")
        for spine in ax.spines.values():
            spine.set_edgecolor(_P["border"])
        plt.tight_layout()
        plt.savefig(output_path, format="svg",
                     facecolor=_P["fig"])
        plt.close(fig)
        return output_path
    except Exception:
        return None

def _generate_plots(de_result: dict, sample_qc: dict,
                    counts_filt, metadata, design_factor: str,
                    output_dir: str, padj_thr: float,
                    lfc_thr: float,
                    title_suffix: str = "") -> dict:
    """Generate volcano, heatmap, and sample PCA plots."""
    plots = {}

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        results_df = de_result.get("results")
        if results_df is not None and not results_df.empty:
            # ── Volcano plot ───────────────────────────────────────────
            fig, ax = plt.subplots(figsize=(8, 6))
            fig.patch.set_facecolor(_P["fig"])
            ax.set_facecolor(_P["ax"])

            lfc = results_df["log2FoldChange"].clip(-8, 8)
            neg_log_p = -np.log10(results_df["padj"].clip(1e-300, 1) + 1e-300)

            # Color by significance
            sig_mask = results_df["padj"] < padj_thr
            up_mask   = sig_mask & (results_df["log2FoldChange"] > 0)
            down_mask = sig_mask & (results_df["log2FoldChange"] < 0)

            ax.scatter(lfc[~sig_mask], neg_log_p[~sig_mask],
                       color=_P["ns"], alpha=0.5, s=8, linewidths=0)
            ax.scatter(lfc[up_mask], neg_log_p[up_mask],
                       color=_P["up"], alpha=0.8, s=12, linewidths=0)
            ax.scatter(lfc[down_mask], neg_log_p[down_mask],
                       color=_P["down"], alpha=0.8, s=12, linewidths=0)

            # Reference lines
            ax.axhline(-np.log10(padj_thr), color=_P["ref"],
                       linestyle="--", alpha=0.6, linewidth=0.8)
            ax.axvline( lfc_thr, color=_P["ref"],
                       linestyle="--", alpha=0.6, linewidth=0.8)
            ax.axvline(-lfc_thr, color=_P["ref"],
                       linestyle="--", alpha=0.6, linewidth=0.8)

            # Labels for top genes
            top_sig = results_df[sig_mask].nsmallest(10, "padj")
            for gene, row in top_sig.iterrows():
                ax.annotate(
                    str(gene)[:12],
                    xy=(float(row["log2FoldChange"]),
                        float(-np.log10(row["padj"] + 1e-300))),
                    xytext=(3, 3), textcoords="offset points",
                    fontsize=6, color=_P["annot"], alpha=0.9,
                )

            ax.set_xlabel("log₂ Fold Change", color=_P["muted"])
            ax.set_ylabel("-log₁₀ adjusted p-value", color=_P["muted"])
            n_up   = de_result.get("n_up",   0)
            n_down = de_result.get("n_down", 0)
            title_text = (
                f"Differential Expression"
                + (f" — {title_suffix}" if title_suffix else "")
                + f"\n{n_up} up (red)  {n_down} down (blue)"
            )
            ax.set_title(title_text, color=_P["title"], fontsize=11,
                         fontweight="bold")
            ax.tick_params(colors=_P["muted"])
            for spine in ax.spines.values():
                spine.set_edgecolor(_P["border"])

            volcano_path = str(Path(output_dir) / "volcano.svg")
            plt.tight_layout()
            plt.savefig(volcano_path, format="svg",
                        facecolor=_P["fig"])
            plt.close()
            plots["volcano"] = volcano_path

        # ── Heatmaps of top DE genes (two complementary views) ──────────
        # (a) Top 50 by padj — most statistically confident genes
        # (b) Top 50 by |log2FC| — largest effect sizes (user-requested v3.8)
        #
        # Input data: VST if available (sample_qc provides it), else
        # log2(counts+1). Then row z-score for visualization.
        # Symbol annotation: use symbol_map to replace ENSG with HGNC.
        de_results_df = de_result.get("results")
        vst_matrix    = sample_qc.get("vst_matrix")   # may be None
        symbol_map    = sample_qc.get("_symbol_map", {})  # passed via QC dict

        if de_results_df is not None and not de_results_df.empty:
            # Intersect with gene matrix (only plot genes we have data for)
            available_all = [g for g in de_results_df.index
                              if (vst_matrix is not None and g in vst_matrix.index)
                              or (counts_filt is not None and g in counts_filt.index)]

            if available_all:
                # (a) Top 50 by padj — already sig_genes-ordered by padj asc
                top_padj = [g for g in de_result.get("sig_genes", [])[:50]
                            if g in available_all]
                hm_a = _plot_heatmap(
                    genes=top_padj,
                    vst_matrix=vst_matrix,
                    counts_filt=counts_filt,
                    symbol_map=symbol_map,
                    title_prefix="Top 50 DE genes by padj",
                    title_suffix=title_suffix,
                    output_path=str(Path(output_dir) / "heatmap_padj.svg"),
                )
                if hm_a:
                    plots["heatmap"]      = hm_a   # backward compat
                    plots["heatmap_padj"] = hm_a

                # (b) Top 50 by |log2FC| (NEW in v3.8)
                lfc_sorted = de_results_df.dropna(subset=["padj","log2FoldChange"])
                lfc_sorted = lfc_sorted[lfc_sorted["padj"] < padj_thr]
                lfc_sorted = lfc_sorted.assign(
                    abs_lfc=lfc_sorted["log2FoldChange"].abs()
                ).sort_values("abs_lfc", ascending=False)
                top_lfc = [g for g in lfc_sorted.index[:50]
                            if g in available_all]
                hm_b = _plot_heatmap(
                    genes=top_lfc,
                    vst_matrix=vst_matrix,
                    counts_filt=counts_filt,
                    symbol_map=symbol_map,
                    title_prefix="Top 50 DE genes by |log2FC|",
                    title_suffix=title_suffix,
                    output_path=str(Path(output_dir) / "heatmap_lfc.svg"),
                )
                if hm_b:
                    plots["heatmap_lfc"] = hm_b

    except Exception as e:
        plots["error"] = str(e)

    # PCA and MDS were already generated in _sample_qc (v3.8: VST-based)
    if sample_qc.get("pca_plot"):
        plots["pca"] = sample_qc["pca_plot"]
    if sample_qc.get("mds_plot"):
        plots["mds"] = sample_qc["mds_plot"]

    return plots

def _plot_heatmap(genes: list, vst_matrix, counts_filt,
                    symbol_map: dict,
                    title_prefix: str, title_suffix: str,
                    output_path: str) -> str | None:
    """
    Plot a row-z-score heatmap of the given genes.

    Data source priority:
      1. VST matrix if provided (homoscedastic, lib-size corrected)
      2. log2(counts_filt+1) as fallback

    Row labels are replaced with HGNC symbols when available.
    Dark theme matching the rest of the report.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd

        if not genes:
            return None

        if vst_matrix is not None:
            available = [g for g in genes if g in vst_matrix.index]
            if not available:
                return None
            hm_data = vst_matrix.loc[available]
        elif counts_filt is not None:
            available = [g for g in genes if g in counts_filt.index]
            if not available:
                return None
            hm_data = np.log2(counts_filt.loc[available].astype(float) + 1)
        else:
            return None

        # Row z-score across samples
        row_mean = hm_data.mean(axis=1)
        row_std  = hm_data.std(axis=1).replace(0, 1)
        hm_z     = ((hm_data.T - row_mean) / row_std).T

        # Gene labels: symbols when available
        row_labels = []
        for g in available:
            clean = str(g).split(".")[0]
            sym   = (symbol_map or {}).get(clean, "")
            row_labels.append(sym if sym else clean)

        n_genes = len(available)
        fig, ax = plt.subplots(
            figsize=(
                max(6, len(hm_data.columns) * 0.6),
                max(5, n_genes * 0.18),
            )
        )
        fig.patch.set_facecolor(_P["fig"])

        im = ax.imshow(hm_z, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
        ax.set_xticks(range(len(hm_data.columns)))
        ax.set_xticklabels(hm_data.columns, rotation=45,
                            ha="right", fontsize=9, color=_P["text"])
        ax.set_yticks(range(n_genes))
        ax.set_yticklabels(row_labels, fontsize=7, color=_P["text"])
        ax.set_facecolor(_P["ax"])

        cbar = plt.colorbar(im, ax=ax, label="Z-score (VST or log2)",
                             fraction=0.025, pad=0.04)
        cbar.ax.tick_params(colors=_P["muted"])
        cbar.set_label("Z-score (VST or log2)", color=_P["text"])

        title = title_prefix + (f" — {title_suffix}" if title_suffix else "")
        ax.set_title(title, color=_P["title"], fontsize=11,
                     fontweight="bold", pad=10)

        plt.tight_layout()
        plt.savefig(output_path, format="svg",
                     facecolor=_P["fig"])
        plt.close(fig)
        return output_path
    except Exception as e:
        log.warning(f"Heatmap failed ({title_prefix}): {e}")
        return None

def _plot_sample_pca(coords, samples, metadata, var_exp: list,
                      output_dir: str) -> str | None:
    """Save sample PCA plot colored by condition."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        fig, ax = plt.subplots(figsize=(7, 5))
        fig.patch.set_facecolor(_P["fig"])
        ax.set_facecolor(_P["ax"])

        COLORS = _P["palette"]

        # Get condition for each sample
        condition_col = None
        for col in metadata.columns:
            if col not in ("sample", "batch", "replicate"):
                condition_col = col
                break

        groups = metadata[condition_col].values if condition_col else \
                 ["unknown"] * len(samples)
        unique = sorted(set(groups))

        for i, grp in enumerate(unique):
            mask = [g == grp for g in groups]
            ax.scatter(
                coords[mask, 0], coords[mask, 1],
                label=grp, color=COLORS[i % len(COLORS)],
                s=80, alpha=0.85, edgecolors=_P["fig"], linewidths=0.5,
            )

        for j, s in enumerate(samples):
            ax.annotate(s[:10], (coords[j, 0], coords[j, 1]),
                        xytext=(4, 4), textcoords="offset points",
                        fontsize=6, color=_P["muted"])

        pct1 = round(float(var_exp[0]) * 100, 1) if len(var_exp) > 0 else 0
        pct2 = round(float(var_exp[1]) * 100, 1) if len(var_exp) > 1 else 0
        ax.set_xlabel(f"PC1 ({pct1}%)", color=_P["muted"])
        ax.set_ylabel(f"PC2 ({pct2}%)", color=_P["muted"])
        ax.set_title("Sample PCA", color=_P["title"], fontweight="bold")
        ax.tick_params(colors=_P["muted"])
        for spine in ax.spines.values():
            spine.set_edgecolor(_P["border"])
        ax.legend(fontsize=8, facecolor=_P["fig"],
                  labelcolor=_P["text"], edgecolor=_P["border"])

        pca_path = str(Path(output_dir) / "sample_pca.svg")
        plt.tight_layout()
        plt.savefig(pca_path, format="svg",
                    facecolor=_P["fig"])
        plt.close()
        return pca_path

    except Exception:
        return None
