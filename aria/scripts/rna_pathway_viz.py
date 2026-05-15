"""
ARIA RNA Pathway Visualization
-------------------------------
Generates publication-grade figures for pathway enrichment results:
  1. ORA dotplot (per database) — log2(OR) vs term, colored by FDR, sized by gene count
  2. GSEA running sum plots (top 3 pathways per contrast) — via blitzgsea
  3. GSEA top_table figure — quick summary of top enriched terms

Patterns ported from the user's production scripts (ORA_mirna.py, gsea_*.py).
Outputs PNG (high-DPI, publication-ready) into figures_dir/.

This module is imported and called directly by rna_bulk_de.py — it does NOT
run as a standalone subprocess (so it shares the rna environment with DESeq2).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

log = logging.getLogger("aria.pathway_viz")


# Paper theme — white background, publication-ready
_DARK_BG    = "white"
_DARK_PANEL = "white"
_DARK_TEXT  = "#1e293b"
_DARK_MUTED = "#475569"


def make_ora_dotplot(pathways_list: list,
                       db_name: str,
                       contrast_name: str,
                       output_path: str,
                       top_n: int = 15) -> Optional[str]:
    """
    Generate ORA dotplot in the style of the user's ORA_mirna.py.

    Args:
        pathways_list: list of dicts with keys: term, padj, overlap, odds_ratio,
                       combined_score, genes
        db_name:       database label (e.g. "GO_BP", "KEGG", "Reactome")
        contrast_name: e.g. "treated vs control"
        output_path:   where to write PNG
        top_n:         max pathways to plot (default 15)

    Returns:
        output_path on success, None on failure.
    """
    if not pathways_list:
        return None

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.lines as mlines
        import seaborn as sns
        import numpy as np
        import pandas as pd
        import re

        df = pd.DataFrame(pathways_list)
        if df.empty:
            return None

        df = df.head(top_n).copy()

        # Derived columns matching ORA_mirna.py
        df["-log10(FDR)"] = -np.log10(df["padj"].clip(lower=1e-300)).round(1)

        # Genes column: in our schema it's already a list; count length
        df["Genes"] = df["genes"].apply(
            lambda g: len(g) if isinstance(g, list) else 1
        )

        # Clean term names (strip GO IDs, Reactome IDs)
        df["Term"] = (df["term"]
                      .str.replace(r"\s*\(GO:\d+\)$", "", regex=True)
                      .str.replace(r"\s+R-HSA-\d+$", "", regex=True)
                      .str.replace(r"\s+WP\d+$", "", regex=True))
        df["Term"] = df["Term"].str.slice(0, 70)  # truncate very long terms

        # log2 of OR (clip to avoid log(0))
        df["log2_OR"] = np.log2(df["odds_ratio"].clip(lower=0.01)).round(1)

        # Sort by Odds Ratio descending (most enriched at top after invert)
        df = df.sort_values("log2_OR", ascending=False)

        # ── Build the figure ─────────────────────────────────────────────
        sns.set(font="DejaVu Sans")
        sns.set_context("paper", font_scale=1.4)
        sns.set_style("whitegrid")

        fig, ax = plt.subplots(figsize=(9, max(6, 0.45 * len(df))))
        fig.patch.set_facecolor(_DARK_BG)
        ax.set_facecolor(_DARK_PANEL)

        scatter = sns.scatterplot(
            data=df,
            x="log2_OR",
            y="Term",
            hue="-log10(FDR)",
            palette="flare_r",
            size="Genes",
            sizes=(80, 500),
            legend=False,
            ax=ax,
        )

        # Color bar for FDR
        norm = plt.Normalize(df["-log10(FDR)"].min(), df["-log10(FDR)"].max())
        sm = plt.cm.ScalarMappable(cmap="flare_r", norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation="vertical",
                            fraction=0.025, pad=0.07, anchor=(1.0, 0.65))
        cbar.set_label("-log10(FDR)", fontsize=12, color=_DARK_TEXT)
        cbar.ax.tick_params(colors=_DARK_MUTED)
        cbar.outline.set_edgecolor(_DARK_MUTED)

        # Size legend
        gene_sizes = sorted(df["Genes"].unique())
        # Cap legend entries to 5 for readability
        if len(gene_sizes) > 5:
            step = max(1, len(gene_sizes) // 5)
            gene_sizes = gene_sizes[::step][:5]
        marker_sizes = np.linspace(80, 500, len(gene_sizes))
        legend_handles = [
            mlines.Line2D([], [], color=_DARK_TEXT, marker="o",
                          markersize=np.sqrt(s), label=f"{g}",
                          linestyle="None")
            for g, s in zip(gene_sizes, marker_sizes)
        ]
        leg = ax.legend(
            handles=legend_handles,
            title="Gene count",
            title_fontsize=11,
            fontsize=10,
            loc="upper left",
            bbox_to_anchor=(1.02, 0.4),
            borderaxespad=0,
            frameon=False,
            labelspacing=0.8,
            labelcolor=_DARK_TEXT,
        )
        leg.get_title().set_color(_DARK_TEXT)

        # Labels and title
        ax.set_xlabel("log₂(Odds Ratio)", fontsize=13, color=_DARK_TEXT)
        ax.set_ylabel("", fontsize=13)
        ax.set_title(f"{db_name} Enrichment — {contrast_name}",
                     fontsize=14, fontweight="bold", color="#0f172a", pad=12)

        ax.tick_params(colors=_DARK_MUTED, labelsize=10)
        for spine in ax.spines.values():
            spine.set_edgecolor("#e2e8f0")

        plt.tight_layout()
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=200, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        plt.close(fig)

        return str(output_path)

    except Exception as e:
        log.warning(f"ORA dotplot failed for {db_name}/{contrast_name}: {e}")
        return None


def make_gsea_running_sums(de_results_df,
                              symbol_map: dict,
                              contrast_name: str,
                              output_dir: str,
                              gene_set_library: str = "GO_Biological_Process_2023",
                              organism: str = "human",
                              top_n_plots: int = 3) -> dict:
    """
    Run GSEA via blitzgsea on a ranked gene list (log2FC), then plot the
    top N enriched pathways as running sum plots.

    Pattern from user's gsea_results.py and gsea_graph.py.

    Args:
        de_results_df: DESeq2 results DataFrame (index = gene_id,
                       must have columns log2FoldChange and padj)
        symbol_map:    {ensembl_id (no version): hgnc_symbol}
        contrast_name: e.g. "treated vs control"
        output_dir:    figures dir for this contrast
        gene_set_library: blitzgsea library name
        organism:      "human" or "mouse" (matches gseapy convention)
        top_n_plots:   how many running sum plots to make

    Returns:
        dict with:
          'gsea_table':    path to GSEA results CSV (sorted by NES)
          'top_table_fig': path to summary table figure
          'running_sums':  list of paths to running sum plots
          'n_pathways':    number of significant pathways (FDR<0.25 standard)
    """
    out = {
        "gsea_table":    None,
        "top_table_fig": None,
        "running_sums":  [],
        "n_pathways":    0,
    }

    try:
        import blitzgsea as blitz
        import pandas as pd
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # ── Build the ranked signature ───────────────────────────────────
        # blitzgsea expects DataFrame with 2 cols: [symbol, score]
        # where score is log2FC (positive = upregulated).
        # Convert IDs → symbols, drop unmapped, dedupe (keep largest |log2FC|)
        df = de_results_df.dropna(subset=["log2FoldChange"]).copy()
        df["clean_id"] = df.index.astype(str).str.split(".").str[0]
        df["symbol"] = df["clean_id"].map(symbol_map)
        df = df.dropna(subset=["symbol"])

        if df.empty:
            log.warning(f"GSEA[{contrast_name}]: no genes after symbol mapping")
            return out

        # Dedupe: if multiple Ensembl IDs map to same symbol, keep the one
        # with largest |log2FC|
        df["abs_lfc"] = df["log2FoldChange"].abs()
        df = df.sort_values("abs_lfc", ascending=False)
        df = df.drop_duplicates(subset=["symbol"], keep="first")

        signature = pd.DataFrame({
            0: df["symbol"].values,
            1: df["log2FoldChange"].values,
        })
        signature.columns = ["gene", "score"]

        log.info(f"GSEA[{contrast_name}]: {len(signature)} ranked genes")

        # ── Fetch gene set library ───────────────────────────────────────
        try:
            library = blitz.enrichr.get_library(gene_set_library)
        except Exception as e:
            log.warning(f"GSEA[{contrast_name}]: library fetch failed: {e}")
            return out

        # ── Run GSEA ─────────────────────────────────────────────────────
        try:
            result = blitz.gsea(signature, library)
        except Exception as e:
            log.warning(f"GSEA[{contrast_name}]: gsea() failed: {e}")
            return out

        if result is None or len(result) == 0:
            return out

        # Save the full GSEA table
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        table_path = Path(output_dir) / "gsea_results.csv"
        result_sorted = result.sort_values("fdr").copy()
        result_sorted.to_csv(table_path)
        out["gsea_table"] = str(table_path)
        out["n_pathways"] = int((result_sorted["fdr"] < 0.25).sum())

        # ── Top-N running sum plots ──────────────────────────────────────
        top = result_sorted.head(top_n_plots)
        for i, (term, _) in enumerate(top.iterrows()):
            try:
                fig = blitz.plot.running_sum(
                    signature, term, library,
                    result=result, compact=False,
                )
                rank_path = Path(output_dir) / f"gsea_running_sum_{i+1}.png"
                fig.savefig(rank_path, dpi=200, bbox_inches="tight",
                            facecolor="white")
                plt.close(fig)
                out["running_sums"].append(str(rank_path))
            except Exception as e:
                log.warning(f"GSEA running_sum {i+1} failed: {e}")
                continue

        # ── Top table figure (summary) ───────────────────────────────────
        try:
            fig_table = blitz.plot.top_table(signature, library, result, n=15)
            tt_path = Path(output_dir) / "gsea_top_table.png"
            fig_table.savefig(tt_path, dpi=200, bbox_inches="tight",
                              facecolor="white")
            plt.close(fig_table)
            out["top_table_fig"] = str(tt_path)
        except Exception as e:
            log.warning(f"GSEA top_table failed: {e}")

        return out

    except ImportError:
        log.warning("blitzgsea not installed — skipping GSEA visualization")
        return out
    except Exception as e:
        log.warning(f"GSEA failed for {contrast_name}: {e}")
        return out


def export_de_table(de_results_df,
                      symbol_map: dict,
                      contrast_name: str,
                      output_path: str,
                      padj_thr: float = 0.05) -> Optional[str]:
    """
    Export significant DE genes as a TSV supplementary table.
    Includes: gene_id, gene_symbol, log2FoldChange, padj, baseMean.
    Sorted by padj ascending.
    """
    try:
        import pandas as pd

        df = de_results_df.dropna(subset=["padj"]).copy()
        df = df[df["padj"] < padj_thr].sort_values("padj")
        if df.empty:
            return None

        df["gene_id"] = df.index.astype(str)
        df["clean_id"] = df["gene_id"].str.split(".").str[0]
        df["gene_symbol"] = df["clean_id"].map(symbol_map).fillna("")

        cols = ["gene_id", "gene_symbol"]
        for c in ["baseMean", "log2FoldChange", "lfcSE",
                   "stat", "pvalue", "padj"]:
            if c in df.columns:
                cols.append(c)

        df = df[cols]
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, sep="\t", index=False)
        return str(output_path)
    except Exception as e:
        log.warning(f"DE table export failed for {contrast_name}: {e}")
        return None


def export_pathways_table(pathways_dict: dict,
                            contrast_name: str,
                            output_path: str) -> Optional[str]:
    """
    Export pathway enrichment results as a TSV (all databases combined).
    Columns: database, term, padj, overlap, odds_ratio, combined_score, genes
    """
    try:
        import pandas as pd

        rows = []
        for db, terms in pathways_dict.items():
            if not isinstance(terms, list):
                continue
            for t in terms:
                rows.append({
                    "database":       db,
                    "term":           t.get("term", ""),
                    "padj":           t.get("padj", 1),
                    "overlap":        t.get("overlap", ""),
                    "odds_ratio":     t.get("odds_ratio", 1),
                    "combined_score": t.get("combined_score", 0),
                    "genes":          ";".join(t.get("genes", []))
                                       if isinstance(t.get("genes"), list)
                                       else str(t.get("genes", "")),
                })
        if not rows:
            return None

        df = pd.DataFrame(rows).sort_values(["database", "padj"])
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, sep="\t", index=False)
        return str(output_path)
    except Exception as e:
        log.warning(f"Pathways table export failed for {contrast_name}: {e}")
        return None
