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

from aria.scripts.rna_bulk.gtf_io import _strip_gene_version

log = logging.getLogger("aria.pathway_viz")


# Paper theme — white background, publication-ready
_DARK_BG    = "white"
_DARK_PANEL = "white"
_DARK_TEXT  = "#1e293b"
_DARK_MUTED = "#475569"


def _ranked_signature_frame(de_results_df, symbol_map: Optional[dict] = None):
    """
    Build the two-column preranked GSEA signature.

    Ensembl IDs require a symbol map. Already-symbolic gene IDs are kept when
    no map exists or when a mapped lookup misses a non-Ensembl token.
    """
    import numpy as np
    import pandas as pd

    symbol_map = symbol_map or {}
    df = de_results_df.copy()
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df = df.loc[np.isfinite(df["log2FoldChange"])].copy()
    df["clean_id"] = df.index.astype(str).map(_strip_gene_version)

    def resolve_symbol(gene_id: str):
        symbol = symbol_map.get(gene_id)
        if symbol:
            return symbol
        if not gene_id.startswith(("ENSG", "ENSMUSG", "ENS")):
            return gene_id
        return None

    df["symbol"] = df["clean_id"].map(resolve_symbol)
    df = df.dropna(subset=["symbol"])
    if df.empty:
        return pd.DataFrame(columns=["gene", "score"])

    # Dedupe: if multiple Ensembl IDs map to the same symbol, keep the one
    # with largest |log2FC|.
    df["abs_lfc"] = df["log2FoldChange"].abs()
    df = df.sort_values("abs_lfc", ascending=False)
    df = df.drop_duplicates(subset=["symbol"], keep="first")

    signature = pd.DataFrame({
        "gene": df["symbol"].values,
        "score": df["log2FoldChange"].values,
    })
    return signature


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

        # Derived columns matching ORA_mirna.py
        padj_num = pd.to_numeric(df["padj"], errors="coerce")
        odds_num = pd.to_numeric(df["odds_ratio"], errors="coerce")
        valid = (np.isfinite(padj_num) & padj_num.between(0, 1)
                 & np.isfinite(odds_num) & (odds_num > 0))
        df = df.loc[valid].head(top_n).copy()
        if df.empty:
            return None
        padj_num = padj_num.loc[df.index]
        positive = padj_num > 0
        neglog = -np.log10(padj_num.where(positive))
        cap = max(1.0, float(neglog[positive].max()) * 1.15) if positive.any() else 10.0
        df["-log10(FDR)"] = neglog.fillna(cap)
        has_underflow = bool((~positive).any())

        # Gene count: prefer the true overlap numerator (e.g. "12/134" -> 12)
        # because the display `genes` list is capped at 10 entries in ORA output.
        def _gene_count(row) -> int:
            overlap = row.get("overlap")
            if isinstance(overlap, str) and "/" in overlap:
                try:
                    k = int(overlap.split("/", 1)[0].strip())
                    if k > 0:
                        return k
                except (TypeError, ValueError):
                    pass
            genes = row.get("genes")
            if isinstance(genes, list):
                return max(1, len(genes))
            if isinstance(genes, str) and genes.strip():
                return max(1, len([g for g in genes.split(";") if g.strip()]))
            return 1

        df["Genes"] = df.apply(_gene_count, axis=1)

        # Clean term names (strip GO IDs, Reactome IDs)
        df["Term"] = (df["term"]
                      .str.replace(r"\s*\(GO:\d+\)$", "", regex=True)
                      .str.replace(r"\s+R-HSA-\d+$", "", regex=True)
                      .str.replace(r"\s+WP\d+$", "", regex=True))
        df["Term"] = df["Term"].str.slice(0, 70)  # truncate very long terms
        # Keep colliding display names on distinct rows, including exact duplicates.
        used_labels = set()
        labels = []
        for label in df["Term"].fillna("Unnamed term"):
            candidate = label
            suffix = 2
            while candidate in used_labels:
                candidate = f"{label} [{suffix}]"
                suffix += 1
            used_labels.add(candidate)
            labels.append(candidate)
        df["Term"] = labels

        df["log2_OR"] = np.log2(odds_num.loc[df.index]).round(1)

        # Sort by Odds Ratio descending (most enriched at top after invert)
        df = df.sort_values("log2_OR", ascending=False)

        # ── Build the figure ─────────────────────────────────────────────
        sns.set(font="DejaVu Sans")
        sns.set_context("paper", font_scale=1.4)
        sns.set_style("whitegrid")

        fig, ax = plt.subplots(figsize=(9, max(6, 0.45 * len(df))))
        fig.patch.set_facecolor(_DARK_BG)
        ax.set_facecolor(_DARK_PANEL)

        # Color normalization shared between scatter and colorbar
        vmin = float(df["-log10(FDR)"].min())
        vmax = float(df["-log10(FDR)"].max())
        if not np.isfinite(vmin) or not np.isfinite(vmax):
            vmin, vmax = 0.0, 1.0
        elif vmin == vmax:
            vmin = max(0.0, vmin - 0.5)
            vmax = vmax + 0.5
        norm = plt.Normalize(vmin, vmax)

        scatter = sns.scatterplot(
            data=df,
            x="log2_OR",
            y="Term",
            hue="-log10(FDR)",
            hue_norm=norm,
            palette="flare_r",
            size="Genes",
            sizes=(80, 500),
            legend=False,
            ax=ax,
        )

        # Color bar for FDR
        sm = plt.cm.ScalarMappable(cmap="flare_r", norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation="vertical",
                            fraction=0.025, pad=0.07, anchor=(1.0, 0.65))
        cbar_label = "-log10(FDR)"
        if has_underflow:
            cbar_label += f" (FDR=0 shown at {cap:.2f}; underflow)"
        cbar.set_label(cbar_label, fontsize=12, color=_DARK_TEXT)
        cbar.ax.tick_params(colors=_DARK_MUTED)
        cbar.outline.set_edgecolor(_DARK_MUTED)

        # Size legend (scaled linearly on exact gene-count value to match scatter)
        gene_sizes = sorted(int(g) for g in df["Genes"].unique())
        # Cap legend entries to 5 for readability while preserving min & max
        if len(gene_sizes) > 5:
            idx = np.linspace(0, len(gene_sizes) - 1, 5).round().astype(int)
            gene_sizes = sorted({gene_sizes[i] for i in idx})
        g_min, g_max = min(gene_sizes), max(gene_sizes)
        if g_max == g_min:
            marker_sizes = [250.0 for _ in gene_sizes]
        else:
            marker_sizes = [
                80.0 + (g - g_min) / (g_max - g_min) * (500.0 - 80.0)
                for g in gene_sizes
            ]
        legend_handles = [
            mlines.Line2D([], [], color=_DARK_TEXT, marker="o",
                          markersize=np.sqrt(s), label=f"{g}",
                          linestyle="None")
            for g, s in zip(gene_sizes, marker_sizes, strict=True)
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
        if Path(output_path).suffix.lower() == ".png":
            fig.savefig(Path(output_path).with_suffix(".svg"), bbox_inches="tight",
                        facecolor=fig.get_facecolor())
        plt.close(fig)

        return str(output_path)

    except Exception as e:
        log.warning(f"ORA dotplot failed for {db_name}/{contrast_name}: {e}")
        return None


def _partition_gsea_results(result):
    """Split a blitzgsea result frame into (stable_sorted, unstable).

    ``blitzgsea`` estimates p-values and NES via parametric Gamma anchor
    interpolation, which can occasionally emit ``NES = ±inf`` (and ``fdr = 0.0``)
    for oversized or tail-saturated gene sets. Keeping non-finite rows out of
    counts, running-sum plots, and ``top_table`` prevents ``NES=inf`` from
    taking rank #1 while preserving those rows at the tail of
    ``gsea_results.csv`` so ``BulkRnaNarrator`` can honestly disclose them.
    """
    import numpy as np
    import pandas as pd

    df = result.copy()
    finite_mask = pd.Series(True, index=df.index)
    for col in ("nes", "fdr", "es", "pval"):
        if col in df.columns:
            vals = pd.to_numeric(df[col], errors="coerce")
            finite_mask &= np.isfinite(vals)

    stable = df.loc[finite_mask].copy()
    unstable = df.loc[~finite_mask].copy()

    if not stable.empty:
        stable["_sort_fdr"] = pd.to_numeric(stable["fdr"], errors="coerce")
        stable["_sort_abs_nes"] = -pd.to_numeric(
            stable["nes"], errors="coerce"
        ).abs()
        stable["_sort_term"] = stable.index.astype(str)
        stable = stable.sort_values(
            ["_sort_fdr", "_sort_abs_nes", "_sort_term"],
            ascending=[True, True, True],
            kind="mergesort",
        ).drop(columns=["_sort_fdr", "_sort_abs_nes", "_sort_term"])

    return stable, unstable


def _load_gsea_library(gene_set_library: str, blitz_mod) -> tuple[Optional[dict], str]:
    """Load a gene-set library for GSEA, honoring local GMTs and air-gap mode.

    Resolution order:
      1. Local versioned GMT matching ``gene_set_library`` (offline, deterministic).
      2. ``blitzgsea.enrichr.get_library(gene_set_library)`` when egress is allowed.
      3. Local ``GO_Biological_Process_2021`` fallback when air-gapped or offline
         and a GO Biological Process library was requested.

    All gene symbols are upper-cased so both human and mouse signatures match.
    """
    try:
        from aria.utils import ora as _ora
        loaded = _ora.load_local_library(gene_set_library)
        if loaded is not None:
            gene_sets, _ = loaded
            return (
                {k: [str(g).strip().upper() for g in v if str(g).strip()]
                 for k, v in gene_sets.items() if v},
                gene_set_library,
            )
    except Exception:
        pass

    egress_ok = True
    try:
        from aria.utils.privacy import egress_allowed
        egress_ok = bool(egress_allowed())
    except Exception:
        egress_ok = True

    if egress_ok:
        try:
            lib = blitz_mod.enrichr.get_library(gene_set_library)
            if lib:
                return (
                    {k: [str(g).strip().upper() for g in v if str(g).strip()]
                     for k, v in lib.items() if v},
                    gene_set_library,
                )
        except Exception as e:
            log.warning(f"GSEA library fetch failed for {gene_set_library}: {e}")

    # Offline / air-gapped fallback to staged local GO_Biological_Process_2021
    if gene_set_library.startswith("GO_Biological_Process"):
        try:
            from aria.utils import ora as _ora
            fallback_name = "GO_Biological_Process_2021"
            loaded = _ora.load_local_library(fallback_name)
            if loaded is not None:
                gene_sets, _ = loaded
                log.info(
                    f"GSEA using local offline fallback library {fallback_name}"
                )
                return (
                    {k: [str(g).strip().upper() for g in v if str(g).strip()]
                     for k, v in gene_sets.items() if v},
                    fallback_name,
                )
        except Exception:
            pass

    return None, gene_set_library


def _prepare_gsea_signature(signature):
    """Return a descending-ranked, deduplicated, upper-cased signature frame
    plus centered score vector and ``{gene: rank_index}`` lookup map matching
    ``blitzgsea.gsea(..., center=True)``."""
    import numpy as np
    import pandas as pd

    sig = signature.copy()
    sig["score"] = pd.to_numeric(sig["score"], errors="coerce")
    sig = sig.loc[np.isfinite(sig["score"])].copy()
    sig["gene"] = sig["gene"].astype(str).str.strip().str.upper()
    sig = sig[sig["gene"] != ""]
    sig = sig.sort_values("score", ascending=False, kind="mergesort")
    sig = sig.drop_duplicates(subset=["gene"], keep="first").reset_index(drop=True)
    scores_raw = sig["score"].to_numpy(dtype=float)
    scores_centered = scores_raw - float(np.mean(scores_raw)) if len(scores_raw) else scores_raw
    gene_to_rank = {g: idx for idx, g in enumerate(sig["gene"].tolist())}
    return sig, scores_raw, scores_centered, gene_to_rank


def _compute_running_sum(abs_scores, gene_to_rank: dict, gene_set) -> tuple:
    """Exact weighted GSEA running-sum curve matching ``blitzgsea.enrichment_score``."""
    import numpy as np

    n_total = len(abs_scores)
    hits = sorted({gene_to_rank[g] for g in gene_set if g in gene_to_rank})
    if not hits or len(hits) >= n_total:
        return np.zeros(n_total, dtype=float), 0.0, hits
    hit_ind = np.zeros(n_total, dtype=float)
    hit_ind[hits] = 1.0
    no_hit_ind = 1.0 - hit_ind
    sum_hits = float(np.sum(abs_scores[hits]))
    if sum_hits <= 0.0:
        return np.zeros(n_total, dtype=float), 0.0, hits
    norm_hit = 1.0 / sum_hits
    norm_miss = 1.0 / float(n_total - len(hits))
    rs = np.cumsum(hit_ind * abs_scores * norm_hit - no_hit_ind * norm_miss)
    peak_idx = int(np.argmax(np.abs(rs)))
    return rs, float(rs[peak_idx]), hits


def _fmt_fdr(val: float) -> str:
    if val < 1e-3:
        return f"{val:.1e}"
    if val < 0.01:
        return f"{val:.4f}"
    return f"{val:.3f}"


def _geneset_size(row, fallback: int) -> int:
    """Use observed hits when an optional gene-set size is missing/non-finite."""
    import math

    try:
        value = float(row.get("geneset_size", fallback))
        return int(value) if math.isfinite(value) and value >= 0 else fallback
    except (TypeError, ValueError, OverflowError):
        return fallback


def _render_gsea_running_sum(
    scores_raw,
    scores_centered,
    gene_to_rank: dict,
    term: str,
    gene_set,
    row,
    contrast_name: str,
    output_png: Path,
) -> str:
    """Render a publication-grade 3-panel GSEA running sum (300-DPI PNG + SVG)."""
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    abs_scores = np.abs(scores_centered)
    rs, es_calc, hits = _compute_running_sum(abs_scores, gene_to_rank, gene_set)
    n_genes = len(scores_raw)
    nes = float(row["nes"])
    fdr = float(row["fdr"])
    es_val = float(row["es"]) if "es" in row and np.isfinite(row["es"]) else es_calc
    gsize = _geneset_size(row, len(hits))

    is_pos = nes >= 0
    accent = "#dc2626" if is_pos else "#2563eb"
    peak_idx = int(np.argmax(np.abs(rs))) if len(rs) else 0

    fig = plt.figure(figsize=(8.2, 5.6), facecolor="white")
    gs = fig.add_gridspec(3, 1, height_ratios=[4.0, 0.75, 1.9], hspace=0.06)
    ax_es = fig.add_subplot(gs[0])
    ax_hit = fig.add_subplot(gs[1], sharex=ax_es)
    ax_rank = fig.add_subplot(gs[2], sharex=ax_es)

    # Panel 1: Running enrichment score
    x = np.arange(n_genes)
    ax_es.axhline(0.0, color="#94a3b8", lw=0.9, ls="-", zorder=1)
    ax_es.plot(x, rs, color=accent, lw=2.3, zorder=3)
    ax_es.fill_between(x, 0.0, rs, color=accent, alpha=0.10, zorder=2)
    ax_es.axvline(peak_idx, color=accent, lw=1.1, ls="--", alpha=0.75, zorder=2)
    ax_es.scatter([peak_idx], [rs[peak_idx]], color=accent, s=36, zorder=4)

    # Place stats badge in the corner opposite the peak so it never covers the curve
    badge_x = 0.98 if peak_idx < n_genes * 0.5 else 0.02
    badge_ha = "right" if badge_x > 0.5 else "left"
    badge_y = 0.94 if is_pos else 0.08
    badge_va = "top" if is_pos else "bottom"
    badge_txt = (
        f"NES = {nes:+.2f}   |   ES = {es_val:+.2f}\n"
        f"FDR = {_fmt_fdr(fdr)}   |   n = {gsize} genes"
    )
    ax_es.text(
        badge_x, badge_y, badge_txt,
        transform=ax_es.transAxes, ha=badge_ha, va=badge_va,
        fontsize=9.5, fontweight="bold", color=_DARK_TEXT,
        bbox={
            "boxstyle": "round,pad=0.45",
            "facecolor": "white",
            "edgecolor": "#cbd5e1",
            "alpha": 0.94,
            "linewidth": 0.9,
        },
        zorder=10,
    )
    ax_es.set_ylabel("Enrichment Score (ES)", fontsize=10.5, color=_DARK_TEXT)
    ax_es.grid(True, axis="y", color="#e2e8f0", lw=0.7)
    ax_es.tick_params(axis="x", labelbottom=False, length=0)
    ax_es.tick_params(axis="y", colors=_DARK_MUTED, labelsize=9.5)
    for sp in ("top", "right"):
        ax_es.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax_es.spines[sp].set_edgecolor("#cbd5e1")

    title_wrapped = textwrap.fill(str(term), width=64)
    ax_es.set_title(
        f"{title_wrapped}\n{contrast_name} · Pre-ranked GSEA (log₂FC)",
        fontsize=11.5, fontweight="bold", color="#0f172a", pad=10,
    )

    # Panel 2: Gene-hit barcode track
    ax_hit.set_facecolor("#f8fafc")
    if hits:
        ax_hit.vlines(hits, ymin=0.08, ymax=0.92, color=accent, lw=0.95, alpha=0.85)
    ax_hit.set_ylim(0, 1)
    ax_hit.set_yticks([])
    ax_hit.tick_params(axis="x", labelbottom=False, length=0)
    ax_hit.set_ylabel("Hits", fontsize=9, color=_DARK_MUTED, rotation=0,
                      labelpad=18, va="center")
    for sp in ax_hit.spines.values():
        sp.set_edgecolor("#cbd5e1")
        sp.set_linewidth(0.8)

    # Panel 3: Ranked list metric (log2FC)
    pos_mask = scores_raw >= 0
    neg_mask = scores_raw < 0
    ax_rank.axhline(0.0, color="#64748b", lw=0.8)
    ax_rank.fill_between(x, 0.0, scores_raw, where=pos_mask,
                         color="#fca5a5", alpha=0.55, interpolate=True)
    ax_rank.fill_between(x, 0.0, scores_raw, where=neg_mask,
                         color="#93c5fd", alpha=0.55, interpolate=True)
    ax_rank.plot(x, scores_raw, color="#334155", lw=0.9)

    zero_idx = int(np.argmin(np.abs(scores_raw))) if n_genes else 0
    ax_rank.axvline(zero_idx, color="#64748b", lw=0.9, ls=":")
    ax_rank.text(
        zero_idx, 0.88, f"Zero cross @ {zero_idx:,}",
        transform=ax_rank.get_xaxis_transform(),
        ha="center", va="top", fontsize=8.5, color=_DARK_MUTED,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.2},
    )
    ax_rank.set_xlim(0, max(1, n_genes - 1))
    ax_rank.set_xlabel("Rank in Ordered Dataset (high → low log₂FC)",
                       fontsize=10.5, color=_DARK_TEXT)
    ax_rank.set_ylabel("log₂FC", fontsize=10, color=_DARK_TEXT)
    ax_rank.grid(True, axis="y", color="#e2e8f0", lw=0.6)
    ax_rank.tick_params(colors=_DARK_MUTED, labelsize=9.5)
    for sp in ("top", "right"):
        ax_rank.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax_rank.spines[sp].set_edgecolor("#cbd5e1")

    fig.subplots_adjust(left=0.12, right=0.97, top=0.86, bottom=0.11)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight", facecolor="white")
    try:
        fig.savefig(
            output_png.with_suffix(".svg"), format="svg",
            bbox_inches="tight", facecolor="white", metadata={"Date": None},
        )
    except Exception:
        pass
    plt.close(fig)
    return str(output_png)


def _render_gsea_top_table(
    scores_raw,
    gene_to_rank: dict,
    library: dict,
    result_stable,
    contrast_name: str,
    library_label: str,
    n_sig_fdr25: int,
    output_png: Path,
    n: int = 15,
) -> str:
    """Render a full-width publication GSEA top-pathways table with hit barcodes
    for both positive (red) and negative (blue) NES (300-DPI PNG + SVG)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np

    top = result_stable.head(n)
    n_rows = len(top)
    n_genes = max(1, len(scores_raw))
    zero_frac = (
        float(np.argmin(np.abs(scores_raw))) / float(n_genes)
        if len(scores_raw) else 0.5
    )

    fig_h = max(3.6, 1.35 + 0.42 * n_rows)
    fig, ax = plt.subplots(figsize=(11.5, fig_h), facecolor="white")
    ax.set_facecolor("white")
    ax.axis("off")

    # Coordinate layout across x in [0.0, 1.0]:
    # [0.00..0.09] NES | [0.09..0.18] FDR | [0.19..0.46] Barcode | [0.48..1.00] Gene Set
    x_nes = 0.08
    x_fdr = 0.17
    bar_x0, bar_x1 = 0.195, 0.465
    bar_w = bar_x1 - bar_x0
    x_term = 0.485

    # Y coordinates: header at y = n_rows + 0.35, rows at y = n_rows-1 .. 0
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.75, n_rows + 1.05)

    # Column headers
    hdr_y = n_rows + 0.35
    ax.text(x_nes, hdr_y, "NES", ha="right", va="center",
            fontsize=10.5, fontweight="bold", color="#0f172a")
    ax.text(x_fdr, hdr_y, "FDR", ha="right", va="center",
            fontsize=10.5, fontweight="bold", color="#0f172a")
    ax.text((bar_x0 + bar_x1) / 2.0, hdr_y, "Ranked Gene Hits",
            ha="center", va="center",
            fontsize=10.5, fontweight="bold", color="#0f172a")
    ax.text(x_term, hdr_y, "Gene Set (size)", ha="left", va="center",
            fontsize=10.5, fontweight="bold", color="#0f172a")
    ax.hlines(n_rows - 0.02, 0.0, 1.0, color="#334155", lw=1.3)

    # Zero-crossing dashed guide down the barcode column
    zx = bar_x0 + zero_frac * bar_w
    ax.vlines(zx, ymin=-0.05, ymax=n_rows - 0.05,
              color="#94a3b8", lw=0.8, ls=":", zorder=2)

    for idx, (term, row) in enumerate(top.iterrows()):
        y_center = (n_rows - 1 - idx) + 0.45
        y_bot = (n_rows - 1 - idx)
        y_top = y_bot + 0.92

        # Subtle alternating zebra shading
        if idx % 2 == 1:
            ax.add_patch(mpatches.Rectangle(
                (0.0, y_bot), 1.0, 0.96,
                facecolor="#f8fafc", edgecolor="none", zorder=0,
            ))

        nes = float(row["nes"])
        fdr = float(row["fdr"])
        is_pos = nes >= 0
        color = "#dc2626" if is_pos else "#2563eb"

        # Barcode track background box
        ax.add_patch(mpatches.Rectangle(
            (bar_x0, y_bot + 0.08), bar_w, 0.76,
            facecolor="white", edgecolor="#cbd5e1", linewidth=0.6, zorder=1,
        ))

        gs = library.get(term, ())
        hits = [gene_to_rank[g] for g in gs if g in gene_to_rank]
        if hits:
            hx = bar_x0 + (np.asarray(hits, dtype=float) / float(n_genes)) * bar_w
            ax.vlines(hx, ymin=y_bot + 0.11, ymax=y_top - 0.11,
                      color=color, lw=0.75, alpha=0.70, zorder=3)

        gsize = _geneset_size(row, len(hits))
        term_str = str(term)
        if len(term_str) > 62:
            term_str = term_str[:59].rstrip() + "…"

        ax.text(x_nes, y_center, f"{nes:+.2f}", ha="right", va="center",
                fontsize=9.5, fontweight="bold", color=color, zorder=4)
        ax.text(x_fdr, y_center, _fmt_fdr(fdr), ha="right", va="center",
                fontsize=9.0, color=_DARK_MUTED, zorder=4)
        ax.text(x_term, y_center, f"{term_str}  (n={gsize})",
                ha="left", va="center", fontsize=9.5, color=_DARK_TEXT, zorder=4)

        ax.hlines(y_bot, 0.0, 1.0, color="#e2e8f0", lw=0.6, zorder=1)

    # Bottom axis annotations under the barcode column
    ax.hlines(0.0, 0.0, 1.0, color="#334155", lw=1.0, zorder=2)
    ax.text(bar_x0, -0.35, "↑ High log₂FC", ha="left", va="center",
            fontsize=8.5, color="#dc2626", fontweight="bold")
    ax.text(zx, -0.35, "0", ha="center", va="center",
            fontsize=8.5, color="#64748b")
    ax.text(bar_x1, -0.35, "Low log₂FC ↓", ha="right", va="center",
            fontsize=8.5, color="#2563eb", fontweight="bold")

    ax.set_title(
        f"GSEA Top {n_rows} Enriched Pathways — {contrast_name}\n"
        f"{library_label} · {n_sig_fdr25} significant pathways (FDR < 0.25) · "
        f"Red: NES > 0   Blue: NES < 0",
        fontsize=11.5, fontweight="bold", color="#0f172a", pad=12, loc="left",
    )

    plt.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight", facecolor="white")
    try:
        fig.savefig(
            output_png.with_suffix(".svg"), format="svg",
            bbox_inches="tight", facecolor="white", metadata={"Date": None},
        )
    except Exception:
        pass
    plt.close(fig)
    return str(output_png)


def make_gsea_running_sums(de_results_df,
                              symbol_map: dict,
                              contrast_name: str,
                              output_dir: str,
                              gene_set_library: str = "GO_Biological_Process_2023",
                              organism: str = "human",
                              top_n_plots: int = 3,
                              min_size: int = 15,
                              max_size: int = 500,
                              seed: int = 0) -> dict:
    """
    Run GSEA via blitzgsea on a ranked gene list (log2FC), then plot the
    top N enriched pathways as running sum plots and a summary top-table
    with hit barcodes.

    Pattern from user's gsea_results.py and gsea_graph.py.

    Args:
        de_results_df: DESeq2 results DataFrame (index = gene_id,
                       must have columns log2FoldChange and padj)
        symbol_map:    {ensembl_id (no version): hgnc_symbol}
        contrast_name: e.g. "treated vs control"
        output_dir:    figures dir for this contrast
        gene_set_library: blitzgsea / GMT library name
        organism:      "human" or "mouse" (matches gseapy convention)
        top_n_plots:   how many running sum plots to make
        min_size:      minimum gene-set size (Broad GSEA standard = 15)
        max_size:      maximum gene-set size (Broad GSEA standard = 500;
                       avoids mega-generic >500-gene GO terms whose Gamma
                       anchor extrapolation destabilizes blitzgsea NES)
        seed:          deterministic RNG seed for blitzgsea permutations

    Returns:
        dict with:
          'gsea_table':         path to GSEA results CSV (stable rows first)
          'top_table_fig':      path to summary table figure
          'running_sums':       list of paths to running sum plots
          'n_pathways':         number of significant finite pathways (FDR<0.25)
          'n_numeric_unstable': count of rows dropped from plots due to non-finite NES/FDR
    """
    out = {
        "gsea_table":         None,
        "top_table_fig":      None,
        "running_sums":       [],
        "n_pathways":         0,
        "n_numeric_unstable": 0,
    }

    try:
        import blitzgsea as blitz
        import pandas as pd

        # ── Build the ranked signature ───────────────────────────────────
        # blitzgsea expects DataFrame with columns [gene, score], where score
        # is log2FC (positive = upregulated). Upper-case symbols for mouse/human
        # library matching while preserving descending rank order.
        signature = _ranked_signature_frame(de_results_df, symbol_map)
        if signature.empty:
            log.warning(f"GSEA[{contrast_name}]: no genes after symbol mapping")
            return out

        sig_gsea, scores_raw, scores_centered, gene_to_rank = (
            _prepare_gsea_signature(signature)
        )
        if sig_gsea.empty:
            return out

        log.info(f"GSEA[{contrast_name}]: {len(sig_gsea)} ranked genes")

        # ── Load gene set library (local GMT first; respects air-gap) ────
        library, library_label = _load_gsea_library(gene_set_library, blitz)
        if not library:
            log.warning(
                f"GSEA[{contrast_name}]: gene set library '{gene_set_library}' unavailable"
            )
            return out

        # ── Run GSEA ─────────────────────────────────────────────────────
        try:
            result = blitz.gsea(
                sig_gsea, library,
                min_size=min_size,
                max_size=max_size,
                seed=seed,
            )
        except TypeError:
            # Compatibility fallback for minimal test doubles / older signatures
            result = blitz.gsea(sig_gsea, library)
        except Exception as e:
            log.warning(f"GSEA[{contrast_name}]: gsea() failed: {e}")
            return out

        if result is None or len(result) == 0:
            return out

        result_stable, result_unstable = _partition_gsea_results(result)
        n_unstable = int(len(result_unstable))
        out["n_numeric_unstable"] = n_unstable
        if n_unstable:
            log.warning(
                f"GSEA[{contrast_name}]: excluded {n_unstable} term(s) with "
                "non-finite NES/FDR from plots and significance counts."
            )

        # Save the full GSEA table (finite rows sorted first; any non-finite
        # rows preserved at the bottom for transparent narrator disclosure).
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        table_path = Path(output_dir) / "gsea_results.csv"
        result_export = (
            pd.concat([result_stable, result_unstable])
            if n_unstable
            else result_stable
        )
        result_export.to_csv(table_path)
        out["gsea_table"] = str(table_path)
        if result_stable.empty:
            return out

        n_sig_fdr25 = int((result_stable["fdr"] < 0.25).sum())
        out["n_pathways"] = n_sig_fdr25

        # ── Top-N running sum plots (strictly finite pathways) ───────────
        top = result_stable.head(top_n_plots)
        for i, (term, row) in enumerate(top.iterrows()):
            try:
                rank_path = Path(output_dir) / f"gsea_running_sum_{i+1}.png"
                saved = _render_gsea_running_sum(
                    scores_raw=scores_raw,
                    scores_centered=scores_centered,
                    gene_to_rank=gene_to_rank,
                    term=str(term),
                    gene_set=library.get(term, ()),
                    row=row,
                    contrast_name=contrast_name,
                    output_png=rank_path,
                )
                out["running_sums"].append(saved)
            except Exception as e:
                log.warning(f"GSEA running_sum {i+1} failed: {e}")
                continue

        # ── Top table figure (summary; strictly finite pathways) ─────────
        try:
            n_table = min(15, len(result_stable))
            tt_path = Path(output_dir) / "gsea_top_table.png"
            out["top_table_fig"] = _render_gsea_top_table(
                scores_raw=scores_raw,
                gene_to_rank=gene_to_rank,
                library=library,
                result_stable=result_stable,
                contrast_name=contrast_name,
                library_label=library_label,
                n_sig_fdr25=n_sig_fdr25,
                output_png=tt_path,
                n=n_table,
            )
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
        df["clean_id"] = df["gene_id"].map(_strip_gene_version)
        df["gene_symbol"] = df["clean_id"].map(symbol_map or {}).fillna("")

        cols = ["gene_id", "gene_symbol"]
        for c in ["baseMean", "log2FoldChange", "log2FoldChange_raw", "lfcSE",
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
