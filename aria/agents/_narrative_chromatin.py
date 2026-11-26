"""scATAC / chromatin figure generation (W0.1, scATAC P0 pre-preprint plan).

Mirrors ``aria.agents._narrative_scrna.generate_figures``: render publication
figures in the chromatin conda stack via ``EnvironmentManager`` and write their
paths into ``findings["figures"]``, which ``ChromatinNarrator.figures()`` surfaces.

No fabrication: a figure is produced ONLY when its underlying data exists (a
clustered ``.h5ad`` with an embedding) and an env manager is available. With
``env_manager=None`` or no ``.h5ad`` the figure set stays empty (honest absence),
exactly like the modality's QC ``None`` discipline.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from aria.agents.narrative.narrators.chromatin import _CHROMATIN_MODALITY_KEYS

log = logging.getLogger(__name__)


def live_findings(agent_result: dict) -> Optional[dict]:
    """Return the LIVE findings dict that ``unwrap_chromatin_findings`` surfaces.

    ``unwrap_chromatin_findings`` rebuilds a fresh dict when the per-modality
    wrapper is present, so a caller that wants a mutation (adding ``figures``) to be
    visible on the narrator's re-unwrap must mutate the SAME nested dict unwrap
    reads from. This returns that object.
    """
    if not isinstance(agent_result, dict):
        return None
    findings = agent_result.get("findings", agent_result)
    if not isinstance(findings, dict):
        return None
    for key, value in findings.items():
        if key in _CHROMATIN_MODALITY_KEYS and isinstance(value, dict):
            nested = value.get("findings")
            return nested if isinstance(nested, dict) else value
    return findings


def _umap_color_keys(findings: dict) -> list:
    """Obs columns to color the UMAP by, in priority order. ``rna_figure_umap``
    silently skips columns that are absent, so listing fallbacks is safe."""
    da = findings.get("differential_accessibility") or {}
    lsi = findings.get("lsi") or findings.get("lsi_clustering") or {}
    groupby = da.get("groupby") or lsi.get("cluster_key") or "leiden"
    keys: list = []
    # `log10_n_fragments` is the per-cell accessibility depth lsi writes to obs
    # (W0.2 leftover): a standard scATAC UMAP QC overlay. rna_figure_umap silently
    # skips obs columns that are absent, so listing it is safe on older outputs.
    for candidate in (groupby, "leiden", "sample", "sample_id", "batch",
                      "condition", "log10_n_fragments"):
        if candidate and candidate not in keys:
            keys.append(candidate)
    return keys


def _safe(value) -> str:
    return (str(value).replace("/", "_").replace(" ", "_")
            .replace(":", "_").replace("|", "_"))


def _save_dual(fig, png_path) -> str:
    """P4.1: save a line-art chromatin figure as raster PNG (inlined in the HTML
    report) AND a publication-grade vector SVG alongside it (`<stem>.svg`), matching
    the RNA side's vector convention. Only the PNG path is returned/surfaced as the
    report figure; the SVG sits next to it for manuscript use, so the narrator's
    one-figure-per-key contract is unchanged. UMAP stays raster (rna_figure_umap)."""
    png_path = Path(png_path)
    fig.savefig(png_path, dpi=160, bbox_inches="tight", facecolor="white")
    fig.savefig(png_path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    return str(png_path)


def _neglog10_capped(p, headroom: float = 1.12, all_zero_floor: float = 10.0):
    """Wave 2A: data-adaptive ``-log10(p)`` for display via canonical ``stats.adaptive_neglog10``."""
    import numpy as np
    from aria.utils.stats import adaptive_neglog10

    arr = np.asarray(p, dtype=float)
    y, cap, has_zero = adaptive_neglog10(
        arr, headroom=headroom, all_zero_floor=all_zero_floor
    )
    zero = np.isfinite(arr) & (arr == 0.0)
    return y, zero, (cap if has_zero else None)


def _finite_unit_mask(p):
    """Boolean mask of finite p-values inside ``[0, 1]``."""
    import numpy as np

    p = np.asarray(p, dtype=float)
    return np.isfinite(p) & (p >= 0.0) & (p <= 1.0)


_BINDING_MODALITIES = frozenset({"ChIP", "CUT_AND_RUN", "CUT_AND_TAG"})


def _signal_terms(data_type: Optional[str]) -> tuple[str, str, str]:
    """Return ``(noun, title_phrase, peak_phrase)`` for binding vs accessibility."""
    if str(data_type or "") in _BINDING_MODALITIES:
        return "binding", "Differential binding", "differential-binding"
    return "accessibility", "Differential accessibility", "differential-accessibility"


def render_da_figures(da_full_csv: str, output_dir: Path,
                      padj_max: float = 0.05,
                      data_type: Optional[str] = None) -> dict:
    """W0.3: pseudobulk-DA volcano + MA per comparison, rendered INLINE (matplotlib
    over the full per-peak CSV from `chromatin_diffacc`). T15 honesty: only peaks in
    the convergence-gated `significant` column are drawn as up/down; everything else
    (incl. non-converged) is background grey, never a real effect. Returns
    {fig_key: path}; empty when the table is missing/empty."""
    out: dict = {}
    if not da_full_csv or not Path(da_full_csv).exists():
        return out
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    noun, title_phrase, _peak_phrase = _signal_terms(data_type)
    is_binding = str(data_type or "") in _BINDING_MODALITIES

    df = pd.read_csv(da_full_csv)
    # Accept both column conventions: the scATAC pseudobulk CSV
    # (`chromatin_diffacc`) writes `log2fc`/`base_mean`, while the bulk ATAC CSV
    # (`chromatin_bulk_diffacc`) writes the native DESeq2 `log2FoldChange`/`baseMean`.
    # Renaming is a no-op when a convention's columns are absent, so neither lane
    # is disturbed.
    df = df.rename(columns={"log2FoldChange": "log2fc", "baseMean": "base_mean"})
    if df.empty or "comparison" not in df.columns:
        return out
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for comp, sub in df.groupby("comparison"):
        comp_s = _safe(comp)
        lfc_all = pd.to_numeric(sub["log2fc"], errors="coerce").to_numpy(dtype=float)
        padj_all = pd.to_numeric(sub["padj"], errors="coerce").to_numpy(dtype=float)
        sig_all = (sub["significant"].fillna(False).to_numpy(dtype=bool)
                   if "significant" in sub.columns
                   else np.zeros(len(sub), dtype=bool))
        # Wave 2A: only finite log2FC and finite padj in [0, 1] are drawable;
        # NaN/inf/negative padj never reach the axes.
        keep = np.isfinite(lfc_all) & _finite_unit_mask(padj_all)
        lfc, padj, sig = lfc_all[keep], padj_all[keep], sig_all[keep]

        # ── Volcano ───────────────────────────────────────────────────────
        fig = None
        if lfc.size:
            try:
                # Data-adaptive display: padj == 0 underflow sits just above the
                # finite cloud (never at y=300), and is disclosed explicitly.
                neglog, zero, cap = _neglog10_capped(padj)
                fig, ax = plt.subplots(figsize=(6.5, 6), dpi=160)
                ns = ~sig & ~zero
                ax.scatter(lfc[ns], neglog[ns], s=4, c="#bdbdbd",
                           linewidths=0, alpha=0.5, label="not significant")
                up, down = sig & (lfc > 0), sig & (lfc < 0)
                ax.scatter(lfc[up & ~zero], neglog[up & ~zero], s=7, c="#c0392b",
                           linewidths=0, alpha=0.85, label="up")
                ax.scatter(lfc[down & ~zero], neglog[down & ~zero], s=7,
                           c="#2c7fb8", linewidths=0, alpha=0.85, label="down")
                if cap is not None:
                    zc = np.where(sig & (lfc > 0), "#c0392b",
                                  np.where(sig & (lfc < 0), "#2c7fb8", "#7f7f7f"))
                    ax.scatter(lfc[zero], neglog[zero], s=16, c=list(zc[zero]),
                               marker="^", linewidths=0.3, edgecolors="black",
                               alpha=0.9,
                               label=f"padj = 0 (underflow, n={int(zero.sum())}) "
                                     f"drawn at cap {cap:.1f}")
                    ax.axhline(cap, color="#7f7f7f", lw=0.4, ls=":")
                ax.axhline(-np.log10(padj_max), color="black", lw=0.5, ls="--")
                ax.axvline(0, color="black", lw=0.4)
                ax.set_xlabel(f"log2 fold change ({noun})", fontsize=9)
                ax.set_ylabel("-log10 padj" + (" (padj=0 capped)"
                                               if cap is not None else ""),
                              fontsize=9)
                ax.set_title(f"{title_phrase} — {comp}", fontsize=10,
                             fontweight="bold")
                ax.legend(loc="upper right", fontsize=7, frameon=False,
                          markerscale=2.0)
                ax.tick_params(labelsize=7)
                path = output_dir / f"da_volcano_{comp_s}.png"
                _save_dual(fig, path)
                plt.close(fig)
                out[f"da_volcano_{comp_s}"] = str(path)
            except Exception as e:
                if fig is not None:
                    plt.close(fig)
                log.warning("DA volcano failed for %s: %s", comp, e)

        # ── MA (needs baseMean) ───────────────────────────────────────────
        fig = None
        if "base_mean" in sub.columns and not sub["base_mean"].isna().all():
            try:
                base_all = pd.to_numeric(sub["base_mean"],
                                         errors="coerce").to_numpy(dtype=float)
                ma_keep = np.isfinite(lfc_all) & np.isfinite(base_all)
                if not ma_keep.any():
                    continue
                base = base_all[ma_keep]
                lfc_ma, sig_ma = lfc_all[ma_keep], sig_all[ma_keep]
                x = np.log10(np.clip(base, 1e-3, None))
                fig, ax = plt.subplots(figsize=(6.5, 5), dpi=160)
                ax.scatter(x[~sig_ma], lfc_ma[~sig_ma], s=4, c="#bdbdbd",
                           linewidths=0, alpha=0.5)
                ax.scatter(x[sig_ma], lfc_ma[sig_ma], s=7, c="#c0392b",
                           linewidths=0, alpha=0.85, label="significant")
                ax.axhline(0, color="black", lw=0.5)
                ax.set_xlabel(f"log10 mean {noun} (baseMean)", fontsize=9)
                ax.set_ylabel("log2 fold change", fontsize=9)
                ax.set_title(
                    f"Differential binding MA plot — {comp}"
                    if is_binding
                    else f"Differential accessibility MA plot — {comp}",
                    fontsize=10,
                    fontweight="bold",
                )
                ax.legend(loc="upper right", fontsize=7, frameon=False,
                          markerscale=2.0)
                ax.tick_params(labelsize=7)
                path = output_dir / f"da_ma_{comp_s}.png"
                _save_dual(fig, path)
                plt.close(fig)
                out[f"da_ma_{comp_s}"] = str(path)
            except Exception as e:
                if fig is not None:
                    plt.close(fig)
                log.warning("DA MA failed for %s: %s", comp, e)
    return out


def render_sample_correlation(counts_tsv, png_path) -> Optional[str]:
    """Bulk ATAC sample-QC: Pearson correlation heatmap across replicate samples,
    computed on log1p peak counts (the standard ATAC sample-similarity diagnostic).
    Reads the replicate peak-count TSV (`peak_id` + one column per replicate sample)
    emitted by `chromatin_bulk_diffacc`. Honest-skip (returns None) when the matrix
    is missing or has fewer than two samples. Writes dual PNG+SVG."""
    if not counts_tsv or not Path(counts_tsv).exists():
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    df = pd.read_csv(counts_tsv, sep="\t")
    peak_col = df.columns[0]
    samples = [c for c in df.columns if c != peak_col]
    if len(samples) < 2:
        return None
    mat = np.log1p(df[samples].to_numpy(dtype=float))
    corr = np.corrcoef(mat, rowvar=False)
    if not np.all(np.isfinite(corr)):
        return None

    n = len(samples)
    fig, ax = plt.subplots(figsize=(max(4.5, 0.6 * n + 2.5),
                                    max(4.0, 0.6 * n + 2.0)), dpi=160)
    vmin = min(0.0, float(corr.min()))
    im = ax.imshow(corr, cmap="viridis", vmin=vmin, vmax=1.0)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(samples, fontsize=7)
    # Wave 2A: text contrast follows the COLORMAP midpoint (viridis is dark below
    # it, bright yellow above), not the data midpoint — otherwise replicates in
    # 0.90-1.00 got white text on bright-yellow cells.
    thresh = 0.5 * (vmin + 1.0)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center",
                    fontsize=6,
                    color="white" if corr[i, j] < thresh else "black")
    ax.set_title("Replicate sample correlation (log1p peak counts)",
                 fontsize=10, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    out = _save_dual(fig, Path(png_path))
    plt.close(fig)
    return out


def render_peak_annotation(annotation: dict, png_path,
                           data_type: Optional[str] = None) -> Optional[str]:
    """Bulk ATAC / epigenomic genomic-context figure: horizontal bar chart of the
    differential-peak feature distribution (Promoter/Exonic/Intronic/Distal
    Intergenic) from `chromatin_peak_annotation`. Dual PNG+SVG. Honest-None when
    no annotation ran or the distribution is empty."""
    dist = (annotation or {}).get("feature_distribution_overall") or {}
    dist = {k: int(v) for k, v in dist.items() if v}
    if not dist:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    eff_dt = data_type or (annotation or {}).get("data_type")
    _noun, _title_phrase, peak_phrase = _signal_terms(eff_dt)
    is_binding = str(eff_dt or "") in _BINDING_MODALITIES

    # Canonical ChIPseeker-style order (most proximal first).
    order = ["Promoter", "Exonic", "Intronic", "Distal Intergenic"]
    feats = [f for f in order if f in dist] + [f for f in dist if f not in order]
    counts = [dist[f] for f in feats]
    total = sum(counts) or 1
    colors = {"Promoter": "#c0392b", "Exonic": "#e67e22",
              "Intronic": "#2c7fb8", "Distal Intergenic": "#7f8c8d"}

    fig, ax = plt.subplots(figsize=(7, max(2.4, 0.6 * len(feats) + 1.4)), dpi=160)
    ypos = range(len(feats))
    ax.barh(list(ypos), counts,
            color=[colors.get(f, "#95a5a6") for f in feats])
    for y, c in zip(ypos, counts):
        ax.text(c, y, f" {c:,} ({100 * c / total:.0f}%)", va="center",
                fontsize=8)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(feats, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Differential binding peaks" if is_binding else "DA peaks",
                  fontsize=9)
    ax.set_title(f"Genomic distribution of {peak_phrase} peaks",
                 fontsize=10, fontweight="bold")
    ax.margins(x=0.18)
    out = _save_dual(fig, Path(png_path))
    plt.close(fig)
    return out


def render_motif_dotplot(motifs: dict, out_path: Path, top_n: int = 8) -> Optional[str]:
    """W0.3: motif-enrichment dotplot (group × motif; size = -log10 padj, color =
    log2 enrichment), rendered INLINE from the structured `motifs["per_group"]`
    (top_motifs carry name/log2_enrichment/padj). Returns the path or None when
    there is nothing enriched to plot. A motif id is a DB identifier, not a claim."""
    import math

    per_group = (motifs or {}).get("per_group") or {}
    raw_padj_max = (motifs or {}).get("padj_max", 0.05)
    try:
        padj_max = float(raw_padj_max) if raw_padj_max is not None else 0.05
    except (TypeError, ValueError):
        padj_max = 0.05
    if not (math.isfinite(padj_max) and 0.0 < padj_max <= 1.0):
        padj_max = 0.05
    # Wave 2A + Wave 5C: keep only genuinely SIGNIFICANT ENRICHED motifs
    # (finite log2_enrichment > 0 and finite 0 <= padj < padj_max).
    best: dict = {}
    for group, info in per_group.items():
        for m in (info or {}).get("top_motifs", []) or []:
            if not isinstance(m, dict):
                continue
            try:
                padj = float(m.get("padj"))
                l2e = float(m.get("log2_enrichment"))
            except (TypeError, ValueError):
                continue
            if not (math.isfinite(padj) and 0.0 <= padj < padj_max):
                continue
            if not (math.isfinite(l2e) and l2e > 0.0):
                continue
            row = {"group": str(group),
                   "motif": str(m.get("name") or m.get("motif_id") or "?"),
                   "log2_enrichment": l2e, "padj": padj}
            key = (row["group"], row["motif"])
            # Duplicate TF names (several JASPAR matrices) overplot at identical
            # (x, y): keep the strongest (lowest padj, then largest enrichment).
            prev = best.get(key)
            if prev is None or (padj, -l2e) < (prev["padj"], -prev["log2_enrichment"]):
                best[key] = row
    by_group: dict = {}
    for row in best.values():
        by_group.setdefault(row["group"], []).append(row)
    rows = []
    for group in by_group:
        ranked = sorted(by_group[group],
                        key=lambda r: (r["padj"], -r["log2_enrichment"], r["motif"]))
        rows.extend(ranked[:top_n])
    if not rows:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from matplotlib.colors import Normalize

    df = pd.DataFrame(rows)
    groups = sorted(df["group"].unique())
    motif_order = list(dict.fromkeys(df["motif"]))
    gx = {g: i for i, g in enumerate(groups)}
    my = {m: i for i, m in enumerate(motif_order)}
    fig = None
    try:
        fig, ax = plt.subplots(
            figsize=(max(4.0, 1.1 * len(groups) + 2),
                     max(3.0, 0.32 * len(motif_order) + 1)), dpi=160)
        neglog, zero, cap = _neglog10_capped(df["padj"].to_numpy(float),
                                             headroom=1.15)
        # Bounded marker AREA: -log10(padj) mapped linearly onto [25, 280] pt²
        # (the old 12 + 40*300 = 12,012 pt² swallowed neighbouring rows).
        s_min, s_max = 25.0, 280.0
        y_top = float(np.max(neglog)) if neglog.size else 0.0

        def _area(v):
            v = np.asarray(v, dtype=float)
            if y_top <= 0:
                return np.full(v.shape, s_min)
            return s_min + (s_max - s_min) * np.clip(v / y_top, 0.0, 1.0)

        sizes = _area(neglog)
        colors = df["log2_enrichment"].to_numpy(dtype=float)
        vmin, vmax = float(np.min(colors)), float(np.max(colors))
        if not vmax > vmin:   # degenerate Normalize(v, v)
            vmin, vmax = vmin - 0.5, vmax + 0.5
        sc = ax.scatter([gx[g] for g in df["group"]],
                        [my[m] for m in df["motif"]],
                        s=sizes, c=colors, cmap="viridis",
                        norm=Normalize(vmin=vmin, vmax=vmax),
                        edgecolors="black", linewidths=0.3)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(motif_order)))
        ax.set_yticklabels(motif_order, fontsize=7)
        ax.set_xlim(-0.6, len(groups) - 0.4)
        ax.set_ylim(len(motif_order) - 0.5, -0.5)
        ax.set_title("Motif enrichment (size = -log10 padj; enriched only)",
                     fontsize=10, fontweight="bold")
        plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02, label="log2 enrichment")
        # Size legend: three representative -log10(padj) levels.
        if y_top > 0:
            raw_levels = [v for v in (y_top / 4.0, y_top / 2.0, y_top) if v > 0]
            if y_top >= 1.0:
                levels = sorted({round(v, 1) for v in raw_levels if round(v, 1) > 0})
            else:
                levels = sorted({float(f"{v:.2g}") for v in raw_levels})
            if not levels:
                levels = [float(f"{y_top:.2g}")]
        else:
            levels = [0.0]
        handles = [ax.scatter([], [], s=float(_area([lv])[0]), c="#bdbdbd",
                              edgecolors="black", linewidths=0.3)
                   for lv in levels]
        fmt_level = (lambda lv: f"{lv:g}" if lv >= 1.0 else f"{lv:.2g}")
        labels = [fmt_level(lv) for lv in levels]
        title = "-log10 padj"
        if cap is not None:
            labels[-1] = f"{fmt_level(levels[-1])} (padj=0 cap)"
            title += f"\n({int(zero.sum())} padj=0 capped at {cap:.1f})"
        ax.legend(handles, labels, title=title, loc="upper left",
                  bbox_to_anchor=(1.28, 1.0), frameon=False, fontsize=6,
                  title_fontsize=6, labelspacing=1.2, borderpad=0.8)
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path = _save_dual(fig, out_path)
        plt.close(fig)
        return out_path
    except Exception as e:
        if fig is not None:
            plt.close(fig)
        log.warning("motif dotplot failed: %s", e)
        return None


# Fragment-size class bands, mirroring the windows `chromatin_qc` measures
# (ATAC: NFR <100, mono 180-247, di 315-473 bp; CUT&RUN/CUT&Tag: sub-nucleosomal
# <120, mono 140-250 bp). ChIP has no nucleosomal banding expectation -> no bands.
_ATAC_FRAGMENT_BANDS = (
    (0, 100, "NFR (<100 bp)", "#2ecc71"),
    (180, 247, "mono-nucleosome (180-247 bp)", "#e67e22"),
    (315, 473, "di-nucleosome (315-473 bp)", "#8e44ad"),
)
_CUTRUN_FRAGMENT_BANDS = (
    (0, 120, "sub-nucleosomal (<120 bp)", "#2ecc71"),
    (140, 250, "mono-nucleosome (140-250 bp)", "#e67e22"),
)
_MODALITY_LABEL = {
    "scATAC": "scATAC", "bulk_ATAC": "bulk ATAC", "ChIP": "ChIP-seq",
    "CUT_AND_RUN": "CUT&RUN", "CUT_AND_TAG": "CUT&Tag",
}
# Read-only mirror of `chromatin_qc._default_bulk_qc_thresholds` / the scATAC
# FRiP gate, used only when the QC payload carries no explicit `thresholds`.
_DEFAULT_FRIP_GATE = {
    "scATAC": 0.20, "bulk_ATAC": 0.20, "ChIP": 0.01,
    "CUT_AND_RUN": 0.05, "CUT_AND_TAG": 0.10,
}


def _fragment_bands(data_type: Optional[str]):
    if data_type in ("CUT_AND_RUN", "CUT_AND_TAG"):
        return _CUTRUN_FRAGMENT_BANDS
    if data_type == "ChIP":
        return ()
    # scATAC / bulk_ATAC / legacy payloads without a data_type (scATAC default).
    return _ATAC_FRAGMENT_BANDS


def render_fragment_size_figure(qc: dict, out_path: Path,
                                data_type: Optional[str] = None) -> Optional[str]:
    """W0.2: fragment-size / nucleosome-banding plot from the QC size histogram
    (`qc["fragment_sizes"]["size_histogram"]`). Rendered INLINE. Returns None when
    no histogram was computed (e.g. a `.h5mu` peak-matrix run with no fragments) —
    honest absence, never a fabricated distribution.

    Wave 2A: drawn as a true step histogram, with modality-appropriate size-class
    bands (ATAC vs CUT&RUN/CUT&Tag; none for ChIP) matching `chromatin_qc`, and the
    subsample (size-sample cap / truncated scan / per-BAM read cap) disclosed."""
    fs = (qc or {}).get("fragment_sizes") or {}
    hist = fs.get("size_histogram") or {}
    edges = hist.get("bin_edges")
    counts = hist.get("counts")
    if edges is None or counts is None or len(edges) < 2 or len(counts) == 0:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    edges = np.asarray(edges, dtype=float)
    counts = np.asarray(counts, dtype=float)
    n = min(len(edges) - 1, len(counts))
    counts = np.where(np.isfinite(counts[:n]), counts[:n], 0.0)
    if n < 1 or float(counts.sum()) <= 0:
        return None
    data_type = data_type or (qc or {}).get("data_type")
    n_plot = int(round(float(counts.sum())))
    n_total = fs.get("n_fragments")
    truncated = bool(fs.get("scan_truncated"))
    if isinstance(n_total, (int, float)) and n_total > n_plot:
        sample_note = (f"subsample: {n_plot:,} of {int(n_total):,}"
                       f"{'+' if truncated else ''} fragments")
    elif truncated:
        sample_note = f"subsample: {n_plot:,} fragments (scan capped)"
    elif fs.get("n_fragments_sampled") is not None:
        sample_note = f"sampled: {n_plot:,} read pairs (per-BAM read cap)"
    else:
        sample_note = f"n={n_plot:,} fragments"
    fig = None
    try:
        fig, ax = plt.subplots(figsize=(6.5, 4), dpi=160)
        # Step histogram over the real bin edges (no piecewise-linear overlay that
        # disagrees with the step fill).
        step_y = np.append(counts, counts[-1])
        ax.fill_between(edges[:n + 1], step_y, step="post", alpha=0.6,
                        color="#34699a", linewidth=0)
        ax.step(edges[:n + 1], step_y, where="post", lw=0.8, color="#1f4068")
        for lo, hi, label, color in _fragment_bands(data_type):
            ax.axvspan(lo, hi, color=color, alpha=0.10, lw=0, label=label)
        ax.set_xlim(float(edges[0]), float(edges[n]))
        ax.set_xlabel("fragment size (bp)", fontsize=9)
        ax.set_ylabel("fragment count (subsample)" if "subsample" in sample_note
                      or "sampled" in sample_note else "fragment count",
                      fontsize=9)
        mod = _MODALITY_LABEL.get(str(data_type), "")
        ax.set_title(f"{mod + ' ' if mod else ''}fragment-size distribution\n"
                     f"({sample_note})", fontsize=10, fontweight="bold")
        ax.tick_params(labelsize=7)
        if _fragment_bands(data_type):
            ax.legend(fontsize=6.5, frameon=False, loc="upper right")
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path = _save_dual(fig, out_path)
        plt.close(fig)
        return out_path
    except Exception as e:
        if fig is not None:
            plt.close(fig)
        log.warning("fragment-size figure failed: %s", e)
        return None


def render_tss_depth_scatter(qc: dict, out_path: Path) -> Optional[str]:
    """P4.1 (2)[A]: the canonical ATAC QC gating panel — per-cell TSS enrichment
    vs log10 fragment depth, with the ENCODE TSSe gate lines. Reads the additive
    ``qc["tss_depth_scatter"]`` arrays (only present on a real fragments run, W0.5
    path). Returns None when TSSe arrays were not computed — honest absence, never
    a fabricated cloud. When the importer recorded no per-cell depth, falls back to
    a 1-D TSSe strip so the gate is still visible."""
    sc = (qc or {}).get("tss_depth_scatter") or {}
    tsse = sc.get("tsse")
    if tsse is None or len(tsse) == 0:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    tsse = np.asarray(tsse, dtype=float)
    depth = sc.get("log10_depth")
    min_tss = sc.get("min_tss")
    warn_tss = sc.get("warn_tss")
    # Wave 2A: `if depth and ...` raised on a NumPy array (silently dropping the
    # figure); test identity + length explicitly and mask non-finite points.
    has_depth = depth is not None and len(depth) == len(tsse)
    if has_depth:
        x = np.asarray(depth, dtype=float)
        keep = np.isfinite(tsse) & np.isfinite(x)
        x, tsse = x[keep], tsse[keep]
        xlabel = "log10 fragment depth"
    else:
        tsse = tsse[np.isfinite(tsse)]
        # No per-cell depth — jitter on x so points are separable.
        rng = np.random.default_rng(0)
        x = rng.normal(0.0, 0.04, size=tsse.shape[0])
        xlabel = "(no per-cell depth recorded)"
    if tsse.size == 0:
        return None
    fig = None
    try:
        fig, ax = plt.subplots(figsize=(6.0, 4.5), dpi=160)
        ax.scatter(x, tsse, s=4, alpha=0.35, color="#34699a", linewidths=0)
        for gate, c, lab in ((min_tss, "#c0392b", "ENCODE min"),
                             (warn_tss, "#e08e0b", "marginal")):
            if gate is not None:
                ax.axhline(float(gate), color=c, lw=0.8, ls="--",
                           label=f"{lab} ({float(gate):.0f})")
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel("TSS enrichment (per cell)", fontsize=9)
        ax.set_title("Per-cell QC: TSS enrichment vs depth",
                     fontsize=10, fontweight="bold")
        ax.tick_params(labelsize=7)
        if min_tss is not None or warn_tss is not None:
            ax.legend(fontsize=7, frameon=False)
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path = _save_dual(fig, out_path)
        plt.close(fig)
        return out_path
    except Exception as e:
        if fig is not None:
            plt.close(fig)
        log.warning("TSS×depth scatter figure failed: %s", e)
        return None


def render_frip_distribution(qc: dict, out_path: Path,
                             data_type: Optional[str] = None,
                             thresholds: Optional[dict] = None) -> Optional[str]:
    """P4.1 (2)[B]: FRiP distribution histogram with the FRiP quality gate. Reads
    the additive ``qc["frip_distribution"]`` list (per-barcode for scATAC; per-sample
    for bulk ATAC/ChIP/CUT&RUN/CUT&Tag). Returns None when no distribution was
    computed — honest absence, never a fabricated histogram.

    Wave 2A: the unit ("barcodes" vs "samples") and the gate follow the modality —
    ``thresholds["min_frip"]`` (``qc["thresholds"]``) when supplied, else the
    modality default (ChIP 0.01, CUT&RUN 0.05, CUT&Tag 0.10, ATAC 0.20)."""
    dist = (qc or {}).get("frip_distribution")
    if dist is None or len(dist) == 0:
        per_sample = (qc or {}).get("frip_per_sample")
        if isinstance(per_sample, dict) and per_sample:
            dist = [v for v in per_sample.values() if v is not None]
    if dist is None or len(dist) == 0:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    vals = np.asarray(dist, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return None
    data_type = data_type or (qc or {}).get("data_type")
    th = thresholds if isinstance(thresholds, dict) else (
        (qc or {}).get("thresholds") if isinstance((qc or {}).get("thresholds"),
                                                     dict) else {})
    gate = th.get("min_frip")
    try:
        gate = float(gate) if gate is not None else None
    except (TypeError, ValueError):
        gate = None
    if gate is None or not np.isfinite(gate):
        gate = _DEFAULT_FRIP_GATE.get(str(data_type), 0.20)
    is_bulk = data_type in _BULK_CHROMATIN_MODALITIES
    unit = "samples" if is_bulk else "barcodes"
    per = "Per-sample" if is_bulk else "Per-barcode"
    mod = _MODALITY_LABEL.get(str(data_type), "")
    fig = None
    try:
        fig, ax = plt.subplots(figsize=(6.0, 4.0), dpi=160)
        ax.hist(vals, bins=40, range=(0, 1), color="#34699a", alpha=0.75)
        ax.axvline(gate, color="#c0392b", lw=0.8, ls="--",
                   label=f"{mod + ' ' if mod else ''}FRiP gate (≥{gate:.2f})")
        ax.set_xlabel("FRiP (fraction of reads in peaks)", fontsize=9)
        ax.set_ylabel(unit, fontsize=9)
        ax.set_title(f"{per} FRiP distribution "
                     f"({mod + ', ' if mod else ''}n={vals.size} {unit})",
                     fontsize=10, fontweight="bold")
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7, frameon=False)
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path = _save_dual(fig, out_path)
        plt.close(fig)
        return out_path
    except Exception as e:
        if fig is not None:
            plt.close(fig)
        log.warning("FRiP distribution figure failed: %s", e)
        return None


def _marker_peaks_by_cluster(findings: dict) -> dict:
    """Extract {cluster_id: [peak, ...]} from the per-cluster DA findings."""
    da = findings.get("differential_accessibility") or {}
    pc = da.get("per_cluster") or {}
    by_cluster = pc.get("da_peaks_by_cluster") or {}
    out: dict = {}
    for cl, recs in by_cluster.items():
        peaks = [str(r.get("peak")) for r in (recs or []) if r.get("peak")]
        if peaks:
            out[str(cl)] = peaks
    return out


_BULK_ATAC_SPINE_KEYS = ("preprocessing", "qc", "peaks", "peak_counts",
                         "differential_accessibility")
_BULK_CHROMATIN_MODALITIES = frozenset({
    "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG",
})


def _is_bulk_atac_findings(findings: dict) -> bool:
    return any(isinstance(findings.get(k), dict)
               and findings[k].get("data_type") in _BULK_CHROMATIN_MODALITIES
               for k in _BULK_ATAC_SPINE_KEYS)


def _findings_data_type(findings: dict) -> Optional[str]:
    """Resolve the chromatin modality from the measured spine (bulk keys first,
    then any spine/QC `data_type`, e.g. "scATAC")."""
    seen = None
    for k in _BULK_ATAC_SPINE_KEYS:
        block = findings.get(k)
        if isinstance(block, dict) and block.get("data_type"):
            dt = str(block["data_type"])
            if dt in _BULK_CHROMATIN_MODALITIES:
                return dt
            seen = seen or dt
    return seen


def _bulk_figure_params(findings: dict, output_dir: Path) -> dict:
    """IPC params for `chromatin_figure_bulk` from the measured findings only.
    Unset optional fields are omitted (the script contract rejects empties)."""
    counts = findings.get("peak_counts") or {}
    da = findings.get("differential_accessibility") or {}
    counts_path = (counts.get("counts_cpm_path")
                   or da.get("replicate_counts_path")
                   or counts.get("counts_matrix_path"))
    params: dict = {"output_dir": str(output_dir)}
    if counts_path:
        params["counts_path"] = str(counts_path)
    # Replicate-aggregated columns are named "<condition>::<replicate>", so the
    # per-library metadata does not apply to them; only pass metadata for the
    # per-library matrices.
    if counts_path and counts_path != da.get("replicate_counts_path"):
        if counts.get("sample_metadata_path"):
            params["sample_metadata"] = str(counts["sample_metadata_path"])
    elif counts_path:
        cond_col = da.get("condition_col") or "condition"
        try:
            import csv
            with open(counts_path, newline="", encoding="utf-8") as handle:
                header = next(csv.reader(handle, delimiter="\t"))
            meta = {c: {cond_col: c.split("::", 1)[0]}
                    for c in header[1:] if "::" in c}
            if meta:
                params["sample_metadata"] = meta
        except (OSError, StopIteration):
            pass
    if da.get("condition_col"):
        params["condition_col"] = str(da["condition_col"])
    qc = findings.get("qc")
    if isinstance(qc, dict):
        params["qc"] = {k: v for k, v in qc.items() if k in (
            "frip", "frip_per_sample", "tss_enrichment",
            "tss_enrichment_per_sample", "dup_rate", "dup_rate_per_sample",
            "mito_frac", "mito_fraction")}
        # Wave 2A: modality-specific gates (ChIP FRiP 0.01, CUT&RUN 0.05, ...).
        if isinstance(qc.get("thresholds"), dict) and qc["thresholds"]:
            params["thresholds"] = dict(qc["thresholds"])
    data_type = _findings_data_type(findings)
    if data_type in _BULK_CHROMATIN_MODALITIES:
        params["data_type"] = data_type
    peaks = findings.get("peaks")
    if isinstance(peaks, dict):
        params["peaks"] = {k: v for k, v in peaks.items() if k in (
            "frip", "frip_per_sample", "n_peaks", "n_consensus_peaks")}
    return params


def _render_bulk_atac_figures(findings: dict, output_dir: Path,
                              env_manager=None) -> dict:
    params = _bulk_figure_params(findings, output_dir)
    if env_manager is not None:
        res = env_manager.run_in_stack(
            stack="rna",
            script_path="aria/scripts/chromatin_figure_bulk.py",
            params=params,
        )
    else:
        from aria.scripts.chromatin_figure_bulk import render_bulk_figures
        res = render_bulk_figures(params)
    if not isinstance(res, dict) or res.get("status") != "success":
        log.warning("bulk ATAC figures failed: %s — %s",
                    (res or {}).get("error_type"),
                    str((res or {}).get("details", ""))[:200])
        return {}
    return {k: v for k, v in (res.get("figures") or {}).items()
            if isinstance(v, str)}


def generate_figures(findings: dict, h5ad_path: Optional[str], output_dir,
                     env_manager=None) -> dict:
    """Generate scATAC figures into ``findings['figures']`` (mutates + returns it).

    - **UMAP** (W0.1): rendered via ``rna_figure_umap.py`` in the ``chromatin``
      stack (needs scanpy on the clustered ``.h5ad``).
    - **DA volcano + MA** and **motif dotplot** (W0.3): rendered INLINE (matplotlib
      over the DA full CSV / structured motif data), like the scRNA inline figures.

    No fabrication: a figure appears only when its data exists.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figs = findings.setdefault("figures", {})
    figs_svg = findings.setdefault("figures_svg", {})

    # 1. UMAP — needs the embedding in the .h5ad, so render in the chromatin stack.
    if h5ad_path and env_manager is not None:
        color_keys = _umap_color_keys(findings)
        try:
            res = env_manager.run_in_stack(
                stack="chromatin",
                script_path="aria/scripts/rna_figure_umap.py",
                params={
                    "h5ad_path":  str(h5ad_path),
                    "color_by":   color_keys,
                    "output_dir": str(output_dir),
                },
            )
            if res.get("status") == "success":
                for key, path in (res.get("figures") or {}).items():
                    figs[f"umap_{key}"] = path
                for key, path in (res.get("figures_svg") or {}).items():
                    if isinstance(path, str):
                        figs_svg[f"umap_{key}"] = path
            else:
                log.warning(
                    "chromatin UMAP figure failed: %s — %s",
                    res.get("error_type"), str(res.get("details", ""))[:200],
                )
        except Exception as e:  # subprocess crash must not abort the report
            log.warning("chromatin UMAP figure subprocess crashed: %s", e)

    # 2. DA volcano + MA (inline) from the convergence-gated full DA table.
    da = findings.get("differential_accessibility") or {}
    padj_max = float(da.get("padj_max", 0.05))
    da_dt = da.get("data_type") or _findings_data_type(findings)
    # scATAC shape: a single pseudobulk table.
    pb = da.get("pseudobulk") or {}
    da_full_csv = pb.get("full_results_csv")
    if da_full_csv:
        try:
            figs.update(render_da_figures(
                da_full_csv, output_dir, padj_max, data_type=da_dt))
        except Exception as e:
            log.warning("chromatin DA figures failed: %s", e)
    # Bulk chromatin shape (bulk_ATAC, ChIP, CUT&RUN, CUT&Tag): one
    # convergence-gated full CSV per condition comparison.
    if da.get("data_type") in _BULK_CHROMATIN_MODALITIES:
        for comp in (da.get("comparisons") or []):
            comp_csv = (comp or {}).get("full_results_csv")
            if not comp_csv:
                continue
            try:
                figs.update(render_da_figures(
                    comp_csv, output_dir, padj_max,
                    data_type=da.get("data_type")))
            except Exception as e:
                log.warning("bulk ATAC DA figures failed for %s: %s",
                            (comp or {}).get("test"), e)
        # Sample-level QC: replicate-correlation heatmap from the peak-count matrix
        # used for DA (honest-skip when the matrix is missing/degenerate).
        try:
            path = render_sample_correlation(
                da.get("replicate_counts_path"),
                output_dir / "bulk_atac_sample_correlation.png")
            if path:
                figs["bulk_atac_sample_correlation"] = path
        except Exception as e:
            log.warning("bulk ATAC sample-correlation figure failed: %s", e)

    # Genomic peak-annotation distribution (B2) + functional ORA dotplots (B3).
    # Wave 2A: `chromatin_peak_annotation` / `chromatin_peak_ora` also run for
    # scATAC, so surface them for every chromatin modality (not bulk-only).
    data_type = _findings_data_type(findings)
    ann = findings.get("peak_annotation")
    ann_dt = (ann.get("data_type") if isinstance(ann, dict) else None) or data_type
    ann_prefix = ("bulk_atac" if ann_dt in _BULK_CHROMATIN_MODALITIES
                  else "scatac")
    try:
        path = render_peak_annotation(
            ann, output_dir / f"{ann_prefix}_peak_annotation.png",
            data_type=ann_dt)
        if path:
            figs[f"{ann_prefix}_peak_annotation"] = path
    except Exception as e:
        log.warning("peak-annotation figure failed: %s", e)
    # ORA dotplots are rendered in the rna stack by the ORA script (seaborn lives
    # there); surface the returned paths (PNG + companion SVG).
    ora = findings.get("peak_ora") or {}
    for key, path in ((ora.get("figures") or {}) if isinstance(ora, dict)
                      else {}).items():
        if isinstance(path, str):
            figs[key] = path
    for key, path in ((ora.get("figures_svg") or {}) if isinstance(ora, dict)
                      else {}).items():
        if isinstance(path, str):
            figs_svg[key] = path

    # 2b. Bulk ATAC sample PCA + ENCODE QC summary (ADR-077), rendered by
    # `chromatin_figure_bulk` in the rna stack (or in-process without an env
    # manager). Each panel self-skips when its data is absent/degenerate.
    if _is_bulk_atac_findings(findings):
        try:
            figs.update(_render_bulk_atac_figures(findings, output_dir,
                                                  env_manager))
        except Exception as e:  # a figure failure must not abort the report
            log.warning("bulk ATAC PCA/QC figures failed: %s", e)

    # 3. Motif enrichment dotplot (inline) from structured motif findings.
    motifs = findings.get("motifs") or {}
    if motifs.get("per_group"):
        path = render_motif_dotplot(motifs, output_dir / "motif_dotplot.png")
        if path:
            figs["motif_dotplot"] = path

    # 4. Fragment-size / nucleosome banding (inline) from the QC size histogram.
    qc = findings.get("qc") or {}
    qc_data_type = qc.get("data_type") or data_type
    path = render_fragment_size_figure(qc, output_dir / "fragment_size.png",
                                       data_type=qc_data_type)
    if path:
        figs["fragment_size"] = path

    # 4b. Per-cell QC gating scatter (TSSe × depth) + per-barcode FRiP distribution
    # (P4.1 (2)). Both honest-skip unless a real fragments/peaks run populated the
    # additive QC arrays.
    path = render_tss_depth_scatter(qc, output_dir / "qc_tss_depth.png")
    if path:
        figs["qc_tss_depth"] = path
    path = render_frip_distribution(
        qc, output_dir / "qc_frip_distribution.png", data_type=qc_data_type,
        thresholds=qc.get("thresholds") if isinstance(qc.get("thresholds"),
                                                      dict) else None)
    if path:
        figs["qc_frip_distribution"] = path

    # 5. Marker-peak heatmap + per-cluster QC depth violin (W0.2). These need the
    # AnnData matrix, so they render in the chromatin stack like the UMAP.
    if h5ad_path and env_manager is not None:
        try:
            res = env_manager.run_in_stack(
                stack="chromatin",
                script_path="aria/scripts/chromatin_figure_clusters.py",
                params={
                    "h5ad_path":    str(h5ad_path),
                    "output_dir":   str(output_dir),
                    "cluster_key":  _umap_color_keys(findings)[0],
                    "marker_peaks": _marker_peaks_by_cluster(findings),
                },
            )
            if res.get("status") == "success":
                for key, path in (res.get("figures") or {}).items():
                    if isinstance(path, str):
                        figs[key] = path
                for key, path in (res.get("figures_svg") or {}).items():
                    if isinstance(path, str):
                        figs_svg[key] = path
            else:
                log.warning(
                    "chromatin cluster figures failed: %s — %s",
                    res.get("error_type"), str(res.get("details", ""))[:200],
                )
        except Exception as e:
            log.warning("chromatin cluster figure subprocess crashed: %s", e)

    # 6. Differential TF footprinting (ADR-078): the volcano + per-TF aggregate
    # footprints are rendered by the TOBIAS backend; surface only files that exist.
    footprinting = findings.get("footprinting") or {}
    if isinstance(footprinting, dict) and footprinting.get("ran"):
        volcano_plot = footprinting.get("volcano_plot") or {}
        volcano_png = volcano_plot.get("png")
        if isinstance(volcano_png, str) and Path(volcano_png).exists():
            figs["footprint_volcano"] = volcano_png
        volcano_svg = volcano_plot.get("svg")
        if isinstance(volcano_svg, str) and Path(volcano_svg).exists():
            figs_svg["footprint_volcano"] = volcano_svg
        for tf, paths in (footprinting.get("aggregate_plots") or {}).items():
            png = (paths or {}).get("png") if isinstance(paths, dict) else None
            svg = (paths or {}).get("svg") if isinstance(paths, dict) else None
            fk = f"footprint_{_safe(tf)}"
            if isinstance(png, str) and Path(png).exists():
                figs[fk] = png
            if isinstance(svg, str) and Path(svg).exists():
                figs_svg[fk] = svg

    # Record companion sibling SVGs emitted by _save_dual / figure scripts.
    for key, png_path in figs.items():
        if key not in figs_svg and isinstance(png_path, str) and png_path:
            svg_cand = Path(png_path).with_suffix(".svg")
            if svg_cand.is_file():
                figs_svg[key] = str(svg_cand)

    return findings
