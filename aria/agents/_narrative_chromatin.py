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


def render_da_figures(da_full_csv: str, output_dir: Path,
                      padj_max: float = 0.05) -> dict:
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
        lfc = sub["log2fc"].to_numpy(dtype=float)
        padj = sub["padj"].to_numpy(dtype=float)
        sig = sub["significant"].fillna(False).to_numpy(dtype=bool)
        neglog = -np.log10(np.clip(padj, 1e-300, None))

        # ── Volcano ───────────────────────────────────────────────────────
        try:
            fig, ax = plt.subplots(figsize=(6.5, 6), dpi=160)
            ax.scatter(lfc[~sig], neglog[~sig], s=4, c="#bdbdbd",
                       linewidths=0, alpha=0.5, label="not significant")
            up, down = sig & (lfc > 0), sig & (lfc < 0)
            ax.scatter(lfc[up], neglog[up], s=7, c="#c0392b", linewidths=0,
                       alpha=0.85, label="up")
            ax.scatter(lfc[down], neglog[down], s=7, c="#2c7fb8", linewidths=0,
                       alpha=0.85, label="down")
            ax.axhline(-np.log10(padj_max), color="black", lw=0.5, ls="--")
            ax.axvline(0, color="black", lw=0.4)
            ax.set_xlabel("log2 fold change (accessibility)", fontsize=9)
            ax.set_ylabel("-log10 padj", fontsize=9)
            ax.set_title(f"Differential accessibility — {comp}", fontsize=10,
                         fontweight="bold")
            ax.legend(loc="upper right", fontsize=7, frameon=False,
                      markerscale=2.0)
            ax.tick_params(labelsize=7)
            path = output_dir / f"da_volcano_{comp_s}.png"
            _save_dual(fig, path)
            plt.close(fig)
            out[f"da_volcano_{comp_s}"] = str(path)
        except Exception as e:
            log.warning("DA volcano failed for %s: %s", comp, e)

        # ── MA (needs baseMean) ───────────────────────────────────────────
        if "base_mean" in sub.columns and not sub["base_mean"].isna().all():
            try:
                base = sub["base_mean"].to_numpy(dtype=float)
                x = np.log10(np.clip(base, 1e-3, None))
                fig, ax = plt.subplots(figsize=(6.5, 5), dpi=160)
                ax.scatter(x[~sig], lfc[~sig], s=4, c="#bdbdbd",
                           linewidths=0, alpha=0.5)
                ax.scatter(x[sig], lfc[sig], s=7, c="#c0392b", linewidths=0,
                           alpha=0.85, label="significant")
                ax.axhline(0, color="black", lw=0.5)
                ax.set_xlabel("log10 mean accessibility (baseMean)", fontsize=9)
                ax.set_ylabel("log2 fold change", fontsize=9)
                ax.set_title(f"MA plot — {comp}", fontsize=10, fontweight="bold")
                ax.legend(loc="upper right", fontsize=7, frameon=False,
                          markerscale=2.0)
                ax.tick_params(labelsize=7)
                path = output_dir / f"da_ma_{comp_s}.png"
                _save_dual(fig, path)
                plt.close(fig)
                out[f"da_ma_{comp_s}"] = str(path)
            except Exception as e:
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
    im = ax.imshow(corr, cmap="viridis", vmin=min(0.0, float(corr.min())),
                   vmax=1.0)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(samples, fontsize=7)
    thresh = (corr.max() + corr.min()) / 2.0
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


def render_peak_annotation(annotation: dict, png_path) -> Optional[str]:
    """Bulk ATAC genomic-context figure: horizontal bar chart of the DA-peak feature
    distribution (Promoter/Exonic/Intronic/Distal Intergenic) from
    `chromatin_peak_annotation`. Dual PNG+SVG. Honest-None when no annotation ran or
    the distribution is empty."""
    dist = (annotation or {}).get("feature_distribution_overall") or {}
    dist = {k: int(v) for k, v in dist.items() if v}
    if not dist:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

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
    ax.set_xlabel("DA peaks", fontsize=9)
    ax.set_title("Genomic distribution of differential-accessibility peaks",
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
    per_group = (motifs or {}).get("per_group") or {}
    rows = []
    for group, info in per_group.items():
        for m in (info or {}).get("top_motifs", [])[:top_n]:
            padj = m.get("padj")
            if padj is None:
                continue
            rows.append({
                "group": str(group),
                "motif": str(m.get("name") or m.get("motif_id") or "?"),
                "log2_enrichment": m.get("log2_enrichment"),
                "padj": float(padj),
            })
    if not rows:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    df = pd.DataFrame(rows)
    groups = sorted(df["group"].unique())
    motif_order = list(dict.fromkeys(df["motif"]))
    gx = {g: i for i, g in enumerate(groups)}
    my = {m: i for i, m in enumerate(motif_order)}
    try:
        fig, ax = plt.subplots(
            figsize=(max(4.0, 1.1 * len(groups) + 2),
                     max(3.0, 0.32 * len(motif_order) + 1)), dpi=160)
        sizes = 12 + 40 * (-np.log10(np.clip(df["padj"].to_numpy(float),
                                             1e-300, 1.0)))
        colors = df["log2_enrichment"].to_numpy(dtype=float)
        sc = ax.scatter([gx[g] for g in df["group"]],
                        [my[m] for m in df["motif"]],
                        s=sizes, c=colors, cmap="viridis",
                        edgecolors="black", linewidths=0.3)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(motif_order)))
        ax.set_yticklabels(motif_order, fontsize=7)
        ax.set_title("Motif enrichment (size=-log10 padj)", fontsize=10,
                     fontweight="bold")
        plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02, label="log2 enrichment")
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path = _save_dual(fig, out_path)
        plt.close(fig)
        return out_path
    except Exception as e:
        log.warning("motif dotplot failed: %s", e)
        return None


def render_fragment_size_figure(qc: dict, out_path: Path) -> Optional[str]:
    """W0.2: fragment-size / nucleosome-banding plot from the QC size histogram
    (`qc["fragment_sizes"]["size_histogram"]`). Rendered INLINE. Returns None when
    no histogram was computed (e.g. a `.h5mu` peak-matrix run with no fragments) —
    honest absence, never a fabricated distribution."""
    hist = ((qc or {}).get("fragment_sizes") or {}).get("size_histogram") or {}
    edges = hist.get("bin_edges") or []
    counts = hist.get("counts") or []
    if len(edges) < 2 or not counts or sum(counts) == 0:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    edges = np.asarray(edges, dtype=float)
    counts = np.asarray(counts, dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    n = min(len(centers), len(counts))
    try:
        fig, ax = plt.subplots(figsize=(6.5, 4), dpi=160)
        ax.fill_between(centers[:n], counts[:n], step="mid", alpha=0.6,
                        color="#34699a")
        ax.plot(centers[:n], counts[:n], lw=0.8, color="#1f4068")
        for x in (147, 294):   # mono- / di-nucleosome guides
            ax.axvline(x, color="grey", lw=0.5, ls="--")
        ax.set_xlabel("fragment size (bp)", fontsize=9)
        ax.set_ylabel("fragment count", fontsize=9)
        ax.set_title("Fragment-size distribution (nucleosome banding)",
                     fontsize=10, fontweight="bold")
        ax.tick_params(labelsize=7)
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path = _save_dual(fig, out_path)
        plt.close(fig)
        return out_path
    except Exception as e:
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
    if not tsse:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    tsse = np.asarray(tsse, dtype=float)
    depth = sc.get("log10_depth")
    min_tss = sc.get("min_tss")
    warn_tss = sc.get("warn_tss")
    try:
        fig, ax = plt.subplots(figsize=(6.0, 4.5), dpi=160)
        if depth and len(depth) == len(tsse):
            x = np.asarray(depth, dtype=float)
            xlabel = "log10 fragment depth"
        else:
            # No per-cell depth — jitter on x so points are separable.
            rng = np.random.default_rng(0)
            x = rng.normal(0.0, 0.04, size=tsse.shape[0])
            xlabel = "(no per-cell depth recorded)"
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
        log.warning("TSS×depth scatter figure failed: %s", e)
        return None


def render_frip_distribution(qc: dict, out_path: Path) -> Optional[str]:
    """P4.1 (2)[B]: per-barcode FRiP distribution histogram with the FRiP quality
    gate. Reads the additive ``qc["frip_distribution"]`` list (only present when a
    called-peak BED was supplied, W0.5 path). Returns None when no distribution was
    computed — honest absence, never a fabricated histogram."""
    dist = (qc or {}).get("frip_distribution")
    if not dist:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    vals = np.asarray(dist, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return None
    try:
        fig, ax = plt.subplots(figsize=(6.0, 4.0), dpi=160)
        ax.hist(vals, bins=40, range=(0, 1), color="#34699a", alpha=0.75)
        ax.axvline(0.2, color="#c0392b", lw=0.8, ls="--", label="FRiP gate (0.20)")
        ax.set_xlabel("FRiP (fraction of reads in peaks)", fontsize=9)
        ax.set_ylabel("barcodes", fontsize=9)
        ax.set_title(f"Per-barcode FRiP distribution (n={vals.size})",
                     fontsize=10, fontweight="bold")
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7, frameon=False)
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path = _save_dual(fig, out_path)
        plt.close(fig)
        return out_path
    except Exception as e:
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
    # scATAC shape: a single pseudobulk table.
    pb = da.get("pseudobulk") or {}
    da_full_csv = pb.get("full_results_csv")
    if da_full_csv:
        try:
            figs.update(render_da_figures(da_full_csv, output_dir, padj_max))
        except Exception as e:
            log.warning("chromatin DA figures failed: %s", e)
    # Bulk ATAC shape: one convergence-gated full CSV per condition comparison.
    if da.get("data_type") == "bulk_ATAC":
        for comp in (da.get("comparisons") or []):
            comp_csv = (comp or {}).get("full_results_csv")
            if not comp_csv:
                continue
            try:
                figs.update(render_da_figures(comp_csv, output_dir, padj_max))
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
        # Genomic peak-annotation distribution (B2): feature class bar chart.
        try:
            path = render_peak_annotation(
                findings.get("peak_annotation"),
                output_dir / "bulk_atac_peak_annotation.png")
            if path:
                figs["bulk_atac_peak_annotation"] = path
        except Exception as e:
            log.warning("bulk ATAC peak-annotation figure failed: %s", e)
        # Functional ORA dotplots (B3): rendered in the rna stack by the ORA script
        # (seaborn lives there); surface the returned paths.
        ora = findings.get("peak_ora") or {}
        for key, path in (ora.get("figures") or {}).items():
            if isinstance(path, str):
                figs[key] = path

    # 3. Motif enrichment dotplot (inline) from structured motif findings.
    motifs = findings.get("motifs") or {}
    if motifs.get("per_group"):
        path = render_motif_dotplot(motifs, output_dir / "motif_dotplot.png")
        if path:
            figs["motif_dotplot"] = path

    # 4. Fragment-size / nucleosome banding (inline) from the QC size histogram.
    qc = findings.get("qc") or {}
    path = render_fragment_size_figure(qc, output_dir / "fragment_size.png")
    if path:
        figs["fragment_size"] = path

    # 4b. Per-cell QC gating scatter (TSSe × depth) + per-barcode FRiP distribution
    # (P4.1 (2)). Both honest-skip unless a real fragments/peaks run populated the
    # additive QC arrays.
    path = render_tss_depth_scatter(qc, output_dir / "qc_tss_depth.png")
    if path:
        figs["qc_tss_depth"] = path
    path = render_frip_distribution(qc, output_dir / "qc_frip_distribution.png")
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
            else:
                log.warning(
                    "chromatin cluster figures failed: %s — %s",
                    res.get("error_type"), str(res.get("details", ""))[:200],
                )
        except Exception as e:
            log.warning("chromatin cluster figure subprocess crashed: %s", e)

    return findings
