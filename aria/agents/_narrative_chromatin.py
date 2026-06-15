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
    for candidate in (groupby, "leiden", "sample", "sample_id", "batch",
                      "condition"):
        if candidate and candidate not in keys:
            keys.append(candidate)
    return keys


def _safe(value) -> str:
    return (str(value).replace("/", "_").replace(" ", "_")
            .replace(":", "_").replace("|", "_"))


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
            fig.savefig(path, dpi=160, bbox_inches="tight", facecolor="white")
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
                fig.savefig(path, dpi=160, bbox_inches="tight", facecolor="white")
                plt.close(fig)
                out[f"da_ma_{comp_s}"] = str(path)
            except Exception as e:
                log.warning("DA MA failed for %s: %s", comp, e)
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
        fig.savefig(out_path, dpi=160, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        return str(out_path)
    except Exception as e:
        log.warning("motif dotplot failed: %s", e)
        return None


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
    pb = da.get("pseudobulk") or {}
    da_full_csv = pb.get("full_results_csv")
    if da_full_csv:
        try:
            figs.update(render_da_figures(
                da_full_csv, output_dir, float(da.get("padj_max", 0.05))))
        except Exception as e:
            log.warning("chromatin DA figures failed: %s", e)

    # 3. Motif enrichment dotplot (inline) from structured motif findings.
    motifs = findings.get("motifs") or {}
    if motifs.get("per_group"):
        path = render_motif_dotplot(motifs, output_dir / "motif_dotplot.png")
        if path:
            figs["motif_dotplot"] = path

    return findings
