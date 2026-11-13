"""ARIA 3D Genome Architecture Publication Figures (ADR-082).

validation_level: production (stack="hic", headless matplotlib).

Renders up to four dual-format (300-DPI ``.png`` + editable ``.svg``) figures
exclusively from measured Hi-C outputs (zero fabrication — any panel lacking
input data is omitted and recorded under ``skipped``):

* ``hic_contact_decay_ps`` — log-log contact probability ``P(s)`` vs genomic
  separation ``s`` + log-derivative polymer slope panel.
* ``hic_compartment_insulation_track`` — stacked chromosome ``E1`` A/B
  eigenvector track (+A coral / -B slate blue) and Crane ``log2`` insulation
  profile with prominent TAD boundaries marked.
* ``hic_contact_map`` — out-of-core region slice of the ``log10`` ICE-balanced
  contact matrix with TAD domain overlays (bounded to <= 400 bins so RAM stays
  strictly < 50 MB).
* ``hic_differential_compartments`` — volcano / effect plot of ``delta_E1`` vs
  ``-log10(padj)`` highlighting ``A_to_B`` and ``B_to_A`` switches when
  multi-replicate ``differential_3d`` ran.

Zero imports from ``aria.agents.*``, ``litellm``, or ``pydantic``.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from aria.scripts._base import run_script

_SURFACE = "#fcfcfb"
_TEXT = "#111111"
_SUBTEXT = "#52514e"
_COLOR_A = "#d6453d"
_COLOR_B = "#2a78d6"
_PALETTE = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6", "#e34948")


def _save_dual(fig, png_path: Path) -> dict[str, str]:
    png_path = Path(png_path)
    svg_path = png_path.with_suffix(".svg")
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor=_SURFACE)
    fig.savefig(svg_path, bbox_inches="tight", facecolor=_SURFACE)
    return {"png": str(png_path), "svg": str(svg_path)}


def _style_ax(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#b5b4ae")
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=_SUBTEXT, labelsize=8)
    ax.set_facecolor(_SURFACE)


def _render_ps_decay(ps_files: list[tuple[str, str]], out_dir: Path) -> dict[str, str] | None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    loaded: list[tuple[str, pd.DataFrame]] = []
    for label, p in ps_files:
        if p and Path(p).is_file():
            df = pd.read_csv(p, sep="\t")
            sub = df[(df["dist_bp"] > 0) & np.isfinite(df["contact_prob"]) & (df["contact_prob"] > 0)]
            if len(sub) >= 3:
                loaded.append((label, sub.sort_values("dist_bp")))
    if not loaded:
        return None

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6.2, 5.4), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.0]}
    )
    for idx, (label, sub) in enumerate(loaded):
        col = _PALETTE[idx % len(_PALETTE)]
        s_kb = sub["dist_bp"].to_numpy(dtype=float) / 1000.0
        prob = sub["contact_prob"].to_numpy(dtype=float)
        ax1.plot(s_kb, prob, lw=1.8, color=col, label=label)
        if len(s_kb) >= 4:
            lx = np.log10(s_kb)
            ly = np.log10(prob)
            deriv = np.gradient(ly, lx)
            ax2.plot(s_kb, deriv, lw=1.3, color=col, alpha=0.85)

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylabel("Contact probability P(s)", fontsize=9, color=_TEXT)
    ax1.set_title("3D Contact Decay P(s) & Polymer Scaling Derivative", fontsize=10, color=_TEXT)
    if len(loaded) <= 8:
        ax1.legend(frameon=False, fontsize=7.5)
    _style_ax(ax1)

    ax2.axhline(-1.0, ls="--", lw=0.8, color="#888888")
    ax2.set_xscale("log")
    ax2.set_xlabel("Genomic separation s (kb)", fontsize=9, color=_TEXT)
    ax2.set_ylabel("d log P(s) / d log s", fontsize=8.5, color=_TEXT)
    ax2.set_ylim(-2.5, 0.5)
    _style_ax(ax2)

    fig.tight_layout()
    paths = _save_dual(fig, out_dir / "hic_contact_decay_ps.png")
    plt.close(fig)
    return paths


def _render_compartment_insulation_track(
    compartments_tsv: str | None,
    insulation_tsv: str | None,
    out_dir: Path,
) -> dict[str, str] | None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    comp_df = (
        pd.read_csv(compartments_tsv, sep="\t")
        if compartments_tsv and Path(compartments_tsv).is_file() else None
    )
    ins_df = (
        pd.read_csv(insulation_tsv, sep="\t")
        if insulation_tsv and Path(insulation_tsv).is_file() else None
    )
    if comp_df is None and ins_df is None:
        return None

    chrom = None
    if comp_df is not None and not comp_df.empty:
        valid_c = comp_df[comp_df["E1"].notna()]
        if not valid_c.empty:
            chrom = str(valid_c["chrom"].iloc[0])
    if chrom is None and ins_df is not None and not ins_df.empty:
        chrom = str(ins_df["chrom"].iloc[0])
    if chrom is None:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(7.4, 4.4), sharex=True)

    # Top panel: E1 eigenvector track
    ax_e1 = axes[0]
    if comp_df is not None and not comp_df.empty:
        sub_c = comp_df[comp_df["chrom"].astype(str) == chrom].copy()
        mid_mb = 0.5 * (sub_c["start"].to_numpy(dtype=float) + sub_c["end"].to_numpy(dtype=float)) / 1e6
        e1 = sub_c["E1"].to_numpy(dtype=float)
        width_mb = max(np.median(sub_c["end"] - sub_c["start"]) / 1e6, 0.01) if len(sub_c) else 0.05
        pos = np.where(np.isfinite(e1) & (e1 > 0), e1, 0.0)
        neg = np.where(np.isfinite(e1) & (e1 < 0), e1, 0.0)
        ax_e1.bar(mid_mb, pos, width=width_mb, color=_COLOR_A, alpha=0.85, label="A (+E1)")
        ax_e1.bar(mid_mb, neg, width=width_mb, color=_COLOR_B, alpha=0.85, label="B (-E1)")
        ax_e1.axhline(0.0, lw=0.7, color="#444444")
        ax_e1.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax_e1.set_ylabel("Eigenvector E1", fontsize=8.5, color=_TEXT)
    ax_e1.set_title(f"Compartment E1 & Insulation Profile ({chrom})", fontsize=10, color=_TEXT)
    _style_ax(ax_e1)

    # Bottom panel: Crane log2 insulation score + boundaries
    ax_ins = axes[1]
    if ins_df is not None and not ins_df.empty:
        sub_i = ins_df[ins_df["chrom"].astype(str) == chrom].copy()
        sc_cols = [c for c in sub_i.columns if c.startswith("log2_insulation_score_")]
        bd_cols = [c for c in sub_i.columns if c.startswith("is_boundary_")]
        if sc_cols and not sub_i.empty:
            sc_col = sc_cols[0]
            mid_mb = 0.5 * (sub_i["start"].to_numpy(dtype=float) + sub_i["end"].to_numpy(dtype=float)) / 1e6
            vals = sub_i[sc_col].to_numpy(dtype=float)
            ax_ins.plot(mid_mb, vals, lw=1.2, color="#222222")
            if bd_cols:
                b_mask = sub_i[bd_cols[0]].fillna(False).astype(bool).to_numpy()
                for bx in mid_mb[b_mask]:
                    ax_ins.axvline(bx, lw=0.6, ls=":", color="#eb6834", alpha=0.75)
    ax_ins.set_xlabel(f"{chrom} position (Mb)", fontsize=9, color=_TEXT)
    ax_ins.set_ylabel("log2 insulation", fontsize=8.5, color=_TEXT)
    _style_ax(ax_ins)

    fig.tight_layout()
    paths = _save_dual(fig, out_dir / "hic_compartment_insulation_track.png")
    plt.close(fig)
    return paths


def _render_contact_map(
    cooler_path: str | None,
    resolution: int,
    tads_bed: str | None,
    out_dir: Path,
    max_bins: int = 300,
) -> dict[str, str] | None:
    if not cooler_path or not Path(cooler_path).exists():
        return None
    import cooler
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    if cooler_path.lower().endswith(".mcool"):
        avail = sorted(
            int(p.split("/")[-1])
            for p in cooler.fileops.list_coolers(cooler_path)
            if p.split("/")[-1].isdigit()
        )
        if not avail:
            return None
        chosen = int(resolution) if int(resolution) in avail else min(
            avail, key=lambda r: abs(r - int(resolution))
        )
        uri = f"{cooler_path}::resolutions/{chosen}"
    else:
        uri = cooler_path
        chosen = int(cooler.Cooler(uri).binsize)

    clr = cooler.Cooler(uri)
    chroms = [c for c in clr.chromnames if str(c) not in {"chrM", "MT", "M"}] or list(clr.chromnames)
    if not chroms:
        return None
    chrom = "chr1" if "chr1" in chroms else chroms[0]
    chrom_len = int(clr.chromsizes[chrom])
    region_end = min(chrom_len, max_bins * chosen)
    has_w = "weight" in list(clr.bins().columns)
    mat = clr.matrix(balance=has_w).fetch((chrom, 0, region_end))
    if mat.size == 0:
        return None

    arr = np.asarray(mat, dtype=float)
    pos = arr[np.isfinite(arr) & (arr > 0)]
    if pos.size == 0 and has_w:
        arr = np.asarray(clr.matrix(balance=False).fetch((chrom, 0, region_end)), dtype=float)
        pos = arr[np.isfinite(arr) & (arr > 0)]
    if pos.size == 0:
        return None

    floor = float(np.percentile(pos, 5)) * 0.5
    log_mat = np.log10(np.where(np.isfinite(arr) & (arr > 0), arr, max(floor, 1e-8)))

    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    extent_mb = [0.0, region_end / 1e6, region_end / 1e6, 0.0]
    im = ax.imshow(
        log_mat,
        cmap="YlOrRd",
        extent=extent_mb,
        interpolation="nearest",
    )
    if tads_bed and Path(tads_bed).is_file() and Path(tads_bed).stat().st_size > 0:
        tdf = pd.read_csv(tads_bed, sep="\t", header=None)
        if tdf.shape[1] >= 3:
            sub_t = tdf[(tdf[0].astype(str) == str(chrom)) & (tdf[2].astype(float) <= region_end)]
            for _, row in sub_t.iterrows():
                s_mb = float(row[1]) / 1e6
                e_mb = float(row[2]) / 1e6
                ax.plot([s_mb, e_mb, e_mb, s_mb, s_mb], [s_mb, s_mb, e_mb, e_mb, s_mb],
                        lw=0.8, color="#111111", alpha=0.75)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("log10 balanced contact frequency", fontsize=8)
    ax.set_xlabel(f"{chrom} (Mb)", fontsize=9, color=_TEXT)
    ax.set_ylabel(f"{chrom} (Mb)", fontsize=9, color=_TEXT)
    ax.set_title(f"ICE-Balanced Contact Map ({chrom}, {chosen // 1000} kb)", fontsize=9.5, color=_TEXT)
    fig.tight_layout()
    paths = _save_dual(fig, out_dir / "hic_contact_map.png")
    plt.close(fig)
    return paths


def _render_differential_compartments(diff_csv: str | None, out_dir: Path) -> dict[str, str] | None:
    if not diff_csv or not Path(diff_csv).is_file():
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    df = pd.read_csv(diff_csv)
    if df.empty or "delta_E1" not in df.columns or "padj" not in df.columns:
        return None
    sub = df[df["delta_E1"].notna() & df["padj"].notna()].copy()
    if sub.empty:
        return None

    neglog = -np.log10(np.clip(sub["padj"].to_numpy(dtype=float), 1e-300, 1.0))
    dx = sub["delta_E1"].to_numpy(dtype=float)
    cls = sub.get("switch_class", pd.Series(["unchanged"] * len(sub))).astype(str).to_numpy()

    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    bg = cls == "unchanged"
    b2a = cls == "B_to_A"
    a2b = cls == "A_to_B"
    qsh = cls == "quantitative_shift"

    ax.scatter(dx[bg], neglog[bg], s=12, color="#b8b7b2", alpha=0.55, label="Unchanged")
    if np.any(qsh):
        ax.scatter(dx[qsh], neglog[qsh], s=16, color="#eda100", alpha=0.8, label="Quantitative E1 shift")
    if np.any(b2a):
        ax.scatter(dx[b2a], neglog[b2a], s=22, color=_COLOR_A, alpha=0.9, label="B -> A switch")
    if np.any(a2b):
        ax.scatter(dx[a2b], neglog[a2b], s=22, color=_COLOR_B, alpha=0.9, label="A -> B switch")

    ax.axhline(-np.log10(0.05), ls="--", lw=0.8, color="#666666")
    ax.axvline(0.20, ls=":", lw=0.8, color="#666666")
    ax.axvline(-0.20, ls=":", lw=0.8, color="#666666")
    ax.set_xlabel("Compartment ΔE1 (test - reference)", fontsize=9, color=_TEXT)
    ax.set_ylabel("-log10(BH padj)", fontsize=9, color=_TEXT)
    ax.set_title("Differential 3D Compartment Switching", fontsize=10, color=_TEXT)
    ax.legend(frameon=False, fontsize=7.5)
    _style_ax(ax)
    fig.tight_layout()
    paths = _save_dual(fig, out_dir / "hic_differential_compartments.png")
    plt.close(fig)
    return paths


def hic_figure(params: dict) -> dict:
    out_dir = Path(str(params.get("output_dir") or "/tmp/aria_hic_figures"))
    out_dir.mkdir(parents=True, exist_ok=True)

    figures: dict[str, str] = {}
    figure_bundles: dict[str, dict[str, str]] = {}
    skipped: dict[str, str] = {}

    # 1. P(s) decay curve
    ps_inputs: list[tuple[str, str]] = []
    for rec in params.get("per_file") or []:
        if isinstance(rec, dict) and rec.get("ps_curve_tsv"):
            ps_inputs.append((str(rec.get("sample_id") or "sample"), str(rec["ps_curve_tsv"])))
    if not ps_inputs and params.get("ps_curve_tsv"):
        ps_inputs.append(("cohort", str(params["ps_curve_tsv"])))
    ps_bundle = _render_ps_decay(ps_inputs, out_dir)
    if ps_bundle:
        figures["hic_contact_decay_ps"] = ps_bundle["png"]
        figure_bundles["hic_contact_decay_ps"] = ps_bundle
    else:
        skipped["hic_contact_decay_ps"] = "no_ps_curve_tsv"

    # 2. Compartment + Insulation track
    track_bundle = _render_compartment_insulation_track(
        params.get("compartments_tsv"),
        params.get("insulation_tsv"),
        out_dir,
    )
    if track_bundle:
        figures["hic_compartment_insulation_track"] = track_bundle["png"]
        figure_bundles["hic_compartment_insulation_track"] = track_bundle
    else:
        skipped["hic_compartment_insulation_track"] = "no_compartment_or_insulation_tsv"

    # 3. Contact map with TAD overlay
    map_bundle = _render_contact_map(
        params.get("cooler_path"),
        int(params.get("resolution") or 40_000),
        params.get("tads_bed"),
        out_dir,
    )
    if map_bundle:
        figures["hic_contact_map"] = map_bundle["png"]
        figure_bundles["hic_contact_map"] = map_bundle
    else:
        skipped["hic_contact_map"] = "no_cooler_contact_matrix"

    # 4. Differential compartments volcano
    diff_bundle = _render_differential_compartments(
        params.get("differential_compartments_csv"),
        out_dir,
    )
    if diff_bundle:
        figures["hic_differential_compartments"] = diff_bundle["png"]
        figure_bundles["hic_differential_compartments"] = diff_bundle
    else:
        skipped["hic_differential_compartments"] = "no_differential_compartments_csv"

    return {
        "status": "success",
        "ran": bool(figures),
        "data_type": "HiC",
        "analysis": "hic_figure",
        "validation_level": "production",
        "figures": figures,
        "figure_bundles": figure_bundles,
        "skipped": skipped,
        "output_dir": str(out_dir),
    }


if __name__ == "__main__":
    run_script(hic_figure)
