"""ARIA 3D Genome Architecture Publication Figures (ADR-082).

validation_level: production (stack="hic", headless matplotlib).

Renders up to four dual-format (300-DPI ``.png`` + editable ``.svg``) figures
exclusively from measured Hi-C outputs (zero fabrication — any panel lacking
input data is omitted and recorded under ``skipped``):

* ``hic_contact_decay_ps`` — log-log contact probability ``P(s)`` vs genomic
  separation ``s`` + smoothed log-derivative polymer slope panel.
* ``hic_compartment_insulation_track`` — stacked chromosome ``E1`` A/B
  eigenvector track (+A coral / -B slate blue, or unphased polarity when
  ``pc1_validated=False``) and Crane ``log2`` insulation profile with TAD
  boundaries marked.
* ``hic_contact_map`` — out-of-core region slice of the first mappable
  ``log10`` contact matrix window with TAD domain overlays (bounded to <= 400
  bins so RAM stays strictly < 50 MB).
* ``hic_differential_compartments`` — volcano / effect plot of ``delta_E1`` vs
  ``-log10(padj)`` highlighting ``A_to_B``, ``B_to_A``, and unphased sign flips
  when multi-replicate ``differential_3d`` ran.

Zero imports from ``aria.agents.*``, ``litellm``, or ``pydantic``.
"""

from __future__ import annotations

import math
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
_COLOR_UNPHASED = "#6250d6"
_PALETTE = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6", "#e34948")
EXCLUDED_CHROMS = frozenset({"chrM", "MT", "chrMT", "M", "chrEBV", "EBV"})


def _chrom_aliases(chrom: str) -> tuple[str, ...]:
    c = str(chrom).strip()
    if not c:
        return ()
    if c.lower().startswith("chr"):
        bare = c[3:]
        if bare.upper() in {"M", "MT"}:
            return (c, "chrM", "chrMT", "MT", "M")
        return (c, bare)
    if c.upper() in {"M", "MT"}:
        return (c, "chrM", "chrMT", "MT", "M")
    return (c, f"chr{c}")


def _pick_display_chrom(
    available_chroms: list[str],
    requested: Any = None,
) -> str | None:
    """Choose display chromosome honoring explicit `chrom`/`chromosomes` and `EXCLUDED_CHROMS`."""
    all_chroms = [str(c) for c in available_chroms if str(c).strip()]
    if not all_chroms:
        return None
    nuclear = [
        c for c in all_chroms
        if c not in EXCLUDED_CHROMS
        and "random" not in c.lower()
        and not c.lower().startswith("chrun")
    ] or all_chroms

    req_list: list[str] = []
    if requested is not None and requested != "all":
        if isinstance(requested, str):
            if requested.strip() and requested.strip().lower() != "all":
                req_list = [requested.strip()]
        elif isinstance(requested, (list, tuple, set)):
            req_list = [str(c).strip() for c in requested if str(c).strip()]
    if req_list:
        for req in req_list:
            aliases = set(_chrom_aliases(req))
            for cand in all_chroms:
                if cand in aliases:
                    return cand

    for preferred in ("chr1", "1"):
        if preferred in nuclear:
            return preferred
    return nuclear[0]


def _format_resolution(res_bp: int) -> str:
    res = int(res_bp)
    if res >= 1_000_000 and res % 1_000_000 == 0:
        return f"{res // 1_000_000} Mb"
    if res >= 1_000 and res % 1_000 == 0:
        return f"{res // 1_000} kb"
    if res >= 1_000:
        return f"{res / 1_000:.1f} kb"
    return f"{res} bp"


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


def _smoothed_log_derivative(s_kb, prob):
    """Compute d log10 P(s) / d log10 s cleanly without duplicate-x div0 or shot-noise spikes."""
    import numpy as np

    lx = np.log10(s_kb)
    ly = np.log10(prob)
    finite = np.isfinite(lx) & np.isfinite(ly)
    lx = lx[finite]
    ly = ly[finite]
    if lx.size < 4:
        return np.array([]), np.array([])

    # Bin into log-spaced intervals when curve is dense to suppress Poisson tail spikes.
    if lx.size >= 12 and (lx[-1] - lx[0]) > 1e-6:
        n_bins = min(60, max(10, lx.size // 2))
        edges = np.linspace(lx[0], lx[-1], n_bins + 1)
        bx, by = [], []
        for idx in range(n_bins):
            if idx == n_bins - 1:
                mask = (lx >= edges[idx]) & (lx <= edges[idx + 1])
            else:
                mask = (lx >= edges[idx]) & (lx < edges[idx + 1])
            if np.any(mask):
                bx.append(float(np.median(lx[mask])))
                by.append(float(np.median(ly[mask])))
        bx_arr = np.asarray(bx, dtype=float)
        by_arr = np.asarray(by, dtype=float)
        if bx_arr.size >= 4:
            # Ensure strictly increasing x before differentiating.
            keep = np.concatenate(([True], np.diff(bx_arr) > 1e-9))
            bx_arr = bx_arr[keep]
            by_arr = by_arr[keep]
            if bx_arr.size >= 4:
                deriv = np.gradient(by_arr, bx_arr)
                ok = np.isfinite(deriv)
                return np.power(10.0, bx_arr[ok]), deriv[ok]

    keep = np.concatenate(([True], np.diff(lx) > 1e-9))
    lx = lx[keep]
    ly = ly[keep]
    if lx.size < 4:
        return np.array([]), np.array([])
    deriv = np.gradient(ly, lx)
    ok = np.isfinite(deriv)
    return np.power(10.0, lx[ok]), deriv[ok]


def _render_ps_decay(ps_files: list[tuple[str, str]], out_dir: Path) -> dict[str, str] | None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    loaded: list[tuple[str, pd.DataFrame]] = []
    for label, p in ps_files:
        if not p:
            continue
        path = Path(p)
        if not path.is_file() or path.stat().st_size == 0:
            continue
        df = pd.read_csv(path, sep="\t")
        if "dist_bp" not in df.columns or "contact_prob" not in df.columns:
            continue
        dist = pd.to_numeric(df["dist_bp"], errors="coerce").to_numpy(dtype=float)
        prob = pd.to_numeric(df["contact_prob"], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(dist) & (dist > 0) & np.isfinite(prob) & (prob > 0)
        if int(valid.sum()) < 3:
            continue
        clean = pd.DataFrame({"dist_bp": dist[valid], "contact_prob": prob[valid]})
        # Deduplicate dist_bp so np.gradient never divides by zero.
        clean = clean.groupby("dist_bp", as_index=False)["contact_prob"].mean().sort_values("dist_bp")
        if len(clean) >= 3:
            loaded.append((label, clean))
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
            dx_kb, deriv = _smoothed_log_derivative(s_kb, prob)
            if dx_kb.size >= 2:
                ax2.plot(dx_kb, deriv, lw=1.3, color=col, alpha=0.85)

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
    pc1_validated: bool | None = None,
    chrom: Any = None,
) -> dict[str, str] | None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    comp_df = None
    if compartments_tsv and Path(compartments_tsv).is_file() and Path(compartments_tsv).stat().st_size > 0:
        cand = pd.read_csv(compartments_tsv, sep="\t")
        if {"chrom", "start", "end", "E1"}.issubset(cand.columns) and not cand.empty:
            e1_num = pd.to_numeric(cand["E1"], errors="coerce").to_numpy(dtype=float)
            if np.any(np.isfinite(e1_num)):
                comp_df = cand

    ins_df = None
    ins_sc_col = None
    ins_bd_col = None
    if insulation_tsv and Path(insulation_tsv).is_file() and Path(insulation_tsv).stat().st_size > 0:
        cand = pd.read_csv(insulation_tsv, sep="\t")
        sc_cols = [c for c in cand.columns if c.startswith("log2_insulation_score_")]
        bd_cols = [c for c in cand.columns if c.startswith("is_boundary_")]
        if {"chrom", "start", "end"}.issubset(cand.columns) and sc_cols and not cand.empty:
            sc_num = pd.to_numeric(cand[sc_cols[0]], errors="coerce").to_numpy(dtype=float)
            if np.any(np.isfinite(sc_num)):
                ins_df = cand
                ins_sc_col = sc_cols[0]
                ins_bd_col = bd_cols[0] if bd_cols else None

    if comp_df is None and ins_df is None:
        return None

    avail_chroms: list[str] = []
    if comp_df is not None:
        e1_vals = pd.to_numeric(comp_df["E1"], errors="coerce").to_numpy(dtype=float)
        valid_c = comp_df[np.isfinite(e1_vals)]
        for ch_val in valid_c["chrom"].astype(str).tolist():
            if ch_val not in avail_chroms:
                avail_chroms.append(ch_val)
    if ins_df is not None and ins_sc_col is not None:
        sc_vals = pd.to_numeric(ins_df[ins_sc_col], errors="coerce").to_numpy(dtype=float)
        valid_i = ins_df[np.isfinite(sc_vals)]
        for ch_val in valid_i["chrom"].astype(str).tolist():
            if ch_val not in avail_chroms:
                avail_chroms.append(ch_val)
    chosen_chrom = _pick_display_chrom(avail_chroms, requested=chrom)
    if chosen_chrom is None:
        return None
    chrom = chosen_chrom

    sub_c = None
    if comp_df is not None:
        tmp_c = comp_df[comp_df["chrom"].astype(str) == chrom].copy()
        if not tmp_c.empty:
            e1_c = pd.to_numeric(tmp_c["E1"], errors="coerce").to_numpy(dtype=float)
            if np.any(np.isfinite(e1_c)):
                sub_c = tmp_c

    sub_i = None
    if ins_df is not None and ins_sc_col is not None:
        tmp_i = ins_df[ins_df["chrom"].astype(str) == chrom].copy()
        if not tmp_i.empty:
            sc_i = pd.to_numeric(tmp_i[ins_sc_col], errors="coerce").to_numpy(dtype=float)
            if np.any(np.isfinite(sc_i)):
                sub_i = tmp_i

    n_panels = int(sub_c is not None) + int(sub_i is not None)
    if n_panels == 0:
        return None

    fig, axes = plt.subplots(
        n_panels, 1, figsize=(7.4, 2.5 * n_panels + 0.2), sharex=(n_panels > 1)
    )
    if n_panels == 1:
        axes = [axes]
    ax_idx = 0

    if sub_c is not None:
        ax_e1 = axes[ax_idx]
        ax_idx += 1
        starts = pd.to_numeric(sub_c["start"], errors="coerce").to_numpy(dtype=float)
        ends = pd.to_numeric(sub_c["end"], errors="coerce").to_numpy(dtype=float)
        mid_mb = 0.5 * (starts + ends) / 1e6
        e1 = pd.to_numeric(sub_c["E1"], errors="coerce").to_numpy(dtype=float)
        ok_coord = np.isfinite(mid_mb)
        mid_mb = mid_mb[ok_coord]
        e1_clean = np.where(np.isfinite(e1[ok_coord]), e1[ok_coord], np.nan)
        order = np.argsort(mid_mb, kind="mergesort")
        mid_mb = mid_mb[order]
        e1_clean = e1_clean[order]

        unphased = pc1_validated is not True
        pos_label = "+E1 (unphased)" if unphased else "A (+E1)"
        neg_label = "-E1 (unphased)" if unphased else "B (-E1)"
        ax_e1.fill_between(
            mid_mb, e1_clean, 0.0,
            where=np.isfinite(e1_clean) & (e1_clean > 0),
            interpolate=True, color=_COLOR_A, alpha=0.82, linewidth=0.0, label=pos_label,
        )
        ax_e1.fill_between(
            mid_mb, e1_clean, 0.0,
            where=np.isfinite(e1_clean) & (e1_clean < 0),
            interpolate=True, color=_COLOR_B, alpha=0.82, linewidth=0.0, label=neg_label,
        )
        ax_e1.axhline(0.0, lw=0.7, color="#444444")
        ax_e1.legend(frameon=False, fontsize=7.5, loc="upper right")
        ax_e1.set_ylabel("Eigenvector E1", fontsize=8.5, color=_TEXT)
        title_suffix = " (unphased polarity)" if unphased else ""
        if sub_i is not None:
            ax_e1.set_title(f"Compartment E1 & Insulation Profile ({chrom}){title_suffix}", fontsize=10, color=_TEXT)
        else:
            ax_e1.set_title(f"Compartment E1 Profile ({chrom}){title_suffix}", fontsize=10, color=_TEXT)
            ax_e1.set_xlabel(f"{chrom} position (Mb)", fontsize=9, color=_TEXT)
        _style_ax(ax_e1)

    if sub_i is not None and ins_sc_col is not None:
        ax_ins = axes[ax_idx]
        starts = pd.to_numeric(sub_i["start"], errors="coerce").to_numpy(dtype=float)
        ends = pd.to_numeric(sub_i["end"], errors="coerce").to_numpy(dtype=float)
        mid_mb = 0.5 * (starts + ends) / 1e6
        vals = pd.to_numeric(sub_i[ins_sc_col], errors="coerce").to_numpy(dtype=float)
        ok_coord = np.isfinite(mid_mb)
        mid_mb = mid_mb[ok_coord]
        vals = np.where(np.isfinite(vals[ok_coord]), vals[ok_coord], np.nan)
        order = np.argsort(mid_mb, kind="mergesort")
        mid_mb = mid_mb[order]
        vals = vals[order]
        ax_ins.plot(mid_mb, vals, lw=1.2, color="#222222")
        if ins_bd_col:
            b_raw = sub_i[ins_bd_col].fillna(False).astype(bool).to_numpy()[ok_coord][order]
            bx_arr = mid_mb[b_raw]
            by_arr = vals[b_raw]
            if bx_arr.size > 0:
                if bx_arr.size <= 80:
                    for bx in bx_arr:
                        ax_ins.axvline(bx, lw=0.6, ls=":", color="#eb6834", alpha=0.65)
                else:
                    # High boundary count (e.g. full chr1): mark on curve + rug ticks to avoid wall.
                    fin_b = np.isfinite(by_arr)
                    if np.any(fin_b):
                        ax_ins.scatter(
                            bx_arr[fin_b], by_arr[fin_b], s=8, color="#eb6834",
                            alpha=0.75, zorder=3, label=f"TAD boundaries (n={int(bx_arr.size)})",
                        )
                        ax_ins.legend(frameon=False, fontsize=7.5, loc="upper right")
        if sub_c is None:
            ax_ins.set_title(f"Crane Insulation Profile ({chrom})", fontsize=10, color=_TEXT)
        ax_ins.set_xlabel(f"{chrom} position (Mb)", fontsize=9, color=_TEXT)
        ax_ins.set_ylabel("log2 insulation", fontsize=8.5, color=_TEXT)
        _style_ax(ax_ins)

    fig.tight_layout()
    paths = _save_dual(fig, out_dir / "hic_compartment_insulation_track.png")
    plt.close(fig)
    return paths


def _pick_mappable_window(clr, chrom: str, chosen: int, max_bins: int, has_w: bool) -> tuple[int, int]:
    """Select first mappable genomic window on ``chrom`` (avoids mm10/mm39 0-3 Mb N-gap)."""
    import numpy as np

    chrom_len = int(clr.chromsizes[chrom])
    span = min(chrom_len, int(max_bins) * int(chosen))
    region_start: int | None = None

    try:
        bins_df = clr.bins().fetch(chrom)
        if has_w and "weight" in bins_df.columns and not bins_df.empty:
            w = np.asarray(bins_df["weight"], dtype=float)
            good = np.flatnonzero(np.isfinite(w) & (w > 0))
            if good.size > 0:
                region_start = int(bins_df["start"].iloc[int(good[0])])
        if region_start is None and not bins_df.empty:
            # Scan strictly cis sparse contacts on chrom (zero genome-wide RAM load).
            mat = clr.matrix(balance=False, sparse=True).fetch(chrom)
            coo = mat.tocoo()
            pos_mask = np.isfinite(coo.data) & (coo.data > 0)
            if np.any(pos_mask):
                first_local = min(
                    int(np.min(coo.row[pos_mask])),
                    int(np.min(coo.col[pos_mask])),
                )
                if 0 <= first_local < len(bins_df):
                    region_start = int(bins_df["start"].iloc[first_local])
    except Exception:
        region_start = 0

    if region_start is None:
        region_start = 0

    if region_start + span > chrom_len:
        region_start = max(0, chrom_len - span)
    region_end = min(chrom_len, region_start + span)
    return int(region_start), int(region_end)


def _render_contact_map(
    cooler_path: str | None,
    resolution: int,
    tads_bed: str | None,
    out_dir: Path,
    max_bins: int = 300,
    chrom: Any = None,
) -> dict[str, str] | None:
    if not cooler_path or not Path(cooler_path).exists():
        return None
    import cooler
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    from aria.utils.hic_constants import list_mcool_resolutions

    if cooler_path.lower().endswith(".mcool"):
        res_map = list_mcool_resolutions(cooler_path)
        if not res_map:
            return None
        avail = sorted(res_map.keys())
        chosen = int(resolution) if int(resolution) in res_map else min(
            avail, key=lambda r: abs(r - int(resolution))
        )
        uri = f"{cooler_path}::{res_map[chosen]}"
    else:
        uri = cooler_path
        chosen = int(cooler.Cooler(uri).binsize)

    clr = cooler.Cooler(uri)
    chosen_chrom = _pick_display_chrom([str(c) for c in clr.chromnames], requested=chrom)
    if not chosen_chrom:
        return None
    chrom = chosen_chrom
    has_w = "weight" in list(clr.bins().columns)
    region_start, region_end = _pick_mappable_window(clr, chrom, chosen, max_bins, has_w)
    if region_end <= region_start:
        return None

    used_balanced = False
    arr = None
    pos = np.array([], dtype=float)
    if has_w:
        mat = clr.matrix(balance=True).fetch((chrom, region_start, region_end))
        if mat.size > 0:
            arr = np.asarray(mat, dtype=float)
            pos = arr[np.isfinite(arr) & (arr > 0)]
            if pos.size > 0:
                used_balanced = True

    if not used_balanced:
        mat = clr.matrix(balance=False).fetch((chrom, region_start, region_end))
        if mat.size == 0:
            return None
        arr = np.asarray(mat, dtype=float)
        pos = arr[np.isfinite(arr) & (arr > 0)]
        if pos.size == 0:
            return None

    # Preserve masked/unmappable NaN bins as NaN (rendered in neutral gray) rather
    # than painting them as low-contact signal; map finite zeros to floor.
    floor = max(float(np.percentile(pos, 5)) * 0.5, 1e-8)
    log_mat = np.full_like(arr, np.nan, dtype=float)
    fin_mask = np.isfinite(arr)
    log_mat[fin_mask] = np.log10(np.where(arr[fin_mask] > 0, arr[fin_mask], floor))

    # Contrast clipping on off-diagonal contacts so extreme diagonal spikes do not wash out TADs.
    n_side = arr.shape[0]
    if n_side >= 3:
        off_mask = (~np.eye(n_side, dtype=bool)) & fin_mask & (arr > 0)
        contrast_vals = log_mat[off_mask] if np.any(off_mask) else np.log10(pos)
    else:
        contrast_vals = np.log10(pos)
    vmin, vmax = np.percentile(contrast_vals, [2, 98])
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        vmin, vmax = float(np.nanmin(log_mat)), float(np.nanmax(log_mat))
    if vmax <= vmin:
        vmax = vmin + 1.0

    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad("#e8e7e3")

    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    extent_mb = [
        region_start / 1e6,
        region_end / 1e6,
        region_end / 1e6,
        region_start / 1e6,
    ]
    im = ax.imshow(
        log_mat,
        cmap=cmap,
        vmin=float(vmin),
        vmax=float(vmax),
        extent=extent_mb,
        interpolation="nearest",
    )
    if tads_bed and Path(tads_bed).is_file() and Path(tads_bed).stat().st_size > 0:
        tdf = pd.read_csv(tads_bed, sep="\t", header=None)
        if tdf.shape[1] >= 3:
            t_start = pd.to_numeric(tdf[1], errors="coerce")
            t_end = pd.to_numeric(tdf[2], errors="coerce")
            sub_t = tdf[
                (tdf[0].astype(str) == str(chrom))
                & (t_start >= region_start)
                & (t_end <= region_end)
            ]
            for _, row in sub_t.iterrows():
                s_mb = float(row[1]) / 1e6
                e_mb = float(row[2]) / 1e6
                ax.plot([s_mb, e_mb, e_mb, s_mb, s_mb], [s_mb, s_mb, e_mb, e_mb, s_mb],
                        lw=0.8, color="#111111", alpha=0.75)

    res_label = _format_resolution(chosen)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    if used_balanced:
        cbar.set_label("log10 balanced contact frequency", fontsize=8)
        title = f"ICE-Balanced Contact Map ({chrom}, {res_label})"
    else:
        cbar.set_label("log10 raw contact count (unbalanced)", fontsize=8)
        title = f"Raw Contact Map — Unbalanced Fallback ({chrom}, {res_label})"
    ax.set_xlabel(f"{chrom} (Mb)", fontsize=9, color=_TEXT)
    ax.set_ylabel(f"{chrom} (Mb)", fontsize=9, color=_TEXT)
    ax.set_title(title, fontsize=9.5, color=_TEXT)
    fig.tight_layout()
    paths = _save_dual(fig, out_dir / "hic_contact_map.png")
    plt.close(fig)
    return paths


def _render_differential_compartments(diff_csv: str | None, out_dir: Path) -> dict[str, str] | None:
    if not diff_csv or not Path(diff_csv).is_file() or Path(diff_csv).stat().st_size == 0:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    df = pd.read_csv(diff_csv)
    if df.empty or "delta_E1" not in df.columns or "padj" not in df.columns:
        return None
    dx_all = pd.to_numeric(df["delta_E1"], errors="coerce").to_numpy(dtype=float)
    p_all = pd.to_numeric(df["padj"], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(dx_all) & np.isfinite(p_all) & (p_all >= 0.0) & (p_all <= 1.0)
    if not np.any(valid):
        return None

    sub = df.loc[valid].copy()
    dx = dx_all[valid]
    padj = p_all[valid]
    cls = sub.get("switch_class", pd.Series(["unchanged"] * len(sub), index=sub.index)).astype(str).to_numpy()

    # Data-adaptive ceiling for padj == 0.0 (never clip to 1e-300 -> y=300 spike).
    pos_p = padj[padj > 0.0]
    if pos_p.size > 0:
        max_fin_y = float(-np.log10(np.min(pos_p)))
        cap_y = max(max_fin_y * 1.15, max_fin_y + 1.0, 5.0)
    else:
        cap_y = 10.0
    zero_mask = padj == 0.0
    neglog = np.where(zero_mask, cap_y, -np.log10(np.where(zero_mask, 1.0, padj)))

    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    bg = cls == "unchanged"
    b2a = cls == "B_to_A"
    a2b = cls == "A_to_B"
    qsh = cls == "quantitative_shift"
    usf = cls == "unphased_sign_flip"

    if np.any(bg):
        ax.scatter(dx[bg], neglog[bg], s=12, color="#b8b7b2", alpha=0.55, label="Unchanged")
    if np.any(qsh):
        ax.scatter(dx[qsh], neglog[qsh], s=16, color="#eda100", alpha=0.8, label="Quantitative E1 shift")
    if np.any(usf):
        ax.scatter(dx[usf], neglog[usf], s=22, color=_COLOR_UNPHASED, alpha=0.9, label="Sign flip (unphased)")
    if np.any(b2a):
        ax.scatter(dx[b2a], neglog[b2a], s=22, color=_COLOR_A, alpha=0.9, label="B -> A switch")
    if np.any(a2b):
        ax.scatter(dx[a2b], neglog[a2b], s=22, color=_COLOR_B, alpha=0.9, label="A -> B switch")

    if np.any(zero_mask):
        ax.scatter(
            dx[zero_mask], neglog[zero_mask], s=28, facecolors="none",
            edgecolors="#111111", marker="^", linewidths=0.7,
            label=f"padj = 0 capped (y={cap_y:.1f})",
        )

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
    figures_svg: dict[str, str] = {}
    figure_bundles: dict[str, dict[str, str]] = {}
    skipped: dict[str, str] = {}

    def _record(key: str, bundle: dict[str, str] | None, empty_reason: str) -> None:
        if bundle:
            figures[key] = bundle["png"]
            if bundle.get("svg"):
                figures_svg[key] = bundle["svg"]
            figure_bundles[key] = bundle
        else:
            skipped[key] = empty_reason

    # 1. P(s) decay curve
    ps_inputs: list[tuple[str, str]] = []
    for rec in params.get("per_file") or []:
        if isinstance(rec, dict) and rec.get("ps_curve_tsv"):
            ps_inputs.append((str(rec.get("sample_id") or "sample"), str(rec["ps_curve_tsv"])))
    if not ps_inputs and params.get("ps_curve_tsv"):
        ps_inputs.append(("cohort", str(params["ps_curve_tsv"])))
    try:
        _record("hic_contact_decay_ps", _render_ps_decay(ps_inputs, out_dir), "no_ps_curve_tsv")
    except Exception as exc:
        skipped["hic_contact_decay_ps"] = f"render_error: {exc}"

    # 2. Compartment + Insulation track
    pc1_validated = params.get("pc1_validated")
    display_chrom = params.get("chrom") or params.get("chromosomes")
    try:
        _record(
            "hic_compartment_insulation_track",
            _render_compartment_insulation_track(
                params.get("compartments_tsv"),
                params.get("insulation_tsv"),
                out_dir,
                pc1_validated=pc1_validated if isinstance(pc1_validated, bool) else None,
                chrom=display_chrom,
            ),
            "no_compartment_or_insulation_tsv",
        )
    except Exception as exc:
        skipped["hic_compartment_insulation_track"] = f"render_error: {exc}"

    # 3. Contact map with TAD overlay
    try:
        _record(
            "hic_contact_map",
            _render_contact_map(
                params.get("cooler_path"),
                int(params.get("resolution") or 40_000),
                params.get("tads_bed"),
                out_dir,
                chrom=display_chrom,
            ),
            "no_cooler_contact_matrix",
        )
    except Exception as exc:
        skipped["hic_contact_map"] = f"render_error: {exc}"

    # 4. Differential compartments volcano
    try:
        _record(
            "hic_differential_compartments",
            _render_differential_compartments(
                params.get("differential_compartments_csv"),
                out_dir,
            ),
            "no_differential_compartments_csv",
        )
    except Exception as exc:
        skipped["hic_differential_compartments"] = f"render_error: {exc}"

    return {
        "status": "success",
        "ran": bool(figures),
        "data_type": "HiC",
        "analysis": "hic_figure",
        "validation_level": "production",
        "figures": figures,
        "figures_svg": figures_svg,
        "figure_bundles": figure_bundles,
        "skipped": skipped,
        "output_dir": str(out_dir),
    }


if __name__ == "__main__":
    run_script(hic_figure)
