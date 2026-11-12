"""Bulk ATAC publication figures (ADR-077): sample PCA + ENCODE QC summary.

validation_level: production (stack="rna", headless matplotlib).

Renders two sample-level panels for the bulk ATAC lane, each saved as a raster
``.png`` (inlined in the HTML report) AND a vector ``.svg`` (manuscript use):

* ``bulk_atac_pca``  — PCA (SVD on the centred ``log2(CPM + 1)`` peak matrix):
  PC1 vs PC2 with percent variance explained on the axes, every sample labelled,
  coloured (and shape-encoded) by condition when sample metadata is supplied.
  Raw integer counts are normalised to counts-per-million first; an already-CPM
  matrix (``*cpm*`` file name or non-integer values) is used as-is.
* ``bulk_atac_qc_summary`` — per-sample (or cohort) bars for the measured QC
  metrics (FRiP, TSS enrichment, duplicate rate) with the ENCODE reference
  thresholds drawn as dashed lines (FRiP 0.20, TSSe 5.0).

Zero fabrication: a panel is produced only when its data exists and is not
degenerate (PCA needs >= 2 samples and >= 2 non-constant peaks; the QC summary
needs at least one measured metric). Absent panels are simply omitted and the
reason is listed under ``skipped``.

Params:
  output_dir       str  (required)
  counts_path      str  peak_id-first TSV/CSV (raw counts or CPM)
  sample_metadata  dict {sample: {condition_col: level}} or a TSV/CSV path
  condition_col    str  default "condition"
  qc               dict chromatin_qc output (frip / tss_enrichment / dup_rate,
                        optionally *_per_sample, mito_frac)
  peaks            dict chromatin_peaks output (frip / frip_per_sample, n_peaks,
                        n_consensus_peaks)
"""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from aria.scripts._base import run_script

# ENCODE ATAC-seq reference thresholds (dashed guides, not pass/fail verdicts).
ENCODE_FRIP_MIN = 0.20
ENCODE_TSSE_MIN = 5.0

# Validated categorical palette, assigned in fixed order to the SORTED condition
# levels so a colour follows the condition, never its rank in a given run.
_CATEGORICAL = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                "#e87ba4", "#008300", "#6250d6", "#e34948")
_MARKERS = ("o", "s", "^", "D", "v", "P", "X", "*")
_SINGLE_SERIES = "#2a78d6"
_TEXT_PRIMARY = "#0b0b0b"
_TEXT_SECONDARY = "#52514e"
_GUIDE = "#52514e"
_SURFACE = "#fcfcfb"


def _save_dual(fig, png_path: Path) -> str:
    """PNG (report) + SVG (manuscript) side by side; returns the PNG path."""
    png_path = Path(png_path)
    fig.savefig(png_path, dpi=200, bbox_inches="tight", facecolor=_SURFACE)
    fig.savefig(png_path.with_suffix(".svg"), bbox_inches="tight",
                facecolor=_SURFACE)
    return str(png_path)


def _style_axis(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#b5b4ae")
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=_TEXT_SECONDARY, labelsize=7, length=3, width=0.6)
    ax.set_facecolor(_SURFACE)


# ── inputs ────────────────────────────────────────────────────────────────

def _read_table(path: str):
    import pandas as pd

    name = Path(path).name.lower()
    sep = "," if name.endswith((".csv", ".csv.gz")) else "\t"
    return pd.read_csv(path, sep=sep)


def _load_counts(path: Optional[str]):
    """Return (peaks x samples float DataFrame, is_cpm) or (None, reason)."""
    if not path or not Path(str(path)).is_file():
        return None, "no counts matrix supplied"
    import numpy as np
    import pandas as pd

    df = _read_table(str(path))
    if df.shape[1] < 2:
        return None, "counts matrix has no sample columns"
    mat = df.set_index(df.columns[0])
    mat = mat.apply(lambda col: pd.to_numeric(col, errors="coerce"))
    mat = mat.dropna(axis=0, how="any")
    if mat.empty:
        return None, "counts matrix has no complete numeric rows"
    values = mat.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        return None, "counts matrix contains negative or non-finite values"
    integer_like = bool(np.allclose(values, np.round(values)))
    is_cpm = ("cpm" in Path(str(path)).name.lower()) or not integer_like
    return mat.astype(float), is_cpm


def _load_conditions(sample_metadata, condition_col: str,
                     samples: list[str]) -> dict[str, str]:
    """Map sample -> condition level; empty when not resolvable (no guessing)."""
    if not sample_metadata:
        return {}
    records: dict[str, Any] = {}
    if isinstance(sample_metadata, dict):
        records = sample_metadata
    elif isinstance(sample_metadata, (str, Path)) and Path(
            str(sample_metadata)).is_file():
        meta = _read_table(str(sample_metadata))
        id_col = "sample_id" if "sample_id" in meta.columns else meta.columns[0]
        if condition_col not in meta.columns:
            return {}
        records = {str(r[id_col]): {condition_col: r[condition_col]}
                   for _, r in meta.iterrows()}
    out: dict[str, str] = {}
    for sample in samples:
        meta = records.get(sample)
        if isinstance(meta, dict) and meta.get(condition_col) not in (None, ""):
            value = meta.get(condition_col)
            if isinstance(value, float) and math.isnan(value):
                continue
            out[sample] = str(value)
    return out


# ── panels ────────────────────────────────────────────────────────────────

def _pca(mat, is_cpm: bool):
    """SVD-PCA on centred log2(CPM+1). Returns (scores, pct_var, samples) or
    (None, reason)."""
    import numpy as np

    values = mat.to_numpy(dtype=float)
    if not is_cpm:
        lib = values.sum(axis=0)
        if (lib <= 0).any():
            return None, "a sample has zero total counts"
        values = values / lib * 1e6
    logv = np.log2(values + 1.0)
    keep = logv.std(axis=1) > 0
    if int(keep.sum()) < 2:
        return None, "fewer than 2 non-constant peaks"
    x = logv[keep].T  # samples x peaks
    x = x - x.mean(axis=0, keepdims=True)
    _u, s, vt = np.linalg.svd(x, full_matrices=False)
    var = s ** 2
    if var.sum() <= 0 or len(s) < 2:
        return None, "no variance across samples"
    scores = x @ vt.T
    pct = 100.0 * var / var.sum()
    return (scores[:, :2], pct[:2], list(mat.columns)), None


def render_pca(mat, is_cpm: bool, conditions: dict[str, str],
               out_path: Path) -> tuple[Optional[str], Optional[str]]:
    if mat is None or mat.shape[1] < 2:
        return None, "fewer than 2 samples"
    res, reason = _pca(mat, is_cpm)
    if res is None:
        return None, reason
    scores, pct, samples = res
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.6, 4.6), dpi=200)
    fig.patch.set_facecolor(_SURFACE)
    _style_axis(ax)
    ax.axhline(0, color="#e2e1dc", lw=0.6, zorder=0)
    ax.axvline(0, color="#e2e1dc", lw=0.6, zorder=0)
    levels = sorted(set(conditions.values()))
    if levels:
        for i, level in enumerate(levels):
            idx = [j for j, s in enumerate(samples) if conditions.get(s) == level]
            ax.scatter(scores[idx, 0], scores[idx, 1], s=64,
                       c=_CATEGORICAL[i % len(_CATEGORICAL)],
                       marker=_MARKERS[i % len(_MARKERS)],
                       edgecolors=_SURFACE, linewidths=1.6, zorder=3,
                       label=level)
        unassigned = [j for j, s in enumerate(samples) if s not in conditions]
        if unassigned:
            ax.scatter(scores[unassigned, 0], scores[unassigned, 1], s=64,
                       c="#9a9992", edgecolors=_SURFACE, linewidths=1.6,
                       zorder=3, label="no condition")
        ax.legend(frameon=False, fontsize=7, labelcolor=_TEXT_SECONDARY,
                  loc="best", title=None)
    else:
        ax.scatter(scores[:, 0], scores[:, 1], s=64, c=_SINGLE_SERIES,
                   edgecolors=_SURFACE, linewidths=1.6, zorder=3)
    for j, sample in enumerate(samples):
        ax.annotate(str(sample), (scores[j, 0], scores[j, 1]),
                    xytext=(5, 4), textcoords="offset points", fontsize=6,
                    color=_TEXT_SECONDARY)
    ax.set_xlabel(f"PC1 ({pct[0]:.1f}% variance)", fontsize=8,
                  color=_TEXT_PRIMARY)
    ax.set_ylabel(f"PC2 ({pct[1]:.1f}% variance)", fontsize=8,
                  color=_TEXT_PRIMARY)
    ax.set_title("Bulk ATAC sample PCA (log2 CPM + 1)", fontsize=9,
                 fontweight="bold", color=_TEXT_PRIMARY, loc="left")
    ax.margins(0.15)
    out = _save_dual(fig, out_path)
    plt.close(fig)
    return out, None


def _metric_series(qc: dict, peaks: dict, per_sample_keys: tuple[str, ...],
                   cohort_keys: tuple[str, ...]) -> dict[str, float]:
    """Per-sample values when measured, else a single cohort value, else {}."""
    for source in (qc, peaks):
        for key in per_sample_keys:
            raw = (source or {}).get(key)
            if isinstance(raw, dict):
                vals = {str(k): float(v) for k, v in raw.items()
                        if isinstance(v, (int, float)) and math.isfinite(v)}
                if vals:
                    return vals
    for source in (qc, peaks):
        for key in cohort_keys:
            raw = (source or {}).get(key)
            if isinstance(raw, (int, float)) and math.isfinite(raw):
                return {"cohort": float(raw)}
    return {}


def render_qc_summary(qc: dict, peaks: dict,
                      out_path: Path) -> tuple[Optional[str], Optional[str]]:
    panels = []
    frip = _metric_series(qc, peaks, ("frip_per_sample",), ("frip",))
    if frip:
        panels.append(("FRiP", frip, ENCODE_FRIP_MIN,
                       f"ENCODE ≥ {ENCODE_FRIP_MIN:.2f}"))
    tsse = _metric_series(qc, peaks, ("tss_enrichment_per_sample",),
                          ("tss_enrichment",))
    if tsse:
        panels.append(("TSS enrichment", tsse, ENCODE_TSSE_MIN,
                       f"ENCODE ≥ {ENCODE_TSSE_MIN:.0f}"))
    dup = _metric_series(qc, peaks, ("dup_rate_per_sample",), ("dup_rate",))
    if dup:
        panels.append(("Duplicate rate", dup, None, None))
    if not panels:
        return None, "no measured FRiP / TSS enrichment / duplicate rate"

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_samples = max(len(series) for _, series, _, _ in panels)
    width = max(2.6, 0.32 * n_samples + 1.4)
    fig, axes = plt.subplots(1, len(panels),
                             figsize=(width * len(panels), 3.4), dpi=200,
                             squeeze=False)
    fig.patch.set_facecolor(_SURFACE)
    for ax, (label, series, ref, ref_label) in zip(axes[0], panels):
        _style_axis(ax)
        names = list(series)
        vals = [series[n] for n in names]
        ax.bar(range(len(names)), vals, width=0.7, color=_SINGLE_SERIES,
               edgecolor=_SURFACE, linewidth=1.0, zorder=2)
        if ref is not None:
            ax.axhline(ref, color=_GUIDE, lw=1.0, ls="--", zorder=3,
                       label=ref_label)
            ax.legend(frameon=False, fontsize=6, labelcolor=_TEXT_SECONDARY,
                      loc="upper right")
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=45 if len(names) > 1 else 0,
                           ha="right" if len(names) > 1 else "center",
                           fontsize=6, color=_TEXT_SECONDARY)
        top = max(vals + ([ref] if ref is not None else [])) or 1.0
        ax.set_ylim(0, top * 1.18)
        ax.set_title(label, fontsize=8, fontweight="bold",
                     color=_TEXT_PRIMARY, loc="left")
        ax.grid(axis="y", color="#ecebe6", lw=0.5, zorder=0)
    fig.suptitle("Bulk ATAC library QC (ENCODE reference guides dashed)",
                 fontsize=9, fontweight="bold", color=_TEXT_PRIMARY,
                 x=0.01, ha="left")
    fig.tight_layout()
    out = _save_dual(fig, out_path)
    plt.close(fig)
    return out, None


# ── entry point ─────────────────────────────────────────────────────────────

def render_bulk_figures(params: dict) -> dict:
    """In-process entry (also used by ``_narrative_chromatin`` when no env
    manager is available). Returns the standard script result dict."""
    output_dir = params.get("output_dir")
    if not output_dir:
        return {"status": "error", "error_type": "MissingParameter",
                "details": "output_dir is required"}
    output_dir = Path(str(output_dir))
    output_dir.mkdir(parents=True, exist_ok=True)
    condition_col = str(params.get("condition_col") or "condition")
    qc = params.get("qc") if isinstance(params.get("qc"), dict) else {}
    peaks = params.get("peaks") if isinstance(params.get("peaks"), dict) else {}
    counts_path = (params.get("counts_path")
                   or params.get("counts_cpm_path")
                   or params.get("replicate_counts_path")
                   or params.get("counts_matrix_path"))

    figures: dict[str, str] = {}
    skipped: dict[str, str] = {}
    warnings: list[str] = []

    mat, is_cpm_or_reason = _load_counts(counts_path)
    if mat is None:
        skipped["bulk_atac_pca"] = str(is_cpm_or_reason)
    else:
        conditions = _load_conditions(params.get("sample_metadata"),
                                      condition_col, list(mat.columns))
        try:
            path, reason = render_pca(mat, bool(is_cpm_or_reason), conditions,
                                      output_dir / "bulk_atac_pca.png")
        except Exception as exc:  # a figure failure never aborts the report
            path, reason = None, f"PCA rendering failed: {str(exc)[:160]}"
            warnings.append(reason)
        if path:
            figures["bulk_atac_pca"] = path
        else:
            skipped["bulk_atac_pca"] = str(reason)

    try:
        path, reason = render_qc_summary(qc, peaks,
                                         output_dir / "bulk_atac_qc_summary.png")
    except Exception as exc:
        path, reason = None, f"QC summary rendering failed: {str(exc)[:160]}"
        warnings.append(reason)
    if path:
        figures["bulk_atac_qc_summary"] = path
    else:
        skipped["bulk_atac_qc_summary"] = str(reason)

    return {
        "status": "success",
        "validation_level": "production",
        "data_type": "bulk_ATAC",
        "figures": figures,
        "skipped": skipped,
        "warnings": warnings,
    }


def chromatin_figure_bulk(params: dict) -> dict:
    return render_bulk_figures(params)


if __name__ == "__main__":
    run_script(chromatin_figure_bulk)
