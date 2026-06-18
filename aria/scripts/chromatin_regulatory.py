"""
ARIA scATAC regulatory layers (P2 primary-lane track)
-----------------------------------------------------
Optional regulatory analyses on a clustered scATAC peak matrix.

This script is deliberately input-gated. It computes only layers whose required
evidence is present and records a structured skip for the rest:

  - chromVAR-style motif activity from an explicit motif->peak map;
  - gene activity scores from peak coordinates + a local GTF;
  - peak-to-gene correlations from paired RNA+ATAC cells + a local GTF;
  - Tn5-bias-corrected footprinting only when fragments, motif sites, and a
    bias model are provided;
  - scRNA label transfer only as a hypothesis column, never an inferential group.

No biological identity is inferred from peak names. Coordinates are used only as
genomic coordinates for overlap/distance calculations (ADR-011 technical
exception). Missing prerequisites are honest `ran:false` results, never mocks.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def _skip(reason: str, **extra) -> dict[str, Any]:
    out = {"ran": False, "reason": reason}
    out.update(extra)
    return out


def _peak_coord(name: str) -> tuple[str, int, int] | None:
    text = str(name)
    try:
        if ":" in text and "-" in text:
            chrom, pos = text.split(":", 1)
            start, end = pos.replace("_", "-").split("-", 1)
            return chrom, int(start), int(end)
        parts = text.split("\t")
        if len(parts) >= 3:
            return parts[0], int(parts[1]), int(parts[2])
    except Exception:
        return None
    return None


def _load_motif_peak_map(path: str | None, peak_universe: set[str]) -> dict[str, list[str]]:
    """Read a motif->peak map CSV/TSV with columns motif_id and peak."""
    if not path or not Path(path).is_file():
        return {}
    sep = "\t" if str(path).endswith((".tsv", ".tsv.gz")) else ","
    out: dict[str, list[str]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        if not reader.fieldnames:
            return {}
        fields = {f.lower(): f for f in reader.fieldnames}
        motif_col = (
            fields.get("motif_id") or fields.get("motif") or fields.get("tf")
        )
        peak_col = fields.get("peak") or fields.get("region")
        if not motif_col or not peak_col:
            return {}
        for row in reader:
            motif = str(row.get(motif_col, "")).strip()
            peak = str(row.get(peak_col, "")).strip()
            if motif and peak in peak_universe:
                out.setdefault(motif, []).append(peak)
    return {m: sorted(set(peaks)) for m, peaks in out.items() if peaks}


def _parse_gtf_tss(gtf_path: str | None, wanted: set[str] | None = None) -> dict[str, tuple[str, int]]:
    """Parse gene TSS coordinates from a local GTF/GTF.GZ."""
    if not gtf_path or not Path(gtf_path).is_file():
        return {}
    import gzip

    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    coords: dict[str, tuple[str, int]] = {}
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            gene = _gtf_attr(attrs, "gene_name") or _gtf_attr(attrs, "gene_id")
            if not gene or (wanted is not None and gene not in wanted):
                continue
            start = int(parts[3])
            end = int(parts[4])
            tss = start if parts[6] != "-" else end
            coords[gene] = (parts[0], tss)
            if wanted is not None and len(coords) >= len(wanted):
                break
    return coords


def _gtf_attr(attrs: str, key: str) -> str | None:
    prefix = key + " "
    for chunk in attrs.split(";"):
        item = chunk.strip()
        if item.startswith(prefix):
            value = item[len(prefix):].strip()
            return value.strip('"') or None
    return None


def _genes_to_peak_indices(
    peak_names,
    gene_tss: dict[str, tuple[str, int]],
    *,
    window_bp: int,
) -> dict[str, list[int]]:
    peaks_by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for idx, peak in enumerate(peak_names):
        coord = _peak_coord(str(peak))
        if coord is None:
            continue
        chrom, start, end = coord
        peaks_by_chrom.setdefault(chrom, []).append((start, end, idx))

    out: dict[str, list[int]] = {}
    for gene, (chrom, tss) in gene_tss.items():
        hits = []
        lo, hi = tss - window_bp, tss + window_bp
        for start, end, idx in peaks_by_chrom.get(chrom, []):
            if end >= lo and start <= hi:
                hits.append(idx)
        if hits:
            out[gene] = hits
    return out


def _safe_mean(values) -> float | None:
    vals = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    if not vals:
        return None
    return round(sum(vals) / len(vals), 6)


def _corr(x, y) -> float | None:
    import numpy as np

    a = np.asarray(x, dtype=float).ravel()
    b = np.asarray(y, dtype=float).ravel()
    if a.shape[0] != b.shape[0] or a.shape[0] < 3:
        return None
    if float(np.std(a)) == 0.0 or float(np.std(b)) == 0.0:
        return None
    r = float(np.corrcoef(a, b)[0, 1])
    return round(r, 6) if math.isfinite(r) else None


def _vectorized_corr(x, cols):
    """Pearson correlation of vector ``x`` (n,) against each column of ``cols``
    (n, k), in one pass. Returns a length-k array, NaN where a column has zero
    variance or n < 3 — the same undefined cases :func:`_corr` returns None for.
    Used to vectorize the peak-to-gene inner loop (P4.2: O(genes*peaks) sparse
    per-pair getcol -> one masked matmul per gene)."""
    import numpy as np

    x = np.asarray(x, dtype=float).ravel()
    n = x.shape[0]
    k = cols.shape[1] if getattr(cols, "ndim", 1) == 2 else 1
    if n < 3:
        return np.full(k, np.nan)
    xc = x - x.mean()
    xs = float(np.sqrt((xc * xc).sum()))
    C = np.asarray(cols, dtype=float)
    if C.ndim == 1:
        C = C.reshape(-1, 1)
    Cc = C - C.mean(axis=0, keepdims=True)
    cs = np.sqrt((Cc * Cc).sum(axis=0))
    denom = xs * cs
    num = xc @ Cc
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(denom > 0, num / denom, np.nan)
    return r


def _param_int(params: dict, key: str, default: int) -> int:
    value = params.get(key)
    if value is None or value == "":
        return default
    return int(value)


def _param_float(params: dict, key: str, default: float) -> float:
    value = params.get(key)
    if value is None or value == "":
        return default
    return float(value)


def _col(mat, idx: int):
    import numpy as np
    from scipy import sparse

    if sparse.issparse(mat):
        return np.asarray(mat.getcol(idx).todense()).ravel()
    return np.asarray(mat[:, idx]).ravel()


def _sum_cols(mat, indices: list[int]):
    import numpy as np
    from scipy import sparse

    if sparse.issparse(mat):
        return np.asarray(mat[:, indices].sum(axis=1)).ravel()
    return np.asarray(mat[:, indices]).sum(axis=1)


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> str | None:
    if not rows:
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return str(path)


def _motif_activity(adata, params: dict, out_dir: Path) -> dict[str, Any]:
    import numpy as np

    motif_map_path = params.get("motif_peak_map") or params.get("motif_peak_map_csv")
    peak_names = [str(p) for p in adata.var_names]
    motif_map = _load_motif_peak_map(motif_map_path, set(peak_names))
    if not motif_map:
        return _skip(
            "motif_peak_map was not provided or did not overlap the peak universe",
            method="motif_peak_accessibility_deviation",
        )

    peak_index = {p: i for i, p in enumerate(peak_names)}
    rows = []
    for motif, peaks in motif_map.items():
        idx = [peak_index[p] for p in peaks if p in peak_index]
        if not idx:
            continue
        score = _sum_cols(adata.X, idx)
        mean_score = _safe_mean(score)
        depth = adata.obs.get("n_fragments")
        depth_corr = _corr(score, depth) if depth is not None else None
        rows.append({
            "motif_id": motif,
            "n_peaks": len(idx),
            "mean_accessibility_score": mean_score,
            "depth_correlation": depth_corr,
        })
    rows.sort(key=lambda r: (-(r["mean_accessibility_score"] or 0), r["motif_id"]))
    csv_path = _write_rows(out_dir / "chromatin_motif_activity.csv", rows)
    return {
        "ran": True,
        "method": "motif_peak_accessibility_deviation",
        "n_motifs": len(rows),
        "output_csv": csv_path,
        "top_motifs": rows[:20],
        "caveat": (
            "Motif activity is an accessibility deviation over motif-annotated "
            "peaks; it is not TF expression or proof of regulator activity."
        ),
    }


def _gene_scores(adata, params: dict, out_dir: Path) -> dict[str, Any]:
    window_bp = _param_int(params, "gene_score_window_bp", 100_000)
    gtf = params.get("gtf_path") or params.get("gene_annotation_gtf")
    gene_tss = _parse_gtf_tss(gtf)
    if not gene_tss:
        return _skip(
            "local GTF gene annotation was not provided; gene scores require real TSS coordinates",
            method="peak_window_gene_activity",
        )
    mapping = _genes_to_peak_indices(adata.var_names, gene_tss, window_bp=window_bp)
    if not mapping:
        return _skip(
            "no peaks overlap gene windows in the supplied GTF",
            method="peak_window_gene_activity",
        )
    rows = []
    for gene, idx in mapping.items():
        score = _sum_cols(adata.X, idx)
        rows.append({
            "gene": gene,
            "n_linked_peaks": len(idx),
            "mean_gene_activity_score": _safe_mean(score),
        })
    rows.sort(key=lambda r: (-(r["mean_gene_activity_score"] or 0), r["gene"]))
    csv_path = _write_rows(out_dir / "chromatin_gene_scores.csv", rows)
    return {
        "ran": True,
        "method": "peak_window_gene_activity",
        "window_bp": window_bp,
        "n_genes_scored": len(rows),
        "output_csv": csv_path,
        "top_genes": rows[:20],
        "caveat": (
            "Gene activity scores summarize nearby accessibility; they are not "
            "RNA expression measurements."
        ),
    }


def _peak_to_gene(adata, params: dict, out_dir: Path) -> dict[str, Any]:
    rna_path = params.get("rna_data_path") or params.get("scrna_data_path")
    if not rna_path or not Path(rna_path).is_file():
        return _skip(
            "paired scRNA h5ad was not provided",
            method="paired_cell_peak_gene_correlation",
        )
    gtf = params.get("gtf_path") or params.get("gene_annotation_gtf")
    min_shared = _param_int(params, "min_shared_cells", 100)
    min_corr = _param_float(params, "min_peak_gene_corr", 0.3)
    distance_bp = _param_int(params, "peak_gene_distance_bp", 500_000)

    import numpy as np
    from aria.utils.safe_h5ad import read_h5ad

    rna = read_h5ad(str(rna_path))
    common = [c for c in adata.obs_names if c in set(rna.obs_names)]
    if len(common) < min_shared:
        return _skip(
            f"only {len(common)} shared cells; peak-to-gene requires >= {min_shared}",
            method="paired_cell_peak_gene_correlation",
        )
    atac_aligned = adata[common, :]
    rna_aligned = rna[common, :]
    genes = {str(g) for g in rna_aligned.var_names}
    gene_tss = _parse_gtf_tss(gtf, wanted=genes)
    if not gene_tss:
        return _skip(
            "local GTF gene annotation was not provided or did not overlap RNA genes",
            method="paired_cell_peak_gene_correlation",
        )

    peak_names = [str(p) for p in atac_aligned.var_names]
    peak_coords = [_peak_coord(p) for p in peak_names]
    gene_index = {str(g): i for i, g in enumerate(rna_aligned.var_names)}

    # Spatial index: candidate peaks per chromosome (skip unparseable coords), so a
    # gene only ever scans its own chromosome instead of the whole peak universe.
    peaks_by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for pi, coord in enumerate(peak_coords):
        if coord is None:
            continue
        peaks_by_chrom.setdefault(coord[0], []).append((coord[1], coord[2], pi))

    rows = []
    tested = 0
    atac_X = atac_aligned.X
    for gene, (chrom, tss) in gene_tss.items():
        gi = gene_index.get(gene)
        if gi is None:
            continue
        expr = _col(rna_aligned.X, gi)
        if float(np.mean(expr)) <= 0.0:
            continue
        local = peaks_by_chrom.get(chrom)
        if not local:
            continue
        # In-window peaks: nearest-endpoint distance to the TSS <= distance_bp
        # (identical predicate to the per-pair version it replaces).
        cand_idx: list[int] = []
        cand_dist: list[int] = []
        for start, end, pi in local:
            dist = min(abs(start - tss), abs(end - tss))
            if dist <= distance_bp:
                cand_idx.append(pi)
                cand_dist.append(int(dist))
        if not cand_idx:
            continue
        tested += len(cand_idx)
        sub = atac_X[:, cand_idx]
        sub = np.asarray(sub.todense()) if hasattr(sub, "todense") else np.asarray(sub)
        rvals = _vectorized_corr(expr, sub)
        for j, pi in enumerate(cand_idx):
            rv = rvals[j]
            if not math.isfinite(float(rv)):  # NaN = undefined (zero variance / n<3)
                continue
            r = round(float(rv), 6)
            if abs(r) >= min_corr:
                rows.append({
                    "gene": gene,
                    "peak": peak_names[pi],
                    "correlation": r,
                    "distance_bp": cand_dist[j],
                    "direction": "positive" if r > 0 else "negative",
                })
    rows.sort(key=lambda r: (-abs(float(r["correlation"])), r["gene"], r["peak"]))
    csv_path = _write_rows(out_dir / "chromatin_peak_to_gene.csv", rows)
    return {
        "ran": True,
        # Peak-to-gene link recovery is the one regulatory layer promoted out of
        # scaffold: it is externally concordance-validated against a canonical
        # single-cell peak-gene linker (evidence lives in docs/ADR/benchmark
        # artifact, never hardcoded here per ADR-011). The other regulatory layers
        # (motif activity, gene scores, footprinting, label transfer) stay
        # scaffold/exploratory, so the script-level umbrella level is unchanged.
        "validation_level": "beta",
        "method": "paired_cell_peak_gene_correlation",
        "n_tested": tested,
        "n_links": len(rows),
        "n_positive": sum(1 for r in rows if r["direction"] == "positive"),
        "n_negative": sum(1 for r in rows if r["direction"] == "negative"),
        "output_csv": csv_path,
        "top_links": rows[:20],
        "caveat": (
            "Peak-to-gene links are correlations across paired cells; sign and "
            "distance are descriptive and do not establish regulatory mechanism."
        ),
    }


def _footprinting(params: dict) -> dict[str, Any]:
    fragments = params.get("fragments_file") or params.get("fragments_path")
    motif_sites = params.get("motif_sites_bed")
    bias_model = params.get("tn5_bias_model")
    if not fragments:
        return _skip("fragments_file was not provided", method="tn5_bias_corrected_footprint")
    if not Path(str(fragments)).is_file():
        return _skip(f"fragments_file not found: {fragments}", method="tn5_bias_corrected_footprint")
    if not motif_sites or not Path(str(motif_sites)).is_file():
        return _skip("motif_sites_bed was not provided", method="tn5_bias_corrected_footprint")
    if not bias_model or not Path(str(bias_model)).is_file():
        return _skip(
            "tn5_bias_model was not provided; ARIA does not report uncorrected footprints",
            method="tn5_bias_corrected_footprint",
        )
    return _skip(
        "Tn5-corrected footprinting runs via the dedicated TOBIAS backend "
        "(aria/scripts/chromatin_footprint_tobias.py in aria-tobias-env: ATACorrect -> "
        "ScoreBigwig -> BINDetect), dispatched separately like the SnapATAC2/Signac "
        "drivers; it is not run inline from the chromatin stack",
        method="tn5_bias_corrected_footprint",
    )


def _label_transfer(adata, params: dict, out_dir: Path) -> dict[str, Any]:
    rna_path = params.get("rna_data_path") or params.get("scrna_data_path")
    label_col = params.get("rna_label_col") or params.get("scrna_label_col")
    if not rna_path or not label_col:
        return _skip(
            "paired scRNA h5ad and rna_label_col are required",
            method="barcode_matched_label_hypothesis",
            trusted_for_inference=False,
        )
    if not Path(str(rna_path)).is_file():
        return _skip(
            f"paired scRNA h5ad not found: {rna_path}",
            method="barcode_matched_label_hypothesis",
            trusted_for_inference=False,
        )
    from aria.utils.safe_h5ad import read_h5ad

    rna = read_h5ad(str(rna_path))
    if label_col not in rna.obs.columns:
        return _skip(
            f"rna_label_col '{label_col}' is not present in scRNA obs",
            method="barcode_matched_label_hypothesis",
            trusted_for_inference=False,
        )
    common = [c for c in adata.obs_names if c in set(rna.obs_names)]
    if not common:
        return _skip(
            "no shared barcodes between scATAC and scRNA",
            method="barcode_matched_label_hypothesis",
            trusted_for_inference=False,
        )
    labels = rna.obs.loc[common, label_col].astype(str)
    rows = []
    for label, count in labels.value_counts().items():
        rows.append({"label": str(label), "n_cells": int(count)})
    csv_path = _write_rows(out_dir / "chromatin_label_transfer_hypotheses.csv", rows)
    return {
        "ran": True,
        "method": "barcode_matched_label_hypothesis",
        "n_labeled_cells": len(common),
        "n_labels": len(rows),
        "label_col": label_col,
        "output_csv": csv_path,
        "trusted_for_inference": False,
        "caveat": (
            "Transferred labels are report hypotheses only and must not be used "
            "as an inferential groupby."
        ),
    }


def chromatin_regulatory(params: dict) -> dict[str, Any]:
    from aria.utils.safe_h5ad import read_h5ad

    data_path = params.get("data_path")
    output_dir = Path(params.get("output_dir") or ".")
    if not data_path or not Path(str(data_path)).is_file():
        return {
            "status": "error",
            "error_type": "FileNotFound",
            "details": f"data_path does not exist: {data_path}",
        }
    try:
        adata = read_h5ad(str(data_path))
    except Exception as e:
        return {
            "status": "error",
            "error_type": "UnreadableInput",
            "details": str(e)[:400],
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "success",
        "validation_level": "alpha",
        "ran": True,
        "motif_activity": _motif_activity(adata, params, output_dir),
        "gene_scores": _gene_scores(adata, params, output_dir),
        "footprinting": _footprinting(params),
        "peak_to_gene": _peak_to_gene(adata, params, output_dir),
        "label_transfer": _label_transfer(adata, params, output_dir),
        "warnings": [],
    }
    for key in ("motif_activity", "gene_scores", "footprinting", "peak_to_gene", "label_transfer"):
        sub = result[key]
        if not sub.get("ran"):
            result["warnings"].append(f"{key} skipped: {sub.get('reason')}")
    return result


if __name__ == "__main__":
    run_script(chromatin_regulatory)
