"""Deterministic SEACR peak calling, Li et al. (2011) Gaussian-copula IDR,
broad megadomain calling, and normalized BigWig export for epigenomic assays.

validation_level: production (ADR-081).

Boundary invariant: this module is imported inside Conda scientific subprocesses
(``aria-chromatin-env``) and must NEVER import ``aria.agents.*``, ``litellm``,
or ``pydantic``.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np


def _open_alignment(path: str | Path, reference_fasta: str | None = None):
    import pysam

    mode = "rc" if str(path).lower().endswith(".cram") else "rb"
    kwargs = {"reference_filename": str(reference_fasta)} if reference_fasta else {}
    return pysam.AlignmentFile(str(path), mode, **kwargs)


def bam_to_fragment_bedgraph(
    bam_path: str | Path,
    output_bedgraph: str | Path,
    *,
    scale_factor: float = 1.0,
    max_fragment_bp: int = 1000,
    single_end_extsize: int = 200,
    min_mapq: int = 0,
    keep_duplicates: bool = False,
    mito_contigs: tuple[str, ...] = ("chrM", "MT", "M"),
    reference_fasta: str | None = None,
) -> dict[str, Any]:
    """Stream coordinate-sorted BAM/CRAM to an exact fragment coverage bedGraph.

    Follows the Meers et al. (2019) SEACR paired-end fragment definition
    (properly paired reads on the same contig with fragment length < 1000 bp)
    without requiring an intermediate name-sorted BAM, and supports single-end
    extension when unpaired alignments are supplied.
    """
    scale = float(scale_factor)
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError(f"scale_factor must be finite and positive, got {scale_factor!r}")

    mito = {c.lower() for c in mito_contigs}
    frags_by_chrom: dict[str, tuple[list[int], list[int]]] = {}
    chrom_lengths: dict[str, int] = {}
    n_fragments = 0
    n_paired = 0
    n_single = 0

    with _open_alignment(bam_path, reference_fasta) as bam:
        for ref, length in zip(bam.header.references, bam.header.lengths):
            chrom_lengths[str(ref)] = int(length)

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not keep_duplicates and read.is_duplicate:
                continue
            if int(read.mapping_quality) < int(min_mapq):
                continue
            chrom = read.reference_name
            if not chrom or chrom.lower() in mito or chrom.lower().startswith("chrm"):
                continue
            chrom_len = chrom_lengths.get(chrom, 0)

            if read.is_paired:
                if (
                    not read.is_proper_pair
                    or not read.is_read1
                    or read.mate_is_unmapped
                    or read.next_reference_name != chrom
                ):
                    continue
                tlen = abs(int(read.template_length))
                if tlen <= 0 or tlen >= int(max_fragment_bp):
                    continue
                start = min(int(read.reference_start), int(read.next_reference_start))
                end = start + tlen
                if chrom_len > 0:
                    start = max(0, start)
                    end = min(chrom_len, end)
                if end <= start:
                    continue
                n_paired += 1
            else:
                if read.reference_end is None:
                    continue
                ext = max(1, int(single_end_extsize))
                if read.is_reverse:
                    end = int(read.reference_end)
                    start = max(0, end - ext)
                else:
                    start = max(0, int(read.reference_start))
                    end = start + ext
                if chrom_len > 0:
                    end = min(chrom_len, end)
                if end <= start:
                    continue
                n_single += 1

            starts, ends = frags_by_chrom.setdefault(chrom, ([], []))
            starts.append(start)
            ends.append(end)
            n_fragments += 1

    out_path = Path(output_bedgraph)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    chrom_order = [c for c in chrom_lengths if c in frags_by_chrom]
    for extra in sorted(set(frags_by_chrom) - set(chrom_order)):
        chrom_order.append(extra)

    n_intervals = 0
    with open(out_path, "w", encoding="utf-8") as out:
        for chrom in chrom_order:
            starts_l, ends_l = frags_by_chrom[chrom]
            if not starts_l:
                continue
            s_arr = np.asarray(starts_l, dtype=np.int64)
            e_arr = np.asarray(ends_l, dtype=np.int64)
            coords = np.unique(np.concatenate([s_arr, e_arr]))
            if coords.size < 2:
                continue
            delta = np.zeros(coords.size, dtype=np.int64)
            np.add.at(delta, np.searchsorted(coords, s_arr), 1)
            np.add.at(delta, np.searchsorted(coords, e_arr), -1)
            cov = np.cumsum(delta)[:-1]
            nz = np.flatnonzero(cov > 0)
            for idx in nz:
                val = float(cov[idx]) * scale
                val_str = (
                    str(int(round(val)))
                    if abs(val - round(val)) < 1e-9
                    else f"{val:.6f}".rstrip("0").rstrip(".")
                )
                out.write(f"{chrom}\t{int(coords[idx])}\t{int(coords[idx + 1])}\t{val_str}\n")
                n_intervals += 1

    return {
        "bedgraph_path": str(out_path),
        "n_fragments": int(n_fragments),
        "n_paired_fragments": int(n_paired),
        "n_single_end_reads": int(n_single),
        "n_bedgraph_intervals": int(n_intervals),
        "scale_factor": float(scale),
        "chrom_lengths": chrom_lengths,
    }


def _integrate_bedgraph_signal_blocks(bedgraph_path: str | Path) -> list[dict[str, Any]]:
    """Collapse contiguous non-zero bedGraph intervals into SEACR signal blocks.

    In SEACR (Meers et al., 2019), contiguous non-zero coverage intervals
    ``(prev_end == start)`` on the same chromosome form one signal block whose
    ``total_signal`` is the integrated area under the curve
    ``sum((end - start) * score)`` and ``max_signal`` is ``max(score)``.
    """
    blocks: list[dict[str, Any]] = []
    cur_chrom: str | None = None
    cur_start = 0
    cur_end = 0
    cur_auc = 0.0
    cur_max = -1.0
    max_start = 0
    max_end = 0

    def _flush() -> None:
        if cur_chrom is not None and cur_end > cur_start and cur_auc > 0:
            blocks.append({
                "chrom": cur_chrom,
                "start": int(cur_start),
                "end": int(cur_end),
                "total_signal": float(cur_auc),
                "max_signal": float(cur_max),
                "max_region": f"{cur_chrom}:{int(max_start)}-{int(max_end)}",
            })

    with open(bedgraph_path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                chrom = parts[0]
                start = int(float(parts[1]))
                end = int(float(parts[2]))
                score = float(parts[3])
            except ValueError:
                continue
            if end <= start or score <= 0:
                continue

            area = (end - start) * score
            if cur_chrom == chrom and start == cur_end:
                cur_end = end
                cur_auc += area
                if score > cur_max + 1e-12:
                    cur_max = score
                    max_start = start
                    max_end = end
                elif abs(score - cur_max) <= 1e-12 and start == max_end:
                    max_end = end
            else:
                _flush()
                cur_chrom = chrom
                cur_start = start
                cur_end = end
                cur_auc = area
                cur_max = score
                max_start = start
                max_end = end

    _flush()
    return blocks


def _estimate_seacr_control_scale(
    target_blocks: list[dict[str, Any]],
    control_blocks: list[dict[str, Any]],
) -> float:
    """Estimate background normalization factor matching SEACR ``norm`` mode.

    Uses the median per-bp background density across lower-half signal blocks
    so enriched target peaks do not inflate the background ratio.
    """
    if not target_blocks or not control_blocks:
        return 1.0
    t_dens = np.array(
        [b["total_signal"] / max(1, b["end"] - b["start"]) for b in target_blocks],
        dtype=float,
    )
    c_dens = np.array(
        [b["total_signal"] / max(1, b["end"] - b["start"]) for b in control_blocks],
        dtype=float,
    )
    t_bg = float(np.median(t_dens[t_dens <= np.quantile(t_dens, 0.5)]))
    c_bg = float(np.median(c_dens[c_dens <= np.quantile(c_dens, 0.5)]))
    if c_bg <= 0 or not math.isfinite(t_bg / c_bg):
        return 1.0
    return float(t_bg / c_bg)


def _select_seacr_empirical_threshold(
    target_auc: np.ndarray,
    control_auc: np.ndarray,
    mode: str = "stringent",
) -> tuple[float, dict[str, Any]]:
    """Select empirical AUC cutoff from target vs IgG/control survival curves.

    Implements SEACR's empirical thresholding principle:
      - Evaluate candidate cutoffs across sorted combined block total signals.
      - Compute target retention ``S_T(t) = #{T >= t}/|T|``, control retention
        ``S_C(t) = #{C >= t}/|C|``, and target share among retained blocks
        ``P(t) = #{T >= t} / (#{T >= t} + #{C >= t})``.
      - ``t_global`` maximizes ``P(t) - S_C(t)`` (empirical precision minus
        background survival). ``stringent`` selects ``t_global``; ``relaxed``
        selects the knee threshold halfway in target-retention space between
        ``t_global`` and the median background block so broader domains are
        retained while guaranteeing ``threshold_relaxed <= threshold_stringent``
        (hence ``n_peaks_relaxed >= n_peaks_stringent``).
    """
    mode_norm = str(mode or "stringent").strip().lower()
    if mode_norm not in {"stringent", "relaxed"}:
        raise ValueError(f"SEACR mode must be 'stringent' or 'relaxed', got {mode!r}")

    t_sorted = np.sort(np.asarray(target_auc, dtype=float))
    c_sorted = np.sort(np.asarray(control_auc, dtype=float))
    if t_sorted.size == 0:
        return float("inf"), {"empirical_fdr": None, "precision": None}
    if c_sorted.size == 0:
        thr = float(np.quantile(t_sorted, 0.99))
        return thr, {"empirical_fdr": 0.0, "precision": 1.0}

    grid = np.unique(
        np.quantile(
            np.concatenate([t_sorted, c_sorted]),
            np.linspace(0.0, 0.999, min(2000, t_sorted.size + c_sorted.size)),
        )
    )
    n_t = float(t_sorted.size)
    n_c = float(c_sorted.size)
    t_cum_above = np.concatenate([[0.0], np.cumsum(t_sorted[::-1])])[::-1]
    c_cum_above = np.concatenate([[0.0], np.cumsum(c_sorted[::-1])])[::-1]
    t_tot = max(1e-12, float(t_cum_above[0]))
    c_tot = max(1e-12, float(c_cum_above[0]))

    t_idx = np.searchsorted(t_sorted, grid, side="left")
    c_idx = np.searchsorted(c_sorted, grid, side="left")
    t_above = (t_sorted.size - t_idx).astype(float)
    c_above = (c_sorted.size - c_idx).astype(float)
    t_sig_frac = t_cum_above[t_idx] / t_tot
    c_sig_frac = c_cum_above[c_idx] / c_tot
    sig_denom = np.maximum(1e-12, t_cum_above[t_idx] + c_cum_above[c_idx])
    sig_precision = t_cum_above[t_idx] / sig_denom
    block_precision = t_above / np.maximum(1.0, t_above + c_above)

    # Meers et al. (2019) SEACR signal-excess objective:
    # relaxed maximizes total target signal fraction minus control signal fraction;
    # stringent weights signal excess by empirical block/signal precision so IgG
    # tail noise is rejected without collapsing onto a tiny extreme-tail subset.
    excess = np.where(t_above >= 1.0, t_sig_frac - c_sig_frac, -np.inf)
    stringent_obj = np.where(
        t_above >= 1.0,
        excess * (0.5 * sig_precision + 0.5 * block_precision) * (1.0 - c_above / n_c),
        -np.inf,
    )

    rel_idx = int(np.argmax(excess))
    thr_relaxed = float(grid[rel_idx])
    # Stringent threshold must be >= relaxed threshold
    str_candidates = np.where(grid >= thr_relaxed)[0]
    if str_candidates.size > 0:
        best_idx = int(str_candidates[int(np.argmax(stringent_obj[str_candidates]))])
    else:
        best_idx = int(np.argmax(stringent_obj))
    thr_stringent = float(max(thr_relaxed, float(grid[best_idx])))

    if mode_norm == "stringent":
        thr = thr_stringent
    else:
        thr = thr_relaxed

    n_t_keep = int(t_sorted.size - np.searchsorted(t_sorted, thr, side="left"))
    n_c_keep = int(c_sorted.size - np.searchsorted(c_sorted, thr, side="left"))
    emp_fdr = float(n_c_keep / max(1, n_t_keep + n_c_keep))
    return thr, {
        "threshold_stringent": thr_stringent,
        "threshold_selected": thr,
        "empirical_fdr": round(emp_fdr, 6),
        "precision": round(float(n_t_keep / max(1, n_t_keep + n_c_keep)), 6),
        "n_target_blocks_total": int(t_sorted.size),
        "n_control_blocks_total": int(c_sorted.size),
        "n_control_blocks_above_threshold": n_c_keep,
    }


def call_seacr_peaks(
    target_inputs: list[str] | str,
    output_bed: str | Path,
    *,
    control_inputs: list[str] | str | None = None,
    mode: str = "stringent",
    norm: str = "norm",
    top_quantile: float = 0.01,
    target_scale_factors: list[float] | None = None,
    control_scale_factors: list[float] | None = None,
    already_spikein_scaled: bool = False,
    max_fragment_bp: int = 1000,
    single_end_extsize: int = 200,
    keep_duplicates: bool = True,
    reference_fasta: str | None = None,
    work_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Call sparse chromatin enrichment peaks with SEACR (Meers et al., 2019).

    Accepts either pre-computed bedGraph files (``.bedgraph`` / ``.bg``) or
    coordinate-sorted BAM/CRAM alignments. When ``already_spikein_scaled=True``
    (or ``target_scale_factors`` are provided), ``norm`` is forced to ``"non"``
    per SEACR best practices so heterologous spike-in global scaling is preserved.
    """
    targets = [target_inputs] if isinstance(target_inputs, (str, Path)) else list(target_inputs or [])
    controls = (
        []
        if control_inputs is None
        else ([control_inputs] if isinstance(control_inputs, (str, Path)) else list(control_inputs))
    )
    if not targets:
        return {"status": "error", "error_type": "NoInputFiles", "details": "SEACR requires target input(s)."}

    out_path = Path(output_bed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = Path(work_dir) if work_dir else out_path.parent / "_seacr_work"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    def _to_bedgraphs(paths: list[str], scales: list[float] | None, prefix: str) -> list[str]:
        bgs: list[str] = []
        for idx, p in enumerate(paths):
            p_str = str(p)
            scale = float(scales[idx]) if (scales and idx < len(scales)) else 1.0
            if not math.isfinite(scale) or scale <= 0.0:
                raise ValueError(f"Invalid SEACR scale factor for {p_str}: {scale!r}")
            if p_str.lower().endswith((".bedgraph", ".bg", ".bdg")):
                if abs(scale - 1.0) <= 1e-12:
                    bgs.append(p_str)
                else:
                    bg_scaled = tmp_dir / f"{prefix}_{idx}_{Path(p_str).stem}.scaled.bedgraph"
                    with open(p_str, encoding="utf-8") as src, open(bg_scaled, "w", encoding="utf-8") as dst:
                        for line in src:
                            if not line.strip() or line.startswith(("#", "track", "browser")):
                                continue
                            cols = line.rstrip("\n").split("\t")
                            if len(cols) >= 4:
                                try:
                                    val = float(cols[3]) * scale
                                    dst.write(f"{cols[0]}\t{cols[1]}\t{cols[2]}\t{val:.6f}\n")
                                except ValueError:
                                    continue
                    bgs.append(str(bg_scaled))
            else:
                bg_out = tmp_dir / f"{prefix}_{idx}_{Path(p_str).stem}.bedgraph"
                bam_to_fragment_bedgraph(
                    p_str,
                    bg_out,
                    scale_factor=scale,
                    max_fragment_bp=max_fragment_bp,
                    single_end_extsize=single_end_extsize,
                    keep_duplicates=keep_duplicates,
                    reference_fasta=reference_fasta,
                )
                bgs.append(str(bg_out))
        return bgs

    def _merge_bedgraphs_mean(bg_paths: list[str], merged_out: Path) -> str:
        if len(bg_paths) == 1:
            return bg_paths[0]
        by_chrom: dict[str, tuple[list[int], list[int], list[float]]] = {}
        for bg in bg_paths:
            with open(bg, encoding="utf-8") as fh:
                for line in fh:
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 4:
                        continue
                    try:
                        c, s, e, v = cols[0], int(float(cols[1])), int(float(cols[2])), float(cols[3])
                    except ValueError:
                        continue
                    if e > s and v > 0:
                        sl, el, vl = by_chrom.setdefault(c, ([], [], []))
                        sl.append(s)
                        el.append(e)
                        vl.append(v / len(bg_paths))
        with open(merged_out, "w", encoding="utf-8") as out:
            for chrom in sorted(by_chrom):
                s_arr = np.asarray(by_chrom[chrom][0], dtype=np.int64)
                e_arr = np.asarray(by_chrom[chrom][1], dtype=np.int64)
                v_arr = np.asarray(by_chrom[chrom][2], dtype=np.float64)
                coords = np.unique(np.concatenate([s_arr, e_arr]))
                if coords.size < 2:
                    continue
                delta = np.zeros(coords.size, dtype=np.float64)
                np.add.at(delta, np.searchsorted(coords, s_arr), v_arr)
                np.add.at(delta, np.searchsorted(coords, e_arr), -v_arr)
                cov = np.cumsum(delta)[:-1]
                for idx in np.flatnonzero(cov > 1e-9):
                    out.write(f"{chrom}\t{int(coords[idx])}\t{int(coords[idx + 1])}\t{cov[idx]:.6f}\n")
        return str(merged_out)

    used_spikein = bool(already_spikein_scaled or target_scale_factors or control_scale_factors)
    effective_norm = "non" if used_spikein else str(norm or "norm").strip().lower()

    target_bgs = _to_bedgraphs(targets, target_scale_factors, "target")
    merged_target_bg = _merge_bedgraphs_mean(target_bgs, tmp_dir / "pooled_target.bedgraph")
    target_blocks = _integrate_bedgraph_signal_blocks(merged_target_bg)
    if not target_blocks:
        out_path.write_text("", encoding="utf-8")
        return {
            "status": "success",
            "peak_caller": "seacr",
            "mode": mode,
            "norm": effective_norm,
            "n_peaks": 0,
            "peaks_path": str(out_path),
            "threshold": None,
            "warnings": ["Target alignment/bedGraph contained zero positive signal blocks."],
        }

    target_auc = np.array([b["total_signal"] for b in target_blocks], dtype=float)
    warnings_out: list[str] = []
    threshold_meta: dict[str, Any] = {}

    if controls:
        ctrl_bgs = _to_bedgraphs(controls, control_scale_factors, "control")
        merged_ctrl_bg = _merge_bedgraphs_mean(ctrl_bgs, tmp_dir / "pooled_control.bedgraph")
        ctrl_blocks = _integrate_bedgraph_signal_blocks(merged_ctrl_bg)
        ctrl_scale = (
            _estimate_seacr_control_scale(target_blocks, ctrl_blocks)
            if effective_norm == "norm"
            else 1.0
        )
        ctrl_auc = np.array([b["total_signal"] * ctrl_scale for b in ctrl_blocks], dtype=float)
        threshold, threshold_meta = _select_seacr_empirical_threshold(
            target_auc, ctrl_auc, mode=mode
        )
        threshold_meta["control_scale_factor"] = round(float(ctrl_scale), 6)
        threshold_meta["threshold_method"] = f"empirical_igg_{mode}"
    else:
        q = float(top_quantile)
        if not (0.0 < q <= 1.0):
            raise ValueError(f"top_quantile must be in (0, 1], got {top_quantile!r}")
        k_keep = max(1, int(math.ceil(q * target_auc.size)))
        sorted_desc = np.sort(target_auc)[::-1]
        threshold = float(sorted_desc[min(k_keep - 1, sorted_desc.size - 1)])
        threshold_meta = {
            "threshold_method": "numeric_top_quantile",
            "top_quantile": q,
            "threshold_selected": threshold,
            "n_target_blocks_total": int(target_auc.size),
        }

    kept = [b for b in target_blocks if b["total_signal"] >= threshold]
    with open(out_path, "w", encoding="utf-8") as fh:
        for b in kept:
            fh.write(
                f"{b['chrom']}\t{b['start']}\t{b['end']}\t"
                f"{b['total_signal']:.4f}\t{b['max_signal']:.4f}\t{b['max_region']}\n"
            )

    return {
        "status": "success",
        "peak_caller": "seacr",
        "mode": str(mode),
        "norm": effective_norm,
        "already_spikein_scaled": used_spikein,
        "n_peaks": len(kept),
        "peaks_path": str(out_path),
        "threshold": round(float(threshold), 6),
        "diagnostics": threshold_meta,
        "warnings": warnings_out,
    }


def _read_scored_peaks(path: str | Path) -> list[tuple[str, int, int, float]]:
    """Read BED/narrowPeak/broadPeak/SEACR peak intervals with signal score.

    Uses MACS3 ``signalValue`` (column 7, 0-based index 6) when >=7 columns are
    present (narrowPeak/broadPeak), SEACR ``total_signal`` (column 4, index 3)
    when numeric in 6-column SEACR output, or column 5 score otherwise.
    """
    peaks: list[tuple[str, int, int, float]] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            try:
                chrom = cols[0]
                start = int(float(cols[1]))
                end = int(float(cols[2]))
            except ValueError:
                continue
            if end <= start:
                continue
            score = 1.0
            if len(cols) >= 7:
                try:
                    score = float(cols[6])
                except ValueError:
                    score = 1.0
            elif len(cols) >= 5:
                for idx in (3, 4):
                    try:
                        score = float(cols[idx])
                        break
                    except ValueError:
                        continue
            elif len(cols) == 4:
                try:
                    score = float(cols[3])
                except ValueError:
                    score = 1.0
            peaks.append((chrom, start, end, float(score)))
    return peaks


def _pair_two_replicate_peaks(
    peaks1: list[tuple[str, int, int, float]],
    peaks2: list[tuple[str, int, int, float]],
) -> list[tuple[str, int, int, float, float]]:
    """Build matched overlapping peak pairs ``(chrom, start, end, s1, s2)``."""
    tagged: list[tuple[str, int, int, int, float]] = []
    for chrom, s, e, sc in peaks1:
        tagged.append((chrom, s, e, 0, sc))
    for chrom, s, e, sc in peaks2:
        tagged.append((chrom, s, e, 1, sc))
    tagged.sort(key=lambda r: (r[0], r[1], r[2], r[3]))

    paired: list[tuple[str, int, int, float, float]] = []
    cur_c: str | None = None
    cur_s = 0
    cur_e = 0
    best0: float | None = None
    best1: float | None = None

    def _emit() -> None:
        if cur_c is not None and best0 is not None and best1 is not None:
            paired.append((cur_c, cur_s, cur_e, float(best0), float(best1)))

    for chrom, s, e, rep, sc in tagged:
        if cur_c != chrom or s >= cur_e:
            _emit()
            cur_c, cur_s, cur_e = chrom, s, e
            best0 = sc if rep == 0 else None
            best1 = sc if rep == 1 else None
        else:
            cur_e = max(cur_e, e)
            if rep == 0:
                best0 = sc if best0 is None else max(best0, sc)
            else:
                best1 = sc if best1 is None else max(best1, sc)
    _emit()
    return paired


def _invert_mixture_cdf(
    u: np.ndarray,
    pi1: float,
    mu: float,
    sigma: float,
) -> np.ndarray:
    """Invert ``G(z) = (1 - pi1)*Phi(z) + pi1*Phi((z - mu)/sigma) = u`` via
    bracketed vectorized Newton-Raphson to machine precision."""
    from scipy.special import ndtr, ndtri

    u_clip = np.clip(u, 1e-10, 1.0 - 1e-10)
    z0 = ndtri(u_clip)
    z1 = mu + sigma * z0
    lo = np.minimum(z0, z1) - 1.0
    hi = np.maximum(z0, z1) + 1.0
    z = (1.0 - pi1) * z0 + pi1 * z1
    inv_sqrt2pi = 1.0 / math.sqrt(2.0 * math.pi)

    for _ in range(16):
        cdf = (1.0 - pi1) * ndtr(z) + pi1 * ndtr((z - mu) / sigma)
        diff = cdf - u_clip
        hi = np.where(diff > 0, z, hi)
        lo = np.where(diff < 0, z, lo)
        pdf = (1.0 - pi1) * inv_sqrt2pi * np.exp(-0.5 * z * z) + (
            pi1 / sigma
        ) * inv_sqrt2pi * np.exp(-0.5 * ((z - mu) / sigma) ** 2)
        newton = z - diff / np.maximum(pdf, 1e-15)
        z = np.where((newton > lo) & (newton < hi), newton, 0.5 * (lo + hi))
    return z


def fit_copula_idr_pair(
    scores1: np.ndarray,
    scores2: np.ndarray,
    *,
    max_iter: int = 80,
    tol: float = 1e-4,
) -> dict[str, Any]:
    """Fit Li et al. (2011) bimodal Gaussian copula EM model on paired scores."""
    from scipy.stats import rankdata

    x1 = np.asarray(scores1, dtype=float)
    x2 = np.asarray(scores2, dtype=float)
    n = int(x1.size)
    if n < 20 or x2.size != n:
        raise ValueError(f"IDR copula requires >=20 paired peaks, got {n}.")

    u1 = rankdata(x1, method="average") / (n + 1.0)
    u2 = rankdata(x2, method="average") / (n + 1.0)

    pi1 = 0.35
    mu = 2.0
    sigma = 1.0
    rho = 0.65
    converged = False
    n_iter = 0

    for it in range(1, int(max_iter) + 1):
        n_iter = it
        z1 = _invert_mixture_cdf(u1, pi1, mu, sigma)
        z2 = _invert_mixture_cdf(u2, pi1, mu, sigma)

        log_f0 = -math.log(2.0 * math.pi) - 0.5 * (z1 * z1 + z2 * z2)
        one_minus_r2 = max(1e-6, 1.0 - rho * rho)
        d1 = z1 - mu
        d2 = z2 - mu
        quad = (d1 * d1 - 2.0 * rho * d1 * d2 + d2 * d2) / (sigma * sigma * one_minus_r2)
        log_f1 = -math.log(2.0 * math.pi * sigma * sigma * math.sqrt(one_minus_r2)) - 0.5 * quad

        a0 = math.log(max(1e-8, 1.0 - pi1)) + log_f0
        a1 = math.log(max(1e-8, pi1)) + log_f1
        m_max = np.maximum(a0, a1)
        log_sum = m_max + np.log(np.exp(a0 - m_max) + np.exp(a1 - m_max))
        gamma = np.clip(np.exp(a1 - log_sum), 1e-9, 1.0 - 1e-9)

        w_sum = float(np.sum(gamma))
        if w_sum <= 1e-6:
            pi1_new, mu_new, sigma_new, rho_new = 0.01, mu, sigma, 0.0
            converged = True
            break

        pi1_new = float(np.clip(w_sum / n, 0.01, 0.99))
        mu_new = float(np.clip(np.sum(gamma * (z1 + z2)) / (2.0 * w_sum), 0.15, 8.0))
        var_new = float(
            np.clip(
                np.sum(gamma * ((z1 - mu_new) ** 2 + (z2 - mu_new) ** 2)) / (2.0 * w_sum),
                0.10,
                16.0,
            )
        )
        sigma_new = math.sqrt(var_new)
        cov_new = float(np.sum(gamma * (z1 - mu_new) * (z2 - mu_new)) / w_sum)
        rho_new = float(np.clip(cov_new / var_new, 0.0, 0.99))

        delta = max(
            abs(pi1_new - pi1),
            abs(mu_new - mu),
            abs(sigma_new - sigma),
            abs(rho_new - rho),
        )
        pi1, mu, sigma, rho = pi1_new, mu_new, sigma_new, rho_new
        if delta < tol:
            converged = True
            break

    # Final E-step at converged parameters
    z1 = _invert_mixture_cdf(u1, pi1, mu, sigma)
    z2 = _invert_mixture_cdf(u2, pi1, mu, sigma)
    log_f0 = -math.log(2.0 * math.pi) - 0.5 * (z1 * z1 + z2 * z2)
    one_minus_r2 = max(1e-6, 1.0 - rho * rho)
    d1 = z1 - mu
    d2 = z2 - mu
    quad = (d1 * d1 - 2.0 * rho * d1 * d2 + d2 * d2) / (sigma * sigma * one_minus_r2)
    log_f1 = -math.log(2.0 * math.pi * sigma * sigma * math.sqrt(one_minus_r2)) - 0.5 * quad
    a0 = math.log(max(1e-8, 1.0 - pi1)) + log_f0
    a1 = math.log(max(1e-8, pi1)) + log_f1
    m_max = np.maximum(a0, a1)
    gamma = np.clip(np.exp(a1 - (m_max + np.log(np.exp(a0 - m_max) + np.exp(a1 - m_max)))), 0.0, 1.0)

    local_idr = 1.0 - gamma
    # Tie-break identical local_idr by descending mean rank so stronger peaks come first
    mean_u = 0.5 * (u1 + u2)
    order = np.lexsort((-mean_u, local_idr))
    sorted_local = local_idr[order]
    sorted_global = np.cumsum(sorted_local) / np.arange(1, n + 1, dtype=float)
    global_idr = np.empty(n, dtype=float)
    global_idr[order] = sorted_global

    return {
        "local_idr": local_idr,
        "global_idr": global_idr,
        "params": {
            "pi1": round(pi1, 5),
            "mu": round(mu, 5),
            "sigma": round(sigma, 5),
            "rho": round(rho, 5),
        },
        "converged": bool(converged),
        "n_iter": int(n_iter),
        "n_paired_peaks": n,
    }


def compute_idr_reproducible_peaks(
    peak_files: list[str],
    output_path: str | Path,
    *,
    idr_threshold: float = 0.05,
    max_pairs: int = 6,
) -> dict[str, Any]:
    """Compute Li et al. (2011) Gaussian-copula IDR consensus across >=2 replicates."""
    if len(peak_files) < 2:
        return {
            "ran": False,
            "reason": "fewer_than_two_replicate_peak_files",
            "n_replicate_peak_files": len(peak_files),
        }

    scored = [_read_scored_peaks(p) for p in peak_files]
    pair_indices: list[tuple[int, int]] = []
    for i in range(len(scored)):
        for j in range(i + 1, len(scored)):
            pair_indices.append((i, j))
    if len(pair_indices) > max_pairs:
        # Adjacent + wraparound pairs keep runtime O(R) on large cohorts
        pair_indices = [(i, (i + 1) % len(scored)) for i in range(len(scored))]

    pair_reports: list[dict[str, Any]] = []
    passing_intervals: list[tuple[str, int, int, float, float]] = []

    for i, j in pair_indices:
        paired = _pair_two_replicate_peaks(scored[i], scored[j])
        pair_key = f"rep{i + 1}_vs_rep{j + 1}"
        if len(paired) < 20:
            pair_reports.append({
                "pair": pair_key,
                "ran": False,
                "n_paired_peaks": len(paired),
                "reason": "insufficient_paired_peaks_for_copula_idr",
            })
            continue

        s1 = np.array([r[3] for r in paired], dtype=float)
        s2 = np.array([r[4] for r in paired], dtype=float)
        fit = fit_copula_idr_pair(s1, s2)
        g_idr = fit["global_idr"]
        l_idr = fit["local_idr"]
        keep_idx = np.flatnonzero(g_idr <= float(idr_threshold))
        for idx in keep_idx:
            c, s, e, v1, v2 = paired[int(idx)]
            passing_intervals.append((c, s, e, 0.5 * (v1 + v2), float(g_idr[int(idx)])))

        pair_reports.append({
            "pair": pair_key,
            "ran": True,
            "n_paired_peaks": len(paired),
            "n_passing_idr": int(keep_idx.size),
            "idr_threshold": float(idr_threshold),
            "converged": fit["converged"],
            "n_iter": fit["n_iter"],
            "copula_params": fit["params"],
            "median_local_idr_passing": (
                round(float(np.median(l_idr[keep_idx])), 6) if keep_idx.size else None
            ),
        })

    ran_pairs = [p for p in pair_reports if p.get("ran")]
    if not ran_pairs:
        return {
            "ran": False,
            "method": "gaussian_copula_em_idr",
            "reason": "insufficient_paired_peaks_for_copula_idr",
            "pair_reports": pair_reports,
        }

    # Merge passing intervals across replicate pairs (ENCODE conservative/optimal union)
    passing_intervals.sort(key=lambda r: (r[0], r[1], r[2]))
    merged: list[tuple[str, int, int, float, float]] = []
    for c, s, e, mean_sig, gidr in passing_intervals:
        if not merged or merged[-1][0] != c or s > merged[-1][2]:
            merged.append((c, s, e, mean_sig, gidr))
        else:
            pc, ps, pe, psig, pidr = merged[-1]
            merged[-1] = (pc, ps, max(pe, e), max(psig, mean_sig), min(pidr, gidr))

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        for idx, (c, s, e, sig, gidr) in enumerate(merged, start=1):
            neg_log10_idr = -math.log10(max(gidr, 1e-12))
            score_int = min(1000, max(0, int(round(neg_log10_idr * 125.0))))
            fh.write(
                f"{c}\t{s}\t{e}\tidr_peak_{idx}\t{score_int}\t.\t"
                f"{sig:.4f}\t-1\t{neg_log10_idr:.4f}\t{(e - s) // 2}\n"
            )

    counts_passing = [p["n_passing_idr"] for p in ran_pairs]
    self_consistency_ratio = (
        round(float(max(counts_passing) / max(1, min(counts_passing))), 4)
        if len(counts_passing) >= 2 and min(counts_passing) > 0
        else None
    )
    return {
        "ran": True,
        "method": "gaussian_copula_em_idr",
        "idr_threshold": float(idr_threshold),
        "n_replicate_peak_files": len(peak_files),
        "n_pairs_evaluated": len(ran_pairs),
        "n_reproducible_peaks": len(merged),
        "self_consistency_ratio": self_consistency_ratio,
        "output_path": str(out),
        "pair_reports": pair_reports,
    }


def call_broad_domains(
    ip_bams: list[str],
    output_bed: str | Path,
    *,
    control_bams: list[str] | None = None,
    bin_size: int = 10_000,
    min_log2_enrichment: float = 0.585,
    max_gap_bp: int = 30_000,
    min_domain_bp: int = 20_000,
    ip_scale_factors: list[float] | None = None,
    control_scale_factors: list[float] | None = None,
    reference_fasta: str | None = None,
) -> dict[str, Any]:
    """Call broad heterochromatin / lamina megadomains (e.g. H3K9me3, Lamin B1).

    Quantifies primary non-duplicate reads in fixed genomic windows (default
    10 kb), normalizes by CPM (or explicit spike-in scale factors when supplied),
    computes ``log2((IP + pseudocount) / (Control + pseudocount))`` (or enrichment
    over chromosome median background when no control BAM is provided), stitches
    enriched bins across gaps ``<= max_gap_bp``, and retains domains ``>= min_domain_bp``.
    """
    if not ip_bams:
        return {"status": "error", "error_type": "NoInputFiles", "details": "No IP BAMs provided."}

    win = max(1000, int(bin_size))
    mito = {"chrm", "mt", "m"}

    def _bin_cohort(paths: list[str], scales: list[float] | None) -> tuple[dict[str, np.ndarray], dict[str, int]]:
        chrom_sizes: dict[str, int] = {}
        accum: dict[str, np.ndarray] = {}
        for idx, p in enumerate(paths):
            sample_bins: dict[str, np.ndarray] = {}
            total_mapped = 0
            with _open_alignment(p, reference_fasta) as bam:
                for ref, length in zip(bam.header.references, bam.header.lengths):
                    if ref.lower() in mito or ref.lower().startswith("chrm"):
                        continue
                    chrom_sizes[str(ref)] = int(length)
                    n_bins = max(1, int(math.ceil(int(length) / win)))
                    sample_bins[str(ref)] = np.zeros(n_bins, dtype=np.float64)
                for read in bam.fetch(until_eof=True):
                    if (
                        read.is_unmapped
                        or read.is_secondary
                        or read.is_supplementary
                        or read.is_duplicate
                    ):
                        continue
                    ref = read.reference_name
                    if not ref or ref not in sample_bins:
                        continue
                    total_mapped += 1
                    pos = int(read.reference_start)
                    b_idx = min(sample_bins[ref].size - 1, max(0, pos // win))
                    sample_bins[ref][b_idx] += 1.0

            scale = (
                float(scales[idx])
                if (scales and idx < len(scales))
                else (1e6 / max(1, total_mapped))
            )
            for ref, arr in sample_bins.items():
                if ref not in accum:
                    accum[ref] = arr * (scale / len(paths))
                else:
                    accum[ref] += arr * (scale / len(paths))
        return accum, chrom_sizes

    ip_norm, chrom_sizes = _bin_cohort(ip_bams, ip_scale_factors)
    ctrl_norm = (
        _bin_cohort(control_bams, control_scale_factors)[0]
        if control_bams
        else {}
    )

    pseudo = 0.5
    domains: list[tuple[str, int, int, float]] = []
    for chrom in sorted(ip_norm):
        ip_arr = ip_norm[chrom]
        chrom_len = chrom_sizes.get(chrom, ip_arr.size * win)
        if chrom in ctrl_norm and ctrl_norm[chrom].size == ip_arr.size:
            bg_arr = ctrl_norm[chrom]
        else:
            nz = ip_arr[ip_arr > 0]
            bg_val = float(np.median(nz)) if nz.size else 1.0
            bg_arr = np.full_like(ip_arr, bg_val)

        l2fc = np.log2((ip_arr + pseudo) / (bg_arr + pseudo))
        enriched_idx = np.flatnonzero(l2fc >= float(min_log2_enrichment))
        if enriched_idx.size == 0:
            continue

        # Stitch adjacent enriched windows across gaps <= max_gap_bp
        cur_s = int(enriched_idx[0]) * win
        cur_e = min(chrom_len, (int(enriched_idx[0]) + 1) * win)
        vals = [float(l2fc[int(enriched_idx[0])])]

        for idx in enriched_idx[1:]:
            s = int(idx) * win
            e = min(chrom_len, (int(idx) + 1) * win)
            if s - cur_e <= int(max_gap_bp):
                cur_e = e
                vals.append(float(l2fc[int(idx)]))
            else:
                if cur_e - cur_s >= int(min_domain_bp):
                    domains.append((chrom, cur_s, cur_e, float(np.mean(vals))))
                cur_s, cur_e, vals = s, e, [float(l2fc[int(idx)])]
        if cur_e - cur_s >= int(min_domain_bp):
            domains.append((chrom, cur_s, cur_e, float(np.mean(vals))))

    out_path = Path(output_bed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        for i, (c, s, e, score) in enumerate(domains, start=1):
            score_int = min(1000, max(0, int(round(score * 200.0))))
            fh.write(f"{c}\t{s}\t{e}\tmegadomain_{i}\t{score_int}\t.\t{score:.4f}\t-1\t-1\n")

    return {
        "status": "success",
        "peak_caller": "broad_domain",
        "bin_size": win,
        "min_log2_enrichment": float(min_log2_enrichment),
        "max_gap_bp": int(max_gap_bp),
        "min_domain_bp": int(min_domain_bp),
        "n_peaks": len(domains),
        "peaks_path": str(out_path),
        "total_domain_bp": int(sum(e - s for _, s, e, _ in domains)),
    }


def write_normalized_bigwig(
    bam_path: str | Path,
    output_bw: str | Path,
    *,
    normalization: str = "CPM",
    spikein_scale_factor: float | None = None,
    effective_genome_size: int | None = None,
    keep_duplicates: bool = False,
    reference_fasta: str | None = None,
) -> dict[str, Any]:
    """Export a normalized BigWig signal track (CPM, RPGC, or heterologous spike-in)."""
    import pyBigWig

    out_bw = Path(output_bw)
    out_bw.parent.mkdir(parents=True, exist_ok=True)
    tmp_bg = out_bw.with_suffix(".tmp.bedgraph")

    # First pass at scale=1.0 to get exact fragment count + chrom header
    raw_info = bam_to_fragment_bedgraph(
        bam_path,
        tmp_bg,
        scale_factor=1.0,
        keep_duplicates=keep_duplicates,
        reference_fasta=reference_fasta,
    )
    chrom_lengths = raw_info["chrom_lengths"]
    n_frags = max(1, int(raw_info["n_fragments"]))

    norm_mode = str(normalization or "CPM").strip().upper()
    if spikein_scale_factor is not None:
        scale = float(spikein_scale_factor)
        norm_mode = "SPIKEIN"
    elif norm_mode == "RPGC" and effective_genome_size and effective_genome_size > 0:
        scale = float(effective_genome_size) / float(n_frags * 200.0)
    else:
        scale = 1e6 / float(n_frags)
        norm_mode = "CPM"

    header = [(c, int(length)) for c, length in chrom_lengths.items() if int(length) > 0]
    bw = pyBigWig.open(str(out_bw), "w")
    try:
        bw.addHeader(header)
        chroms: list[str] = []
        starts: list[int] = []
        ends: list[int] = []
        vals: list[float] = []
        chrom_set = {c: int(length) for c, length in header}
        with open(tmp_bg, encoding="utf-8") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 4 or cols[0] not in chrom_set:
                    continue
                c = cols[0]
                s = max(0, int(cols[1]))
                e = min(chrom_set[c], int(cols[2]))
                if e <= s:
                    continue
                chroms.append(c)
                starts.append(s)
                ends.append(e)
                vals.append(float(cols[3]) * scale)
        if chroms:
            bw.addEntries(chroms, starts, ends=ends, values=vals)
    finally:
        bw.close()
        tmp_bg.unlink(missing_ok=True)

    return {
        "bigwig_path": str(out_bw),
        "normalization": norm_mode,
        "scale_factor": round(scale, 8),
        "n_fragments": int(raw_info["n_fragments"]),
    }
