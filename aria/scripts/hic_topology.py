"""ARIA 3D Genome Topology & Multi-Replicate Differential Engine (ADR-082).

validation_level: production (stack="hic").

Computes chromosome-scoped, out-of-core 3D chromatin organization using
``cooler``, ``cooltools``, and ``bioframe``:

* ``compartments`` — cis eigenvector decomposition (``cooltools.api.eigdecomp.eigs_cis``)
  phased against local GC fraction (from ``genome_fasta``) or TSS gene density
  (from ``gtf``) so positive ``E1`` corresponds to active ``A`` chromatin;
  explicitly flags ``pc1_validated=False`` when neither biological phasing track
  is staged locally.
* ``tads`` — Crane et al. (2015) multi-scale insulation score and Li-thresholded
  prominent boundary calling (``cooltools.api.insulation.insulation``).
* ``insulation_calibration`` — deterministic single-pass multi-window boundary
  prominence calibration on ``chr1`` (no fabricated random numbers).
* ``loops`` — ``cooltools.api.expected.expected_cis`` + donut-kernel ``dots``
  calling (``<= 25 kb`` resolution) with optional GTF promoter-enhancer anchor
  annotation.
* ``differential_3d`` — replicate-aware moderated linear model + BH-FDR over
  phased ``E1`` compartment tracks (``A_to_B`` / ``B_to_A`` switches and
  quantitative shifts) and Crane ``log2_insulation_score`` tracks
  (``boundary_strengthened`` / ``boundary_weakened``). Descriptive-only without
  p-values when biological replicates fall below ``min_replicates_per_condition``.

Zero imports from ``aria.agents.*``, ``litellm``, or ``pydantic``.
"""

from __future__ import annotations

import gzip
import math
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from aria.scripts._base import mocks_allowed, run_script
from aria.scripts.hic_preprocess import _check_hic_stage_cache, _save_hic_stage_cache

STAGE_VERSION = 1
EXCLUDED_CHROMS = frozenset({"chrM", "MT", "chrMT", "M", "chrEBV", "EBV"})
_ATTR_GENE_NAME_RE = re.compile(r'gene_name\s+"?([^";]+)"?')
_ATTR_GENE_ID_RE = re.compile(r'gene_id\s+"?([^";]+)"?')


# ── Shared helpers ───────────────────────────────────────────────────────────


def _resolve_uri(fpath: str, resolution: int) -> tuple[str, int]:
    import cooler

    if str(fpath).lower().endswith(".mcool"):
        avail = sorted(
            int(p.split("/")[-1])
            for p in cooler.fileops.list_coolers(fpath)
            if p.split("/")[-1].isdigit()
        )
        if not avail:
            raise ValueError(f"No resolutions inside .mcool: {fpath}")
        chosen = int(resolution) if int(resolution) in avail else min(
            avail, key=lambda r: abs(r - int(resolution))
        )
        return f"{fpath}::resolutions/{chosen}", chosen
    clr = cooler.Cooler(fpath)
    return fpath, int(clr.binsize)


def _has_valid_weights(clr) -> bool:
    import numpy as np

    if "weight" not in list(clr.bins().columns):
        return False
    w = np.asarray(clr.bins()[:]["weight"], dtype=float)
    return bool(np.isfinite(w).any() and np.nanmax(w) > 0)


def _ensure_balanced_cooler(fpath: str, resolution: int, out_dir: Path):
    """Open ``fpath`` at ``resolution``; if unbalanced, stage a writable copy in
    ``out_dir`` and run ICE ``balance_cooler`` so ``cooltools`` has valid weights."""
    import cooler

    uri, chosen_res = _resolve_uri(fpath, resolution)
    clr = cooler.Cooler(uri)
    if _has_valid_weights(clr):
        return clr, chosen_res

    staged_dir = out_dir / "_staged_balanced"
    staged_dir.mkdir(parents=True, exist_ok=True)
    ext = ".mcool" if str(fpath).lower().endswith(".mcool") else ".cool"
    stem = Path(fpath).stem.replace(".mcool", "").replace(".cool", "")
    staged = staged_dir / f"{stem}.{chosen_res}{ext}"
    if not staged.exists():
        shutil.copy2(fpath, staged)
        staged.chmod(0o644)
    s_uri = f"{staged}::resolutions/{chosen_res}" if ext == ".mcool" else str(staged)
    s_clr = cooler.Cooler(s_uri)
    if not _has_valid_weights(s_clr):
        nbins = int(s_clr.info.get("nbins", 0))
        min_nnz = 10 if nbins >= 500 else max(1, min(3, nbins // 10))
        ignore_diags = 2 if nbins >= 20 else 1
        cooler.balance_cooler(
            s_clr,
            cis_only=False,
            ignore_diags=ignore_diags,
            min_nnz=min_nnz,
            mad_max=5,
            max_iters=300,
            store=True,
        )
        s_clr = cooler.Cooler(s_uri)
    return s_clr, chosen_res


def _get_chromosomes(clr, chromosomes) -> list[str]:
    all_chroms = [str(c) for c in clr.chromnames]
    auto_chroms = [
        c for c in all_chroms
        if c not in EXCLUDED_CHROMS
        and "random" not in c.lower()
        and not c.lower().startswith("chrun")
        and "_" not in c
    ]
    if not auto_chroms:
        auto_chroms = [c for c in all_chroms if c not in EXCLUDED_CHROMS] or all_chroms
    if chromosomes == "all" or chromosomes is None:
        return auto_chroms
    if isinstance(chromosomes, (list, tuple)):
        picked = [str(c) for c in chromosomes if str(c) in all_chroms]
        return picked or auto_chroms[:1]
    return auto_chroms


def _build_view_df(clr, chroms: list[str]):
    import pandas as pd

    rows = []
    for c in chroms:
        if c in clr.chromsizes:
            length = int(clr.chromsizes[c])
            if length > 0:
                rows.append({"chrom": c, "start": 0, "end": length, "name": c})
    return pd.DataFrame(rows)


def _parse_gtf_tss(gtf_path: str | None) -> list[tuple[str, int, str]]:
    """Extract unique ``(chrom, tss_0based, gene_symbol)`` records from GTF/GFF."""
    if not gtf_path or not Path(str(gtf_path)).is_file():
        return []
    p = Path(str(gtf_path))
    opener = gzip.open if p.name.lower().endswith(".gz") else open
    genes: list[tuple[str, int, str]] = []
    seen: set[tuple[str, str]] = set()
    with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat = parts[2].lower()
            if feat not in {"gene", "transcript", "mrna"}:
                continue
            chrom = parts[0]
            try:
                start_1 = int(parts[3])
                end_1 = int(parts[4])
            except ValueError:
                continue
            strand = parts[6]
            tss = max(0, (start_1 - 1) if strand != "-" else (end_1 - 1))
            attrs = parts[8]
            m_name = _ATTR_GENE_NAME_RE.search(attrs)
            m_id = _ATTR_GENE_ID_RE.search(attrs)
            gene = (m_name.group(1) if m_name else (m_id.group(1) if m_id else "")).strip()
            if not gene:
                continue
            key = (chrom, gene)
            if key in seen:
                continue
            seen.add(key)
            genes.append((chrom, tss, gene))
    return genes


def _build_phasing_track(clr, chroms: list[str], genome_fasta: str | None, gtf_path: str | None):
    """Build a ``(chrom, start, end, phasing_score)`` DataFrame for ``eigs_cis``.

    Returns ``(phasing_df | None, phasing_method, pc1_validated)``."""
    import numpy as np
    import pandas as pd

    bins_df = clr.bins()[:][["chrom", "start", "end"]].copy()
    bins_df["chrom"] = bins_df["chrom"].astype(str)
    mask = bins_df["chrom"].isin(set(chroms))

    # 1. Local genome FASTA GC fraction (gold standard for cooltools eigs_cis)
    if genome_fasta and Path(str(genome_fasta)).is_file():
        try:
            import pysam

            with pysam.FastaFile(str(genome_fasta)) as fa:
                refs = set(fa.references)
                shared = [c for c in chroms if c in refs]
                if shared:
                    gc_vals = np.full(len(bins_df), np.nan, dtype=float)
                    for ch in shared:
                        ch_len = int(fa.get_reference_length(ch))
                        seq = fa.fetch(ch, 0, ch_len).upper()
                        ch_idx = np.flatnonzero(bins_df["chrom"].to_numpy() == ch)
                        starts = bins_df["start"].to_numpy(dtype=int)[ch_idx]
                        ends = bins_df["end"].to_numpy(dtype=int)[ch_idx]
                        for idx_pos, s, e in zip(ch_idx, starts, ends):
                            sub = seq[max(0, s):min(ch_len, e)]
                            if not sub:
                                continue
                            atgc = (
                                sub.count("A") + sub.count("T")
                                + sub.count("G") + sub.count("C")
                            )
                            if atgc >= 20:
                                gc_vals[idx_pos] = (sub.count("G") + sub.count("C")) / float(atgc)
                    if np.isfinite(gc_vals[mask]).sum() >= 5 and np.nanstd(gc_vals[mask]) > 0:
                        bins_df["GC"] = gc_vals
                        return bins_df[["chrom", "start", "end", "GC"]], "gc_content", True
        except Exception:
            pass

    # 2. Local GTF TSS / gene density per bin
    tss_records = _parse_gtf_tss(gtf_path)
    if tss_records:
        binsize = int(clr.binsize)
        counts = np.zeros(len(bins_df), dtype=float)
        offsets: dict[str, int] = {}
        lengths: dict[str, int] = {}
        for ch, grp in bins_df.groupby("chrom", sort=False, observed=False):
            offsets[str(ch)] = int(grp.index[0])
            lengths[str(ch)] = len(grp)
        for ch, tss, _gene in tss_records:
            if ch not in offsets:
                continue
            b = tss // max(binsize, 1)
            if 0 <= b < lengths[ch]:
                counts[offsets[ch] + b] += 1.0
        if counts[mask].sum() > 0 and np.std(counts[mask]) > 0:
            bins_df["gene_density"] = counts
            return bins_df[["chrom", "start", "end", "gene_density"]], "gene_density", True

    return None, "unphased_arbitrary_sign", False


def _bh_fdr(pvals):
    import numpy as np

    p = np.asarray(pvals, dtype=float)
    n = p.size
    out = np.full(n, np.nan, dtype=float)
    valid = np.isfinite(p)
    m = int(valid.sum())
    if m == 0:
        return out
    pv = np.clip(p[valid], 0.0, 1.0)
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order]
    adj = ranked * float(m) / np.arange(1, m + 1, dtype=float)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    inv = np.empty(m, dtype=int)
    inv[order] = np.arange(m)
    out[valid] = adj[inv]
    return out


# ── 1. Compartments A/B ──────────────────────────────────────────────────────


def _compute_compartments(
    fpath: str,
    genome: str,
    organism: str,
    resolution: int,
    chromosomes,
    out_dir: str,
    genome_fasta: str | None = None,
    gtf_path: str | None = None,
    allow_mock: bool = False,
) -> dict:
    try:
        import cooler  # noqa: F401
        import cooltools.api.eigdecomp as eig
        import numpy as np
    except ImportError as exc:
        if allow_mock:
            return _mock_compartments(resolution, str(exc))
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": f"Compartment analysis requires cooler/cooltools: {exc}",
        }

    try:
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        clr, actual_res = _ensure_balanced_cooler(fpath, resolution, out_path)
        chroms = _get_chromosomes(clr, chromosomes)
        view_df = _build_view_df(clr, chroms)
        if view_df.empty:
            return {
                "status": "error",
                "error_type": "NoValidChromosomes",
                "details": "No analyzable nuclear chromosomes found in cooler.",
            }

        phasing_df, phasing_method, pc1_validated = _build_phasing_track(
            clr, chroms, genome_fasta, gtf_path
        )
        nbins_min = min(int(clr.chromsizes[c]) // max(actual_res, 1) for c in view_df["chrom"])
        n_eigs = min(3, max(1, nbins_min - 2))
        ignore_diags = 2 if nbins_min >= 12 else 1

        eigvals_df, eigvecs_df = eig.eigs_cis(
            clr,
            phasing_track=phasing_df,
            view_df=view_df,
            n_eigs=n_eigs,
            ignore_diags=ignore_diags,
        )

        compartment_tracks: dict[str, list[float | None]] = {}
        a_bins = 0
        b_bins = 0
        n_transitions = 0

        for ch in view_df["chrom"]:
            sub = eigvecs_df[eigvecs_df["chrom"].astype(str) == str(ch)]
            if sub.empty or "E1" not in sub.columns:
                continue
            vals = sub["E1"].to_numpy(dtype=float)
            finite = vals[np.isfinite(vals)]
            if finite.size == 0:
                continue
            a_bins += int((finite > 0).sum())
            b_bins += int((finite < 0).sum())
            signs = np.sign(finite[finite != 0])
            if signs.size >= 2:
                n_transitions += int(np.sum(signs[1:] != signs[:-1]))
            compartment_tracks[str(ch)] = [
                round(float(v), 6) if math.isfinite(float(v)) else None
                for v in vals
            ]

        total_bins = a_bins + b_bins
        if total_bins == 0:
            return {
                "status": "error",
                "error_type": "CompartmentFailed",
                "details": "Eigenvector decomposition yielded zero finite E1 bins across all chromosomes.",
            }

        pct_A = a_bins / total_bins * 100.0
        pct_B = b_bins / total_bins * 100.0

        export_df = eigvecs_df[["chrom", "start", "end", "E1"]].copy()
        export_df["compartment"] = np.where(
            export_df["E1"] > 0, "A",
            np.where(export_df["E1"] < 0, "B", "unassigned"),
        )
        comp_tsv = out_path / "compartments.tsv"
        export_df.to_csv(comp_tsv, sep="\t", index=False)

        bg_df = export_df[np.isfinite(export_df["E1"])][["chrom", "start", "end", "E1"]]
        comp_bg = out_path / "compartments_E1.bedGraph"
        bg_df.to_csv(comp_bg, sep="\t", header=False, index=False)

        warnings: list[str] = []
        note = (
            f"E1 phased via {phasing_method} (positive E1 = active A compartment)."
            if pc1_validated else
            "PC1 sign is unphased (neither genome_fasta nor GTF staged); "
            "provide genome_fasta or gtf to validate A/B polarity."
        )
        if not pc1_validated:
            warnings.append("unphased_E1_sign_arbitrary")

        return {
            "status": "success",
            "ran": True,
            "data_type": "HiC",
            "analysis": "compartments",
            "validation_level": "production",
            "pct_A": round(float(pct_A), 2),
            "pct_B": round(float(pct_B), 2),
            "n_valid_bins": int(total_bins),
            "n_ab_transitions": int(n_transitions),
            "n_ab_switches": 0,  # single-sample run; cross-condition switches live in differential_3d
            "pc1_validated": bool(pc1_validated),
            "phasing_method": phasing_method,
            "compartments_tsv": str(comp_tsv),
            "compartments_bedgraph": str(comp_bg),
            "compartment_tracks": compartment_tracks,
            "resolution": int(actual_res),
            "warnings": warnings,
            "note": note,
        }
    except Exception as exc:
        return {
            "status": "error",
            "error_type": "CompartmentFailed",
            "details": str(exc)[:500],
        }


# ── 2. TAD calling — Insulation Score (Crane 2015) ───────────────────────────


def _compute_tads(
    fpath: str,
    genome: str,
    resolution: int,
    window_size: int,
    chromosomes,
    out_dir: str,
    allow_mock: bool = False,
) -> dict:
    try:
        import cooler  # noqa: F401
        import cooltools.api.insulation as ins
        import numpy as np
    except ImportError as exc:
        if allow_mock:
            return _mock_tads(resolution, window_size, str(exc))
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": f"TAD analysis requires cooler/cooltools: {exc}",
        }

    try:
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        clr, actual_res = _ensure_balanced_cooler(fpath, resolution, out_path)
        window_bp = max(int(window_size), int(actual_res) * 3)
        chroms = _get_chromosomes(clr, chromosomes)
        view_df = _build_view_df(clr, chroms)
        if view_df.empty:
            return {
                "status": "error",
                "error_type": "NoValidChromosomes",
                "details": "No analyzable chromosomes found for TAD calling.",
            }

        nbins_min = min(int(clr.chromsizes[c]) // max(actual_res, 1) for c in view_df["chrom"])
        ignore_diags = 2 if nbins_min >= 12 else 1
        ins_table = ins.insulation(
            clr,
            [window_bp],
            view_df=view_df,
            ignore_diags=ignore_diags,
            nproc=1,
        )

        score_col = f"log2_insulation_score_{window_bp}"
        bound_col = f"is_boundary_{window_bp}"
        strength_col = f"boundary_strength_{window_bp}"

        if score_col not in ins_table.columns or ins_table[score_col].notna().sum() == 0:
            return {
                "status": "error",
                "error_type": "TADFailed",
                "details": f"Insulation score returned zero valid bins at window={window_bp}bp.",
            }

        ins_tsv = out_path / "insulation.tsv"
        ins_table.to_csv(ins_tsv, sep="\t", index=False)

        all_boundaries: dict[str, list[int]] = {}
        tad_rows: list[tuple[str, int, int, str, float]] = []
        tad_sizes_kb: list[float] = []
        strengths: list[float] = []

        for ch in view_df["chrom"]:
            ch_df = ins_table[ins_table["chrom"].astype(str) == str(ch)]
            b_df = ch_df[ch_df[bound_col].fillna(False).astype(bool)]
            positions = [int(x) for x in b_df["start"].tolist()]
            all_boundaries[str(ch)] = positions
            if strength_col in b_df.columns:
                st_vals = b_df[strength_col].dropna().to_numpy(dtype=float)
                strengths.extend(float(x) for x in st_vals if math.isfinite(float(x)))
            st_map = (
                dict(zip(b_df["start"].astype(int), b_df[strength_col].fillna(0.0).astype(float)))
                if strength_col in b_df.columns else {}
            )
            if len(positions) >= 2:
                for idx in range(len(positions) - 1):
                    s_bp = positions[idx]
                    e_bp = positions[idx + 1]
                    if e_bp > s_bp:
                        size_kb = (e_bp - s_bp) / 1000.0
                        tad_sizes_kb.append(size_kb)
                        mean_st = 0.5 * (st_map.get(s_bp, 0.0) + st_map.get(e_bp, 0.0))
                        tad_rows.append((str(ch), s_bp, e_bp, f"TAD_{ch}_{idx+1}", round(mean_st, 4)))

        tads_bed = out_path / "tads.bed"
        with open(tads_bed, "w", encoding="utf-8") as fh:
            for ch, s_bp, e_bp, name, st in tad_rows:
                fh.write(f"{ch}\t{s_bp}\t{e_bp}\t{name}\t{st}\n")

        n_tads = len(tad_sizes_kb)
        n_boundaries = sum(len(v) for v in all_boundaries.values())
        median_size = float(np.median(tad_sizes_kb)) if tad_sizes_kb else 0.0
        mean_strength = round(float(np.mean(strengths)), 4) if strengths else 0.0

        return {
            "status": "success",
            "ran": True,
            "data_type": "HiC",
            "analysis": "tads",
            "validation_level": "production",
            "n_tads": int(n_tads),
            "n_boundaries": int(n_boundaries),
            "median_size_kb": round(median_size, 1),
            "mean_boundary_strength": mean_strength,
            "tad_boundaries": all_boundaries,
            "tads_bed": str(tads_bed),
            "insulation_tsv": str(ins_tsv),
            "window_size": int(window_bp),
            "resolution": int(actual_res),
            "algorithm": "Insulation Score (Crane 2015, Li prominence threshold)",
        }
    except Exception as exc:
        return {
            "status": "error",
            "error_type": "TADFailed",
            "details": str(exc)[:500],
        }


# ── 3. Insulation Score window calibration ───────────────────────────────────


def _calibrate_insulation(
    fpath: str,
    resolution: int,
    windows: list[int],
    chromosomes,
    out_dir: str,
    allow_mock: bool = False,
) -> dict:
    try:
        import cooler  # noqa: F401
        import cooltools.api.insulation as ins
        import numpy as np
    except ImportError as exc:
        if allow_mock:
            return {
                "status": "success",
                "boundary_strengths": {str(w): 0.5 for w in windows},
                "calibration_chrom": "chr1",
                "resolution": resolution,
                "validation_level": "scaffold",
                "note": f"Mock calibration — {exc}",
            }
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": f"Insulation calibration requires cooler/cooltools: {exc}",
        }

    try:
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        clr, actual_res = _ensure_balanced_cooler(fpath, resolution, out_path)
        chroms = _get_chromosomes(clr, chromosomes)
        cal_chrom = "chr1" if "chr1" in chroms else chroms[0]
        view_df = _build_view_df(clr, [cal_chrom])
        valid_windows = [max(int(w), int(actual_res) * 3) for w in windows]
        nbins = int(clr.chromsizes[cal_chrom]) // max(actual_res, 1)
        ignore_diags = 2 if nbins >= 12 else 1

        table = ins.insulation(
            clr,
            sorted(set(valid_windows)),
            view_df=view_df,
            ignore_diags=ignore_diags,
            nproc=1,
        )

        boundary_strengths: dict[str, float] = {}
        for orig_w, eff_w in zip(windows, valid_windows):
            bound_col = f"is_boundary_{eff_w}"
            st_col = f"boundary_strength_{eff_w}"
            sc_col = f"log2_insulation_score_{eff_w}"
            b_rows = (
                table[table[bound_col].fillna(False).astype(bool)]
                if bound_col in table.columns else table.iloc[0:0]
            )
            if not b_rows.empty and st_col in b_rows.columns and b_rows[st_col].notna().any():
                val = float(np.nanmean(b_rows[st_col].to_numpy(dtype=float)))
            elif not b_rows.empty and sc_col in b_rows.columns and b_rows[sc_col].notna().any():
                val = float(np.nanmean(np.abs(b_rows[sc_col].to_numpy(dtype=float))))
            else:
                val = 0.0
            boundary_strengths[str(orig_w)] = round(val if math.isfinite(val) else 0.0, 4)

        return {
            "status": "success",
            "ran": True,
            "data_type": "HiC",
            "analysis": "insulation_calibration",
            "validation_level": "production",
            "boundary_strengths": boundary_strengths,
            "calibration_chrom": cal_chrom,
            "resolution": int(actual_res),
            "windows_tested": [int(w) for w in windows],
        }
    except Exception as exc:
        return {
            "status": "error",
            "error_type": "CalibrationFailed",
            "details": str(exc)[:500],
        }


# ── 4. Chromatin loop calling (cooltools expected_cis + dots) ────────────────


def _annotate_loops_with_gtf(loops_df, gtf_path: str | None, promoter_window_bp: int = 5000):
    tss_list = _parse_gtf_tss(gtf_path)
    if not tss_list or loops_df.empty:
        return loops_df, {}
    by_chrom: dict[str, list[tuple[int, str]]] = {}
    for ch, tss, gene in tss_list:
        by_chrom.setdefault(str(ch), []).append((int(tss), gene))

    def _hit_genes(chrom: str, start: int, end: int) -> list[str]:
        lo = int(start) - promoter_window_bp
        hi = int(end) + promoter_window_bp
        hits = [g for tss, g in by_chrom.get(str(chrom), ()) if lo <= tss <= hi]
        return sorted(set(hits))

    a1_genes: list[str] = []
    a2_genes: list[str] = []
    classes: list[str] = []
    counts = {"promoter_promoter": 0, "promoter_enhancer": 0, "distal_distal": 0}
    for row in loops_df.itertuples(index=False):
        g1 = _hit_genes(str(row.chrom1), int(row.start1), int(row.end1))
        g2 = _hit_genes(str(row.chrom2), int(row.start2), int(row.end2))
        a1_genes.append(",".join(g1))
        a2_genes.append(",".join(g2))
        if g1 and g2:
            cls = "promoter_promoter"
        elif g1 or g2:
            cls = "promoter_enhancer"
        else:
            cls = "distal_distal"
        classes.append(cls)
        counts[cls] += 1

    out = loops_df.copy()
    out["anchor1_genes"] = a1_genes
    out["anchor2_genes"] = a2_genes
    out["loop_class"] = classes
    return out, counts


def _ensure_toeplitz_shift_invariant(dotfinder_mod) -> None:
    """Fix upstream cooltools 0.7.1 + cooler 0.10.x global-bin offset bug on chr2+.

    ``dotfinder.score_tile`` slices ``LazyToeplitz(expected_chr)`` (shaped
    ``(n_bins_chr, n_bins_chr)``) with global cooler bin spans ``(start_i, end_i)``.
    For any chromosome after ``chr1`` where ``start_i >= n_bins_chr``,
    ``cooler._IndexingMixin._process_slice`` clamps both bounds to ``n_bins_chr``
    and yields ``(0, 0)``. Because any Toeplitz matrix depends solely on
    ``j - i``, shifting ``(slice_i, slice_j)`` by ``-min(start_i, start_j)`` is
    mathematically exact on all cis tiles.
    """
    lz_cls = getattr(dotfinder_mod, "LazyToeplitz", None)
    if lz_cls is None or getattr(lz_cls, "_aria_shift_patched", False):
        return
    orig_getitem = lz_cls.__getitem__

    def _shifted_getitem(self, key):
        if isinstance(key, tuple) and len(key) == 2:
            s0, s1 = key
            if (
                isinstance(s0, slice)
                and isinstance(s1, slice)
                and s0.step is None
                and s1.step is None
                and s0.start is not None
                and s1.start is not None
            ):
                m = min(int(s0.start), int(s1.start))
                if m > 0:
                    s0 = slice(int(s0.start) - m, None if s0.stop is None else int(s0.stop) - m)
                    s1 = slice(int(s1.start) - m, None if s1.stop is None else int(s1.stop) - m)
                    key = (s0, s1)
        return orig_getitem(self, key)

    lz_cls.__getitem__ = _shifted_getitem
    lz_cls._aria_shift_patched = True


def _compute_loops(
    fpath: str,
    genome: str,
    resolution: int,
    chromosomes,
    out_dir: str,
    gtf_path: str | None = None,
    fdr_threshold: float = 0.1,
    max_loci_separation: int = 2_000_000,
    allow_mock: bool = False,
) -> dict:
    try:
        import cooler  # noqa: F401
        from cooltools.api import dotfinder, expected
    except ImportError as exc:
        if allow_mock:
            return _mock_loops(resolution, str(exc))
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": f"Loop calling requires cooler/cooltools: {exc}",
        }

    try:
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        clr, actual_res = _ensure_balanced_cooler(fpath, resolution, out_path)
        if actual_res > 25_000:
            return {
                "status": "skipped",
                "ran": False,
                "data_type": "HiC",
                "analysis": "loops",
                "validation_level": "production",
                "resolution": int(actual_res),
                "n_loops": 0,
                "reason": "resolution_too_coarse_for_dot_calling",
                "message": (
                    f"Loop dot calling requires <= 25 kb resolution; "
                    f"finest available resolution was {actual_res:,} bp."
                ),
            }

        chroms = _get_chromosomes(clr, chromosomes)
        view_df = _build_view_df(clr, chroms)
        if view_df.empty:
            return {
                "status": "error",
                "error_type": "NoValidChromosomes",
                "details": "No analyzable chromosomes found for loop calling.",
            }

        exp_df = expected.expected_cis(
            clr,
            view_df=view_df,
            smooth=False,
            aggregate_smoothed=False,
            ignore_diags=2,
            nproc=1,
        )
        max_sep = min(
            int(max_loci_separation),
            max(int(clr.chromsizes[c]) for c in view_df["chrom"]),
        )
        max_sep = max(max_sep, actual_res * 10)
        tile_sz = max(max_sep * 2, 2_000_000)

        _ensure_toeplitz_shift_invariant(dotfinder)
        dots_df = dotfinder.dots(
            clr,
            expected=exp_df,
            expected_value_col="balanced.avg",
            view_df=view_df,
            max_loci_separation=max_sep,
            lambda_bin_fdr=float(fdr_threshold),
            tile_size=tile_sz,
            nproc=1,
        )

        loops_path = out_path / "loops.bedpe"
        dots_df.to_csv(loops_path, sep="\t", index=False)

        annotated_path: str | None = None
        loop_class_counts: dict[str, int] = {}
        if gtf_path and not dots_df.empty:
            ann_df, loop_class_counts = _annotate_loops_with_gtf(dots_df, gtf_path)
            if loop_class_counts:
                ann_file = out_path / "annotated_loops.tsv"
                ann_df.to_csv(ann_file, sep="\t", index=False)
                annotated_path = str(ann_file)
                dots_df = ann_df

        top_loops: list[dict[str, Any]] = []
        if not dots_df.empty:
            sort_cols = [c for c in ("la_exp.donut.qval", "count") if c in dots_df.columns]
            asc = [c != "count" for c in sort_cols]
            ranked = dots_df.sort_values(sort_cols, ascending=asc) if sort_cols else dots_df
            for _, row in ranked.head(20).iterrows():
                item: dict[str, Any] = {
                    "chrom1": str(row["chrom1"]),
                    "start1": int(row["start1"]),
                    "end1": int(row["end1"]),
                    "chrom2": str(row["chrom2"]),
                    "start2": int(row["start2"]),
                    "end2": int(row["end2"]),
                }
                if "count" in row and row["count"] is not None:
                    item["count"] = int(row["count"])
                if "la_exp.donut.qval" in row and math.isfinite(float(row["la_exp.donut.qval"])):
                    item["qval"] = round(float(row["la_exp.donut.qval"]), 6)
                if "anchor1_genes" in row and str(row["anchor1_genes"]):
                    item["anchor1_genes"] = str(row["anchor1_genes"])
                if "anchor2_genes" in row and str(row["anchor2_genes"]):
                    item["anchor2_genes"] = str(row["anchor2_genes"])
                if "loop_class" in row:
                    item["loop_class"] = str(row["loop_class"])
                top_loops.append(item)

        return {
            "status": "success",
            "ran": True,
            "data_type": "HiC",
            "analysis": "loops",
            "validation_level": "production",
            "n_loops": int(len(dots_df)),
            "top_loops": top_loops,
            "loops_path": str(loops_path),
            "annotated_loops_path": annotated_path,
            "loop_class_counts": loop_class_counts,
            "fdr_threshold": float(fdr_threshold),
            "algorithm": "cooltools expected_cis + donut dots",
            "resolution": int(actual_res),
            "ctcf_validated": False,
            "atac_corroborated": False,
        }
    except Exception as exc:
        if allow_mock:
            return _mock_loops(resolution, str(exc))
        return {
            "status": "error",
            "error_type": "LoopCallingFailed",
            "details": str(exc)[:500],
        }


# ── 5. Multi-Replicate Differential 3D Architecture (ΔE1 & ΔInsulation) ──────


def _moderated_linear_contrast(
    mat,  # shape (n_bins, n_replicates)
    conditions: list[str],
    test_group: str,
    ref_group: str,
    covariate_rows: list[dict[str, str]] | None = None,
    prior_df: float = 4.0,
):
    """Empirical-Bayes variance-moderated OLS contrast across genomic bins."""
    import numpy as np
    import pandas as pd
    from scipy import stats

    y = np.asarray(mat, dtype=float)
    n_bins, n_samples = y.shape
    cond_arr = np.asarray(conditions, dtype=str)
    keep_samples = (cond_arr == str(test_group)) | (cond_arr == str(ref_group))
    y = y[:, keep_samples]
    cond_sub = cond_arr[keep_samples]
    is_test = (cond_sub == str(test_group)).astype(float)

    df_design = pd.DataFrame({"intercept": 1.0, "is_test": is_test})
    dropped_covariates: list[str] = []
    if covariate_rows:
        cov_sub = [row for row, kp in zip(covariate_rows, keep_samples) if kp]
        cov_df = pd.DataFrame(cov_sub)
        for col in cov_df.columns:
            vals = cov_df[col].astype(str)
            if vals.nunique() < 2:
                dropped_covariates.append(f"{col}:constant_within_contrast")
                continue
            dummies = pd.get_dummies(vals, prefix=str(col), drop_first=True, dtype=float)
            cand = pd.concat([df_design, dummies], axis=1)
            cand_rank = int(np.linalg.matrix_rank(cand.to_numpy(dtype=float)))
            if cand_rank < cand.shape[1]:
                raise ValueError(
                    f"Covariate {col!r} is collinear/confounded with contrast "
                    f"{test_group}_vs_{ref_group}; refusing unadjusted inference."
                )
            if cand.shape[1] >= len(cond_sub):
                raise ValueError(
                    f"Covariate {col!r} saturates residual degrees of freedom "
                    f"for contrast {test_group}_vs_{ref_group}."
                )
            df_design = cand

    X = df_design.to_numpy(dtype=float)
    n_obs, n_params = X.shape
    df_resid = n_obs - n_params
    if df_resid < 1 or np.linalg.matrix_rank(X) < n_params:
        raise ValueError("Design matrix is rank-deficient or has zero residual degrees of freedom.")

    # Complete-case across replicates for clean QR vectorization, with fallback for sparse NaNs
    finite_mask = np.isfinite(y)
    n_test_valid = finite_mask[:, is_test == 1.0].sum(axis=1)
    n_ref_valid = finite_mask[:, is_test == 0.0].sum(axis=1)
    testable = (n_test_valid >= 2) & (n_ref_valid >= 2)

    betas = np.full(n_bins, np.nan, dtype=float)
    s2 = np.full(n_bins, np.nan, dtype=float)
    v_unscaled = np.full(n_bins, np.nan, dtype=float)
    d_vec = np.full(n_bins, np.nan, dtype=float)
    mean_test = np.full(n_bins, np.nan, dtype=float)
    mean_ref = np.full(n_bins, np.nan, dtype=float)

    all_complete = testable & finite_mask.all(axis=1)
    if np.any(all_complete):
        yc = y[all_complete]
        xtx_inv = np.linalg.inv(X.T @ X)
        coef = yc @ X @ xtx_inv.T  # (n_complete, n_params)
        fitted = coef @ X.T
        resid = yc - fitted
        rss = np.sum(resid * resid, axis=1)
        betas[all_complete] = coef[:, 1]
        s2[all_complete] = rss / float(df_resid)
        v_unscaled[all_complete] = float(xtx_inv[1, 1])
        d_vec[all_complete] = float(df_resid)
        mean_test[all_complete] = np.mean(yc[:, is_test == 1.0], axis=1)
        mean_ref[all_complete] = np.mean(yc[:, is_test == 0.0], axis=1)

    partial_idx = np.flatnonzero(testable & (~all_complete))
    for b in partial_idx:
        ok = finite_mask[b]
        xb = X[ok]
        yb = y[b, ok]
        if xb.shape[0] <= xb.shape[1] or np.linalg.matrix_rank(xb) < xb.shape[1]:
            continue
        xtx_inv_b = np.linalg.inv(xb.T @ xb)
        cb = xtx_inv_b @ (xb.T @ yb)
        rb = yb - xb @ cb
        dfb = xb.shape[0] - xb.shape[1]
        betas[b] = float(cb[1])
        s2[b] = float(np.dot(rb, rb) / dfb)
        v_unscaled[b] = float(xtx_inv_b[1, 1])
        d_vec[b] = float(dfb)
        mean_test[b] = float(np.mean(yb[is_test[ok] == 1.0]))
        mean_ref[b] = float(np.mean(yb[is_test[ok] == 0.0]))

    valid_s2 = s2[np.isfinite(s2) & (s2 > 0)]
    s0_sq = float(np.median(valid_s2)) if valid_s2.size > 0 else 1e-4
    s0_sq = max(s0_sq, 1e-6)

    s2_post = (prior_df * s0_sq + d_vec * s2) / (prior_df + d_vec)
    df_total = prior_df + d_vec
    se = np.sqrt(np.maximum(s2_post * v_unscaled, 1e-16))
    t_stat = betas / se
    pvals = 2.0 * stats.t.sf(np.abs(t_stat), df=df_total)
    padj = _bh_fdr(pvals)
    return {
        "delta": betas,
        "mean_test": mean_test,
        "mean_ref": mean_ref,
        "t_stat": t_stat,
        "pvalue": pvals,
        "padj": padj,
        "prior_s0_sq": s0_sq,
        "dropped_covariates": dropped_covariates,
    }


def _parse_comparisons(raw) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    if isinstance(raw, (list, tuple)):
        if len(raw) == 2 and all(isinstance(x, str) for x in raw):
            return [(str(raw[0]), str(raw[1]))]
        for item in raw:
            if isinstance(item, dict):
                num = item.get("test") or item.get("numerator")
                den = item.get("reference") or item.get("denominator")
                if num and den:
                    out.append((str(num), str(den)))
            elif isinstance(item, (list, tuple)) and len(item) == 2:
                out.append((str(item[0]), str(item[1])))
    elif isinstance(raw, dict):
        num = raw.get("test") or raw.get("numerator")
        den = raw.get("reference") or raw.get("denominator")
        if num and den:
            out.append((str(num), str(den)))
    return out


def _compute_differential_3d(params: dict, out_dir: str) -> dict:
    import cooler  # noqa: F401
    import cooltools.api.eigdecomp as eig
    import cooltools.api.insulation as ins
    import numpy as np
    import pandas as pd

    files = [str(f) for f in (params.get("files") or []) if Path(str(f)).exists()]
    if len(files) < 2:
        return {
            "status": "skipped",
            "ran": False,
            "data_type": "HiC",
            "analysis": "differential_3d",
            "validation_level": "production",
            "reason": "fewer_than_two_hic_samples",
        }

    sample_metadata = params.get("sample_metadata")
    meta_by_sid: dict[str, dict] = {}
    if isinstance(sample_metadata, (str, Path)) and Path(str(sample_metadata)).is_file():
        sep = "\t" if str(sample_metadata).lower().endswith((".tsv", ".txt")) else ","
        df_meta = pd.read_csv(sample_metadata, sep=sep)
        sid_col = next(
            (c for c in ("sample_id", "sample", "id") if c in df_meta.columns),
            df_meta.columns[0],
        )
        for row_dict in df_meta.to_dict(orient="records"):
            sid_val = row_dict.get(sid_col)
            if sid_val is not None:
                meta_by_sid[str(sid_val)] = row_dict
    elif isinstance(sample_metadata, dict):
        for sid, m in sample_metadata.items():
            if isinstance(m, dict):
                meta_by_sid[str(sid)] = m
    elif isinstance(sample_metadata, list):
        for m in sample_metadata:
            if isinstance(m, dict):
                sid_val = m.get("sample_id") or m.get("sample") or m.get("id")
                if sid_val is not None:
                    meta_by_sid[str(sid_val)] = m

    declared_sids = params.get("sample_ids")
    if isinstance(declared_sids, list) and len(declared_sids) == len(files):
        sample_ids = [str(s) for s in declared_sids]
    else:
        sample_ids = [Path(f).stem.replace(".balanced", "").replace(".mcool", "").replace(".cool", "") for f in files]

    comparisons = _parse_comparisons(params.get("comparisons"))
    if not meta_by_sid or not comparisons:
        return {
            "status": "skipped",
            "ran": False,
            "data_type": "HiC",
            "analysis": "differential_3d",
            "validation_level": "production",
            "reason": "missing_sample_metadata_or_comparisons",
        }

    cond_col = str(params.get("condition_col") or "condition")
    rep_col = str(params.get("replicate_col") or "replicate")
    covariates = [str(c) for c in (params.get("covariates") or []) if c]
    min_reps = max(2, int(params.get("min_replicates_per_condition") or 3))
    comp_res = int(params.get("comp_resolution") or params.get("resolution") or 100_000)
    tad_res = int(params.get("tad_resolution") or params.get("resolution") or 40_000)
    window_bp = int(params.get("window_size") or (tad_res * 5))
    delta_e1_min = float(params.get("delta_e1_min", 0.20))
    delta_ins_min = float(params.get("delta_insulation_min", 0.15))
    padj_max = float(params.get("padj_max", 0.05))
    genome_fasta = params.get("genome_fasta")
    gtf_path = params.get("gtf") or params.get("gtf_path")
    chromosomes = params.get("chromosomes", "all")

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Compute per-sample E1 and insulation tracks
    e1_by_sid: dict[str, pd.DataFrame] = {}
    ins_by_sid: dict[str, pd.DataFrame] = {}
    phasing_df = None
    pc1_validated = False
    phasing_method = "unphased_arbitrary_sign"

    for idx, (sid, fpath) in enumerate(zip(sample_ids, files)):
        if sid not in meta_by_sid:
            continue
        clr_c, actual_comp_res = _ensure_balanced_cooler(fpath, comp_res, out_path / sid)
        chroms = _get_chromosomes(clr_c, chromosomes)
        vdf_c = _build_view_df(clr_c, chroms)
        if idx == 0:
            phasing_df, phasing_method, pc1_validated = _build_phasing_track(
                clr_c, chroms, genome_fasta, gtf_path
            )
        nb_min = min(int(clr_c.chromsizes[c]) // max(actual_comp_res, 1) for c in vdf_c["chrom"])
        _, ev_df = eig.eigs_cis(
            clr_c,
            phasing_track=phasing_df,
            view_df=vdf_c,
            n_eigs=min(3, max(1, nb_min - 2)),
            ignore_diags=2 if nb_min >= 12 else 1,
        )
        # If unphased, align sign of sample i >= 1 to sample 0 PER CHROMOSOME so arbitrary
        # per-chromosome eigenvector flips do not create artificial cross-replicate deltas!
        if not pc1_validated and e1_by_sid:
            ref_df = next(iter(e1_by_sid.values()))
            ev_df = ev_df.copy()
            for ch_name in vdf_c["chrom"]:
                r_mask = ref_df["chrom"].astype(str) == str(ch_name)
                c_mask = ev_df["chrom"].astype(str) == str(ch_name)
                ref_ch = ref_df.loc[r_mask, "E1"].to_numpy(dtype=float)
                cur_ch = ev_df.loc[c_mask, "E1"].to_numpy(dtype=float)
                if ref_ch.shape == cur_ch.shape:
                    ok = np.isfinite(ref_ch) & np.isfinite(cur_ch)
                    if ok.sum() >= 5 and float(np.dot(ref_ch[ok], cur_ch[ok])) < 0:
                        ev_df.loc[c_mask, "E1"] = -cur_ch
        e1_by_sid[sid] = ev_df[["chrom", "start", "end", "E1"]].reset_index(drop=True)

        clr_t, actual_tad_res = _ensure_balanced_cooler(fpath, tad_res, out_path / sid)
        vdf_t = _build_view_df(clr_t, _get_chromosomes(clr_t, chromosomes))
        eff_w = max(window_bp, actual_tad_res * 3)
        nb_t = min(int(clr_t.chromsizes[c]) // max(actual_tad_res, 1) for c in vdf_t["chrom"])
        ins_df = ins.insulation(
            clr_t, [eff_w], view_df=vdf_t, ignore_diags=2 if nb_t >= 12 else 1, nproc=1
        )
        sc_col = f"log2_insulation_score_{eff_w}"
        ins_by_sid[sid] = ins_df[["chrom", "start", "end", sc_col]].rename(
            columns={sc_col: "insulation"}
        ).reset_index(drop=True)

    # Collapse technical replicates to biological replicates (mean track per replicate unit)
    rep_units: dict[tuple[str, str], list[str]] = {}
    rep_covs: dict[tuple[str, str], dict[str, str]] = {}
    for sid in e1_by_sid:
        m = meta_by_sid[sid]
        cond = m.get(cond_col) or m.get("condition")
        rep = m.get(rep_col) if m.get(rep_col) is not None else sid
        if cond is None:
            continue
        key = (str(cond), str(rep))
        rep_units.setdefault(key, []).append(sid)
        rep_covs.setdefault(key, {c: str(m.get(c, "")) for c in covariates})

    contrast_results: list[dict[str, Any]] = []
    warnings: list[str] = []
    if min_reps == 2:
        warnings.append("low_power_warning: min_replicates_per_condition=2")

    base_e1_coords = next(iter(e1_by_sid.values()))[["chrom", "start", "end"]]
    base_ins_coords = next(iter(ins_by_sid.values()))[["chrom", "start", "end"]]
    for sid, df_e1 in e1_by_sid.items():
        if len(df_e1) != len(base_e1_coords) or not df_e1[["chrom", "start", "end"]].equals(base_e1_coords):
            return {
                "status": "error",
                "error_type": "CoordinateMismatch",
                "details": (
                    f"Sample {sid!r} compartment bin grid does not match cohort reference "
                    "(mismatched resolution or chromosome assembly across Hi-C inputs)."
                ),
            }
    for sid, df_in in ins_by_sid.items():
        if len(df_in) != len(base_ins_coords) or not df_in[["chrom", "start", "end"]].equals(base_ins_coords):
            return {
                "status": "error",
                "error_type": "CoordinateMismatch",
                "details": (
                    f"Sample {sid!r} insulation bin grid does not match cohort reference "
                    "(mismatched resolution or chromosome assembly across Hi-C inputs)."
                ),
            }

    rep_keys = sorted(rep_units.keys())
    rep_conditions = [k[0] for k in rep_keys]
    rep_cov_list = [rep_covs[k] for k in rep_keys]

    e1_mat = np.column_stack([
        np.nanmean(
            np.column_stack([e1_by_sid[s]["E1"].to_numpy(dtype=float) for s in rep_units[k]]),
            axis=1,
        )
        for k in rep_keys
    ])
    ins_mat = np.column_stack([
        np.nanmean(
            np.column_stack([ins_by_sid[s]["insulation"].to_numpy(dtype=float) for s in rep_units[k]]),
            axis=1,
        )
        for k in rep_keys
    ])

    for test_grp, ref_grp in comparisons:
        n_test = sum(1 for c in rep_conditions if c == test_grp)
        n_ref = sum(1 for c in rep_conditions if c == ref_grp)
        tag = f"{test_grp}_vs_{ref_grp}"

        if n_test < min_reps or n_ref < min_reps:
            # Descriptive-only delta (never fabricates p-values from n < min_reps)
            t_idx = [i for i, c in enumerate(rep_conditions) if c == test_grp]
            r_idx = [i for i, c in enumerate(rep_conditions) if c == ref_grp]
            desc_csv = None
            if t_idx and r_idx:
                mt = np.nanmean(e1_mat[:, t_idx], axis=1)
                mr = np.nanmean(e1_mat[:, r_idx], axis=1)
                desc_df = base_e1_coords.copy()
                desc_df["mean_E1_test"] = mt
                desc_df["mean_E1_ref"] = mr
                desc_df["delta_E1"] = mt - mr
                desc_path = out_path / f"descriptive_compartments_{tag}.csv"
                desc_df.to_csv(desc_path, index=False)
                desc_csv = str(desc_path)
            contrast_results.append({
                "contrast": tag,
                "test": test_grp,
                "reference": ref_grp,
                "status": "skipped",
                "ran": False,
                "descriptive_only": True,
                "n_test_replicates": n_test,
                "n_ref_replicates": n_ref,
                "reason": (
                    f"insufficient_biological_replicates ({n_test} vs {n_ref}; "
                    f"requires >= {min_reps} per condition)"
                ),
                "descriptive_csv": desc_csv,
            })
            continue

        try:
            e1_fit = _moderated_linear_contrast(
                e1_mat, rep_conditions, test_grp, ref_grp, covariate_rows=rep_cov_list
            )
            ins_fit = _moderated_linear_contrast(
                ins_mat, rep_conditions, test_grp, ref_grp, covariate_rows=rep_cov_list
            )
        except ValueError as exc:
            return {
                "status": "error",
                "error_type": "ConfoundedOrRankDeficientDesign",
                "details": str(exc),
            }
        for dropped in e1_fit.get("dropped_covariates") or []:
            warnings.append(f"{tag}: dropped_covariate={dropped}")
        comp_df = base_e1_coords.copy()
        comp_df["mean_E1_test"] = e1_fit["mean_test"]
        comp_df["mean_E1_ref"] = e1_fit["mean_ref"]
        comp_df["delta_E1"] = e1_fit["delta"]
        comp_df["t_stat"] = e1_fit["t_stat"]
        comp_df["pvalue"] = e1_fit["pvalue"]
        comp_df["padj"] = e1_fit["padj"]
        sig_e1 = (
            np.isfinite(comp_df["padj"])
            & (comp_df["padj"] <= padj_max)
            & (np.abs(comp_df["delta_E1"]) >= delta_e1_min)
        )
        comp_df["significant"] = sig_e1
        if pc1_validated:
            switch_cls = np.where(
                sig_e1 & (comp_df["mean_E1_test"] > 0) & (comp_df["mean_E1_ref"] < 0),
                "B_to_A",
                np.where(
                    sig_e1 & (comp_df["mean_E1_test"] < 0) & (comp_df["mean_E1_ref"] > 0),
                    "A_to_B",
                    np.where(sig_e1, "quantitative_shift", "unchanged"),
                ),
            )
        else:
            warnings.append(
                f"{tag}: unphased_E1_directional_A_B_labels_suppressed (provide genome_fasta or gtf)"
            )
            switch_cls = np.where(
                sig_e1 & ((comp_df["mean_E1_test"] * comp_df["mean_E1_ref"]) < 0),
                "unphased_sign_flip",
                np.where(sig_e1, "quantitative_shift", "unchanged"),
            )
        comp_df["switch_class"] = switch_cls
        comp_csv = out_path / f"differential_compartments_{tag}.csv"
        comp_df.to_csv(comp_csv, index=False)
        ins_out = base_ins_coords.copy()
        ins_out["mean_insulation_test"] = ins_fit["mean_test"]
        ins_out["mean_insulation_ref"] = ins_fit["mean_ref"]
        ins_out["delta_insulation"] = ins_fit["delta"]
        ins_out["t_stat"] = ins_fit["t_stat"]
        ins_out["pvalue"] = ins_fit["pvalue"]
        ins_out["padj"] = ins_fit["padj"]
        sig_ins = (
            np.isfinite(ins_out["padj"])
            & (ins_out["padj"] <= padj_max)
            & (np.abs(ins_out["delta_insulation"]) >= delta_ins_min)
        )
        ins_out["significant"] = sig_ins
        # Lower log2 insulation score = stronger insulation dip at boundary
        ins_out["boundary_change"] = np.where(
            sig_ins & (ins_out["delta_insulation"] < 0),
            "boundary_strengthened",
            np.where(
                sig_ins & (ins_out["delta_insulation"] > 0),
                "boundary_weakened",
                "unchanged",
            ),
        )
        ins_csv = out_path / f"differential_insulation_{tag}.csv"
        ins_out.to_csv(ins_csv, index=False)

        n_a_to_b = int((switch_cls == "A_to_B").sum())
        n_b_to_a = int((switch_cls == "B_to_A").sum())
        n_switches = n_a_to_b + n_b_to_a
        switch_rows = comp_df[comp_df["switch_class"].isin(["A_to_B", "B_to_A"])].sort_values("padj")
        top_switches = [
            {
                "chrom": str(r.chrom),
                "start": int(r.start),
                "end": int(r.end),
                "delta_E1": round(float(r.delta_E1), 4),
                "padj": float(r.padj),
                "switch_class": str(r.switch_class),
            }
            for r in switch_rows.head(15).itertuples(index=False)
        ]

        contrast_results.append({
            "contrast": tag,
            "test": test_grp,
            "reference": ref_grp,
            "status": "success",
            "ran": True,
            "descriptive_only": False,
            "n_test_replicates": n_test,
            "n_ref_replicates": n_ref,
            "n_bins_tested_compartments": int(comp_df["padj"].notna().sum()),
            "n_differential_compartment_bins": int(sig_e1.sum()),
            "n_ab_switches": int(n_switches),
            "n_a_to_b": n_a_to_b,
            "n_b_to_a": n_b_to_a,
            "top_switches": top_switches,
            "differential_compartments_csv": str(comp_csv),
            "n_bins_tested_insulation": int(ins_out["padj"].notna().sum()),
            "n_differential_insulation_bins": int(sig_ins.sum()),
            "n_boundaries_strengthened": int((ins_out["boundary_change"] == "boundary_strengthened").sum()),
            "n_boundaries_weakened": int((ins_out["boundary_change"] == "boundary_weakened").sum()),
            "differential_insulation_csv": str(ins_csv),
        })

    any_ran = any(c.get("ran") for c in contrast_results)
    total_switches = sum(int(c.get("n_ab_switches", 0)) for c in contrast_results if c.get("ran"))
    total_diff_ins = sum(int(c.get("n_differential_insulation_bins", 0)) for c in contrast_results if c.get("ran"))
    primary_comp_csv = next(
        (c.get("differential_compartments_csv") for c in contrast_results if c.get("differential_compartments_csv")),
        None,
    )
    primary_ins_csv = next(
        (c.get("differential_insulation_csv") for c in contrast_results if c.get("differential_insulation_csv")),
        None,
    )

    return {
        "status": "success" if any_ran else "skipped",
        "ran": bool(any_ran),
        "descriptive_only": not any_ran,
        "data_type": "HiC",
        "analysis": "differential_3d",
        "validation_level": "production",
        "method": "empirical_bayes_moderated_ols_bh_fdr",
        "pc1_validated": bool(pc1_validated),
        "phasing_method": phasing_method,
        "padj_max": padj_max,
        "delta_e1_min": delta_e1_min,
        "delta_insulation_min": delta_ins_min,
        "n_ab_switches": int(total_switches),
        "n_differential_insulation_bins": int(total_diff_ins),
        "differential_compartments_csv": primary_comp_csv,
        "differential_insulation_csv": primary_ins_csv,
        "comparisons": contrast_results,
        "warnings": warnings,
    }


# ── Main entrypoint ──────────────────────────────────────────────────────────


def hic_topology(params: dict) -> dict:
    files = [str(f) for f in (params.get("files") or []) if f]
    genome = str(params.get("genome") or "hg38")
    organism = str(params.get("organism") or "Homo sapiens")
    analysis = str(params.get("analysis") or "compartments")
    resolution = int(params.get("resolution") or 100_000)
    chromosomes = params.get("chromosomes", "all")
    allow_mock = mocks_allowed(params)

    valid_files = [f for f in files if Path(f).exists()]
    if not valid_files:
        return {
            "status": "error",
            "error_type": "FileNotFound",
            "details": "No valid Hi-C files found for topology analysis.",
        }

    default_out = Path(valid_files[0]).parent / "aria_hic_topology" / analysis
    out_dir = str(params.get("output_dir") or default_out)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    manifest_params = {
        "analysis": analysis,
        "genome": genome,
        "organism": organism,
        "resolution": resolution,
        "comp_resolution": params.get("comp_resolution"),
        "tad_resolution": params.get("tad_resolution"),
        "window_size": params.get("window_size"),
        "windows": params.get("windows"),
        "chromosomes": chromosomes,
        "sample_ids": params.get("sample_ids"),
        "sample_metadata": (
            str(params.get("sample_metadata"))
            if isinstance(params.get("sample_metadata"), (str, Path))
            else params.get("sample_metadata")
        ),
        "condition_col": params.get("condition_col"),
        "replicate_col": params.get("replicate_col"),
        "covariates": params.get("covariates"),
        "comparisons": params.get("comparisons"),
        "min_replicates_per_condition": params.get("min_replicates_per_condition"),
        "padj_max": params.get("padj_max"),
        "delta_e1_min": params.get("delta_e1_min"),
        "delta_insulation_min": params.get("delta_insulation_min"),
        "fdr_threshold": params.get("fdr_threshold"),
        "max_loci_separation": params.get("max_loci_separation"),
        "genome_fasta": str(params.get("genome_fasta") or ""),
        "gtf": str(params.get("gtf") or params.get("gtf_path") or ""),
    }
    cache_inputs = list(valid_files)
    for aux_key in ("genome_fasta", "gtf", "gtf_path", "sample_metadata"):
        aux_val = params.get(aux_key)
        if isinstance(aux_val, (str, Path)) and Path(str(aux_val)).is_file():
            cache_inputs.append(str(aux_val))
    if bool(params.get("resume", True)) and not allow_mock:
        cached = _check_hic_stage_cache(
            out_path,
            stage=f"hic_topology_{analysis}",
            version=STAGE_VERSION,
            params=manifest_params,
            inputs=cache_inputs,
            manifest_name=f"stage_manifest_{analysis}.json",
        )
        if cached is not None and isinstance(cached.get("result"), dict):
            resumed = dict(cached["result"])
            resumed["resumed"] = True
            return resumed

    fpath = valid_files[0]
    genome_fasta = params.get("genome_fasta")
    gtf_path = params.get("gtf") or params.get("gtf_path")

    if analysis == "compartments":
        result = _compute_compartments(
            fpath, genome, organism, resolution, chromosomes, out_dir,
            genome_fasta=genome_fasta, gtf_path=gtf_path, allow_mock=allow_mock,
        )
    elif analysis == "tads":
        window = int(params.get("window_size") or (resolution * 5))
        result = _compute_tads(
            fpath, genome, resolution, window, chromosomes, out_dir,
            allow_mock=allow_mock,
        )
    elif analysis == "loops":
        result = _compute_loops(
            fpath, genome, resolution, chromosomes, out_dir,
            gtf_path=gtf_path,
            fdr_threshold=float(params.get("fdr_threshold", 0.1)),
            max_loci_separation=int(params.get("max_loci_separation", 2_000_000)),
            allow_mock=allow_mock,
        )
    elif analysis == "insulation_calibration":
        windows = [
            int(w) for w in (params.get("windows") or [resolution * 3, resolution * 5, resolution * 10])
        ]
        result = _calibrate_insulation(
            fpath, resolution, windows, chromosomes, out_dir,
            allow_mock=allow_mock,
        )
    elif analysis == "differential_3d":
        result = _compute_differential_3d(params, out_dir)
    else:
        return {
            "status": "error",
            "error_type": "UnknownAnalysis",
            "details": f"Unknown Hi-C topology analysis: {analysis}",
        }

    if result.get("status") == "success" and not allow_mock:
        outs = [
            str(result[k]) for k in (
                "compartments_tsv", "compartments_bedgraph",
                "tads_bed", "insulation_tsv", "loops_path",
                "annotated_loops_path", "differential_compartments_csv",
                "differential_insulation_csv",
            )
            if result.get(k) and Path(str(result[k])).exists()
        ]
        _save_hic_stage_cache(
            out_path,
            stage=f"hic_topology_{analysis}",
            version=STAGE_VERSION,
            params=manifest_params,
            inputs=cache_inputs,
            outputs=outs,
            extra={"result": result},
            manifest_name=f"stage_manifest_{analysis}.json",
        )
    return result


# ── Explicit mock helpers (only when ARIA_ALLOW_MOCKS=1) ─────────────────────


def _mock_compartments(resolution: int, reason: str) -> dict:
    return {
        "status": "success",
        "pct_A": 45.2,
        "pct_B": 54.8,
        "n_ab_switches": 0,
        "pc1_validated": False,
        "compartment_tracks": {},
        "resolution": resolution,
        "validation_level": "scaffold",
        "note": f"Mock compartments — ({reason})",
    }


def _mock_tads(resolution: int, window: int, reason: str) -> dict:
    return {
        "status": "success",
        "n_tads": 3200,
        "median_size_kb": 185.0,
        "tad_boundaries": {},
        "window_size": window,
        "resolution": resolution,
        "validation_level": "scaffold",
        "algorithm": "Insulation Score (mock)",
        "note": f"Mock TADs — ({reason})",
    }


def _mock_loops(resolution: int, reason: str) -> dict:
    return {
        "status": "success",
        "n_loops": 8500,
        "loops_path": None,
        "algorithm": "mock",
        "resolution": resolution,
        "validation_level": "scaffold",
        "note": f"Mock loops — ({reason})",
    }


if __name__ == "__main__":
    run_script(hic_topology)
