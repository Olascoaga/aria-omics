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
from aria.utils.hic_constants import list_mcool_resolutions

STAGE_VERSION = 1
EXCLUDED_CHROMS = frozenset({"chrM", "MT", "chrMT", "M", "chrEBV", "EBV"})
_ATTR_GENE_NAME_RE = re.compile(r'(?:gene_name|Name)[=\s]+"?([^";]+)"?')
_ATTR_GENE_ID_RE = re.compile(r'(?:gene_id|ID)[=\s]+"?([^";]+)"?')


# ── Shared helpers ───────────────────────────────────────────────────────────


def _snap_window(window_bp: int, resolution_bp: int) -> int:
    """Snap insulation window to an exact integer multiple of ``resolution_bp`` (>= 3 bins)."""
    res = max(int(resolution_bp), 1)
    n_bins = max(3, int(round(int(window_bp) / float(res))))
    return n_bins * res


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


def _safe_ignore_diags(clr, view_df, actual_res: int) -> int:
    import numpy as np

    if view_df is None or view_df.empty:
        return 2
    bins_per_chr = [
        int(clr.chromsizes[c]) // max(int(actual_res), 1)
        for c in view_df["chrom"]
        if c in clr.chromsizes
    ]
    if not bins_per_chr:
        return 2
    med_bins = float(np.median(bins_per_chr))
    return 2 if med_bins >= 12.0 else 1


def _resolve_uri(fpath: str, resolution: int) -> tuple[str, int]:
    import cooler

    if str(fpath).lower().endswith(".mcool"):
        res_map = list_mcool_resolutions(fpath)
        if not res_map:
            raise ValueError(f"No resolutions inside .mcool: {fpath}")
        avail = sorted(res_map.keys())
        chosen = int(resolution) if int(resolution) in res_map else min(
            avail, key=lambda r: abs(r - int(resolution))
        )
        return f"{fpath}::{res_map[chosen]}", chosen
    clr = cooler.Cooler(fpath)
    return fpath, int(clr.binsize)


def _coords_match(df_a, df_b) -> bool:
    """Compare genomic bin coordinates without pandas Categorical category-set pitfalls."""
    import numpy as np

    if len(df_a) != len(df_b):
        return False
    return bool(
        np.array_equal(
            df_a["chrom"].astype(str).to_numpy(),
            df_b["chrom"].astype(str).to_numpy(),
        )
        and np.array_equal(
            df_a["start"].to_numpy(dtype=np.int64),
            df_b["start"].to_numpy(dtype=np.int64),
        )
        and np.array_equal(
            df_a["end"].to_numpy(dtype=np.int64),
            df_b["end"].to_numpy(dtype=np.int64),
        )
    )


def _has_valid_weights(clr) -> bool:
    import numpy as np

    if "weight" not in list(clr.bins().columns):
        return False
    w = np.asarray(clr.bins()[:]["weight"], dtype=float)
    finite = w[np.isfinite(w)]
    return bool(finite.size > 0 and float(np.max(finite)) > 0.0)


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
    s_uri, _ = _resolve_uri(str(staged), chosen_res)
    s_clr = cooler.Cooler(s_uri)
    if not _has_valid_weights(s_clr):
        nbins = int(s_clr.info.get("nbins", 0))
        min_nnz = 10 if nbins >= 500 else max(1, min(3, nbins // 10))
        ignore_diags = 2 if nbins >= 20 else 1
        bal_ret = cooler.balance_cooler(
            s_clr,
            cis_only=False,
            ignore_diags=ignore_diags,
            min_nnz=min_nnz,
            mad_max=5,
            max_iters=300,
            store=True,
        )
        stats = bal_ret[1] if isinstance(bal_ret, tuple) and len(bal_ret) >= 2 else None
        s_clr = cooler.Cooler(s_uri)
        if not _has_valid_weights(s_clr):
            raise RuntimeError(
                f"ICE balancing produced zero valid finite weights for {fpath} at {chosen_res:,}bp."
            )
        if isinstance(stats, dict) and stats.get("converged") is False:
            raise RuntimeError(
                f"ICE balancing did not converge for {fpath} at {chosen_res:,}bp."
            )
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
    if isinstance(chromosomes, str):
        req_list = [chromosomes]
    elif isinstance(chromosomes, (list, tuple, set)):
        req_list = [str(c) for c in chromosomes]
    else:
        return auto_chroms
    req_aliases = {alias for c in req_list for alias in _chrom_aliases(c)}
    # Preserve cooler's native chromosome ordering so cooltools viewframe sorting checks pass.
    picked = [c for c in all_chroms if c in req_aliases]
    return picked or auto_chroms[:1]


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


def _phasing_valid_chroms(phasing_df, chroms: list[str], clr=None) -> set[str]:
    """Return chromosomes whose phasing track has finite non-zero variance over valid ICE bins."""
    import numpy as np

    if phasing_df is None or phasing_df.empty:
        return set()
    val_col = [c for c in phasing_df.columns if c not in {"chrom", "start", "end"}]
    if not val_col:
        return set()
    vcol = val_col[0]

    ice_valid_by_chrom: dict[str, np.ndarray] | None = None
    ice_valid_full: np.ndarray | None = None
    if clr is not None and "weight" in list(clr.bins().columns):
        try:
            bins_w = clr.bins()[:][["chrom", "start", "end", "weight"]]
            w_arr = np.asarray(bins_w["weight"], dtype=float)
            w_ok = np.isfinite(w_arr) & (w_arr > 0.0)
            if len(w_ok) == len(phasing_df) and _coords_match(bins_w, phasing_df):
                ice_valid_full = w_ok
            else:
                ice_valid_by_chrom = {}
                b_ch = bins_w["chrom"].astype(str).to_numpy()
                b_st = bins_w["start"].to_numpy(dtype=np.int64)
                for ch in chroms:
                    cm = b_ch == str(ch)
                    ice_valid_by_chrom[str(ch)] = b_st[cm & w_ok]
        except (OSError, ValueError, KeyError, TypeError):
            return set()

    valid: set[str] = set()
    p_ch = phasing_df["chrom"].astype(str).to_numpy()
    p_st = phasing_df["start"].to_numpy(dtype=np.int64)
    p_val = phasing_df[vcol].to_numpy(dtype=float)
    for ch in chroms:
        ch_s = str(ch)
        mask = p_ch == ch_s
        if ice_valid_full is not None:
            mask = mask & ice_valid_full
        elif ice_valid_by_chrom is not None:
            valid_starts = ice_valid_by_chrom.get(ch_s)
            if valid_starts is None or valid_starts.size == 0:
                continue
            mask = mask & np.isin(p_st, valid_starts)
        if not np.any(mask):
            continue
        vals = p_val[mask]
        fin = vals[np.isfinite(vals)]
        if fin.size >= 5 and float(np.std(fin)) > 1e-12:
            valid.add(ch_s)
    return valid


def _run_eigs_per_chrom_safe(
    eig_mod,
    clr,
    *,
    phasing_df,
    view_df,
    n_eigs: int,
    ignore_diags: int,
    seed: int = 0,
):
    """Run deterministic eigs_cis without wiping chromosomes that lack phasing variance."""
    import pandas as pd

    chroms = [str(c) for c in view_df["chrom"].tolist()]
    if phasing_df is None or phasing_df.empty:
        return _run_eigs_cis_deterministic(
            eig_mod, clr, phasing_track=None, view_df=view_df,
            n_eigs=n_eigs, ignore_diags=ignore_diags, seed=seed,
        )

    valid_phased = _phasing_valid_chroms(phasing_df, chroms, clr=clr)
    if len(valid_phased) == len(chroms):
        return _run_eigs_cis_deterministic(
            eig_mod, clr, phasing_track=phasing_df, view_df=view_df,
            n_eigs=n_eigs, ignore_diags=ignore_diags, seed=seed,
        )
    if not valid_phased:
        return _run_eigs_cis_deterministic(
            eig_mod, clr, phasing_track=None, view_df=view_df,
            n_eigs=n_eigs, ignore_diags=ignore_diags, seed=seed,
        )

    phased_vdf = view_df[view_df["chrom"].astype(str).isin(valid_phased)].reset_index(drop=True)
    unphased_vdf = view_df[~view_df["chrom"].astype(str).isin(valid_phased)].reset_index(drop=True)
    ev_parts = []
    vec_parts = []
    if not phased_vdf.empty:
        ev_p, vec_p = _run_eigs_cis_deterministic(
            eig_mod, clr, phasing_track=phasing_df, view_df=phased_vdf,
            n_eigs=n_eigs, ignore_diags=ignore_diags, seed=seed,
        )
        ev_parts.append(ev_p)
        vec_parts.append(vec_p)
    if not unphased_vdf.empty:
        ev_u, vec_u = _run_eigs_cis_deterministic(
            eig_mod, clr, phasing_track=None, view_df=unphased_vdf,
            n_eigs=n_eigs, ignore_diags=ignore_diags, seed=seed,
        )
        ev_parts.append(ev_u)
        vec_parts.append(vec_u)
    chrom_rank = {c: i for i, c in enumerate(chroms)}
    ev_out = pd.concat(ev_parts, ignore_index=True)
    vec_out = pd.concat(vec_parts, ignore_index=True)
    vec_out["_ord"] = vec_out["chrom"].astype(str).map(lambda c: chrom_rank.get(c, 999999))
    vec_out = vec_out.sort_values(["_ord", "start"], kind="mergesort").drop(columns=["_ord"]).reset_index(drop=True)
    return ev_out, vec_out


def _build_phasing_track(
    clr,
    chroms: list[str],
    genome_fasta: str | None,
    gtf_path: str | None,
    fallback_reasons: list[str] | None = None,
):
    """Build a ``(chrom, start, end, phasing_score)`` DataFrame for ``eigs_cis``.

    Returns ``(phasing_df | None, phasing_method, pc1_validated, phased_chroms)``.
    If ``fallback_reasons`` is supplied, records any FASTA open/index failure that
    triggered fallback to GTF or unphased mode.
    """
    import numpy as np
    import pandas as pd

    bins_df = clr.bins()[:][["chrom", "start", "end"]].copy()
    bins_df["chrom"] = bins_df["chrom"].astype(str)
    chrom_strs = [str(c) for c in chroms]
    fasta_open_error: str | None = None

    # 1. Local genome FASTA GC fraction (gold standard for cooltools eigs_cis)
    if genome_fasta and Path(str(genome_fasta)).is_file():
        try:
            import pysam
        except ImportError:
            pysam = None
            fasta_open_error = (
                f"FASTA GC phasing unavailable ({Path(str(genome_fasta)).name}: "
                "pysam is not installed)"
            )
        if pysam is not None:
            fa = None
            try:
                fa = pysam.FastaFile(str(genome_fasta))
            except OSError as exc:
                # Valid plain-gzip FASTAs (.fa.gz) and FASTAs in read-only dirs
                # without a pre-built .fai raise OSError at open/index time in
                # htslib; record the reason and fall back to GTF phasing.
                # Errors after a successful open (corrupt sequence during fetch)
                # propagate uncaught below.
                fasta_open_error = (
                    f"FASTA open/index failed ({Path(str(genome_fasta)).name}: "
                    f"{str(exc)[:180]})"
                )
            if fa is not None:
                with fa:
                    refs = set(fa.references)
                    alias_to_ref: dict[str, str] = {}
                    for r in refs:
                        for al in _chrom_aliases(r):
                            alias_to_ref.setdefault(al, r)
                    matched = [(c, alias_to_ref[c]) for c in chrom_strs if c in alias_to_ref]
                    if matched:
                        gc_vals = np.full(len(bins_df), np.nan, dtype=float)
                        for ch, ref_ch in matched:
                            ch_len = int(fa.get_reference_length(ref_ch))
                            seq = fa.fetch(ref_ch, 0, ch_len).upper()
                            ch_idx = np.flatnonzero(bins_df["chrom"].to_numpy() == ch)
                            starts = bins_df["start"].to_numpy(dtype=int)[ch_idx]
                            ends = bins_df["end"].to_numpy(dtype=int)[ch_idx]
                            for idx_pos, s, e in zip(ch_idx, starts, ends, strict=True):
                                sub = seq[max(0, s):min(ch_len, e)]
                                if not sub:
                                    continue
                                atgc = (
                                    sub.count("A") + sub.count("T")
                                    + sub.count("G") + sub.count("C")
                                )
                                if atgc >= 20:
                                    gc_vals[idx_pos] = (sub.count("G") + sub.count("C")) / float(atgc)
                        cand_df = bins_df[["chrom", "start", "end"]].copy()
                        cand_df["GC"] = gc_vals
                        valid_ch = _phasing_valid_chroms(cand_df, chrom_strs, clr=clr)
                        if valid_ch:
                            all_phased = bool(chrom_strs) and (set(chrom_strs) <= valid_ch)
                            return cand_df, "gc_content", all_phased, valid_ch

    # 2. Local GTF TSS / gene density per bin
    tss_records = _parse_gtf_tss(gtf_path)
    if tss_records:
        binsize = int(clr.binsize)
        counts = np.zeros(len(bins_df), dtype=float)
        offsets: dict[str, int] = {}
        lengths: dict[str, int] = {}
        alias_to_cooler_ch: dict[str, str] = {}
        for ch, grp in bins_df.groupby("chrom", sort=False, observed=False):
            ch_s = str(ch)
            offsets[ch_s] = int(grp.index[0])
            lengths[ch_s] = len(grp)
            for al in _chrom_aliases(ch_s):
                alias_to_cooler_ch.setdefault(al, ch_s)
        for ch, tss, _gene in tss_records:
            target_ch = alias_to_cooler_ch.get(str(ch))
            if target_ch is None:
                continue
            b = tss // max(binsize, 1)
            if 0 <= b < lengths[target_ch]:
                counts[offsets[target_ch] + b] += 1.0
        cand_df = bins_df[["chrom", "start", "end"]].copy()
        cand_df["gene_density"] = counts
        valid_ch = _phasing_valid_chroms(cand_df, chrom_strs, clr=clr)
        if valid_ch:
            if fasta_open_error and fallback_reasons is not None:
                fallback_reasons.append(
                    f"{fasta_open_error}; fell back to GTF gene_density phasing"
                )
            all_phased = bool(chrom_strs) and (set(chrom_strs) <= valid_ch)
            return cand_df, "gene_density", all_phased, valid_ch

    if fasta_open_error and fallback_reasons is not None:
        fallback_reasons.append(
            f"{fasta_open_error}; no valid GTF phasing track available"
        )
    return None, "unphased_arbitrary_sign", False, set()


from aria.utils.stats import bh_correct as _bh_fdr


# ── 1. Compartments A/B ──────────────────────────────────────────────────────


def _run_eigs_cis_deterministic(eig_mod, clr, *, phasing_track, view_df, n_eigs: int, ignore_diags: int, seed: int = 0):
    """Run ``cooltools.api.eigdecomp.eigs_cis`` with deterministic ARPACK ``v0``."""
    import numpy as np
    import scipy.sparse.linalg as spla

    rng_state = np.random.get_state()
    orig_eigsh = spla.eigsh
    orig_eigs = spla.eigs

    def _det_eigsh(A, *args, **kwargs):
        if kwargs.get("v0") is None and len(args) < 5:
            n = int(A.shape[0])
            kwargs = dict(kwargs)
            kwargs["v0"] = np.ones(n, dtype=np.float64) / np.sqrt(max(1, n))
        return orig_eigsh(A, *args, **kwargs)

    def _det_eigs(A, *args, **kwargs):
        if kwargs.get("v0") is None and len(args) < 5:
            n = int(A.shape[0])
            kwargs = dict(kwargs)
            kwargs["v0"] = np.ones(n, dtype=np.float64) / np.sqrt(max(1, n))
        return orig_eigs(A, *args, **kwargs)

    np.random.seed(int(seed))
    spla.eigsh = _det_eigsh
    spla.eigs = _det_eigs
    try:
        return eig_mod.eigs_cis(
            clr,
            phasing_track=phasing_track,
            view_df=view_df,
            n_eigs=n_eigs,
            ignore_diags=ignore_diags,
        )
    finally:
        spla.eigsh = orig_eigsh
        spla.eigs = orig_eigs
        np.random.set_state(rng_state)


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

        # Filter out tiny scaffolds that would degrade genome-wide ignore_diags or n_eigs.
        long_chroms = [
            c for c in view_df["chrom"]
            if (int(clr.chromsizes[c]) // max(actual_res, 1)) >= 6
        ]
        if long_chroms and len(long_chroms) < len(view_df):
            view_df = view_df[view_df["chrom"].isin(long_chroms)].reset_index(drop=True)
        analyzed_chroms = [str(c) for c in view_df["chrom"].tolist()]
        phasing_fallback_reasons: list[str] = []
        phasing_df, phasing_method, _, phased_chroms = _build_phasing_track(
            clr, analyzed_chroms, genome_fasta, gtf_path,
            fallback_reasons=phasing_fallback_reasons,
        )
        phasing_fallback_reason = (
            "; ".join(phasing_fallback_reasons) if phasing_fallback_reasons else None
        )
        phased_chroms_analyzed = [c for c in analyzed_chroms if c in phased_chroms]
        unphased_chroms = [c for c in analyzed_chroms if c not in phased_chroms]
        pc1_validated = bool(analyzed_chroms) and (set(analyzed_chroms) <= phased_chroms)
        nbins_min = min(int(clr.chromsizes[c]) // max(actual_res, 1) for c in view_df["chrom"])
        n_eigs = min(3, max(1, nbins_min - 2))
        ignore_diags = _safe_ignore_diags(clr, view_df, actual_res)

        eigvals_df, eigvecs_df = _run_eigs_per_chrom_safe(
            eig,
            clr,
            phasing_df=phasing_df,
            view_df=view_df,
            n_eigs=n_eigs,
            ignore_diags=ignore_diags,
            seed=0,
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
        row_phased = export_df["chrom"].astype(str).isin(phased_chroms).to_numpy(dtype=bool)
        export_df["compartment"] = np.where(
            row_phased & (export_df["E1"] > 0), "A",
            np.where(
                row_phased & (export_df["E1"] < 0), "B",
                np.where(
                    (~row_phased) & (export_df["E1"] > 0), "+E1_unphased",
                    np.where((~row_phased) & (export_df["E1"] < 0), "-E1_unphased", "unassigned"),
                ),
            ),
        )
        comp_tsv = out_path / "compartments.tsv"
        export_df.to_csv(comp_tsv, sep="\t", index=False)

        bg_df = export_df[np.isfinite(export_df["E1"])][["chrom", "start", "end", "E1"]]
        comp_bg = out_path / "compartments_E1.bedGraph"
        bg_df.to_csv(comp_bg, sep="\t", header=False, index=False)

        warnings: list[str] = []
        if phasing_fallback_reason:
            warnings.append(phasing_fallback_reason)
        if pc1_validated:
            note = f"E1 phased via {phasing_method} (positive E1 = active A compartment)."
        elif phased_chroms_analyzed:
            note = (
                f"E1 phased via {phasing_method} on {len(phased_chroms_analyzed)}/{len(analyzed_chroms)} "
                f"chromosome(s); unphased chromosomes ({', '.join(unphased_chroms[:5])}) retain arbitrary E1 sign."
            )
            warnings.append("unphased_E1_sign_arbitrary")
            warnings.append(f"partial_unphased_chroms:{','.join(unphased_chroms)}")
        else:
            note = (
                "PC1 sign is unphased (neither genome_fasta nor GTF staged); "
                "provide genome_fasta or gtf to validate A/B polarity."
            )
            warnings.append("unphased_E1_sign_arbitrary")

        comp_result = {
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
            "phased_chroms": phased_chroms_analyzed,
            "unphased_chroms": unphased_chroms,
            "compartments_tsv": str(comp_tsv),
            "compartments_bedgraph": str(comp_bg),
            "compartment_tracks": compartment_tracks,
            "resolution": int(actual_res),
            "warnings": warnings,
            "note": note,
        }
        if phasing_fallback_reason:
            comp_result["phasing_fallback_reason"] = phasing_fallback_reason
        return comp_result
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
        window_bp = _snap_window(window_size, actual_res)
        chroms = _get_chromosomes(clr, chromosomes)
        view_df = _build_view_df(clr, chroms)
        if view_df.empty:
            return {
                "status": "error",
                "error_type": "NoValidChromosomes",
                "details": "No analyzable chromosomes found for TAD calling.",
            }

        ignore_diags = _safe_ignore_diags(clr, view_df, actual_res)
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

        if (
            score_col not in ins_table.columns
            or int(np.isfinite(ins_table[score_col].to_numpy(dtype=float)).sum()) == 0
        ):
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
            st_map: dict[int, float] = {}
            if strength_col in b_df.columns:
                raw_st = b_df[strength_col].to_numpy(dtype=float)
                for pos_i, sval in zip(positions, raw_st, strict=True):
                    if math.isfinite(float(sval)):
                        fv = float(sval)
                        strengths.append(fv)
                        st_map[pos_i] = fv
                    else:
                        st_map[pos_i] = 0.0
            if len(positions) >= 2:
                for idx in range(len(positions) - 1):
                    s_bp = positions[idx]
                    e_bp = positions[idx + 1]
                    if e_bp > s_bp:
                        size_kb = (e_bp - s_bp) / 1000.0
                        tad_sizes_kb.append(size_kb)
                        mean_st = 0.5 * (st_map.get(s_bp, 0.0) + st_map.get(e_bp, 0.0))
                        if not math.isfinite(mean_st):
                            mean_st = 0.0
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
        if not chroms:
            return {
                "status": "error",
                "error_type": "NoValidChromosomes",
                "details": "No analyzable chromosomes found for insulation calibration.",
            }
        cal_chrom = "chr1" if "chr1" in chroms else chroms[0]
        view_df = _build_view_df(clr, [cal_chrom])
        if view_df.empty:
            return {
                "status": "error",
                "error_type": "NoValidChromosomes",
                "details": f"Calibration chromosome {cal_chrom!r} has zero valid length.",
            }
        valid_windows = [_snap_window(w, actual_res) for w in windows]
        ignore_diags = _safe_ignore_diags(clr, view_df, actual_res)

        table = ins.insulation(
            clr,
            sorted(set(valid_windows)),
            view_df=view_df,
            ignore_diags=ignore_diags,
            nproc=1,
        )

        boundary_strengths: dict[str, float] = {}
        for orig_w, eff_w in zip(windows, valid_windows, strict=True):
            bound_col = f"is_boundary_{eff_w}"
            st_col = f"boundary_strength_{eff_w}"
            sc_col = f"log2_insulation_score_{eff_w}"
            b_rows = (
                table[table[bound_col].fillna(False).astype(bool)]
                if bound_col in table.columns else table.iloc[0:0]
            )
            val = 0.0
            if not b_rows.empty and st_col in b_rows.columns:
                arr_st = b_rows[st_col].to_numpy(dtype=float)
                fin_st = arr_st[np.isfinite(arr_st)]
                if fin_st.size > 0:
                    val = float(np.mean(fin_st))
            if val == 0.0 and not b_rows.empty and sc_col in b_rows.columns:
                arr_sc = b_rows[sc_col].to_numpy(dtype=float)
                fin_sc = arr_sc[np.isfinite(arr_sc)]
                if fin_sc.size > 0:
                    val = float(np.mean(np.abs(fin_sc)))
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
        for al in _chrom_aliases(ch):
            by_chrom.setdefault(al, []).append((int(tss), gene))

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

        # Check [4_000, 25_000] resolution bounds (and .mcool fallback) BEFORE staging/balancing.
        target_res = int(resolution)
        if str(fpath).lower().endswith(".mcool"):
            res_map = list_mcool_resolutions(fpath)
            if not res_map:
                raise ValueError(f"No resolutions inside .mcool: {fpath}")
            avail = sorted(res_map.keys())
            resolved_res = target_res if target_res in res_map else min(
                avail, key=lambda r: abs(r - target_res)
            )
            if not (4_000 <= resolved_res <= 25_000):
                avail_ok = [r for r in avail if 4_000 <= r <= 25_000]
                if avail_ok:
                    resolved_res = avail_ok[0]
            target_res = resolved_res
        else:
            target_res = int(cooler.Cooler(fpath).binsize)

        if target_res < 4_000:
            return {
                "status": "skipped",
                "ran": False,
                "data_type": "HiC",
                "analysis": "loops",
                "validation_level": "production",
                "resolution": int(target_res),
                "n_loops": 0,
                "reason": "resolution_too_fine_for_default_donut_kernels",
                "message": (
                    f"Loop dot calling with cooltools default donut kernels requires "
                    f"4 kb - 25 kb resolution; available resolution was {target_res:,} bp."
                ),
            }
        if target_res > 25_000:
            return {
                "status": "skipped",
                "ran": False,
                "data_type": "HiC",
                "analysis": "loops",
                "validation_level": "production",
                "resolution": int(target_res),
                "n_loops": 0,
                "reason": "resolution_too_coarse_for_dot_calling",
                "message": (
                    f"Loop dot calling requires <= 25 kb resolution; "
                    f"finest available resolution was {target_res:,} bp."
                ),
            }

        clr, actual_res = _ensure_balanced_cooler(fpath, target_res, out_path)

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
                sort_cols_ann = [
                    c for c in ("la_exp.donut.qval", "count") if c in ann_df.columns
                ]
                asc_ann = [c != "count" for c in sort_cols_ann]
                if sort_cols_ann:
                    ann_df = ann_df.sort_values(
                        sort_cols_ann, ascending=asc_ann, na_position="last"
                    )
                ann_file = out_path / "annotated_loops.tsv"
                ann_df.to_csv(ann_file, sep="\t", index=False)
                annotated_path = str(ann_file)
                dots_df = ann_df

        top_loops: list[dict[str, Any]] = []
        if not dots_df.empty:
            sort_cols = [c for c in ("la_exp.donut.qval", "count") if c in dots_df.columns]
            asc = [c != "count" for c in sort_cols]
            ranked = (
                dots_df.sort_values(sort_cols, ascending=asc, na_position="last")
                if sort_cols
                else dots_df
            )
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
                    try:
                        c_f = float(row["count"])
                        if math.isfinite(c_f):
                            item["count"] = int(round(c_f))
                    except (TypeError, ValueError):
                        pass
                if "la_exp.donut.qval" in row and row["la_exp.donut.qval"] is not None:
                    try:
                        q_f = float(row["la_exp.donut.qval"])
                        if math.isfinite(q_f):
                            item["qval"] = q_f
                    except (TypeError, ValueError):
                        pass
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
        cov_sub = [row for row, kp in zip(covariate_rows, keep_samples, strict=True) if kp]
        cov_df = pd.DataFrame(cov_sub)
        for col in cov_df.columns:
            raw_vals = list(cov_df[col])
            missing_any = False
            for v in raw_vals:
                if v is None:
                    missing_any = True
                    break
                try:
                    if pd.isna(v):
                        missing_any = True
                        break
                except (TypeError, ValueError):
                    pass
                if isinstance(v, str) and not v.strip():
                    missing_any = True
                    break
            if missing_any:
                raise ValueError(
                    f"Covariate {col!r} contains missing/NaN values for contrast "
                    f"{test_group}_vs_{ref_group}; refusing phantom-category encoding."
                )

            # Treat finite numeric covariates as continuous numeric columns (1 df)
            # instead of one-hot encoding unique numbers as categorical dummies.
            parsed_nums: list[float] = []
            all_numeric = True
            for v in raw_vals:
                if isinstance(v, bool):
                    all_numeric = False
                    break
                try:
                    fv = float(v)
                except (TypeError, ValueError):
                    all_numeric = False
                    break
                if not math.isfinite(fv):
                    all_numeric = False
                    break
                parsed_nums.append(fv)

            if all_numeric:
                num_arr = np.asarray(parsed_nums, dtype=float)
                std_val = float(np.std(num_arr))
                if std_val <= 1e-12:
                    dropped_covariates.append(f"{col}:constant_within_contrast")
                    continue
                cov_block = pd.DataFrame(
                    {str(col): (num_arr - float(np.mean(num_arr))) / std_val}
                )
            else:
                vals = pd.Series([str(v).strip() for v in raw_vals])
                if vals.nunique() < 2:
                    dropped_covariates.append(f"{col}:constant_within_contrast")
                    continue
                cov_block = pd.get_dummies(
                    vals, prefix=str(col), drop_first=True, dtype=float
                )
            cand = pd.concat([df_design, cov_block], axis=1)
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

    raw_files = [str(f) for f in (params.get("files") or []) if f]
    missing_files = [f for f in raw_files if not Path(f).exists()]
    if missing_files:
        return {
            "status": "error",
            "error_type": "FileNotFound",
            "details": f"Declared Hi-C input file(s) not found: {missing_files}",
        }
    files = raw_files
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
        try:
            df_meta = pd.read_csv(sample_metadata, sep=sep)
        except pd.errors.EmptyDataError:
            df_meta = pd.DataFrame()
        if not df_meta.empty and len(df_meta.columns) > 0:
            sid_col = next(
                (c for c in ("sample_id", "sample", "id") if c in df_meta.columns),
                df_meta.columns[0],
            )
            for row_dict in df_meta.to_dict(orient="records"):
                sid_val = row_dict.get(sid_col)
                if sid_val is not None and not pd.isna(sid_val):
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
    if declared_sids is not None:
        if not isinstance(declared_sids, list) or len(declared_sids) != len(files):
            return {
                "status": "error",
                "error_type": "SampleIdLengthMismatch",
                "details": (
                    f"Declared sample_ids length "
                    f"({len(declared_sids) if isinstance(declared_sids, list) else type(declared_sids).__name__}) "
                    f"does not match files length ({len(files)})."
                ),
            }
        raw_sids = [str(s) for s in declared_sids]
    else:
        raw_sids = [Path(f).stem.replace(".balanced", "").replace(".mcool", "").replace(".cool", "") for f in files]
    sid_counts: dict[str, int] = {}
    for s in raw_sids:
        sid_counts[s] = sid_counts.get(s, 0) + 1
    seen_sids: set[str] = set()
    sample_ids: list[str] = []
    for idx_s, (s, fpath_s) in enumerate(zip(raw_sids, files, strict=True)):
        cand = s
        if sid_counts.get(s, 0) > 1:
            parent_tag = Path(fpath_s).parent.name.strip()
            cand = f"{parent_tag}_{s}" if parent_tag else f"{s}_{idx_s + 1}"
            if cand in seen_sids:
                cand = f"{s}_{idx_s + 1}"
        while cand in seen_sids:
            cand = f"{cand}_{idx_s + 1}"
        seen_sids.add(cand)
        sample_ids.append(cand)

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
    ref_e1_by_chrom: dict[str, np.ndarray] = {}
    phasing_df = None
    phasing_built = False
    phasing_coords = None
    phased_chroms: set[str] = set()
    pc1_validated = False
    phasing_method = "unphased_arbitrary_sign"
    phasing_fallback_reasons: list[str] = []

    for sid, fpath in zip(sample_ids, files, strict=True):
        if sid not in meta_by_sid:
            continue
        clr_c, actual_comp_res = _ensure_balanced_cooler(fpath, comp_res, out_path / sid)
        chroms = _get_chromosomes(clr_c, chromosomes)
        vdf_c = _build_view_df(clr_c, chroms)
        if vdf_c.empty:
            return {
                "status": "error",
                "error_type": "NoValidChromosomes",
                "details": f"Sample {sid!r} has no analyzable chromosomes for compartment testing.",
            }
        long_c = [
            c for c in vdf_c["chrom"]
            if (int(clr_c.chromsizes[c]) // max(actual_comp_res, 1)) >= 6
        ]
        if long_c and len(long_c) < len(vdf_c):
            vdf_c = vdf_c[vdf_c["chrom"].isin(long_c)].reset_index(drop=True)
        analyzed_c = [str(c) for c in vdf_c["chrom"].tolist()]
        vdf_chrom_set = set(analyzed_c)
        all_bins_c = clr_c.bins()[:][["chrom", "start", "end"]]
        cur_bins = all_bins_c[all_bins_c["chrom"].astype(str).isin(vdf_chrom_set)].reset_index(drop=True)
        if not phasing_built:
            phasing_df, phasing_method, _, phased_chroms = _build_phasing_track(
                clr_c, analyzed_c, genome_fasta, gtf_path,
                fallback_reasons=phasing_fallback_reasons,
            )
            phasing_coords = cur_bins
            phasing_built = True
        elif phasing_coords is not None and not _coords_match(cur_bins, phasing_coords):
            return {
                "status": "error",
                "error_type": "CoordinateMismatch",
                "details": (
                    f"Sample {sid!r} compartment bin grid does not match cohort reference "
                    "(mismatched resolution or chromosome assembly across Hi-C inputs)."
                ),
            }
        else:
            sample_valid_phased = _phasing_valid_chroms(phasing_df, analyzed_c, clr=clr_c)
            phased_chroms &= sample_valid_phased

        nb_min = min(int(clr_c.chromsizes[c]) // max(actual_comp_res, 1) for c in vdf_c["chrom"])
        # Pass phasing track restricted to the analyzed chromosomes so extra contigs do not misalign.
        phasing_for_sample = (
            phasing_df[phasing_df["chrom"].astype(str).isin(vdf_chrom_set)].reset_index(drop=True)
            if phasing_df is not None
            else None
        )
        _, ev_df = _run_eigs_per_chrom_safe(
            eig,
            clr_c,
            phasing_df=phasing_for_sample,
            view_df=vdf_c,
            n_eigs=min(3, max(1, nb_min - 2)),
            ignore_diags=_safe_ignore_diags(clr_c, vdf_c, actual_comp_res),
            seed=0,
        )
        # Always align cross-replicate E1 sign per chromosome for chromosomes NOT in
        # phased_chroms using the first sample with >=5 finite bins on that chromosome.
        unphased_for_align = [c for c in analyzed_c if c not in phased_chroms]
        if unphased_for_align:
            ev_df = ev_df.copy()
            for ch_name in unphased_for_align:
                c_mask = ev_df["chrom"].astype(str) == ch_name
                cur_ch = ev_df.loc[c_mask, "E1"].to_numpy(dtype=float)
                ref_ch = ref_e1_by_chrom.get(ch_name)
                if ref_ch is None:
                    if int(np.isfinite(cur_ch).sum()) >= 5:
                        ref_e1_by_chrom[ch_name] = cur_ch.copy()
                elif ref_ch.shape == cur_ch.shape:
                    ok = np.isfinite(ref_ch) & np.isfinite(cur_ch)
                    if int(ok.sum()) >= 5 and float(np.dot(ref_ch[ok], cur_ch[ok])) < 0:
                        ev_df.loc[c_mask, "E1"] = -cur_ch
        e1_by_sid[sid] = ev_df[["chrom", "start", "end", "E1"]].reset_index(drop=True)

        clr_t, actual_tad_res = _ensure_balanced_cooler(fpath, tad_res, out_path / sid)
        vdf_t = _build_view_df(clr_t, _get_chromosomes(clr_t, chromosomes))
        if vdf_t.empty:
            return {
                "status": "error",
                "error_type": "NoValidChromosomes",
                "details": f"Sample {sid!r} has no analyzable chromosomes for insulation testing.",
            }
        eff_w = _snap_window(window_bp, actual_tad_res)
        ins_df = ins.insulation(
            clr_t, [eff_w], view_df=vdf_t, ignore_diags=_safe_ignore_diags(clr_t, vdf_t, actual_tad_res), nproc=1
        )
        sc_col = f"log2_insulation_score_{eff_w}"
        ins_by_sid[sid] = ins_df[["chrom", "start", "end", sc_col]].rename(
            columns={sc_col: "insulation"}
        ).reset_index(drop=True)

    # Collapse technical replicates to biological replicates (mean track per replicate unit)
    def _clean_cov_val(val: Any) -> Any:
        if val is None:
            return None
        try:
            if pd.isna(val):
                return None
        except (TypeError, ValueError):
            pass
        if isinstance(val, bool):
            return str(val)
        if isinstance(val, (int, float, np.integer, np.floating)):
            fv = float(val)
            return fv if math.isfinite(fv) else None
        s = str(val).strip()
        return s if s != "" else None

    rep_units: dict[tuple[str, str], list[str]] = {}
    rep_covs: dict[tuple[str, str], dict[str, Any]] = {}
    for sid in e1_by_sid:
        m = meta_by_sid[sid]
        raw_cond = m.get(cond_col) if cond_col in m else m.get("condition")
        if raw_cond is None or pd.isna(raw_cond) or str(raw_cond).strip() == "":
            continue
        raw_rep = m.get(rep_col)
        if raw_rep is None or pd.isna(raw_rep) or str(raw_rep).strip() == "":
            rep = sid
        else:
            rep = str(raw_rep)
        key = (str(raw_cond), str(rep))
        rep_units.setdefault(key, []).append(sid)
        rep_covs.setdefault(key, {c: _clean_cov_val(m.get(c)) for c in covariates})

    if not e1_by_sid or not rep_units:
        return {
            "status": "skipped",
            "ran": False,
            "data_type": "HiC",
            "analysis": "differential_3d",
            "validation_level": "production",
            "reason": "no_matched_samples_in_metadata",
        }

    contrast_results: list[dict[str, Any]] = []
    warnings: list[str] = []
    if min_reps == 2:
        warnings.append("low_power_warning: min_replicates_per_condition=2")

    base_e1_coords = next(iter(e1_by_sid.values()))[["chrom", "start", "end"]]
    base_ins_coords = next(iter(ins_by_sid.values()))[["chrom", "start", "end"]]
    for sid, df_e1 in e1_by_sid.items():
        if not _coords_match(df_e1, base_e1_coords):
            return {
                "status": "error",
                "error_type": "CoordinateMismatch",
                "details": (
                    f"Sample {sid!r} compartment bin grid does not match cohort reference "
                    "(mismatched resolution or chromosome assembly across Hi-C inputs)."
                ),
            }
    for sid, df_in in ins_by_sid.items():
        if not _coords_match(df_in, base_ins_coords):
            return {
                "status": "error",
                "error_type": "CoordinateMismatch",
                "details": (
                    f"Sample {sid!r} insulation bin grid does not match cohort reference "
                    "(mismatched resolution or chromosome assembly across Hi-C inputs)."
                ),
            }

    analyzed_chroms = list(dict.fromkeys(base_e1_coords["chrom"].astype(str).tolist()))
    phased_chroms_analyzed = [c for c in analyzed_chroms if c in phased_chroms]
    unphased_chroms = [c for c in analyzed_chroms if c not in phased_chroms]
    pc1_validated = bool(analyzed_chroms) and (set(analyzed_chroms) <= phased_chroms)

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
        row_phased = comp_df["chrom"].astype(str).isin(phased_chroms).to_numpy(dtype=bool)
        if not pc1_validated:
            warnings.append(
                f"{tag}: unphased_E1_directional_A_B_labels_suppressed (provide genome_fasta or gtf)"
            )
        switch_cls = np.where(
            sig_e1 & row_phased & (comp_df["mean_E1_test"] > 0) & (comp_df["mean_E1_ref"] < 0),
            "B_to_A",
            np.where(
                sig_e1 & row_phased & (comp_df["mean_E1_test"] < 0) & (comp_df["mean_E1_ref"] > 0),
                "A_to_B",
                np.where(
                    sig_e1 & (~row_phased) & ((comp_df["mean_E1_test"] * comp_df["mean_E1_ref"]) < 0),
                    "unphased_sign_flip",
                    np.where(sig_e1, "quantitative_shift", "unchanged"),
                ),
            ),
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
        n_unphased_flips = int((switch_cls == "unphased_sign_flip").sum())
        switch_classes_to_rank = (
            ["A_to_B", "B_to_A", "unphased_sign_flip"]
            if n_unphased_flips > 0
            else ["A_to_B", "B_to_A"]
        )
        switch_rows = comp_df[comp_df["switch_class"].isin(switch_classes_to_rank)].sort_values("padj")
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
            "n_unphased_sign_flips": n_unphased_flips,
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
    total_unphased_flips = sum(
        int(c.get("n_unphased_sign_flips", 0)) for c in contrast_results if c.get("ran")
    )
    total_diff_ins = sum(int(c.get("n_differential_insulation_bins", 0)) for c in contrast_results if c.get("ran"))
    primary_comp_csv = next(
        (c.get("differential_compartments_csv") for c in contrast_results if c.get("differential_compartments_csv")),
        None,
    )
    primary_ins_csv = next(
        (c.get("differential_insulation_csv") for c in contrast_results if c.get("differential_insulation_csv")),
        None,
    )

    phasing_fallback_reason = (
        "; ".join(phasing_fallback_reasons) if phasing_fallback_reasons else None
    )
    if phasing_fallback_reason:
        warnings.append(phasing_fallback_reason)
    diff_result = {
        "status": "success" if any_ran else "skipped",
        "ran": bool(any_ran),
        "descriptive_only": not any_ran,
        "data_type": "HiC",
        "analysis": "differential_3d",
        "validation_level": "production",
        "method": "empirical_bayes_moderated_ols_bh_fdr",
        "pc1_validated": bool(pc1_validated),
        "phasing_method": phasing_method,
        "phased_chroms": phased_chroms_analyzed,
        "unphased_chroms": unphased_chroms,
        "padj_max": padj_max,
        "delta_e1_min": delta_e1_min,
        "delta_insulation_min": delta_ins_min,
        "n_ab_switches": int(total_switches),
        "n_unphased_sign_flips": int(total_unphased_flips),
        "n_differential_insulation_bins": int(total_diff_ins),
        "differential_compartments_csv": primary_comp_csv,
        "differential_insulation_csv": primary_ins_csv,
        "comparisons": contrast_results,
        "warnings": warnings,
    }
    if phasing_fallback_reason:
        diff_result["phasing_fallback_reason"] = phasing_fallback_reason
    return diff_result


# ── Main entrypoint ──────────────────────────────────────────────────────────


def hic_topology(params: dict) -> dict:
    files = [str(f) for f in (params.get("files") or []) if f]
    genome = str(params.get("genome") or "hg38")
    organism = str(params.get("organism") or "Homo sapiens")
    analysis = str(params.get("analysis") or "compartments")
    resolution = int(params.get("resolution") or 100_000)
    chromosomes = params.get("chromosomes", "all")
    allow_mock = mocks_allowed(params)

    if not files:
        return {
            "status": "error",
            "error_type": "FileNotFound",
            "details": "No valid Hi-C files found for topology analysis.",
        }
    missing_files = [f for f in files if not Path(f).exists()]
    if missing_files:
        return {
            "status": "error",
            "error_type": "FileNotFound",
            "details": f"Declared Hi-C input file(s) not found: {missing_files}",
        }
    declared_sids = params.get("sample_ids")
    if declared_sids is not None:
        if not isinstance(declared_sids, list) or len(declared_sids) != len(files):
            return {
                "status": "error",
                "error_type": "SampleIdLengthMismatch",
                "details": (
                    f"Declared sample_ids length "
                    f"({len(declared_sids) if isinstance(declared_sids, list) else type(declared_sids).__name__}) "
                    f"does not match files length ({len(files)})."
                ),
            }
    valid_files = list(files)

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
        out_candidates = [
            str(result[k]) for k in (
                "compartments_tsv", "compartments_bedgraph",
                "tads_bed", "insulation_tsv", "loops_path",
                "annotated_loops_path", "differential_compartments_csv",
                "differential_insulation_csv", "descriptive_csv",
            )
            if result.get(k)
        ]
        for comp in result.get("comparisons") or []:
            if not isinstance(comp, dict):
                continue
            for ck in (
                "differential_compartments_csv",
                "differential_insulation_csv",
                "descriptive_csv",
            ):
                if comp.get(ck):
                    out_candidates.append(str(comp[ck]))
        outs = [p for p in dict.fromkeys(out_candidates) if Path(p).is_file()]
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
