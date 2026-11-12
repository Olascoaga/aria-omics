"""Unified paired peak-to-gene correlation & background-null engine (V48 Slice I5).

Shared pure-Python / NumPy / SciPy implementation backing BOTH:
  1. `aria/scripts/chromatin_regulatory.py` (`_peak_to_gene`, ADR-050 beta-scoped
     scATAC regulatory sub-layer, preserving `preprocess="as_is"` and
     `background_null=False` for 100% backward compatibility).
  2. `aria/scripts/integration_peak2gene.py` (V48 Slice I5 paired multiome Tier 2b
     regulatory layer, implementing `Signac::LinkPeaks`-faithful `log1p(CP10k)` RNA
     vs raw-count ATAC correlation, deterministic `(GC, accessibility, length)`
     binned trans-chromosomal background z-scores, two-sided + one-sided normal
     tail p-values, and pre-registered Benjamini-Hochberg FDR across all tested
     in-window cis pairs).

Governance invariants (ADR-002, ADR-005, ADR-011, ADR-050, ADR-064):
  - Strictly associative: correlations, z-scores, and distances describe paired
    accessibility-expression co-variation and NEVER establish enhancer, silencer,
    or poised regulatory mechanism without independent orthogonal assays.
  - Zero hardcoded gene lists or dataset-specific heuristics; coordinates are
    parsed strictly as genomic intervals.
  - Uses ARIA's shared Benjamini-Hochberg implementation (`aria.utils.stats.bh_correct`)
    across all tested cis pairs prior to effect-size thresholding.
"""

from __future__ import annotations

import gzip
import math
import os
from pathlib import Path
import re
from typing import Any, Sequence

import numpy as np

from aria.utils.stats import bh_correct

_PEAK_SEP_RE = re.compile(r"^([^:\-_\t]+)[:\-_](\d+)[\-_](\d+)$")
_STANDARD_AUTOSOMES = {str(i) for i in range(1, 23)}


def canonical_contig(chrom: str | None) -> str | None:
    """Normalize UCSC (`chr1`, `chrM`) and Ensembl (`1`, `MT`) contig aliases.

    Standard human/mouse nuclear chromosomes (`1..22`, `X`, `Y`) and mitochondrial
    aliases (`M`, `MT`, `chrM`, `chrMT`) map to canonical `chr*` keys so UCSC
    peaks (`chr1:100-200`) and Ensembl GTFs (`1`) interoperate seamlessly while
    preserving non-standard contig identifiers in normalized form.
    """
    if chrom is None:
        return None
    raw = str(chrom).strip()
    if not raw:
        return None
    low = raw.lower()
    if low in ("m", "mt", "chrm", "chrmt"):
        return "chrM"
    stripped = raw[3:] if low.startswith("chr") else raw
    up = stripped.upper()
    if stripped in _STANDARD_AUTOSOMES or up in ("X", "Y"):
        return f"chr{up}"
    # Preserve mouse autosomes 1..19 (already in 1..22) or rat/other numeric contigs
    if stripped.isdigit() and 1 <= int(stripped) <= 99:
        return f"chr{int(stripped)}"
    return raw


def is_standard_contig(chrom: str | None) -> bool:
    """Return True if `chrom` maps to `chr1..chr22`, `chrX`, or `chrY`."""
    canon = canonical_contig(chrom)
    if canon is None or not canon.startswith("chr"):
        return False
    body = canon[3:]
    return body in _STANDARD_AUTOSOMES or body in ("X", "Y")


def parse_peak_coord(name: str) -> tuple[str, int, int] | None:
    """Parse a peak identifier (`chr1:100-200`, `chr1-100-200`, or BED tab line)."""
    text = str(name).strip()
    if not text:
        return None
    try:
        if "\t" in text:
            parts = text.split("\t")
            if len(parts) >= 3:
                s, e = int(parts[1]), int(parts[2])
                return parts[0].strip(), min(s, e), max(s, e)
        if ":" in text and "-" in text:
            chrom, pos = text.split(":", 1)
            s_txt, e_txt = pos.replace("_", "-").split("-", 1)
            s, e = int(s_txt), int(e_txt)
            return chrom.strip(), min(s, e), max(s, e)
        m = _PEAK_SEP_RE.match(text)
        if m:
            s, e = int(m.group(2)), int(m.group(3))
            return m.group(1).strip(), min(s, e), max(s, e)
    except Exception:
        return None
    return None


def _gtf_attr(attrs: str, key: str) -> str | None:
    prefix_eq = f"{key}="
    prefix_sp = f"{key} "
    for chunk in attrs.split(";"):
        item = chunk.strip()
        if not item:
            continue
        if item.startswith(prefix_sp):
            val = item[len(prefix_sp) :].strip().strip('"')
            return val or None
        if item.startswith(prefix_eq):
            val = item[len(prefix_eq) :].strip().strip('"')
            return val or None
    return None


def _gtf_has_tag(attrs: str, tag_value: str) -> bool:
    target = tag_value.lower()
    for chunk in attrs.split(";"):
        item = chunk.strip()
        if item.lower().startswith("tag ") or item.lower().startswith("tag="):
            val = item.split(None, 1)[-1] if " " in item else item.split("=", 1)[-1]
            if val.strip().strip('"').lower() == target:
                return True
    return False


def parse_gtf_tss_detailed(
    gtf_path: str | None,
    wanted: set[str] | None = None,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    """Parse canonical per-gene TSS records from a local `.gtf` or `.gtf.gz`.

    Resolves duplicate `gene_name` entries deterministically:
      1. Drops pseudoautosomal Y duplicates (`_PAR_Y` / `PAR` tag on Y).
      2. Prefers `MANE_Select` / `Ensembl_canonical` tags.
      3. Prefers `protein_coding` biotype over pseudogenes/fragments.
      4. Prefers `basic` tag.
      5. Prefers longest genomic gene span (`end - start`).
      6. Breaks remaining ties by lexicographically smallest `gene_id`.
    """
    stats: dict[str, Any] = {
        "n_gene_rows_scanned": 0,
        "n_par_y_dropped": 0,
        "n_duplicate_symbols_resolved": 0,
        "n_unique_genes": 0,
    }
    if not gtf_path or not Path(str(gtf_path)).is_file():
        return {}, stats

    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    best_by_gene: dict[str, tuple[tuple[int, int, int, int], str, dict[str, Any]]] = {}

    def _scan_feature(target_feature: str) -> None:
        needle = f"\t{target_feature}\t"
        with opener(str(gtf_path), "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line or line[0] == "#" or needle not in line:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != target_feature:
                    continue
                stats["n_gene_rows_scanned"] += 1
                attrs = parts[8]
                gene_id = _gtf_attr(attrs, "gene_id") or ""
                gene_name = _gtf_attr(attrs, "gene_name") or gene_id
                if not gene_name:
                    continue
                if wanted is not None and gene_name not in wanted:
                    if not gene_id or gene_id not in wanted:
                        continue
                    key_symbol = gene_id
                else:
                    key_symbol = gene_name

                if gene_id.upper().endswith("_PAR_Y") or (
                    _gtf_has_tag(attrs, "PAR")
                    and canonical_contig(parts[0]) == "chrY"
                ):
                    stats["n_par_y_dropped"] += 1
                    continue

                try:
                    start = int(parts[3])
                    end = int(parts[4])
                except ValueError:
                    continue
                strand = "-" if parts[6] == "-" else "+"
                tss = start if strand == "+" else end
                span = max(0, end - start)

                biotype = (
                    _gtf_attr(attrs, "gene_type")
                    or _gtf_attr(attrs, "gene_biotype")
                    or ""
                ).lower()
                is_mane = int(
                    _gtf_has_tag(attrs, "MANE_Select")
                    or _gtf_has_tag(attrs, "Ensembl_canonical")
                )
                is_pc = int(biotype == "protein_coding")
                is_basic = int(_gtf_has_tag(attrs, "basic"))
                prio = (is_mane, is_pc, is_basic, span)

                rec = {
                    "gene": key_symbol,
                    "gene_id": gene_id or key_symbol,
                    "chrom": parts[0].strip(),
                    "canon_chrom": canonical_contig(parts[0]),
                    "start": start,
                    "end": end,
                    "tss": tss,
                    "strand": strand,
                    "gene_biotype": biotype or None,
                }

                prev = best_by_gene.get(key_symbol)
                if prev is None:
                    best_by_gene[key_symbol] = (prio, rec["gene_id"], rec)
                else:
                    stats["n_duplicate_symbols_resolved"] += 1
                    prev_prio, prev_gid, _ = prev
                    if prio > prev_prio or (
                        prio == prev_prio and rec["gene_id"] < prev_gid
                    ):
                        best_by_gene[key_symbol] = (prio, rec["gene_id"], rec)

    _scan_feature("gene")
    if not best_by_gene:
        _scan_feature("transcript")

    out = {k: v[2] for k, v in best_by_gene.items()}
    stats["n_unique_genes"] = len(out)
    return out, stats


def parse_gtf_tss(
    gtf_path: str | None,
    wanted: set[str] | None = None,
) -> dict[str, tuple[str, int]]:
    """Backward-compatible `{gene: (chrom, tss)}` wrapper over `parse_gtf_tss_detailed`."""
    detailed, _ = parse_gtf_tss_detailed(gtf_path, wanted=wanted)
    return {g: (rec["chrom"], int(rec["tss"])) for g, rec in detailed.items()}


def resolve_gtf_path(
    explicit_gtf: str | None = None,
    genome: str = "hg38",
) -> str | None:
    """Resolve local GTF path deterministically without network access."""
    if explicit_gtf is not None and str(explicit_gtf).strip():
        cand = Path(str(explicit_gtf)).expanduser()
        if not cand.is_file():
            raise FileNotFoundError(f"Explicit GTF file does not exist: {explicit_gtf}")
        return str(cand)

    env_gtf = os.environ.get("ARIA_GTF_PATH", "").strip()
    if env_gtf:
        env_p = Path(env_gtf).expanduser()
        if env_p.is_file():
            return str(env_p)

    genome_dir = Path.home() / ".aria" / "genomes" / str(genome)
    for fname in (
        "annotation.gtf",
        "annotation.gtf.gz",
        f"{genome}.gtf",
        f"{genome}.gtf.gz",
    ):
        cand = genome_dir / fname
        if cand.is_file():
            return str(cand)
    return None


def resolve_indexed_fasta_path(
    explicit_fasta: str | None = None,
    genome: str = "hg38",
) -> str | None:
    """Resolve an indexed local genome FASTA (`.fa` with `.fa.fai`) when present."""
    if explicit_fasta is not None and str(explicit_fasta).strip():
        cand = Path(str(explicit_fasta)).expanduser()
        if not cand.is_file():
            raise FileNotFoundError(
                f"Explicit genome_fasta does not exist: {explicit_fasta}"
            )
        fai = Path(f"{cand}.fai")
        return str(cand) if fai.is_file() else None

    genome_dir = Path.home() / ".aria" / "genomes" / str(genome)
    for fname in (
        "genome.fa",
        "GRCh38.primary_assembly.genome.fa",
        f"{genome}.fa",
    ):
        cand = genome_dir / fname
        if cand.is_file() and Path(f"{cand}.fai").is_file():
            return str(cand)
    return None


def compute_peak_gc_from_indexed_fasta(
    fasta_path: str | None,
    peak_coords: Sequence[tuple[str, int, int] | None],
) -> tuple[np.ndarray | None, dict[str, Any]]:
    """Compute exact per-peak GC fraction `[0, 1]` via `.fai` byte-offset seeking."""
    if not fasta_path:
        return None, {"available": False, "reason": "no_genome_fasta"}
    fa_p = Path(str(fasta_path))
    fai_p = Path(f"{fa_p}.fai")
    if not fa_p.is_file() or not fai_p.is_file():
        return None, {"available": False, "reason": "missing_fasta_or_fai"}

    # Parse FASTA index: contig -> (seq_len, byte_offset, line_bases, line_bytes)
    fai_by_canon: dict[str, tuple[int, int, int, int]] = {}
    fai_exact: dict[str, tuple[int, int, int, int]] = {}
    with open(fai_p, "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            cname = parts[0].strip()
            entry = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))
            fai_exact[cname] = entry
            canon = canonical_contig(cname)
            if canon is not None and canon not in fai_by_canon:
                fai_by_canon[canon] = entry

    n_peaks = len(peak_coords)
    gc = np.full(n_peaks, np.nan, dtype=np.float64)
    n_computed = 0

    with open(fa_p, "rb") as fh:
        for i, coord in enumerate(peak_coords):
            if coord is None:
                continue
            chrom, start, end = coord
            rec = fai_exact.get(chrom) or fai_by_canon.get(canonical_contig(chrom) or "")
            if rec is None:
                continue
            seq_len, offset, line_bases, line_bytes = rec
            # Support both 0-based and 1-based start safely clamped to [0, seq_len).
            s0 = max(0, min(seq_len, int(start) - 1 if int(start) >= 1 else 0))
            e0 = max(s0, min(seq_len, int(end)))
            if e0 <= s0 or line_bases <= 0:
                continue

            start_line, start_col = divmod(s0, line_bases)
            end_line, end_col = divmod(e0, line_bases)
            byte_start = offset + start_line * line_bytes + start_col
            byte_end = offset + end_line * line_bytes + end_col
            fh.seek(byte_start)
            raw_bytes = fh.read(max(0, byte_end - byte_start))
            seq = raw_bytes.translate(None, b"\r\n\t ").upper()
            a_cnt = seq.count(b"A")
            c_cnt = seq.count(b"C")
            g_cnt = seq.count(b"G")
            t_cnt = seq.count(b"T")
            atgc = a_cnt + c_cnt + g_cnt + t_cnt
            if atgc > 0:
                gc[i] = float(c_cnt + g_cnt) / float(atgc)
                n_computed += 1

    valid_coord_count = sum(1 for c in peak_coords if c is not None)
    if valid_coord_count == 0 or n_computed < max(1, int(0.5 * valid_coord_count)):
        return None, {
            "available": False,
            "reason": "insufficient_fasta_contig_overlap",
            "n_computed": int(n_computed),
            "n_valid_peaks": int(valid_coord_count),
        }

    med_gc = float(np.nanmedian(gc))
    gc = np.where(np.isfinite(gc), gc, med_gc)
    return gc, {
        "available": True,
        "n_computed": int(n_computed),
        "n_valid_peaks": int(valid_coord_count),
        "median_gc": round(med_gc, 4),
    }


def genes_to_peak_indices(
    peak_names: Sequence[Any],
    gene_tss: dict[str, tuple[str, int]],
    *,
    window_bp: int,
) -> dict[str, list[int]]:
    """Map each gene to in-window peak indices using canonical contig matching."""
    peaks_by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for idx, peak in enumerate(peak_names):
        coord = parse_peak_coord(str(peak))
        if coord is None:
            continue
        chrom, start, end = coord
        ckey = canonical_contig(chrom) or chrom
        peaks_by_chrom.setdefault(ckey, []).append((start, end, idx))

    out: dict[str, list[int]] = {}
    for gene, (chrom, tss) in gene_tss.items():
        ckey = canonical_contig(chrom) or chrom
        hits = []
        lo, hi = int(tss) - int(window_bp), int(tss) + int(window_bp)
        for start, end, idx in peaks_by_chrom.get(ckey, []):
            if end >= lo and start <= hi:
                hits.append(idx)
        if hits:
            out[gene] = hits
    return out


def scalar_pearson_corr(x: Any, y: Any) -> float | None:
    """Scalar Pearson correlation rounded to 6 decimals; None when undefined."""
    a = np.asarray(x, dtype=np.float64).ravel()
    b = np.asarray(y, dtype=np.float64).ravel()
    if a.shape[0] != b.shape[0] or a.shape[0] < 3:
        return None
    if float(np.std(a)) == 0.0 or float(np.std(b)) == 0.0:
        return None
    r = float(np.corrcoef(a, b)[0, 1])
    return round(r, 6) if math.isfinite(r) else None


def vectorized_pearson_corr(x: Any, cols: Any) -> np.ndarray:
    """Pearson correlation of 1D vector `x` `(n,)` against each column of `cols` `(n, k)`."""
    xv = np.asarray(x, dtype=np.float64).ravel()
    n = xv.shape[0]
    k = cols.shape[1] if getattr(cols, "ndim", 1) == 2 else 1
    if n < 3:
        return np.full(k, np.nan, dtype=np.float64)
    xc = xv - xv.mean()
    xs = float(np.sqrt(np.dot(xc, xc)))
    c_mat = np.asarray(cols, dtype=np.float64)
    if c_mat.ndim == 1:
        c_mat = c_mat.reshape(-1, 1)
    cc = c_mat - c_mat.mean(axis=0, keepdims=True)
    cs = np.sqrt(np.sum(cc * cc, axis=0))
    denom = xs * cs
    num = xc @ cc
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(denom > 0.0, num / denom, np.nan)
    return np.clip(r, -1.0, 1.0)


def _quantile_bin_labels(values: np.ndarray, n_bins: int) -> np.ndarray:
    """Assign finite 1D values into `0..n_bins-1` rank-quantile bins deterministically."""
    n = values.size
    if n == 0 or n_bins <= 1:
        return np.zeros(n, dtype=int)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(n, dtype=np.float64)
    ranks[order] = np.arange(n, dtype=np.float64)
    bins = np.floor(ranks * float(n_bins) / float(n)).astype(int)
    return np.clip(bins, 0, n_bins - 1)


def _log1p_cp10k_csc(mat: Any, target_sum: float = 1e4):
    """Normalize cell-by-gene raw counts to `log1p(target_sum * x / cell_sum)` in CSC."""
    import scipy.sparse as sp

    if sp.issparse(mat):
        csr = sp.csr_matrix(mat, dtype=np.float64, copy=True)
        row_sums = np.asarray(csr.sum(axis=1)).ravel()
        scale = np.zeros_like(row_sums, dtype=np.float64)
        nz = row_sums > 0.0
        scale[nz] = float(target_sum) / row_sums[nz]
        csr = sp.diags(scale) @ csr
        csr.data = np.log1p(csr.data)
        return csr.tocsc()

    arr = np.asarray(mat, dtype=np.float64)
    row_sums = np.sum(arr, axis=1, keepdims=True)
    scale = np.where(row_sums > 0.0, float(target_sum) / np.maximum(row_sums, 1e-12), 0.0)
    return np.log1p(arr * scale)


def _signac_tfidf_csc(mat: Any, scale_factor: float = 1e4):
    """Signac `RunTFIDF(method=1)`: `log1p(scale_factor * (x / cell_sum) * (n_cells / peak_sum))`."""
    import scipy.sparse as sp

    n_cells = mat.shape[0]
    if sp.issparse(mat):
        csr = sp.csr_matrix(mat, dtype=np.float64, copy=True)
        cell_sums = np.asarray(csr.sum(axis=1)).ravel()
        peak_sums = np.asarray(csr.sum(axis=0)).ravel()
        inv_cell = np.zeros_like(cell_sums)
        nz_c = cell_sums > 0.0
        inv_cell[nz_c] = 1.0 / cell_sums[nz_c]
        idf = np.zeros_like(peak_sums)
        nz_p = peak_sums > 0.0
        idf[nz_p] = float(n_cells) / peak_sums[nz_p]
        scaled = sp.diags(inv_cell * float(scale_factor)) @ csr @ sp.diags(idf)
        scaled.data = np.log1p(scaled.data)
        return scaled.tocsc()

    arr = np.asarray(mat, dtype=np.float64)
    cell_sums = np.sum(arr, axis=1, keepdims=True)
    peak_sums = np.sum(arr, axis=0, keepdims=True)
    tf = np.where(cell_sums > 0.0, arr / np.maximum(cell_sums, 1e-12), 0.0)
    idf = np.where(peak_sums > 0.0, float(n_cells) / np.maximum(peak_sums, 1e-12), 0.0)
    return np.log1p(tf * idf * float(scale_factor))


def _col_dense(csc_or_dense: Any, idx: int) -> np.ndarray:
    import scipy.sparse as sp

    if sp.issparse(csc_or_dense):
        return np.asarray(csc_or_dense.getcol(idx).todense(), dtype=np.float64).ravel()
    return np.asarray(csc_or_dense[:, idx], dtype=np.float64).ravel()


def _cols_dense(csc_or_dense: Any, indices: Sequence[int]) -> np.ndarray:
    import scipy.sparse as sp

    sub = csc_or_dense[:, indices]
    if sp.issparse(sub):
        return np.asarray(sub.todense(), dtype=np.float64)
    return np.asarray(sub, dtype=np.float64)


def _nonzero_cells_per_col(mat: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return `(nnz_cells_per_col, col_sums)` without densifying sparse inputs."""
    import scipy.sparse as sp

    if sp.issparse(mat):
        csc = mat if sp.isspmatrix_csc(mat) else mat.tocsc()
        pos = csc.copy()
        pos.data = (pos.data > 0).astype(np.int32)
        nnz = np.asarray(pos.sum(axis=0)).ravel().astype(int)
        sums = np.asarray(csc.sum(axis=0)).ravel().astype(np.float64)
        return nnz, sums
    arr = np.asarray(mat, dtype=np.float64)
    return np.sum(arr > 0.0, axis=0).astype(int), np.sum(arr, axis=0).astype(np.float64)


def compute_peak_to_gene_links(
    *,
    rna_matrix: Any,
    rna_var_names: Sequence[str],
    atac_matrix: Any,
    atac_var_names: Sequence[str],
    gene_records: dict[str, dict[str, Any]],
    distance_bp: int = 500_000,
    min_distance_bp: int = 0,
    min_corr: float = 0.3,
    max_padj: float | None = None,
    max_pvalue_signac: float | None = None,
    rna_is_raw_counts: bool = False,
    atac_is_raw_counts: bool = False,
    preprocess: str = "as_is",
    atac_transform: str = "raw_counts",
    min_cells_detected: int = 1,
    standard_chroms_only: bool = False,
    background_null: bool = False,
    genome_fasta: str | None = None,
    n_bg_bins_per_axis: int = 5,
    n_bg_per_bin: int = 200,
    min_bg_peaks: int = 30,
    promoter_bp: int = 2000,
    distance_mode: str = "interval_overlap_zero",
    seed: int = 0,
) -> dict[str, Any]:
    """Core vectorized cis peak-to-gene correlation & trans-background linker.

    Parameters
    ----------
    preprocess:
        - `"as_is"`: correlate `rna_matrix` and `atac_matrix` directly without
          transforming values (used by `chromatin_regulatory._peak_to_gene` for
          exact ADR-050 backward compatibility).
        - `"paired"`: normalize raw RNA counts via `log1p(CP10k)` (matching
          `Signac::NormalizeData`) and apply `atac_transform` (`"raw_counts"`
          matching `Signac::LinkPeaks` `peak.slot="counts"`, or `"tfidf"`).
    """
    import scipy.sparse as sp

    n_cells = int(rna_matrix.shape[0])
    if int(atac_matrix.shape[0]) != n_cells:
        raise ValueError(
            f"RNA and ATAC cell counts differ ({n_cells} vs {atac_matrix.shape[0]})"
        )
    if n_cells < 3:
        return {
            "status": "error",
            "error_type": "InsufficientSharedCells",
            "details": f"Peak-to-gene requires >= 3 shared cells; got {n_cells}.",
        }

    warnings: list[str] = []

    # 1. Prepare RNA & ATAC matrices for correlation.
    if preprocess == "paired":
        if rna_is_raw_counts:
            rna_work = _log1p_cp10k_csc(rna_matrix)
            rna_norm_applied = "log1p_cp10k"
        else:
            rna_work = rna_matrix.tocsc() if sp.issparse(rna_matrix) else np.asarray(rna_matrix, dtype=np.float64)
            rna_norm_applied = "pre_normalized_input"
        if atac_transform == "tfidf" and atac_is_raw_counts:
            atac_work = _signac_tfidf_csc(atac_matrix)
            atac_norm_applied = "signac_tfidf"
        else:
            atac_work = atac_matrix.tocsc() if sp.issparse(atac_matrix) else np.asarray(atac_matrix, dtype=np.float64)
            atac_norm_applied = "raw_counts" if atac_is_raw_counts else "pre_normalized_input"
    else:
        rna_work = rna_matrix.tocsc() if sp.issparse(rna_matrix) else np.asarray(rna_matrix, dtype=np.float64)
        atac_work = atac_matrix.tocsc() if sp.issparse(atac_matrix) else np.asarray(atac_matrix, dtype=np.float64)
        rna_norm_applied = "as_is"
        atac_norm_applied = "as_is"

    # Clamp min_cells_detected safely on small unit-test matrices (< 50 cells).
    if int(min_cells_detected) <= 1:
        eff_min_cells = 1
    elif n_cells < 50:
        eff_min_cells = max(2, min(int(min_cells_detected), n_cells // 5))
    else:
        eff_min_cells = int(min_cells_detected)

    gene_nnz, gene_sums = _nonzero_cells_per_col(rna_matrix)
    peak_nnz, peak_raw_sums = _nonzero_cells_per_col(atac_matrix)

    peak_names = [str(p) for p in atac_var_names]
    peak_coords_all = [parse_peak_coord(p) for p in peak_names]

    # Build testable peak universe.
    valid_peak_indices: list[int] = []
    valid_peak_canons: list[str] = []
    valid_peak_coords: list[tuple[str, int, int]] = []
    n_unparseable_peaks = 0
    n_nonstandard_dropped = 0

    for pi, coord in enumerate(peak_coords_all):
        if coord is None:
            n_unparseable_peaks += 1
            continue
        if standard_chroms_only and not is_standard_contig(coord[0]):
            n_nonstandard_dropped += 1
            continue
        if eff_min_cells > 1:
            # Match Signac LinkPeaks `rowSums(x > 0) > min.cells` strictly (> eff_min_cells)
            # when eff_min_cells >= 10, or >= eff_min_cells on tiny cohorts.
            ok_cells = (
                peak_nnz[pi] > eff_min_cells
                if n_cells >= 50
                else peak_nnz[pi] >= eff_min_cells
            )
            if not ok_cells:
                continue
        elif peak_raw_sums[pi] <= 0.0 and preprocess == "paired":
            continue
        ckey = canonical_contig(coord[0]) or coord[0]
        valid_peak_indices.append(pi)
        valid_peak_canons.append(ckey)
        valid_peak_coords.append(coord)

    # Contig overlap sanity check against GTF.
    gtf_contigs = {
        rec["canon_chrom"] or rec["chrom"]
        for rec in gene_records.values()
    }
    if len(valid_peak_indices) >= 20 and gtf_contigs:
        n_contig_overlap = sum(1 for c in valid_peak_canons if c in gtf_contigs)
        frac_overlap = float(n_contig_overlap) / float(len(valid_peak_indices))
        if frac_overlap < 0.50:
            return {
                "status": "error",
                "error_type": "ContigMismatch",
                "details": (
                    f"Only {n_contig_overlap}/{len(valid_peak_indices)} "
                    f"({frac_overlap:.1%}) peaks map to contigs present in GTF."
                ),
            }

    peaks_by_chrom: dict[str, list[tuple[int, int, int, int]]] = {}
    for pos_in_valid, pi in enumerate(valid_peak_indices):
        _, s_bp, e_bp = valid_peak_coords[pos_in_valid]
        ckey = valid_peak_canons[pos_in_valid]
        peaks_by_chrom.setdefault(ckey, []).append((s_bp, e_bp, pi, pos_in_valid))

    # Map gene symbols to RNA column index (first occurrence if duplicates exist).
    gene_index: dict[str, int] = {}
    for gi, gname in enumerate(rna_var_names):
        gs = str(gname)
        if gs not in gene_index:
            gene_index[gs] = gi

    # Identify candidate genes that are expressed, annotated, and have >=1 cis peak.
    cand_genes: list[tuple[str, int, dict[str, Any], list[int], list[int], list[int], list[int]]] = []
    for gene_sym, rec in gene_records.items():
        gi = gene_index.get(gene_sym)
        if gi is None:
            continue
        if eff_min_cells > 1:
            ok_g = (
                gene_nnz[gi] > eff_min_cells
                if n_cells >= 50
                else gene_nnz[gi] >= eff_min_cells
            )
            if not ok_g:
                continue
        else:
            if gene_sums[gi] <= 0.0:
                expr_probe = _col_dense(rna_work, gi)
                if float(np.mean(expr_probe)) <= 0.0:
                    continue

        ckey = rec["canon_chrom"] or rec["chrom"]
        local_peaks = peaks_by_chrom.get(ckey)
        if not local_peaks:
            continue

        tss = int(rec["tss"])
        strand_sign = 1 if rec.get("strand", "+") != "-" else -1
        c_pidx: list[int] = []
        c_vpos: list[int] = []
        c_dist: list[int] = []
        c_off: list[int] = []

        for s_bp, e_bp, pi, vpos in local_peaks:
            if distance_mode == "endpoint_min":
                dist = min(abs(s_bp - tss), abs(e_bp - tss))
            else:
                dist = (
                    0
                    if (s_bp <= tss <= e_bp)
                    else min(abs(s_bp - tss), abs(e_bp - tss))
                )
            if int(min_distance_bp) <= dist <= int(distance_bp):
                mid = (s_bp + e_bp) // 2
                c_pidx.append(pi)
                c_vpos.append(vpos)
                c_dist.append(int(dist))
                c_off.append(int((mid - tss) * strand_sign))

        if c_pidx:
            cand_genes.append((gene_sym, gi, rec, c_pidx, c_vpos, c_dist, c_off))

    # 2. Optional BLAS-vectorized matched background null (Signac LinkPeaks style).
    bg_method = "none"
    bg_match = "none"
    gc_info: dict[str, Any] = {"available": False, "reason": "not_requested"}
    valid_peak_stratum = np.zeros(len(valid_peak_indices), dtype=int)
    bg_corrs_all: np.ndarray | None = None
    bg_canons: list[str] = []
    bg_strata: np.ndarray = np.zeros(0, dtype=int)
    bg_midpoints: np.ndarray = np.zeros(0, dtype=int)
    bg_acc_bin: np.ndarray = np.zeros(0, dtype=int)
    valid_acc_bin: np.ndarray = np.zeros(len(valid_peak_indices), dtype=int)
    n_bins_ax = max(1, int(n_bg_bins_per_axis))

    if background_null and len(valid_peak_indices) >= 5 and len(cand_genes) > 0:
        log_counts = np.log10(
            np.maximum(0.0, peak_raw_sums[valid_peak_indices]) + 1.0
        )
        lengths = np.asarray(
            [max(1, e - s) for _, s, e in valid_peak_coords], dtype=np.float64
        )
        log_len = np.log10(lengths + 1.0)

        gc_vals, gc_info = compute_peak_gc_from_indexed_fasta(
            genome_fasta, valid_peak_coords
        )
        eff_bins = min(n_bins_ax, max(1, len(valid_peak_indices) // 10))
        valid_acc_bin = _quantile_bin_labels(log_counts, eff_bins)
        if gc_vals is not None and gc_info.get("available"):
            sec_bin = _quantile_bin_labels(gc_vals, eff_bins)
            bg_match = "gc_and_accessibility"
        else:
            sec_bin = _quantile_bin_labels(log_len, eff_bins)
            bg_match = "accessibility_and_length"
        bg_method = "binned_matched_background"
        valid_peak_stratum = valid_acc_bin * eff_bins + sec_bin

        rng = np.random.default_rng(int(seed))
        n_strata = eff_bins * eff_bins
        sampled_vpos: list[int] = []
        for st_id in range(n_strata):
            members = np.flatnonzero(valid_peak_stratum == st_id)
            if members.size == 0:
                continue
            # Sort by original peak name for order-invariant deterministic sampling.
            members_sorted = sorted(
                members.tolist(),
                key=lambda vp: peak_names[valid_peak_indices[vp]],
            )
            if len(members_sorted) <= int(n_bg_per_bin):
                chosen = members_sorted
            else:
                pick = rng.choice(
                    len(members_sorted), size=int(n_bg_per_bin), replace=False
                )
                chosen = [members_sorted[int(p)] for p in np.sort(pick)]
            sampled_vpos.extend(chosen)

        sampled_vpos = sorted(set(sampled_vpos))
        bg_orig_pidx = [valid_peak_indices[vp] for vp in sampled_vpos]
        bg_mat = _cols_dense(atac_work, bg_orig_pidx)  # (n_cells, n_bg)
        bg_centered = bg_mat - np.mean(bg_mat, axis=0, keepdims=True)
        bg_norms = np.sqrt(np.sum(bg_centered * bg_centered, axis=0))
        keep_bg = bg_norms > 1e-12
        if int(np.sum(keep_bg)) >= 5:
            sampled_kept = [
                sampled_vpos[i] for i in range(len(sampled_vpos)) if keep_bg[i]
            ]
            a_bg_tilde = bg_centered[:, keep_bg] / bg_norms[keep_bg][None, :]
            bg_canons = [valid_peak_canons[vp] for vp in sampled_kept]
            bg_strata = np.asarray(
                [valid_peak_stratum[vp] for vp in sampled_kept], dtype=int
            )
            bg_acc_bin = np.asarray(
                [valid_acc_bin[vp] for vp in sampled_kept], dtype=int
            )
            bg_midpoints = np.asarray(
                [
                    (valid_peak_coords[vp][1] + valid_peak_coords[vp][2]) // 2
                    for vp in sampled_kept
                ],
                dtype=int,
            )

            # Standardize all candidate gene expression columns in chunks and
            # compute exact background correlations via BLAS DGEMM!
            gene_cols = [g_item[1] for g_item in cand_genes]
            e_mat = _cols_dense(rna_work, gene_cols).T  # (n_cand_genes, n_cells)
            e_centered = e_mat - np.mean(e_mat, axis=1, keepdims=True)
            e_norms = np.sqrt(np.sum(e_centered * e_centered, axis=1, keepdims=True))
            e_tilde = np.where(e_norms > 1e-12, e_centered / np.maximum(e_norms, 1e-12), 0.0)
            bg_corrs_all = np.clip(e_tilde @ a_bg_tilde, -1.0, 1.0)

    # 3. Evaluate all in-window cis `(gene, peak)` pairs.
    bg_canon_arr = np.asarray(bg_canons, dtype=object) if bg_canons else np.empty(0, dtype=object)
    eff_min_bg = min(int(min_bg_peaks), max(5, len(bg_canons) // 3)) if bg_canons else 5
    used_distal_fallback = False

    tested_records: list[dict[str, Any]] = []
    n_tested = 0

    for cg_idx, (gene_sym, gi, rec, c_pidx, c_vpos, c_dist, c_off) in enumerate(
        cand_genes
    ):
        expr = _col_dense(rna_work, gi)
        if float(np.std(expr)) == 0.0:
            if preprocess != "as_is":
                continue
        sub_atac = _cols_dense(atac_work, c_pidx)
        rvals = vectorized_pearson_corr(expr, sub_atac)

        ckey = rec["canon_chrom"] or rec["chrom"]
        tss = int(rec["tss"])

        # Precompute trans mask for this gene's chromosome if background_null active.
        trans_mask: np.ndarray | None = None
        bg_row: np.ndarray | None = None
        if bg_corrs_all is not None:
            bg_row = bg_corrs_all[cg_idx]
            trans_mask = bg_canon_arr != ckey

        for j, pi in enumerate(c_pidx):
            n_tested += 1
            rv = float(rvals[j])
            if not math.isfinite(rv):
                continue
            r_round = round(rv, 6)

            z_val: float | None = None
            p_two: float | None = None
            p_one: float | None = None
            bg_scope = "not_run"

            if bg_row is not None and trans_mask is not None:
                vpos = c_vpos[j]
                st_id = int(valid_peak_stratum[vpos])
                acc_b = int(valid_acc_bin[vpos])

                # Tier 1: exact trans-chromosomal matched stratum
                idx_pool = np.flatnonzero(trans_mask & (bg_strata == st_id))
                if idx_pool.size < eff_min_bg:
                    # Expand within trans chromosomes to same accessibility bin, then all trans
                    idx_acc = np.flatnonzero(trans_mask & (bg_acc_bin == acc_b))
                    idx_pool = (
                        idx_acc
                        if idx_acc.size >= eff_min_bg
                        else np.flatnonzero(trans_mask)
                    )

                if idx_pool.size >= eff_min_bg:
                    bg_scope = "trans"
                else:
                    # Tier 2: single-chromosome synthetic fixture fallback -> distal peaks outside cis window
                    distal_mask = (bg_canon_arr == ckey) & (
                        np.abs(bg_midpoints - tss) > int(distance_bp)
                    )
                    idx_distal = np.flatnonzero(distal_mask)
                    if idx_distal.size >= 5:
                        idx_pool = idx_distal
                        bg_scope = "same_chrom_distal"
                        used_distal_fallback = True
                    else:
                        idx_pool = np.empty(0, dtype=int)
                        bg_scope = "unavailable"

                if idx_pool.size >= 5:
                    vals_bg = bg_row[idx_pool]
                    mu_bg = float(np.mean(vals_bg))
                    sd_bg = float(np.std(vals_bg, ddof=1)) if idx_pool.size > 1 else 0.0
                    if sd_bg > 1e-12 and math.isfinite(sd_bg):
                        z_raw = (rv - mu_bg) / sd_bg
                        z_val = round(float(z_raw), 6)
                        # Exact Gaussian tail via math.erfc(|z| / sqrt(2))
                        p_two_raw = float(math.erfc(abs(z_raw) / math.sqrt(2.0)))
                        p_two = max(0.0, min(1.0, p_two_raw))
                        p_one = 0.5 * p_two

            tested_records.append(
                {
                    "gene": gene_sym,
                    "gene_id": rec.get("gene_id", gene_sym),
                    "peak": peak_names[pi],
                    "correlation": r_round,
                    "zscore": z_val,
                    "p_value": p_two,
                    "p_value_signac": p_one,
                    "distance_bp": c_dist[j],
                    "distance_kb": round(float(c_dist[j]) / 1000.0, 3),
                    "tss_offset_bp": c_off[j],
                    "direction": "positive" if r_round > 0.0 else "negative",
                    "region_class": (
                        "promoter_proximal"
                        if c_dist[j] <= int(promoter_bp)
                        else "distal_cis"
                    ),
                    "background_scope": bg_scope,
                }
            )

    if used_distal_fallback:
        warnings.append(
            "Background trans-chromosomal peak pool had fewer than minimum peaks "
            "for one or more genes; used same_chrom_distal fallback outside cis window."
        )

    # 4. Pre-registered Benjamini-Hochberg FDR across ALL tested valid cis pairs
    # BEFORE applying `min_corr` / `max_padj` filters (anti-selection-bias invariant).
    n_fdr_sig_all = 0
    if background_null and tested_records:
        finite_idx = [
            i for i, r in enumerate(tested_records) if r["p_value"] is not None
        ]
        if finite_idx:
            p_arr = np.asarray(
                [float(tested_records[i]["p_value"]) for i in finite_idx],
                dtype=np.float64,
            )
            padj_arr = bh_correct(p_arr)
            for pos, rec_idx in enumerate(finite_idx):
                padj_val = float(padj_arr[pos])
                tested_records[rec_idx]["padj"] = padj_val
                if padj_val <= 0.05:
                    n_fdr_sig_all += 1
        for r in tested_records:
            r.setdefault("padj", None)

    # 5. Apply correlation & statistical cutoffs.
    surviving: list[dict[str, Any]] = []
    for r in tested_records:
        if abs(float(r["correlation"])) < float(min_corr):
            continue
        if background_null and max_padj is not None:
            if r.get("padj") is None or float(r["padj"]) > float(max_padj):
                continue
        if background_null and max_pvalue_signac is not None:
            if (
                r.get("p_value_signac") is None
                or float(r["p_value_signac"]) > float(max_pvalue_signac)
            ):
                continue
        surviving.append(r)

    surviving.sort(
        key=lambda r: (
            -abs(float(r["correlation"])),
            -(abs(float(r["zscore"])) if r.get("zscore") is not None else 0.0),
            r["gene"],
            r["peak"],
        )
    )

    n_pos = sum(1 for r in surviving if r["direction"] == "positive")
    n_neg = sum(1 for r in surviving if r["direction"] == "negative")
    n_prom = sum(1 for r in surviving if r["region_class"] == "promoter_proximal")
    n_dist = sum(1 for r in surviving if r["region_class"] == "distal_cis")

    return {
        "status": "success",
        "n_cells": int(n_cells),
        "n_genes_tested": int(len(cand_genes)),
        "n_peaks_tested": int(len(valid_peak_indices)),
        "n_unparseable_peaks": int(n_unparseable_peaks),
        "n_nonstandard_peaks_dropped": int(n_nonstandard_dropped),
        "n_tested": int(n_tested),
        "n_fdr_significant_all_tested": int(n_fdr_sig_all),
        "n_links": int(len(surviving)),
        "n_positive": int(n_pos),
        "n_negative": int(n_neg),
        "n_promoter_proximal": int(n_prom),
        "n_distal_cis": int(n_dist),
        "rna_normalization": rna_norm_applied,
        "atac_normalization": atac_norm_applied,
        "min_cells_detected_effective": int(eff_min_cells),
        "background_method": bg_method,
        "background_match": bg_match,
        "gc_matching": gc_info,
        "links": surviving,
        "warnings": warnings,
    }
