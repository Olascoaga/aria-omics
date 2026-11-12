"""
Unified Paired Candidate eRegulon Engine (V48 Slice I6 — Tier 2b Regulatory Layer)
---------------------------------------------------------------------------------
Assembles deterministic, strictly associative candidate enhancer-driven regulons
(`TF -> cis peak -> target gene`) across paired same-cell RNA + ATAC profiles by
combining three independently controlled statistical legs:

  1. Sequence-specific motif cistrome (`TF -> Peak`):
     - Offline versioned MEME PWM parsing (`JASPAR2024_CORE_vertebrates` default).
     - Strict 1-to-1 symbol mapping to expressed dataset genes (`>= min_cells_detected`),
       stripping technical `(var.N)` suffixes, requiring ALL subunits of `::`
       heterodimers to be expressed (`motif_source="heterodimer"`), logging
       unmapped/fusion/species-specific symbols without any hardcoded gene or
       ortholog dictionary (`ADR-005`, `ADR-011`).
     - Fixed-width peak-center DNA window extraction (`scan_window_bp=500`,
       `center +/- 250 bp`) via pure-Python `.fai` byte-offset seeking, eliminating
       variable-peak-length `2*L*p` bias.
     - Dual-backend motif scanning (`snapatac2_rust` parallel `PyDNAMotifScanner`
       primary; `aria_numpy_dp` discrete dynamic-programming log-odds fallback for
       environments without `snapatac2`; never mixed within a run).
     - GC-quintile Cochran-Mantel-Haenszel (CMH) stratified cistrome enrichment in
       cis-linked peaks vs tested peak background.

  2. Cis chromatin-to-gene association (`Peak -> Gene`):
     - Reuses `aria.utils.peak2gene_core.compute_peak_to_gene_links`
       (`preprocess="paired"`, `(GC, accessibility)` trans-chromosomal matched
       background null, shared `bh_correct` `padj`), defaulting to `p2g_min_corr=0.05`
       and `p2g_max_padj=0.05` so single-cell sparsity does not starve regulon
       assembly while maintaining genome-wide FDR control.

  3. Paired TF co-variation & covariate-matched module calibration (`TF -> Gene/Peak`):
     - Genome-wide paired Pearson `r_tf_gene` (`log1p(CP10k)` RNA) and `r_tf_peak`,
       excluding `g == TF` self-loops (recorded separately as `self_linked_peaks`)
       and zero-variance features.
     - Library-depth-partial correlations (`r_tf_gene_partial_depth`,
       `r_tf_peak_partial_depth`) conditioning on `log10(1 + total_counts)` to flag
       sequencing-depth-confounded candidates (`depth_confounded`).
     - Descriptive non-causal sign-class partitioning (`pos_tf_pos_p2g`,
       `neg_tf_pos_p2g`, `pos_tf_neg_p2g`, `neg_tf_neg_p2g` — NEVER labeled
       "activator" or "repressor" per `ADR-064`).
     - Module-level significance calibrated against a `(log1p_mean_expr,
       log10_cells_detected)` 2D bin-stratified finite-population null (analytic
       variance + `n_null` seeded Monte Carlo draws), `bh_correct` across ALL
       tested `(motif_id, tf, sign_class)` modules BEFORE selecting the best motif
       per `(tf, sign_class)` (preventing variant-union inflation), and weighted
       GSEA leading-edge pruning (`|r_tf_gene| >= min_tf_gene_corr`,
       `n_target_genes >= min_target_genes`).
     - Paralog motif-sharing ambiguity disclosure (`motif_family_ambiguity`,
       `ambiguous_paralog_tfs`) using cis-peak Jaccard and aligned PWM similarity,
       excluding physical `::` heterodimer subunit pairs.
     - Deterministic rank-based AUCell per-cell regulon activity scoring on both
       RNA target genes and ATAC cis peaks.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re
from typing import Any, Sequence

import numpy as np
from scipy import sparse, stats

from aria.utils.peak2gene_core import (
    canonical_contig,
    compute_peak_to_gene_links,
    parse_peak_coord,
)
from aria.utils.stats import bh_correct


class ContigMismatchError(ValueError):
    """Raised when peak contigs fail to match GTF or indexed FASTA contigs."""

_VARIANT_SUFFIX_RE = re.compile(r"\s*\(var\.(\d+)\)\s*$", flags=re.IGNORECASE)
_ALPHABET = b"ACGT"
_SIGN_CLASSES = (
    "pos_tf_pos_p2g",
    "neg_tf_pos_p2g",
    "pos_tf_neg_p2g",
    "neg_tf_neg_p2g",
)


@dataclass(frozen=True)
class MemeMotif:
    """Parsed MEME v4 nucleotide position frequency matrix."""

    motif_id: str
    raw_name: str
    pfm: np.ndarray  # shape (w, 4) over A, C, G, T
    bg_probs: tuple[float, float, float, float] = (0.25, 0.25, 0.25, 0.25)

    @property
    def width(self) -> int:
        return int(self.pfm.shape[0])


@dataclass(frozen=True)
class ResolvedMotifTF:
    """Mapping from one MEME motif to one expressed dataset TF gene."""

    motif_idx: int
    motif_id: str
    raw_motif_name: str
    tf_gene: str
    motif_source: str  # "monomer" | "heterodimer"
    match_mode: str  # "exact" | "case_folded"
    heterodimer_partners: tuple[str, ...]
    variant_tag: str | None


# --------------------------------------------------------------------------- #
# 1. MEME parsing & strict dataset-backed TF symbol resolution                #
# --------------------------------------------------------------------------- #


def parse_meme_motifs(meme_path: str | Path) -> list[MemeMotif]:
    """Parse all nucleotide PWMs/PFMs from a local MEME v4 file."""
    p = Path(meme_path)
    if not p.is_file():
        raise FileNotFoundError(f"MEME motif file not found: {meme_path}")

    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    bg = [0.25, 0.25, 0.25, 0.25]
    motifs: list[MemeMotif] = []

    i = 0
    n = len(lines)
    while i < n:
        s = lines[i].strip()
        if s.lower().startswith("background letter frequencies"):
            if i + 1 < n:
                toks = lines[i + 1].strip().split()
                parsed_bg: dict[str, float] = {}
                for k in range(0, len(toks) - 1, 2):
                    try:
                        parsed_bg[toks[k].upper()] = float(toks[k + 1])
                    except ValueError:
                        pass
                if all(b in parsed_bg for b in ("A", "C", "G", "T")):
                    tot = sum(parsed_bg[b] for b in ("A", "C", "G", "T"))
                    if tot > 0:
                        bg = [parsed_bg[b] / tot for b in ("A", "C", "G", "T")]
            i += 2
            continue

        if s.startswith("MOTIF "):
            rest = s[6:].strip()
            parts = rest.split(None, 1)
            motif_id = parts[0]
            raw_name = parts[1].strip() if len(parts) > 1 else motif_id
            i += 1
            while i < n and not lines[i].strip().lower().startswith(
                "letter-probability matrix"
            ):
                if lines[i].strip().startswith("MOTIF "):
                    break
                i += 1
            if i >= n or not lines[i].strip().lower().startswith(
                "letter-probability matrix"
            ):
                continue
            hdr = lines[i].strip()
            m_w = re.search(r"\bw\s*=\s*(\d+)", hdr)
            i += 1
            rows: list[list[float]] = []
            if m_w:
                w = int(m_w.group(1))
                for _ in range(w):
                    if i >= n:
                        break
                    vals = [float(x) for x in lines[i].strip().split()[:4]]
                    if len(vals) == 4:
                        rows.append(vals)
                    i += 1
            else:
                while i < n:
                    row_s = lines[i].strip()
                    if not row_s or row_s.startswith(("URL", "MOTIF")):
                        break
                    try:
                        vals = [float(x) for x in row_s.split()[:4]]
                    except ValueError:
                        break
                    if len(vals) != 4:
                        break
                    rows.append(vals)
                    i += 1
            if rows:
                arr = np.asarray(rows, dtype=np.float64)
                row_sums = arr.sum(axis=1, keepdims=True)
                row_sums = np.where(row_sums > 0, row_sums, 1.0)
                arr = arr / row_sums
                motifs.append(
                    MemeMotif(
                        motif_id=motif_id,
                        raw_name=raw_name,
                        pfm=arr,
                        bg_probs=(bg[0], bg[1], bg[2], bg[3]),
                    )
                )
            continue
        i += 1
    return motifs


def resolve_motifs_to_expressed_tfs(
    motifs: Sequence[MemeMotif],
    gene_names: Sequence[str],
    gene_cells_detected: np.ndarray,
    *,
    min_cells_detected: int = 10,
) -> tuple[list[ResolvedMotifTF], list[dict[str, str]]]:
    """Map MEME motifs strictly to expressed dataset genes (zero hardcoded panels).

    Rules (`ADR-005`, `ADR-011`, `Concilio d9842f`):
      - Strip JASPAR technical variant suffix `(var.N)` (recorded in `variant_tag`).
      - Split `::` heterodimers; require ALL subunits to map uniquely to expressed
        dataset genes (`>= min_cells_detected` cells) before assigning the dimer
        motif to each subunit (`motif_source="heterodimer"`).
      - Prefer exact symbol match (`match_mode="exact"`), falling back to unique
        1-to-1 case-folded symbol match (`match_mode="case_folded"` for conserved
        vertebrate JASPAR titles like `Arnt` -> `ARNT`, `Sox2` -> `SOX2`).
      - Ambiguous case collisions, non-expressed TFs, fusion symbols (`EWSR1-FLI1`),
        and non-matching species symbols (`Zfp809`) are recorded in `unmapped`.
    """
    expressed_exact: set[str] = set()
    upper_buckets: dict[str, list[str]] = {}
    all_upper_buckets: dict[str, list[str]] = {}
    for idx, g_raw in enumerate(gene_names):
        g = str(g_raw).strip()
        if not g:
            continue
        all_upper_buckets.setdefault(g.upper(), []).append(g)
        if int(gene_cells_detected[idx]) >= min_cells_detected:
            expressed_exact.add(g)
            upper_buckets.setdefault(g.upper(), []).append(g)

    def _lookup_one(token: str) -> tuple[str | None, str | None, str | None]:
        tok = token.strip()
        if not tok:
            return None, None, "empty_symbol"
        if tok in expressed_exact:
            return tok, "exact", None
        cands = upper_buckets.get(tok.upper(), [])
        if len(cands) == 1:
            return cands[0], "case_folded", None
        if len(cands) > 1:
            return None, None, "ambiguous_case_collision"
        if tok.upper() in all_upper_buckets:
            return None, None, "tf_below_min_cells_detected"
        return None, None, "symbol_not_in_rna_universe"

    resolved: list[ResolvedMotifTF] = []
    unmapped: list[dict[str, str]] = []

    for m_idx, m in enumerate(motifs):
        base_name = m.raw_name.strip()
        var_match = _VARIANT_SUFFIX_RE.search(base_name)
        variant_tag = f"var.{var_match.group(1)}" if var_match else None
        cleaned = _VARIANT_SUFFIX_RE.sub("", base_name).strip()

        raw_parts = [p.strip() for p in cleaned.split("::") if p.strip()]
        if not raw_parts:
            unmapped.append(
                {
                    "motif_id": m.motif_id,
                    "raw_name": m.raw_name,
                    "reason": "empty_symbol",
                }
            )
            continue

        if len(raw_parts) == 1:
            gene_sym, mode, err = _lookup_one(raw_parts[0])
            if gene_sym is None or mode is None:
                unmapped.append(
                    {
                        "motif_id": m.motif_id,
                        "raw_name": m.raw_name,
                        "reason": err or "unmapped",
                    }
                )
                continue
            resolved.append(
                ResolvedMotifTF(
                    motif_idx=m_idx,
                    motif_id=m.motif_id,
                    raw_motif_name=m.raw_name,
                    tf_gene=gene_sym,
                    motif_source="monomer",
                    match_mode=mode,
                    heterodimer_partners=(),
                    variant_tag=variant_tag,
                )
            )
        else:
            mapped_parts: list[tuple[str, str]] = []
            fail_reason: str | None = None
            for part in raw_parts:
                gene_sym, mode, err = _lookup_one(part)
                if gene_sym is None or mode is None:
                    fail_reason = f"heterodimer_subunit_{err}:{part}"
                    break
                mapped_parts.append((gene_sym, mode))
            if fail_reason is not None:
                unmapped.append(
                    {
                        "motif_id": m.motif_id,
                        "raw_name": m.raw_name,
                        "reason": fail_reason,
                    }
                )
                continue
            partner_tuple = tuple(dict.fromkeys(g for g, _ in mapped_parts))
            seen_subunits: set[str] = set()
            for gene_sym, mode in mapped_parts:
                if gene_sym in seen_subunits:
                    continue
                seen_subunits.add(gene_sym)
                resolved.append(
                    ResolvedMotifTF(
                        motif_idx=m_idx,
                        motif_id=m.motif_id,
                        raw_motif_name=m.raw_name,
                        tf_gene=gene_sym,
                        motif_source="heterodimer",
                        match_mode=mode,
                        heterodimer_partners=partner_tuple,
                        variant_tag=variant_tag,
                    )
                )

    return resolved, unmapped


# --------------------------------------------------------------------------- #
# 2. Fixed-width `.fai` sequence extraction & dual-backend motif scanning     #
# --------------------------------------------------------------------------- #


def _read_fai_index(fai_path: Path) -> dict[str, tuple[int, int, int, int]]:
    out: dict[str, tuple[int, int, int, int]] = {}
    with open(fai_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            try:
                out[parts[0]] = (
                    int(parts[1]),
                    int(parts[2]),
                    int(parts[3]),
                    int(parts[4]),
                )
            except ValueError:
                continue
    return out


def fetch_peak_center_sequences_from_indexed_fasta(
    peak_names: Sequence[str],
    fasta_path: str | Path,
    *,
    scan_window_bp: int = 500,
    max_n_fraction: float = 0.50,
) -> dict[str, Any]:
    """Extract fixed-width `scan_window_bp` DNA sequences around each peak center.

    Uses pure-Python `.fai` byte-offset seeking (zero `pysam` requirement) and
    resolves UCSC (`chr1`) <-> Ensembl (`1`) contig prefixes transparently.
    """
    fa = Path(fasta_path)
    fai = Path(f"{fa}.fai")
    if not fa.is_file() or not fai.is_file():
        raise FileNotFoundError(
            f"Indexed reference FASTA (.fa + .fai) required: {fasta_path}"
        )
    idx = _read_fai_index(fai)
    if not idx:
        raise ValueError(f"Empty or malformed FASTA index: {fai}")

    canon_to_raw: dict[str, str] = {}
    for raw_c in idx:
        canon_to_raw.setdefault(canonical_contig(raw_c), raw_c)

    half_w = max(1, int(scan_window_bp) // 2)
    width = 2 * half_w

    sequences: list[str | None] = [None] * len(peak_names)
    gc_values = np.full(len(peak_names), np.nan, dtype=np.float64)
    valid_mask = np.zeros(len(peak_names), dtype=bool)
    n_missing_contig = 0
    n_high_n = 0

    by_contig: dict[str, list[tuple[int, int, int]]] = {}
    for i, pname in enumerate(peak_names):
        parsed = parse_peak_coord(pname)
        if parsed is None:
            continue
        chrom, s, e = parsed
        raw_c = chrom if chrom in idx else canon_to_raw.get(chrom)
        if raw_c is None:
            n_missing_contig += 1
            continue
        contig_len = idx[raw_c][0]
        mid = (int(s) + int(e)) // 2
        w_start = max(0, mid - half_w)
        w_end = min(contig_len, w_start + width)
        if w_end - w_start < width and w_end == contig_len:
            w_start = max(0, contig_len - width)
        if w_end <= w_start:
            continue
        by_contig.setdefault(raw_c, []).append((i, w_start, w_end))

    with open(fa, "rb") as fh:
        for raw_c, intervals in by_contig.items():
            _clen, offset, line_bases, line_width = idx[raw_c]
            intervals.sort(key=lambda t: t[1])
            for peak_idx, s0, e0 in intervals:
                byte_start = offset + (s0 // line_bases) * line_width + (s0 % line_bases)
                end_excl = e0 - 1
                byte_end = (
                    offset
                    + (end_excl // line_bases) * line_width
                    + (end_excl % line_bases)
                    + 1
                )
                fh.seek(byte_start)
                chunk = fh.read(byte_end - byte_start)
                seq_bytes = (
                    chunk.replace(b"\n", b"").replace(b"\r", b"").upper()
                )
                seq_len = len(seq_bytes)
                if seq_len == 0:
                    continue
                n_cnt = seq_bytes.count(b"N")
                if (n_cnt / float(seq_len)) > max_n_fraction:
                    n_high_n += 1
                    continue
                gc = seq_bytes.count(b"G") + seq_bytes.count(b"C")
                at = seq_bytes.count(b"A") + seq_bytes.count(b"T")
                denom = gc + at
                if denom == 0:
                    n_high_n += 1
                    continue
                sequences[peak_idx] = seq_bytes.decode("ascii", errors="replace")
                gc_values[peak_idx] = float(gc) / float(denom)
                valid_mask[peak_idx] = True

    return {
        "sequences": sequences,
        "gc_fraction": gc_values,
        "valid_mask": valid_mask,
        "scan_window_bp": width,
        "n_valid_sequences": int(valid_mask.sum()),
        "n_missing_contig": int(n_missing_contig),
        "n_masked_high_n": int(n_high_n),
    }


def _numpy_dp_threshold_and_ipwm(
    motif: MemeMotif,
    *,
    pvalue: float,
    pseudocount: float = 1e-4,
    scale: float = 1000.0,
) -> tuple[np.ndarray, int]:
    """Compute integer-scaled log-odds PWM and exact discrete-DP score threshold.

    Matches MOODS / `snapatac2` Rust `PyDNAMotifScanner` parameterization:
      - `pfm_adj = (pfm + pseudocount) / (1 + 4 * pseudocount)`
      - natural-log odds against `(a, c, g, t)` background probabilities
      - exact discrete convolution over motif positions
      - clamps threshold to the maximum attainable score when `pvalue` is smaller
        than the minimum possible tail probability of a short (`w <= 6`) motif.
    """
    bg = np.asarray(motif.bg_probs, dtype=np.float64)
    bg = bg / bg.sum()
    pfm_adj = (motif.pfm + pseudocount) / (1.0 + 4.0 * pseudocount)
    pwm = np.log(pfm_adj / bg[None, :])
    ipwm = np.rint(pwm * scale).astype(np.int32)

    dp = np.array([1.0], dtype=np.float64)
    offset = 0
    for r in ipwm:
        rmin = int(r.min())
        rmax = int(r.max())
        ndp = np.zeros(len(dp) + (rmax - rmin), dtype=np.float64)
        for b in range(4):
            sh = int(r[b]) - rmin
            ndp[sh : sh + len(dp)] += bg[b] * dp
        dp = ndp
        offset += rmin

    sf = np.cumsum(dp[::-1])[::-1]
    idx = np.where(sf <= float(pvalue))[0]
    if len(idx) > 0:
        ithr = int(idx[0]) + offset
    else:
        ithr = (len(sf) - 1) + offset
    return ipwm, ithr


def scan_motifs_numpy_dp(
    motifs: Sequence[MemeMotif],
    sequences: Sequence[str],
    *,
    pvalue: float = 5e-5,
    rc: bool = True,
) -> np.ndarray:
    """Pure-NumPy discrete-DP PWM scanner returning `(n_motifs, n_seqs)` bool matrix."""
    n_m = len(motifs)
    n_s = len(sequences)
    hits = np.zeros((n_m, n_s), dtype=bool)
    if n_m == 0 or n_s == 0:
        return hits

    # Encode sequences once: A=0, C=1, G=2, T=3, any ambiguous/N=4 (neutral 0 score)
    bmap = np.full(256, 4, dtype=np.int8)
    bmap[ord("A")] = 0
    bmap[ord("C")] = 1
    bmap[ord("G")] = 2
    bmap[ord("T")] = 3
    bmap[ord("a")] = 0
    bmap[ord("c")] = 1
    bmap[ord("g")] = 2
    bmap[ord("t")] = 3

    # Group equal-length sequences (typically 99%+ are exact `scan_window_bp=500`!)
    by_len: dict[int, list[int]] = {}
    encoded_list: list[np.ndarray] = []
    for idx, seq in enumerate(sequences):
        arr = bmap[np.frombuffer(seq.encode("ascii", errors="replace"), dtype=np.uint8)]
        encoded_list.append(arr)
        by_len.setdefault(len(arr), []).append(idx)

    len_batches: list[tuple[int, np.ndarray, np.ndarray]] = []
    for seq_len, idxs in by_len.items():
        if seq_len == 0:
            continue
        mat = np.vstack([encoded_list[i] for i in idxs])  # (batch, seq_len)
        len_batches.append((seq_len, np.asarray(idxs, dtype=np.int64), mat))

    for m_idx, motif in enumerate(motifs):
        w = motif.width
        ipwm4, ithr = _numpy_dp_threshold_and_ipwm(motif, pvalue=pvalue)
        # Pad 5th column = 0 for N/ambiguous bases
        ipwm5 = np.zeros((w, 5), dtype=np.int32)
        ipwm5[:, :4] = ipwm4
        rc_ipwm5 = np.zeros((w, 5), dtype=np.int32)
        rc_ipwm5[:, :4] = ipwm4[::-1, ::-1]
        pos_idx = np.arange(w)[:, None, None]

        for seq_len, idxs, batch_mat in len_batches:
            if seq_len < w:
                continue
            # windows shape: (batch, seq_len - w + 1, w) -> transpose to (w, batch, n_win)
            win = np.lib.stride_tricks.sliding_window_view(
                batch_mat, window_shape=w, axis=1
            )
            win_t = np.transpose(win, (2, 0, 1))
            max_fwd = ipwm5[pos_idx, win_t].sum(axis=0).max(axis=1)
            if rc:
                max_rev = rc_ipwm5[pos_idx, win_t].sum(axis=0).max(axis=1)
                passed = np.maximum(max_fwd, max_rev) >= ithr
            else:
                passed = max_fwd >= ithr
            hits[m_idx, idxs] = passed

    return hits


def scan_motifs_in_sequences(
    meme_path: str | Path,
    motifs: Sequence[MemeMotif],
    sequences: Sequence[str],
    *,
    motif_Indices: Sequence[int] | None = None,
    pvalue: float = 5e-5,
    rc: bool = True,
    backend: str = "auto",
) -> tuple[np.ndarray, str]:
    """Scan selected `motif_Indices` across `sequences`, returning `(hits_bool, backend_used)`.

    Never mixes backends within a call (`Concilio d9842f`).
    """
    if not (0.0 < float(pvalue) < 1.0):
        raise ValueError(f"motif pvalue must be strictly in (0, 1), got {pvalue}")

    sel_indices = (
        list(range(len(motifs)))
        if motif_Indices is None
        else [int(i) for i in motif_Indices]
    )
    if not sel_indices or not sequences:
        used = "aria_numpy_dp" if backend == "numpy" else "snapatac2_rust"
        return np.zeros((len(sel_indices), len(sequences)), dtype=bool), used

    req = backend.strip().lower()
    if req not in ("auto", "snapatac2", "snapatac2_rust", "numpy", "aria_numpy_dp"):
        raise ValueError(f"Unsupported motif scanner backend: {backend}")

    if req in ("auto", "snapatac2", "snapatac2_rust"):
        try:
            import snapatac2 as snap  # type: ignore

            snap_motifs = snap.datasets.read_motifs(str(meme_path))
            id_map = {}
            for sm in snap_motifs:
                sm_id = str(sm.id).split()[0].strip()
                id_map[sm_id] = sm

            seqs_list = [str(s) for s in sequences]
            rows: list[np.ndarray] = []
            all_found = True
            for m_idx in sel_indices:
                mm = motifs[m_idx]
                sm = id_map.get(mm.motif_id)
                if sm is None:
                    all_found = False
                    break
                bg = mm.bg_probs
                scanner = sm.with_nucl_prob(a=bg[0], c=bg[1], g=bg[2], t=bg[3])
                rows.append(
                    np.asarray(
                        scanner.exists(seqs_list, pvalue=float(pvalue), rc=bool(rc)),
                        dtype=bool,
                    )
                )
            if all_found:
                return np.vstack(rows), "snapatac2_rust"
            if req in ("snapatac2", "snapatac2_rust"):
                raise RuntimeError(
                    "Requested snapatac2_rust scanner but motif IDs mismatched MEME file."
                )
        except ImportError:
            if req in ("snapatac2", "snapatac2_rust"):
                raise

    sub_motifs = [motifs[i] for i in sel_indices]
    return (
        scan_motifs_numpy_dp(sub_motifs, sequences, pvalue=pvalue, rc=rc),
        "aria_numpy_dp",
    )


# --------------------------------------------------------------------------- #
# 3. GC-stratified CMH motif enrichment & PWM similarity                      #
# --------------------------------------------------------------------------- #


def compute_gc_stratified_motif_enrichment(
    motif_hit_matrix: np.ndarray,
    foreground_mask: np.ndarray,
    peak_gc: np.ndarray,
    *,
    n_gc_bins: int = 5,
) -> dict[str, np.ndarray]:
    """One-sided Cochran-Mantel-Haenszel GC-stratified motif enrichment across peaks.

    Controls for CpG-island / GC-composition differences between cis-linked peaks
    and unlinked background peaks (`Concilio d9842f` Blocking #1).
    """
    n_motifs, n_peaks = motif_hit_matrix.shape
    fg = np.asarray(foreground_mask, dtype=bool)
    if n_motifs == 0 or n_peaks == 0 or fg.sum() == 0 or (~fg).sum() == 0:
        return {
            "cmh_z": np.zeros(n_motifs, dtype=np.float64),
            "pvalue": np.ones(n_motifs, dtype=np.float64),
            "log2_fold_enrichment": np.zeros(n_motifs, dtype=np.float64),
        }

    qs = np.linspace(0.0, 1.0, n_gc_bins + 1)
    edges = np.quantile(peak_gc, qs)
    edges[0] -= 1e-9
    edges[-1] += 1e-9
    bins = np.clip(np.digitize(peak_gc, edges[1:-1]), 0, n_gc_bins - 1)

    hits_f = motif_hit_matrix.astype(np.float64)
    num = np.zeros(n_motifs, dtype=np.float64)
    var = np.zeros(n_motifs, dtype=np.float64)
    exp_hits = np.zeros(n_motifs, dtype=np.float64)
    obs_hits = hits_f[:, fg].sum(axis=1)

    for b in range(n_gc_bins):
        in_b = bins == b
        n_b = int(in_b.sum())
        if n_b < 2:
            continue
        fg_b = in_b & fg
        n1_b = int(fg_b.sum())
        n0_b = n_b - n1_b
        if n1_b == 0 or n0_b == 0:
            continue
        m1_b = hits_f[:, in_b].sum(axis=1)  # total motif hits in GC stratum b
        m0_b = n_b - m1_b
        a_b = hits_f[:, fg_b].sum(axis=1)  # observed fg motif hits in stratum b
        e_b = (n1_b * m1_b) / float(n_b)
        v_b = (n1_b * n0_b * m1_b * m0_b) / (float(n_b * n_b) * float(n_b - 1))
        num += a_b - e_b
        var += v_b
        exp_hits += e_b

    z = np.where(var > 1e-12, num / np.sqrt(np.maximum(var, 1e-12)), 0.0)
    pvals = stats.norm.sf(z)
    log2_fe = np.log2((obs_hits + 0.5) / (exp_hits + 0.5))
    return {
        "cmh_z": z,
        "pvalue": pvals,
        "log2_fold_enrichment": log2_fe,
    }


def compute_max_aligned_pwm_similarity(pfm_a: np.ndarray, pfm_b: np.ndarray) -> float:
    """Maximum column-wise Pearson correlation over all >=4 bp alignments (fwd + RC)."""
    a = np.asarray(pfm_a, dtype=np.float64)
    b = np.asarray(pfm_b, dtype=np.float64)
    b_rc = b[::-1, ::-1]
    wa, wb = a.shape[0], b.shape[0]
    min_ov = min(4, wa, wb)
    best = -1.0

    for cand in (b, b_rc):
        for shift in range(-(wb - min_ov), wa - min_ov + 1):
            a_s = max(0, shift)
            b_s = max(0, -shift)
            ov = min(wa - a_s, wb - b_s)
            if ov < min_ov:
                continue
            sub_a = a[a_s : a_s + ov].ravel()
            sub_b = cand[b_s : b_s + ov].ravel()
            sa = float(sub_a.std())
            sb = float(sub_b.std())
            if sa > 1e-12 and sb > 1e-12:
                r = float(np.corrcoef(sub_a, sub_b)[0, 1])
                if r > best:
                    best = r
    return round(max(0.0, best), 6)


# --------------------------------------------------------------------------- #
# 4. Weighted GSEA leading-edge extraction & deterministic AUCell scoring     #
# --------------------------------------------------------------------------- #


def _extract_weighted_gsea_leading_edge(
    gene_scores: np.ndarray,
    target_indices: Sequence[int],
    *,
    positive_direction: bool,
    min_abs_corr: float = 0.05,
) -> tuple[float, list[int]]:
    """Deterministic weighted-KS leading-edge extraction along signed `r_tf_gene`.

    Used strictly to prune trailing low-coexpression targets inside a candidate
    module (`Concilio d9842f` Blocking #2); module significance comes from the
    `(mean_expr, cells_detected)` bin-stratified null.
    """
    signed = gene_scores if positive_direction else -gene_scores
    order = np.argsort(-signed, kind="mergesort")
    n_total = len(signed)
    in_set = np.zeros(n_total, dtype=bool)
    in_set[list(target_indices)] = True
    hits = in_set[order]
    n_hits = int(hits.sum())
    n_miss = n_total - n_hits
    if n_hits == 0 or n_miss == 0:
        return 0.0, []

    weights = np.abs(signed[order])
    hit_w = np.where(hits, weights, 0.0)
    sum_w = float(hit_w.sum())
    if sum_w <= 1e-12:
        hit_step = hits.astype(np.float64) / float(n_hits)
    else:
        hit_step = hit_w / sum_w
    miss_step = (~hits).astype(np.float64) / float(n_miss)
    running = np.cumsum(hit_step - miss_step)
    peak_pos = int(np.argmax(running))
    es = float(running[peak_pos])

    le_mask = np.zeros(n_total, dtype=bool)
    le_mask[: peak_pos + 1] = True
    le_genes = [
        int(order[k])
        for k in range(peak_pos + 1)
        if hits[k] and signed[order[k]] >= float(min_abs_corr)
    ]
    return es, le_genes


def compute_aucell_scores(
    matrix: Any,
    feature_names: Sequence[str],
    feature_sets: Sequence[Sequence[str]],
    *,
    top_frac: float = 0.05,
    min_max_rank: int = 10,
) -> np.ndarray:
    """Deterministic tie-aware AUCell recovery-curve scores `(n_cells, n_sets)`.

    Computes the normalized area under the cumulative recovery curve in the top
    `max(min_max_rank, round(top_frac * n_features))` highest-valued features of
    each cell, handling sparse zero ties analytically so unexpressed features
    cannot arbitrarily enter the top-k window by column index!
    """
    n_cells = int(matrix.shape[0])
    n_feats = int(matrix.shape[1])
    n_sets = len(feature_sets)
    if n_cells == 0 or n_feats == 0 or n_sets == 0:
        return np.zeros((n_cells, n_sets), dtype=np.float32)

    name_to_idx = {str(f): i for i, f in enumerate(feature_names)}
    indexed_sets: list[np.ndarray] = []
    for fset in feature_sets:
        idxs = sorted({name_to_idx[str(g)] for g in fset if str(g) in name_to_idx})
        indexed_sets.append(np.asarray(idxs, dtype=np.int64))

    max_rank = min(n_feats, max(int(min_max_rank), int(round(top_frac * n_feats))))
    out = np.zeros((n_cells, n_sets), dtype=np.float32)
    csr = matrix.tocsr() if sparse.issparse(matrix) else None

    for c in range(n_cells):
        if csr is not None:
            start, end = csr.indptr[c], csr.indptr[c + 1]
            nz_cols = csr.indices[start:end]
            nz_vals = np.asarray(csr.data[start:end], dtype=np.float64)
            pos_mask = nz_vals > 0
            nz_cols = nz_cols[pos_mask]
            nz_vals = nz_vals[pos_mask]
        else:
            row = np.asarray(matrix[c], dtype=np.float64).ravel()
            nz_cols = np.flatnonzero(row > 0)
            nz_vals = row[nz_cols]

        n_pos = len(nz_cols)
        if n_pos == 0:
            continue

        # Sort positive features descending; average ranks across exact value ties
        order = np.argsort(-nz_vals, kind="mergesort")
        s_cols = nz_cols[order]
        s_vals = nz_vals[order]
        ranks = np.arange(n_pos, dtype=np.float64)  # 0-based (0 = highest)
        i = 0
        while i < n_pos:
            j = i
            while j + 1 < n_pos and s_vals[j + 1] == s_vals[i]:
                j += 1
            if j > i:
                ranks[i : j + 1] = 0.5 * (i + j)
            i = j + 1

        effective_k = min(max_rank, n_pos)
        keep = ranks < float(max_rank)
        if not np.any(keep):
            continue
        top_cols = s_cols[keep]
        top_ranks = ranks[keep]
        # Contribution of a hit at 0-based rank r in [0, max_rank): (max_rank - r)
        # Normalized by maximum attainable area: sum_{r=0..min(|S|, max_rank)-1} (max_rank - r)
        weights_by_col = np.zeros(n_feats, dtype=np.float64)
        weights_by_col[top_cols] = np.maximum(0.0, float(max_rank) - top_ranks)

        for s_idx, f_idxs in enumerate(indexed_sets):
            m = len(f_idxs)
            if m == 0:
                continue
            m_cap = min(m, max_rank)
            max_auc = m_cap * float(max_rank) - 0.5 * m_cap * (m_cap - 1)
            if max_auc <= 0:
                continue
            auc_val = float(weights_by_col[f_idxs].sum()) / max_auc
            # Penalize modestly if cell has fewer positive features than effective_k? No:
            # standard AUCell uses fixed `max_rank` ceiling and excludes zero-count ties.
            _ = effective_k
            out[c, s_idx] = np.float32(min(1.0, max(0.0, auc_val)))

    return out


# --------------------------------------------------------------------------- #
# 5. End-to-End Candidate eRegulon Computation                                #
# --------------------------------------------------------------------------- #


def _zscore_columns_and_depth_residual(
    dense_mat: np.ndarray,
    log_depth: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return `(Z_raw, Z_depth_resid, valid_nonzero_var_mask)` for `(n_cells, n_features)`.

    `Z_raw[:, j]` has zero mean and unit L2 norm (`Z_raw.T @ Z_raw` = Pearson `r`).
    `Z_depth_resid[:, j]` orthogonalizes each column against centered `log_depth`
    and re-normalizes to unit L2 norm (`Z_a.T @ Z_b` = partial Pearson `r | depth`).
    """
    n_cells = dense_mat.shape[0]
    centered = dense_mat - dense_mat.mean(axis=0, keepdims=True)
    norms = np.linalg.norm(centered, axis=0)
    valid = norms > 1e-12
    safe_norms = np.where(valid, norms, 1.0)
    z_raw = centered / safe_norms[None, :]
    z_raw[:, ~valid] = 0.0

    d_c = np.asarray(log_depth, dtype=np.float64).ravel()
    if len(d_c) == n_cells:
        d_c = d_c - d_c.mean()
        d_norm = float(np.linalg.norm(d_c))
    else:
        d_norm = 0.0

    if d_norm > 1e-12:
        u_d = d_c / d_norm
        proj = u_d @ centered  # (n_features,)
        resid = centered - np.outer(u_d, proj)
        r_norms = np.linalg.norm(resid, axis=0)
        r_valid = valid & (r_norms > 1e-12)
        safe_rn = np.where(r_valid, r_norms, 1.0)
        z_resid = resid / safe_rn[None, :]
        z_resid[:, ~r_valid] = 0.0
    else:
        z_resid = z_raw.copy()

    return z_raw, z_resid, valid


def _assign_2d_quantile_bins(
    cov1: np.ndarray,
    cov2: np.ndarray,
    *,
    n_bins_1: int = 5,
    n_bins_2: int = 5,
) -> np.ndarray:
    q1 = np.quantile(cov1, np.linspace(0.0, 1.0, n_bins_1 + 1))
    q2 = np.quantile(cov2, np.linspace(0.0, 1.0, n_bins_2 + 1))
    q1[0] -= 1e-9
    q1[-1] += 1e-9
    q2[0] -= 1e-9
    q2[-1] += 1e-9
    b1 = np.clip(np.digitize(cov1, q1[1:-1]), 0, n_bins_1 - 1)
    b2 = np.clip(np.digitize(cov2, q2[1:-1]), 0, n_bins_2 - 1)
    return (b1 * n_bins_2 + b2).astype(np.int32)


def compute_candidate_eregulons(
    rna_matrix: Any,
    rna_gene_names: Sequence[str],
    atac_matrix: Any,
    atac_peak_names: Sequence[str],
    gene_tss_records: dict[str, Any],
    *,
    meme_path: str | Path,
    fasta_path: str | Path,
    distance_bp: int = 500_000,
    min_distance_bp: int = 0,
    promoter_bp: int = 2_000,
    p2g_min_corr: float = 0.05,
    p2g_max_padj: float = 0.05,
    min_cells_detected: int = 10,
    scan_window_bp: int = 500,
    motif_pvalue: float = 5e-5,
    scanner_backend: str = "auto",
    min_tf_gene_corr: float = 0.05,
    min_target_genes: int = 5,
    max_regulon_padj: float = 0.05,
    n_null_permutations: int = 500,
    sign_classes: Sequence[str] = ("pos_tf_pos_p2g", "neg_tf_pos_p2g"),
    compute_cell_activity: bool = True,
    aucell_top_frac: float = 0.05,
    standard_chroms_only: bool = True,
    seed: int = 0,
    precomputed_p2g: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble candidate eRegulons (`TF -> cis peak -> target gene`) deterministically."""
    n_cells = int(rna_matrix.shape[0])
    if int(atac_matrix.shape[0]) != n_cells:
        raise ValueError(
            f"RNA and ATAC cell counts must match: {n_cells} vs {atac_matrix.shape[0]}"
        )

    # ── Step 1: Cis peak-to-gene links (reusing Slice I5 `peak2gene_core`) ──
    if precomputed_p2g is not None:
        p2g_res = precomputed_p2g
    else:
        p2g_res = compute_peak_to_gene_links(
            rna_matrix=rna_matrix,
            rna_var_names=rna_gene_names,
            atac_matrix=atac_matrix,
            atac_var_names=atac_peak_names,
            gene_records=gene_tss_records,
            distance_bp=int(distance_bp),
            min_distance_bp=int(min_distance_bp),
            promoter_bp=int(promoter_bp),
            min_corr=float(p2g_min_corr),
            max_padj=float(p2g_max_padj),
            rna_is_raw_counts=True,
            atac_is_raw_counts=True,
            preprocess="paired",
            atac_transform="raw_counts",
            min_cells_detected=int(min_cells_detected),
            background_null=True,
            n_bg_per_bin=200,
            seed=int(seed),
            genome_fasta=str(fasta_path),
            standard_chroms_only=bool(standard_chroms_only),
        )
        if p2g_res.get("status") == "error":
            if p2g_res.get("error_type") == "ContigMismatch":
                raise ContigMismatchError(
                    str(p2g_res.get("details") or "Contig mismatch between ATAC peaks and GTF.")
                )
            raise ValueError(
                str(p2g_res.get("details") or f"Peak-to-gene step failed: {p2g_res}")
            )

    cis_links: list[dict[str, Any]] = list(p2g_res.get("links") or [])

    # ── Step 2: Normalize RNA (`log1p(CP10k)`) & detect expressed genes/TFs ──
    if sparse.issparse(rna_matrix):
        rna_csc = rna_matrix.tocsc()
        gene_detected = np.diff(rna_csc.indptr)
        rna_lib = np.asarray(rna_csc.sum(axis=1)).ravel().astype(np.float64)
    else:
        rna_arr_in = np.asarray(rna_matrix, dtype=np.float64)
        gene_detected = np.sum(rna_arr_in > 0, axis=0)
        rna_lib = rna_arr_in.sum(axis=1)

    if sparse.issparse(atac_matrix):
        atac_csc = atac_matrix.tocsc()
        peak_detected = np.diff(atac_csc.indptr)
        atac_lib = np.asarray(atac_csc.sum(axis=1)).ravel().astype(np.float64)
    else:
        atac_arr_in = np.asarray(atac_matrix, dtype=np.float64)
        peak_detected = np.sum(atac_arr_in > 0, axis=0)
        atac_lib = atac_arr_in.sum(axis=1)

    log_rna_depth = np.log10(1.0 + np.maximum(rna_lib, 0.0))
    log_atac_depth = np.log10(1.0 + np.maximum(atac_lib, 0.0))

    motifs = parse_meme_motifs(meme_path)
    resolved_tfs, unmapped_motifs = resolve_motifs_to_expressed_tfs(
        motifs,
        rna_gene_names,
        gene_detected,
        min_cells_detected=int(min_cells_detected),
    )

    if not cis_links or not resolved_tfs:
        return {
            "eregulons": [],
            "triplets": [],
            "n_eregulons": 0,
            "n_triplets": 0,
            "n_unique_tfs": 0,
            "n_unique_enhancers": 0,
            "n_unique_target_genes": 0,
            "n_cis_links_input": len(cis_links),
            "n_motifs_total": len(motifs),
            "n_motif_tf_mappings": len(resolved_tfs),
            "n_unmapped_motifs": len(unmapped_motifs),
            "unmapped_motifs_sample": unmapped_motifs[:25],
            "n_candidate_modules_tested": 0,
            "scanner_backend": "not_run",
            "scan_window_bp": int(scan_window_bp),
            "motif_pvalue": float(motif_pvalue),
            "p2g_summary": {
                "n_tested": p2g_res.get("n_tested", 0),
                "n_fdr_significant": p2g_res.get(
                    "n_fdr_significant_all_tested",
                    p2g_res.get("n_fdr_significant", 0),
                ),
                "n_links_used": len(cis_links),
                "background_match": p2g_res.get("background_match"),
            },
            "rna_aucell": np.zeros((n_cells, 0), dtype=np.float32),
            "atac_aucell": np.zeros((n_cells, 0), dtype=np.float32),
        }

    # ── Step 3: Fixed-width 500bp peak-center sequences & motif scanning ──
    # Tested peak universe = detected standard peaks (>= min_cells_detected)
    tested_peak_indices: list[int] = []
    tested_peak_names: list[str] = []
    for p_idx, pname in enumerate(atac_peak_names):
        if int(peak_detected[p_idx]) < int(min_cells_detected):
            continue
        parsed = parse_peak_coord(pname)
        if parsed is None:
            continue
        if standard_chroms_only:
            from aria.utils.peak2gene_core import is_standard_contig

            if not is_standard_contig(parsed[0]):
                continue
        tested_peak_indices.append(p_idx)
        tested_peak_names.append(str(pname))

    seq_bundle = fetch_peak_center_sequences_from_indexed_fasta(
        tested_peak_names,
        fasta_path,
        scan_window_bp=int(scan_window_bp),
    )
    valid_seq_mask = seq_bundle["valid_mask"]
    if int(valid_seq_mask.sum()) == 0:
        raise ContigMismatchError(
            "Zero accessible peaks yielded valid sequences from reference FASTA."
        )

    valid_tested_positions = np.flatnonzero(valid_seq_mask)
    scanned_peak_names = [tested_peak_names[i] for i in valid_tested_positions]
    scanned_atac_col_indices = [tested_peak_indices[i] for i in valid_tested_positions]
    scanned_seqs = [str(seq_bundle["sequences"][i]) for i in valid_tested_positions]
    scanned_gc = seq_bundle["gc_fraction"][valid_tested_positions]
    scanned_peak_to_local = {p: i for i, p in enumerate(scanned_peak_names)}

    unique_motif_indices = sorted({r.motif_idx for r in resolved_tfs})
    motif_idx_to_scan_row = {m_i: r_i for r_i, m_i in enumerate(unique_motif_indices)}

    motif_hits_bool, backend_used = scan_motifs_in_sequences(
        meme_path=meme_path,
        motifs=motifs,
        sequences=scanned_seqs,
        motif_Indices=unique_motif_indices,
        pvalue=float(motif_pvalue),
        rc=True,
        backend=scanner_backend,
    )

    # GC-stratified CMH enrichment of motifs in cis-linked peaks vs tested universe
    cis_peak_set = {str(lk["peak"]) for lk in cis_links}
    fg_peak_mask = np.asarray(
        [p in cis_peak_set for p in scanned_peak_names], dtype=bool
    )
    cistrome_enr = compute_gc_stratified_motif_enrichment(
        motif_hits_bool, fg_peak_mask, scanned_gc
    )
    cistrome_padj = np.asarray(
        bh_correct(cistrome_enr["pvalue"].tolist()), dtype=np.float64
    )

    # ── Step 4: Expressed RNA standardization, depth residuals & 2D null bins ──
    expr_gene_col_indices = np.flatnonzero(gene_detected >= int(min_cells_detected))
    if sparse.issparse(rna_matrix):
        sub_rna = rna_csc[:, expr_gene_col_indices].toarray().astype(np.float64)
    else:
        sub_rna = np.asarray(
            rna_matrix[:, expr_gene_col_indices], dtype=np.float64
        )
    safe_lib = np.where(rna_lib > 0, rna_lib, 1.0)
    log_cp10k = np.log1p((sub_rna / safe_lib[:, None]) * 1e4)
    z_rna, z_rna_resid, valid_var_mask = _zscore_columns_and_depth_residual(
        log_cp10k, log_rna_depth
    )

    # Restrict expressed gene universe strictly to non-zero-variance genes
    keep_expr = np.flatnonzero(valid_var_mask)
    expr_gene_col_indices = expr_gene_col_indices[keep_expr]
    z_rna = z_rna[:, keep_expr]
    z_rna_resid = z_rna_resid[:, keep_expr]
    log_cp10k = log_cp10k[:, keep_expr]
    expr_gene_names = [str(rna_gene_names[i]) for i in expr_gene_col_indices]
    expr_gene_to_local = {g: i for i, g in enumerate(expr_gene_names)}

    gene_mean_expr = log_cp10k.mean(axis=0)
    gene_log_det = np.log10(
        1.0 + gene_detected[expr_gene_col_indices].astype(np.float64)
    )
    gene_null_bins = _assign_2d_quantile_bins(
        gene_mean_expr, gene_log_det, n_bins_1=5, n_bins_2=5
    )

    # Index cis-linked peaks for fast `r_tf_peak` & `r_tf_peak_partial_depth`
    cis_local_scan_indices = sorted(
        {
            scanned_peak_to_local[str(lk["peak"])]
            for lk in cis_links
            if str(lk["peak"]) in scanned_peak_to_local
        }
    )
    cis_scan_to_sub = {s_i: k for k, s_i in enumerate(cis_local_scan_indices)}
    cis_atac_cols = [scanned_atac_col_indices[s_i] for s_i in cis_local_scan_indices]
    if cis_atac_cols:
        if sparse.issparse(atac_matrix):
            sub_atac = atac_csc[:, cis_atac_cols].toarray().astype(np.float64)
        else:
            sub_atac = np.asarray(atac_matrix[:, cis_atac_cols], dtype=np.float64)
        # Condition ATAC peaks on ATAC depth (and RNA depth for cross-modal partial r)
        z_atac, z_atac_resid, _ = _zscore_columns_and_depth_residual(
            sub_atac, log_atac_depth
        )
    else:
        z_atac = np.zeros((n_cells, 0), dtype=np.float64)
        z_atac_resid = np.zeros((n_cells, 0), dtype=np.float64)

    # Index cis links by scanned peak local index
    links_by_scan_peak: dict[int, list[dict[str, Any]]] = {}
    for lk in cis_links:
        p_str = str(lk["peak"])
        g_str = str(lk["gene"])
        if p_str not in scanned_peak_to_local or g_str not in expr_gene_to_local:
            continue
        s_idx = scanned_peak_to_local[p_str]
        links_by_scan_peak.setdefault(s_idx, []).append(lk)

    # Precompute TF correlation vectors once per unique expressed TF
    unique_tfs = sorted(
        {r.tf_gene for r in resolved_tfs if r.tf_gene in expr_gene_to_local}
    )
    tf_corr_cache: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    if unique_tfs:
        tf_locs = [expr_gene_to_local[tf] for tf in unique_tfs]
        z_tf = z_rna[:, tf_locs]  # (n_cells, n_tfs)
        z_tf_res = z_rna_resid[:, tf_locs]
        r_tf_gene_all = np.clip(z_tf.T @ z_rna, -1.0, 1.0)  # (n_tfs, n_expr_genes)
        r_tf_gene_part = np.clip(z_tf_res.T @ z_rna_resid, -1.0, 1.0)
        if z_atac.shape[1] > 0:
            r_tf_peak_all = np.clip(z_tf.T @ z_atac, -1.0, 1.0)
            r_tf_peak_part = np.clip(z_tf_res.T @ z_atac_resid, -1.0, 1.0)
        else:
            r_tf_peak_all = np.zeros((len(unique_tfs), 0), dtype=np.float64)
            r_tf_peak_part = np.zeros((len(unique_tfs), 0), dtype=np.float64)
        for k, tf in enumerate(unique_tfs):
            tf_corr_cache[tf] = (
                r_tf_gene_all[k],
                r_tf_gene_part[k],
                r_tf_peak_all[k],
                r_tf_peak_part[k],
            )

    # ── Step 5: Evaluate EVERY `(motif_id, tf, sign_class)` candidate module ──
    allowed_signs = tuple(s for s in sign_classes if s in _SIGN_CLASSES)
    rng = np.random.default_rng(int(seed))
    candidate_modules: list[dict[str, Any]] = []

    for r_tf in resolved_tfs:
        if r_tf.tf_gene not in tf_corr_cache:
            continue
        scan_row = motif_idx_to_scan_row[r_tf.motif_idx]
        hit_scan_peaks = np.flatnonzero(motif_hits_bool[scan_row])
        if len(hit_scan_peaks) == 0:
            continue

        r_g, r_g_part, r_p, r_p_part = tf_corr_cache[r_tf.tf_gene]
        tf_loc = expr_gene_to_local[r_tf.tf_gene]
        tf_upper = r_tf.tf_gene.upper()

        # Partition motif-positive cis links into positive-p2g vs negative-p2g
        pos_p2g_by_gene: dict[int, list[tuple[int, dict[str, Any]]]] = {}
        neg_p2g_by_gene: dict[int, list[tuple[int, dict[str, Any]]]] = {}
        self_peaks: list[str] = []

        for s_idx in hit_scan_peaks:
            lks = links_by_scan_peak.get(int(s_idx))
            if not lks:
                continue
            for lk in lks:
                g_name = str(lk["gene"])
                if g_name.upper() == tf_upper:
                    self_peaks.append(str(lk["peak"]))
                    continue
                g_loc = expr_gene_to_local[g_name]
                if float(lk["correlation"]) >= 0.0:
                    pos_p2g_by_gene.setdefault(g_loc, []).append((int(s_idx), lk))
                else:
                    neg_p2g_by_gene.setdefault(g_loc, []).append((int(s_idx), lk))

        # Background pool strictly excludes self (`g == TF`)
        bg_mask = np.ones(len(expr_gene_names), dtype=bool)
        bg_mask[tf_loc] = False
        bg_indices = np.flatnonzero(bg_mask)
        if len(bg_indices) < max(20, int(min_target_genes) * 2):
            continue

        bg_scores = r_g[bg_indices]
        bg_bins = gene_null_bins[bg_indices]
        orig_to_bg = np.full(len(expr_gene_names), -1, dtype=np.int64)
        orig_to_bg[bg_indices] = np.arange(len(bg_indices), dtype=np.int64)

        # Precompute per-bin pools & moments for this TF's `r_tf_gene` vector
        bin_pools: dict[int, np.ndarray] = {}
        bin_mean: dict[int, float] = {}
        bin_var: dict[int, float] = {}
        for b_id in np.unique(bg_bins):
            b_idx = np.flatnonzero(bg_bins == b_id)
            vals = bg_scores[b_idx]
            bin_pools[int(b_id)] = vals
            bin_mean[int(b_id)] = float(vals.mean())
            bin_var[int(b_id)] = (
                float(vals.var(ddof=1)) if len(vals) > 1 else 0.0
            )
        global_var = float(bg_scores.var(ddof=1)) if len(bg_scores) > 1 else 1e-6

        for s_class in allowed_signs:
            use_pos_p2g = s_class.endswith("_pos_p2g")
            pos_tf_dir = s_class.startswith("pos_tf_")
            gene_map = pos_p2g_by_gene if use_pos_p2g else neg_p2g_by_gene
            if len(gene_map) < int(min_target_genes):
                continue

            cand_orig_locs = sorted(gene_map.keys())
            cand_bg_locs = [int(orig_to_bg[g_i]) for g_i in cand_orig_locs if orig_to_bg[g_i] >= 0]
            if len(cand_bg_locs) < int(min_target_genes):
                continue

            signed_bg = bg_scores if pos_tf_dir else -bg_scores
            obs_sum = float(signed_bg[cand_bg_locs].sum())

            # Covariate-matched finite-population stratified null (`(mean_expr, log_det)`)
            cand_bins = bg_bins[cand_bg_locs]
            uniq_b, counts_b = np.unique(cand_bins, return_counts=True)
            exp_sum = 0.0
            var_sum = 0.0
            mc_sums = np.zeros(int(n_null_permutations), dtype=np.float64)

            for b_id_raw, k_raw in zip(uniq_b, counts_b):
                b_id = int(b_id_raw)
                k_b = int(k_raw)
                pool = bin_pools[b_id]
                signed_pool = pool if pos_tf_dir else -pool
                n_b = len(signed_pool)
                mu_b = float(signed_pool.mean())
                exp_sum += k_b * mu_b
                if n_b > k_b and n_b > 1:
                    fpc = 1.0 - (float(k_b) / float(n_b))
                    v_b = float(signed_pool.var(ddof=1))
                    var_sum += k_b * fpc * v_b
                    if n_null_permutations > 0:
                        n_mc = int(n_null_permutations)
                        if k_b == 1:
                            mc_sums += signed_pool[rng.integers(0, n_b, size=n_mc)]
                        else:
                            # Exact independent without-replacement subsets via random keys
                            pick_idx = np.argpartition(
                                rng.random((n_mc, n_b)), k_b - 1, axis=1
                            )[:, :k_b]
                            mc_sums += signed_pool[pick_idx].sum(axis=1)
                else:
                    # Tiny stratum fallback: draw from full background
                    var_sum += k_b * global_var
                    if n_null_permutations > 0:
                        n_mc = int(n_null_permutations)
                        n_bg_tot = len(signed_bg)
                        if k_b == 1:
                            mc_sums += signed_bg[rng.integers(0, n_bg_tot, size=n_mc)]
                        elif k_b < n_bg_tot:
                            pick_idx = np.argpartition(
                                rng.random((n_mc, n_bg_tot)), k_b - 1, axis=1
                            )[:, :k_b]
                            mc_sums += signed_bg[pick_idx].sum(axis=1)

            z_mod = (
                (obs_sum - exp_sum) / math.sqrt(max(var_sum, 1e-12))
                if var_sum > 1e-12
                else 0.0
            )
            p_analytic = float(stats.norm.sf(z_mod))
            if n_null_permutations > 0:
                n_ge = int(np.sum(mc_sums >= obs_sum - 1e-12))
                p_emp = float(n_ge + 1) / float(int(n_null_permutations) + 1)
                # Combine analytic tail (captures extreme z > 4 beyond MC resolution)
                # while preventing anti-conservative normal tail when MC disagrees:
                if n_ge == 0:
                    p_mod = p_analytic
                else:
                    p_mod = max(p_analytic, p_emp)
            else:
                p_mod = p_analytic

            # Weighted GSEA leading-edge extraction (pruning only, NOT a second test)
            gsea_es, le_bg_locs = _extract_weighted_gsea_leading_edge(
                bg_scores,
                cand_bg_locs,
                positive_direction=pos_tf_dir,
                min_abs_corr=float(min_tf_gene_corr),
            )
            le_orig_locs = [int(bg_indices[b_l]) for b_l in le_bg_locs]

            candidate_modules.append(
                {
                    "resolved": r_tf,
                    "scan_row": scan_row,
                    "sign_class": s_class,
                    "n_candidate_genes": len(cand_bg_locs),
                    "module_zscore": round(float(z_mod), 6),
                    "module_pvalue": float(p_mod),
                    "gsea_es": round(float(gsea_es), 6),
                    "le_orig_locs": le_orig_locs,
                    "gene_map": gene_map,
                    "self_peaks": sorted(set(self_peaks)),
                    "r_g": r_g,
                    "r_g_part": r_g_part,
                    "r_p": r_p,
                    "r_p_part": r_p_part,
                }
            )

    if not candidate_modules:
        return {
            "eregulons": [],
            "triplets": [],
            "n_eregulons": 0,
            "n_triplets": 0,
            "n_unique_tfs": 0,
            "n_unique_enhancers": 0,
            "n_unique_target_genes": 0,
            "n_cis_links_input": len(cis_links),
            "n_motifs_total": len(motifs),
            "n_motif_tf_mappings": len(resolved_tfs),
            "n_unmapped_motifs": len(unmapped_motifs),
            "unmapped_motifs_sample": unmapped_motifs[:25],
            "n_candidate_modules_tested": 0,
            "scanner_backend": backend_used,
            "scan_window_bp": int(seq_bundle["scan_window_bp"]),
            "motif_pvalue": float(motif_pvalue),
            "p2g_summary": {
                "n_tested": p2g_res.get("n_tested", 0),
                "n_fdr_significant": p2g_res.get(
                    "n_fdr_significant_all_tested",
                    p2g_res.get("n_fdr_significant", 0),
                ),
                "n_links_used": len(cis_links),
                "background_match": p2g_res.get("background_match"),
            },
            "rna_aucell": np.zeros((n_cells, 0), dtype=np.float32),
            "atac_aucell": np.zeros((n_cells, 0), dtype=np.float32),
        }

    # ── Step 6: Module-level BH-FDR BEFORE collapsing multiple motifs per TF ──
    raw_pvals = [m["module_pvalue"] for m in candidate_modules]
    mod_padjs = bh_correct(raw_pvals)
    for m, padj in zip(candidate_modules, mod_padjs):
        m["regulon_padj"] = float(padj)

    # Count how many motifs were tested per `(tf_gene, sign_class)`
    motifs_per_tf_sign: dict[tuple[str, str], list[str]] = {}
    for m in candidate_modules:
        key = (m["resolved"].tf_gene, m["sign_class"])
        motifs_per_tf_sign.setdefault(key, []).append(m["resolved"].motif_id)

    # Filter passing modules (`regulon_padj <= max_regulon_padj` & `len(LE) >= min_target_genes`)
    # and keep the single strongest motif per `(tf_gene, sign_class)` (no union inflation!)
    best_by_tf_sign: dict[tuple[str, str], dict[str, Any]] = {}
    passing_motifs_by_tf_sign: dict[tuple[str, str], list[str]] = {}

    for m in candidate_modules:
        if m["regulon_padj"] > float(max_regulon_padj):
            continue
        if len(m["le_orig_locs"]) < int(min_target_genes):
            continue
        key = (m["resolved"].tf_gene, m["sign_class"])
        passing_motifs_by_tf_sign.setdefault(key, []).append(m["resolved"].motif_id)
        prev = best_by_tf_sign.get(key)
        if prev is None or (
            m["regulon_padj"],
            -m["module_zscore"],
            m["resolved"].motif_id,
        ) < (
            prev["regulon_padj"],
            -prev["module_zscore"],
            prev["resolved"].motif_id,
        ):
            best_by_tf_sign[key] = m

    # Sort surviving eRegulons deterministically
    selected_keys = sorted(
        best_by_tf_sign.keys(),
        key=lambda k: (
            best_by_tf_sign[k]["regulon_padj"],
            -best_by_tf_sign[k]["module_zscore"],
            k[0],
            k[1],
        ),
    )

    eregulons: list[dict[str, Any]] = []
    all_triplets: list[dict[str, Any]] = []
    peak_sets_for_ambiguity: list[set[str]] = []
    pfms_for_ambiguity: list[np.ndarray] = []

    for key in selected_keys:
        m = best_by_tf_sign[key]
        r_tf: ResolvedMotifTF = m["resolved"]
        s_class: str = m["sign_class"]
        scan_row: int = m["scan_row"]
        regulon_id = f"{r_tf.tf_gene}__{s_class}"

        r_g = m["r_g"]
        r_g_part = m["r_g_part"]
        r_p = m["r_p"]
        r_p_part = m["r_p_part"]

        reg_triplets: list[dict[str, Any]] = []
        reg_genes: list[str] = []
        reg_peaks_set: set[str] = set()

        # Order leading-edge genes by `|r_tf_gene|` descending
        le_sorted = sorted(
            m["le_orig_locs"],
            key=lambda g_i: (-abs(float(r_g[g_i])), expr_gene_names[g_i]),
        )
        for g_loc in le_sorted:
            g_sym = expr_gene_names[g_loc]
            reg_genes.append(g_sym)
            peak_entries = sorted(
                m["gene_map"][g_loc],
                key=lambda pair: (
                    float(pair[1].get("padj") or 1.0),
                    -abs(float(pair[1]["correlation"])),
                    str(pair[1]["peak"]),
                ),
            )
            for s_idx, lk in peak_entries:
                p_name = str(lk["peak"])
                reg_peaks_set.add(p_name)
                sub_p_idx = cis_scan_to_sub.get(s_idx)
                tf_pk_r = float(r_p[sub_p_idx]) if sub_p_idx is not None else 0.0
                tf_pk_part = (
                    float(r_p_part[sub_p_idx]) if sub_p_idx is not None else 0.0
                )
                trip = {
                    "regulon_id": regulon_id,
                    "tf": r_tf.tf_gene,
                    "motif_id": r_tf.motif_id,
                    "motif_source": r_tf.motif_source,
                    "sign_class": s_class,
                    "peak": p_name,
                    "gene": g_sym,
                    "r_tf_gene": round(float(r_g[g_loc]), 6),
                    "r_tf_gene_partial_depth": round(float(r_g_part[g_loc]), 6),
                    "r_tf_peak": round(tf_pk_r, 6),
                    "r_tf_peak_partial_depth": round(tf_pk_part, 6),
                    "r_peak_gene": round(float(lk["correlation"]), 6),
                    "padj_peak_gene": (
                        float(lk["padj"]) if lk.get("padj") is not None else None
                    ),
                    "distance_bp": int(lk["distance_bp"]),
                    "region_class": str(lk.get("region_class") or "cis"),
                    "regulon_padj": float(m["regulon_padj"]),
                }
                reg_triplets.append(trip)
                all_triplets.append(trip)

        mean_raw_r = float(np.mean([r_g[g_i] for g_i in le_sorted]))
        mean_part_r = float(np.mean([r_g_part[g_i] for g_i in le_sorted]))
        depth_confounded = bool(
            abs(mean_part_r) < 0.5 * abs(mean_raw_r)
            and abs(mean_part_r) < float(min_tf_gene_corr)
        )

        eregulons.append(
            {
                "regulon_id": regulon_id,
                "tf": r_tf.tf_gene,
                "sign_class": s_class,
                "best_motif_id": r_tf.motif_id,
                "raw_motif_name": r_tf.raw_motif_name,
                "motif_source": r_tf.motif_source,
                "tf_symbol_match_mode": r_tf.match_mode,
                "heterodimer_partners": list(r_tf.heterodimer_partners),
                "n_motifs_tested_for_tf": len(motifs_per_tf_sign.get(key, [])),
                "supporting_motif_ids": sorted(
                    set(passing_motifs_by_tf_sign.get(key, []))
                ),
                "n_candidate_genes": int(m["n_candidate_genes"]),
                "n_target_genes": len(reg_genes),
                "n_enhancer_peaks": len(reg_peaks_set),
                "n_triplets": len(reg_triplets),
                "module_zscore": m["module_zscore"],
                "module_pvalue": m["module_pvalue"],
                "regulon_padj": float(m["regulon_padj"]),
                "regulon_padj_scope": "module_level",
                "gsea_leading_edge_es": m["gsea_es"],
                "mean_r_tf_gene": round(mean_raw_r, 6),
                "mean_r_tf_gene_partial_depth": round(mean_part_r, 6),
                "depth_confounded": depth_confounded,
                "cistrome_cmh_z": round(float(cistrome_enr["cmh_z"][scan_row]), 6),
                "cistrome_cmh_padj": float(cistrome_padj[scan_row]),
                "cistrome_log2_fe": round(
                    float(cistrome_enr["log2_fold_enrichment"][scan_row]), 6
                ),
                "self_linked_peaks": m["self_peaks"],
                "target_genes": reg_genes,
                "enhancer_peaks": sorted(reg_peaks_set),
                "motif_family_ambiguity": False,
                "ambiguous_paralog_tfs": [],
                "mean_rna_aucell": None,
                "mean_atac_aucell": None,
                "rna_atac_activity_corr": None,
            }
        )
        peak_sets_for_ambiguity.append(reg_peaks_set)
        pfms_for_ambiguity.append(motifs[r_tf.motif_idx].pfm)

    # ── Step 7: Paralog motif-sharing ambiguity screen (excluding dimer pairs) ──
    n_reg = len(eregulons)
    for i in range(n_reg):
        amb_list: list[dict[str, Any]] = []
        tf_i = eregulons[i]["tf"]
        partners_i = set(eregulons[i]["heterodimer_partners"])
        ps_i = peak_sets_for_ambiguity[i]
        if not ps_i:
            continue
        for j in range(n_reg):
            if i == j:
                continue
            tf_j = eregulons[j]["tf"]
            if tf_i == tf_j:
                continue
            # Exclude physical heterodimer partners assigned from the same dimer motif
            partners_j = set(eregulons[j]["heterodimer_partners"])
            if (tf_j in partners_i and tf_i in partners_j) and (
                eregulons[i]["best_motif_id"] == eregulons[j]["best_motif_id"]
            ):
                continue
            ps_j = peak_sets_for_ambiguity[j]
            union_n = len(ps_i | ps_j)
            if union_n == 0:
                continue
            p_jacc = float(len(ps_i & ps_j)) / float(union_n)
            if p_jacc < 0.30:
                continue
            pwm_sim = compute_max_aligned_pwm_similarity(
                pfms_for_ambiguity[i], pfms_for_ambiguity[j]
            )
            if p_jacc >= 0.50 or (p_jacc >= 0.30 and pwm_sim >= 0.80):
                amb_list.append(
                    {
                        "tf": tf_j,
                        "sign_class": eregulons[j]["sign_class"],
                        "peak_jaccard": round(p_jacc, 4),
                        "pwm_similarity": round(pwm_sim, 4),
                    }
                )
        if amb_list:
            amb_list.sort(key=lambda d: (-d["peak_jaccard"], d["tf"]))
            eregulons[i]["motif_family_ambiguity"] = True
            eregulons[i]["ambiguous_paralog_tfs"] = amb_list

    # ── Step 8: Per-cell deterministic AUCell RNA & ATAC regulon activities ──
    if compute_cell_activity and eregulons:
        rna_sets = [r["target_genes"] for r in eregulons]
        atac_sets = [r["enhancer_peaks"] for r in eregulons]
        rna_aucell = compute_aucell_scores(
            rna_matrix,
            rna_gene_names,
            rna_sets,
            top_frac=float(aucell_top_frac),
        )
        atac_aucell = compute_aucell_scores(
            atac_matrix,
            atac_peak_names,
            atac_sets,
            top_frac=float(aucell_top_frac),
        )
        for idx_r, reg in enumerate(eregulons):
            ra = rna_aucell[:, idx_r].astype(np.float64)
            aa = atac_aucell[:, idx_r].astype(np.float64)
            reg["mean_rna_aucell"] = round(float(ra.mean()), 6)
            reg["mean_atac_aucell"] = round(float(aa.mean()), 6)
            if ra.std() > 1e-12 and aa.std() > 1e-12:
                reg["rna_atac_activity_corr"] = round(
                    float(np.corrcoef(ra, aa)[0, 1]), 6
                )
            else:
                reg["rna_atac_activity_corr"] = None
    else:
        rna_aucell = np.zeros((n_cells, len(eregulons)), dtype=np.float32)
        atac_aucell = np.zeros((n_cells, len(eregulons)), dtype=np.float32)

    return {
        "eregulons": eregulons,
        "triplets": all_triplets,
        "n_eregulons": len(eregulons),
        "n_triplets": len(all_triplets),
        "n_unique_tfs": len({r["tf"] for r in eregulons}),
        "n_unique_enhancers": len({t["peak"] for t in all_triplets}),
        "n_unique_target_genes": len({t["gene"] for t in all_triplets}),
        "n_cis_links_input": len(cis_links),
        "n_motifs_total": len(motifs),
        "n_motif_tf_mappings": len(resolved_tfs),
        "n_unmapped_motifs": len(unmapped_motifs),
        "unmapped_motifs_sample": unmapped_motifs[:25],
        "n_candidate_modules_tested": len(candidate_modules),
        "scanner_backend": backend_used,
        "scan_window_bp": int(seq_bundle["scan_window_bp"]),
        "motif_pvalue": float(motif_pvalue),
        "p2g_summary": {
            "n_tested": p2g_res.get("n_tested", 0),
            "n_fdr_significant": p2g_res.get(
                "n_fdr_significant_all_tested",
                p2g_res.get("n_fdr_significant", 0),
            ),
            "n_links_used": len(cis_links),
            "background_match": p2g_res.get("background_match"),
        },
        "rna_aucell": rna_aucell,
        "atac_aucell": atac_aucell,
    }
