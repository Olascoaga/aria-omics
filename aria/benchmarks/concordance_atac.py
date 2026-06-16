"""P3a — scATAC external-concordance scoring (ARIA vs SnapATAC2/ArchR/Signac).

This is the PURE, dependency-light (numpy only) scoring layer for the external
concordance benchmark (master-plan W3.1 / scATAC P3a). It scores ARIA's scATAC
outputs against a reference tool's outputs on the metrics a reviewer expects:

- cluster concordance: Adjusted Rand Index (ARI) and Normalized Mutual
  Information (NMI) over the shared cells;
- DA-peak overlap: exact-string Jaccard AND genomic-overlap matching (the
  ADR-043 lesson — boundary-shifted peaks that overlap in genomic space must
  count as the same peak, which exact-string overlap misses);
- motif concordance: top-k Jaccard + rank-biased overlap (RBO);
- seed stability: mean pairwise ARI across repeated runs.

No-fabrication contract (ADR-002): ``score_atac_concordance`` returns
``status="not_run"`` with a reason when the external tool's outputs are absent —
it never emits a perfect-but-fake concordance. ARI/NMI are implemented in pure
numpy so this scores anywhere, without sklearn.

P3a is buildable now: it CONSUMES the external tool's outputs (clusters/DA
peaks/motifs as plain lists). The infra-gated P3b drop-in is the actual
SnapATAC2/ArchR/Signac driver that PRODUCES those outputs inside
``aria-bench-env``; it requires installing those tools + public datasets, which
is Samael's call and is not done here.

Neutral labels only (ADR-011): this module makes no biological claim from any
cluster id, peak coordinate, or motif name.
"""
from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from math import comb
from typing import Any, Sequence

import numpy as np

from aria.benchmarks.synthetic_de import _rank_biased_overlap


# --------------------------------------------------------------------------- #
# Cluster concordance: ARI / NMI (pure numpy)                                  #
# --------------------------------------------------------------------------- #

def _contingency(a: Sequence[Any], b: Sequence[Any]) -> np.ndarray:
    """Contingency (confusion) matrix of two label assignments of equal length."""
    a = np.asarray(list(a))
    b = np.asarray(list(b))
    if a.shape[0] != b.shape[0]:
        raise ValueError(
            f"label vectors must be the same length: {a.shape[0]} != {b.shape[0]}")
    _, ia = np.unique(a, return_inverse=True)
    _, ib = np.unique(b, return_inverse=True)
    table = np.zeros((ia.max() + 1, ib.max() + 1), dtype=np.int64)
    np.add.at(table, (ia, ib), 1)
    return table


def adjusted_rand_index(a: Sequence[Any], b: Sequence[Any]) -> float:
    """Adjusted Rand Index — chance-corrected partition agreement.

    1.0 == identical partitions (invariant to label renaming), ~0.0 == agreement
    no better than random, can be negative for worse-than-random.
    """
    table = _contingency(a, b)
    n = int(table.sum())
    if n < 2:
        return 1.0
    sum_comb_rows = sum(comb(int(x), 2) for x in table.sum(axis=1))
    sum_comb_cols = sum(comb(int(x), 2) for x in table.sum(axis=0))
    sum_comb_cells = sum(comb(int(x), 2) for x in table.ravel())
    total = comb(n, 2)
    expected = (sum_comb_rows * sum_comb_cols) / total
    max_index = 0.5 * (sum_comb_rows + sum_comb_cols)
    if max_index == expected:
        return 1.0
    return float((sum_comb_cells - expected) / (max_index - expected))


def normalized_mutual_info(a: Sequence[Any], b: Sequence[Any]) -> float:
    """Normalized Mutual Information (arithmetic-mean normalization), in [0, 1].

    1.0 == one partition fully determines the other; 0.0 == independent (one
    labeling carries no information about the other).
    """
    table = _contingency(a, b).astype(float)
    n = table.sum()
    if n == 0:
        return 1.0
    pij = table / n
    pi = pij.sum(axis=1)
    pj = pij.sum(axis=0)
    nz = table > 0
    outer = np.outer(pi, pj)
    mi = float(np.sum(pij[nz] * np.log(pij[nz] / outer[nz])))
    hi = -float(np.sum(pi[pi > 0] * np.log(pi[pi > 0])))
    hj = -float(np.sum(pj[pj > 0] * np.log(pj[pj > 0])))
    denom = 0.5 * (hi + hj)
    if denom == 0:
        # both labelings are single-cluster (no structure either way) -> agree
        return 1.0
    return max(0.0, min(1.0, mi / denom))


# --------------------------------------------------------------------------- #
# DA-peak overlap (exact + genomic)                                           #
# --------------------------------------------------------------------------- #

_PEAK_RE = re.compile(r"^(chr[\w]+)[:\-_](\d+)[\-_](\d+)$")


def _parse_peak(name: Any) -> tuple[str, int, int] | None:
    """Parse ``chr1:100-200`` / ``chr1-100-200`` / ``chr1_100_200`` -> interval."""
    m = _PEAK_RE.match(str(name).strip())
    if not m:
        return None
    start, end = int(m.group(2)), int(m.group(3))
    if end < start:
        start, end = end, start
    return (m.group(1), start, end)


def _count_genomic_matches(
    peaks_a: Sequence[Any], peaks_b: Sequence[Any], slack: int = 0
) -> tuple[int, int, int]:
    """Return (matched_a, matched_b, n_parsed_b) under interval overlap.

    A peak in A is "matched" if it overlaps (within ``slack`` bp) at least one
    peak in B, and vice versa. Boundary-shifted peaks that overlap are matched —
    this is exactly what exact-string comparison misses (ADR-043).
    """
    by_chrom_b: dict[str, list[tuple[int, int]]] = {}
    for p in peaks_b:
        parsed = _parse_peak(p)
        if parsed is None:
            continue
        chrom, s, e = parsed
        by_chrom_b.setdefault(chrom, []).append((s, e))
    for chrom in by_chrom_b:
        by_chrom_b[chrom].sort()
    n_parsed_b = sum(len(v) for v in by_chrom_b.values())

    def overlaps_any(chrom: str, s: int, e: int) -> bool:
        for bs, be in by_chrom_b.get(chrom, ()):  # small per-chrom lists in practice
            if bs - slack > e:
                break
            if be + slack >= s:
                return True
        return False

    matched_a = 0
    matched_b_set: set[tuple[str, int, int]] = set()
    for p in peaks_a:
        parsed = _parse_peak(p)
        if parsed is None:
            continue
        chrom, s, e = parsed
        if overlaps_any(chrom, s, e):
            matched_a += 1
    # symmetric pass for matched_b
    by_chrom_a: dict[str, list[tuple[int, int]]] = {}
    for p in peaks_a:
        parsed = _parse_peak(p)
        if parsed is None:
            continue
        chrom, s, e = parsed
        by_chrom_a.setdefault(chrom, []).append((s, e))
    for chrom in by_chrom_a:
        by_chrom_a[chrom].sort()

    def overlaps_any_a(chrom: str, s: int, e: int) -> bool:
        for as_, ae in by_chrom_a.get(chrom, ()):
            if as_ - slack > e:
                break
            if ae + slack >= s:
                return True
        return False

    for chrom, intervals in by_chrom_b.items():
        for s, e in intervals:
            if overlaps_any_a(chrom, s, e):
                matched_b_set.add((chrom, s, e))
    return matched_a, len(matched_b_set), n_parsed_b


def peak_set_overlap(
    peaks_a: Sequence[Any],
    peaks_b: Sequence[Any],
    *,
    mode: str = "genomic",
    slack: int = 0,
) -> dict[str, Any]:
    """Overlap between two DA-peak sets.

    ``mode="exact"``: Jaccard of the raw peak-id strings.
    ``mode="genomic"``: interval-overlap matching (boundary-shifted peaks count),
    reporting directional recovery (``matched_a``/``matched_b``) and a symmetric
    ``overlap_fraction = (matched_a + matched_b) / (n_a + n_b)``.
    """
    a = list(peaks_a)
    b = list(peaks_b)
    n_a, n_b = len(a), len(b)
    out: dict[str, Any] = {"mode": mode, "n_a": n_a, "n_b": n_b}

    set_a, set_b = set(map(str, a)), set(map(str, b))
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    out["jaccard"] = (inter / union) if union else 1.0

    if mode == "exact":
        out["intersection"] = inter
        out["union"] = union
        return out

    matched_a, matched_b, _ = _count_genomic_matches(a, b, slack=slack)
    out["slack"] = slack
    out["matched_a"] = matched_a
    out["matched_b"] = matched_b
    out["recall_a_in_b"] = (matched_a / n_a) if n_a else 1.0
    out["recall_b_in_a"] = (matched_b / n_b) if n_b else 1.0
    denom = n_a + n_b
    out["overlap_fraction"] = ((matched_a + matched_b) / denom) if denom else 1.0
    return out


# --------------------------------------------------------------------------- #
# Motif concordance                                                           #
# --------------------------------------------------------------------------- #

def motif_concordance(
    motifs_a: Sequence[Any],
    motifs_b: Sequence[Any],
    *,
    top_k: int | None = None,
) -> dict[str, Any]:
    """Concordance of two RANKED motif lists (most-enriched first).

    Reports top-k Jaccard and rank-biased overlap (RBO, reusing the DE benchmark
    helper). ``top_k=None`` uses the shorter of the two lists.
    """
    a = [str(m) for m in motifs_a]
    b = [str(m) for m in motifs_b]
    k = top_k if top_k is not None else min(len(a), len(b))
    k = max(k, 0)
    top_a, top_b = set(a[:k]), set(b[:k])
    union = len(top_a | top_b)
    jaccard = (len(top_a & top_b) / union) if union else 1.0
    rbo = _rank_biased_overlap(a[:k], b[:k], p=0.9) if k else 1.0
    return {
        "top_k": k,
        "n_a": len(a),
        "n_b": len(b),
        "jaccard_top_k": jaccard,
        "rank_biased_overlap": float(rbo),
    }


# --------------------------------------------------------------------------- #
# Seed stability                                                              #
# --------------------------------------------------------------------------- #

def seed_stability(label_sets: Sequence[Sequence[Any]]) -> dict[str, Any]:
    """Mean / min pairwise ARI across repeated clusterings (seed robustness)."""
    sets = [list(s) for s in label_sets]
    if len(sets) < 2:
        return {"n_seeds": len(sets), "mean_pairwise_ari": 1.0,
                "min_pairwise_ari": 1.0, "n_pairs": 0}
    aris: list[float] = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            aris.append(adjusted_rand_index(sets[i], sets[j]))
    return {
        "n_seeds": len(sets),
        "n_pairs": len(aris),
        "mean_pairwise_ari": float(np.mean(aris)),
        "min_pairwise_ari": float(np.min(aris)),
    }


# --------------------------------------------------------------------------- #
# Aggregator                                                                  #
# --------------------------------------------------------------------------- #

@dataclass
class ConcordanceResult:
    status: str
    tool: str
    cluster_concordance: dict[str, Any] | None = None
    da_peak_concordance: dict[str, Any] | None = None
    motif_concordance: dict[str, Any] | None = None
    seed_stability: dict[str, Any] | None = None
    reason: str | None = None
    messages: list[str] = field(default_factory=list)

    def to_manifest(self) -> dict[str, Any]:
        out = asdict(self)
        return {k: v for k, v in out.items() if v is not None or k in
                ("cluster_concordance", "da_peak_concordance", "motif_concordance")}


def score_atac_concordance(
    aria: dict[str, Any],
    external: dict[str, Any] | None,
    *,
    tool: str = "external",
    top_k_motifs: int | None = None,
    slack: int = 0,
) -> dict[str, Any]:
    """Score ARIA's scATAC outputs against a reference tool's outputs.

    Each output dict may carry ``clusters`` (per-cell labels, ALIGNED cell order
    between ARIA and external), ``da_peaks`` (peak ids/coords), and ``motifs``
    (ranked motif names). Any missing field is scored honestly as ``None`` for
    that axis — never fabricated.

    Returns ``status="not_run"`` with a reason when ``external`` is absent, so a
    missing comparator never reads as perfect agreement (ADR-002).
    """
    if not external:
        return ConcordanceResult(
            status="not_run", tool=tool,
            reason=(f"no {tool} outputs provided; the external comparator "
                    "(SnapATAC2/ArchR/Signac) must run in aria-bench-env (P3b) "
                    "before concordance can be scored"),
        ).to_manifest()

    res = ConcordanceResult(status="success", tool=tool)

    a_clusters = aria.get("clusters")
    b_clusters = external.get("clusters")
    if a_clusters is not None and b_clusters is not None and len(a_clusters):
        if len(a_clusters) != len(b_clusters):
            res.cluster_concordance = None
            res.messages.append(
                f"cluster vectors not aligned ({len(a_clusters)} vs "
                f"{len(b_clusters)} cells); concordance not scored")
        else:
            res.cluster_concordance = {
                "ari": adjusted_rand_index(a_clusters, b_clusters),
                "nmi": normalized_mutual_info(a_clusters, b_clusters),
                "n_cells": len(a_clusters),
            }

    a_peaks = aria.get("da_peaks")
    b_peaks = external.get("da_peaks")
    if a_peaks is not None and b_peaks is not None:
        res.da_peak_concordance = peak_set_overlap(
            a_peaks, b_peaks, mode="genomic", slack=slack)

    a_motifs = aria.get("motifs")
    b_motifs = external.get("motifs")
    if a_motifs is not None and b_motifs is not None:
        res.motif_concordance = motif_concordance(
            a_motifs, b_motifs, top_k=top_k_motifs)

    aria_seeds = aria.get("cluster_seeds")
    if aria_seeds:
        res.seed_stability = seed_stability(aria_seeds)

    return res.to_manifest()
