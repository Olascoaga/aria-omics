"""Consensus-by-genomic-overlap pseudobulk assembly for multi-sample scATAC DA.

When per-donor scATAC peak matrices are called independently they do NOT share a
peak feature space, and without fragment files they cannot be re-quantified onto a
common tile/peak set. This module builds a shared feature space the only way the
processed matrices allow: union-merge the donors' peak intervals into a consensus
set (the ADR-043 lesson — overlapping peaks are the same peak), then re-map each
donor's per-peak PSEUDOBULK counts (already summed over cells) onto the consensus
by peak-midpoint containment.

This is an approximation versus fragment-level re-quantification (boundary effects
in the count re-mapping) and is labelled as such wherever its outputs are scored.
Pure numpy; the peak parser is the single ``concordance_atac._parse_peak``.
"""
from __future__ import annotations

import bisect
from typing import Mapping, Sequence

import numpy as np

from aria.benchmarks.concordance_atac import _parse_peak


def merge_consensus_peaks(peak_lists: Sequence[Sequence[str]]) -> list[str]:
    """Union-merge overlapping peak intervals across donors into a consensus set.

    Returns sorted ``chrom:start-end`` ids. Unparseable peak names are skipped.
    """
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for peaks in peak_lists:
        for p in peaks:
            parsed = _parse_peak(p)
            if parsed is None:
                continue
            chrom, s, e = parsed
            by_chrom.setdefault(chrom, []).append((s, e))
    consensus: list[tuple[str, int, int]] = []
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        cs, ce = ivs[0]
        for s, e in ivs[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                consensus.append((chrom, cs, ce))
                cs, ce = s, e
        consensus.append((chrom, cs, ce))
    consensus.sort()
    return [f"{c}:{s}-{e}" for c, s, e in consensus]


def _consensus_index(consensus_ids: Sequence[str]) -> dict[str, tuple[list[int], list[int], list[int]]]:
    """Per-chrom sorted (starts, ends, row_index) for midpoint-containment lookup."""
    index: dict[str, list[tuple[int, int, int]]] = {}
    for i, cid in enumerate(consensus_ids):
        parsed = _parse_peak(cid)
        if parsed is None:
            continue
        chrom, s, e = parsed
        index.setdefault(chrom, []).append((s, e, i))
    out: dict[str, tuple[list[int], list[int], list[int]]] = {}
    for chrom, rows in index.items():
        rows.sort()
        out[chrom] = ([r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows])
    return out


def remap_to_consensus(
    peaks: Sequence[str],
    counts: Sequence[float],
    consensus_index: Mapping[str, tuple[list[int], list[int], list[int]]],
    n_consensus: int,
) -> np.ndarray:
    """Re-map one donor's per-peak counts onto consensus bins by peak midpoint.

    A donor peak's counts go to the consensus interval that CONTAINS its midpoint.
    Because the consensus is a union-merge of all donor peaks, each donor peak's
    midpoint falls inside exactly one consensus interval (or none, if the peak was
    unparseable and excluded from the merge).
    """
    vec = np.zeros(n_consensus, dtype=np.float64)
    for p, c in zip(peaks, counts):
        parsed = _parse_peak(p)
        if parsed is None:
            continue
        chrom, s, e = parsed
        entry = consensus_index.get(chrom)
        if entry is None:
            continue
        starts, ends, idx = entry
        mid = (s + e) // 2
        k = bisect.bisect_right(starts, mid) - 1
        if k >= 0 and starts[k] <= mid <= ends[k]:
            vec[idx[k]] += c
    return vec


def build_consensus_pseudobulk(
    donor_pseudobulks: Mapping[str, tuple[Sequence[str], Sequence[float]]],
) -> tuple[np.ndarray, list[str], list[str]]:
    """Assemble a (consensus_peaks x donors) pseudobulk count matrix.

    ``donor_pseudobulks`` maps donor -> (peak_ids, per_peak_pseudobulk_counts).
    Returns ``(matrix, consensus_ids, donor_order)``. ``donor_order`` is the input
    key order so callers can align condition labels.
    """
    donor_order = list(donor_pseudobulks.keys())
    consensus_ids = merge_consensus_peaks(
        [donor_pseudobulks[d][0] for d in donor_order])
    index = _consensus_index(consensus_ids)
    mat = np.zeros((len(consensus_ids), len(donor_order)), dtype=np.float64)
    for j, d in enumerate(donor_order):
        peaks, counts = donor_pseudobulks[d]
        mat[:, j] = remap_to_consensus(peaks, counts, index, len(consensus_ids))
    return mat, consensus_ids, donor_order
