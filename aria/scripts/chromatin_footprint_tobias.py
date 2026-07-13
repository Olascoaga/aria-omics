#!/usr/bin/env python3
"""scATAC P4.3 — Tn5-bias-corrected footprinting + differential TF binding (TOBIAS).

Runs in the dedicated ``aria-tobias-env`` (``envs/aria-tobias-env.yml``). This is the
backend the honest stub ``chromatin_regulatory._footprinting`` always pointed to:
TOBIAS ATACorrect (Tn5 bias correction) -> ScoreBigwig (footprint scores) -> BINDetect
(differential TF binding between two cell-type groups or conditions).

Pipeline:
  1. Build one BAM per explicit biological replicate/donor (scATAC barcode design or
     bulk replicate map), merging only technical BAMs within a replicate.
  2. ``TOBIAS ATACorrect`` + ``ScoreBigwig`` independently per replicate.
  3. ``TOBIAS BINDetect`` over all replicate signals to obtain motif mean scores.
  4. Test replicate mean scores with Welch tests, apply BH across motifs, and report
     condition-label permutations as null diagnostics. If replicate identity/counts are
     insufficient, retain only the legacy descriptive condition/group ranking.

No fabrication (ADR-002 / W2.2): TOBIAS, the genome FASTA, and the motif collection are
all required; any missing one is an honest ``ran: false`` with a concrete reason, never
an uncorrected footprint or an invented score. Site-level BINDetect p-values never enter
replicate inference. Differential TF binding is ASSOCIATIVE, not causal regulation.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Iterator

from aria.scripts._base import run_script

# Tn5 insertion offsets (the canonical +4 / -5 shift); applied by TOBIAS via
# --read_shift, kept here as the documented constant for the cut-site reads.
TN5_SHIFT = (4, -5)


def _skip(reason: str, **extra: Any) -> dict[str, Any]:
    """Honest not-run result (mirrors chromatin_regulatory._skip)."""
    return {"ran": False, "reason": reason, **extra}


# --------------------------------------------------------------------------- #
# Pure helpers (no pysam / no TOBIAS) — unit-testable in aria-env              #
# --------------------------------------------------------------------------- #

def _chrom_sizes_from_fai(fai_path: str) -> dict[str, int]:
    """Parse a samtools ``.fai`` index into ``{chrom: length}``. The BAM header and
    the fragment chrom filter both come from this, so a fragment on a contig absent
    from the genome is dropped rather than crashing the writer."""
    sizes: dict[str, int] = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    sizes[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    return sizes


def _load_barcode_groups(tsv_path: str) -> dict[str, str]:
    """Load a ``barcode<TAB>group`` TSV (header optional) into ``{barcode: group}``.
    Groups are the cell-type labels (e.g. transferred from the paired RNA)."""
    groups: dict[str, str] = {}
    with open(tsv_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            bc, grp = parts[0].strip(), parts[1].strip()
            if not bc or not grp or bc.lower() in ("barcode", "cell", "cell_id"):
                continue
            groups[bc] = grp
    return groups


def _load_barcode_design(tsv_path: str) -> dict[str, dict[str, str | None]]:
    """Load an explicit barcode/group/replicate design.

    The first two columns remain compatible with the historical
    ``barcode<TAB>group`` contract. A third ``replicate``/``donor`` column is
    required before scATAC footprint scores can be used inferentially; missing
    identity is represented as ``None`` and never guessed from a barcode or file
    name.
    """
    design: dict[str, dict[str, str | None]] = {}
    with open(tsv_path, encoding="utf-8") as fh:
        rows = csv.reader(fh, delimiter="\t")
        try:
            first = next(rows)
        except StopIteration:
            return design
        lowered = [cell.strip().lower() for cell in first]
        has_header = bool(lowered and lowered[0] in {"barcode", "cell", "cell_id"})
        if has_header:
            barcode_idx = 0
            group_idx = next((i for i, value in enumerate(lowered)
                              if value in {"group", "condition", "cell_type"}), 1)
            replicate_idx = next((i for i, value in enumerate(lowered)
                                  if value in {"replicate", "replicate_id", "donor",
                                               "donor_id", "sample", "sample_id"}), None)
            data_rows = rows
        else:
            barcode_idx, group_idx = 0, 1
            replicate_idx = 2 if len(first) >= 3 else None
            data_rows = itertools.chain([first], rows)
        for parts in data_rows:
            if len(parts) <= max(barcode_idx, group_idx):
                continue
            barcode = parts[barcode_idx].strip()
            group = parts[group_idx].strip()
            replicate = None
            if replicate_idx is not None and len(parts) > replicate_idx:
                replicate = parts[replicate_idx].strip() or None
            if barcode and group:
                design[barcode] = {"group": group, "replicate": replicate}
    return design


def _bh_adjust(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values, preserving NaN positions."""
    adjusted = [math.nan] * len(pvalues)
    valid = [(idx, float(value)) for idx, value in enumerate(pvalues)
             if math.isfinite(float(value))]
    valid.sort(key=lambda item: item[1])
    m = len(valid)
    running = 1.0
    for rank_from_end in range(m - 1, -1, -1):
        idx, pvalue = valid[rank_from_end]
        rank = rank_from_end + 1
        running = min(running, max(0.0, min(1.0, pvalue)) * m / rank)
        adjusted[idx] = min(1.0, running)
    return adjusted


def parse_bindetect_replicate_scores(
    results_txt: str,
    replicate_groups: dict[str, str],
) -> dict[str, dict[str, Any]]:
    """Read TOBIAS per-signal mean scores for explicit biological replicates.

    BINDetect site-level p-values are intentionally ignored. Only
    ``<signal>_mean_score`` columns enter the replicate model.
    """
    table: dict[str, dict[str, Any]] = {}
    with open(results_txt, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = set(reader.fieldnames or [])
        score_columns = {
            replicate: f"{replicate}_mean_score" for replicate in replicate_groups
        }
        if not score_columns or any(column not in fields
                                    for column in score_columns.values()):
            return table
        for row_number, row in enumerate(reader, start=1):
            tf_name = str(row.get("name") or "").strip()
            if not tf_name:
                continue
            motif_id = str(row.get("motif_id") or row.get("id") or tf_name).strip()
            feature_id = motif_id
            if feature_id in table:
                feature_id = f"{motif_id}__row_{row_number}"
            try:
                scores = {replicate: float(row[column])
                          for replicate, column in score_columns.items()}
                n_sites = int(float(row.get("total_tfbs") or 0))
            except (TypeError, ValueError):
                continue
            if all(math.isfinite(value) for value in scores.values()):
                table[feature_id] = {
                    "tf": tf_name,
                    "motif_id": motif_id,
                    "scores": scores,
                    "n_sites": n_sites,
                }
    return table


def _welch_pvalues(matrix: Any, idx_a: list[int], idx_b: list[int]) -> list[float]:
    """Vectorized two-sided Welch tests over rows of a score matrix."""
    import numpy as np
    from scipy.stats import ttest_ind

    with np.errstate(all="ignore"):
        result = ttest_ind(
            matrix[:, idx_a], matrix[:, idx_b], axis=1, equal_var=False,
            nan_policy="omit",
        )
    values = np.asarray(result.pvalue, dtype=float)
    # Identical constant groups yield NaN (0/0) and represent no difference.
    return [float(value) if math.isfinite(float(value)) else 1.0 for value in values]


def _paired_pvalues(matrix: Any, idx_a: list[int], idx_b: list[int]) -> list[float]:
    """Vectorized paired t-tests over rows of a score matrix."""
    import numpy as np
    from scipy.stats import ttest_rel

    with np.errstate(all="ignore"):
        result = ttest_rel(
            matrix[:, idx_a], matrix[:, idx_b], axis=1, nan_policy="omit")
    values = np.asarray(result.pvalue, dtype=float)
    return [float(value) if math.isfinite(float(value)) else 1.0 for value in values]


def _null_label_controls(
    matrix: Any,
    n_a: int,
    alpha: float,
    max_permutations: int,
    random_seed: int,
    paired: bool = False,
) -> dict[str, Any]:
    """Stress-test BH discoveries under label permutations.

    The observed assignment and its exact label complement are excluded so a
    real group effect is not mislabeled as a null control. These are diagnostics,
    not post-hoc thresholds and not a replacement for the primary BH correction.
    """
    n_total = int(matrix.shape[1])
    if paired:
        n_pairs = n_a
        masks = list(range(1, (2 ** n_pairs) - 1))
        if len(masks) > max_permutations:
            masks = random.Random(random_seed).sample(masks, max_permutations)
        fractions: list[float] = []
        for mask in masks:
            idx_a = [n_pairs + idx if mask & (1 << idx) else idx
                     for idx in range(n_pairs)]
            idx_b = [idx if mask & (1 << idx) else n_pairs + idx
                     for idx in range(n_pairs)]
            adjusted = _bh_adjust(_paired_pvalues(matrix, idx_a, idx_b))
            n_discoveries = sum(1 for value in adjusted
                                if math.isfinite(value) and value <= alpha)
            fractions.append(n_discoveries / int(matrix.shape[0]))
        return {
            "method": "within_donor_condition_label_swap_with_bh",
            "n_permutations": len(fractions),
            "observed_labels_excluded": True,
            "mean_discovery_fraction": (sum(fractions) / len(fractions)
                                        if fractions else 0.0),
            "max_discovery_fraction": max(fractions, default=0.0),
            "fraction_with_any_discoveries": (
                sum(value > 0 for value in fractions) / len(fractions)
                if fractions else 0.0
            ),
            "role": "diagnostic_only_not_a_post_hoc_significance_gate",
        }
    observed = tuple(range(n_a))
    complement = tuple(range(n_a, n_total)) if n_a == n_total - n_a else None
    total_assignments = math.comb(n_total, n_a)
    candidates: list[tuple[int, ...]] = []
    if total_assignments <= max_permutations + 2:
        source = itertools.combinations(range(n_total), n_a)
        candidates = [combo for combo in source
                      if combo != observed and combo != complement]
    else:
        rng = random.Random(random_seed)
        seen = {observed}
        if complement is not None:
            seen.add(complement)
        attempts = 0
        while len(candidates) < max_permutations and attempts < max_permutations * 50:
            combo = tuple(sorted(rng.sample(range(n_total), n_a)))
            attempts += 1
            if combo not in seen:
                seen.add(combo)
                candidates.append(combo)

    if len(candidates) > max_permutations:
        candidates = random.Random(random_seed).sample(candidates, max_permutations)
    fractions: list[float] = []
    all_idx = set(range(n_total))
    for idx_a_tuple in candidates:
        idx_a = list(idx_a_tuple)
        idx_b = sorted(all_idx.difference(idx_a))
        adjusted = _bh_adjust(_welch_pvalues(matrix, idx_a, idx_b))
        n_discoveries = sum(1 for value in adjusted
                            if math.isfinite(value) and value <= alpha)
        fractions.append(n_discoveries / int(matrix.shape[0]))
    return {
        "method": "condition_label_permutation_with_bh",
        "n_permutations": len(fractions),
        "observed_labels_excluded": True,
        "mean_discovery_fraction": (sum(fractions) / len(fractions)
                                    if fractions else 0.0),
        "max_discovery_fraction": max(fractions, default=0.0),
        "fraction_with_any_discoveries": (
            sum(value > 0 for value in fractions) / len(fractions)
            if fractions else 0.0
        ),
        "role": "diagnostic_only_not_a_post_hoc_significance_gate",
    }


def infer_replicate_footprints(
    score_table: dict[str, dict[str, Any]],
    replicate_groups: dict[str, str],
    group_a: str,
    group_b: str,
    *,
    alpha: float = 0.05,
    min_replicates_per_condition: int = 3,
    top_n: int = 15,
    max_label_permutations: int = 100,
    random_seed: int = 0,
    replicate_ids: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Test TF footprint mean scores with biological replicates as observations."""
    import numpy as np

    reps_a = [rep for rep, group in replicate_groups.items() if group == group_a]
    reps_b = [rep for rep, group in replicate_groups.items() if group == group_b]
    pairing_status = "independent"
    excluded_unpaired_replicates: list[str] = []
    if replicate_ids:
        ids_a = {replicate_ids.get(rep, rep): rep for rep in reps_a}
        ids_b = {replicate_ids.get(rep, rep): rep for rep in reps_b}
        common = set(ids_a).intersection(ids_b)
        if common:
            ordered_ids = sorted(common)
            excluded_unpaired_replicates = sorted(
                (set(ids_a).symmetric_difference(set(ids_b))))
            reps_a = [ids_a[replicate_id] for replicate_id in ordered_ids]
            reps_b = [ids_b[replicate_id] for replicate_id in ordered_ids]
            pairing_status = ("complete" if not excluded_unpaired_replicates
                              else "partial_complete_pair_subset")
    replicates = reps_a + reps_b
    complete = []
    for motif, record in sorted(score_table.items()):
        scores = record.get("scores") or {}
        try:
            values = [float(scores[rep]) for rep in replicates]
        except (KeyError, TypeError, ValueError):
            continue
        if all(math.isfinite(value) for value in values):
            complete.append((motif, record, values))

    base_inference = {
        "inferential_unit": "biological_replicate_or_donor",
        "replicates_per_condition": {group_a: len(reps_a), group_b: len(reps_b)},
        "minimum_replicates_per_condition": int(min_replicates_per_condition),
        "pairing_status": pairing_status,
        "excluded_unpaired_replicates": excluded_unpaired_replicates,
    }
    if (len(reps_a) < min_replicates_per_condition
            or len(reps_b) < min_replicates_per_condition):
        ranked = []
        for motif, record, values in complete:
            mean_a = sum(values[:len(reps_a)]) / len(reps_a) if reps_a else math.nan
            mean_b = sum(values[len(reps_a):]) / len(reps_b) if reps_b else math.nan
            ranked.append({
                "tf": str(record.get("tf") or motif),
                "motif_id": str(record.get("motif_id") or motif),
                "change": round(mean_a - mean_b, 4),
                "mean_score_a": mean_a, "mean_score_b": mean_b,
                "n_sites": int(record.get("n_sites") or 0),
            })
        return {
            "parsed": bool(complete),
            "n_motifs_tested": len(complete),
            "n_ranked_candidates": len(complete),
            "ranking_basis": {
                "fdr_controlled": False,
                "replicate_inference": False,
                "pseudobulk": False,
                "filter": "none; ranked by absolute replicate-mean difference",
                "note": "descriptive only because biological replication is insufficient",
            },
            "inference": {
                **base_inference,
                "status": "descriptive_only",
                "reason": "insufficient_biological_replicates",
            },
            f"top_toward_{group_a}": sorted(
                ranked, key=lambda row: -row["change"])[:top_n],
            f"top_toward_{group_b}": sorted(
                ranked, key=lambda row: row["change"])[:top_n],
        }
    if not complete:
        return {
            "parsed": False,
            "reason": "no motifs with complete finite replicate mean scores",
            "inference": {**base_inference, "status": "not_run"},
        }

    matrix = np.asarray([values for _, _, values in complete], dtype=float)
    idx_a = list(range(len(reps_a)))
    idx_b = list(range(len(reps_a), len(replicates)))
    paired = pairing_status != "independent"
    pvalues = (_paired_pvalues(matrix, idx_a, idx_b)
               if paired
               else _welch_pvalues(matrix, idx_a, idx_b))
    adjusted = _bh_adjust(pvalues)
    results = []
    for row_idx, ((motif, record, values), pvalue, padj) in enumerate(
            zip(complete, pvalues, adjusted)):
        mean_a = float(np.mean(matrix[row_idx, :len(reps_a)]))
        mean_b = float(np.mean(matrix[row_idx, len(reps_a):]))
        results.append({
            "tf": str(record.get("tf") or motif),
            "motif_id": str(record.get("motif_id") or motif),
            "change": round(mean_a - mean_b, 4),
            "mean_score_a": mean_a,
            "mean_score_b": mean_b,
            "pvalue": pvalue,
            "padj": padj,
            "n_sites": int(record.get("n_sites") or 0),
        })
    significant = [row for row in results if row["padj"] <= alpha]
    null_controls = _null_label_controls(
        matrix, len(reps_a), alpha, max_label_permutations, random_seed,
        paired=paired,
    )
    return {
        "parsed": True,
        "n_motifs_tested": len(results),
        "n_significant": len(significant),
        "alpha": alpha,
        "inference": {
            **base_inference,
            "status": "success",
            "test": ("paired_t_test_on_within_donor_mean_score_differences"
                     if paired
                     else "welch_t_test_on_replicate_mean_scores"),
            "multiple_testing": "benjamini_hochberg_across_motifs",
            "low_power_warning": min(len(reps_a), len(reps_b)) < 3,
            "null_label_controls": null_controls,
        },
        f"top_toward_{group_a}": sorted(
            significant, key=lambda row: -row["change"])[:top_n],
        f"top_toward_{group_b}": sorted(
            significant, key=lambda row: row["change"])[:top_n],
        "all_results": results,
    }


def _fragment_cut_reads(chrom: str, start: int, end: int, name: str,
                        read_len: int = 50) -> Iterator[dict[str, Any]]:
    """Yield the two Tn5 cut-site reads for one fragment as pure read specs (no pysam):
    a forward read at the left cut site and a reverse read at the right cut site. The
    pysam writer materializes these; keeping it pure makes the coordinate logic
    unit-testable. Degenerate fragments (end <= start) yield nothing."""
    if end <= start:
        return
    rlen = max(1, min(read_len, end - start))
    # Forward read anchored at the left Tn5 insertion.
    yield {"name": name, "chrom": chrom, "pos": start, "is_reverse": False,
           "length": rlen}
    # Reverse read anchored at the right Tn5 insertion.
    yield {"name": name, "chrom": chrom, "pos": max(start, end - rlen),
           "is_reverse": True, "length": rlen}


def _iter_fragments(frag_file: str, scan_limit: int | None = None
                    ) -> Iterator[tuple[str, int, int, str]]:
    """Stream ``(chrom, start, end, barcode)`` from a (gzipped) fragments TSV."""
    import gzip
    opener = gzip.open if str(frag_file).endswith(".gz") else open
    with opener(frag_file, "rt") as fh:
        for i, line in enumerate(fh):
            if scan_limit is not None and i >= scan_limit:
                break
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                try:
                    yield parts[0], int(parts[1]), int(parts[2]), parts[3]
                except ValueError:
                    continue


# --------------------------------------------------------------------------- #
# pysam BAM writer (gated)                                                     #
# --------------------------------------------------------------------------- #

def _write_group_bams(frag_file: str, barcode_groups: dict[str, str],
                      target_groups: Iterable[str], chrom_sizes: dict[str, int],
                      out_dir: Path, read_len: int = 50,
                      reuse: bool = True) -> dict[str, dict[str, Any]]:
    """Write one sorted+indexed pseudobulk BAM per target group from the fragments.
    Returns ``{group: {"bam": path, "n_fragments": int}}``. Requires pysam. When
    ``reuse`` and a group's sorted BAM + index already exist, that group is reused
    (counts recovered from the index) so a re-run does not rewrite multi-GB BAMs."""
    import pysam

    targets = {g for g in target_groups}
    reused: dict[str, dict[str, Any]] = {}
    if reuse:
        for g in list(targets):
            bam = out_dir / f"{g}.bam"
            if bam.is_file() and (out_dir / f"{g}.bam.bai").is_file():
                mapped = pysam.AlignmentFile(str(bam)).mapped
                reused[g] = {"bam": str(bam), "n_fragments": int(mapped // 2),
                             "reused": True}
                targets.discard(g)
    if not targets:
        return reused
    refs = list(chrom_sizes.keys())
    lens = [chrom_sizes[r] for r in refs]
    tid = {r: i for i, r in enumerate(refs)}
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": r, "LN": chrom_sizes[r]} for r in refs]}

    unsorted = {g: out_dir / f"{g}.unsorted.bam" for g in targets}
    writers = {g: pysam.AlignmentFile(str(p), "wb", header=header)
               for g, p in unsorted.items()}
    counts = {g: 0 for g in targets}
    try:
        for chrom, start, end, bc in _iter_fragments(frag_file):
            grp = barcode_groups.get(bc)
            if grp not in targets or chrom not in tid:
                continue
            wrote = False
            for r in _fragment_cut_reads(chrom, start, end, f"{bc}:{counts[grp]}",
                                         read_len):
                seg = pysam.AlignedSegment()
                seg.query_name = r["name"]
                seg.reference_id = tid[r["chrom"]]
                seg.reference_start = max(0, r["pos"])
                seg.flag = 16 if r["is_reverse"] else 0
                seg.mapping_quality = 60
                seg.cigarstring = f"{r['length']}M"
                seg.query_sequence = "N" * r["length"]
                seg.query_qualities = pysam.qualitystring_to_array("I" * r["length"])
                writers[grp].write(seg)
                wrote = True
            if wrote:
                counts[grp] += 1
    finally:
        for w in writers.values():
            w.close()

    out: dict[str, dict[str, Any]] = dict(reused)
    for g in targets:
        sorted_bam = out_dir / f"{g}.bam"
        pysam.sort("-o", str(sorted_bam), str(unsorted[g]))
        pysam.index(str(sorted_bam))
        unsorted[g].unlink(missing_ok=True)
        out[g] = {"bam": str(sorted_bam), "n_fragments": counts[g]}
    return out


# --------------------------------------------------------------------------- #
# TOBIAS orchestration (gated; CLI subprocess)                                 #
# --------------------------------------------------------------------------- #

def summarize_bindetect(results_txt: str, group_a: str, group_b: str,
                        top_n: int = 15, pvalue_max: float = 0.05,
                        min_sites: int = 20) -> dict[str, Any]:
    """Pure summary of a TOBIAS ``bindetect_results.txt``: the top differential-binding
    TFs toward each group (by the ``<A>_<B>_change`` score), plus counts.

    INTERIM (preprint audit B7, 2026-07-09): this is a DESCRIPTIVE candidate ranking,
    NOT FDR-controlled significance. TOBIAS pools every cell/BAM of a condition into a
    single pseudobulk, so the per-site ``<A>_<B>_pvalue`` ignores biological replication
    and cannot support a significance claim; the ``pvalue_max``/``min_sites`` gate only
    prunes the ranked candidate list by an uncorrected p and a site-count floor. Real
    inference (per-replicate/donor model + BH multiplicity + null-label controls) is
    deferred to the full B7 fix. Until then the count is exported as
    ``n_ranked_candidates`` with an explicit ``ranking_basis`` disclosure — never as
    ``n_significant``. Differential TF binding is ASSOCIATIVE, never causal. Returns an
    honest ``parsed: false`` when the expected columns are absent."""
    import csv

    change_col = f"{group_a}_{group_b}_change"
    pval_col = f"{group_a}_{group_b}_pvalue"
    rows: list[tuple[str, float, float, int]] = []
    with open(results_txt) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or change_col not in reader.fieldnames:
            return {"parsed": False, "reason": f"missing column {change_col}"}
        for r in reader:
            try:
                rows.append((r["name"], float(r[change_col]),
                             float(r.get(pval_col, "nan")), int(float(r["total_tfbs"]))))
            except (ValueError, KeyError):
                continue

    def _ranked(rs):
        # Descriptive candidate filter (see docstring): an uncorrected p<pvalue_max
        # gate + site-count floor over a pseudobulk-per-condition contrast. NOT
        # FDR-controlled significance; it only prunes the ranked candidate list.
        return [x for x in rs if (x[2] == x[2]) and x[2] < pvalue_max and x[3] >= min_sites]

    ranked = _ranked(rows)
    fmt = lambda x: {"tf": x[0], "change": round(x[1], 4), "pvalue": x[2], "n_sites": x[3]}
    return {
        "parsed": True,
        "n_motifs_tested": len(rows),
        "n_ranked_candidates": len(ranked),
        "ranking_basis": {
            "fdr_controlled": False,
            "replicate_inference": False,
            "pseudobulk": True,
            "filter": f"uncorrected p<{pvalue_max} and n_sites>={min_sites}",
            "note": ("descriptive candidate ranking by |change|; NOT FDR-controlled "
                     "significance. TOBIAS pools each condition into one pseudobulk, so "
                     "per-site p-values ignore biological replication."),
        },
        f"top_toward_{group_a}": [fmt(x) for x in sorted(ranked, key=lambda x: -x[1])[:top_n]],
        f"top_toward_{group_b}": [fmt(x) for x in sorted(ranked, key=lambda x: x[1])[:top_n]],
    }


def _top_tfs(summary: dict, group_a: str, group_b: str, n: int = 5) -> list[str]:
    """The TFs to draw aggregate footprints for: the top differential toward EACH group
    (both directions, so no biology is dropped). Deduplicated, order preserved."""
    tfs: list[str] = []
    for grp in (group_a, group_b):
        for x in (summary.get(f"top_toward_{grp}") or [])[:n]:
            name = str(x.get("tf", "")).strip()
            if name and name not in tfs:
                tfs.append(name)
    return tfs


def _aggregate_plots(bindetect_dir: Path, signals: dict[str, str], tfs: list[str],
                     out_dir: Path, cores: int = 8) -> dict[str, dict[str, str]]:
    """Per-TF aggregate footprint plots via TOBIAS PlotAggregate (the canonical tool —
    ARIA orchestrates it, does not re-implement footprint aggregation). For each TF the
    corrected signal of BOTH groups is averaged around the motif's TFBS; dual PNG+SVG
    (publication, P4.1 discipline). A TF whose TFBS bed is absent is skipped honestly."""
    agg_dir = out_dir / "aggregate"
    agg_dir.mkdir(parents=True, exist_ok=True)
    log = out_dir / "tobias.log"
    sig_args: list[str] = []
    for g in signals:
        sig_args += [str(signals[g])]
    out: dict[str, dict[str, str]] = {}
    for tf in tfs:
        sanitized = tf.replace("::", "")
        beds = sorted(bindetect_dir.glob(f"{sanitized}_*/beds/*_all.bed"))
        if not beds:
            continue
        made: dict[str, str] = {}
        for ext in ("png", "svg"):
            dest = agg_dir / f"{sanitized}_aggregate.{ext}"
            try:
                _run(["TOBIAS", "PlotAggregate", "--TFBS", str(beds[0]),
                      "--signals", *sig_args, "--output", str(dest),
                      "--title", f"{tf} footprint ({' vs '.join(signals)})",
                      "--share-y", "both", "--plot-boundaries"], log)
                made[ext] = str(dest)
            except Exception:
                continue
        if made.get("png"):
            out[tf] = made
    return out


def _run(cmd: list[str], log: Path) -> None:
    with open(log, "a") as fh:
        fh.write("\n$ " + " ".join(cmd) + "\n")
        subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)


def _tobias_pipeline(group_bams: dict[str, dict[str, Any]], genome_fasta: str,
                     peaks_bed: str, motif_meme: str, group_a: str, group_b: str,
                     out_dir: Path, cores: int = 8) -> dict[str, Any]:
    """ATACorrect (per group) -> ScoreBigwig (per group) -> BINDetect (A vs B).
    ``cores`` is passed to every TOBIAS step (single-threaded is pathologically slow
    genome-wide)."""
    log = out_dir / "tobias.log"
    c = ["--cores", str(cores)]
    signals = {}
    for g in (group_a, group_b):
        bam = group_bams[g]["bam"]
        _run(["TOBIAS", "ATACorrect", "--bam", bam, "--genome", genome_fasta,
              "--peaks", peaks_bed, "--outdir", str(out_dir), "--prefix", g] + c, log)
        corrected = out_dir / f"{g}_corrected.bw"
        score = out_dir / f"{g}_footprints.bw"
        _run(["TOBIAS", "ScoreBigwig", "--signal", str(corrected),
              "--regions", peaks_bed, "--output", str(score)] + c, log)
        signals[g] = str(score)

    bindetect_out = out_dir / "bindetect"
    _run(["TOBIAS", "BINDetect", "--motifs", motif_meme,
          "--signals", signals[group_a], signals[group_b],
          "--genome", genome_fasta, "--peaks", peaks_bed,
          "--cond_names", group_a, group_b, "--outdir", str(bindetect_out)] + c, log)

    results = bindetect_out / "bindetect_results.txt"
    return {"bindetect_results": str(results) if results.is_file() else None,
            "footprint_signals": signals, "bindetect_outdir": str(bindetect_out)}


def _safe_signal_label(group: str, replicate: str, ordinal: int) -> str:
    """Stable TOBIAS column label without treating names as biological metadata."""
    group_part = re.sub(r"[^A-Za-z0-9]+", "_", str(group)).strip("_") or "group"
    rep_part = re.sub(r"[^A-Za-z0-9]+", "_", str(replicate)).strip("_") or "replicate"
    return f"{group_part}_rep{ordinal}_{rep_part}"


def _tobias_replicate_pipeline(
    replicate_bams: dict[str, dict[str, Any]],
    replicate_groups: dict[str, str],
    genome_fasta: str,
    peaks_bed: str,
    motif_meme: str,
    out_dir: Path,
    cores: int = 8,
) -> dict[str, Any]:
    """Run TOBIAS per biological replicate, then one multi-signal BINDetect."""
    log = out_dir / "tobias.log"
    c = ["--cores", str(cores)]
    footprint_signals: dict[str, str] = {}
    corrected_signals: dict[str, str] = {}
    for signal_name in replicate_groups:
        bam = replicate_bams[signal_name]["bam"]
        _run(["TOBIAS", "ATACorrect", "--bam", bam, "--genome", genome_fasta,
              "--peaks", peaks_bed, "--outdir", str(out_dir),
              "--prefix", signal_name] + c, log)
        corrected = out_dir / f"{signal_name}_corrected.bw"
        score = out_dir / f"{signal_name}_footprints.bw"
        _run(["TOBIAS", "ScoreBigwig", "--signal", str(corrected),
              "--regions", peaks_bed, "--output", str(score)] + c, log)
        corrected_signals[signal_name] = str(corrected)
        footprint_signals[signal_name] = str(score)

    bindetect_out = out_dir / "bindetect_replicates"
    ordered = list(replicate_groups)
    _run(["TOBIAS", "BINDetect", "--motifs", motif_meme,
          "--signals", *[footprint_signals[name] for name in ordered],
          "--genome", genome_fasta, "--peaks", peaks_bed,
          "--cond_names", *ordered, "--outdir", str(bindetect_out)] + c, log)
    results = bindetect_out / "bindetect_results.txt"
    return {
        "bindetect_results": str(results) if results.is_file() else None,
        "footprint_signals": footprint_signals,
        "corrected_signals": corrected_signals,
        "bindetect_outdir": str(bindetect_out),
        "replicate_groups": replicate_groups,
    }


def _write_replicate_inference_table(summary: dict[str, Any], out_dir: Path) -> str | None:
    """Persist the complete BH table and keep the JSON summary compact."""
    rows = summary.pop("all_results", None)
    if not rows:
        return None
    path = out_dir / "replicate_footprint_inference.tsv"
    fields = ["motif_id", "tf", "mean_score_a", "mean_score_b", "change",
              "pvalue", "padj", "n_sites"]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({key: row.get(key) for key in fields} for row in rows)
    summary["results_table"] = str(path)
    return str(path)


def _replicate_summary(
    tobias: dict[str, Any],
    group_a: str,
    group_b: str,
    out_dir: Path,
    params: dict[str, Any],
) -> dict[str, Any]:
    results_path = tobias.get("bindetect_results")
    if not results_path:
        return {"parsed": False, "reason": "BINDetect produced no results table"}
    score_table = parse_bindetect_replicate_scores(
        str(results_path), tobias["replicate_groups"])
    summary = infer_replicate_footprints(
        score_table, tobias["replicate_groups"], group_a, group_b,
        alpha=float(params.get("footprint_fdr", 0.05)),
        min_replicates_per_condition=int(
            params.get("min_replicates_per_condition", 3)),
        max_label_permutations=int(params.get("max_label_permutations", 100)),
        random_seed=int(params.get("permutation_seed", 0)),
        replicate_ids=tobias.get("replicate_ids"),
    )
    _write_replicate_inference_table(summary, out_dir)
    return summary


# --------------------------------------------------------------------------- #
# Entry point                                                                 #
# --------------------------------------------------------------------------- #

def chromatin_footprint_tobias(params: dict) -> dict[str, Any]:
    """Run the TOBIAS footprinting + differential-binding pipeline, or skip honestly.

    Required params: ``fragments_file``, ``genome_fasta``, ``peaks_bed``,
    ``motif_meme``, ``barcode_groups`` (TSV), ``genome_fai`` (or sibling ``.fai``),
    ``group_a``, ``group_b``, ``output_dir``.
    """
    frag = params.get("fragments_file")
    genome = params.get("genome_fasta")
    peaks = params.get("peaks_bed")
    motifs = params.get("motif_meme")
    groups_tsv = params.get("barcode_groups")
    group_a = params.get("group_a")
    group_b = params.get("group_b")
    out_dir = Path(params.get("output_dir", "."))

    if shutil.which("TOBIAS") is None:
        return _skip("TOBIAS not installed (needs aria-tobias-env); footprinting "
                     "refuses uncorrected output", method="tobias")
    for label, val in (("fragments_file", frag), ("genome_fasta", genome),
                       ("peaks_bed", peaks), ("motif_meme", motifs),
                       ("barcode_groups", groups_tsv)):
        if not val or not Path(str(val)).is_file():
            return _skip(f"{label} missing or not a file: {val}", method="tobias")
    if not group_a or not group_b:
        return _skip("group_a and group_b (the two cell-type groups to contrast) "
                     "are required", method="tobias")

    fai = params.get("genome_fai") or f"{genome}.fai"
    if not Path(fai).is_file():
        return _skip(f"genome .fai index not found: {fai} (samtools faidx)",
                     method="tobias")

    out_dir.mkdir(parents=True, exist_ok=True)
    chrom_sizes = _chrom_sizes_from_fai(fai)
    barcode_design = _load_barcode_design(str(groups_tsv))
    barcode_groups = {barcode: str(row["group"])
                      for barcode, row in barcode_design.items()}
    present = set(barcode_groups.values())
    missing = [g for g in (group_a, group_b) if g not in present]
    if missing:
        return _skip(f"groups absent from barcode_groups: {missing}", method="tobias")

    min_replicates = int(params.get("min_replicates_per_condition", 3))
    explicit_replicates = bool(barcode_design) and all(
        row.get("replicate") for row in barcode_design.values()
        if row.get("group") in {group_a, group_b}
    )
    reps_by_group = {
        group: sorted({str(row["replicate"]) for row in barcode_design.values()
                       if row.get("group") == group and row.get("replicate")})
        for group in (group_a, group_b)
    }
    replicate_ready = explicit_replicates and all(
        len(reps_by_group[group]) >= min_replicates for group in (group_a, group_b)
    )

    if replicate_ready:
        signal_for_pair: dict[tuple[str, str], str] = {}
        replicate_groups: dict[str, str] = {}
        replicate_ids: dict[str, str] = {}
        for group in (group_a, group_b):
            for ordinal, replicate in enumerate(reps_by_group[group], start=1):
                signal = _safe_signal_label(str(group), replicate, ordinal)
                signal_for_pair[(str(group), replicate)] = signal
                replicate_groups[signal] = str(group)
                replicate_ids[signal] = replicate
        barcode_signals = {
            barcode: signal_for_pair[(str(row["group"]), str(row["replicate"]))]
            for barcode, row in barcode_design.items()
            if (str(row.get("group")), str(row.get("replicate"))) in signal_for_pair
        }
        replicate_bams = _write_group_bams(
            str(frag), barcode_signals, replicate_groups, chrom_sizes, out_dir)
        empty = [name for name, info in replicate_bams.items()
                 if info["n_fragments"] == 0]
        if empty:
            return _skip(
                f"biological replicate(s) have no fragments after barcode mapping: {empty}",
                method="tobias", replicate_bams=replicate_bams)
        tobias = _tobias_replicate_pipeline(
            replicate_bams, replicate_groups, str(genome), str(peaks), str(motifs),
            out_dir, cores=int(params.get("cores", 8)))
        tobias["replicate_ids"] = replicate_ids
        summary = _replicate_summary(tobias, str(group_a), str(group_b), out_dir, params)
        aggregate = {}
        if (params.get("make_aggregate_plots", True) and summary.get("parsed")
                and tobias.get("bindetect_outdir")):
            aggregate = _aggregate_plots(
                Path(tobias["bindetect_outdir"]), tobias["corrected_signals"],
                _top_tfs(summary, str(group_a), str(group_b)), out_dir,
                cores=int(params.get("cores", 8)))
        return {
            "ran": True,
            "method": "tobias_replicate_atacorrect_scorebigwig_bindetect_welch_bh",
            "group_a": group_a, "group_b": group_b,
            "replicate_bams": replicate_bams,
            "differential_summary": summary,
            "aggregate_plots": aggregate,
            **tobias,
            "caveat": (
                "TF footprint inference uses biological replicate/donor mean scores "
                "with Welch tests and BH across motifs. Label permutations are null "
                "diagnostics. Footprint differences remain associative, not causal."
            ),
        }

    group_bams = _write_group_bams(str(frag), barcode_groups, (group_a, group_b),
                                   chrom_sizes, out_dir)
    for g in (group_a, group_b):
        if group_bams[g]["n_fragments"] == 0:
            return _skip(f"group '{g}' has no fragments after barcode mapping",
                         method="tobias", group_bams=group_bams)

    tobias = _tobias_pipeline(group_bams, str(genome), str(peaks), str(motifs),
                              group_a, group_b, out_dir,
                              cores=int(params.get("cores", 8)))
    summary = (summarize_bindetect(tobias["bindetect_results"], group_a, group_b)
               if tobias.get("bindetect_results") else {"parsed": False,
               "reason": "BINDetect produced no results table"})
    summary["inference"] = {
        "status": "descriptive_only",
        "reason": ("missing_explicit_replicate_or_donor_identity"
                   if not explicit_replicates else "insufficient_biological_replicates"),
        "inferential_unit": "biological_replicate_or_donor",
        "replicates_per_condition": {group: len(reps_by_group[group])
                                     for group in (group_a, group_b)},
        "minimum_replicates_per_condition": min_replicates,
    }
    aggregate = {}
    if (params.get("make_aggregate_plots", True) and summary.get("parsed")
            and tobias.get("bindetect_outdir")):
        aggregate = _aggregate_plots(
            Path(tobias["bindetect_outdir"]),
            {group_a: str(out_dir / f"{group_a}_corrected.bw"),
             group_b: str(out_dir / f"{group_b}_corrected.bw")},
            _top_tfs(summary, group_a, group_b), out_dir,
            cores=int(params.get("cores", 8)))
    return {
        "ran": True,
        "method": "tobias_atacorrect_scorebigwig_bindetect",
        "group_a": group_a, "group_b": group_b,
        "group_bams": group_bams,
        "differential_summary": summary,
        "aggregate_plots": aggregate,
        **tobias,
        "caveat": (
            "Differential TF binding (BINDetect) is an associative footprint-signal "
            "difference between cell-type groups, not causal regulation. Footprints "
            "are Tn5-bias-corrected (ATACorrect). Without explicit sufficient donor/"
            "replicate identity, cell groups remain descriptive pseudobulks."
        ),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Bulk ATAC footprinting (replicate inference or descriptive condition pools)
# ─────────────────────────────────────────────────────────────────────────────

def _merge_condition_bams(condition_bams: dict[str, list],
                          out_dir: Path, reuse: bool = True
                          ) -> dict[str, dict[str, Any]]:
    """Merge the per-sample BAMs of each condition into one sorted+indexed
    condition BAM (the bulk analogue of the scATAC per-group pseudobulk BAM).
    Returns ``{condition: {"bam": path|None, "n_fragments": int, "n_bams": int}}``.
    Requires pysam; honest (``bam=None``/``0``) when a condition has no readable BAM.
    Coordinate-sorted inputs are assumed (post-alignment); ``pysam.merge`` keeps the
    output sorted. ``reuse`` recovers an existing merged+indexed BAM without rewriting."""
    import pysam

    out: dict[str, dict[str, Any]] = {}
    for cond, bams in condition_bams.items():
        merged = out_dir / f"{cond}.bam"
        index = out_dir / f"{cond}.bam.bai"
        valid = [str(b) for b in (bams or []) if b and Path(str(b)).is_file()]
        if reuse and merged.is_file() and index.is_file():
            out[cond] = {"bam": str(merged),
                         "n_fragments": int(pysam.AlignmentFile(str(merged)).mapped),
                         "n_bams": len(valid)}
            continue
        if not valid:
            out[cond] = {"bam": None, "n_fragments": 0, "n_bams": 0}
            continue
        pysam.merge("-f", str(merged), *valid)
        pysam.index(str(merged))
        out[cond] = {"bam": str(merged),
                     "n_fragments": int(pysam.AlignmentFile(str(merged)).mapped),
                     "n_bams": len(valid)}
    return out


def _prepare_bulk_replicate_bams(
    replicate_bams: dict[str, dict[str, Any]],
    group_a: str,
    group_b: str,
    out_dir: Path,
) -> tuple[dict[str, dict[str, Any]], dict[str, str], dict[str, str]]:
    """Merge technical BAMs within each explicitly named biological replicate."""
    merge_inputs: dict[str, list[Any]] = {}
    replicate_groups: dict[str, str] = {}
    original_ids: dict[str, str] = {}
    for group in (group_a, group_b):
        group_replicates = replicate_bams.get(group) or {}
        if not isinstance(group_replicates, dict):
            continue
        for ordinal, (replicate_id, bams) in enumerate(group_replicates.items(), start=1):
            signal = _safe_signal_label(group, str(replicate_id), ordinal)
            paths = bams if isinstance(bams, list) else [bams]
            merge_inputs[signal] = paths
            replicate_groups[signal] = group
            original_ids[signal] = str(replicate_id)
    replicate_dir = out_dir / "replicate_bams"
    replicate_dir.mkdir(parents=True, exist_ok=True)
    prepared = _merge_condition_bams(merge_inputs, replicate_dir)
    return prepared, replicate_groups, original_ids


def chromatin_footprint_tobias_bulk(params: dict) -> dict[str, Any]:
    """Condition-level differential TF footprinting for bulk ATAC.

    An explicit ``replicate_bams`` map runs ATACorrect/ScoreBigwig independently per
    biological replicate, then models replicate mean scores with Welch+BH. The legacy
    ``condition_bams`` shape remains a descriptive compatibility path only.

    Required params: ``replicate_bams`` (``{condition: {replicate: [bam, ...]}}``) or
    descriptive ``condition_bams`` (``{condition: [bam, ...]}``), ``genome_fasta``,
    ``peaks_bed``, ``motif_meme``, ``group_a``, ``group_b`` (the two conditions),
    ``genome_fai`` (or sibling ``.fai``), ``output_dir``. Any missing one is an honest
    ``ran: false`` with a concrete reason — footprinting never reports uncorrected
    output, and a missing condition BAM is never fabricated."""
    condition_bams = params.get("condition_bams") or {}
    replicate_design = params.get("replicate_bams") or {}
    genome = params.get("genome_fasta")
    peaks = params.get("peaks_bed")
    motifs = params.get("motif_meme")
    group_a = params.get("group_a")
    group_b = params.get("group_b")
    out_dir = Path(params.get("output_dir", "."))

    if shutil.which("TOBIAS") is None:
        return _skip("TOBIAS not installed (needs aria-tobias-env); footprinting "
                     "refuses uncorrected output", method="tobias")
    for label, val in (("genome_fasta", genome), ("peaks_bed", peaks),
                       ("motif_meme", motifs)):
        if not val or not Path(str(val)).is_file():
            return _skip(f"{label} missing or not a file: {val}", method="tobias")
    if not group_a or not group_b:
        return _skip("group_a and group_b (the two conditions to contrast) are "
                     "required", method="tobias")
    if (group_a not in condition_bams or group_b not in condition_bams) and replicate_design:
        condition_bams = {
            group: [bam for bams in (replicate_design.get(group) or {}).values()
                    for bam in (bams if isinstance(bams, list) else [bams])]
            for group in (group_a, group_b)
        }
    if group_a not in condition_bams or group_b not in condition_bams:
        return _skip("condition_bams must contain BAM lists for both group_a and "
                     "group_b", method="tobias")
    fai = params.get("genome_fai") or f"{genome}.fai"
    if not Path(str(fai)).is_file():
        return _skip(f"genome .fai index not found: {fai} (samtools faidx)",
                     method="tobias")

    out_dir.mkdir(parents=True, exist_ok=True)
    min_replicates = int(params.get("min_replicates_per_condition", 3))
    replicate_counts = {
        group: (len(replicate_design.get(group) or {})
                if isinstance(replicate_design.get(group), dict) else 0)
        for group in (group_a, group_b)
    }
    replicate_ready = all(
        replicate_counts[group] >= min_replicates for group in (group_a, group_b)
    )
    if replicate_ready:
        replicate_bam_info, replicate_groups, original_ids = (
            _prepare_bulk_replicate_bams(
                replicate_design, str(group_a), str(group_b), out_dir))
        empty = [name for name, info in replicate_bam_info.items()
                 if not info.get("bam") or info.get("n_fragments", 0) == 0]
        if empty:
            return _skip(
                f"biological replicate(s) have no readable BAM reads: {empty}",
                method="tobias", replicate_bams=replicate_bam_info)
        tobias = _tobias_replicate_pipeline(
            replicate_bam_info, replicate_groups, str(genome), str(peaks),
            str(motifs), out_dir, cores=int(params.get("cores", 8)))
        tobias["replicate_ids"] = original_ids
        summary = _replicate_summary(tobias, str(group_a), str(group_b), out_dir, params)
        aggregate = {}
        if (params.get("make_aggregate_plots", True) and summary.get("parsed")
                and tobias.get("bindetect_outdir")):
            aggregate = _aggregate_plots(
                Path(tobias["bindetect_outdir"]), tobias["corrected_signals"],
                _top_tfs(summary, str(group_a), str(group_b)), out_dir,
                cores=int(params.get("cores", 8)))
        return {
            "ran": True,
            "method": "tobias_replicate_atacorrect_scorebigwig_bindetect_welch_bh",
            "data_type": "bulk_ATAC",
            "group_a": group_a, "group_b": group_b,
            "group_label": "Conditions", "group_kind": "conditions",
            "replicate_bams": replicate_bam_info,
            "replicate_ids": original_ids,
            "differential_summary": summary,
            "aggregate_plots": aggregate,
            **tobias,
            "caveat": (
                "TF footprint inference uses biological-replicate mean scores with "
                "Welch tests and BH across motifs; label permutations are null "
                "diagnostics. Footprint differences remain associative, not causal."
            ),
        }

    group_bams = _merge_condition_bams(
        {group_a: condition_bams[group_a], group_b: condition_bams[group_b]},
        out_dir)
    for g in (group_a, group_b):
        if not group_bams[g]["bam"] or group_bams[g]["n_fragments"] == 0:
            return _skip(f"condition '{g}' has no readable BAM reads after merge",
                         method="tobias", group_bams=group_bams)

    tobias = _tobias_pipeline(group_bams, str(genome), str(peaks), str(motifs),
                              group_a, group_b, out_dir,
                              cores=int(params.get("cores", 8)))
    summary = (summarize_bindetect(tobias["bindetect_results"], group_a, group_b)
               if tobias.get("bindetect_results") else {"parsed": False,
               "reason": "BINDetect produced no results table"})
    summary["inference"] = {
        "status": "descriptive_only",
        "reason": ("missing_explicit_biological_replicate_design"
                   if not replicate_design else "insufficient_biological_replicates"),
        "inferential_unit": "biological_replicate_or_donor",
        "replicates_per_condition": replicate_counts,
        "minimum_replicates_per_condition": min_replicates,
    }
    aggregate = {}
    if (params.get("make_aggregate_plots", True) and summary.get("parsed")
            and tobias.get("bindetect_outdir")):
        aggregate = _aggregate_plots(
            Path(tobias["bindetect_outdir"]),
            {group_a: str(out_dir / f"{group_a}_corrected.bw"),
             group_b: str(out_dir / f"{group_b}_corrected.bw")},
            _top_tfs(summary, group_a, group_b), out_dir,
            cores=int(params.get("cores", 8)))
    return {
        "ran": True,
        "method": "tobias_atacorrect_scorebigwig_bindetect",
        "data_type": "bulk_ATAC",
        "group_a": group_a, "group_b": group_b,
        "group_label": "Conditions", "group_kind": "conditions",
        "group_bams": group_bams,
        "differential_summary": summary,
        "aggregate_plots": aggregate,
        **tobias,
        "caveat": (
            "Differential TF binding (BINDetect) is an associative footprint-signal "
            "difference between conditions, not causal regulation. Footprints are "
            "Tn5-bias-corrected (ATACorrect). Without an explicit sufficient "
            "biological-replicate design, condition BAMs remain descriptive pools."
        ),
    }


def chromatin_footprint_tobias_dispatch(params: dict) -> dict[str, Any]:
    """EnvironmentManager IPC entrypoint.

    The CLI remains available for benchmark/repro runs; live agent dispatch uses
    this JSON-IPC wrapper so all inputs and honest skips are captured in the
    normal ARIA subprocess contract.
    """
    mode = params.get("mode", "scatac")
    if mode == "bulk":
        return chromatin_footprint_tobias_bulk(params)
    return chromatin_footprint_tobias(params)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("scatac", "bulk"), default="scatac",
                   help="scatac: split fragments by barcode group; "
                        "bulk: merge per-sample BAMs per condition")
    # scATAC inputs (fragments + barcode->group).
    p.add_argument("--fragments-file")
    p.add_argument("--barcode-groups",
                   help="barcode<TAB>group[<TAB>replicate/donor] TSV (scatac mode)")
    # bulk inputs (per-condition BAM lists).
    p.add_argument("--condition-bams",
                   help="JSON {condition: [bam, ...]} (bulk mode)")
    p.add_argument("--replicate-bams",
                   help="JSON {condition: {biological_replicate: [bam, ...]}}")
    # shared.
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--peaks-bed", required=True)
    p.add_argument("--motif-meme", required=True)
    p.add_argument("--group-a", required=True)
    p.add_argument("--group-b", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--output-json", default=None)
    p.add_argument("--cores", type=int, default=8)
    p.add_argument("--min-replicates-per-condition", type=int, default=3)
    p.add_argument("--footprint-fdr", type=float, default=0.05)
    p.add_argument("--max-label-permutations", type=int, default=100)
    p.add_argument("--permutation-seed", type=int, default=0)
    p.add_argument("--skip-aggregate-plots", action="store_true")
    args = p.parse_args(argv)

    shared = {
        "genome_fasta": args.genome_fasta, "peaks_bed": args.peaks_bed,
        "motif_meme": args.motif_meme, "group_a": args.group_a,
        "group_b": args.group_b, "output_dir": args.output_dir,
        "cores": args.cores,
        "min_replicates_per_condition": args.min_replicates_per_condition,
        "footprint_fdr": args.footprint_fdr,
        "max_label_permutations": args.max_label_permutations,
        "permutation_seed": args.permutation_seed,
        "make_aggregate_plots": not args.skip_aggregate_plots,
    }
    if args.mode == "bulk":
        if not args.condition_bams and not args.replicate_bams:
            p.error("--condition-bams or --replicate-bams is required in bulk mode")
        condition_bams = (json.loads(Path(args.condition_bams).read_text())
                          if args.condition_bams else {})
        replicate_bams = (json.loads(Path(args.replicate_bams).read_text())
                          if args.replicate_bams else {})
        res = chromatin_footprint_tobias_bulk(
            {**shared, "condition_bams": condition_bams,
             "replicate_bams": replicate_bams})
    else:
        if not args.fragments_file or not args.barcode_groups:
            p.error("--fragments-file and --barcode-groups are required in scatac mode")
        res = chromatin_footprint_tobias({
            **shared, "fragments_file": args.fragments_file,
            "barcode_groups": args.barcode_groups,
        })
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(res, indent=2, sort_keys=True))
    print(json.dumps({k: v for k, v in res.items()
                      if k not in {"group_bams", "replicate_bams"}}, indent=2))
    return 0


if __name__ == "__main__":
    if len(sys.argv) == 3:
        run_script(chromatin_footprint_tobias_dispatch)
    else:
        raise SystemExit(main())
