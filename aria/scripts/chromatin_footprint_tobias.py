#!/usr/bin/env python3
"""scATAC P4.3 — Tn5-bias-corrected footprinting + differential TF binding (TOBIAS).

Runs in the dedicated ``aria-tobias-env`` (``envs/aria-tobias-env.yml``). This is the
backend the honest stub ``chromatin_regulatory._footprinting`` always pointed to:
TOBIAS ATACorrect (Tn5 bias correction) -> ScoreBigwig (footprint scores) -> BINDetect
(differential TF binding between two cell-type groups).

Pipeline:
  1. Split the ATAC fragments into pseudobulk-per-cell-type BAMs (barcode -> group from
     a label TSV, e.g. transferred from the paired RNA). Fragments -> BAM via pysam +
     chrom sizes from the genome ``.fai``.
  2. ``TOBIAS ATACorrect`` per group (genome FASTA + peaks BED) -> bias-corrected signal.
  3. ``TOBIAS ScoreBigwig`` per group -> footprint-score bigwig.
  4. ``TOBIAS BINDetect`` (the two groups, JASPAR2024 MEME motifs, genome, peaks) ->
     differential TF binding table.

No fabrication (ADR-002 / W2.2): TOBIAS, the genome FASTA, and the motif collection are
all required; any missing one is an honest ``ran: false`` with a concrete reason, never
an uncorrected footprint or an invented score. Differential TF binding is reported as an
ASSOCIATIVE signal (a TF whose footprint differs between groups), not causal regulation.
"""
from __future__ import annotations

import argparse
import json
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
    barcode_groups = _load_barcode_groups(str(groups_tsv))
    present = set(barcode_groups.values())
    missing = [g for g in (group_a, group_b) if g not in present]
    if missing:
        return _skip(f"groups absent from barcode_groups: {missing}", method="tobias")

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
    aggregate = {}
    if summary.get("parsed") and tobias.get("bindetect_outdir"):
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
            "are Tn5-bias-corrected (ATACorrect); single-cell footprinting is noisy "
            "so groups are pseudobulk."
        ),
    }


# ─────────────────────────────────────────────────────────────────────────────
# B4: bulk ATAC condition-level footprinting (per-sample BAMs merged per condition)
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


def chromatin_footprint_tobias_bulk(params: dict) -> dict[str, Any]:
    """Condition-level differential TF footprinting for bulk ATAC.

    Mirrors the scATAC entry but builds the two pseudobulk BAMs by MERGING each
    condition's per-sample BAMs (not by splitting fragments on a barcode label), then
    reuses the SAME TOBIAS core (``_tobias_pipeline`` → ATACorrect/ScoreBigwig/
    BINDetect) + ``summarize_bindetect`` + ``_aggregate_plots``.

    Required params: ``condition_bams`` (``{condition: [bam, ...]}``), ``genome_fasta``,
    ``peaks_bed``, ``motif_meme``, ``group_a``, ``group_b`` (the two conditions),
    ``genome_fai`` (or sibling ``.fai``), ``output_dir``. Any missing one is an honest
    ``ran: false`` with a concrete reason — footprinting never reports uncorrected
    output, and a missing condition BAM is never fabricated."""
    condition_bams = params.get("condition_bams") or {}
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
    if group_a not in condition_bams or group_b not in condition_bams:
        return _skip("condition_bams must contain BAM lists for both group_a and "
                     "group_b", method="tobias")
    fai = params.get("genome_fai") or f"{genome}.fai"
    if not Path(str(fai)).is_file():
        return _skip(f"genome .fai index not found: {fai} (samtools faidx)",
                     method="tobias")

    out_dir.mkdir(parents=True, exist_ok=True)
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
    aggregate = {}
    if summary.get("parsed") and tobias.get("bindetect_outdir"):
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
            "Tn5-bias-corrected (ATACorrect) over per-condition merged replicate BAMs."
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
    p.add_argument("--barcode-groups", help="barcode<TAB>group TSV (scatac mode)")
    # bulk inputs (per-condition BAM lists).
    p.add_argument("--condition-bams",
                   help="JSON {condition: [bam, ...]} (bulk mode)")
    # shared.
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--peaks-bed", required=True)
    p.add_argument("--motif-meme", required=True)
    p.add_argument("--group-a", required=True)
    p.add_argument("--group-b", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--output-json", default=None)
    p.add_argument("--cores", type=int, default=8)
    args = p.parse_args(argv)

    shared = {
        "genome_fasta": args.genome_fasta, "peaks_bed": args.peaks_bed,
        "motif_meme": args.motif_meme, "group_a": args.group_a,
        "group_b": args.group_b, "output_dir": args.output_dir,
        "cores": args.cores,
    }
    if args.mode == "bulk":
        if not args.condition_bams:
            p.error("--condition-bams is required in bulk mode")
        condition_bams = json.loads(Path(args.condition_bams).read_text())
        res = chromatin_footprint_tobias_bulk(
            {**shared, "condition_bams": condition_bams})
    else:
        if not args.fragments_file or not args.barcode_groups:
            p.error("--fragments-file and --barcode-groups are required in scatac mode")
        res = chromatin_footprint_tobias({
            **shared, "fragments_file": args.fragments_file,
            "barcode_groups": args.barcode_groups,
        })
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(res, indent=2, sort_keys=True))
    print(json.dumps({k: v for k, v in res.items() if k != "group_bams"}, indent=2))
    return 0


if __name__ == "__main__":
    if len(sys.argv) == 3:
        run_script(chromatin_footprint_tobias_dispatch)
    else:
        raise SystemExit(main())
