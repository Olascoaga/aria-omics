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
from pathlib import Path
from typing import Any, Iterable, Iterator

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
    TFs toward each group (by the ``<A>_<B>_change`` score, gated on pvalue + site
    count) plus counts. Differential TF binding is ASSOCIATIVE, never causal. Returns
    an honest ``parsed: false`` when the expected columns are absent."""
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

    def _sig(rs):
        return [x for x in rs if (x[2] == x[2]) and x[2] < pvalue_max and x[3] >= min_sites]

    sig = _sig(rows)
    fmt = lambda x: {"tf": x[0], "change": round(x[1], 4), "pvalue": x[2], "n_sites": x[3]}
    return {
        "parsed": True,
        "n_motifs_tested": len(rows),
        "n_significant": len(sig),
        f"top_toward_{group_a}": [fmt(x) for x in sorted(sig, key=lambda x: -x[1])[:top_n]],
        f"top_toward_{group_b}": [fmt(x) for x in sorted(sig, key=lambda x: x[1])[:top_n]],
    }


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
    return {
        "ran": True,
        "method": "tobias_atacorrect_scorebigwig_bindetect",
        "group_a": group_a, "group_b": group_b,
        "group_bams": group_bams,
        "differential_summary": summary,
        **tobias,
        "caveat": (
            "Differential TF binding (BINDetect) is an associative footprint-signal "
            "difference between cell-type groups, not causal regulation. Footprints "
            "are Tn5-bias-corrected (ATACorrect); single-cell footprinting is noisy "
            "so groups are pseudobulk."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fragments-file", required=True)
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--peaks-bed", required=True)
    p.add_argument("--motif-meme", required=True)
    p.add_argument("--barcode-groups", required=True, help="barcode<TAB>group TSV")
    p.add_argument("--group-a", required=True)
    p.add_argument("--group-b", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--output-json", default=None)
    p.add_argument("--cores", type=int, default=8)
    args = p.parse_args(argv)

    res = chromatin_footprint_tobias({
        "fragments_file": args.fragments_file, "genome_fasta": args.genome_fasta,
        "peaks_bed": args.peaks_bed, "motif_meme": args.motif_meme,
        "barcode_groups": args.barcode_groups, "group_a": args.group_a,
        "group_b": args.group_b, "output_dir": args.output_dir, "cores": args.cores,
    })
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(res, indent=2, sort_keys=True))
    print(json.dumps({k: v for k, v in res.items() if k != "group_bams"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
