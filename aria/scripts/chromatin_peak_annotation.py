"""B2 (bulk ATAC pre-print): genomic annotation of differential-accessibility peaks.

State-of-the-art ATAC interpretation annotates each peak with its genomic context
(ChIPseeker / HOMER `annotatePeaks` style): the nearest gene by TSS, the signed
distance to that TSS (strand-aware; negative = upstream), and a genomic feature
class — Promoter / Exonic / Intronic / Distal Intergenic.

Design:
- Pure, deterministic core (`parse_gtf_gene_models`, `annotate_peaks`); no network,
  no fabrication. The feature is assigned relative to the nearest-TSS gene within a
  TSS-sorted neighborhood (the standard nearest-gene convention), so a peak whose
  nearest TSS has no locatable gene model is reported as Distal Intergenic, never
  imputed.
- The dispatchable entry (`chromatin_peak_annotation`) auto-resolves a GTF from the
  standard ARIA locations (`~/.aria/genomes/*/annotation.gtf{,.gz}`) or an explicit
  hint / `ARIA_GTF`, and returns an HONEST skip when none is found — annotation is a
  positional/structural fact, so without a gene model there is nothing to assert.

Annotation is associative (a peak's proximity to a gene), NOT evidence that the gene
is regulated by that peak.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional

from aria.benchmarks.bulk_atac_marker_mapping import parse_peak_interval
from aria.scripts.rna_bulk.gtf_io import _locate_gtf
from aria.utils.csv_artifacts import persist_csv

_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_GENE_NAME = re.compile(r'gene_name "([^"]+)"')

# Feature priority (highest first) when more than one candidate gene applies.
FEATURE_PROMOTER = "Promoter"
FEATURE_EXONIC = "Exonic"
FEATURE_INTRONIC = "Intronic"
FEATURE_INTERGENIC = "Distal Intergenic"

DEFAULT_PROMOTER_UPSTREAM = 2000
DEFAULT_PROMOTER_DOWNSTREAM = 2000


def _norm_chrom(chrom) -> str:
    """Reconcile UCSC (`chr1`/`chrX`/`chrM`) and Ensembl (`1`/`X`/`MT`) contig naming
    so chr-prefixed ATAC peaks match an Ensembl GTF (and vice versa)."""
    c = str(chrom)
    if c[:3].lower() == "chr":
        c = c[3:]
    return "MT" if c in ("M", "MT") else c
# How many TSS-sorted neighbors to consider as feature candidates around a peak.
# Bounds the work per peak while still catching a long gene body whose TSS is not
# the single closest one.
_NEIGHBORS = 25


def parse_gtf_gene_models(gtf_path) -> dict:
    """Stream a GTF into per-chromosome gene models (low memory).

    Returns ``{chrom: [gene_dict, ...]}`` sorted by TSS, where each ``gene_dict`` has
    ``gene`` (symbol or id), ``strand``, ``start``, ``end``, ``tss`` (1-based), and
    ``exons`` (sorted ``[(start, end), ...]``). Mirrors the parsing conventions in
    ``rna_bulk/gtf_io``. Returns ``{}`` on a read/parse failure (honest, never partial
    fabrication).
    """
    genes: dict = {}
    exons = defaultdict(list)
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    try:
        with opener(gtf_path, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                feature = fields[2]
                if feature not in ("gene", "exon"):
                    continue
                attrs = fields[8]
                gid_m = _GENE_ID.search(attrs)
                if not gid_m:
                    continue
                gid = gid_m.group(1)
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError:
                    continue
                if feature == "gene":
                    strand = fields[6]
                    name_m = _GENE_NAME.search(attrs)
                    genes[gid] = {
                        "gene": name_m.group(1) if name_m else gid.split(".")[0],
                        "chrom": _norm_chrom(fields[0]),
                        "strand": strand,
                        "start": start,
                        "end": end,
                        "tss": start if strand != "-" else end,
                    }
                else:  # exon
                    exons[gid].append((start, end))
    except Exception:
        return {}

    by_chrom: dict = defaultdict(list)
    for gid, gene in genes.items():
        gene["exons"] = sorted(exons.get(gid, []))
        by_chrom[gene["chrom"]].append(gene)
    out: dict = {}
    for chrom, glist in by_chrom.items():
        glist.sort(key=lambda g: g["tss"])
        out[chrom] = glist
    return out


def _signed_tss_distance(center: int, gene: dict) -> int:
    """Signed distance from a peak center to a gene TSS (negative = upstream)."""
    if gene["strand"] == "-":
        return gene["tss"] - center
    return center - gene["tss"]


def _in_exon(center: int, exons: list) -> bool:
    for s, e in exons:
        if s <= center <= e:
            return True
        if s > center:
            break
    return False


def annotate_one_peak(chrom: str, start: int, end: int, models_by_chrom: dict,
                      tss_index: dict, promoter_upstream: int,
                      promoter_downstream: int) -> dict:
    """Annotate a single peak against the gene models. Pure + deterministic.

    Returns ``{nearest_gene, distance_to_tss, feature}``. ``nearest_gene`` / distance
    are ``None`` only when the chromosome has no gene model (honest intergenic)."""
    chrom = _norm_chrom(chrom)
    glist = models_by_chrom.get(chrom)
    if not glist:
        return {"nearest_gene": None, "distance_to_tss": None,
                "feature": FEATURE_INTERGENIC}
    center = (start + end) // 2
    pos = bisect.bisect_left(tss_index[chrom], center)
    lo = max(0, pos - _NEIGHBORS)
    hi = min(len(glist), pos + _NEIGHBORS)
    cands = glist[lo:hi] or glist

    nearest = min(cands, key=lambda g: abs(center - g["tss"]))

    # Priority: Promoter > Exonic > Intronic > Distal Intergenic.
    for g in cands:
        sdist = _signed_tss_distance(center, g)
        if -promoter_upstream <= sdist <= promoter_downstream:
            return {"nearest_gene": g["gene"],
                    "distance_to_tss": _signed_tss_distance(center, g),
                    "feature": FEATURE_PROMOTER}
    for g in cands:
        if g["start"] <= center <= g["end"]:
            feature = FEATURE_EXONIC if _in_exon(center, g["exons"]) \
                else FEATURE_INTRONIC
            return {"nearest_gene": g["gene"],
                    "distance_to_tss": _signed_tss_distance(center, g),
                    "feature": feature}
    return {"nearest_gene": nearest["gene"],
            "distance_to_tss": _signed_tss_distance(center, nearest),
            "feature": FEATURE_INTERGENIC}


def annotate_peaks(peaks, models_by_chrom: dict,
                   promoter_upstream: int = DEFAULT_PROMOTER_UPSTREAM,
                   promoter_downstream: int = DEFAULT_PROMOTER_DOWNSTREAM) -> dict:
    """Annotate an iterable of peak ids. Returns
    ``{rows: [{peak, nearest_gene, distance_to_tss, feature}], feature_distribution,
    n_annotated, n_unparsed}``. Unparseable peak ids are skipped + counted (never
    fabricated)."""
    tss_index = {c: [g["tss"] for g in glist]
                 for c, glist in models_by_chrom.items()}
    rows: list = []
    dist: dict = defaultdict(int)
    n_unparsed = 0
    for peak in peaks:
        parsed = parse_peak_interval(peak)
        if parsed is None:
            n_unparsed += 1
            continue
        chrom, start, end = parsed
        ann = annotate_one_peak(chrom, start, end, models_by_chrom, tss_index,
                                promoter_upstream, promoter_downstream)
        rows.append({"peak": peak, **ann})
        dist[ann["feature"]] += 1
    return {"rows": rows, "feature_distribution": dict(dist),
            "n_annotated": len(rows), "n_unparsed": n_unparsed}


def _resolve_gtf(params: dict, out_dir: Path) -> Optional[Path]:
    """Auto-resolve a GTF: explicit ``gtf`` hint > ``ARIA_GTF`` env >
    assembly-specific ``~/.aria/genomes/<genome>/annotation.gtf{,.gz}`` > generic
    standard-location search. The assembly-specific step matters: the generic
    `_locate_gtf` returns whatever genome dir it finds first, so without it a human
    dataset can be annotated against a mouse GTF. Power-user override stays a
    fallback, per the auto-resolution principle."""
    hint = params.get("gtf") or params.get("gtf_path") or os.environ.get("ARIA_GTF")
    if hint and Path(hint).exists():
        return Path(hint)
    # Prefer the GTF for THIS experiment's assembly before the generic search.
    genome = params.get("genome") or params.get("assembly")
    if genome:
        genome_dir = Path.home() / ".aria" / "genomes" / str(genome)
        for name in ("annotation.gtf", "annotation.gtf.gz"):
            cand = genome_dir / name
            if cand.exists() and cand.stat().st_size > 100_000:
                return cand
    return _locate_gtf(out_dir / "peaks", hint)


def _significant_peaks(csv_path: str):
    """Yield significant peak ids from a bulk ATAC DA full CSV (`significant==True`)."""
    import pandas as pd

    df = pd.read_csv(csv_path)
    peak_col = "peak" if "peak" in df.columns else df.columns[0]
    if "significant" in df.columns:
        df = df[df["significant"].fillna(False).astype(bool)]
    return [str(p) for p in df[peak_col].tolist()]


def chromatin_peak_annotation(params: dict) -> dict:
    """Dispatchable entry: annotate the significant DA peaks of each bulk ATAC
    comparison against an auto-resolved GTF. Honest-skip without a GTF or peaks."""
    out_dir = Path(params.get("output_dir", "."))
    out_dir.mkdir(parents=True, exist_ok=True)
    _up = params.get("promoter_upstream")
    _down = params.get("promoter_downstream")
    promoter_up = int(_up) if _up is not None else DEFAULT_PROMOTER_UPSTREAM
    promoter_down = int(_down) if _down is not None else DEFAULT_PROMOTER_DOWNSTREAM

    gtf_path = _resolve_gtf(params, out_dir)
    if not gtf_path:
        return {
            "status": "skipped", "ran": False,
            "analysis": "peak_annotation", "data_type": "bulk_ATAC",
            "validation_level": "beta",
            "reason": "no_gtf_for_peak_annotation",
            "details": ("No GTF found. Place one at "
                        "~/.aria/genomes/<assembly>/annotation.gtf[.gz], pass "
                        "gtf=<path>, or set ARIA_GTF=<path>."),
        }

    models = parse_gtf_gene_models(str(gtf_path))
    if not models:
        return {
            "status": "skipped", "ran": False,
            "analysis": "peak_annotation", "data_type": "bulk_ATAC",
            "validation_level": "beta",
            "reason": "gtf_parsed_zero_gene_models",
            "gtf": Path(gtf_path).name,
        }

    comparisons_in = params.get("comparisons") or []
    warnings: list = []
    comp_results: list = []
    overall: dict = defaultdict(int)
    for comp in comparisons_in:
        csv_path = (comp or {}).get("full_results_csv")
        if not csv_path or not Path(csv_path).exists():
            continue
        peaks = _significant_peaks(csv_path)
        if not peaks:
            continue
        ann = annotate_peaks(peaks, models, promoter_up, promoter_down)
        comp_key = f"{comp.get('test')}_vs_{comp.get('reference')}"
        ann_csv = out_dir / f"bulk_atac_peak_annotation_{comp_key}.csv"

        # persist_csv calls write_fn() with no args; bind the per-iteration path +
        # rows via defaults to avoid late-binding over the loop variables.
        def _write(_path=ann_csv, _rows=ann["rows"]):
            with open(_path, "w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["peak", "nearest_gene", "distance_to_tss",
                                "feature"])
                writer.writeheader()
                writer.writerows(_rows)

        written = persist_csv(str(ann_csv), f"peak annotation [{comp_key}]",
                              _write, warnings)
        for feat, n in ann["feature_distribution"].items():
            overall[feat] += n
        comp_results.append({
            "test": comp.get("test"), "reference": comp.get("reference"),
            "n_annotated": ann["n_annotated"], "n_unparsed": ann["n_unparsed"],
            "feature_distribution": ann["feature_distribution"],
            "annotation_csv": written,
        })

    if not comp_results:
        return {
            "status": "skipped", "ran": False,
            "analysis": "peak_annotation", "data_type": "bulk_ATAC",
            "validation_level": "beta",
            "reason": "no_significant_peaks_to_annotate",
            "gtf": Path(gtf_path).name, "warnings": warnings,
        }

    return {
        "status": "success", "ran": True,
        "analysis": "peak_annotation", "data_type": "bulk_ATAC",
        "validation_level": "beta",
        "method": ("nearest-TSS genomic annotation (Promoter/Exonic/Intronic/"
                   "Distal Intergenic) over GTF gene models"),
        "gtf": Path(gtf_path).name,
        "promoter_upstream": promoter_up,
        "promoter_downstream": promoter_down,
        "comparisons": comp_results,
        "feature_distribution_overall": dict(overall),
        "warnings": warnings,
    }


if __name__ == "__main__":
    import json
    import sys

    print(json.dumps(chromatin_peak_annotation(json.loads(sys.argv[1]))))
