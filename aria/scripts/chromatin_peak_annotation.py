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
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

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


def _exon_overlaps_interval(start: int, end: int, exons: list) -> bool:
    if end <= start:
        return False
    for s, e in exons:
        if start < e and end > s:
            return True
        if s >= end:
            break
    return False


def _build_body_index(models_by_chrom: dict) -> dict:
    """Build per-chromosome start-sorted gene lists with prefix max end so
    finding every gene overlapping a peak interval/midpoint is exact even for
    multi-Mb genes spanning hundreds of nested shorter genes."""
    index: dict = {}
    for chrom, glist in models_by_chrom.items():
        by_start = sorted(glist, key=lambda g: (g["start"], g["end"]))
        starts = [g["start"] for g in by_start]
        prefix_max_end: list[int] = []
        running_max = -1
        for g in by_start:
            if g["end"] > running_max:
                running_max = g["end"]
            prefix_max_end.append(running_max)
        index[chrom] = (by_start, starts, prefix_max_end)
    return index


def _overlapping_genes(
    chrom: str,
    start: int,
    end: int,
    center: int,
    models_by_chrom: dict,
    body_index: Optional[dict] = None,
) -> list[dict]:
    """Return every gene on ``chrom`` whose body contains ``center`` or
    physically overlaps ``[start, end)``."""
    hi_coord = max(center, end)
    lo_coord = min(center, start)
    out: list[dict] = []
    if body_index is not None and chrom in body_index:
        by_start, starts, prefix_max_end = body_index[chrom]
        idx = bisect.bisect_right(starts, hi_coord) - 1
        while idx >= 0 and prefix_max_end[idx] >= lo_coord:
            g = by_start[idx]
            if (g["start"] <= center <= g["end"]) or (
                end > start and start < g["end"] and end > g["start"]
            ):
                out.append(g)
            idx -= 1
        return out
    for g in models_by_chrom.get(chrom) or []:
        if (g["start"] <= center <= g["end"]) or (
            end > start and start < g["end"] and end > g["start"]
        ):
            out.append(g)
    return out


def annotate_one_peak(chrom: str, start: int, end: int, models_by_chrom: dict,
                      tss_index: dict, promoter_upstream: int,
                      promoter_downstream: int,
                      body_index: Optional[dict] = None) -> dict:
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
    tss_cands = glist[lo:hi] or glist

    def _tss_key(g: dict):
        return (abs(center - g["tss"]), g["tss"], str(g["gene"]))

    # Candidates ordered by |center - TSS| (closest first; ties broken by TSS then
    # gene name for determinism), NOT by genomic TSS order -- otherwise a
    # bidirectional/overlapping promoter resolves to the leftmost TSS.
    tss_cands = sorted(tss_cands, key=_tss_key)
    nearest = tss_cands[0]

    # Physically overlapping genes across the whole chromosome (catches long
    # genes whose TSS is > _NEIGHBORS away as well as interval-overlapping exons).
    body_all = _overlapping_genes(
        chrom, start, end, center, models_by_chrom, body_index=body_index
    )
    body_all.sort(key=_tss_key)

    prom_cands = list(tss_cands)
    if body_all:
        seen_ids = {id(g) for g in prom_cands}
        for g in body_all:
            if id(g) not in seen_ids:
                prom_cands.append(g)
                seen_ids.add(id(g))
        prom_cands.sort(key=_tss_key)

    # Strict priority across ALL candidates: Promoter > Exonic > Intronic >
    # Distal Intergenic (an exon or gene body physically overlapping the peak
    # always beats a closer non-overlapping gene's Distal Intergenic call).
    for g in prom_cands:
        sdist = _signed_tss_distance(center, g)
        if -promoter_upstream <= sdist <= promoter_downstream:
            return {"nearest_gene": g["gene"],
                    "distance_to_tss": sdist,
                    "feature": FEATURE_PROMOTER}
    for g in body_all:
        if _in_exon(center, g["exons"]):
            return {"nearest_gene": g["gene"],
                    "distance_to_tss": _signed_tss_distance(center, g),
                    "feature": FEATURE_EXONIC}
    for g in body_all:
        if _exon_overlaps_interval(start, end, g["exons"]):
            return {"nearest_gene": g["gene"],
                    "distance_to_tss": _signed_tss_distance(center, g),
                    "feature": FEATURE_EXONIC}
    body_center = [g for g in body_all if g["start"] <= center <= g["end"]]
    if body_center:
        g = body_center[0]
        return {"nearest_gene": g["gene"],
                "distance_to_tss": _signed_tss_distance(center, g),
                "feature": FEATURE_INTRONIC}
    if body_all:
        g = body_all[0]
        return {"nearest_gene": g["gene"],
                "distance_to_tss": _signed_tss_distance(center, g),
                "feature": FEATURE_INTRONIC}
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
    body_index = _build_body_index(models_by_chrom)
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
                                promoter_upstream, promoter_downstream,
                                body_index=body_index)
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
        # An explicit assembly NEVER falls through to the generic search: that
        # returns the first GTF under any genome dir and can silently annotate
        # mouse peaks with a human GTF (or vice versa). Honest None instead.
        return _assembly_gtf(str(genome))
    return _locate_gtf(out_dir / "peaks", hint)


# UCSC <-> GENCODE/Ensembl assembly aliases (same assembly, different names).
_ASSEMBLY_ALIASES = {
    "hg38": {"hg38", "grch38"}, "grch38": {"hg38", "grch38"},
    "hg19": {"hg19", "grch37"}, "grch37": {"hg19", "grch37"},
    "mm10": {"mm10", "grcm38"}, "grcm38": {"mm10", "grcm38"},
    "mm39": {"mm39", "grcm39"}, "grcm39": {"mm39", "grcm39"},
}


def _assembly_gtf(genome: str) -> Optional[Path]:
    """GTF for exactly this assembly (alias-aware, case-insensitive) under the
    managed genome store(s): ``ARIA_GENOME_DIR`` / ``$ARIA_HOME/genomes`` and
    ``~/.aria/genomes``. Returns None when absent (no cross-assembly fallback)."""
    key = str(genome).strip().lower()
    names = _ASSEMBLY_ALIASES.get(key, {key}) if key else set()
    if not names:
        return None
    bases = []
    override = os.environ.get("ARIA_GENOME_DIR")
    if override:
        bases.append(Path(override).expanduser())
    home = os.environ.get("ARIA_HOME")
    if home:
        bases.append(Path(home).expanduser() / "genomes")
    bases.append(Path.home() / ".aria" / "genomes")
    seen = set()
    for base in bases:
        try:
            subdirs = [d for d in base.iterdir() if d.is_dir()] if base.is_dir() else []
        except OSError:
            continue
        for d in sorted(subdirs, key=lambda x: x.name):
            if d.name.lower() not in names or d in seen:
                continue
            seen.add(d)
            for name in ("annotation.gtf", "annotation.gtf.gz"):
                cand = d / name
                if cand.exists() and cand.stat().st_size > 100_000:
                    return cand
    return None


_PRODUCTION_MODALITIES = {
    "scATAC": "scatac",
    "bulk_ATAC": "bulk_atac",
    "ChIP": "chip",
    "CUT_AND_RUN": "cut_and_run",
    "CUT_AND_TAG": "cut_and_tag",
}


def _modality_stamp(params: dict) -> tuple:
    """(data_type, validation_level, filename_prefix) for this dispatch.

    Defaults to the bulk ATAC stamp so every pre-existing bulk call (which never
    passes ``data_type``) keeps byte-identical output names and stamps."""
    data_type = str(params.get("data_type") or "bulk_ATAC")
    if data_type in _PRODUCTION_MODALITIES:
        return data_type, "production", _PRODUCTION_MODALITIES[data_type]
    safe = re.sub(r"[^A-Za-z0-9]+", "_", data_type).strip("_").lower() or "bulk_atac"
    return data_type, "beta", safe


def _select_comparison_rows(df, comp: Optional[dict]):
    """Restrict a DA table to one contrast when it is a combined multi-contrast
    table (a ``comparison`` column, as in the scATAC pseudobulk full CSV).

    A multi-contrast table (>1 distinct ``comparison``) is ALWAYS filtered to
    ``f"{test}_vs_{reference}"`` -- possibly to zero rows; it never falls back to
    the other contrasts' peaks. A single-contrast table (no ``comparison``
    column, or one distinct value) is returned unchanged when the key does not
    match, so per-comparison CSVs are never emptied."""
    comp = comp or {}
    test, ref = comp.get("test"), comp.get("reference")
    if "comparison" not in df.columns or test is None or ref is None:
        return df
    comp_key = f"{test}_vs_{ref}"
    labels = df["comparison"].astype(str)
    mask = labels == comp_key
    if bool(mask.any()) or labels.nunique(dropna=True) > 1:
        return df[mask]
    return df


def _significant_peaks(csv_path: str, comp: Optional[dict] = None):
    """Yield significant peak ids from a DA full CSV (`significant==True`),
    restricted to ``comp``'s contrast when the table is multi-contrast."""
    import pandas as pd

    df = _select_comparison_rows(pd.read_csv(csv_path), comp)
    peak_col = "peak" if "peak" in df.columns else df.columns[0]
    if "significant" in df.columns:
        df = df[df["significant"].fillna(False).astype(bool)]
    return [str(p) for p in df[peak_col].tolist()]


def chromatin_peak_annotation(params: dict) -> dict:
    """Dispatchable entry: annotate the significant DA peaks of each ATAC
    comparison (bulk ATAC DA, or scATAC pseudobulk DA when ``data_type`` is
    ``"scATAC"``) against an auto-resolved GTF. Honest-skip without a GTF or
    peaks."""
    data_type, validation_level, prefix = _modality_stamp(params)
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
            "analysis": "peak_annotation", "data_type": data_type,
            "validation_level": validation_level,
            "reason": "no_gtf_for_peak_annotation",
            "details": ("No GTF found. Place one at "
                        "~/.aria/genomes/<assembly>/annotation.gtf[.gz], pass "
                        "gtf=<path>, or set ARIA_GTF=<path>."),
        }

    models = parse_gtf_gene_models(str(gtf_path))
    if not models:
        return {
            "status": "skipped", "ran": False,
            "analysis": "peak_annotation", "data_type": data_type,
            "validation_level": validation_level,
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
        peaks = _significant_peaks(csv_path, comp)
        if not peaks:
            continue
        ann = annotate_peaks(peaks, models, promoter_up, promoter_down)
        comp_key = f"{comp.get('test')}_vs_{comp.get('reference')}"
        safe_key = re.sub(r"[^A-Za-z0-9._-]+", "_", comp_key)
        ann_csv = out_dir / f"{prefix}_peak_annotation_{safe_key}.csv"

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
            "analysis": "peak_annotation", "data_type": data_type,
            "validation_level": validation_level,
            "reason": "no_significant_peaks_to_annotate",
            "gtf": Path(gtf_path).name, "warnings": warnings,
        }

    return {
        "status": "success", "ran": True,
        "analysis": "peak_annotation", "data_type": data_type,
        "validation_level": validation_level,
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

    if len(sys.argv) == 2 and sys.argv[1].lstrip().startswith("{"):
        print(json.dumps(chromatin_peak_annotation(json.loads(sys.argv[1]))))
    else:
        from aria.scripts._base import run_script

        run_script(chromatin_peak_annotation)
