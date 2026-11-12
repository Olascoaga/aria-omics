"""Production bulk ATAC peak-by-sample count matrix.

validation_level: production (stack="chromatin").

Counts aligned bulk ATAC reads/fragments overlapping an explicit peak universe.
This is a technical prerequisite for replicate-aware differential accessibility;
it does not run DA or infer biological groups from filenames.
"""

from __future__ import annotations

import csv
import os
import sys
from bisect import bisect_left, bisect_right
import shutil
import subprocess
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script

# Compatibility with the original A4 API; keep all manifests content-addressed.
try:
    from aria.utils.stage_manifest import (
        file_record, load_valid_stage_result, stage_manifest_path,
        write_stage_manifest,
    )
except ImportError:
    from aria.utils.stage_manifest import (
        build_stage_manifest, stage_is_current, write_stage_manifest,
    )

    def file_record(path):
        return str(Path(path).expanduser().resolve())

    def stage_manifest_path(output_dir, stage):
        return Path(output_dir) / f"{stage}.stage_manifest.json"

    def load_valid_stage_result(path, *, expected_inputs, expected_parameters):
        inputs = [(str(i), f) for i, f in enumerate(expected_inputs)]
        try:
            current, _ = stage_is_current(
                path, inputs=inputs,
                params={**expected_parameters, "_input_paths": expected_inputs},
                tool_version="ADR-077-track1-v1",
            )
            if not current:
                return None
            import json
            result = json.loads(Path(path).read_text()).get("extra", {}).get("result")
            return result if isinstance(result, dict) else None
        except (OSError, ValueError, TypeError, AttributeError):
            return None


def _save_stage_result(path, stage, inputs, parameters, result):
    from aria.utils.stage_manifest import build_stage_manifest
    manifest = build_stage_manifest(
        stage=stage, inputs=[(str(i), f) for i, f in enumerate(inputs)],
        params={**parameters, "_input_paths": inputs}, tool_version="ADR-077-track1-v1",
        extra={"result": result},
    )
    write_stage_manifest(path, manifest)


def _open_alignment(path, reference_fasta=None):
    import pysam
    return pysam.AlignmentFile(
        str(path), "rc" if str(path).lower().endswith(".cram") else "rb",
        **({"reference_filename": str(reference_fasta)} if reference_fasta else {}),
    )


def _interval_index(rows):
    """Merge intervals so a predecessor bisect handles nested/overlapping BEDs."""
    by_chrom = {}
    for chrom, start, end in rows:
        by_chrom.setdefault(chrom, []).append((start, end))
    index = {}
    for chrom, intervals in by_chrom.items():
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        index[chrom] = ([s for s, e in merged], [e for s, e in merged])
    return index


def _overlaps(index, chrom, start, end):
    starts, ends = index.get(chrom, ([], []))
    i = bisect_left(starts, end) - 1
    return i >= 0 and ends[i] > start


def _primary_mapped(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary)


def _count_cram(peaks_path, file_path, reference_fasta=None):
    """Count alignment overlaps, matching bedtools coverage -counts semantics.

    A prefix maximum of peak ends bounds the search even for nested intervals.
    """
    rows = []
    with open(peaks_path) as handle:
        for line in handle:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                cols = line.split("\t")
                rows.append((cols[0], int(cols[1]), int(cols[2])))
    by_chrom = {}
    for i, (chrom, start, end) in enumerate(rows):
        by_chrom.setdefault(chrom, []).append((start, end, i))
    index = {}
    for chrom, intervals in by_chrom.items():
        intervals.sort()
        max_end = []
        for start, end, i in intervals:
            max_end.append(max(end, max_end[-1] if max_end else end))
        index[chrom] = (intervals, [v[0] for v in intervals], max_end)
    counts = [0] * len(rows)
    with _open_alignment(file_path, reference_fasta) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_end is None:
                continue
            intervals, starts, max_end = index.get(read.reference_name, ([], [], []))
            lo = bisect_right(max_end, read.reference_start)
            hi = bisect_left(starts, read.reference_end)
            for start, end, i in intervals[lo:hi]:
                if end > read.reference_start:
                    counts[i] += 1
    return counts


def _error(error_type: str, details: str, **extra) -> dict:
    out = {"status": "error", "error_type": error_type, "details": details}
    out.update(extra)
    return out


def _strip_known_suffixes(path: str | Path) -> str:
    name = Path(path).name
    for suffix in (
        ".filtered.bam", ".markdup.bam", ".coordsort.bam",
        ".sorted.bam", ".dedup.bam", ".filtered.cram",
        ".bam", ".cram", ".bed", ".bed.gz",
        ".fragments.tsv.gz", ".fragments.tsv",
    ):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def _resolve_sample_ids(files: list[str],
                        sample_ids: list[str] | None = None) -> list[str]:
    if sample_ids is not None:
        ids = [str(s).strip() for s in sample_ids]
        if len(ids) != len(files):
            raise ValueError(
                "sample_ids length must match files length "
                f"({len(ids)} != {len(files)})."
            )
    else:
        ids = [_strip_known_suffixes(f) for f in files]
    if any(not s for s in ids):
        raise ValueError("sample IDs must be non-empty.")
    if len(set(ids)) != len(ids):
        raise ValueError("sample IDs must be unique; pass explicit sample_ids.")
    return ids


def _read_peak_ids(peaks_path: str | Path) -> list[str]:
    ids: list[str] = []
    seen: dict[str, int] = {}
    with open(peaks_path, "r", newline="") as handle:
        for lineno, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                raise ValueError(
                    f"Peak file line {lineno} has fewer than 3 BED columns."
                )
            chrom = cols[0]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError as exc:
                raise ValueError(
                    f"Peak file line {lineno} has non-integer coordinates."
                ) from exc
            if end <= start:
                raise ValueError(
                    f"Peak file line {lineno} has end <= start."
                )
            base = f"{chrom}:{start}-{end}"
            n = seen.get(base, 0) + 1
            seen[base] = n
            ids.append(base if n == 1 else f"{base}#{n}")
    if not ids:
        raise ValueError("Peak file contains no usable BED intervals.")
    return ids


def _parse_bedtools_counts(stdout: str, peak_ids: list[str]) -> list[int]:
    counts: list[int] = []
    for line in stdout.splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        cols = line.rstrip("\n").split("\t")
        try:
            counts.append(int(float(cols[-1])))
        except (ValueError, IndexError) as exc:
            raise ValueError(
                "bedtools coverage output ended with a non-numeric count."
            ) from exc
    if len(counts) != len(peak_ids):
        raise ValueError(
            "bedtools coverage row count did not match peak file "
            f"({len(counts)} != {len(peak_ids)})."
        )
    return counts


def _write_counts_matrix(path: Path, peak_ids: list[str],
                         sample_counts: dict[str, list[int]]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        samples = list(sample_counts)
        writer.writerow(["peak_id", *samples])
        for idx, peak_id in enumerate(peak_ids):
            writer.writerow([peak_id, *[sample_counts[s][idx] for s in samples]])


def _metadata_rows(files: list[str], sample_ids: list[str],
                   sample_metadata: Any = None) -> list[dict[str, Any]]:
    rows = [
        {"sample_id": sample_id, "file": str(file)}
        for sample_id, file in zip(sample_ids, files)
    ]
    if sample_metadata is None:
        return rows
    if isinstance(sample_metadata, dict):
        meta_by_sample = sample_metadata
    elif isinstance(sample_metadata, list):
        meta_by_sample = {
            str(row.get("sample_id")): row
            for row in sample_metadata
            if isinstance(row, dict) and row.get("sample_id") is not None
        }
    else:
        raise ValueError("sample_metadata must be a dict or a list of dicts.")

    for row in rows:
        meta = meta_by_sample.get(row["sample_id"]) or {}
        if not isinstance(meta, dict):
            raise ValueError("sample_metadata entries must be dicts.")
        for key, value in meta.items():
            if key not in {"sample_id", "file"}:
                row[str(key)] = value
    return rows


def _write_metadata(path: Path, rows: list[dict[str, Any]]) -> None:
    extra = sorted({k for row in rows for k in row} - {"sample_id", "file"})
    fields = ["sample_id", "file", *extra]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _genome_sizes_from_bam(file_path: str | Path,
                           reference_fasta: str | None = None) -> list[tuple[str, int]] | None:
    """Read (chrom, length) in header order from a BAM/CRAM via pysam.

    Used to drive `bedtools coverage -sorted` (low-memory sweeping algorithm).
    Returns None when the header cannot be read (e.g. pysam unavailable, fragment
    TSV input, or an unreadable file) so the caller falls back to the plain call.
    """
    try:
        with _open_alignment(file_path, reference_fasta) as bam:
            refs = list(bam.references or [])
            lens = list(bam.lengths or [])
    except Exception:
        return None
    if refs and lens and len(refs) == len(lens):
        return [(str(c), int(l)) for c, l in zip(refs, lens)]
    return None


def _write_genome_file(path: Path, genome_sizes: list[tuple[str, int]]) -> None:
    with open(path, "w", newline="") as handle:
        for chrom, length in genome_sizes:
            handle.write(f"{chrom}\t{int(length)}\n")


def _sort_peaks_by_genome(peaks_path: str | Path, genome_order: list[str],
                          out_path: Path) -> None:
    """Sort a peak file by the alignment's chromosome order, then start/end.

    `bedtools coverage -sorted` requires `-a` to share the chromosome ordering of
    `-g`. Malformed lines are preserved (not dropped) so `_read_peak_ids` still
    raises on them honestly rather than silently shrinking the peak universe.
    """
    order = {chrom: idx for idx, chrom in enumerate(genome_order)}
    rows = []
    with open(peaks_path, "r", newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            chrom = cols[0] if cols else ""
            try:
                start, end = int(cols[1]), int(cols[2])
            except (IndexError, ValueError):
                start, end = 0, 0
            rows.append((order.get(chrom, len(order)), start, end,
                         line.rstrip("\n")))
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    with open(out_path, "w", newline="") as handle:
        for row in rows:
            handle.write(row[3] + "\n")


def _count_one_sample(bedtools: str, peaks_path: Path, file_path: Path,
                      peak_ids: list[str],
                      genome_file: Path | None = None) -> list[int]:
    cmd = [bedtools, "coverage"]
    if genome_file is not None:
        # Sweeping algorithm: streams one chromosome at a time instead of loading
        # the whole BAM into memory (which OOMs on large bulk ATAC libraries).
        cmd += ["-sorted", "-g", str(genome_file)]
    cmd += ["-counts", "-a", str(peaks_path), "-b", str(file_path)]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "bedtools coverage failed for "
            f"{file_path}: {(proc.stderr or proc.stdout).strip()[:500]}"
        )
    return _parse_bedtools_counts(proc.stdout, peak_ids)


def chromatin_peak_counts(params: dict) -> dict:
    data_type = params.get("data_type", "bulk_ATAC")
    if data_type != "bulk_ATAC":
        return _error(
            "UnsupportedDataType",
            "chromatin_peak_counts currently supports bulk_ATAC only.",
            data_type=data_type,
        )

    files = [str(f) for f in (params.get("files") or [])]
    peaks_value = (
        params.get("peaks_path")
        or params.get("consensus_peaks_path")
        or params.get("peak_universe_path")
    )
    if not files:
        return _error("MissingInput", "files must contain at least one BAM/CRAM.")
    if not peaks_value:
        return _error("MissingInput", "peaks_path is required.")

    peaks_path = Path(str(peaks_value)).expanduser()
    if not peaks_path.exists():
        return _error("FileNotFound", f"Peak file not found: {peaks_path}")
    missing = [f for f in files if not Path(f).expanduser().exists()]
    if missing:
        return _error("FileNotFound", f"Input alignment file(s) not found: {missing}")

    output_dir = Path(str(params.get("output_dir") or
                          peaks_path.parent / "bulk_atac_counts"))
    output_dir.mkdir(parents=True, exist_ok=True)
    warnings = []

    try:
        sample_ids = _resolve_sample_ids(files, params.get("sample_ids"))
        reference_fasta = params.get("reference_fasta")
        manifest_path = stage_manifest_path(output_dir, "chromatin_peak_counts")
        inputs = [file_record(f) for f in [*files, str(peaks_path)]]
        if reference_fasta:
            inputs.append(file_record(reference_fasta))
        parameters = {"files": files, "peaks_path": str(peaks_path),
                      "sample_ids": sample_ids,
                      "sample_metadata": params.get("sample_metadata"),
                      "reference_fasta": reference_fasta}
        cached = load_valid_stage_result(
            manifest_path, expected_inputs=inputs, expected_parameters=parameters)
        if cached and all(Path(cached.get(k, "")).is_file() for k in (
                "counts_matrix_path", "counts_cpm_path", "sample_metadata_path")):
            return {**cached, "reused_from_stage_manifest": True}
        Path(manifest_path).unlink(missing_ok=True)
        bedtools = shutil.which(str(params.get("bedtools") or "bedtools"))
        if not bedtools and any(not f.lower().endswith(".cram") for f in files):
            return _error("MissingDependency",
                          "bedtools is required for bulk ATAC BAM peak counting.")
        # Prefer the low-memory sweeping algorithm. Without `-sorted`, bedtools
        # loads the whole BAM into memory and OOMs on large bulk ATAC libraries.
        genome_sizes = (_genome_sizes_from_bam(Path(files[0]).expanduser(), reference_fasta)
                        if files[0].lower().endswith(".cram") else
                        _genome_sizes_from_bam(Path(files[0]).expanduser()))
        genome_file = None
        count_peaks_path = peaks_path
        if genome_sizes:
            genome_file = output_dir / "genome.txt"
            _write_genome_file(genome_file, genome_sizes)
            count_peaks_path = output_dir / "peaks.genome_sorted.bed"
            _sort_peaks_by_genome(peaks_path, [c for c, _ in genome_sizes],
                                  count_peaks_path)
        else:
            warnings.append(
                "Could not read a genome from the alignment header; counting "
                "with the in-memory bedtools algorithm (higher memory use)."
            )
        peak_ids = _read_peak_ids(count_peaks_path)
        sample_counts = {}
        for sample_id, file_path in zip(sample_ids, files):
            if file_path.lower().endswith(".cram"):
                sample_counts[sample_id] = _count_cram(
                    count_peaks_path, file_path, reference_fasta)
                continue
            sample_counts[sample_id] = _count_one_sample(
                bedtools, count_peaks_path, Path(file_path).expanduser(),
                peak_ids, genome_file=genome_file)
        metadata = _metadata_rows(files, sample_ids,
                                  params.get("sample_metadata"))
    except Exception as exc:
        return _error(type(exc).__name__, str(exc))

    matrix_path = output_dir / "bulk_atac_peak_counts.tsv"
    metadata_path = output_dir / "bulk_atac_samples.tsv"
    cpm_path = output_dir / "bulk_atac_peak_cpm.tsv"
    cpm = {}
    for sid, counts in sample_counts.items():
        total = max(sum(counts), 1)
        cpm[sid] = [format(c / total * 1e6, ".10g") for c in counts]
    _write_counts_matrix(cpm_path, peak_ids, cpm)
    _write_counts_matrix(matrix_path, peak_ids, sample_counts)
    _write_metadata(metadata_path, metadata)

    has_condition = any("condition" in row for row in metadata)
    has_replicate = any("replicate" in row or "replicate_id" in row
                        for row in metadata)
    if not (has_condition and has_replicate):
        warnings.append(
            "Peak-count matrix was built, but condition and replicate metadata "
            "were not both supplied; differential accessibility remains gated."
        )

    result = {
        "status": "success",
        "data_type": "bulk_ATAC",
        "validation_level": "production",
        "analysis": "peak_count_matrix",
        "counting_method": ("pysam CRAM / bedtools coverage -counts"
                            if any(f.lower().endswith(".cram") for f in files)
                            else "bedtools coverage -sorted -counts"
                            if genome_file is not None
                            else "bedtools coverage -counts"),
        "counts_matrix_path": str(matrix_path),
        "counts_cpm_path": str(cpm_path),
        "reused_from_stage_manifest": False,
        "sample_metadata_path": str(metadata_path),
        "peaks_path": str(peaks_path),
        "n_peaks": len(peak_ids),
        "n_samples": len(sample_ids),
        "sample_ids": sample_ids,
        "warnings": warnings,
    }

    _save_stage_result(manifest_path, "chromatin_peak_counts", inputs, parameters, result)
    return result


if __name__ == "__main__":
    run_script(chromatin_peak_counts)
