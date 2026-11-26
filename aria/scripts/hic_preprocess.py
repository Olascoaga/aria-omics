"""Universal Hi-C / Micro-C preprocessing & cooler intake (ADR-082).

validation_level: production (stack="hic").

Normalizes heterogeneous 3D genome inputs into canonical Cooler containers
(``.cool`` / ``.mcool``) using strictly out-of-core streaming:

* ``.mcool`` / ``.cool``: metadata-validated via ``cooler``; single-resolution
  ``.cool`` files can optionally be zoomified to ``.mcool`` when coarser
  compartment resolutions are requested.
* ``.pairs`` / ``.pairs.gz`` (4DN pairs): requires explicit chromosome sizes
  (``chrom_sizes`` dict/TSV or ``genome_fasta`` ``.fai`` — never guessed from
  assembly names); binned via ``cooler cload pairs`` + optional ``cooler zoomify``.
* ``.bam`` / ``.sam`` (chimeric ``bwa-mem2 -SP5M`` alignments): parsed, sorted,
  and deduplicated via ``pairtools parse -> sort -> dedup`` into ``.pairs.gz``
  with measured pairtools summary statistics, then binned via ``cooler cload pairs``.
* ``.hic`` (Juicer): streamed chromosome-pair by chromosome-pair via ``hicstraw``
  into an out-of-core ``.cool`` (+ optional ``.mcool``) without loading the
  genome-wide dense matrix into RAM.

Zero imports from ``aria.agents.*``, ``litellm``, or ``pydantic``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from aria.scripts._base import run_script
from aria.utils.hic_constants import list_mcool_resolutions
from aria.utils.provenance import hash_file
from aria.utils.stage_manifest import (
    build_stage_manifest,
    stage_is_current,
    write_stage_manifest,
)

STAGE_VERSION = 1


def _check_hic_stage_cache(
    out_dir: Path,
    *,
    stage: str,
    version: int,
    params: dict[str, Any],
    inputs: list[str],
    manifest_name: str = "stage_manifest.json",
) -> dict[str, Any] | None:
    manifest_path = Path(out_dir) / manifest_name
    input_tuples = [(f"in_{i}", str(p)) for i, p in enumerate(inputs)]
    ok, _ = stage_is_current(
        manifest_path,
        inputs=input_tuples,
        params=params,
        tool_version=f"{stage}:v{version}",
    )
    if not ok:
        return None
    try:
        import json as _json

        payload = _json.loads(manifest_path.read_text(encoding="utf-8"))
        extra = payload.get("extra") or {}
        outputs = extra.get("outputs") or []
        if any(not Path(str(o)).is_file() for o in outputs if o):
            return None
        output_sha256 = extra.get("output_sha256")
        if outputs and not isinstance(output_sha256, dict):
            return None
        for o in outputs:
            if not o:
                continue
            expected_sha = (output_sha256 or {}).get(str(o))
            if not expected_sha or hash_file(o) != expected_sha:
                return None
        return extra
    except Exception:
        return None


def _save_hic_stage_cache(
    out_dir: Path,
    *,
    stage: str,
    version: int,
    params: dict[str, Any],
    inputs: list[str],
    outputs: list[str],
    extra: dict[str, Any],
    manifest_name: str = "stage_manifest.json",
) -> None:
    manifest_path = Path(out_dir) / manifest_name
    input_tuples = [(f"in_{i}", str(p)) for i, p in enumerate(inputs)]
    clean_outputs = [str(o) for o in outputs if o]
    output_sha256 = {
        o: hash_file(o) for o in clean_outputs if Path(o).is_file()
    }
    manifest = build_stage_manifest(
        stage=stage,
        inputs=input_tuples,
        params=params,
        tool_version=f"{stage}:v{version}",
        extra={
            **extra,
            "outputs": clean_outputs,
            "output_sha256": output_sha256,
        },
    )
    write_stage_manifest(manifest_path, manifest)
DEFAULT_BASE_RESOLUTION = 10_000
DEFAULT_ZOOM_RESOLUTIONS = (10_000, 40_000, 100_000, 500_000, 1_000_000)


def _strip_hic_suffixes(path_str: str) -> str:
    name = Path(path_str).name
    for suf in (
        ".pairs.gz", ".pairs", ".mcool", ".cool", ".hic",
        ".bam", ".sam", ".cram",
    ):
        if name.lower().endswith(suf):
            return name[: -len(suf)] or "sample"
    return Path(path_str).stem or "sample"


def _resolve_chromsizes_file(params: dict, out_dir: Path) -> Path | None:
    """Materialize a 2-column ``chrom<TAB>length`` file from explicit inputs."""
    raw = params.get("chrom_sizes") or params.get("chromsizes")
    dest = out_dir / "chrom.sizes"
    if isinstance(raw, dict) and raw:
        lines = [
            f"{str(k)}\t{int(v)}"
            for k, v in raw.items()
            if str(k).strip() and int(v) > 0
        ]
        if lines:
            dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
            return dest
    if isinstance(raw, (str, Path)) and Path(str(raw)).is_file():
        src = Path(str(raw))
        rows: list[str] = []
        for line in src.read_text(encoding="utf-8").splitlines():
            parts = line.strip().split()
            if len(parts) >= 2 and not parts[0].startswith("#"):
                try:
                    length = int(parts[1])
                except ValueError:
                    continue
                if length > 0:
                    rows.append(f"{parts[0]}\t{length}")
        if rows:
            dest.write_text("\n".join(rows) + "\n", encoding="utf-8")
            return dest

    for key in ("genome_fai", "genome_fasta", "reference_fasta"):
        cand = params.get(key)
        if not cand:
            continue
        p = Path(str(cand))
        fai = p if p.suffix == ".fai" else Path(f"{p}.fai")
        if fai.is_file():
            rows = []
            for line in fai.read_text(encoding="utf-8").splitlines():
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    try:
                        length = int(parts[1])
                    except ValueError:
                        continue
                    if length > 0:
                        rows.append(f"{parts[0]}\t{length}")
            if rows:
                dest.write_text("\n".join(rows) + "\n", encoding="utf-8")
                return dest
        if p.is_file() and p.suffix in {".fa", ".fasta", ".fna"}:
            import pysam

            with pysam.FastaFile(str(p)) as fa:
                rows = [
                    f"{ref}\t{int(fa.get_reference_length(ref))}"
                    for ref in fa.references
                    if int(fa.get_reference_length(ref)) > 0
                ]
            if rows:
                dest.write_text("\n".join(rows) + "\n", encoding="utf-8")
                return dest
    return None


def _stack_bin(name: str) -> str | None:
    sib = Path(sys.executable).resolve().parent / name
    if sib.is_file() and os.access(sib, os.X_OK):
        return str(sib)
    return shutil.which(name)


def _zoomify_cooler(
    cool_path: Path,
    mcool_path: Path,
    base_res: int,
    zoom_resolutions: list[int],
    threads: int,
) -> Path:
    valid_zooms = sorted({int(base_res), *[int(r) for r in zoom_resolutions if int(r) >= int(base_res) and int(r) % int(base_res) == 0]})
    if len(valid_zooms) <= 1:
        return cool_path
    res_arg = ",".join(str(r) for r in valid_zooms)
    cooler_exe = _stack_bin("cooler") or "cooler"
    cmd = [
        cooler_exe, "zoomify",
        "-p", str(max(1, int(threads))),
        "-r", res_arg,
        "-o", str(mcool_path),
        str(cool_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    if proc.returncode != 0 or not mcool_path.is_file():
        raise RuntimeError(
            f"cooler zoomify failed (rc={proc.returncode}): {proc.stderr[-400:]}"
        )
    return mcool_path


def _cload_pairs(
    pairs_path: Path,
    chromsizes_path: Path,
    cool_path: Path,
    base_res: int,
    assembly: str,
) -> None:
    cooler_exe = _stack_bin("cooler") or "cooler"
    cmd = [
        cooler_exe, "cload", "pairs",
        "--assembly", str(assembly),
        "-c1", "2", "-p1", "3", "-c2", "4", "-p2", "5",
        f"{chromsizes_path}:{int(base_res)}",
        str(pairs_path),
        str(cool_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    if proc.returncode != 0 or not cool_path.is_file():
        raise RuntimeError(
            f"cooler cload pairs failed (rc={proc.returncode}): {proc.stderr[-400:]}"
        )


def _parse_pairtools_stats(stats_path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    if not stats_path.is_file():
        return out
    for line in stats_path.read_text(encoding="utf-8").splitlines():
        parts = line.strip().split("\t")
        if len(parts) == 2 and "/" not in parts[0]:
            try:
                out[parts[0]] = int(float(parts[1]))
            except ValueError:
                continue
    return out


def _convert_bam_to_pairs(
    bam_path: Path,
    chromsizes_path: Path,
    sample_dir: Path,
    sample_id: str,
    mapq: int,
    threads: int,
) -> tuple[Path, dict[str, Any]]:
    pairtools_exe = _stack_bin("pairtools")
    if pairtools_exe is None:
        raise RuntimeError("pairtools binary not found in aria-hic-env.")
    parsed_pairs = sample_dir / f"{sample_id}.parsed.pairs"
    sorted_pairs = sample_dir / f"{sample_id}.sorted.pairs"
    dedup_pairs = sample_dir / f"{sample_id}.dedup.pairs.gz"
    stats_path = sample_dir / f"{sample_id}.pairtools.stats.tsv"
    tmpdir = sample_dir / "sort_tmp"
    tmpdir.mkdir(parents=True, exist_ok=True)

    parse_cmd = [
        pairtools_exe, "parse",
        "--min-mapq", str(int(mapq)),
        "--walks-policy", "5unique",
        "--max-inter-align-gap", "30",
        "--chroms-path", str(chromsizes_path),
        "-o", str(parsed_pairs),
        str(bam_path),
    ]
    p1 = subprocess.run(parse_cmd, capture_output=True, text=True, timeout=7200)
    if p1.returncode != 0:
        raise RuntimeError(f"pairtools parse failed: {p1.stderr[-400:]}")

    sort_cmd = [
        pairtools_exe, "sort",
        "--nproc", str(max(1, int(threads))),
        "--tmpdir", str(tmpdir),
        "-o", str(sorted_pairs),
        str(parsed_pairs),
    ]
    p2 = subprocess.run(sort_cmd, capture_output=True, text=True, timeout=7200)
    if p2.returncode != 0:
        raise RuntimeError(f"pairtools sort failed: {p2.stderr[-400:]}")

    stats_path.unlink(missing_ok=True)
    dedup_cmd = [
        pairtools_exe, "dedup",
        "--mark-dups",
        "--output-stats", str(stats_path),
        "-o", str(dedup_pairs),
        str(sorted_pairs),
    ]
    p3 = subprocess.run(dedup_cmd, capture_output=True, text=True, timeout=7200)
    if p3.returncode != 0:
        raise RuntimeError(f"pairtools dedup failed: {p3.stderr[-400:]}")

    for tmp_f in (parsed_pairs, sorted_pairs):
        tmp_f.unlink(missing_ok=True)
    shutil.rmtree(tmpdir, ignore_errors=True)
    return dedup_pairs, _parse_pairtools_stats(stats_path)


def _convert_hic_to_cooler(
    hic_path: Path,
    cool_path: Path,
    target_res: int,
    assembly: str,
) -> tuple[int, list[int]]:
    import bioframe
    import cooler
    import hicstraw
    import numpy as np
    import pandas as pd

    hic = hicstraw.HiCFile(str(hic_path))
    avail = sorted(int(r) for r in hic.getResolutions())
    if not avail:
        raise RuntimeError(f"No resolutions reported inside .hic file: {hic_path}")
    chosen_res = int(target_res) if int(target_res) in avail else min(
        avail, key=lambda r: abs(r - int(target_res))
    )
    chrom_objs = [
        c for c in hic.getChromosomes()
        if str(c.name) not in {"All", "ALL", "all"} and int(c.length) > 0
    ]
    if not chrom_objs:
        raise RuntimeError(f"No valid chromosomes in .hic file: {hic_path}")

    chrom_series = pd.Series(
        {str(c.name): int(c.length) for c in chrom_objs}
    )
    bins_df = bioframe.binnify(chrom_series, chosen_res)
    chrom_offsets: dict[str, int] = {}
    chrom_nbins: dict[str, int] = {}
    for chrom_name, grp in bins_df.groupby("chrom", sort=False, observed=False):
        chrom_offsets[str(chrom_name)] = int(grp.index[0])
        chrom_nbins[str(chrom_name)] = int(len(grp))

    def _pixel_chunks():
        import math

        names = [str(c.name) for c in chrom_objs]
        lengths = {str(c.name): int(c.length) for c in chrom_objs}
        for i, c1 in enumerate(names):
            off1 = chrom_offsets[c1]
            nb1 = chrom_nbins[c1]
            len1 = lengths[c1]
            for j in range(i, len(names)):
                c2 = names[j]
                off2 = chrom_offsets[c2]
                nb2 = chrom_nbins[c2]
                len2 = lengths[c2]
                mzd = hic.getMatrixZoomData(
                    c1, c2, "observed", "NONE", "BP", chosen_res
                )
                if mzd is None:
                    continue
                recs = mzd.getRecords(0, max(0, len1 - 1), 0, max(0, len2 - 1))
                if not recs:
                    continue
                b1_list: list[int] = []
                b2_list: list[int] = []
                cnt_list: list[int] = []
                for r in recs:
                    try:
                        c_raw = float(r.counts)
                    except (TypeError, ValueError):
                        continue
                    if not math.isfinite(c_raw) or c_raw <= 0.0:
                        continue
                    c_int = int(round(c_raw))
                    if c_int <= 0:
                        continue
                    lb1 = int(r.binX) // chosen_res
                    lb2 = int(r.binY) // chosen_res
                    if not (0 <= lb1 < nb1 and 0 <= lb2 < nb2):
                        continue
                    g1 = off1 + lb1
                    g2 = off2 + lb2
                    if g1 > g2:
                        g1, g2 = g2, g1
                    b1_list.append(g1)
                    b2_list.append(g2)
                    cnt_list.append(c_int)
                if not cnt_list:
                    continue
                df_chunk = pd.DataFrame({
                    "bin1_id": np.asarray(b1_list, dtype=np.int64),
                    "bin2_id": np.asarray(b2_list, dtype=np.int64),
                    "count": np.asarray(cnt_list, dtype=np.int64),
                })
                df_chunk = (
                    df_chunk.groupby(["bin1_id", "bin2_id"], as_index=False)["count"]
                    .sum()
                    .sort_values(["bin1_id", "bin2_id"], kind="mergesort")
                )
                yield df_chunk

    cooler.create_cooler(
        str(cool_path),
        bins_df,
        _pixel_chunks(),
        assembly=str(assembly),
        ordered=False,
        symmetric_upper=True,
    )
    return chosen_res, avail


def hic_preprocess(params: dict) -> dict:
    import cooler

    files = [str(f) for f in (params.get("files") or []) if f]
    if not files:
        return {
            "status": "error",
            "error_type": "NoInputFiles",
            "details": "hic_preprocess requires non-empty 'files'.",
        }
    missing = [f for f in files if not Path(f).exists()]
    if missing:
        return {
            "status": "error",
            "error_type": "FileNotFound",
            "details": f"Hi-C input files not found: {missing}",
        }

    genome = str(params.get("genome") or "hg38")
    base_res = int(params.get("resolution") or DEFAULT_BASE_RESOLUTION)
    zoom_res = [
        int(r) for r in (params.get("zoom_resolutions") or DEFAULT_ZOOM_RESOLUTIONS)
    ]
    build_mcool = bool(params.get("zoomify", True))
    threads = max(1, int(params.get("threads", 4)))
    mapq = int(params.get("mapq", 30))
    out_dir = Path(
        str(params.get("output_dir") or (Path(files[0]).parent / "aria_hic_preprocessed"))
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    declared_sids = params.get("sample_ids")
    if isinstance(declared_sids, list) and len(declared_sids) == len(files):
        raw_sids = [str(s) for s in declared_sids]
    else:
        raw_sids = [_strip_hic_suffixes(f) for f in files]
    # Deduplicate colliding derived sample_ids (e.g. rep1/matrix.mcool, rep2/matrix.mcool).
    sid_counts: dict[str, int] = {}
    for s in raw_sids:
        sid_counts[s] = sid_counts.get(s, 0) + 1
    seen_sids: set[str] = set()
    sample_ids: list[str] = []
    for idx_s, (s, fpath_s) in enumerate(zip(raw_sids, files, strict=True)):
        cand = s
        if sid_counts.get(s, 0) > 1:
            parent_tag = Path(fpath_s).parent.name.strip()
            cand = f"{parent_tag}_{s}" if parent_tag else f"{s}_{idx_s + 1}"
            if cand in seen_sids:
                cand = f"{s}_{idx_s + 1}"
        while cand in seen_sids:
            cand = f"{cand}_{idx_s + 1}"
        seen_sids.add(cand)
        sample_ids.append(cand)

    chrom_sizes_param = params.get("chrom_sizes") or params.get("chromsizes")
    genome_fasta_param = params.get("genome_fasta") or params.get("reference_fasta")
    genome_fai_param = params.get("genome_fai")
    cache_inputs = list(files)
    for rv in (chrom_sizes_param, genome_fasta_param, genome_fai_param):
        if isinstance(rv, (str, Path)) and str(rv) and Path(str(rv)).is_file():
            cache_inputs.append(str(rv))
    manifest_params = {
        "genome": genome,
        "base_res": base_res,
        "zoom_res": sorted(zoom_res),
        "zoomify": build_mcool,
        "mapq": mapq,
        "sample_ids": sample_ids,
        "chrom_sizes": (
            chrom_sizes_param
            if isinstance(chrom_sizes_param, dict)
            else str(chrom_sizes_param or "")
        ),
        "genome_fasta": str(genome_fasta_param or ""),
        "genome_fai": str(genome_fai_param or ""),
    }
    if bool(params.get("resume", True)):
        cached = _check_hic_stage_cache(
            out_dir,
            stage="hic_preprocess",
            version=STAGE_VERSION,
            params=manifest_params,
            inputs=cache_inputs,
        )
        if cached is not None and isinstance(cached.get("result"), dict):
            resumed = dict(cached["result"])
            resumed["resumed"] = True
            return resumed

    needs_chromsizes = any(
        str(f).lower().endswith((".pairs", ".pairs.gz", ".bam", ".sam"))
        for f in files
    )
    chromsizes_path = _resolve_chromsizes_file(params, out_dir) if needs_chromsizes else None
    if needs_chromsizes and chromsizes_path is None:
        return {
            "status": "error",
            "error_type": "MissingChromSizes",
            "details": (
                "Converting .pairs or .bam Hi-C inputs to .cool requires explicit "
                "chromosome sizes (chrom_sizes dict/TSV or genome_fasta/.fai)."
            ),
        }

    processed_samples: list[dict[str, Any]] = []
    cooler_files: list[str] = []

    for sid, src_str in zip(sample_ids, files, strict=True):
        src = Path(src_str)
        low = src.name.lower()
        sdir = out_dir / sid
        sdir.mkdir(parents=True, exist_ok=True)

        try:
            if low.endswith(".mcool"):
                res_vals = sorted(list_mcool_resolutions(str(src)).keys())
                cooler_files.append(str(src))
                processed_samples.append({
                    "sample_id": sid,
                    "input_path": str(src),
                    "input_format": "mcool",
                    "cooler_path": str(src),
                    "resolutions": res_vals,
                    "status": "success",
                })
            elif low.endswith(".cool"):
                clr = cooler.Cooler(str(src))
                src_res = int(clr.binsize)
                out_path = src
                res_vals = [src_res]
                if build_mcool and any(r > src_res and r % src_res == 0 for r in zoom_res):
                    mcool_out = sdir / f"{sid}.mcool"
                    out_path = _zoomify_cooler(src, mcool_out, src_res, zoom_res, threads)
                    res_vals = sorted(list_mcool_resolutions(str(out_path)).keys()) or [src_res]
                cooler_files.append(str(out_path))
                processed_samples.append({
                    "sample_id": sid,
                    "input_path": str(src),
                    "input_format": "cool",
                    "cooler_path": str(out_path),
                    "resolutions": res_vals,
                    "status": "success",
                })
            elif low.endswith((".pairs", ".pairs.gz")):
                assert chromsizes_path is not None
                cool_out = sdir / f"{sid}.{base_res}.cool"
                _cload_pairs(src, chromsizes_path, cool_out, base_res, genome)
                out_path = cool_out
                res_vals = [base_res]
                if build_mcool:
                    mcool_out = sdir / f"{sid}.mcool"
                    out_path = _zoomify_cooler(cool_out, mcool_out, base_res, zoom_res, threads)
                    res_vals = sorted(list_mcool_resolutions(str(out_path)).keys()) or [base_res]
                cooler_files.append(str(out_path))
                processed_samples.append({
                    "sample_id": sid,
                    "input_path": str(src),
                    "input_format": "pairs",
                    "cooler_path": str(out_path),
                    "resolutions": res_vals,
                    "status": "success",
                })
            elif low.endswith((".bam", ".sam")):
                assert chromsizes_path is not None
                pairs_out, pt_stats = _convert_bam_to_pairs(
                    src, chromsizes_path, sdir, sid, mapq=mapq, threads=threads
                )
                cool_out = sdir / f"{sid}.{base_res}.cool"
                _cload_pairs(pairs_out, chromsizes_path, cool_out, base_res, genome)
                out_path = cool_out
                res_vals = [base_res]
                if build_mcool:
                    mcool_out = sdir / f"{sid}.mcool"
                    out_path = _zoomify_cooler(cool_out, mcool_out, base_res, zoom_res, threads)
                    res_vals = sorted(list_mcool_resolutions(str(out_path)).keys()) or [base_res]
                cooler_files.append(str(out_path))
                processed_samples.append({
                    "sample_id": sid,
                    "input_path": str(src),
                    "input_format": "bam",
                    "pairs_path": str(pairs_out),
                    "pairtools_stats": pt_stats,
                    "cooler_path": str(out_path),
                    "resolutions": res_vals,
                    "status": "success",
                })
            elif low.endswith(".hic"):
                cool_out = sdir / f"{sid}.cool"
                chosen_res, avail_hic_res = _convert_hic_to_cooler(
                    src, cool_out, base_res, genome
                )
                out_path = cool_out
                res_vals = [chosen_res]
                if build_mcool and any(r > chosen_res and r % chosen_res == 0 for r in zoom_res):
                    mcool_out = sdir / f"{sid}.mcool"
                    out_path = _zoomify_cooler(
                        cool_out, mcool_out, chosen_res, zoom_res, threads
                    )
                    res_vals = sorted(list_mcool_resolutions(str(out_path)).keys()) or [chosen_res]
                cooler_files.append(str(out_path))
                processed_samples.append({
                    "sample_id": sid,
                    "input_path": str(src),
                    "input_format": "hic",
                    "hic_available_resolutions": avail_hic_res,
                    "cooler_path": str(out_path),
                    "resolutions": res_vals,
                    "status": "success",
                })
            else:
                return {
                    "status": "error",
                    "error_type": "UnsupportedHiCFormat",
                    "details": f"Unsupported Hi-C file format: {src}",
                }
        except Exception as exc:
            return {
                "status": "error",
                "error_type": "HiCPreprocessFailed",
                "sample_id": sid,
                "details": str(exc)[:500],
            }

    result = {
        "status": "success",
        "ran": True,
        "analysis": "hic_preprocess",
        "data_type": "HiC",
        "validation_level": "production",
        "genome": genome,
        "base_resolution": int(base_res),
        "n_samples": len(cooler_files),
        "sample_ids": sample_ids,
        "cooler_files": cooler_files,
        "samples": processed_samples,
        "output_dir": str(out_dir),
    }
    _save_hic_stage_cache(
        out_dir,
        stage="hic_preprocess",
        version=STAGE_VERSION,
        params=manifest_params,
        inputs=cache_inputs,
        outputs=cooler_files,
        extra={"result": result},
    )
    return result


if __name__ == "__main__":
    run_script(hic_preprocess)
