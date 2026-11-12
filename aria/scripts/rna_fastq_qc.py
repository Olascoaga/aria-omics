"""
ARIA RNA-seq FASTQ QC Script
------------------------------
Runs fastp (trimming + QC) on all samples, then MultiQC for summary.
Executed inside aria-rnaseq-env via EnvironmentManager.

Detects paired-end vs single-end automatically from file naming.
Naming patterns supported:
  sample_1.fq.gz / sample_2.fq.gz
  sample_R1.fastq.gz / sample_R2.fastq.gz
  sample_R1_001.fastq.gz / sample_R2_001.fastq.gz

Input params:
  fastq_dir:      str   — directory containing FASTQ files
  fastq_files:    list  — optional explicit selection within/alongside fastq_dir
  output_dir:     str   — where to write trimmed FASTQs + QC reports
  threads:        int   — threads per sample (default: 4)
  min_length:     int   — minimum read length after trimming (default: 36)
  quality:        int   — phred quality threshold (default: 20)

Output:
  {
    "status":        "success",
    "samples":       [{"name", "read_layout", "lanes", "r1_trimmed_files",
                       "r2_trimmed_files", "paired", "n_reads_raw",
                       "n_reads_trimmed", "pct_passed"}],
    "multiqc_report": str,   — path to MultiQC HTML
    "output_dir":    str,
    "warnings":      [str]
  }
"""

from __future__ import annotations
import sys, os, subprocess, json, re
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import mocks_allowed, run_script
from aria.utils.stage_manifest import (
    build_stage_manifest, write_stage_manifest, stage_is_current,
)


def _fastp_version() -> str:
    """Best-effort fastp version string for the A4 stage manifest (empty when
    the tool is unavailable)."""
    try:
        proc = subprocess.run(["fastp", "--version"],
                              capture_output=True, text=True, timeout=30)
        return (proc.stderr or proc.stdout or "").strip()
    except Exception:
        return ""


def rna_fastq_qc(params: dict) -> dict:
    fastq_dir  = Path(params["fastq_dir"])
    output_dir = Path(params.get("output_dir",
                                  str(fastq_dir.parent / "aria_qc")))
    threads    = int(params.get("threads", 4))
    library_protocol = str(params.get("library_protocol") or "full_length").strip().lower()
    default_min_len = 20 if library_protocol in {"tag_3prime", "tag_5prime", "umi_bulk"} else 36
    min_len    = int(params.get("min_length", default_min_len))
    quality    = int(params.get("quality", 20))
    umi_config = params.get("umi_config") if isinstance(params.get("umi_config"), dict) else None
    trim_poly_x = bool(
        params.get(
            "trim_poly_x",
            library_protocol in {"tag_3prime", "tag_5prime"}
            or bool((umi_config or {}).get("trim_poly_x", False)),
        )
    )
    trim_poly_g = bool(params.get("trim_poly_g", True))
    trim_front1 = int(params.get("trim_front1", 0) or 0)
    r1_only = bool(params.get("r1_only", False))
    allow_mock = mocks_allowed(params)
    warnings   = []

    trimmed_dir = output_dir / "trimmed"
    qc_dir      = output_dir / "fastqc"
    fastp_dir   = output_dir / "fastp_reports"
    for d in [trimmed_dir, qc_dir, fastp_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # ── Detect samples ────────────────────────────────────────────────────
    samples = _detect_samples(
        fastq_dir,
        warnings,
        fastq_files=params.get("fastq_files"),
        r1_only=r1_only,
        known_samples=params.get("known_samples"),
    )
    if not samples:
        return {
            "status":     "error",
            "error_type": ("FastqManifestInvalid" if warnings
                           else "NoFASTQFound"),
            "details":    ("; ".join(warnings) if warnings else
                           f"No FASTQ files found in {fastq_dir}"),
            "warnings":   warnings,
        }

    # ── Run fastp per sample ──────────────────────────────────────────────
    processed = []
    for sample in samples:
        result = _run_fastp(
            sample=sample,
            trimmed_dir=trimmed_dir,
            fastp_dir=fastp_dir,
            threads=threads,
            min_len=min_len,
            quality=quality,
            warnings=warnings,
            allow_mock=allow_mock,
            library_protocol=library_protocol,
            trim_poly_x=trim_poly_x,
            trim_poly_g=trim_poly_g,
            trim_front1=trim_front1,
            umi_config=umi_config,
        )
        processed.append(result)

    # ── MultiQC summary ───────────────────────────────────────────────────
    multiqc_report = _run_multiqc(fastp_dir, output_dir, warnings)

    failed_samples = [
        sample.get("name", "unknown")
        for sample in processed
        if sample.get("status") != "success"
    ]
    if failed_samples:
        return {
            "status":         "error",
            "error_type":     "PartialFastqQCFailure",
            "details":        (
                "FASTQ QC is all-or-fail; one or more samples or lanes did "
                "not complete successfully."
            ),
            "failed_samples": failed_samples,
            "samples":        processed,
            "n_samples":      len(processed),
            "library_protocol": library_protocol,
            "multiqc_report": str(multiqc_report) if multiqc_report else None,
            "output_dir":     str(output_dir),
            "trimmed_dir":    str(trimmed_dir),
            "warnings":       warnings,
        }

    return {
        "status":         "success",
        "samples":        processed,
        "n_samples":      len(processed),
        "library_protocol": library_protocol,
        "umi_extracted":  bool(umi_config),
        "trim_poly_x":    trim_poly_x,
        "multiqc_report": str(multiqc_report) if multiqc_report else None,
        "output_dir":     str(output_dir),
        "trimmed_dir":    str(trimmed_dir),
        "warnings":       warnings,
    }


def _strip_fastq_ext(filename: str) -> str:
    stem = filename
    for extension in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if stem.lower().endswith(extension):
            return stem[:-len(extension)]
    return stem


def _detect_samples(fastq_dir: Path, warnings: list,
                    fastq_files: list[str] | None = None,
                    *,
                    r1_only: bool = False,
                    known_samples: list[str] | None = None) -> list[dict]:
    """
    Build an explicit sample/lane/read-layout manifest.

    Illumina lane names such as ``sample_S1_L001_R1_001.fastq.gz`` are
    grouped under one biological sample instead of becoming separate samples.
    Single-end files without an ``_R1`` suffix (e.g. ``SRR123456.fastq.gz`` or
    ``ctrl_1.fastq.gz``, ``ctrl_2.fastq.gz``, ``ctrl_3.fastq.gz``) are
    supported as single-end ``R1`` lanes, while SRA-style paired mates
    (``WT1_1.fq.gz`` + ``WT1_2.fq.gz`` with numeric suffixes strictly equal to
    ``{1, 2}`` and not declared as standalone sample IDs) are paired, and
    barcode/index slots (``_R3``, ``_I1``, ``_I2``) remain rejected.
    """
    # Find all FASTQ files
    if fastq_files is None:
        candidates = list(fastq_dir.iterdir())
    else:
        candidates = [Path(path) for path in fastq_files]
    fq_files = sorted(
        f for f in candidates
        if f.name.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))
    )

    if not fq_files:
        return []

    bare_numeric_suffixes: dict[str, set[int]] = {}
    for path in fq_files:
        stem = _strip_fastq_ext(path.name)
        m_num = re.match(r"^(?P<prefix>.+?)_(?P<num>\d+)$", stem)
        if m_num:
            bare_numeric_suffixes.setdefault(
                m_num.group("prefix"), set()
            ).add(int(m_num.group("num")))

    known_set = {str(s) for s in (known_samples or ()) if s}

    grouped: dict[str, dict[str, dict[str, Path]]] = {}
    duplicate_slots: list[str] = []
    unmatched: list[str] = []
    for path in fq_files:
        stem = _strip_fastq_ext(path.name)
        m_num = re.match(r"^(?P<prefix>.+?)_(?P<num>\d+)$", stem)
        allow_bare_12_pair = bool(
            m_num
            and bare_numeric_suffixes.get(m_num.group("prefix")) == {1, 2}
            and stem not in known_set
        )
        parsed = _parse_fastq_name(
            path.name, allow_bare_12_pair=allow_bare_12_pair
        )
        if parsed is None:
            unmatched.append(path.name)
            continue
        sample_name, lane, read = parsed
        if r1_only and read == "R2":
            continue
        slot = grouped.setdefault(sample_name, {}).setdefault(lane, {})
        if read in slot:
            duplicate_slots.append(f"{sample_name}/{lane}/{read}")
        else:
            slot[read] = path

    if unmatched:
        warnings.append(
            "Unrecognized FASTQ names (unexpected barcode/index read slot): "
            + ", ".join(unmatched)
        )
    if duplicate_slots:
        warnings.append(
            "Duplicate FASTQ read slots detected: " + ", ".join(duplicate_slots)
        )
        return []

    manifest_invalid = bool(unmatched)

    samples = []
    for name in sorted(grouped):
        lane_records = []
        invalid = False
        for lane, reads in sorted(grouped[name].items()):
            if "R1" not in reads:
                warnings.append(f"Sample '{name}' lane {lane}: R2 has no R1.")
                invalid = True
                manifest_invalid = True
                continue
            lane_records.append({
                "lane": lane,
                "r1": str(reads["R1"]),
                "r2": None if r1_only else (str(reads["R2"]) if "R2" in reads else None),
            })
        if invalid or not lane_records:
            continue

        paired_lanes = [lane["r2"] is not None for lane in lane_records]
        if any(paired_lanes) and not all(paired_lanes):
            warnings.append(
                f"Sample '{name}' mixes paired and single-end lanes; refusing "
                "an ambiguous read layout."
            )
            manifest_invalid = True
            continue

        paired = all(paired_lanes)
        samples.append({
            "name": name,
            "read_layout": "paired-end" if paired else "single-end",
            "paired": paired,
            "lanes": lane_records,
            # Backward-compatible aliases for single-lane consumers.
            "r1": lane_records[0]["r1"],
            "r2": lane_records[0]["r2"],
        })

    if manifest_invalid:
        return []

    if not samples:
        warnings.append(
            "Could not resolve FASTQ sample manifest. "
            "Expected paired R1/R2 or single-end FASTQ files."
        )

    return samples


def _parse_fastq_name(
    filename: str, *, allow_bare_12_pair: bool = True
) -> tuple[str, str, str] | None:
    """Return ``(sample, lane, R1|R2)`` for supported paired or single-end FASTQ names."""
    stem = _strip_fastq_ext(filename)

    # Explicitly reject non-bulk barcode/index slots (R3, R4, I1, I2).
    if re.search(r"_(?:R[34]|I[12])(?:_\d+)?$", stem, re.IGNORECASE):
        return None

    patterns: list[str] = [
        # Standard bcl2fastq / DRAGEN naming. The optional sample-number token
        # is technical metadata and is not part of the biological sample ID.
        r"^(?P<sample>.+?)(?:_S\d+)?_(?P<lane>L\d{3})_"
        r"(?P<read>R[12])(?:_\d+)?$",
        r"^(?P<sample>.+?)_(?P<read>R[12])(?:_\d+)?$",
    ]
    if allow_bare_12_pair:
        patterns.append(r"^(?P<sample>.+?)_(?P<read>[12])$")
    for pattern in patterns:
        match = re.match(pattern, stem, flags=re.IGNORECASE)
        if not match:
            continue
        read = match.group("read").upper()
        if read in {"1", "2"}:
            read = f"R{read}"
        lane = (match.groupdict().get("lane") or "L001").upper()
        return match.group("sample"), lane, read

    # Single-end lane naming without _R1 (e.g. sample_S1_L002_001.fastq.gz)
    se_lane = re.match(
        r"^(?P<sample>.+?)(?:_S\d+)?_(?P<lane>L\d{3})(?:_\d+)?$",
        stem,
        flags=re.IGNORECASE,
    )
    if se_lane:
        return se_lane.group("sample"), se_lane.group("lane").upper(), "R1"

    # Generic single-end FASTQ (e.g. SRR123456.fastq.gz, ctrl_rep1.fq.gz)
    if stem:
        return stem, "L001", "R1"
    return None


def _fastp_outputs_valid(r1_out: str, r2_out: str, json_out: str,
                          paired: bool) -> bool:
    """
    Check if fastp outputs are present and valid for skipping.
    Requires:
      - R1 trimmed file > 64 bytes (non-empty gzip payload)
      - R2 trimmed file > 64 bytes if paired
      - JSON report > 100 bytes (non-empty report)
    """
    try:
        r1 = Path(r1_out)
        if not r1.exists() or r1.stat().st_size < 64:
            return False
        if paired:
            if not r2_out:
                return False
            r2 = Path(r2_out)
            if not r2.exists() or r2.stat().st_size < 64:
                return False
        j = Path(json_out)
        if not j.exists() or j.stat().st_size < 100:
            return False
        return True
    except Exception:
        return False


def _run_fastp(sample: dict, trimmed_dir: Path, fastp_dir: Path,
               threads: int, min_len: int, quality: int,
               warnings: list,
               allow_mock: bool = False,
               *,
               library_protocol: str = "full_length",
               trim_poly_x: bool = False,
               trim_poly_g: bool = True,
               trim_front1: int = 0,
               umi_config: dict | None = None) -> dict:
    """Run fastp for every lane and aggregate one sample-level result."""
    name = sample["name"]
    read_layout = sample.get(
        "read_layout",
        "paired-end" if sample.get("paired") else "single-end",
    )
    lanes = sample.get("lanes") or [{
        "lane": "L001",
        "r1": sample["r1"],
        "r2": sample.get("r2"),
    }]

    lane_results = []
    multiple_lanes = len(lanes) > 1
    for lane in lanes:
        lane_name = f"{name}_{lane['lane']}" if multiple_lanes else name
        lane_result = _run_fastp_lane(
            sample={
                "name": lane_name,
                "r1": lane["r1"],
                "r2": lane.get("r2"),
                "paired": read_layout == "paired-end",
            },
            trimmed_dir=trimmed_dir,
            fastp_dir=fastp_dir,
            threads=threads,
            min_len=min_len,
            quality=quality,
            warnings=warnings,
            allow_mock=allow_mock,
            library_protocol=library_protocol,
            trim_poly_x=trim_poly_x,
            trim_poly_g=trim_poly_g,
            trim_front1=trim_front1,
            umi_config=umi_config,
        )
        lane_result["lane"] = lane["lane"]
        lane_results.append(lane_result)

    failed_lanes = [
        lane["lane"] for lane in lane_results
        if lane.get("status") != "success"
    ]
    r1_trimmed_files = [lane["r1_trimmed"] for lane in lane_results]
    r2_trimmed_files = [
        lane["r2_trimmed"] for lane in lane_results
        if lane.get("r2_trimmed")
    ]
    n_raw = sum(int(lane.get("n_reads_raw", 0)) for lane in lane_results)
    n_trimmed = sum(
        int(lane.get("n_reads_trimmed", 0)) for lane in lane_results
    )
    q30_numerator = sum(
        float(lane.get("q30_rate", 0)) * int(lane.get("n_reads_trimmed", 0))
        for lane in lane_results
    )

    return {
        "name": name,
        "status": "failed" if failed_lanes else "success",
        "read_layout": read_layout,
        "paired": read_layout == "paired-end",
        "library_protocol": library_protocol,
        "umi_extracted": bool(umi_config),
        "lanes": lane_results,
        "failed_lanes": failed_lanes,
        "r1_raw": lanes[0]["r1"],
        "r2_raw": lanes[0].get("r2"),
        "r1_trimmed": ",".join(r1_trimmed_files),
        "r2_trimmed": (",".join(r2_trimmed_files)
                       if r2_trimmed_files else None),
        "r1_trimmed_files": r1_trimmed_files,
        "r2_trimmed_files": r2_trimmed_files,
        "n_reads_raw": n_raw,
        "n_reads_trimmed": n_trimmed,
        "pct_passed": round(n_trimmed / max(n_raw, 1) * 100, 1),
        "q30_rate": round(q30_numerator / max(n_trimmed, 1), 1),
        "resumed": bool(lane_results) and all(
            lane.get("resumed") is True for lane in lane_results
        ),
    }


def _run_fastp_lane(sample: dict, trimmed_dir: Path, fastp_dir: Path,
                    threads: int, min_len: int, quality: int,
                    warnings: list,
                    allow_mock: bool = False,
                    *,
                    library_protocol: str = "full_length",
                    trim_poly_x: bool = False,
                    trim_poly_g: bool = True,
                    trim_front1: int = 0,
                    umi_config: dict | None = None) -> dict:
    """Run fastp on one lane. Idempotent when its stage manifest is current."""
    name   = sample["name"]
    r1_in  = sample["r1"]
    r2_in  = sample.get("r2")
    paired = sample["paired"]

    r1_out = str(trimmed_dir / f"{name}_R1_trimmed.fq.gz")
    r2_out = str(trimmed_dir / f"{name}_R2_trimmed.fq.gz") if paired else None
    json_out = str(fastp_dir / f"{name}_fastp.json")
    html_out = str(fastp_dir / f"{name}_fastp.html")

    # A4: content-addressed resume manifest for this fastp stage.
    # Preserve exact default keys when protocol/options are at standard defaults
    # so existing stage manifests remain valid; record non-default protocol
    # options so switching protocol/UMI/polyX invalidates stale trimmed reads.
    manifest_path = str(fastp_dir / f"{name}.fastp.stage.json")
    stage_inputs = [("r1", r1_in), ("r2", r2_in)]
    stage_params: dict = {"min_len": min_len, "quality": quality, "paired": paired}
    if library_protocol != "full_length":
        stage_params["library_protocol"] = library_protocol
    if trim_poly_x:
        stage_params["trim_poly_x"] = True
    if not trim_poly_g:
        stage_params["trim_poly_g"] = False
    if trim_front1 > 0:
        stage_params["trim_front1"] = trim_front1
    if umi_config:
        stage_params["umi_config"] = umi_config
    tool_version = _fastp_version()

    result = {
        "name":       name,
        "r1_raw":     r1_in,
        "r2_raw":     r2_in,
        "r1_trimmed": r1_out,
        "r2_trimmed": r2_out,
        "paired":     paired,
    }

    # ── Resume check: outputs valid AND the stage manifest still matches? ──
    manifest_current, manifest_reason = stage_is_current(
        manifest_path, inputs=stage_inputs, params=stage_params,
        tool_version=tool_version,
    )
    if _fastp_outputs_valid(r1_out, r2_out, json_out, paired) and manifest_current:
        try:
            with open(json_out) as f:
                fastp_data = json.load(f)
            summary = fastp_data.get("summary", {})
            before  = summary.get("before_filtering", {})
            after   = summary.get("after_filtering", {})
            n_raw     = int(before.get("total_reads", 0))
            n_trimmed = int(after.get("total_reads",  0))
            pct_pass  = round(n_trimmed / max(n_raw, 1) * 100, 1)
            result.update({
                "status":          "success",
                "n_reads_raw":     n_raw,
                "n_reads_trimmed": n_trimmed,
                "pct_passed":      pct_pass,
                "q30_rate":        round(float(after.get("q30_rate", 0))
                                          * 100, 1),
                "resumed":         True,
            })
            warnings.append(
                f"[resume] fastp skipped for {name} "
                f"(valid outputs exist: {pct_pass}% passed previously)"
            )
            return result
        except Exception:
            # JSON corrupted — fall through to re-run fastp
            pass

    # ADR-075: both single-end and paired-end reads undergo adapter trimming and
    # Phred quality filtering (previously single-end erroneously passed
    # --disable_adapter_trimming and --disable_quality_filtering).
    cmd = [
        "fastp",
        "--in1",  r1_in,
        "--out1", r1_out,
        "--json", json_out,
        "--html", html_out,
        "--thread",      str(threads),
        "--length_required", str(min_len),
        "--qualified_quality_phred", str(quality),
        "--overrepresentation_analysis",
    ]
    if paired and r2_in and r2_out:
        cmd += [
            "--in2", r2_in,
            "--out2", r2_out,
            "--detect_adapter_for_pe",
            "--correction",
        ]
    if trim_poly_g:
        cmd.append("--trim_poly_g")
    if trim_poly_x:
        cmd += ["--trim_poly_x", "--poly_x_min_len", "10"]
    if trim_front1 > 0:
        cmd += ["--trim_front1", str(trim_front1)]
    if umi_config:
        cmd += [
            "--umi",
            "--umi_loc", str(umi_config["umi_loc"]),
            "--umi_len", str(int(umi_config["umi_len"])),
        ]
        if int(umi_config.get("umi_skip", 0) or 0) > 0:
            cmd += ["--umi_skip", str(int(umi_config["umi_skip"]))]
    cmd = [c for c in cmd if c]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600
        )
        if proc.returncode != 0:
            warnings.append(
                f"fastp failed for {name}: {proc.stderr[-200:]}"
            )
            result["status"] = "failed"
            return result

        # Parse fastp JSON for QC metrics
        if Path(json_out).exists():
            with open(json_out) as f:
                fastp_data = json.load(f)
            summary = fastp_data.get("summary", {})
            before  = summary.get("before_filtering", {})
            after   = summary.get("after_filtering", {})

            n_raw     = int(before.get("total_reads", 0))
            n_trimmed = int(after.get("total_reads",  0))
            pct_pass  = round(n_trimmed / max(n_raw, 1) * 100, 1)

            result.update({
                "status":          "success",
                "n_reads_raw":     n_raw,
                "n_reads_trimmed": n_trimmed,
                "pct_passed":      pct_pass,
                "q30_rate":        round(float(
                    after.get("q30_rate", 0)) * 100, 1),
            })

            if pct_pass < 70:
                warnings.append(
                    f"Sample {name}: only {pct_pass}% reads passed QC. "
                    f"Check library quality."
                )
        else:
            result["status"] = "success"

    except subprocess.TimeoutExpired:
        warnings.append(f"fastp timed out for {name}")
        result["status"] = "timeout"
    except FileNotFoundError:
        if allow_mock:
            warnings.append(
                "fastp not found — using mock QC result (explicit mock mode)"
            )
            result["status"] = "missing_tool"
            result.update(_mock_fastp_result(name, r1_in, r2_in,
                                              r1_out, r2_out, paired))
        else:
            warnings.append("fastp not found — FASTQ QC failed.")
            result["status"] = "failed"
            result["error_type"] = "MissingDependency"
            result["details"] = "fastp is required for raw FASTQ preprocessing."

    # A4: record the content-addressed manifest for a genuine successful run so a
    # later resume can tell whether the inputs/params/version still match.
    if result.get("status") == "success":
        try:
            write_stage_manifest(manifest_path, build_stage_manifest(
                stage="fastp", inputs=stage_inputs, params=stage_params,
                tool_version=tool_version,
            ))
        except Exception as exc:
            warnings.append(f"[A4] could not write fastp stage manifest for "
                            f"{name}: {exc}")

    return result


def _run_multiqc(fastp_dir: Path, output_dir: Path,
                  warnings: list) -> Path | None:
    """Run MultiQC over fastp reports."""
    try:
        proc = subprocess.run(
            ["multiqc", str(fastp_dir),
             "--outdir", str(output_dir),
             "--filename", "multiqc_report",
             "--force", "--quiet"],
            capture_output=True, text=True, timeout=300
        )
        report = output_dir / "multiqc_report.html"
        if report.exists():
            return report
        warnings.append(f"MultiQC ran but report not found: {proc.stderr[-100:]}")
    except FileNotFoundError:
        warnings.append("multiqc not found — skipping summary report")
    except Exception as e:
        warnings.append(f"MultiQC failed: {e}")
    return None


def _mock_fastp_result(name, r1_in, r2_in, r1_out, r2_out, paired) -> dict:
    """Mock result when fastp not installed."""
    import random
    import zlib
    # A2: deterministic across processes (builtin hash() is PYTHONHASHSEED-salted).
    random.seed(zlib.crc32(str(name).encode()))
    n_raw = random.randint(20_000_000, 50_000_000)
    pct   = random.uniform(92, 98)
    return {
        "status":          "success",
        "n_reads_raw":     n_raw,
        "n_reads_trimmed": int(n_raw * pct / 100),
        "pct_passed":      round(pct, 1),
        "q30_rate":        round(random.uniform(88, 95), 1),
        "note":            "mock — fastp not installed",
    }


if __name__ == "__main__":
    run_script(rna_fastq_qc)
