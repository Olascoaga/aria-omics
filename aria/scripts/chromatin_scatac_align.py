"""
ARIA scATAC FASTQ -> fragments alignment (chromap)
---------------------------------------------------
Aligns raw single-cell ATAC FASTQ (10x-style: genomic R1/R3 + a barcode read)
to a reference with chromap and emits the barcoded, position-sorted
`fragments.tsv.gz` (+ tabix index) that the fragments->matrix bridge
(chromatin_fragments_to_matrix.py) consumes. Together they take scATAC from raw
sequencer FASTQ all the way to the validated matrix pipeline.

Executed inside aria-atacseq-env (chromap + samtools) via EnvironmentManager.

chromap is the open-source single-cell ATAC aligner (it handles the cell barcode
+ whitelist correction + Tn5 fragment generation in one pass); this is the
single-cell analogue of atac_align.py (bulk, bwa-mem2).

Input params:
  libraries:       list — E5 typed multi-library manifest rows; when present,
                          each row supplies R1/R2/R3 + whitelist + identities
                          and the single-library keys below are ignored
  r1_fastq:        str  — genomic read 1 (10x ATAC R1)
  r3_fastq:        str  — genomic read 2 (10x ATAC R3; "r2_fastq" also accepted
                          for the genomic mate when no separate barcode key)
  barcode_fastq:   str  — cell-barcode read (10x ATAC R2)
  genome_fasta:    str  — reference FASTA
  chromap_index:   str  — (optional) prebuilt chromap index; default <fasta>.index
  barcode_whitelist: str — 10x barcode whitelist (plain or .gz)
  output_dir:      str
  threads:         int  (default 8)

Output:
  {
    "status": "success",
    "fragments_file": str,    — bgzipped, tabix-indexed fragments.tsv.gz
    "aligner": "chromap",
    "n_fragments": int,
    "n_libraries": int,
    "library_manifest": list, — barcode prefix → library/sample/donor metadata
    "warnings": [str]
  }
  Missing chromap/samtools or required inputs -> structured not-run (never a
  fabricated fragments file).
"""

from __future__ import annotations
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def _normalize_fastq_list(val) -> list[str]:
    """Normalize lanes without splitting a literal filename containing commas."""
    if val is None:
        return []
    if isinstance(val, str):
        lanes = [val] if Path(val).is_file() else val.split(",")
    elif isinstance(val, list):
        lanes = val
    else:
        raise ValueError("FASTQs must be a path string or a non-empty list of strings")
    if not lanes:
        raise ValueError("FASTQ lane list must not be empty")
    for lane in lanes:
        if not isinstance(lane, str) or not lane.strip() or not Path(lane).is_file():
            raise ValueError(f"FASTQ lane not found: {lane!r}")
    return lanes


def _open_text(path):
    """Open plain text or gzip by content, tolerating plain test placeholders."""
    with open(path, "rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if compressed else open(path, "rt")


def _whitelist_sequences(path: str) -> list[str]:
    with _open_text(path) as handle:
        return [seq for line in handle
                if (seq := line.strip().upper()) and set(seq) <= set("ACGTN")]


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def detect_barcode_whitelist_orientation(
    barcode_fastqs: list[str], whitelist_path: str, *, max_reads: int = 10000,
    bc_slice: tuple[int, int | None] | None = None, min_match_rate: float = 0.05,
) -> dict:
    """Measure exact direct/RC whitelist matches on valid sampled FASTQ records.

    Low match rates are reported here; the auto-mode caller enforces the
    fail-closed policy. Ties and insufficient evidence retain direct orientation.
    """
    fallback = {"orientation": "direct", "n_sampled": 0, "direct_rate": None,
                "rc_rate": None, "status": "insufficient_reads_fallback"}
    direct_set = set(_whitelist_sequences(whitelist_path))
    if not direct_set:
        return fallback
    rc_set = {_reverse_complement(seq) for seq in direct_set}
    sampled = direct_matches = rc_matches = 0
    for path in barcode_fastqs:
        with _open_text(path) as handle:
            while sampled < max_reads:
                header, sequence, separator, quality = [handle.readline() for _ in range(4)]
                if not header or not quality:
                    break
                seq = sequence.strip().upper()
                if (not header.startswith("@") or not separator.startswith("+")
                        or not seq or len(seq) != len(quality.rstrip("\r\n"))
                        or not set(seq) <= set("ACGTN")):
                    continue
                if bc_slice is not None:
                    seq = seq[slice(*bc_slice)]
                sampled += 1
                direct_matches += seq in direct_set
                rc_matches += seq in rc_set
        if sampled >= max_reads:
            break
    if not sampled:
        return fallback
    direct_rate, rc_rate = direct_matches / sampled, rc_matches / sampled
    return {"orientation": "reverse_complement" if rc_rate > direct_rate else "direct",
            "n_sampled": sampled, "direct_rate": direct_rate, "rc_rate": rc_rate,
            "status": ("whitelist_mismatch" if sampled >= 20
                       and max(direct_rate, rc_rate) < min_match_rate else "sampled")}


def _chromap_version() -> str:
    """Return tool identity without letting a missing executable mask input gates."""
    try:
        proc = subprocess.run(["chromap", "--version"], capture_output=True,
                              text=True, timeout=30)
        return (proc.stdout or proc.stderr or "").strip() if proc.returncode == 0 else ""
    except (OSError, subprocess.TimeoutExpired):
        return ""


def _barcode_source(read_format: str, r1, r3, barcode):
    """Resolve chromap bc:start:length against the available FASTQ streams."""
    for segment in read_format.split(","):
        fields = segment.strip().split(":")
        if len(fields) == 3 and fields[0] == "bc":
            try:
                start, length = int(fields[1]), int(fields[2])
            except ValueError:
                break
            return barcode or r1, (start, None if length < 0 else start + length)
    return barcode, None


def chromatin_scatac_align(params: dict) -> dict:
    if params.get("libraries") is not None:
        return _align_libraries(params)

    from aria.utils.scatac_fastq_manifest import PLATFORM_ROLES, ORIENTATION_ALIASES
    from aria.utils.stage_manifest import build_stage_manifest, write_stage_manifest, stage_is_current

    platform = params.get("platform", "10x_atac")
    read_format = params.get("read_format") or ""
    orientation = ORIENTATION_ALIASES.get(str(params.get("barcode_orientation", "auto")).lower())
    genome_fasta = str(params.get("genome_fasta") or "")
    chromap_index = str(params.get("chromap_index") or (
        f"{genome_fasta}.index" if genome_fasta else ""))
    whitelist = str(params.get("barcode_whitelist") or "")
    output_dir = Path(params.get("output_dir") or "aria_scatac_aligned")
    threads = int(params.get("threads", 8))
    warnings: list[str] = []

    def not_run(reason: str, **extra) -> dict:
        return {"status": "skipped", "ran": False,
                "analysis": "scatac_fastq_to_fragments", "reason": reason,
                "warnings": warnings, **extra}

    if platform not in PLATFORM_ROLES:
        return not_run("invalid_platform", message=f"unknown scATAC platform: {platform!r}")
    if orientation is None:
        return not_run("invalid_barcode_orientation")
    if not isinstance(read_format, str) or (platform == "custom" and not read_format.strip()):
        return not_run("missing_read_format")
    try:
        r1 = _normalize_fastq_list(params.get("r1_fastq"))
        r3 = _normalize_fastq_list(params.get("r3_fastq") or params.get("r2_fastq"))
        barcode = _normalize_fastq_list(params.get("barcode_fastq"))
        if not r1 or (platform != "custom" and not r3):
            raise ValueError("required genomic FASTQ is missing")
        if platform in {"10x_atac", "10x_multiome_atac", "share_seq"} and not barcode:
            raise ValueError("required barcode FASTQ is missing")
        if not barcode and not read_format:
            raise ValueError("inline barcodes require an explicit read_format")
        if (r3 and len(r1) != len(r3)) or (barcode and len(r1) != len(barcode)):
            raise ValueError("genomic and barcode FASTQs must have matching lane counts")
    except ValueError as exc:
        return not_run("missing_input", message=str(exc))
    if not genome_fasta or not Path(genome_fasta).is_file():
        return not_run("missing_genome_fasta", message=f"genome_fasta not found: {genome_fasta}")
    if ((not whitelist and platform in {"10x_atac", "10x_multiome_atac"})
            or (whitelist and not Path(whitelist).is_file())):
        return not_run("missing_barcode_whitelist",
                       message="an existing barcode whitelist is required for this run")

    output_dir.mkdir(parents=True, exist_ok=True)
    fragments_raw = output_dir / "fragments.tsv"
    fragments_gz = output_dir / "fragments.tsv.gz"
    manifest_path = output_dir / "chromatin_scatac_align.stage_manifest.json"
    inputs = [(f"{role}_{i}", path) for role, lanes in (("r1", r1), ("r3", r3), ("bc", barcode))
              for i, path in enumerate(lanes)]
    inputs += [("genome_fasta", genome_fasta), ("barcode_whitelist", whitelist or None)]
    stage_params = {"platform": platform, "read_format": read_format,
                    "barcode_orientation": orientation, "threads": threads,
                    "chromap_index": chromap_index}
    tool_version = _chromap_version() if fragments_gz.exists() else ""

    evidence = {"orientation": "direct" if orientation == "auto" else orientation,
                "n_sampled": 0, "direct_rate": None, "rc_rate": None}
    if whitelist:
        barcode_sources, bc_slice = _barcode_source(read_format, r1, r3, barcode)
        evidence = detect_barcode_whitelist_orientation(barcode_sources, whitelist, bc_slice=bc_slice)
        if orientation == "auto" and evidence.get("status") == "whitelist_mismatch":
            return not_run("barcode_whitelist_mismatch", direct_match_rate=evidence["direct_rate"],
                           rc_match_rate=evidence["rc_rate"], n_sampled=evidence["n_sampled"])
        if orientation != "auto":
            evidence["orientation"] = orientation
    chosen = evidence["orientation"]
    orientation_result = {"barcode_orientation": chosen, "direct_match_rate": evidence["direct_rate"],
                          "rc_match_rate": evidence["rc_rate"], "n_sampled": evidence["n_sampled"]}
    if (_fragments_valid(fragments_gz) and stage_is_current(
            manifest_path, inputs=inputs, params=stage_params, tool_version=tool_version)[0]):
        warnings.append("[resume] reused existing fragments.tsv.gz")
        return {"status": "success", "fragments_file": str(fragments_gz),
                "aligner": "chromap", "resumed": True, "warnings": warnings,
                "validation_level": "production", **orientation_result}

    active_whitelist = whitelist
    if whitelist and chosen == "reverse_complement":
        rc_path = output_dir / "barcode_whitelist.rc.txt"
        rc_path.write_text("".join(_reverse_complement(seq) + "\n"
                                   for seq in _whitelist_sequences(whitelist)))
        active_whitelist = str(rc_path)
        if orientation == "auto":
            warnings.append("Auto-selected reverse-complement barcode whitelist "
                            f"(direct={evidence['direct_rate']}, RC={evidence['rc_rate']})")

    if not Path(chromap_index).is_file():
        idx = _run(["chromap", "-i", "-r", genome_fasta, "-o", chromap_index], warnings, "index")
        if idx == "missing_tool":
            return not_run("chromap_unavailable", message="chromap not found (install aria-atacseq-env)")
        if idx != "ok":
            return not_run("index_build_failed", message=idx)

    cmd = ["chromap", "--preset", "atac", "-x", chromap_index, "-r", genome_fasta]
    for flag, lanes in (("-1", r1), ("-2", r3), ("-b", barcode)):
        for lane in lanes:
            cmd.extend([flag, lane])
    if active_whitelist:
        cmd.extend(["--barcode-whitelist", active_whitelist])
    if read_format:
        cmd.extend(["--read-format", read_format])
    cmd.extend(["-o", str(fragments_raw), "-t", str(threads)])
    # Remove stale output before a rerun so failed tools cannot reuse old fragments.
    fragments_raw.unlink(missing_ok=True)
    res = _run(cmd, warnings, "align")
    if res == "missing_tool":
        return not_run("chromap_unavailable", message="chromap not found (install aria-atacseq-env)")
    if res != "ok":
        return not_run("alignment_failed", message=res)
    if not fragments_raw.is_file() or fragments_raw.stat().st_size == 0:
        return not_run("no_fragments_produced", message="chromap produced no fragments")

    bg = _bgzip_tabix(fragments_raw, fragments_gz, warnings)
    compressed = bg == "ok"
    if not compressed:
        warnings.append(f"bgzip/tabix unavailable ({bg}); leaving plain TSV")
    artifact = fragments_gz if compressed else fragments_raw
    tool_version = _chromap_version()
    write_stage_manifest(manifest_path, build_stage_manifest(
        stage="chromatin_scatac_align", inputs=inputs, params=stage_params, tool_version=tool_version))
    return {"status": "success", "fragments_file": str(artifact), "aligner": "chromap",
            "n_fragments": _count_gz_lines(artifact) if compressed else _count_lines(artifact),
            "compressed": compressed,
            "tabix_indexed": compressed and Path(str(fragments_gz) + ".tbi").exists(),
            "genome_fasta": genome_fasta, "validation_level": "production",
            "warnings": warnings, **orientation_result}


def _align_libraries(params: dict) -> dict:
    """Align every typed library, then build one identity-preserving union."""
    from aria.utils.scatac_fastq_manifest import resolve_scatac_fastq_manifest

    validation = resolve_scatac_fastq_manifest(
        {
            "data_dir": params.get("output_dir") or ".",
            "scatac_fastq_manifest": {
                "schema_version": "1",
                "platform": params.get("platform", "10x_atac"),
                "barcode_orientation": params.get("barcode_orientation", "auto"),
                **({"read_format": params["read_format"]} if params.get("read_format") else {}),
                "libraries": params.get("libraries"),
            },
        },
        require_paths=True,
    )
    if validation.get("status") != "valid":
        return {
            "status": "skipped",
            "ran": False,
            "analysis": "scatac_fastq_to_fragments",
            "reason": "invalid_scatac_fastq_manifest",
            "errors": list(validation.get("errors") or []),
            "warnings": [],
        }

    libraries = validation["manifest"]["libraries"]
    genome_fasta = params.get("genome_fasta")
    output_dir = Path(params.get("output_dir") or "aria_scatac_aligned")
    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for row in libraries:
        roles = row.get("fastq_lanes") or row["fastqs"]
        result = chromatin_scatac_align({
            "r1_fastq": roles["R1"],
            "barcode_fastq": roles.get("R2") if roles.get("R3") else None,
            "r3_fastq": roles.get("R3") or roles.get("R2"),
            "platform": row["platform"],
            "read_format": row.get("read_format"),
            "barcode_orientation": row["barcode_orientation"],
            "genome_fasta": genome_fasta,
            "chromap_index": params.get("chromap_index"),
            "barcode_whitelist": row["barcode_whitelist"],
            "output_dir": str(output_dir / row["library_id"]),
            "threads": params.get("threads", 8),
        })
        results.append({
            **row,
            "status": result.get("status"),
            "fragments_file": result.get("fragments_file"),
            "alignment": result,
        })

    failed = [
        row["library_id"] for row in results
        if row.get("status") != "success" or not row.get("fragments_file")
    ]
    if failed:
        return {
            "status": "skipped",
            "ran": False,
            "analysis": "scatac_fastq_to_fragments",
            "reason": "partial_library_alignment_failure",
            "failed_libraries": failed,
            "library_results": results,
            "warnings": [],
        }

    library_manifest = [
        {
            "library_id": row["library_id"],
            "sample_id": row["sample_id"],
            "donor_id": row["donor_id"],
            "fastqs": row["fastqs"],
            "barcode_whitelist": row["barcode_whitelist"],
            "barcode_orientation": (
                (row.get("alignment") or {}).get("barcode_orientation")
                or row.get("barcode_orientation")
            ),
            "metadata": row["metadata"],
            "barcode_prefix": f"{row['library_id']}#",
            "fragments_file": row["fragments_file"],
        }
        for row in results
    ]
    all_resumed = bool(results) and all(
        (row.get("alignment") or {}).get("resumed") is True for row in results
    )
    combined_gz = output_dir / "combined.fragments.tsv.gz"
    combined_plain = output_dir / "combined.fragments.tsv"
    manifest_path = output_dir / "library_manifest.json"
    if all_resumed and (_fragments_valid(combined_gz) or combined_plain.is_file()):
        fragments_file = combined_gz if _fragments_valid(combined_gz) else combined_plain
        n_fragments = (
            _count_gz_lines(combined_gz)
            if fragments_file == combined_gz
            else _count_lines(combined_plain)
        )
        return {
            "status": "success",
            "ran": True,
            "resumed": True,
            "analysis": "scatac_fastq_to_fragments",
            "fragments_file": str(fragments_file),
            "aligner": "chromap",
            "n_fragments": n_fragments,
            "n_libraries": len(library_manifest),
            "library_manifest": library_manifest,
            "library_manifest_path": str(manifest_path),
            "library_results": results,
            "compressed": fragments_file == combined_gz,
            "warnings": ["[resume] reused existing combined fragments"],
            "validation_level": "production",
        }

    merged = _merge_library_fragments(results, combined_plain)
    if merged.get("status") != "success":
        return {
            "status": "skipped",
            "ran": False,
            "analysis": "scatac_fastq_to_fragments",
            "reason": "fragment_union_failed",
            "message": merged.get("message"),
            "library_results": results,
            "warnings": [],
        }

    manifest_path.write_text(json.dumps({
        "schema_version": "1",
        "libraries": library_manifest,
    }, indent=2, sort_keys=True), encoding="utf-8")

    warnings = []
    compressed = _bgzip_tabix(combined_plain, combined_gz, warnings)
    if compressed == "ok":
        fragments_file = combined_gz
        n_fragments = _count_gz_lines(combined_gz)
    else:
        warnings.append(
            f"bgzip/tabix unavailable ({compressed}); leaving plain TSV"
        )
        fragments_file = combined_plain
        n_fragments = _count_lines(combined_plain)

    return {
        "status": "success",
        "ran": True,
        "resumed": False,
        "analysis": "scatac_fastq_to_fragments",
        "fragments_file": str(fragments_file),
        "aligner": "chromap",
        "n_fragments": n_fragments,
        "n_libraries": len(library_manifest),
        "library_manifest": library_manifest,
        "library_manifest_path": str(manifest_path),
        "library_results": results,
        "compressed": compressed == "ok",
        "warnings": warnings,
        "validation_level": "production",
    }


def _merge_library_fragments(records: list[dict], output_path: Path) -> dict:
    """Union fragments and prefix barcodes so library identity cannot collide."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    unsorted = output_path.with_suffix(output_path.suffix + ".unsorted")
    prefixes = {
        str(row["library_id"]): f"{row['library_id']}#" for row in records
    }
    try:
        with open(unsorted, "w", encoding="utf-8") as destination:
            for row in records:
                source = Path(row["fragments_file"])
                opener = gzip.open if source.suffix == ".gz" else open
                with opener(source, "rt", encoding="utf-8") as handle:
                    for line in handle:
                        if not line.strip() or line.startswith("#"):
                            continue
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) < 4:
                            return {
                                "status": "error",
                                "message": f"invalid fragments row in {source}",
                            }
                        prefix = prefixes[str(row["library_id"])]
                        if not fields[3].startswith(prefix):
                            fields[3] = prefix + fields[3]
                        destination.write("\t".join(fields) + "\n")
        try:
            proc = subprocess.run(
                ["sort", "-k1,1", "-k2,2n", "-k3,3n", "-o",
                 str(output_path), str(unsorted)],
                capture_output=True,
                text=True,
                timeout=86400,
                env={**os.environ, "LC_ALL": "C"},
            )
            if proc.returncode != 0:
                _python_sort_fragments(unsorted, output_path)
        except (OSError, subprocess.TimeoutExpired):
            _python_sort_fragments(unsorted, output_path)
        return {
            "status": "success",
            "output_path": str(output_path),
            "barcode_prefixes": prefixes,
        }
    except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
        return {"status": "error", "message": str(exc)}
    finally:
        try:
            unsorted.unlink()
        except OSError:
            pass


def _python_sort_fragments(source: Path, destination: Path) -> None:
    """Sort coordinates locally when the external sort executable is unavailable."""
    with source.open() as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle]
    rows.sort(key=lambda row: (row[0], int(row[1]), int(row[2])))
    with destination.open("w") as handle:
        for row in rows:
            handle.write("\t".join(row) + "\n")


def _bgzip_tabix(plain: Path, gz: Path, warnings: list) -> str:
    try:
        with open(gz, "wb") as out:
            p = subprocess.run(["bgzip", "-c", str(plain)], stdout=out,
                               stderr=subprocess.PIPE, timeout=3600)
        if p.returncode != 0:
            return p.stderr.decode(errors="replace")[-200:]
        t = subprocess.run(["tabix", "-f", "-p", "bed", str(gz)],
                           capture_output=True, text=True, timeout=1800)
        if t.returncode != 0:
            return t.stderr[-200:]
        try:
            plain.unlink()
        except OSError:
            pass
        return "ok"
    except FileNotFoundError:
        return "missing_tool"
    except subprocess.TimeoutExpired:
        return "bgzip/tabix timed out"


def _run(cmd: list, warnings: list, step: str) -> str:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
        if p.returncode != 0:
            warnings.append(f"{step} failed: {p.stderr[-200:]}")
            return p.stderr[-300:]
        return "ok"
    except FileNotFoundError:
        return "missing_tool"
    except subprocess.TimeoutExpired:
        warnings.append(f"{step} timed out")
        return f"{step} timed out"


def _fragments_valid(gz: Path) -> bool:
    """Check for a readable genomic fragment rather than using file size as QC."""
    try:
        with gzip.open(gz, "rt") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                return (len(fields) >= 4 and bool(fields[0]) and bool(fields[3])
                        and 0 <= int(fields[1]) < int(fields[2]))
        return False
    except (OSError, EOFError, ValueError, UnicodeError):
        return False


def _count_lines(path: Path) -> int:
    try:
        with open(path) as f:
            return sum(1 for _ in f)
    except OSError:
        return 0


def _count_gz_lines(path: Path) -> int:
    import gzip
    try:
        with gzip.open(path, "rt") as f:
            return sum(1 for _ in f)
    except OSError:
        return 0


if __name__ == "__main__":
    run_script(chromatin_scatac_align)
