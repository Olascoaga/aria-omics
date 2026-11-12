"""Universal raw FASTQ protocol & read-geometry contracts for RNA-seq (ADR-075).

Supports five raw-FASTQ transcriptomic families while failing closed on
out-of-scope library architectures:

1. ``full_length``       — TruSeq Stranded, NEBNext Ultra II, Ribo-Zero, SMART-Seq bulk
2. ``tag_3prime``        — Lexogen QuantSeq 3' FWD/REV, Takara SMART-3Seq, TruSeq 3'
3. ``tag_5prime``        — 5' CAGE, nanoCAGE, RAMPAGE
4. ``umi_bulk``          — Per-sample inline-UMI bulk (QuantSeq-UMI, QIAseq UPX, SMARTer UMI)
5. ``pooled_plate_bulk`` — Multiplexed well-barcode + UMI plate bulk (BRB-seq, DRUG-seq,
                           Prime-seq, QuantSeq-Pool, CEL-Seq2 bulk plate)

Out-of-scope architectures rejected fail-closed at DataAudit / AssayContract:
- Long-read transcriptomics (Oxford Nanopore / PacBio Iso-Seq, median read length >= 600 bp)
- Immune-receptor V(D)J TCR/BCR repertoire sequencing
- Unmanifested asymmetric barcode+UMI paired FASTQs (R1 = 14-36 bp, R2 >= 45 bp) routed to
  per-sample bulk alignment without ``pooled_plate_bulk`` or ``scRNA`` ``kb count`` config
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import statistics
from pathlib import Path
from typing import Any


BULK_RNA_PROTOCOLS: tuple[str, ...] = (
    "full_length",
    "tag_3prime",
    "tag_5prime",
    "umi_bulk",
    "pooled_plate_bulk",
)

PROTOCOL_ALIASES: dict[str, str] = {
    # Family 1: Full-length bulk RNA-seq
    "full_length": "full_length",
    "standard": "full_length",
    "truseq": "full_length",
    "truseq_stranded": "full_length",
    "nebnext": "full_length",
    "ribozero": "full_length",
    "ribo_zero": "full_length",
    "smartseq_bulk": "full_length",
    "polya_full": "full_length",
    # Family 2a: 3' end-tag bulk (no UMI)
    "tag_3prime": "tag_3prime",
    "3prime": "tag_3prime",
    "3p": "tag_3prime",
    "3prime_tag": "tag_3prime",
    "quantseq": "tag_3prime",
    "quantseq_fwd": "tag_3prime",
    "quantseq_rev": "tag_3prime",
    "smart3seq": "tag_3prime",
    "smart_3seq": "tag_3prime",
    "truseq_3p": "tag_3prime",
    # Family 2b: 5' end-tag bulk (no UMI)
    "tag_5prime": "tag_5prime",
    "5prime": "tag_5prime",
    "5p": "tag_5prime",
    "5prime_tag": "tag_5prime",
    "cage": "tag_5prime",
    "nanocage": "tag_5prime",
    "rampage": "tag_5prime",
    # Family 3: Per-sample inline-UMI bulk
    "umi_bulk": "umi_bulk",
    "bulk_umi": "umi_bulk",
    "quantseq_umi": "umi_bulk",
    "qiaseq_upx": "umi_bulk",
    "smarter_umi": "umi_bulk",
    "nanocage_umi": "umi_bulk",
    # Family 4: Multiplexed well-barcode + UMI plate bulk (STARsolo)
    "pooled_plate_bulk": "pooled_plate_bulk",
    "brbseq": "pooled_plate_bulk",
    "brb_seq": "pooled_plate_bulk",
    "drugseq": "pooled_plate_bulk",
    "drug_seq": "pooled_plate_bulk",
    "primeseq": "pooled_plate_bulk",
    "prime_seq": "pooled_plate_bulk",
    "quantseq_pool": "pooled_plate_bulk",
    "celseq2_bulk": "pooled_plate_bulk",
    "strt_bulk": "pooled_plate_bulk",
}

POOLED_PLATE_PRESETS: dict[str, dict[str, Any]] = {
    "brbseq": {
        "preset": "brbseq",
        "cb_start": 1,
        "cb_len": 14,
        "umi_start": 15,
        "umi_len": 14,
        "cdna_read": "R2",
        "strand": "Forward",
    },
    "drugseq": {
        "preset": "drugseq",
        "cb_start": 1,
        "cb_len": 10,
        "umi_start": 11,
        "umi_len": 10,
        "cdna_read": "R2",
        "strand": "Forward",
    },
    "primeseq": {
        "preset": "primeseq",
        "cb_start": 1,
        "cb_len": 12,
        "umi_start": 13,
        "umi_len": 16,
        "cdna_read": "R2",
        "strand": "Forward",
    },
    "quantseq_pool": {
        "preset": "quantseq_pool",
        "cb_start": 1,
        "cb_len": 12,
        "umi_start": 13,
        "umi_len": 10,
        "cdna_read": "R2",
        "strand": "Forward",
    },
    "celseq2_bulk": {
        "preset": "celseq2_bulk",
        "cb_start": 7,
        "cb_len": 6,
        "umi_start": 1,
        "umi_len": 6,
        "cdna_read": "R2",
        "strand": "Forward",
    },
}

VALID_FASTP_UMI_LOCS: frozenset[str] = frozenset({
    "read1",
    "read2",
    "per_read",
    "index1",
    "index2",
    "per_index",
})

_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
_BARCODE_SEQ_RE = re.compile(r"^[ACGTN]+$")
_VDJ_TOKEN_RE = re.compile(
    r"(?:^|[_\-./])(?:vdj|vdj[_\-]?[tb]|tcr[_\-]?seq|bcr[_\-]?seq|"
    r"ig[_\-]?repertoire|immune[_\-]?repertoire)(?:[_\-./]|$)",
    re.IGNORECASE,
)
_POLY_A_RE = re.compile(r"A{10,}")
_POLY_T_START_RE = re.compile(r"^T{10,}")
_POLY_G_END_RE = re.compile(r"G{10,}$")


def _normalize_protocol_token(raw: Any) -> str:
    text = str(raw or "").strip().lower()
    text = re.sub(r"[\s\-]+", "_", text)
    return text


def is_end_tag_or_umi_bias_protocol(
    library_protocol: str | None,
    umi_config: dict[str, Any] | None = None,
) -> bool:
    """Return True when transcript counts must NOT be length-normalized to TPM.

    3'-tag, 5'-tag, and pooled plate UMI protocols generate ~1 tag/molecule per
    transcript regardless of transcript length; dividing by exon-union kb biases
    against long transcripts. ``umi_bulk`` defaults to tag-based unless the
    caller explicitly declares ``umi_config={"is_tag_assay": False}``.
    """
    proto = PROTOCOL_ALIASES.get(
        _normalize_protocol_token(library_protocol or "full_length"),
        str(library_protocol or "full_length"),
    )
    if proto in {"tag_3prime", "tag_5prime", "pooled_plate_bulk"}:
        return True
    if proto == "umi_bulk":
        cfg = umi_config or {}
        return bool(cfg.get("is_tag_assay", True))
    return False


def resolve_bulk_rna_protocol(exp_context: dict[str, Any] | None) -> dict[str, Any]:
    """Resolve and validate the bulk RNA raw-FASTQ protocol configuration."""
    ctx = exp_context or {}
    genome_cfg = ctx.get("genome_config") or {}
    design = ctx.get("design") or {}

    raw_proto = (
        ctx.get("rna_library_protocol")
        or genome_cfg.get("library_protocol")
        or design.get("rna_library_protocol")
    )

    # Auto-detect presence of pooled_plate_manifest.json or explicit dict in context
    manifest_in = (
        ctx.get("pooled_plate_manifest")
        or genome_cfg.get("pooled_plate_manifest")
        or _load_sidecar_pooled_manifest(ctx)
    )
    umi_in = ctx.get("umi_config") or genome_cfg.get("umi_config")

    if not raw_proto:
        if manifest_in:
            raw_proto = "pooled_plate_bulk"
            explicit = True
        elif isinstance(umi_in, dict) and umi_in:
            raw_proto = "umi_bulk"
            explicit = True
        else:
            return {
                "status": "valid",
                "protocol": "full_length",
                "raw_token": None,
                "explicit": False,
                "suppress_tpm": False,
                "umi_config": None,
                "pooled_plate_manifest": None,
                "errors": [],
            }
    else:
        explicit = True

    token = _normalize_protocol_token(raw_proto)
    canonical = PROTOCOL_ALIASES.get(token)
    if canonical is None:
        return {
            "status": "invalid",
            "protocol": None,
            "raw_token": str(raw_proto),
            "explicit": True,
            "suppress_tpm": False,
            "umi_config": None,
            "pooled_plate_manifest": None,
            "errors": [
                f"Unsupported rna_library_protocol {raw_proto!r}. Supported "
                f"protocols: {', '.join(BULK_RNA_PROTOCOLS)} (plus aliases such "
                "as 'quantseq', 'cage', 'quantseq_umi', 'brbseq', 'drugseq', "
                "'primeseq')."
            ],
        }

    errors: list[str] = []
    resolved_umi: dict[str, Any] | None = None
    resolved_manifest: dict[str, Any] | None = None

    if canonical == "umi_bulk":
        resolved_umi, umi_errors = validate_umi_bulk_config(umi_in, raw_token=token)
        errors.extend(umi_errors)
    elif canonical == "pooled_plate_bulk":
        preset_hint = token if token in POOLED_PLATE_PRESETS else None
        resolved_manifest, manifest_errors = validate_pooled_plate_manifest(
            manifest_in,
            preset_hint=preset_hint,
        )
        errors.extend(manifest_errors)

    return {
        "status": "invalid" if errors else "valid",
        "protocol": canonical,
        "raw_token": token,
        "explicit": explicit,
        "suppress_tpm": is_end_tag_or_umi_bias_protocol(canonical, resolved_umi),
        "umi_config": resolved_umi,
        "pooled_plate_manifest": resolved_manifest,
        "errors": errors,
    }


def _load_sidecar_pooled_manifest(exp_context: dict[str, Any]) -> dict[str, Any] | None:
    data_dir = exp_context.get("data_dir")
    if not data_dir:
        return None
    candidate = Path(str(data_dir)) / "pooled_plate_manifest.json"
    if not candidate.is_file():
        return None
    try:
        loaded = json.loads(candidate.read_text(encoding="utf-8"))
        return loaded if isinstance(loaded, dict) else None
    except Exception:
        return {"_parse_error": f"Invalid JSON in {candidate}"}


def validate_umi_bulk_config(
    umi_config: Any,
    *,
    raw_token: str = "umi_bulk",
) -> tuple[dict[str, Any] | None, list[str]]:
    """Validate inline-UMI configuration for ``umi_bulk``."""
    preset_defaults: dict[str, dict[str, Any]] = {
        "quantseq_umi": {
            "umi_loc": "read1",
            "umi_len": 6,
            "umi_skip": 4,
            "is_tag_assay": True,
        },
        "nanocage_umi": {
            "umi_loc": "read1",
            "umi_len": 8,
            "umi_skip": 0,
            "is_tag_assay": True,
        },
    }
    cfg: dict[str, Any] = {}
    if raw_token in preset_defaults:
        cfg.update(preset_defaults[raw_token])
    if isinstance(umi_config, dict):
        cfg.update(umi_config)
    elif umi_config is not None:
        return None, ["umi_config must be a dictionary specifying umi_loc and umi_len."]

    if not cfg:
        return None, [
            "Protocol 'umi_bulk' requires explicit umi_config={'umi_loc': "
            "'read1'|'read2'|'per_read'|'index1'|'index2', 'umi_len': <int>, "
            "'umi_skip': <int>} (or preset 'quantseq_umi' / 'nanocage_umi')."
        ]

    errors: list[str] = []
    umi_loc = str(cfg.get("umi_loc") or "read1").strip().lower()
    if umi_loc not in VALID_FASTP_UMI_LOCS:
        errors.append(
            f"umi_config.umi_loc must be one of {sorted(VALID_FASTP_UMI_LOCS)}; "
            f"got {umi_loc!r}."
        )
    try:
        umi_len = int(cfg.get("umi_len", 0))
    except (TypeError, ValueError):
        umi_len = 0
    if umi_len < 4 or umi_len > 32:
        errors.append(f"umi_config.umi_len must be an integer in [4, 32]; got {cfg.get('umi_len')!r}.")

    try:
        umi_skip = int(cfg.get("umi_skip", 0))
    except (TypeError, ValueError):
        umi_skip = -1
    if umi_skip < 0 or umi_skip > 32:
        errors.append(f"umi_config.umi_skip must be an integer in [0, 32]; got {cfg.get('umi_skip')!r}.")

    is_tag = bool(cfg.get("is_tag_assay", True))
    trim_poly_x = bool(cfg.get("trim_poly_x", is_tag))

    if errors:
        return None, errors
    return {
        "umi_loc": umi_loc,
        "umi_len": umi_len,
        "umi_skip": umi_skip,
        "is_tag_assay": is_tag,
        "trim_poly_x": trim_poly_x,
    }, []


def _load_barcode_map_file(path: Path) -> tuple[dict[str, str], list[str]]:
    if not path.is_file():
        return {}, [f"barcode_to_sample file does not exist: {path}"]
    mapping: dict[str, str] = {}
    errors: list[str] = []
    try:
        lines = [
            ln.strip()
            for ln in path.read_text(encoding="utf-8").splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")
        ]
        if not lines:
            return {}, [f"barcode_to_sample file is empty: {path}"]
        delim = "\t" if "\t" in lines[0] else ","
        reader = csv.reader(lines, delimiter=delim)
        for idx, row in enumerate(reader):
            if len(row) < 2:
                errors.append(f"Row {idx + 1} in {path} must have barcode and sample_id columns.")
                continue
            bc, sample = row[0].strip(), row[1].strip()
            if idx == 0 and bc.lower() in {"barcode", "well_barcode", "cb"}:
                continue
            mapping[bc] = sample
    except Exception as exc:
        return {}, [f"Could not parse barcode_to_sample file {path}: {exc}"]
    return mapping, errors


def validate_pooled_plate_manifest(
    manifest: Any,
    *,
    preset_hint: str | None = None,
) -> tuple[dict[str, Any] | None, list[str]]:
    """Validate a ``pooled_plate_bulk`` well-barcode + UMI manifest."""
    if not isinstance(manifest, dict) or not manifest:
        return None, [
            "Protocol 'pooled_plate_bulk' (BRB-seq / DRUG-seq / Prime-seq / "
            "QuantSeq-Pool) requires a pooled_plate_manifest with "
            "'barcode_to_sample' (mapping well barcode -> sample_id) and CB/UMI "
            "coordinates (or 'preset': 'brbseq'|'drugseq'|'primeseq'|'quantseq_pool')."
        ]
    if "_parse_error" in manifest:
        return None, [str(manifest["_parse_error"])]

    errors: list[str] = []
    preset_name = _normalize_protocol_token(manifest.get("preset") or preset_hint or "")
    merged: dict[str, Any] = {}
    if preset_name:
        if preset_name not in POOLED_PLATE_PRESETS:
            errors.append(
                f"Unknown pooled_plate_manifest preset {preset_name!r}; "
                f"supported presets: {sorted(POOLED_PLATE_PRESETS)}."
            )
        else:
            merged.update(POOLED_PLATE_PRESETS[preset_name])
    merged.update(manifest)

    raw_map = merged.get("barcode_to_sample")
    barcode_to_sample: dict[str, str] = {}
    if isinstance(raw_map, (str, Path)):
        barcode_to_sample, map_errors = _load_barcode_map_file(Path(str(raw_map)))
        errors.extend(map_errors)
    elif isinstance(raw_map, dict):
        barcode_to_sample = {
            str(k).strip().upper(): str(v).strip()
            for k, v in raw_map.items()
            if str(k).strip()
        }
    else:
        errors.append(
            "pooled_plate_manifest.barcode_to_sample is required (dict mapping "
            "barcode -> sample_id, or path to a 2-column TSV/CSV)."
        )

    try:
        cb_start = int(merged.get("cb_start", 0))
        cb_len = int(merged.get("cb_len", 0))
        umi_start = int(merged.get("umi_start", 0))
        umi_len = int(merged.get("umi_len", 0))
    except (TypeError, ValueError):
        cb_start = cb_len = umi_start = umi_len = 0

    if cb_start < 1 or cb_len < 4:
        errors.append(
            f"Invalid cell/well barcode coordinates: cb_start={cb_start}, cb_len={cb_len} "
            "(1-based coordinates with cb_len >= 4 required)."
        )
    if umi_start < 1 or umi_len < 4:
        errors.append(
            f"Invalid UMI coordinates: umi_start={umi_start}, umi_len={umi_len} "
            "(1-based coordinates with umi_len >= 4 required)."
        )
    if cb_start >= 1 and cb_len >= 1 and umi_start >= 1 and umi_len >= 1:
        cb_end = cb_start + cb_len - 1
        umi_end = umi_start + umi_len - 1
        if not (cb_end < umi_start or umi_end < cb_start):
            errors.append(
                f"Overlapping barcode [{cb_start}, {cb_end}] and UMI "
                f"[{umi_start}, {umi_end}] intervals on the barcode read."
            )

    cdna_read = str(merged.get("cdna_read") or "R2").strip().upper()
    if cdna_read not in {"R1", "R2"}:
        errors.append(f"pooled_plate_manifest.cdna_read must be 'R1' or 'R2'; got {cdna_read!r}.")

    strand = str(merged.get("strand") or "Forward").strip().capitalize()
    if strand not in {"Forward", "Reverse"}:
        errors.append(
            f"pooled_plate_manifest.strand must be 'Forward' or 'Reverse'; got {strand!r}."
        )

    if barcode_to_sample:
        seen_samples: set[str] = set()
        for bc, sample_id in sorted(barcode_to_sample.items()):
            if not sample_id:
                errors.append(f"Barcode {bc!r} maps to an empty sample_id.")
                continue
            if sample_id in seen_samples:
                errors.append(
                    f"Duplicate sample_id {sample_id!r} in barcode_to_sample; "
                    "each well barcode must map to a distinct library/sample column."
                )
            seen_samples.add(sample_id)
            if not _BARCODE_SEQ_RE.fullmatch(bc):
                errors.append(f"Invalid nucleotide barcode sequence {bc!r}.")
            elif cb_len >= 4 and len(bc) != cb_len:
                errors.append(
                    f"Barcode {bc!r} has length {len(bc)} != declared cb_len={cb_len}."
                )
        if len(barcode_to_sample) < 2:
            errors.append(
                "pooled_plate_manifest.barcode_to_sample must declare at least 2 barcoded samples."
            )

    whitelist_path = merged.get("barcode_whitelist")
    if whitelist_path is not None and str(whitelist_path).strip():
        whitelist_path = str(whitelist_path).strip()
    else:
        whitelist_path = None

    if errors:
        return None, errors

    return {
        "preset": preset_name or None,
        "cb_start": cb_start,
        "cb_len": cb_len,
        "umi_start": umi_start,
        "umi_len": umi_len,
        "cdna_read": cdna_read,
        "barcode_read": "R1" if cdna_read == "R2" else "R2",
        "strand": strand,
        "barcode_to_sample": dict(sorted(barcode_to_sample.items())),
        "barcode_whitelist": whitelist_path,
        "r1_fastqs": [str(p) for p in (merged.get("r1_fastqs") or [])],
        "r2_fastqs": [str(p) for p in (merged.get("r2_fastqs") or [])],
    }, []


def _read_fastq_lengths_and_tails(
    path: Path,
    max_reads: int = 1000,
) -> dict[str, Any]:
    """Stream up to ``max_reads`` valid 4-line FASTQ records without loading the file."""
    if not path.is_file():
        return {"n_reads": 0, "lengths": [], "poly_a_hits": 0, "poly_g_hits": 0}
    try:
        if path.stat().st_size == 0:
            return {"n_reads": 0, "lengths": [], "poly_a_hits": 0, "poly_g_hits": 0}
    except OSError:
        return {"n_reads": 0, "lengths": [], "poly_a_hits": 0, "poly_g_hits": 0}

    opener = gzip.open if path.name.lower().endswith(".gz") else open
    lengths: list[int] = []
    poly_a_hits = 0
    poly_g_hits = 0
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            while len(lengths) < max_reads:
                hdr = fh.readline()
                if not hdr:
                    break
                seq = fh.readline()
                plus = fh.readline()
                qual = fh.readline()
                if not (seq and plus and qual):
                    break
                hdr_s = hdr.strip()
                seq_s = seq.strip().upper()
                plus_s = plus.strip()
                qual_s = qual.strip()
                if (
                    not hdr_s.startswith("@")
                    or not plus_s.startswith("+")
                    or not seq_s
                    or len(seq_s) != len(qual_s)
                ):
                    # Synthetic placeholder file (e.g. b"fastq" in unit tests)
                    return {"n_reads": 0, "lengths": [], "poly_a_hits": 0, "poly_g_hits": 0}
                lengths.append(len(seq_s))
                if _POLY_A_RE.search(seq_s[-35:]) or _POLY_T_START_RE.search(seq_s[:35]):
                    poly_a_hits += 1
                if _POLY_G_END_RE.search(seq_s[-25:]):
                    poly_g_hits += 1
    except Exception:
        return {"n_reads": 0, "lengths": [], "poly_a_hits": 0, "poly_g_hits": 0}

    return {
        "n_reads": len(lengths),
        "lengths": lengths,
        "poly_a_hits": poly_a_hits,
        "poly_g_hits": poly_g_hits,
    }


def _classify_read_slot(filename: str) -> str:
    stem = filename
    for ext in _FASTQ_SUFFIXES:
        if stem.lower().endswith(ext):
            stem = stem[:-len(ext)]
            break
    if re.search(r"_(?:R2|2)(?:_\d+)?$", stem, re.IGNORECASE):
        return "R2"
    if re.search(r"_(?:I[12]|R3)(?:_\d+)?$", stem, re.IGNORECASE):
        return "INDEX"
    return "R1"


def inspect_fastq_geometry(
    files: list[str] | tuple[str, ...],
    *,
    max_reads_per_slot: int = 1000,
) -> dict[str, Any]:
    """Inspect real FASTQ headers/sequences to classify library read geometry.

    Returns ``status="uninspected"`` when files are missing, empty, or synthetic
    placeholders so lightweight unit tests are unaffected.
    """
    fq_paths = [
        Path(str(p))
        for p in (files or [])
        if str(p).lower().endswith(_FASTQ_SUFFIXES)
    ]
    vdj_files = [
        str(p) for p in fq_paths if _VDJ_TOKEN_RE.search(str(p))
    ]
    if not fq_paths:
        return {
            "status": "not_fastq",
            "geometry": "not_fastq",
            "n_reads_sampled": 0,
            "vdj_repertoire_detected": False,
            "vdj_files": [],
        }

    r1_lengths: list[int] = []
    r2_lengths: list[int] = []
    poly_a_hits = 0
    poly_g_hits = 0

    r1_candidates = [p for p in fq_paths if _classify_read_slot(p.name) == "R1"][:2]
    r2_candidates = [p for p in fq_paths if _classify_read_slot(p.name) == "R2"][:2]

    per_file_budget = max(200, max_reads_per_slot // max(1, len(r1_candidates)))
    for p in r1_candidates:
        stats = _read_fastq_lengths_and_tails(p, max_reads=per_file_budget)
        r1_lengths.extend(stats["lengths"])
        poly_a_hits += stats["poly_a_hits"]
        poly_g_hits += stats["poly_g_hits"]

    for p in r2_candidates:
        stats = _read_fastq_lengths_and_tails(p, max_reads=per_file_budget)
        r2_lengths.extend(stats["lengths"])

    n_r1 = len(r1_lengths)
    n_r2 = len(r2_lengths)
    total_sampled = n_r1 + n_r2
    if total_sampled == 0:
        return {
            "status": "uninspected",
            "geometry": "uninspected",
            "n_reads_sampled": 0,
            "r1_median_len": None,
            "r2_median_len": None,
            "r1_poly_a_frac": 0.0,
            "r1_poly_g_frac": 0.0,
            "vdj_repertoire_detected": bool(vdj_files),
            "vdj_files": vdj_files,
        }

    r1_med = float(statistics.median(r1_lengths)) if r1_lengths else None
    r2_med = float(statistics.median(r2_lengths)) if r2_lengths else None
    poly_a_frac = round(poly_a_hits / max(n_r1, 1), 4)
    poly_g_frac = round(poly_g_hits / max(n_r1, 1), 4)

    max_med = max(v for v in (r1_med, r2_med) if v is not None)
    geometry = "standard_short_read"
    if max_med >= 600.0:
        geometry = "long_read"
    elif (
        r1_med is not None
        and r2_med is not None
        and n_r1 >= 10
        and n_r2 >= 10
        and (
            (14.0 <= r1_med <= 36.0 and r2_med >= 45.0)
            or (14.0 <= r2_med <= 36.0 and r1_med >= 45.0)
        )
    ):
        geometry = "asymmetric_barcode_umi"
    elif n_r1 >= 10 and poly_a_frac >= 0.10:
        geometry = "3prime_poly_a_enriched"

    return {
        "status": "inspected",
        "geometry": geometry,
        "n_reads_sampled": total_sampled,
        "n_r1_sampled": n_r1,
        "n_r2_sampled": n_r2,
        "r1_median_len": r1_med,
        "r2_median_len": r2_med,
        "r1_poly_a_frac": poly_a_frac,
        "r1_poly_g_frac": poly_g_frac,
        "vdj_repertoire_detected": bool(vdj_files),
        "vdj_files": vdj_files,
    }


def validate_rna_fastq_contract(
    exp_context: dict[str, Any],
    modality: str,
    files: list[str] | tuple[str, ...],
) -> dict[str, Any]:
    """Validate RNA FASTQ protocol & read geometry before expensive dispatch."""
    fq_files = [
        str(p) for p in (files or [])
        if str(p).lower().endswith(_FASTQ_SUFFIXES)
    ]
    if modality not in {"bulk_RNA", "bulk_RNA_raw", "scRNA"} or not fq_files:
        return {"status": "not_applicable", "findings": []}

    geom = inspect_fastq_geometry(fq_files)
    findings: list[dict[str, Any]] = []

    if geom.get("vdj_repertoire_detected"):
        findings.append({
            "severity": "blocking",
            "check": "assay_contract_vdj_repertoire_unsupported",
            "message": (
                f"{modality} FASTQ input matches V(D)J TCR/BCR immune-repertoire "
                f"naming ({Path(geom['vdj_files'][0]).name}). ARIA analyzes "
                "gene-expression transcriptomics, not clonotype receptor assembly."
            ),
            "recommendation": (
                "Analyze V(D)J repertoire libraries with cellranger vdj, MiXCR, "
                "or TRUST4, and pass only gene-expression (GEX) libraries to ARIA."
            ),
            "modality": modality,
        })

    if geom.get("geometry") == "long_read":
        max_bp = max(
            v for v in (geom.get("r1_median_len"), geom.get("r2_median_len"))
            if v is not None
        )
        findings.append({
            "severity": "blocking",
            "check": "assay_contract_long_read_rna_unsupported",
            "message": (
                f"{modality} FASTQ has long-read geometry (median read length "
                f"{max_bp:.0f} bp >= 600 bp, typical of Oxford Nanopore or "
                "PacBio Iso-Seq). Short-read STAR/kallisto pipelines are not "
                "valid for long-read transcriptomics."
            ),
            "recommendation": (
                "Align long reads with minimap2 (-ax splice) and quantify with "
                "IsoQuant, Bambu, or FLAMES, then supply the resulting gene "
                "count matrix to ARIA."
            ),
            "modality": modality,
        })

    if modality == "scRNA":
        return {
            "status": "blocked" if findings else "pass",
            "geometry": geom,
            "protocol": "scrna_kb_count",
            "findings": findings,
        }

    proto = resolve_bulk_rna_protocol(exp_context)
    if proto["status"] == "invalid":
        findings.append({
            "severity": "blocking",
            "check": "assay_contract_rna_protocol_invalid",
            "message": "; ".join(proto["errors"]),
            "recommendation": (
                "Set `rna_library_protocol` to one of "
                f"{list(BULK_RNA_PROTOCOLS)} and supply required UMI or "
                "pooled-plate barcode manifest fields."
            ),
            "modality": modality,
        })

    if (
        geom.get("geometry") == "asymmetric_barcode_umi"
        and proto.get("protocol") not in {"pooled_plate_bulk", "umi_bulk"}
    ):
        findings.append({
            "severity": "blocking",
            "check": "assay_contract_asymmetric_barcode_umi_fastq",
            "message": (
                f"{modality} paired FASTQs have asymmetric barcode+UMI read "
                f"geometry (R1 median={geom.get('r1_median_len'):.0f} bp, "
                f"R2 median={geom.get('r2_median_len'):.0f} bp). Treating a "
                "short barcode+UMI read as a genomic paired-end mate in STAR "
                "corrupts alignment."
            ),
            "recommendation": (
                "If this is single-cell 3'/5' RNA-seq (e.g. 10x Chromium, "
                "Drop-seq, CEL-Seq2), confirm modality `scRNA` with "
                "`raw_ingestion_kb` chemistry (`10XV3`, `10X_5P`, `CELSEQ2`, "
                "etc.). If this is multiplexed well-barcoded bulk RNA-seq "
                "(BRB-seq, DRUG-seq, Prime-seq, QuantSeq-Pool), set "
                "`rna_library_protocol='pooled_plate_bulk'` and provide "
                "`pooled_plate_manifest`."
            ),
            "modality": modality,
        })

    advisory = None
    if (
        geom.get("geometry") == "3prime_poly_a_enriched"
        and not proto.get("explicit")
    ):
        advisory = (
            f"Sampled R1 reads show {geom.get('r1_poly_a_frac', 0.0) * 100:.1f}% "
            "poly(A)/poly(T) stretches, consistent with a 3'-tag RNA-seq library "
            "(e.g. QuantSeq 3' or SMART-3Seq). Consider setting "
            "`rna_library_protocol='tag_3prime'` to enable poly(X) tail trimming "
            "and suppress length-biased TPM export."
        )

    return {
        "status": "blocked" if findings else "pass",
        "protocol": proto.get("protocol"),
        "explicit_protocol": proto.get("explicit"),
        "suppress_tpm": proto.get("suppress_tpm"),
        "umi_config": proto.get("umi_config"),
        "pooled_plate_manifest_summary": (
            {
                "preset": proto["pooled_plate_manifest"].get("preset"),
                "n_barcoded_samples": len(
                    proto["pooled_plate_manifest"]["barcode_to_sample"]
                ),
                "cdna_read": proto["pooled_plate_manifest"]["cdna_read"],
                "strand": proto["pooled_plate_manifest"]["strand"],
            }
            if proto.get("pooled_plate_manifest")
            else None
        ),
        "geometry": geom,
        "advisory": advisory,
        "findings": findings,
    }
