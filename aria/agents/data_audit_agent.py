"""
ARIA DataAuditAgent
-------------------
The gatekeeper. Always runs FIRST, before any analysis.

Responsibilities:
  1. Scan input directory, detect all data types automatically
  2. Infer organism, genome version, experimental design
  3. Validate completeness (pairs, replicates, barcodes)
  4. Trigger CHECKPOINT #1 — "This is what I found, confirm?"
  5. Homogenize to canonical formats
  6. Build ExperimentContext for all downstream agents

Detects:
  - scRNA-seq / bulk RNA-seq
  - scATAC-seq / bulk ATAC-seq
  - HiC / micro-C
  - ChIP-seq
  - CUT&RUN / CUT&TAG
  - Mixed multimodal experiments
"""

from __future__ import annotations
import os
import re
import uuid
import hashlib
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence, CavemanMode
from aria.memory.memory import ARIAMemory
from aria.utils.assay_detector import AssayDetector
from aria.utils.multiome_contracts import infer_multiome_contract
from aria.utils.provenance import hash_file


# ── File signature patterns ──────────────────────────────────────────────────

SIGNATURES = {
    # Single-cell RNA
    "scRNA": [
        r"barcodes\.tsv(\.gz)?$",
        r"features\.tsv(\.gz)?$",
        r"genes\.tsv(\.gz)?$",        # 10x v2 MEX component
        r"matrix\.mtx(\.gz)?$",
        r".*\.h5$",
        r".*\.h5ad$",
        r".*filtered_feature_bc_matrix.*",
        r".*cellranger.*",
    ],
    # Bulk RNA — counts matrix (ready for DESeq2)
    "bulk_RNA": [
        r".*counts?\.(txt|tsv|csv)$",
        r".*\.counts$",
        r".*featurecounts.*\.(txt|tsv)$",
        r".*htseq.*\.(txt|tsv)$",
        r".*salmon.*/quant\.sf$",
        r".*kallisto.*/abundance\.tsv$",
    ],
    # Bulk RNA — raw FASTQs (need preprocessing: fastp → STAR → featureCounts)
    "bulk_RNA_raw": [
        r".*_R[12]_.*\.fastq(\.gz)?$",
        r".*_R[12]\.fastq(\.gz)?$",
        r".*_[12]\.fastq(\.gz)?$",
        r".*_[12]\.fq(\.gz)?$",
    ],
    # Single-cell ATAC (incl. same-cell paired RNA+ATAC MuData, the v4.6 entry
    # path — C8, audit 2026-05-29). A `.h5mu` carries both modalities; it is
    # routed here because the chromatin/scATAC pipeline owns its ingestion.
    "scATAC": [
        r".*\.h5mu$",
        r"fragments\.tsv(\.gz)?$",
        r".*singlecell\.csv$",
        r".*atac.*barcodes.*",
        r".*scatac.*",
    ],
    # Bulk ATAC
    "bulk_ATAC": [
        r".*atac.*\.fastq(\.gz)?$",
        r".*atac.*\.bam$",
        r".*atac.*peaks?.*\.(bed|narrowPeak|broadPeak)$",
    ],
    # HiC / Micro-C
    "HiC": [
        r".*\.hic$",
        r".*\.cool$",
        r".*\.mcool$",
        r".*\.pairs(\.gz)?$",
        r".*hic.*\.fastq(\.gz)?$",
        r".*micro.?c.*",
    ],
    # ChIP-seq
    "ChIP": [
        r".*chip.*\.fastq(\.gz)?$",
        r".*chip.*\.bam$",
        r".*chip.*peaks?.*\.(bed|narrowPeak|broadPeak)$",
        r".*H[34]K\d+.*\.(fastq|bam)(\.gz)?$",  # histone marks
        r".*input.*\.bam$",
        r".*IgG.*\.bam$",
    ],
    # CUT&RUN
    "CUT_AND_RUN": [
        r".*cut.?and.?run.*\.(fastq|bam)(\.gz)?$",
        r".*cutnrun.*\.(fastq|bam)(\.gz)?$",
        r".*cut.?run.*\.(fastq|bam)(\.gz)?$",
    ],
    # CUT&TAG
    "CUT_AND_TAG": [
        r".*cut.?and.?tag.*\.(fastq|bam)(\.gz)?$",
        r".*cuttag.*\.(fastq|bam)(\.gz)?$",
        r".*cut.?tag.*\.(fastq|bam)(\.gz)?$",
    ],
}

GENOME_HINTS = {
    "hg38": ["hg38", "GRCh38", "human", "homo sapiens"],
    "hg19": ["hg19", "GRCh37"],
    "mm10": ["mm10", "GRCm38", "mouse"],
    "mm39": ["mm39", "GRCm39"],
    "dm6":  ["dm6", "drosophila"],
    "ce11": ["ce11", "worm", "elegans"],
    "danRer11": ["danRer11", "zebrafish"],
    "sacCer3": ["sacCer3", "yeast"],
}

ORGANISM_HINTS = {
    "Homo sapiens":        ["hg38", "hg19", "human", "sapiens"],
    "Mus musculus":        ["mm10", "mm39", "mouse", "musculus"],
    "Drosophila melanogaster": ["dm6", "drosophila"],
    "C. elegans":          ["ce11", "elegans", "worm"],
    "Danio rerio":         ["danRer11", "zebrafish"],
    "S. cerevisiae":       ["sacCer3", "yeast", "cerevisiae"],
}

# Modalities a user may pick when correcting the data audit at CHECKPOINT 1.
# Sourced from the recognized SIGNATURES so the menu always matches what ARIA
# can route (dispatch gating, e.g. scATAC/HiC, is enforced downstream).
SUPPORTED_MODALITIES = list(SIGNATURES.keys())

# Canonical reference genome per organism (the default an alignment-based
# pipeline uses); the first assembly in each organism's hint list.
_DEFAULT_GENOME = {
    "Homo sapiens": "hg38",
    "Mus musculus": "mm10",
    "Drosophila melanogaster": "dm6",
    "C. elegans": "ce11",
    "Danio rerio": "danRer11",
    "S. cerevisiae": "sacCer3",
}


_MEX_MATRIX_RE = re.compile(r"matrix\.mtx(\.gz)?$", re.IGNORECASE)
_MEX_COMPONENT_RE = re.compile(
    r"(matrix\.mtx|barcodes\.tsv|genes\.tsv|features\.tsv)(\.gz)?$", re.IGNORECASE)


_MEX_BARCODES_RE = re.compile(r"barcodes\.tsv(\.gz)?$", re.IGNORECASE)
_MEX_FEATURES_RE = re.compile(r"(genes|features)\.tsv(\.gz)?$", re.IGNORECASE)


def _mex_dir_components(paths: list) -> dict:
    """Map each directory that holds a ``matrix.mtx`` to the MEX components it
    contains: ``{dir: {"matrix", "barcodes", "features"}}`` (set membership).

    Pure helper over path strings so completeness can be reasoned about (and
    unit-tested) without touching disk. Shared by ``_collapse_mex_directories``
    (T5: collapse only complete triplets) and ``_incomplete_mex_warnings``.
    """
    from pathlib import Path

    dirs: dict = {}
    for p in (paths or []):
        pp = Path(str(p))
        name = pp.name
        if not _MEX_COMPONENT_RE.search(name):
            continue
        parent = str(pp.parent)
        comps = dirs.setdefault(parent, set())
        if _MEX_MATRIX_RE.search(name):
            comps.add("matrix")
        elif _MEX_BARCODES_RE.search(name):
            comps.add("barcodes")
        elif _MEX_FEATURES_RE.search(name):
            comps.add("features")
    # Only directories that actually contain a matrix.mtx are MEX candidates.
    return {d: c for d, c in dirs.items() if "matrix" in c}


def _collapse_mex_directories(paths: list) -> list:
    """E2E-1 + T5: collapse a COMPLETE 10x MEX triplet into one sample directory.

    A 10x MEX sample is a directory holding ``matrix.mtx`` plus ``barcodes.tsv``
    and ``genes.tsv``/``features.tsv`` (optionally gzipped). Classifying each
    component as its own scRNA "file" makes the per-sample QC path treat
    ``barcodes.tsv`` as a sample and abort (PerSampleQCFailed). Here, a directory
    that holds the COMPLETE triplet collapses its component files into a single
    entry = that directory, which ``rna_qc`` loads as one MEX sample.

    T5: a directory with ``matrix.mtx`` but MISSING ``barcodes`` or
    ``features``/``genes`` is NOT a loadable sample, so it is left untouched (its
    component files stay) and surfaced via ``_incomplete_mex_warnings`` at CP1,
    rather than collapsed and failing late in ``sc.read_10x_mtx``.

    Non-MEX inputs (``.h5ad``/``.h5``) and loose files are left untouched. Order
    is preserved and each collapsed MEX directory is emitted once.
    """
    from pathlib import Path

    paths = [str(p) for p in (paths or [])]
    components = _mex_dir_components(paths)
    complete_dirs = {
        d for d, comps in components.items()
        if {"matrix", "barcodes", "features"} <= comps
    }
    if not complete_dirs:
        return paths
    out: list = []
    seen: set = set()
    for p in paths:
        pp = Path(p)
        if str(pp.parent) in complete_dirs and _MEX_COMPONENT_RE.search(pp.name):
            if str(pp.parent) not in seen:
                out.append(str(pp.parent))
                seen.add(str(pp.parent))
            continue
        out.append(p)
    return out


def _incomplete_mex_warnings(paths: list) -> list:
    """T5: warn at CP1 for every directory that has a ``matrix.mtx`` but is
    missing the ``barcodes`` and/or ``features``/``genes`` component, naming what
    is missing. An incomplete MEX cannot be loaded as a sample, so this surfaces
    the problem early instead of failing late inside ``rna_qc``.
    """
    warnings: list = []
    for d in sorted(_mex_dir_components(paths)):
        comps = _mex_dir_components(paths)[d]
        missing = []
        if "barcodes" not in comps:
            missing.append("barcodes.tsv")
        if "features" not in comps:
            missing.append("features.tsv/genes.tsv")
        if not missing:
            continue
        warnings.append(
            f"Incomplete 10x MEX directory '{d}': found matrix.mtx but missing "
            f"{', '.join(missing)}. It cannot be loaded as a scRNA sample; "
            f"provide the missing component(s) or remove the directory."
        )
    return warnings


def default_genome_for_organism(organism: str | None) -> str | None:
    """Return the canonical reference assembly for an organism, or None."""
    return _DEFAULT_GENOME.get(str(organism or "").strip())


def _geo_bucket_map(geo_data_type: str) -> dict[str, str]:
    """Map GEOConnector file buckets -> ARIA modality for a given data_type.

    ATAC studies must reach the ChromatinAgent (scATAC/bulk_ATAC) rather than
    being collapsed into the RNA lane. Only buckets a lane can actually consume
    are mapped: ATAC count tables are not ingested directly (peaks/counts are
    recalled internally from BAM/fragments), so they are intentionally omitted.
    """
    if geo_data_type == "scATAC":
        return {"fragments": "scATAC", "peaks": "scATAC", "bam": "scATAC",
                "h5": "scATAC", "h5ad": "scATAC", "mtx": "scATAC",
                "fastq": "scATAC"}
    if geo_data_type == "bulk_ATAC":
        return {"bam": "bulk_ATAC", "peaks": "bulk_ATAC",
                "fastq": "bulk_ATAC"}
    rna_mod = "scRNA" if geo_data_type == "scRNA" else "bulk_RNA"
    mapping = {"counts": rna_mod, "h5ad": "scRNA", "h5": "scRNA", "mtx": "scRNA"}
    mapping["fastq"] = "scRNA" if geo_data_type == "scRNA" else "bulk_RNA_raw"
    return mapping


def apply_metadata_corrections(exp_context: dict, corrections: dict | None) -> dict:
    """Apply a user's CHECKPOINT-1 metadata corrections in place.

    ``corrections`` may carry ``modality`` (re-assign ALL audited files to that
    single modality — the user is asserting "this data is actually X"),
    ``organism``, ``genome``, and an edited ``scatac_fastq_manifest`` (inline or
    path). Returns the same ``exp_context``. A falsy ``corrections`` is a no-op.
    """
    if not corrections:
        return exp_context
    modality = corrections.get("modality")
    if modality:
        all_files = [
            f
            for files in (exp_context.get("modalities") or {}).values()
            for f in (files or [])
        ]
        exp_context["modalities"] = {modality: all_files}
    if corrections.get("organism"):
        exp_context["organism"] = corrections["organism"]
    if corrections.get("genome"):
        exp_context["genome"] = corrections["genome"]
    manifest_edit = (
        corrections.get("scatac_fastq_manifest")
        or corrections.get("scatac_fastq_manifest_path")
    )
    if manifest_edit:
        from aria.utils.scatac_fastq_manifest import (
            manifest_fastq_files,
            resolve_scatac_fastq_manifest,
        )
        candidate_context = {
            "data_dir": exp_context.get("data_dir"),
            "scatac_fastq_manifest": manifest_edit,
        }
        validation = resolve_scatac_fastq_manifest(
            candidate_context, require_paths=True
        )
        exp_context["scatac_fastq_manifest_validation"] = validation
        if validation.get("status") == "valid":
            manifest = validation["manifest"]
            exp_context["scatac_fastq_manifest"] = manifest
            manifest_files = set(manifest_fastq_files(manifest))
            modalities = exp_context.get("modalities") or {}
            updated = {}
            for key, paths in modalities.items():
                if key == "scATAC":
                    continue
                kept = [path for path in paths or [] if path not in manifest_files]
                if kept:
                    updated[key] = kept
            updated["scATAC"] = sorted(manifest_files)
            exp_context["modalities"] = updated
    return exp_context


@dataclass(frozen=True)
class DataAuditScanLimits:
    """Bounds for the initial directory walk.

    Size is optional because DataAudit only inspects file names here; very large
    BAM/H5 files can be valid inputs and should not be skipped by default.
    """

    max_files: int = 5000
    max_entries: int = 20000
    max_depth: int = 8
    max_seconds: float = 10.0
    max_file_size_bytes: Optional[int] = None
    follow_symlinks: bool = False

    @classmethod
    def from_config(cls, config: Optional[dict] = None) -> "DataAuditScanLimits":
        values = {
            "max_files": _env_int("ARIA_DATA_AUDIT_MAX_FILES", cls.max_files),
            "max_entries": _env_int("ARIA_DATA_AUDIT_MAX_ENTRIES", cls.max_entries),
            "max_depth": _env_int("ARIA_DATA_AUDIT_MAX_DEPTH", cls.max_depth),
            "max_seconds": _env_float("ARIA_DATA_AUDIT_MAX_SECONDS", cls.max_seconds),
            "max_file_size_bytes": _env_optional_int(
                "ARIA_DATA_AUDIT_MAX_FILE_SIZE_BYTES",
                cls.max_file_size_bytes,
            ),
            "follow_symlinks": _env_bool(
                "ARIA_DATA_AUDIT_FOLLOW_SYMLINKS",
                cls.follow_symlinks,
            ),
        }
        if config:
            values.update({k: v for k, v in config.items() if k in values})
        return cls(
            max_files=max(1, int(values["max_files"])),
            max_entries=max(1, int(values["max_entries"])),
            max_depth=max(0, int(values["max_depth"])),
            max_seconds=max(0.0, float(values["max_seconds"])),
            max_file_size_bytes=(
                None if values["max_file_size_bytes"] in (None, "")
                else max(0, int(values["max_file_size_bytes"]))
            ),
            follow_symlinks=bool(values["follow_symlinks"]),
        )

    def as_dict(self) -> dict:
        return {
            "max_files": self.max_files,
            "max_entries": self.max_entries,
            "max_depth": self.max_depth,
            "max_seconds": self.max_seconds,
            "max_file_size_bytes": self.max_file_size_bytes,
            "follow_symlinks": self.follow_symlinks,
        }


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw in (None, ""):
        return default
    try:
        return int(raw)
    except ValueError:
        return default


def _env_optional_int(name: str, default: Optional[int]) -> Optional[int]:
    raw = os.environ.get(name)
    if raw in (None, ""):
        return default
    try:
        value = int(raw)
    except ValueError:
        return default
    return None if value < 0 else value


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name)
    if raw in (None, ""):
        return default
    try:
        return float(raw)
    except ValueError:
        return default


def _env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw in (None, ""):
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


class DataAuditAgent(BaseAgent):

    name = "data_audit_agent"
    description = (
        "Scans input directory, detects all omics data types, "
        "validates experimental design, triggers Checkpoint 1."
    )

    def __init__(self, memory: ARIAMemory,
                 llm=None,
                 api_key: str = None):
        super().__init__(memory, llm=llm, api_key=api_key)

    def run(self, experiment_id: str, context: dict) -> dict:
        """
        Main audit pipeline.

        context must contain:
          - data_dir: path to the raw data directory
          - user_question: what the user wants to analyze (optional)
          - geo_metadata: (optional) dict from GEOConnector.fetch()
        """
        data_dir      = Path(context["data_dir"])
        user_question = context.get("user_question", "")
        geo_metadata  = context.get("geo_metadata")
        reproducible_mode = bool(context.get("reproducible_mode"))
        e6_analyzable_paths: set[str] | None = None

        # E6 connector results explicitly declare whether the accession was
        # fully validated and atomically published. Never audit a known partial
        # generation. Legacy GEO metadata without this field remains readable.
        retrieval_status = (
            geo_metadata.get("retrieval_status") if geo_metadata else None
        )
        if retrieval_status and retrieval_status != "complete":
            return {
                "status": "failed",
                "error": (
                    "GEO/SRA retrieval is incomplete and cannot enter data audit "
                    f"(status={retrieval_status!r})"
                ),
            }
        if retrieval_status == "complete":
            from aria.utils.atomic_retrieval import validate_retrieval_manifest

            manifest_path = Path(geo_metadata.get("retrieval_manifest") or "")
            manifest_root = manifest_path.parent.resolve()
            validation = (
                validate_retrieval_manifest(
                    manifest_root,
                    expected_accession=geo_metadata.get("accession"),
                )
                if manifest_path.name == "retrieval_manifest.json"
                else {"status": "invalid", "errors": ["retrieval manifest missing"]}
            )
            if validation.get("status") != "valid":
                return {
                    "status": "failed",
                    "error": (
                        "GEO/SRA retrieval manifest is invalid and cannot enter "
                        "data audit: "
                        + "; ".join(validation.get("errors") or ["unknown error"])
                    ),
                }
            if data_dir.resolve() != manifest_root:
                return {
                    "status": "failed",
                    "error": (
                        "GEO/SRA data directory is outside the validated retrieval "
                        f"manifest root: {data_dir}"
                    ),
                }
            manifested = {
                str(row.get("path") or "")
                for row in validation["manifest"].get("files", [])
            }
            declared_paths = [
                path
                for paths in (geo_metadata.get("files") or {}).values()
                if isinstance(paths, list)
                for path in paths
            ]
            declared_paths.extend(
                (geo_metadata.get("file_modalities") or {}).keys()
            )
            e6_analyzable_paths = {
                str(Path(raw_path).resolve()) for raw_path in declared_paths
            }
            for raw_path in declared_paths:
                path = Path(raw_path).resolve()
                try:
                    relative = path.relative_to(manifest_root).as_posix()
                except ValueError:
                    return {
                        "status": "failed",
                        "error": (
                            "GEO/SRA payload is outside the validated retrieval "
                            f"manifest: {path}"
                        ),
                    }
                if relative not in manifested:
                    return {
                        "status": "failed",
                        "error": (
                            "GEO/SRA payload is absent from the validated retrieval "
                            f"manifest: {relative}"
                        ),
                    }

        # E5: an explicit 10x/scATAC manifest is the authoritative source for
        # library identity and R1/R2/R3 roles. It may be inline, referenced by
        # path, or discovered as scatac_fastq_manifest.json in data_dir, so the
        # same contract works from TUI and headless entrypoints.
        from aria.utils.scatac_fastq_manifest import (
            manifest_library_types,
            resolve_scatac_fastq_manifest,
        )
        scatac_manifest = resolve_scatac_fastq_manifest(
            context, require_paths=True
        )

        self.publish_status(experiment_id,
            f"Scanning {data_dir}...", progress=0.0)

        # 1. Scan files
        all_files = self._scan_directory(
            data_dir,
            limits=context.get("scan_limits"),
        )
        if not all_files and not geo_metadata:
            return {
                "status": "failed",
                "error": f"No files found in {data_dir}"
            }

        # The manifest JSON and barcode whitelists are assay metadata, not
        # modality payloads. Keep them in the audit scan/provenance but do not
        # create a spurious ``unknown`` modality hall for them.
        assay_metadata_paths: set[str] = set()
        if scatac_manifest.get("status") == "valid":
            assay_metadata_paths.update(
                str(row["barcode_whitelist"])
                for row in scatac_manifest["manifest"]["libraries"]
            )
        if scatac_manifest.get("source") not in (None, "inline"):
            assay_metadata_paths.add(str(scatac_manifest["source"]))
        modality_files = [
            path for path in all_files if str(path) not in assay_metadata_paths
        ]
        if e6_analyzable_paths is not None:
            modality_files = [
                path for path in modality_files
                if str(path.resolve()) in e6_analyzable_paths
            ]

        # 2. Classify files by modality. E2: an explicit declared library type
        # (assay manifest / library_type) is authoritative; filenames are a hint.
        per_file_types = dict(context.get("library_types") or {})
        if scatac_manifest.get("status") == "valid":
            per_file_types.update(
                manifest_library_types(scatac_manifest["manifest"])
            )
        declared_library_types = {
            "global": context.get("library_type"),
            "per_file": per_file_types,
        }
        classified = self._classify_files(modality_files, declared_library_types)
        if scatac_manifest.get("status") == "invalid":
            # An explicit-but-invalid scATAC manifest signals assay intent but
            # must never let generic R1/R2 names fall through to bulk RNA. Route
            # all scanned FASTQs to the scATAC readiness card, which then blocks
            # with the manifest errors and remains editable at CP1.
            fastq_paths = {
                str(path) for path in modality_files
                if path.name.lower().endswith(
                    (".fastq", ".fastq.gz", ".fq", ".fq.gz")
                )
            }
            for modality, paths in list(classified.items()):
                kept = [path for path in paths if path not in fastq_paths]
                if kept:
                    classified[modality] = kept
                else:
                    classified.pop(modality, None)
            classified["scATAC"] = sorted(fastq_paths)

        # 2b. When GEO metadata is present, enrich classification with
        #     already-typed files and remove them from "unknown".
        if geo_metadata:
            geo_files     = geo_metadata.get("files", {})
            geo_data_type = geo_metadata.get("data_type", "bulk_RNA")

            # Map each GEO file bucket to an ARIA modality, keyed by the
            # connector's inferred data_type, so ATAC studies reach the
            # ChromatinAgent instead of being collapsed into the RNA lane.
            bucket_map = _geo_bucket_map(geo_data_type)

            for ftype, bucket in bucket_map.items():
                for fpath in geo_files.get(ftype, []):
                    if Path(fpath).exists():
                        bucket_paths = classified.setdefault(bucket, [])
                        if fpath not in bucket_paths:
                            bucket_paths.append(fpath)
                        unknown = classified.get("unknown", [])
                        if fpath in unknown:
                            unknown.remove(fpath)

            # Mixed SRA studies carry a per-FASTQ authoritative modality map.
            # Remove each declared path from every filename-derived bucket before
            # assigning it, otherwise generic R1/R2 hints collapse ATAC into RNA.
            for fpath, bucket in (geo_metadata.get("file_modalities") or {}).items():
                if not Path(fpath).exists():
                    continue
                for existing_bucket, paths in list(classified.items()):
                    classified[existing_bucket] = [
                        path for path in paths if str(path) != str(fpath)
                    ]
                    if not classified[existing_bucket]:
                        classified.pop(existing_bucket, None)
                classified.setdefault(bucket, []).append(str(fpath))

        classified, ignored_intermediates = self._filter_aria_intermediate_outputs(
            classified
        )

        # 2c. Preprocessed h5ad files often already contain the experimental
        # design in obs. Inspect it before DesignAgent so user checkpoints are
        # seeded from data, not filename guesses.
        h5ad_design = self._infer_h5ad_design(
            classified.get("scRNA", []), user_question
        )

        # 3. Infer organism and genome (GEO metadata takes precedence)
        if geo_metadata:
            inferred = geo_metadata.get("inferred_design", {})
            organism = inferred.get("organism", "") or geo_metadata.get("organism", "unknown")
            genome   = inferred.get("genome",   "") or geo_metadata.get("genome",   "unknown")
        else:
            genome, organism = self._infer_genome_organism(
                all_files, data_dir, user_question
            )
            if genome == "unknown" or organism == "unknown":
                h5ad_genome, h5ad_organism = self._infer_h5ad_genome_organism(
                    classified.get("scRNA", [])
                )
                if genome == "unknown" and h5ad_genome != "unknown":
                    genome = h5ad_genome
                if organism == "unknown" and h5ad_organism != "unknown":
                    organism = h5ad_organism

        # 4. Validate design (replicates, pairs, etc.)
        warnings = self._validate_design(classified)
        warnings.extend(self._scan_report_warnings())
        warnings.extend(self._assay_detection_warnings())
        warnings.extend(self._ambiguous_library_type_warnings())
        warnings.extend(getattr(self, "_mex_warnings", []))
        if ignored_intermediates:
            warnings.append(
                "Ignored ARIA intermediate h5ad output(s) during data audit: "
                + ", ".join(Path(p).name for p in ignored_intermediates[:5])
                + ("..." if len(ignored_intermediates) > 5 else "")
            )
        warnings.extend(h5ad_design.get("warnings", []))

        # 5. Build ExperimentContext
        exp_context = self._build_context(
            experiment_id, data_dir, classified,
            genome, organism, warnings, user_question
        )
        exp_context["reproducible_mode"] = reproducible_mode
        exp_context["scatac_fastq_manifest_validation"] = scatac_manifest
        if scatac_manifest.get("status") == "valid":
            exp_context["scatac_fastq_manifest"] = scatac_manifest["manifest"]
            if scatac_manifest.get("source") != "inline":
                exp_context["scatac_fastq_manifest_path"] = scatac_manifest["source"]
        # ADR-057: opt-in SPECULATIVE hypotheses section. Mirrors the
        # reproducible_mode plumbing — surfaced explicitly via the entrypoints
        # (run_headless param / TUI --hypotheses) with an
        # ARIA_ENABLE_HYPOTHESES env fallback for power users. Off by default,
        # never inferred from the data or question.
        env_hyp = os.environ.get("ARIA_ENABLE_HYPOTHESES", "").strip().lower()
        exp_context["enable_hypotheses"] = (
            bool(context.get("enable_hypotheses"))
            or env_hyp in ("1", "true", "yes", "on")
        )
        if context.get("provenance"):
            exp_context["provenance"] = context["provenance"]

        # Propagate GEO metadata into exp_context so downstream agents can use it
        if geo_metadata:
            exp_context["geo_metadata"]    = geo_metadata
            exp_context["inferred_design"] = geo_metadata.get("inferred_design", {})
        elif h5ad_design.get("groups"):
            exp_context["inferred_design"] = h5ad_design

        # 6. Store in memory
        self.memory.create_wing(
            experiment_id,
            name=data_dir.name,
            organism=organism,
            genome=genome
        )

        for modality in classified:
            hall_id = f"{experiment_id}_{modality}"
            self.memory.create_hall(hall_id, experiment_id, modality)

        self.publish_status(experiment_id,
            "Audit complete. Preparing checkpoint.", progress=0.9)

        # 7. CHECKPOINT #1 — show the user what was found
        checkpoint_msg = self._build_checkpoint_summary(
            classified, genome, organism, warnings,
            inferred_design=exp_context.get("inferred_design", {}),
        )

        # P1-8a / W-PRIV: classify input sensitivity and surface it at CP1 so the
        # user can opt into air-gapped mode (block ALL egress). ARIA never flips
        # air-gapped on by itself — it classifies, recommends, and offers a choice.
        from aria.utils.sensitivity import (
            annotate_checkpoint_question,
            checkpoint_options,
            classify_sensitivity,
        )
        sensitivity = classify_sensitivity(
            organism=organism,
            field_names=self._collect_sensitivity_fields(exp_context),
            path_hints=[data_dir.name] + [p.name for p in all_files[:50]],
        )
        exp_context["sensitivity"] = sensitivity

        # 'Confirm and continue' stays first so the default (incl. headless) is
        # unchanged — ARIA never flips air-gapped on without the user's explicit
        # choice; the air-gapped option is always available and recommended in the
        # text when the input looks sensitive.
        checkpoint_msg = annotate_checkpoint_question(checkpoint_msg, sensitivity)
        options = checkpoint_options(sensitivity)

        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=1,
            question=checkpoint_msg,
            options=options,
            context={
                "classified":  classified,
                "genome":      genome,
                "organism":    organism,
                "warnings":    warnings,
                "file_count":  len(all_files),
                "sensitivity": sensitivity,
                "exp_context": exp_context,
            }
        )

        return {
            "status":      "awaiting_checkpoint",
            "checkpoint":  1,
            "exp_context": exp_context
        }

    # ── PRIVATE METHODS ──────────────────────────────────────────────────

    @staticmethod
    def _collect_sensitivity_fields(exp_context: dict) -> list[str]:
        """Gather obs/design field names that feed the sensitivity classifier
        (P1-8a). Uses metadata column NAMES only — never cell-level values."""
        design = exp_context.get("inferred_design", {}) or {}
        fields: list[str] = []
        obs_cols = design.get("obs_columns") or {}
        if isinstance(obs_cols, dict):
            for cols in obs_cols.values():
                fields.extend(str(c) for c in (cols or []))
        elif isinstance(obs_cols, list):
            fields.extend(str(c) for c in obs_cols)
        for key in ("condition_col", "replicate_col", "groupby_col", "main_factor"):
            if design.get(key):
                fields.append(str(design[key]))
        covs = design.get("covariates") or []
        if isinstance(covs, (list, tuple)):
            fields.extend(str(c) for c in covs)
        return fields

    @staticmethod
    def _supported_input_file(path: Path) -> bool:
        extensions = {
            ".fastq", ".gz", ".bam", ".bai", ".sam",
            ".bed", ".narrowPeak", ".broadPeak", ".bigWig", ".bw",
            ".hic", ".cool", ".mcool", ".pairs",
            ".h5", ".h5ad", ".h5mu", ".loom",
            ".mtx", ".tsv", ".csv", ".txt",
        }
        return (
            path.suffix in extensions
            or "".join(path.suffixes) in {".fastq.gz", ".pairs.gz",
                                          ".tsv.gz", ".bed.gz"}
            or AssayDetector().is_supported_file(path)
        )

    def _scan_directory(
        self,
        data_dir: Path,
        limits: Optional[dict | DataAuditScanLimits] = None,
    ) -> list[Path]:
        """Recursively scan directory for input files with explicit bounds."""
        scan_limits = (
            limits if isinstance(limits, DataAuditScanLimits)
            else DataAuditScanLimits.from_config(limits)
        )
        files = []
        report = {
            "limits": scan_limits.as_dict(),
            "entries_seen": 0,
            "files_seen": 0,
            "matched_files": 0,
            "skipped_symlinks": 0,
            "skipped_large_files": 0,
            "skipped_deep_dirs": 0,
            "errors": [],
            "truncated": False,
            "truncated_reason": None,
        }
        self._last_scan_report = report

        start = time.monotonic()
        deadline = start + scan_limits.max_seconds
        stack: list[tuple[Path, int]] = [(data_dir, 0)]
        visited_dirs: set[tuple[int, int]] = set()

        def stop(reason: str) -> bool:
            report["truncated"] = True
            report["truncated_reason"] = reason
            return True

        while stack:
            if time.monotonic() >= deadline:
                stop("scan_timeout")
                break
            current, depth = stack.pop()
            try:
                stat = current.stat()
                marker = (stat.st_dev, stat.st_ino)
                if marker in visited_dirs:
                    continue
                visited_dirs.add(marker)
            except OSError as e:
                report["errors"].append(f"{current}: {e}")
                continue

            try:
                with os.scandir(current) as entries:
                    for entry in entries:
                        if time.monotonic() >= deadline:
                            stop("scan_timeout")
                            break
                        report["entries_seen"] += 1
                        if report["entries_seen"] > scan_limits.max_entries:
                            stop("max_entries")
                            break

                        entry_path = Path(entry.path)
                        try:
                            is_link = entry.is_symlink()
                        except OSError as e:
                            report["errors"].append(f"{entry_path}: {e}")
                            continue
                        if is_link and not scan_limits.follow_symlinks:
                            report["skipped_symlinks"] += 1
                            continue

                        try:
                            if entry.is_dir(follow_symlinks=scan_limits.follow_symlinks):
                                if depth >= scan_limits.max_depth:
                                    report["skipped_deep_dirs"] += 1
                                else:
                                    stack.append((entry_path, depth + 1))
                                continue
                            if not entry.is_file(follow_symlinks=scan_limits.follow_symlinks):
                                continue
                            size = entry.stat(
                                follow_symlinks=scan_limits.follow_symlinks
                            ).st_size
                        except OSError as e:
                            report["errors"].append(f"{entry_path}: {e}")
                            continue

                        report["files_seen"] += 1
                        max_size = scan_limits.max_file_size_bytes
                        if max_size is not None and size > max_size:
                            report["skipped_large_files"] += 1
                            continue
                        if self._supported_input_file(entry_path):
                            files.append(entry_path)
                            report["matched_files"] = len(files)
                            if len(files) >= scan_limits.max_files:
                                stop("max_files")
                                break
                    if report["truncated"]:
                        break
            except PermissionError as e:
                report["errors"].append(f"{current}: {e}")
                self.publish_status("", f"Permission error: {e}")
            except OSError as e:
                report["errors"].append(f"{current}: {e}")
        return files

    def _scan_report_warnings(self) -> list[str]:
        report = getattr(self, "_last_scan_report", None)
        if not report:
            return []
        warnings = []
        if report.get("truncated"):
            limits = report.get("limits", {})
            warnings.append(
                "Data audit directory scan was truncated "
                f"({report.get('truncated_reason')}; "
                f"files={report.get('matched_files')}, "
                f"entries_seen={report.get('entries_seen')}, "
                f"limits={limits})."
            )
        if report.get("skipped_symlinks"):
            warnings.append(
                "Data audit skipped symlinked paths "
                f"({report['skipped_symlinks']}); set "
                "ARIA_DATA_AUDIT_FOLLOW_SYMLINKS=1 or scan_limits.follow_symlinks "
                "to opt in."
            )
        if report.get("skipped_large_files"):
            warnings.append(
                "Data audit skipped files over the configured size limit "
                f"({report['skipped_large_files']})."
            )
        if report.get("skipped_deep_dirs"):
            warnings.append(
                "Data audit skipped directories deeper than the configured limit "
                f"({report['skipped_deep_dirs']})."
            )
        if report.get("errors"):
            warnings.append(
                "Data audit encountered filesystem errors while scanning: "
                + "; ".join(report["errors"][:3])
                + ("..." if len(report["errors"]) > 3 else "")
            )
        return warnings

    _FASTQ_NAME_RE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)

    @classmethod
    def _fastq_modality_hint(cls, fname: str, fpath: str) -> tuple[str | None, bool]:
        """E2: filename modality HINT for a FASTQ (never binding).

        A modality-keyword signature (atac/scatac/hic/chip/cut&run/cut&tag, histone
        marks) beats the generic paired-end ``bulk_RNA_raw`` rule — fixing the
        order-dependent misroute. Returns ``(modality, ambiguous)``: a keyword hit
        is a confident hint (``ambiguous=False``); a plain ``_R[12]`` FASTQ with no
        modality signal falls to ``bulk_RNA_raw`` as an ambiguous hint
        (``ambiguous=True``); no FASTQ signature match → ``(None, False)``.
        """
        for modality, patterns in SIGNATURES.items():
            if modality == "bulk_RNA_raw":
                continue  # generic paired-end rule is the last resort for FASTQ
            for pat in patterns:
                if re.search(pat, fname, re.IGNORECASE) or \
                   re.search(pat, fpath, re.IGNORECASE):
                    return modality, False
        for pat in SIGNATURES["bulk_RNA_raw"]:
            if re.search(pat, fname, re.IGNORECASE) or \
               re.search(pat, fpath, re.IGNORECASE):
                return "bulk_RNA_raw", True
        return None, False

    @staticmethod
    def _declared_library_type(path: Path, declared: dict | None) -> str | None:
        """E2: an explicit declared library type is authoritative over filenames.

        ``declared`` may carry ``per_file`` (basename or full-path → modality) and a
        ``global`` modality applied to every otherwise-undetected file.
        """
        if not declared:
            return None
        per_file = declared.get("per_file") or {}
        for key in (str(path), path.name):
            if key in per_file and per_file[key]:
                return str(per_file[key])
        global_type = declared.get("global")
        return str(global_type) if global_type else None

    def _classify_files(
        self,
        files: list[Path],
        declared_library_types: dict | None = None,
    ) -> dict[str, list[str]]:
        """Map files to their omics modality."""
        classified: dict[str, list[str]] = {}
        detector = AssayDetector()
        self._last_assay_detections = []
        self._ambiguous_library_types = []

        for f in files:
            fname = f.name.lower()
            fpath = str(f).lower()

            detection = detector.detect_file(f)
            if detection is not None:
                classified.setdefault(detection.modality, []).append(str(f))
                self._record_assay_detection(f, detection)
                continue

            # E2: an explicit declared library type is authoritative.
            declared = self._declared_library_type(f, declared_library_types)
            if declared:
                classified.setdefault(declared, []).append(str(f))
                continue

            matched = False
            if self._FASTQ_NAME_RE.search(fname):
                # FASTQ: filenames are only a hint; keyword modality beats the
                # generic paired-end rule, and a signalless R1/R2 is ambiguous.
                modality, ambiguous = self._fastq_modality_hint(fname, fpath)
                if modality is not None:
                    classified.setdefault(modality, []).append(str(f))
                    if ambiguous:
                        self._ambiguous_library_types.append(str(f))
                    matched = True
            else:
                for modality, patterns in SIGNATURES.items():
                    for pat in patterns:
                        if re.search(pat, fname, re.IGNORECASE) or \
                           re.search(pat, fpath, re.IGNORECASE):
                            classified.setdefault(modality, []).append(str(f))
                            matched = True
                            break
                    if matched:
                        break

            if not matched:
                classified.setdefault("unknown", []).append(str(f))

        # E2E-1: a 10x MEX triplet is ONE sample (its directory), not 3 files.
        # T5: an incomplete MEX (matrix.mtx without barcodes/features) is NOT
        # collapsed and is reported at CP1 with the missing component named.
        if classified.get("scRNA"):
            classified = self._demote_non_tenx_mtx_sidecars(classified)
            if classified.get("scRNA"):
                self._mex_warnings = _incomplete_mex_warnings(classified["scRNA"])
                classified["scRNA"] = _collapse_mex_directories(classified["scRNA"])
            else:
                self._mex_warnings = []
        else:
            self._mex_warnings = []

        return classified

    def _demote_non_tenx_mtx_sidecars(
        self, classified: dict[str, list[str]]
    ) -> dict[str, list[str]]:
        records = getattr(self, "_last_assay_detections", []) or []
        ambiguous_dirs = {
            str(Path(rec.get("path", "")).parent)
            for rec in records
            if "tenx_cell_barcode_evidence_missing"
            in set(rec.get("blocking_issues") or [])
        }
        if not ambiguous_dirs:
            return classified

        out = {key: list(value) for key, value in classified.items()}
        kept_scrna = []
        demoted = []
        for path in out.get("scRNA", []):
            p = Path(path)
            if str(p.parent) in ambiguous_dirs and _MEX_COMPONENT_RE.search(p.name):
                demoted.append(path)
            else:
                kept_scrna.append(path)
        if demoted:
            out["scRNA"] = kept_scrna
            out.setdefault("unknown", []).extend(demoted)
        if not out.get("scRNA"):
            out.pop("scRNA", None)
        return out

    def _record_assay_detection(self, path: Path, detection) -> None:
        records = getattr(self, "_last_assay_detections", None)
        if records is None:
            records = []
            self._last_assay_detections = records
        records.append({
            "path": str(path),
            "modality": detection.modality,
            "confidence": detection.confidence,
            "reason": detection.reason,
            "evidence": detection.evidence,
            "possible_alternatives": list(detection.possible_alternatives),
            "blocking_issues": list(detection.blocking_issues),
        })

    def _assay_detection_warnings(self) -> list[str]:
        records = getattr(self, "_last_assay_detections", []) or []
        warnings = []
        for rec in records:
            if rec.get("confidence") != "low" and not rec.get("blocking_issues"):
                continue
            issues = rec.get("blocking_issues") or []
            alt = rec.get("possible_alternatives") or []
            parts = [
                f"AssayDetector classified {Path(rec['path']).name} as "
                f"{rec.get('modality')} with {rec.get('confidence')} confidence"
            ]
            if alt:
                parts.append("alternatives: " + ", ".join(map(str, alt)))
            if issues:
                parts.append("issues: " + ", ".join(map(str, issues)))
            warnings.append("; ".join(parts) + ".")
        return warnings

    def _ambiguous_library_type_warnings(self) -> list[str]:
        """E2: surface generic R1/R2 FASTQ whose library type is only a hint, so
        CHECKPOINT 1 confirms the modality instead of a silent bulk-RNA dispatch."""
        ambiguous = getattr(self, "_ambiguous_library_types", []) or []
        if not ambiguous:
            return []
        names = ", ".join(Path(p).name for p in ambiguous[:5])
        if len(ambiguous) > 5:
            names += ", ..."
        return [
            "Ambiguous library type: generic paired-end FASTQ with no modality "
            f"signal ({names}) — hinted as bulk_RNA_raw but NOT bound. Declare the "
            "library type (assay manifest / library_type) or confirm the modality "
            "at CHECKPOINT 1; filenames are only a hint."
        ]

    @staticmethod
    def _is_aria_generated_output(path: str) -> bool:
        p = Path(path)
        lower_name = p.name.lower()
        if p.parent.name == "pseudobulk":
            return True
        if lower_name in {
            "de_per_cluster.csv",
            "cellcomm_interactions.csv",
            "pathways_per_cluster.csv",
            "pseudobulk_de.csv",
            "differential_abundance.tsv",
        }:
            return True
        if lower_name.endswith(".summary.json"):
            return True
        if p.suffix.lower() != ".h5ad":
            return False
        stem = p.stem
        if stem.startswith("qc_filtered_"):
            return True
        return stem in {
            "qc_filtered",
            "concatenated",
            "integrated",
            "annotated",
            "annotated_marker",
            "clustered",
            "clustered_sketch",
            "trajectory",
            "with_condition",
        }

    @staticmethod
    def _is_aria_intermediate_h5ad(path: str) -> bool:
        return DataAuditAgent._is_aria_generated_output(path)

    def _filter_aria_intermediate_outputs(
        self, classified: dict[str, list[str]]
    ) -> tuple[dict[str, list[str]], list[str]]:
        """
        Remove ARIA-generated intermediates from modality inputs. This prevents
        raw-data directories from being misread because a previous failed or
        partial run left `qc_filtered.h5ad`, `clustered.h5ad`, summary JSONs,
        DE CSVs, or a `pseudobulk/` folder next to the real input.
        """
        ignored: list[str] = []
        filtered: dict[str, list[str]] = {}
        for modality, files in classified.items():
            kept = []
            for f in files:
                if self._is_aria_generated_output(f):
                    ignored.append(f)
                else:
                    kept.append(f)
            if kept:
                filtered[modality] = kept

        return filtered, ignored

    def _infer_genome_organism(self, files: list[Path],
                                data_dir: Path,
                                user_question: str) -> tuple[str, str]:
        """Infer genome and organism from filenames, paths, and question."""
        search_text = " ".join([
            str(data_dir),
            user_question,
            " ".join(f.name for f in files[:50])
        ]).lower()

        genome = "unknown"
        organism = "unknown"

        for g, hints in GENOME_HINTS.items():
            if any(h.lower() in search_text for h in hints):
                genome = g
                break

        for org, hints in ORGANISM_HINTS.items():
            if any(h.lower() in search_text for h in hints):
                organism = org
                break

        # If genome found but not organism, infer organism from genome
        if genome != "unknown" and organism == "unknown":
            genome_to_org = {
                "hg38": "Homo sapiens", "hg19": "Homo sapiens",
                "mm10": "Mus musculus", "mm39": "Mus musculus",
                "dm6":  "Drosophila melanogaster",
                "ce11": "C. elegans",
                "danRer11": "Danio rerio",
                "sacCer3": "S. cerevisiae",
            }
            organism = genome_to_org.get(genome, "unknown")

        return genome, organism

    @staticmethod
    def _infer_h5ad_genome_organism(files: list[str]) -> tuple[str, str]:
        """Infer organism from h5ad feature names when metadata is absent."""
        h5ads = [f for f in files if str(f).lower().endswith(".h5ad")]
        if not h5ads:
            return "unknown", "unknown"

        try:
            import anndata as ad
        except ImportError:
            return "unknown", "unknown"

        # ADR-011 exception: this compact feature-name set is used only for
        # technical organism inference when metadata is absent, not for runtime
        # biological claims or cell-type assignment.
        human_markers = {
            "SAMD11", "ISG15", "TMEM88B", "PRDM16", "MEGF6", "C1QA",
            "C1QB", "C1QC", "NCMAP", "C1orf141", "TNFRSF1B", "PIK3CD",
        }
        mouse_patterns = ("Gm", "Rik")
        for path in h5ads[:2]:
            try:
                adata = ad.read_h5ad(path, backed="r")
                genes = [str(g) for g in list(adata.var_names[:3000])]
                backing_file = getattr(adata, "file", None)
                if backing_file is not None:
                    backing_file.close()
            except Exception:
                continue

            human_score = sum(g in human_markers for g in genes)
            human_score += sum(
                1 for g in genes
                if re.match(r"^C\d+orf\d+", g) or re.match(r"^LINC\d+", g)
            )
            mouse_score = sum(
                1 for g in genes
                if any(g.startswith(prefix) for prefix in mouse_patterns)
            )
            mouse_score += sum(1 for g in genes if re.match(r"^[A-Z][a-z]+$", g))

            if human_score >= 5 and mouse_score == 0:
                return "hg38", "Homo sapiens"

        return "unknown", "unknown"

    def _validate_design(self, classified: dict[str, list]) -> list[str]:
        """Check for common experimental design issues."""
        warnings = []

        for modality, files in classified.items():
            if modality == "unknown":
                continue

            # Check for paired-end completeness
            r1_files = [f for f in files if "_R1_" in f or "_1.fastq" in f]
            r2_files = [f for f in files if "_R2_" in f or "_2.fastq" in f]
            if r1_files and len(r1_files) != len(r2_files):
                warnings.append(
                    f"{modality}: Unequal R1/R2 pairs "
                    f"(R1={len(r1_files)}, R2={len(r2_files)})"
                )

            # Check minimum replicates for differential analysis
            if modality in ("bulk_RNA", "bulk_RNA_raw", "bulk_ATAC", "ChIP") and len(files) < 4:
                warnings.append(
                    f"{modality}: Only {len(files)} files detected. "
                    f"Differential analysis needs ≥2 replicates per condition."
                )

        return warnings

    def _infer_h5ad_design(self, files: list[str],
                           user_question: str = "") -> dict:
        """
        Inspect .h5ad obs metadata and infer design hints for scRNA pseudobulk.

        This is deliberately conservative: it only proposes a design when it can
        identify a condition column and a biological replicate column with at
        least two replicates per condition. The user still confirms the design in
        DesignAgent checkpoints.
        """
        h5ads = [f for f in files if str(f).lower().endswith(".h5ad")]
        if not h5ads:
            return {}

        warnings = []
        try:
            import anndata as ad
            import pandas as pd
        except ImportError:
            return {
                "warnings": [
                    "h5ad design inference skipped: anndata/pandas not available."
                ]
            }

        frames = []
        obs_columns = {}
        n_cells = 0
        for path in h5ads[:3]:
            try:
                adata = ad.read_h5ad(path, backed="r")
                obs = adata.obs.copy()
                backing_file = getattr(adata, "file", None)
                if backing_file is not None:
                    backing_file.close()
            except Exception as e:
                warnings.append(f"h5ad obs inspection failed for {Path(path).name}: {e}")
                continue
            if obs.empty:
                continue
            obs["_aria_source_file"] = Path(path).stem
            frames.append(obs)
            n_cells += int(obs.shape[0])
            obs_columns[Path(path).name] = list(map(str, obs.columns))

        if not frames:
            return {"warnings": warnings}

        obs_all = pd.concat(frames, axis=0, join="outer", sort=False)
        condition_col = self._pick_h5ad_condition_col(obs_all, user_question)
        # A replicate column must distinguish biological replicates within a
        # condition. The condition column itself trivially passes the laxer
        # replicate check (n=2 levels), and so does any duplicate of it
        # (e.g. Seurat's `orig.ident` often equals the condition label in
        # published h5ads). Forbid both so we never report the condition
        # column back as its own replicate.
        replicate_col = self._pick_h5ad_replicate_col(
            obs_all, exclude=[condition_col] if condition_col else None
        )
        groupby_col = self._pick_h5ad_groupby_col(obs_all)
        covariates = self._pick_h5ad_covariates(obs_all, condition_col, replicate_col)

        if not condition_col or not replicate_col:
            return {
                "source": "h5ad_obs",
                "obs_columns": obs_columns,
                "n_cells_inspected": n_cells,
                "warnings": warnings + [
                    "h5ad obs inspected but no reliable condition + replicate "
                    "design could be inferred."
                ],
            }

        obs_design = obs_all[[condition_col, replicate_col]].dropna().copy()
        obs_design[condition_col] = obs_design[condition_col].astype(str)
        obs_design[replicate_col] = obs_design[replicate_col].astype(str)

        groups = {}
        for level, sub in obs_design.groupby(condition_col):
            reps = sorted(r for r in sub[replicate_col].unique() if r and r != "nan")
            if reps:
                groups[str(level)] = reps

        groups = {g: reps for g, reps in groups.items() if len(reps) >= 1}
        if len(groups) < 2:
            return {
                "source": "h5ad_obs",
                "obs_columns": obs_columns,
                "n_cells_inspected": n_cells,
                "warnings": warnings + [
                    f"h5ad obs condition column '{condition_col}' has <2 usable levels."
                ],
            }

        replicate_counts = {g: len(reps) for g, reps in groups.items()}
        if not all(n >= 2 for n in replicate_counts.values()):
            warnings.append(
                "h5ad obs design has condition levels with <2 biological "
                f"replicates: {replicate_counts}. Pseudobulk may skip those contrasts."
            )

        comparisons = self._default_comparisons_from_groups(groups)
        return {
            "source": "h5ad_obs",
            "groups": groups,
            "main_factor": condition_col,
            "condition_col": condition_col,
            "replicate_col": replicate_col,
            "groupby_col": groupby_col,
            "covariates": covariates,
            "comparisons": comparisons,
            "pseudobulk": {
                "from_obs": True,
                "condition_col": condition_col,
                "replicate_col": replicate_col,
                "groupby_col": groupby_col,
                "covariates": covariates,
                "comparisons": comparisons,
            },
            "confidence": "high" if all(n >= 2 for n in replicate_counts.values()) else "medium",
            "reasoning": (
                f"Inferred from h5ad obs: condition='{condition_col}', "
                f"replicate='{replicate_col}', groupby='{groupby_col}'."
            ),
            "obs_columns": obs_columns,
            "n_cells_inspected": n_cells,
            "warnings": warnings,
        }

    @staticmethod
    def _pick_h5ad_condition_col(obs, user_question: str = "") -> str | None:
        # Generic experimental-design vocabulary. NOT dataset-specific terms
        # (no gene names, no disease names, no perturbation names): these are
        # patterns that appear across stimulation, perturbation, treatment,
        # genotype, disease, and age-group studies.
        priority = [
            "age_group", "age_bin", "condition", "treatment", "stim",
            "stimulation", "stim_status", "state", "perturbation",
            "perturb", "intervention", "diagnosis", "disease", "group",
            "genotype", "phenotype", "status", "label", "experimental_group",
        ]
        sustring_keys = (
            "age", "condition", "treatment", "stim", "perturb",
            "intervention", "state", "diagnos", "disease", "group",
            "genotype", "phenotype",
        )
        q = (user_question or "").lower()
        # User question can bias the priority order. Keep the bias generic:
        # we look at the verb the question uses ("vs", "compare", "between",
        # "stimulation", "treatment", "perturbation", "aging") and bring the
        # matching family of column names to the front. We do NOT inject
        # disease- or molecule-specific tokens (no "interferon", no "APOE").
        if any(k in q for k in ("aging", "age group", " age ", "young", "old")):
            priority = ["age_group", "age_bin", "age", *priority]
        if any(k in q for k in ("stimulat", "stim", "vehicle", "perturb",
                                 "agonist", "antagonist", "ligand")):
            priority = ["stim", "stimulation", "stim_status", "state",
                        "perturbation", "perturb", *priority]
        if any(k in q for k in ("treat", "treatment", "drug", "dose",
                                 "compound", "vehicle vs")):
            priority = ["treatment", "condition", *priority]
        # Build a case-insensitive lookup once.
        col_by_lower = {str(c).lower(): str(c) for c in obs.columns}
        seen = set()
        for key in priority:
            actual = col_by_lower.get(key.lower())
            if actual and actual not in seen \
                    and _usable_design_col(obs[actual]):
                return actual
            seen.add(actual)
        candidates = []
        for col in obs.columns:
            name = str(col).lower()
            if any(k in name for k in sustring_keys) \
                    and _usable_design_col(obs[col]):
                candidates.append(str(col))
        return candidates[0] if candidates else None

    @staticmethod
    def _pick_h5ad_replicate_col(obs, exclude: list[str] | None = None) -> str | None:
        # orig.ident is Seurat's default; orig_ident (with underscore) is
        # what scanpy emits after round-tripping. Both need to match.
        priority = [
            "sample_id", "orig.ident", "orig_ident", "donor_id", "donor",
            "subject_id", "subject", "individual", "patient_id", "patient",
            "sample", "batch", "library_id", "library",
        ]
        blocked = {c for c in (exclude or []) if c}
        col_by_lower = {str(c).lower(): str(c) for c in obs.columns
                        if str(c) not in blocked}
        for key in priority:
            actual = col_by_lower.get(key.lower())
            if actual and _usable_replicate_col(obs[actual]):
                return actual
        for col in obs.columns:
            if str(col) in blocked:
                continue
            name = str(col).lower()
            if any(k in name for k in ("donor", "subject", "sample",
                                       "individual", "patient", "library")) \
                    and _usable_replicate_col(obs[col]):
                return str(col)
        return None

    @staticmethod
    def _pick_h5ad_groupby_col(obs) -> str | None:
        priority = [
            "cell_type_celltypist", "cell_type", "celltype", "subclass",
            "class", "annotation", "predicted_labels", "cluster",
            "cluster_label", "clusters", "leiden", "louvain",
        ]
        col_by_lower = {str(c).lower(): str(c) for c in obs.columns}
        for key in priority:
            actual = col_by_lower.get(key.lower())
            if actual and _usable_groupby_col(obs[actual]):
                return actual
        return None

    @staticmethod
    def _pick_h5ad_covariates(obs, condition_col: str | None,
                              replicate_col: str | None) -> list[str]:
        covariates = []
        for col in ("Gender", "gender", "sex", "Sex", "batch", "Batch"):
            if col in obs.columns and col not in {condition_col, replicate_col} \
                    and _usable_design_col(obs[col], min_levels=2, max_levels=8):
                covariates.append(col)
        return covariates[:2]

    @staticmethod
    def _default_comparisons_from_groups(groups: dict) -> list[list[str]]:
        names = sorted(groups)
        lower = {g.lower(): g for g in names}
        for ref_key in ("young", "control", "ctrl", "wt", "healthy", "untreated"):
            if ref_key in lower:
                ref = lower[ref_key]
                return [[g, ref] for g in names if g != ref]
        return [[b, a] for i, a in enumerate(names) for b in names[i + 1:]]

    def _build_context(self, experiment_id: str, data_dir: Path,
                       classified: dict, genome: str,
                       organism: str, warnings: list,
                       user_question: str) -> dict:
        # Preserve the bulk_RNA / bulk_RNA_raw distinction here so the
        # NarrativeAgent can report whether the experiment started from raw
        # FASTQs or pre-quantified counts. OrchestratorAgent._dispatch_agents
        # is the right place to merge them before routing to BulkRNAAgent.
        modalities = {k: v for k, v in classified.items() if k != "unknown"}
        input_files = []
        for modality, files in modalities.items():
            for path in files:
                try:
                    p = Path(path)
                    input_files.append({
                        "modality": modality,
                        "path": str(p),
                        "size_bytes": int(p.stat().st_size),
                        "sha256": hash_file(p),
                    })
                except Exception:
                    input_files.append({
                        "modality": modality,
                        "path": str(path),
                        "size_bytes": None,
                        "sha256": "unavailable",
                    })

        exp_context = {
            "experiment_id":  experiment_id,
            "data_dir":       str(data_dir),
            "user_question":  user_question,
            "modalities":     modalities,
            "input_files":    input_files,
            "unknown_files":  classified.get("unknown", []),
            "assay_detections": getattr(self, "_last_assay_detections", []),
            "genome":         genome,
            "organism":       organism,
            "warnings":       warnings,
            "is_multimodal":  len(modalities) > 1,
        }
        exp_context["multiome_contract"] = infer_multiome_contract(exp_context)
        return exp_context

    def _build_checkpoint_summary(self, classified: dict,
                                   genome: str, organism: str,
                                   warnings: list,
                                   inferred_design: dict = None) -> str:
        modalities = {k: v for k, v in classified.items() if k != "unknown"}
        unknown    = classified.get("unknown", [])

        lines = ["📁 ARIA found the following data:\n"]

        ICONS = {
            "scRNA":        "🧬", "bulk_RNA":     "🧬", "bulk_RNA_raw": "🧬",
            "scATAC":       "🔓", "bulk_ATAC":    "🔓",
            "HiC":          "🧵", "ChIP":         "📌",
            "CUT_AND_RUN":  "✂️", "CUT_AND_TAG":  "🏷️",
        }

        for modality, files in modalities.items():
            icon = ICONS.get(modality, "📄")
            lines.append(f"  {icon} {modality}: {len(files)} file(s)")

        if unknown:
            lines.append(f"  ❓ Unrecognized: {len(unknown)} file(s)")

        genome_flag = "⚠️ (inferred)" if genome == "unknown" else "✓"
        lines.append(f"\n  🌍 Organism: {organism}")
        lines.append(f"  🗺️  Genome:   {genome} {genome_flag}")

        if inferred_design and inferred_design.get("source") == "h5ad_obs":
            lines.append("\n  🧾 h5ad obs design hints:")
            lines.append(
                f"     • condition: {inferred_design.get('condition_col', '?')}"
            )
            lines.append(
                f"     • replicate: {inferred_design.get('replicate_col', '?')}"
            )
            if inferred_design.get("groupby_col"):
                lines.append(
                    f"     • cell type/groupby: {inferred_design.get('groupby_col')}"
                )
            covariates = inferred_design.get("covariates") or []
            lines.append(
                "     • covariates: "
                + (", ".join(map(str, covariates)) if covariates else "none")
            )
            groups = inferred_design.get("groups", {})
            if groups:
                compact = ", ".join(f"{g}={len(v)} reps" for g, v in groups.items())
                lines.append(f"     • groups: {compact}")

        if warnings:
            lines.append("\n  ⚠️  Warnings:")
            for w in warnings:
                lines.append(f"     • {w}")

        lines.append("\nIs this correct?")
        return "\n".join(lines)


def _usable_design_col(series, min_levels: int = 2, max_levels: int = 12) -> bool:
    vals = series.dropna().astype(str)
    if vals.empty:
        return False
    levels = [v for v in vals.unique() if v and v.lower() != "nan"]
    if not (min_levels <= len(levels) <= max_levels):
        return False
    counts = vals.value_counts()
    return bool((counts >= 5).sum() >= min_levels)


def _usable_replicate_col(series) -> bool:
    # A replicate column distinguishes biological replicates within
    # conditions. n_levels == 2 is almost always the condition column itself
    # (e.g. Seurat's `orig.ident` mirroring stim/ctrl). Require at least 3
    # distinct levels so the inference does not silently report the
    # condition column as its own replicate. A real n=2 vs n=2 experiment
    # will still have 4 donor IDs in the replicate column.
    vals = series.dropna().astype(str)
    if vals.empty:
        return False
    n_levels = vals.nunique()
    return 3 <= n_levels <= max(200, len(vals) // 20)


def _usable_groupby_col(series) -> bool:
    vals = series.dropna().astype(str)
    if vals.empty:
        return False
    n_levels = vals.nunique()
    return 2 <= n_levels <= 100
