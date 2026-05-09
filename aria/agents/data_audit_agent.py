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
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence, CavemanMode
from aria.memory.memory import ARIAMemory


# ── File signature patterns ──────────────────────────────────────────────────

SIGNATURES = {
    # Single-cell RNA
    "scRNA": [
        r"barcodes\.tsv(\.gz)?$",
        r"features\.tsv(\.gz)?$",
        r"matrix\.mtx(\.gz)?$",
        r".*\.h5$",
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
    # Single-cell ATAC
    "scATAC": [
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
    "hg38": ["hg38", "GRCh38", "human"],
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


class DataAuditAgent(BaseAgent):

    name = "data_audit_agent"
    description = (
        "Scans input directory, detects all omics data types, "
        "validates experimental design, triggers Checkpoint 1."
    )

    def __init__(self, memory: ARIAMemory, api_key: str = None):
        super().__init__(memory, api_key)

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

        self.publish_status(experiment_id,
            f"Scanning {data_dir}...", progress=0.0)

        # 1. Scan files
        all_files = self._scan_directory(data_dir)
        if not all_files and not geo_metadata:
            return {
                "status": "failed",
                "error": f"No files found in {data_dir}"
            }

        # 2. Classify files by modality
        classified = self._classify_files(all_files)

        # 2b. When GEO metadata is present, enrich classification with
        #     already-typed files and remove them from "unknown".
        if geo_metadata:
            geo_files     = geo_metadata.get("files", {})
            geo_data_type = geo_metadata.get("data_type", "bulk_RNA")
            modality      = "scRNA" if geo_data_type == "scRNA" else "bulk_RNA"

            for ftype, bucket in (("counts", modality), ("h5ad", "scRNA"),
                                  ("h5", "scRNA"), ("mtx", "scRNA")):
                for fpath in geo_files.get(ftype, []):
                    if Path(fpath).exists():
                        classified.setdefault(bucket, []).append(fpath)
                        unknown = classified.get("unknown", [])
                        if fpath in unknown:
                            unknown.remove(fpath)

        # 3. Infer organism and genome (GEO metadata takes precedence)
        if geo_metadata:
            inferred = geo_metadata.get("inferred_design", {})
            organism = inferred.get("organism", "") or geo_metadata.get("organism", "unknown")
            genome   = inferred.get("genome",   "") or geo_metadata.get("genome",   "unknown")
        else:
            genome, organism = self._infer_genome_organism(
                all_files, data_dir, user_question
            )

        # 4. Validate design (replicates, pairs, etc.)
        warnings = self._validate_design(classified)

        # 5. Build ExperimentContext
        exp_context = self._build_context(
            experiment_id, data_dir, classified,
            genome, organism, warnings, user_question
        )

        # Propagate GEO metadata into exp_context so downstream agents can use it
        if geo_metadata:
            exp_context["geo_metadata"]    = geo_metadata
            exp_context["inferred_design"] = geo_metadata.get("inferred_design", {})

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
            classified, genome, organism, warnings
        )

        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=1,
            question=checkpoint_msg,
            options=["Confirm and continue", "Correct metadata", "Cancel"],
            context={
                "classified":  classified,
                "genome":      genome,
                "organism":    organism,
                "warnings":    warnings,
                "file_count":  len(all_files),
                "exp_context": exp_context,
            }
        )

        return {
            "status":      "awaiting_checkpoint",
            "checkpoint":  1,
            "exp_context": exp_context
        }

    # ── PRIVATE METHODS ──────────────────────────────────────────────────

    def _scan_directory(self, data_dir: Path) -> list[Path]:
        """Recursively scan directory for all files."""
        extensions = {
            ".fastq", ".gz", ".bam", ".bai", ".sam",
            ".bed", ".narrowPeak", ".broadPeak", ".bigWig", ".bw",
            ".hic", ".cool", ".mcool", ".pairs",
            ".h5", ".h5ad", ".loom",
            ".mtx", ".tsv", ".csv", ".txt",
        }
        files = []
        try:
            for f in data_dir.rglob("*"):
                if f.is_file():
                    # Check by extension or compound extension (.fastq.gz)
                    if (f.suffix in extensions or
                            "".join(f.suffixes) in {".fastq.gz", ".pairs.gz",
                                                     ".tsv.gz", ".bed.gz"}):
                        files.append(f)
        except PermissionError as e:
            self.publish_status("", f"Permission error: {e}")
        return files

    def _classify_files(self, files: list[Path]) -> dict[str, list[str]]:
        """Map files to their omics modality."""
        classified: dict[str, list[str]] = {}
        
        for f in files:
            fname = f.name.lower()
            fpath = str(f).lower()
            
            matched = False
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
        
        return classified

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

    def _build_context(self, experiment_id: str, data_dir: Path,
                       classified: dict, genome: str,
                       organism: str, warnings: list,
                       user_question: str) -> dict:
        # Merge bulk_RNA_raw into bulk_RNA so the Orchestrator
        # routes to BulkRNAAgent which handles raw FASTQs.
        # BulkRNAAgent detects FASTQs internally via _is_fastq().
        merged = dict(classified)
        if "bulk_RNA_raw" in merged:
            raw_files = merged.pop("bulk_RNA_raw")
            merged.setdefault("bulk_RNA", []).extend(raw_files)

        modalities = {k: v for k, v in merged.items() if k != "unknown"}

        return {
            "experiment_id":  experiment_id,
            "data_dir":       str(data_dir),
            "user_question":  user_question,
            "modalities":     modalities,
            "unknown_files":  merged.get("unknown", []),
            "genome":         genome,
            "organism":       organism,
            "warnings":       warnings,
            "is_multimodal":  len(modalities) > 1,
        }

    def _build_checkpoint_summary(self, classified: dict,
                                   genome: str, organism: str,
                                   warnings: list) -> str:
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

        if warnings:
            lines.append("\n  ⚠️  Warnings:")
            for w in warnings:
                lines.append(f"     • {w}")

        lines.append("\nIs this correct?")
        return "\n".join(lines)
