"""
ARIA DataAuditAgent
-------------------
The gatekeeper. Always runs FIRST, before any analysis.

Responsibilities:
  1. Scan input directory, detect all data types automatically
  2. Infer organism, genome version, experimental design
  3. Validate completeness (pairs, replicates, barcodes)
  4. Trigger CHECKPOINT 1 — "This is what I found, confirm?"
  5. Build ExperimentContext for all downstream agents

Detects automatically:
  scRNA-seq, bulk RNA-seq, scATAC-seq, bulk ATAC-seq,
  HiC / Micro-C, ChIP-seq, CUT&RUN, CUT&TAG,
  Spatial transcriptomics (Visium, Xenium, MERFISH)
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


SIGNATURES = {
    "scRNA": [
        r"barcodes\.tsv(\.gz)?$",
        r"features\.tsv(\.gz)?$",
        r"matrix\.mtx(\.gz)?$",
        r".*\.h5$",
        r".*filtered_feature_bc_matrix.*",
        r".*cellranger.*",
    ],
    "bulk_RNA": [
        r".*counts?\.(txt|tsv|csv)$",
        r".*_R[12]_.*\.fastq(\.gz)?$",
        r".*_[12]\.fastq(\.gz)?$",
    ],
    "scATAC": [
        r"fragments\.tsv(\.gz)?$",
        r".*singlecell\.csv$",
        r".*atac.*barcodes.*",
        r".*scatac.*",
    ],
    "bulk_ATAC": [
        r".*atac.*\.fastq(\.gz)?$",
        r".*atac.*\.bam$",
        r".*atac.*peaks?.*\.(bed|narrowPeak|broadPeak)$",
    ],
    "HiC": [
        r".*\.hic$",
        r".*\.cool$",
        r".*\.mcool$",
        r".*\.pairs(\.gz)?$",
        r".*hic.*\.fastq(\.gz)?$",
        r".*micro.?c.*",
    ],
    "ChIP": [
        r".*chip.*\.fastq(\.gz)?$",
        r".*chip.*\.bam$",
        r".*chip.*peaks?.*\.(bed|narrowPeak|broadPeak)$",
        r".*H[34]K\d+.*\.(fastq|bam)(\.gz)?$",
        r".*input.*\.bam$",
        r".*IgG.*\.bam$",
    ],
    "CUT_AND_RUN": [
        r".*cut.?and.?run.*\.(fastq|bam)(\.gz)?$",
        r".*cutnrun.*\.(fastq|bam)(\.gz)?$",
    ],
    "CUT_AND_TAG": [
        r".*cut.?and.?tag.*\.(fastq|bam)(\.gz)?$",
        r".*cuttag.*\.(fastq|bam)(\.gz)?$",
    ],
    "spatial": [
        r".*tissue_positions.*\.csv$",
        r".*scalefactors_json\.json$",
        r".*tissue_hires_image\.png$",
        r".*visium.*",
        r".*xenium.*",
        r".*merfish.*",
        r".*slide.?seq.*",
        r".*cosmx.*",
        r".*spatial.*\.h5ad$",
    ],
}

GENOME_HINTS = {
    "hg38": ["hg38", "GRCh38", "human"],
    "hg19": ["hg19", "GRCh37"],
    "mm10": ["mm10", "GRCm38", "mouse"],
    "mm39": ["mm39", "GRCm39"],
    "dm6":  ["dm6", "drosophila"],
    "ce11": ["ce11", "worm", "elegans"],
}

ORGANISM_HINTS = {
    "Homo sapiens":            ["hg38", "hg19", "human", "sapiens"],
    "Mus musculus":            ["mm10", "mm39", "mouse", "musculus"],
    "Drosophila melanogaster": ["dm6", "drosophila"],
    "C. elegans":              ["ce11", "elegans", "worm"],
}


class DataAuditAgent(BaseAgent):

    name        = "data_audit_agent"
    description = (
        "Scans input directory, detects all omics data types, "
        "validates experimental design, triggers Checkpoint 1."
    )

    def run(self, experiment_id: str, context: dict) -> dict:
        data_dir      = Path(context["data_dir"])
        user_question = context.get("user_question", "")

        self.publish_status(experiment_id, f"Scanning {data_dir}...", 0.0)

        all_files = self._scan_directory(data_dir)
        if not all_files:
            return {"status": "failed", "error": f"No files found in {data_dir}"}

        classified            = self._classify_files(all_files)
        genome, organism      = self._infer_genome_organism(all_files, data_dir, user_question)
        warnings              = self._validate_design(classified)
        exp_context           = self._build_context(
            experiment_id, data_dir, classified,
            genome, organism, warnings, user_question
        )

        self.memory.create_wing(experiment_id, name=data_dir.name,
                                organism=organism, genome=genome)
        for modality in classified:
            self.memory.create_hall(f"{experiment_id}_{modality}",
                                    experiment_id, modality)

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
        return {"status": "awaiting_checkpoint", "checkpoint": 1,
                "exp_context": exp_context}

    def _scan_directory(self, data_dir: Path) -> list[Path]:
        extensions = {
            ".fastq", ".gz", ".bam", ".bai", ".sam",
            ".bed", ".narrowPeak", ".broadPeak", ".bigWig", ".bw",
            ".hic", ".cool", ".mcool", ".pairs",
            ".h5", ".h5ad", ".loom",
            ".mtx", ".tsv", ".csv", ".txt", ".json", ".png",
        }
        files = []
        try:
            for f in data_dir.rglob("*"):
                if f.is_file():
                    if (f.suffix in extensions or
                            "".join(f.suffixes) in
                            {".fastq.gz", ".pairs.gz", ".tsv.gz", ".bed.gz"}):
                        files.append(f)
        except PermissionError:
            pass
        return files

    def _classify_files(self, files: list[Path]) -> dict[str, list[str]]:
        classified: dict[str, list[str]] = {}
        for f in files:
            matched = False
            for modality, patterns in SIGNATURES.items():
                for pat in patterns:
                    if re.search(pat, f.name, re.IGNORECASE) or \
                       re.search(pat, str(f), re.IGNORECASE):
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
        text = " ".join([
            str(data_dir), user_question,
            " ".join(f.name for f in files[:50])
        ]).lower()

        genome   = "unknown"
        organism = "unknown"

        for g, hints in GENOME_HINTS.items():
            if any(h.lower() in text for h in hints):
                genome = g
                break

        for org, hints in ORGANISM_HINTS.items():
            if any(h.lower() in text for h in hints):
                organism = org
                break

        if genome != "unknown" and organism == "unknown":
            genome_to_org = {
                "hg38": "Homo sapiens", "hg19": "Homo sapiens",
                "mm10": "Mus musculus", "mm39": "Mus musculus",
                "dm6":  "Drosophila melanogaster",
                "ce11": "C. elegans",
            }
            organism = genome_to_org.get(genome, "unknown")

        return genome, organism

    def _validate_design(self, classified: dict) -> list[str]:
        warnings = []
        for modality, files in classified.items():
            if modality == "unknown":
                continue
            r1 = [f for f in files if "_R1_" in f or "_1.fastq" in f]
            r2 = [f for f in files if "_R2_" in f or "_2.fastq" in f]
            if r1 and len(r1) != len(r2):
                warnings.append(
                    f"{modality}: Unequal R1/R2 pairs "
                    f"(R1={len(r1)}, R2={len(r2)})"
                )
            if modality in ("bulk_RNA", "bulk_ATAC", "ChIP") and len(files) < 4:
                warnings.append(
                    f"{modality}: Only {len(files)} files detected. "
                    f"Differential analysis needs >= 2 replicates per condition."
                )
        return warnings

    def _build_context(self, experiment_id, data_dir, classified,
                       genome, organism, warnings, user_question) -> dict:
        return {
            "experiment_id":  experiment_id,
            "data_dir":       str(data_dir),
            "user_question":  user_question,
            "modalities":     {k: v for k, v in classified.items()
                               if k != "unknown"},
            "unknown_files":  classified.get("unknown", []),
            "genome":         genome,
            "organism":       organism,
            "warnings":       warnings,
            "is_multimodal":  len([k for k in classified
                                   if k != "unknown"]) > 1,
        }

    def _build_checkpoint_summary(self, classified, genome,
                                   organism, warnings) -> str:
        modalities = {k: v for k, v in classified.items() if k != "unknown"}
        unknown    = classified.get("unknown", [])
        ICONS = {
            "scRNA": "scRNA-seq", "bulk_RNA": "bulk RNA-seq",
            "scATAC": "scATAC-seq", "bulk_ATAC": "bulk ATAC-seq",
            "HiC": "HiC/Micro-C", "ChIP": "ChIP-seq",
            "CUT_AND_RUN": "CUT&RUN", "CUT_AND_TAG": "CUT&TAG",
            "spatial": "Spatial transcriptomics",
        }
        lines = ["ARIA found the following data:\n"]
        for modality, files in modalities.items():
            label = ICONS.get(modality, modality)
            lines.append(f"  [+] {label}: {len(files)} file(s)")
        if unknown:
            lines.append(f"  [?] Unrecognized: {len(unknown)} file(s)")
        genome_flag = "(inferred)" if genome == "unknown" else ""
        lines.append(f"\n  Organism: {organism}")
        lines.append(f"  Genome:   {genome} {genome_flag}")
        if warnings:
            lines.append("\n  Warnings:")
            for w in warnings:
                lines.append(f"    * {w}")
        lines.append("\nIs this correct?")
        return "\n".join(lines)
