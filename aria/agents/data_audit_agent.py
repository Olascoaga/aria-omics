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

    @staticmethod
    def _is_aria_intermediate_h5ad(path: str) -> bool:
        p = Path(path)
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

    def _filter_aria_intermediate_outputs(
        self, classified: dict[str, list[str]]
    ) -> tuple[dict[str, list[str]], list[str]]:
        """
        Remove ARIA-generated h5ad intermediates from modality inputs when
        source h5ads are present in the same audit. This prevents raw-data
        directories from being misread as multi-sample experiments because a
        previous failed/partial run left `qc_filtered.h5ad` or `clustered.h5ad`
        next to the real input.
        """
        ignored: list[str] = []
        filtered: dict[str, list[str]] = {}
        for modality, files in classified.items():
            if modality != "scRNA":
                filtered[modality] = files
                continue

            intermediates = [
                f for f in files if self._is_aria_intermediate_h5ad(f)
            ]
            sources = [
                f for f in files if not self._is_aria_intermediate_h5ad(f)
            ]
            if sources and intermediates:
                filtered[modality] = sources
                ignored.extend(intermediates)
            else:
                filtered[modality] = files

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
        replicate_col = self._pick_h5ad_replicate_col(obs_all)
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
        priority = [
            "age_group", "age_bin", "condition", "treatment", "diagnosis",
            "disease", "group", "genotype", "phenotype", "status",
        ]
        q = (user_question or "").lower()
        if any(k in q for k in ("age", "aging", "young", "old")):
            priority = ["age_group", "age_bin", "age", *priority]
        for col in priority:
            if col in obs.columns and _usable_design_col(obs[col]):
                return col
        candidates = []
        for col in obs.columns:
            name = str(col).lower()
            if any(k in name for k in ("age", "condition", "treatment", "diagnos",
                                       "disease", "group", "genotype")) \
                    and _usable_design_col(obs[col]):
                candidates.append(str(col))
        return candidates[0] if candidates else None

    @staticmethod
    def _pick_h5ad_replicate_col(obs) -> str | None:
        priority = [
            "sample_id", "orig.ident", "donor_id", "donor", "subject_id",
            "subject", "individual", "sample", "batch",
        ]
        for col in priority:
            if col in obs.columns and _usable_replicate_col(obs[col]):
                return col
        for col in obs.columns:
            name = str(col).lower()
            if any(k in name for k in ("donor", "subject", "sample", "individual")) \
                    and _usable_replicate_col(obs[col]):
                return str(col)
        return None

    @staticmethod
    def _pick_h5ad_groupby_col(obs) -> str | None:
        priority = [
            "cell_type_celltypist", "cell_type", "celltype", "subclass",
            "class", "annotation", "predicted_labels", "leiden",
        ]
        for col in priority:
            if col in obs.columns and _usable_groupby_col(obs[col]):
                return col
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

        return {
            "experiment_id":  experiment_id,
            "data_dir":       str(data_dir),
            "user_question":  user_question,
            "modalities":     modalities,
            "unknown_files":  classified.get("unknown", []),
            "genome":         genome,
            "organism":       organism,
            "warnings":       warnings,
            "is_multimodal":  len(modalities) > 1,
        }

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
    vals = series.dropna().astype(str)
    if vals.empty:
        return False
    n_levels = vals.nunique()
    return 2 <= n_levels <= max(200, len(vals) // 20)


def _usable_groupby_col(series) -> bool:
    vals = series.dropna().astype(str)
    if vals.empty:
        return False
    n_levels = vals.nunique()
    return 2 <= n_levels <= 100
