"""
ARIA ChromatinAgent
-------------------
Handles all chromatin accessibility and protein-DNA interaction data:
  - scATAC-seq / bulk ATAC-seq
  - ChIP-seq (histone marks and TF binding)
  - CUT&RUN (high-resolution TF and histone profiling)
  - CUT&TAG (in situ tagmentation version of CUT&RUN)

Key differences from RNAAgent:
  - Data is inherently sparse (scATAC: ~1-5% of peaks accessible per cell)
  - Dimensionality reduction uses LSI (TF-IDF + SVD), NOT PCA
  - First SVD component MUST be discarded (correlates with sequencing depth,
    not biology) — ParameterAdvisor handles this decision explicitly
  - QC metrics are different: TSS enrichment, fragment size distribution,
    FRiP (Fraction of Reads in Peaks)
  - Peak calling requires MACS3 with ATAC-specific parameters
  - Tn5 insertion bias must be accounted for in footprinting analysis

All heavy computation is delegated to EnvironmentManager (aria-chromatin-env).
DebateCouncil is invoked for:
  - TF motif enrichment claims (Tn5 bias check mandatory)
  - Differential accessibility interpretations
  - Peak annotations near TSS (promoter vs enhancer assignment)
"""

from __future__ import annotations

import json
import logging
import uuid
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence, MessageType, CavemanMode
from aria.llm.provider import LLMProvider, TaskTier
from aria.llm.parameter_advisor import ParameterAdvisor
from aria.memory.memory import ARIAMemory
from aria.agents.chromatin_helpers import (
    bulk_da_motif_regions as _bulk_da_motif_regions,
    is_fastq as _is_fastq,
    pick_read as _pick_read,
    positive_int as _positive_int,
    rank_bulk_da_motif_peaks as _rank_bulk_da_motif_peaks,
)

log = logging.getLogger("aria.chromatin")


def _is_missing(value) -> bool:
    """True for None / NaN / pandas NA (obs cells that carry no label)."""
    if value is None:
        return True
    if isinstance(value, float):
        return value != value
    try:
        import pandas as pd

        return bool(pd.api.types.is_scalar(value) and pd.isna(value))
    except Exception:
        return False


def _drop_none(params: dict) -> dict:
    """Omit unset optional params so script input contracts (which reject an
    explicit None for non-empty fields) only see values that actually exist."""
    return {k: v for k, v in params.items() if v is not None}


CHROMATIN_SYSTEM = """
You are ARIA's ChromatinAgent — a specialist in chromatin accessibility
and protein-DNA interaction analysis.

Your expertise:
- ATAC-seq: nucleosome positioning, TF footprinting, peak calling
- ChIP-seq: histone mark profiling, peak calling, IDR analysis
- CUT&RUN / CUT&TAG: high-resolution TF binding, low background noise
- Chromatin remodeling: accessibility changes between conditions

Critical knowledge:
- scATAC data is extremely sparse (~1-5% peaks per cell): standard
  normalization and PCA are inappropriate. Use TF-IDF + LSI (SVD).
- The FIRST SVD component in LSI almost always correlates with
  sequencing depth (a technical artifact), not biology. It MUST be
  removed before clustering.
- Tn5 transposase has strong sequence insertion bias (prefers
  certain GC contexts). TF footprinting results must account for
  this bias before concluding a TF is active.
- FRiP (Fraction of Reads in Peaks) < 0.2 indicates poor ATAC quality.
- TSS enrichment score < 4 indicates failed ATAC library prep.
- For CUT&RUN/CUT&TAG: expect much lower background than ChIP-seq;
  use --nomodel --extsize 200 for MACS3 peak calling.

Always think about the biological question. Chromatin accessibility
is a proxy for regulatory activity, not direct proof of transcription.
""".strip()


class ChromatinAgent(BaseAgent):

    name        = "chromatin_agent"
    description = (
        "Chromatin accessibility and protein-DNA interaction analysis. "
        "Handles scATAC-seq, bulk ATAC-seq, ChIP-seq, CUT&RUN, CUT&TAG."
    )
    validation_level = "beta"  # de-alpha 2026-06-15 (ADR-048); requires_ack retained
    dispatch_enabled = True
    REQUIRED_SCRIPTS = (
        "aria/scripts/chromatin_qc.py",
        "aria/scripts/chromatin_peaks.py",
        "aria/scripts/chromatin_peak_counts.py",
        "aria/scripts/chromatin_bulk_diffacc.py",
        "aria/scripts/chromatin_lsi_clustering.py",
        "aria/scripts/chromatin_diffacc.py",
        "aria/scripts/chromatin_motifs.py",
        "aria/scripts/chromatin_regulatory.py",
    )
    PLANNED_SCRIPTS = (
        "aria/scripts/chromatin_footprinting.py",
        "aria/scripts/chromatin_accessibility_score.py",
    )

    # Assay-specific MACS3 parameters
    MACS3_PARAMS = {
        "scATAC":      {"format": "BAMPE", "nomodel": True,  "extsize": 200,
                        "keep_dup": "all",  "nolambda": False},
        "bulk_ATAC":   {"format": "BAMPE", "nomodel": True,  "extsize": 200,
                        "keep_dup": "all",  "nolambda": False},
        "ChIP":        {"format": "BAM",   "nomodel": False, "extsize": 147,
                        "keep_dup": "1",    "nolambda": False},
        "CUT_AND_RUN": {"format": "BAMPE", "nomodel": True,  "extsize": 200,
                        "keep_dup": "all",  "nolambda": True},
        "CUT_AND_TAG": {"format": "BAMPE", "nomodel": True,  "extsize": 200,
                        "keep_dup": "all",  "nolambda": True},
    }

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm, api_key)
        self.advisor = ParameterAdvisor(memory, llm)

        # Import here to avoid loading EnvironmentManager at module level
        from aria.utils.environment_manager import env_manager
        self.env = env_manager

    # ── Main entry point ──────────────────────────────────────────────────

    def run(self, experiment_id: str, context: dict) -> dict:
        """
        Run chromatin analysis based on available modalities.

        context must contain:
          - exp_context: dict from DataAuditAgent
          - biological_intent: parsed intent from OrchestratorAgent
        """
        exp_ctx  = context.get("exp_context", {})
        intent   = context.get("biological_intent", {})
        mods     = exp_ctx.get("modalities", {})

        self.publish_status(experiment_id,
                            "ChromatinAgent starting...", 0.0)

        results = {}

        # Dispatch to modality-specific pipelines
        if "scATAC" in mods:
            self.publish_status(experiment_id,
                                "Processing scATAC-seq...", 0.1)
            results["scATAC"] = self._run_scatac(
                experiment_id, exp_ctx, intent, mods["scATAC"]
            )

        if "bulk_ATAC" in mods:
            self.publish_status(experiment_id,
                                "Processing bulk ATAC-seq...", 0.3)
            results["bulk_ATAC"] = self._run_bulk_atac(
                experiment_id, exp_ctx, intent, mods["bulk_ATAC"]
            )

        if "ChIP" in mods:
            self.publish_status(experiment_id,
                                "Processing ChIP-seq...", 0.5)
            results["ChIP"] = self._run_chip(
                experiment_id, exp_ctx, intent, mods["ChIP"]
            )

        if "CUT_AND_RUN" in mods:
            self.publish_status(experiment_id,
                                "Processing CUT&RUN...", 0.7)
            results["CUT_AND_RUN"] = self._run_cut_and_run(
                experiment_id, exp_ctx, intent, mods["CUT_AND_RUN"]
            )

        if "CUT_AND_TAG" in mods:
            self.publish_status(experiment_id,
                                "Processing CUT&TAG...", 0.8)
            results["CUT_AND_TAG"] = self._run_cut_and_run(
                experiment_id, exp_ctx, intent, mods["CUT_AND_TAG"],
                assay_type="CUT_AND_TAG"
            )

        if not results:
            self.publish_finding(
                experiment_id,
                {"error": "No chromatin modalities detected"},
                Confidence.INSUFFICIENT,
            )
            return {"status": "failed", "reason": "no_chromatin_data"}

        self.publish_status(experiment_id,
                            "ChromatinAgent complete.", 1.0)
        return {"status": "done", "findings": results}

    # ── scATAC-seq pipeline ───────────────────────────────────────────────

    def _run_scatac(self, experiment_id: str, exp_ctx: dict,
                   intent: dict, files: list) -> dict:
        """
        Single-cell ATAC-seq pipeline.
        Key: LSI dimensionality reduction, discard SVD component 1.
        """
        # v4.6 peak-matrix path: a 10x ARC `.h5mu` carries pre-called peaks, so
        # the pipeline is QC -> LSI/clustering -> differential accessibility ->
        # motif enrichment (no MACS3 peak calling). Fragments/BAM inputs fall
        # through to the legacy peak-calling path below.
        # ADR-076: one pre-built matrix (`.h5mu` or cell x peak `.h5ad`) routes
        # straight to the matrix pipeline; several matrices (a multi-sample
        # cohort) are first unioned onto one consensus peak set.
        if not _is_fastq(files):
            matrix_files = [f for f in files
                            if str(f).lower().endswith((".h5mu", ".h5ad"))]
            # ADR-079: Cell Ranger ARC `.h5` (`10x_multiome_h5`) and ARC dual
            # MEX directories (`10x_multiome_mtx`) carry pre-called paired
            # RNA+ATAC matrices; materialize them into `.h5mu` so the peak-matrix
            # pipeline runs cleanly instead of falling through to MACS3.
            paired_arc_formats = {"10x_multiome_h5", "10x_multiome_mtx"}
            detected_arc_paths = {
                str(Path(str(rec.get("path"))).resolve())
                for rec in (exp_ctx.get("assay_detections") or [])
                if isinstance(rec, dict)
                and rec.get("path")
                and str((rec.get("evidence") or {}).get("format") or "")
                in paired_arc_formats
            }
            for idx, f in enumerate(files):
                f_str = str(f)
                if f_str in matrix_files:
                    continue
                p_obj = Path(f_str)
                is_arc = (
                    (p_obj.exists() and str(p_obj.resolve()) in detected_arc_paths)
                    or (
                        p_obj.exists()
                        and (p_obj.is_dir() or f_str.lower().endswith(".h5"))
                        and not detected_arc_paths
                    )
                )
                if not is_arc or not p_obj.exists():
                    continue
                asm_dir = (
                    Path(str(exp_ctx.get("data_dir") or p_obj.parent))
                    / "aria_chromatin"
                    / f"assembled_arc_{idx}"
                )
                asm_res = self.env.run_in_stack(
                    stack="chromatin",
                    script_path="aria/scripts/integration_assemble_multiome.py",
                    params={
                        "mode": "cellranger_arc",
                        "arc_path": f_str,
                        "output_dir": str(asm_dir),
                        "min_shared_cells": 1,
                    },
                )
                if asm_res.get("status") == "success" and asm_res.get("output_path"):
                    matrix_files.append(str(asm_res["output_path"]))
            if len(matrix_files) == 1:
                return self._run_scatac_matrix(
                    experiment_id, exp_ctx, intent, matrix_files[0])
            if len(matrix_files) >= 2:
                return self._run_scatac_cohort_matrices(
                    experiment_id, exp_ctx, intent, matrix_files)

        findings = {}

        # 0. FASTQ -> fragments (chromap): raw single-cell reads must be aligned
        # to a barcoded fragments file before QC/peaks. Only runs when the
        # barcode read + whitelist are explicitly provided (the barcode read
        # cannot be guessed); otherwise the raw FASTQ falls through honestly.
        if _is_fastq(files):
            frag_attempt = self._run_scatac_fastq_to_fragments(
                experiment_id, exp_ctx, files)
            findings["fastq_to_fragments"] = frag_attempt
            if (isinstance(frag_attempt, dict)
                    and frag_attempt.get("status") == "success"
                    and frag_attempt.get("fragments_file")):
                files = [frag_attempt["fragments_file"]]
                if frag_attempt.get("library_manifest"):
                    exp_ctx = {
                        **exp_ctx,
                        "scatac_fastq_manifest": {
                            "schema_version": "1",
                            "libraries": frag_attempt["library_manifest"],
                        },
                    }
            else:
                return {"status": "done", "findings": findings,
                        "reason": frag_attempt.get("reason")
                        if isinstance(frag_attempt, dict) else "scatac_fastq_unaligned"}

        # 1. QC via EnvironmentManager
        self.publish_status(experiment_id, "scATAC QC...", 0.1)
        qc_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_qc.py",
            params={
                "data_type":    "scATAC",
                "files":        files,
                "genome":       exp_ctx.get("genome", "hg38"),
                "organism":     exp_ctx.get("organism", "Homo sapiens"),
            },
        )

        if qc_result.get("status") == "error":
            log.warning(f"scATAC QC failed: {qc_result.get('details','')}")
            # Non-fatal: continue with available data
            findings["qc"] = qc_result
        else:
            findings["qc"] = qc_result
            self._publish_qc_finding(experiment_id, qc_result, "scATAC")

            # Checkpoint 3 if QC warns
            if qc_result.get("warnings"):
                self.publish_escalation(
                    experiment_id=experiment_id,
                    checkpoint=3,
                    question=self._format_qc_checkpoint(qc_result, "scATAC"),
                    options=[
                        "Continue with warnings noted",
                        "Exclude low-quality samples",
                        "Abort analysis",
                    ],
                    context={"qc": qc_result, "modality": "scATAC"},
                )

        # 2. Peak calling. Raw fragments/BAM inputs do not yet produce the
        # validated cell x peak matrix that the LSI/DA/regulatory pipeline needs.
        self.publish_status(experiment_id, "Calling peaks (MACS3)...", 0.5)
        peaks_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_peaks.py",
            params={
                "data_type":   "scATAC",
                "files":       files,
                "genome":      exp_ctx.get("genome", "hg38"),
                "macs3_params": self.MACS3_PARAMS["scATAC"],
                "output_dir":  str(Path(files[0]).parent / "peaks"),
            },
        )
        findings["peaks"] = peaks_result

        if peaks_result.get("status") == "success":
            self._publish_peaks_finding(
                experiment_id, peaks_result, "scATAC"
            )

        # 2b. Bridge raw fragments -> cell x peak matrix (snapatac2), reusing the
        # MACS3 peaks just called as the peak set. On success this unlocks the
        # full validated matrix pipeline (LSI -> DA -> motif -> regulatory ->
        # footprint) directly from raw fragments; on failure ARIA falls back to
        # the honest QC + peaks (+ optional motif) C4 path.
        matrix_attempt = self._run_fragments_to_matrix(
            experiment_id, exp_ctx, files, peaks_result)
        findings["fragments_to_matrix"] = matrix_attempt
        if (isinstance(matrix_attempt, dict)
                and matrix_attempt.get("status") == "success"
                and matrix_attempt.get("output_path")):
            self.publish_status(
                experiment_id, "scATAC matrix pipeline from fragments...", 0.6)
            # Downstream regulatory / TOBIAS layers need the fragments (every
            # cohort file) and the MACS3 peak universe just called (ADR-078).
            frags = [f for f in files if str(f).lower().endswith(
                (".tsv", ".tsv.gz", ".bed", ".bed.gz"))]
            if not frags and len(files) == 1:
                frags = [files[0]]
            matrix_ctx = dict(exp_ctx)
            if frags:
                matrix_ctx.setdefault("fragments_file", frags[0])
            if len(frags) > 1:
                matrix_ctx.setdefault("fragments_files", frags)
            if (isinstance(peaks_result, dict)
                    and peaks_result.get("status") == "success"):
                called_bed = (peaks_result.get("consensus_peaks_path")
                              or peaks_result.get("peaks_path"))
                if called_bed:
                    matrix_ctx.setdefault("peaks_bed", called_bed)
            matrix_result = self._run_scatac_matrix(
                experiment_id, matrix_ctx, intent,
                matrix_attempt["output_path"], precomputed_qc=qc_result)
            findings.update(matrix_result.get("findings", {}))
            return {"status": "done", "findings": findings}

        # Fallback: matrix not built -> stay on the honest raw path (C4).
        findings["matrix_pipeline"] = {
            "status": "skipped",
            "ran": False,
            "data_type": "scATAC",
            "analysis": "matrix_pipeline",
            "validation_level": "beta",
            "reason": (matrix_attempt.get("reason")
                       if isinstance(matrix_attempt, dict)
                       else "raw_scatac_matrix_bridge_unavailable"),
            "message": (
                "Raw scATAC fragments support measured QC, MACS3 peak calling, "
                "and optional motif enrichment. The fragments->matrix bridge "
                "(snapatac2) could not build a cell x peak matrix, so LSI "
                "clustering, differential accessibility, and regulatory layers "
                "were skipped (they need a peak matrix)."
            ),
            "required_input": "fragments.tsv(.gz) the snapatac2 bridge can read",
            "completed_steps": [
                "qc",
                "peak_calling",
                "motif_enrichment_if_requested",
            ],
            "skipped_steps": [
                "lsi_clustering",
                "differential_accessibility",
                "regulatory_layers",
            ],
        }

        # 3. TF motif enrichment (if requested)
        if self._needs_tf_analysis(intent):
            self.publish_status(experiment_id,
                                "TF motif enrichment...", 0.7)
            motif_result = self._run_motif_enrichment(
                experiment_id, exp_ctx, intent, peaks_result, "scATAC"
            )
            findings["motifs"] = motif_result

        return {"status": "done", "findings": findings}

    def _run_fragments_to_matrix(self, experiment_id: str, exp_ctx: dict,
                                 files: list, peaks_result: dict) -> dict:
        """Dispatch the snapatac2 fragments->cell x peak matrix bridge.

        Detects the fragments file among the inputs and reuses the MACS3 peaks
        already called (peaks_result) as the peak set when available (the
        deterministic 'provided' mode); otherwise the bridge calls peaks de novo.
        Returns the bridge result (honest skip dict when no fragments file is
        present or snapatac2/genome are unavailable).
        """
        frags = [f for f in files
                 if "fragment" in str(f).lower()
                 and str(f).lower().endswith(
                     (".tsv", ".tsv.gz", ".bed", ".bed.gz"))]
        if not frags:
            # A single tsv/bed.gz with no "fragment" token: accept only when it
            # is the lone input (avoid grabbing peak/annotation BEDs).
            tsvs = [f for f in files
                    if str(f).lower().endswith((".tsv", ".tsv.gz"))]
            frags = tsvs if len(tsvs) == 1 else []
        if not frags:
            return {"status": "skipped", "ran": False,
                    "analysis": "fragments_to_peak_matrix",
                    "reason": "no_fragments_file"}
        frag = frags[0]

        peak_file = None
        if isinstance(peaks_result, dict) and peaks_result.get("status") == "success":
            peak_file = (peaks_result.get("consensus_peaks_path")
                         or peaks_result.get("peaks_path"))
        params = {
            "fragments_file": frag,
            "genome": exp_ctx.get("genome", "hg38"),
            "peak_file": peak_file,
            "output_dir": str(Path(frag).parent / "aria_chromatin"),
            "min_fragments": int(exp_ctx.get("min_fragments", 1000)),
            "resolution": exp_ctx.get("leiden_resolution", 1.0),
        }
        if len(frags) > 1:
            # Multi-sample cohort: the bridge unions every fragments file
            # (sample-prefixed barcodes) onto one peak set.
            params["fragments_files"] = frags
            for key in ("sample_ids", "sample_metadata"):
                if exp_ctx.get(key) is not None:
                    params[key] = exp_ctx[key]
        genome_fasta = exp_ctx.get("genome_fasta")
        genome_fai = exp_ctx.get("genome_fai")
        if not genome_fasta and not genome_fai:
            # Air-gapped: chrom sizes cannot be downloaded, so reuse a locally
            # staged FASTA (+ its .fai) when one exists.
            from aria.utils import genomes, privacy
            if not privacy.egress_allowed():
                local, _src = genomes.resolve_local_genome_fasta(
                    exp_ctx.get("genome"))
                genome_fasta = local
        if genome_fasta:
            params["genome_fasta"] = genome_fasta
        if genome_fai:
            params["genome_fai"] = genome_fai
        library_manifest = (
            (exp_ctx.get("scatac_fastq_manifest") or {}).get("libraries")
        )
        if library_manifest:
            params["library_manifest"] = library_manifest
        return self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_fragments_to_matrix.py",
            params=params,
        )

    def _run_scatac_fastq_to_fragments(self, experiment_id: str, exp_ctx: dict,
                                       files: list) -> dict:
        """Dispatch the chromap scATAC FASTQ->fragments front-end (atacseq stack).

        E5: a typed manifest is authoritative and may carry multiple libraries.
        The legacy single-library keys remain accepted for compatibility, but
        filenames never decide library/sample/donor identity.
        """
        fastqs = [f for f in files if _is_fastq([f])]
        from aria.utils.scatac_fastq_manifest import (
            resolve_scatac_fastq_manifest,
        )
        manifest_result = exp_ctx.get("scatac_fastq_manifest_validation")
        if not isinstance(manifest_result, dict) or manifest_result.get("status") == "missing":
            manifest_result = resolve_scatac_fastq_manifest(
                exp_ctx, require_paths=True
            )

        libraries = None
        if manifest_result.get("status") == "invalid":
            return {
                "status": "skipped",
                "ran": False,
                "analysis": "scatac_fastq_to_fragments",
                "reason": "invalid_scatac_fastq_manifest",
                "errors": list(manifest_result.get("errors") or []),
                "message": "The typed scATAC FASTQ manifest is invalid.",
            }
        if manifest_result.get("status") == "valid":
            libraries = manifest_result["manifest"]["libraries"]

        r1 = exp_ctx.get("r1_fastq") or _pick_read(fastqs, ("_r1", "_R1"))
        r3 = (exp_ctx.get("r3_fastq") or exp_ctx.get("r2_fastq")
              or _pick_read(fastqs, ("_r3", "_R3")))
        barcode = (exp_ctx.get("barcode_fastq") or exp_ctx.get("cell_barcode_fastq"))
        whitelist = (exp_ctx.get("barcode_whitelist")
                     or exp_ctx.get("cell_barcode_whitelist"))

        if libraries is None and (not barcode or not whitelist):
            return {"status": "skipped", "ran": False,
                    "analysis": "scatac_fastq_to_fragments",
                    "reason": "missing_barcode_inputs",
                    "message": ("scATAC FASTQ alignment needs an explicit cell-"
                                "barcode read (barcode_fastq) and 10x whitelist "
                                "(barcode_whitelist); the barcode read cannot be "
                                "inferred from the genomic FASTQs.")}

        from aria.utils import genomes
        assembly = exp_ctx.get("genome", "hg38")
        genome_fasta = exp_ctx.get("genome_fasta")
        if not genome_fasta:
            local, _src = genomes.resolve_local_genome_fasta(assembly)
            genome_fasta = local
        if not genome_fasta:
            return {"status": "skipped", "ran": False,
                    "analysis": "scatac_fastq_to_fragments",
                    "reason": "no_reference_genome_fasta",
                    "message": (f"scATAC FASTQ alignment needs a reference FASTA "
                                f"for assembly '{assembly}'.")}

        self.publish_status(
            experiment_id, "scATAC FASTQ alignment (chromap)...", 0.05)
        params = {
            "genome_fasta": genome_fasta,
            "chromap_index": exp_ctx.get("chromap_index"),
            "threads": int(exp_ctx.get("threads", 8)),
        }
        if libraries is not None:
            raw_r1 = libraries[0]["fastqs"]["R1"]
            # Multi-lane libraries carry a list of R1 files per library.
            first_r1 = raw_r1[0] if isinstance(raw_r1, list) else raw_r1
            params.update({
                "libraries": libraries,
                "output_dir": str(Path(first_r1).parent
                                  / "aria_scatac_aligned"),
            })
        else:
            params.update({
                "r1_fastq": r1,
                "r3_fastq": r3,
                "barcode_fastq": barcode,
                "barcode_whitelist": whitelist,
                "output_dir": str(Path((r1 or barcode)).parent
                                  / "aria_scatac_aligned"),
            })
            for key in ("barcode_orientation", "platform", "read_format"):
                if exp_ctx.get(key):
                    params[key] = exp_ctx[key]
        return self.env.run_in_stack(
            stack="atacseq",
            script_path="aria/scripts/chromatin_scatac_align.py",
            params=params,
        )

    def _run_scatac_cohort_matrices(self, experiment_id: str, exp_ctx: dict,
                                    intent: dict, matrix_files: list) -> dict:
        """Union several per-sample scATAC matrices onto one consensus peak set
        (chromatin_merge_matrices.py), then run the matrix pipeline on the
        merged cohort. A failed union is reported honestly; nothing downstream
        runs on a partial cohort."""
        self.publish_status(
            experiment_id, "scATAC cohort matrix union...", 0.05)
        merge_res = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_merge_matrices.py",
            params=_drop_none({
                "files": matrix_files,
                "output_dir": str(Path(matrix_files[0]).parent
                                  / "aria_chromatin_cohort"),
                "sample_ids": exp_ctx.get("sample_ids"),
                "sample_metadata": exp_ctx.get("sample_metadata"),
            }),
        )
        if (not isinstance(merge_res, dict)
                or merge_res.get("status") != "success"
                or not merge_res.get("output_path")):
            merge_res = merge_res if isinstance(merge_res, dict) else {}
            return {"status": "done",
                    "findings": {"cohort_matrix_union": merge_res},
                    "reason": (merge_res.get("error_type")
                               or "cohort_matrix_union_failed")}
        merged_ctx = {
            **exp_ctx,
            "peak_provenance": merge_res.get(
                "peak_provenance", "cohort_interval_union_consensus"),
        }
        matrix_res = self._run_scatac_matrix(
            experiment_id, merged_ctx, intent, merge_res["output_path"])
        return {"status": matrix_res.get("status", "done"),
                "findings": {"cohort_matrix_union": merge_res,
                             **(matrix_res.get("findings") or {})}}

    def _run_scatac_pseudobulk_interpretation(self, experiment_id: str,
                                              exp_ctx: dict, pb: dict,
                                              out_dir: str,
                                              findings: dict) -> None:
        """Genomic annotation + functional ORA of the scATAC pseudobulk DA peaks,
        reusing the bulk ATAC interpretation cores (ADR-076 parity). Runs only
        when the pseudobulk DA actually persisted per-peak result tables."""
        pb_comparisons = [
            {**c, "full_results_csv": (c.get("full_results_csv")
                                       or pb.get("full_results_csv"))}
            for c in pb.get("comparisons") or []
            if isinstance(c, dict) and c.get("status") == "success"
        ]
        if not any(c.get("full_results_csv") for c in pb_comparisons):
            return
        gtf = (exp_ctx.get("gtf") or exp_ctx.get("gtf_path")
               or exp_ctx.get("gene_annotation_gtf"))
        self.publish_status(
            experiment_id, "scATAC DA peak annotation...", 0.65)
        findings["peak_annotation"] = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/chromatin_peak_annotation.py",
            params=_drop_none({
                "data_type": "scATAC",
                "comparisons": pb_comparisons,
                "genome": exp_ctx.get("genome"),
                "gtf": gtf,
                "promoter_upstream": exp_ctx.get("promoter_upstream"),
                "promoter_downstream": exp_ctx.get("promoter_downstream"),
                "output_dir": str(Path(out_dir) / "scatac_annotation"),
            }),
        )
        self.publish_status(
            experiment_id, "scATAC DA peak functional enrichment...", 0.7)
        findings["peak_ora"] = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/chromatin_peak_ora.py",
            params=_drop_none({
                "data_type": "scATAC",
                "comparisons": pb_comparisons,
                "genome": exp_ctx.get("genome"),
                "organism": exp_ctx.get("organism"),
                "gtf": gtf,
                "output_dir": str(Path(out_dir) / "scatac_ora"),
            }),
        )

    def _resolve_da_design(self, exp_ctx: dict) -> dict:
        """Resolve the differential-accessibility design through the single shared
        RNA/ATAC design contract (E7).

        Both ATAC lanes (scATAC and bulk) derive condition/replicate/covariates
        from the confirmed ``exp_ctx["design"]`` here, instead of re-deriving them
        from scattered exp_ctx/design/pseudobulk keys with divergent fallbacks and
        hard-coded ``"condition"``/``"replicate"`` literal defaults. The confirmed
        batch covariate is carried through as a modelled covariate exactly as the
        bulk RNA lane does.
        """
        from aria.utils.design_matrix import resolve_design_contract

        return resolve_design_contract(exp_ctx.get("design"), exp_ctx)

    def _run_scatac_matrix(self, experiment_id: str, exp_ctx: dict,
                           intent: dict, h5mu_path: str,
                           precomputed_qc: Optional[dict] = None) -> dict:
        """v4.6 scATAC peak-matrix pipeline for a peak matrix (`.h5mu` with
        pre-called peaks, or a cell x peak `.h5ad` built from fragments by the
        fragments->matrix bridge): QC -> LSI/clustering -> differential
        accessibility -> motif enrichment.

        Every stage is honest: each result carries its own status, downstream
        stages only run when the upstream produced a usable output, and motif
        enrichment self-skips unless a local genome FASTA + a versioned motif
        collection are present (it never fabricates). Findings are stored under
        the keys the ChromatinNarrator and run-ledger read (`qc`, `lsi`,
        `differential_accessibility`, `motifs`).

        ``precomputed_qc``: when the caller already ran QC (e.g. the raw-fragments
        path computed real TSSe/FRiP on the fragments before building the matrix),
        it is reused instead of re-running QC on the derived matrix.
        """
        findings: dict = {}
        out_dir = str(Path(h5mu_path).parent / "aria_chromatin")

        # 1. QC (chromatin_qc routes a .h5mu through the MuData reader). Reuse the
        # caller's QC when provided (the fragments path already measured it).
        if precomputed_qc is not None:
            findings["qc"] = precomputed_qc
        else:
            self.publish_status(experiment_id, "scATAC QC (.h5mu)...", 0.1)
            qc_result = self.env.run_in_stack(
                stack="chromatin",
                script_path="aria/scripts/chromatin_qc.py",
                params={
                    "data_type": "scATAC",
                    "files": [h5mu_path],
                    "genome": exp_ctx.get("genome", "hg38"),
                    "organism": exp_ctx.get("organism", "Homo sapiens"),
                },
            )
            findings["qc"] = qc_result
            if qc_result.get("status") != "error":
                self._publish_qc_finding(experiment_id, qc_result, "scATAC")

        # 2. LSI + Leiden clustering on the peak matrix
        self.publish_status(experiment_id, "scATAC LSI clustering...", 0.35)
        # E7: one shared RNA/ATAC design contract (condition/replicate/covariates).
        contract = self._resolve_da_design(exp_ctx)
        condition_col = contract["condition_col"]
        replicate_col = contract["replicate_col"]
        batch_col = contract["batch_covariate"]
        design = exp_ctx.get("design") or {}
        peak_provenance = (
            exp_ctx.get("peak_provenance")
            or exp_ctx.get("consensus_peak_provenance")
            or design.get("peak_provenance")
            or design.get("consensus_peak_provenance")
        )
        lsi_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_lsi_clustering.py",
            params={
                "data_path": h5mu_path,
                "resolution": exp_ctx.get("leiden_resolution", 1.0),
                "output_dir": out_dir,
                "condition_col": condition_col,
                "replicate_col": replicate_col,
                "batch_col": batch_col,
                "peak_provenance": peak_provenance,
            },
        )
        findings["lsi"] = lsi_result
        clustered = (lsi_result.get("output_path")
                     if lsi_result.get("status") == "success" else None)

        # 3. Differential accessibility (needs the clustered .h5ad)
        if clustered:
            self.publish_status(
                experiment_id, "Differential accessibility...", 0.6)
            da_params = {
                "data_path": clustered,
                "output_dir": out_dir,
                "exp_context": exp_ctx,  # confirmed thresholds (P0-7)
            }
            # Forward the confirmed cross-condition design via the shared contract
            # (E7): condition/replicate/covariates/comparisons resolved once so the
            # confirmed batch reaches the pseudobulk DESeq2 fit, never silently
            # dropped. Absent metadata makes the pseudobulk lane honestly skip
            # (P0-5), so only forward fields that actually resolved.
            comparisons = contract["comparisons"] or exp_ctx.get("comparisons")
            for key, value in (
                ("condition_col", condition_col),
                ("replicate_col", replicate_col),
                ("comparisons", comparisons),
                ("covariates", contract["covariates"] or None),
            ):
                if value is not None:
                    da_params[key] = value
            # DA runs in the rna stack: pydeseq2 0.5.4 requires numpy>=2 and the
            # chromatin env is pinned numpy<2 (snapatac2/episcanpy/muon). aria-rna-env
            # already ships scanpy + anndata + pydeseq2 on numpy 2, so chromatin_diffacc
            # (per-cluster wilcoxon + pseudobulk DESeq2) runs there with the same
            # pydeseq2 the rest of ARIA's DE/DA core uses. The downstream motif step
            # reads the DA output CSV by path, so it is unaffected by the stack.
            da_result = self.env.run_in_stack(
                stack="rna",
                script_path="aria/scripts/chromatin_diffacc.py",
                params=da_params,
            )
            findings["differential_accessibility"] = da_result

            pb: dict = {}
            if da_result.get("status") == "success":
                pb = da_result.get("pseudobulk") or {}
                if pb.get("ran") and pb.get("comparisons"):
                    self._run_scatac_pseudobulk_interpretation(
                        experiment_id, exp_ctx, pb, out_dir, findings)

            # 4. TF motif enrichment in the DA peak sets (offline, self-gating)
            da_csv = ((da_result.get("per_cluster") or {}).get("output_csv")
                      if da_result.get("status") == "success" else None)
            if da_csv is None and pb.get("output_csv"):
                da_csv = pb.get("output_csv")

            # Resolve the reference genome automatically (no env var required).
            # If nothing is staged locally, guide the user with a checkpoint
            # rather than silently skipping or demanding an env var: ARIA can
            # auto-download the assembly (heavy) only with explicit approval.
            from aria.utils import genomes, privacy
            assembly = exp_ctx.get("genome")
            genome_fasta = exp_ctx.get("genome_fasta")
            local, _src = genomes.resolve_local_genome_fasta(assembly)
            allow_fetch = bool(exp_ctx.get("allow_genome_fetch", False))
            if not genome_fasta and not local and not allow_fetch:
                attr = genomes.snapatac2_attr(assembly)
                opts = []
                if attr and privacy.egress_allowed():
                    opts.append("Download the reference genome now "
                                "(~hundreds of MB)")
                opts += ["Provide a local genome FASTA path",
                         "Skip TF motif enrichment"]
                self.publish_escalation(
                    experiment_id=experiment_id,
                    checkpoint=3,
                    question=(
                        f"TF motif enrichment needs a reference genome for "
                        f"assembly '{assembly or 'unknown'}', which is not staged "
                        f"locally. How should ARIA obtain it?"),
                    options=opts,
                    context={"modality": "scATAC", "assembly": assembly,
                             "need": "genome_fasta"},
                )

            self.publish_status(experiment_id, "TF motif enrichment...", 0.8)
            motif_result = self.env.run_in_stack(
                stack="chromatin",
                script_path="aria/scripts/chromatin_motifs.py",
                params={
                    "data_path": clustered,
                    "da_csv": da_csv,
                    "genome_fasta": genome_fasta,
                    "assembly": assembly,
                    "allow_genome_fetch": allow_fetch,
                    "motif_collection": exp_ctx.get("motif_collection"),
                    "output_dir": out_dir,
                    "exp_context": exp_ctx,
                },
            )
            findings["motifs"] = motif_result

            # 5. Optional P2 regulatory layers. Each sub-layer self-gates on
            # explicit inputs and returns ran:false with a concrete reason when
            # prerequisites are absent; no labels are trusted for inference.
            self.publish_status(experiment_id, "scATAC regulatory layers...", 0.9)
            fragments_file = (
                exp_ctx.get("fragments_file")
                or exp_ctx.get("fragments_path")
                or exp_ctx.get("fragment_file")
            )
            regulatory_result = self.env.run_in_stack(
                stack="chromatin",
                script_path="aria/scripts/chromatin_regulatory.py",
                params={
                    "data_path": clustered,
                    "output_dir": out_dir,
                    "motif_peak_map": exp_ctx.get("motif_peak_map"),
                    "gtf_path": (
                        exp_ctx.get("gtf_path")
                        or exp_ctx.get("gene_annotation_gtf")
                    ),
                    "rna_data_path": (
                        exp_ctx.get("rna_data_path")
                        or exp_ctx.get("scrna_data_path")
                    ),
                    "rna_label_col": (
                        exp_ctx.get("rna_label_col")
                        or exp_ctx.get("scrna_label_col")
                    ),
                    "fragments_file": fragments_file,
                    "motif_sites_bed": exp_ctx.get("motif_sites_bed"),
                    "tn5_bias_model": exp_ctx.get("tn5_bias_model"),
                    "gene_score_window_bp": exp_ctx.get("gene_score_window_bp"),
                    "peak_gene_distance_bp": exp_ctx.get("peak_gene_distance_bp"),
                    "min_peak_gene_corr": exp_ctx.get("min_peak_gene_corr"),
                    "min_shared_cells": exp_ctx.get("min_shared_cells"),
                    **{k: v for k, v in (
                        ("genome_fasta", genome_fasta or local),
                        ("genome_fai", exp_ctx.get("genome_fai")),
                    ) if v is not None},
                },
            )
            findings["regulatory"] = regulatory_result

            # ADR-078: production TOBIAS footprinting runs by default and
            # self-gates on assets; ``enable_footprinting=False`` opts out.
            if exp_ctx.get("enable_footprinting") is not False:
                footprint_result = self._run_scatac_footprinting(
                    experiment_id, exp_ctx, clustered, out_dir)
                if footprint_result is not None:
                    findings["footprinting"] = footprint_result

        return {"status": "done", "findings": findings}

    # ── Bulk ATAC-seq pipeline ────────────────────────────────────────────

    def _run_bulk_atac_preprocessing(self, experiment_id: str,
                                     fastq_files: list,
                                     exp_ctx: dict) -> tuple:
        """Align raw bulk ATAC FASTQs to filtered BAMs (bwa-mem2).

        Two stages, mirroring the bulk RNA prepend: fastp trimming reuses the
        validated ``rna_fastq_qc`` (assay-agnostic) in the rnaseq stack, then
        ``atac_align`` (bwa-mem2 + ATAC filtering) in the atacseq stack. Returns
        ``(bam_paths | None, preprocessing_findings)``; None signals an honest
        failure the caller surfaces as a skip.

        The explicit FASTQ list (plus the audited ``bulk_fastq_manifest`` /
        ``sample_ids`` when present) is forwarded so the trimmer pairs reads from
        the declared files instead of re-globbing the directory.
        """
        fastq_dir = str(Path(fastq_files[0]).parent)
        out_dir = Path(fastq_files[0]).parent.parent / "aria_atac_processing"

        # Resolve the reference FASTA (explicit > auto-resolved local).
        from aria.utils import genomes
        assembly = exp_ctx.get("genome", "hg38")
        genome_fasta = exp_ctx.get("genome_fasta")
        if not genome_fasta:
            local, _src = genomes.resolve_local_genome_fasta(assembly)
            genome_fasta = local
        if not genome_fasta:
            return None, {
                "status": "skipped", "ran": False,
                "analysis": "bulk_atac_alignment", "data_type": "bulk_ATAC",
                "validation_level": "production",
                "reason": "no_reference_genome_fasta",
                "fastq_files": [str(f) for f in fastq_files],
                "message": (
                    f"Bulk ATAC FASTQ alignment needs a reference genome FASTA "
                    f"for assembly '{assembly}', which is not staged locally. "
                    f"Provide exp_context['genome_fasta'] or stage "
                    f"~/.aria/genomes/{assembly}/genome.fa."),
            }

        # 1. fastp trim/QC (reuse the validated bulk RNA trimmer).
        self.publish_status(experiment_id, "bulk ATAC FASTQ QC (fastp)...", 0.10)
        qc_params = {
            "fastq_dir": fastq_dir,
            "fastq_files": [str(f) for f in fastq_files],
            "output_dir": str(out_dir / "qc"),
            "threads": 8,
        }
        for key in ("bulk_fastq_manifest", "sample_ids"):
            if exp_ctx.get(key):
                qc_params[key] = exp_ctx[key]
        qc = self.env.run_in_stack(
            stack="rnaseq",
            script_path="aria/scripts/rna_fastq_qc.py",
            params=qc_params,
        )
        if qc.get("status") != "success" or not qc.get("samples"):
            return None, {"status": "failed", "ran": False,
                          "analysis": "bulk_atac_alignment",
                          "data_type": "bulk_ATAC",
                          "stage": "fastq_qc", "details": qc}

        # 2. bwa-mem2 alignment + ATAC filtering.
        self.publish_status(experiment_id, "bulk ATAC alignment (bwa-mem2)...",
                            0.20)
        align = self.env.run_in_stack(
            stack="atacseq",
            script_path="aria/scripts/atac_align.py",
            params={
                "samples": qc.get("samples", []),
                "genome_fasta": genome_fasta,
                "bwa_index": exp_ctx.get("bwa_index"),
                "output_dir": str(out_dir / "aligned"),
                "threads": int(exp_ctx.get("threads", 8)),
                "mapq": int(exp_ctx.get("atac_mapq", 30)),
                "remove_mito": bool(exp_ctx.get("atac_remove_mito", True)),
            },
        )
        if align.get("status") != "success":
            return None, {"status": "failed", "ran": False,
                          "analysis": "bulk_atac_alignment",
                          "data_type": "bulk_ATAC", "stage": "align",
                          "details": align}

        bams = [b.get("bam") for b in align.get("bam_files", [])
                if b.get("status") == "success" and b.get("bam")]
        if not bams:
            return None, {"status": "failed", "ran": False,
                          "analysis": "bulk_atac_alignment",
                          "data_type": "bulk_ATAC", "stage": "align",
                          "reason": "no_bam_produced", "details": align}

        preprocessing = {
            "status": "success", "ran": True,
            "analysis": "bulk_atac_alignment",
            "aligner": align.get("aligner", "bwa-mem2"),
            "genome_fasta": genome_fasta, "assembly": assembly,
            "filter": align.get("filter"),
            "filters": align.get("filters"),
            "read_layout": align.get("read_layout"),
            "n_samples_aligned": len(bams),
            "fastq_files": [str(f) for f in fastq_files],
            "fastq_qc": {"n_samples": qc.get("n_samples"),
                         "multiqc_report": qc.get("multiqc_report")},
            "bam_files": align.get("bam_files"),
        }
        preprocessing = {k: v for k, v in preprocessing.items()
                         if v is not None}
        preprocessing["data_type"] = "bulk_ATAC"
        preprocessing["validation_level"] = "production"
        return bams, preprocessing

    # Direct peak-count matrix intake (ADR-077). Tabular files that are NOT
    # fragments / BED / narrowPeak / broadPeak are treated as a peak x sample
    # count matrix; the BAM-level QC / MACS3 / counting stages are skipped.
    _MATRIX_SUFFIXES = (".tsv", ".csv", ".txt", ".tsv.gz", ".csv.gz",
                        ".txt.gz")
    _NON_MATRIX_MARKERS = ("fragments.tsv", ".bed", ".narrowpeak",
                           ".broadpeak")
    _ALIGNMENT_SUFFIXES = (".bam", ".cram", ".sam")
    _PEAK_ID_COLUMNS = ("peak_id", "peak", "feature_id")

    @classmethod
    def _is_peak_count_matrix_path(cls, path) -> bool:
        name = Path(str(path)).name.lower()
        if not name.endswith(cls._MATRIX_SUFFIXES):
            return False
        return not any(marker in name for marker in cls._NON_MATRIX_MARKERS)

    @classmethod
    def _detect_precomputed_counts_matrix(cls, exp_ctx: dict,
                                          files: list) -> Optional[str]:
        """Return the peak-count matrix path when this bulk ATAC run starts from
        a pre-computed matrix, else ``None`` (alignment-level lane)."""
        names = [Path(str(f)).name.lower() for f in (files or [])]
        has_alignment_or_interval = any(
            n.endswith(cls._ALIGNMENT_SUFFIXES)
            or any(m in n for m in cls._NON_MATRIX_MARKERS)
            for n in names
        )
        if has_alignment_or_interval:
            return None
        matrices = [str(f) for f in (files or [])
                    if cls._is_peak_count_matrix_path(f)]
        if matrices:
            return matrices[0]
        explicit = exp_ctx.get("counts_matrix_path")
        if explicit:
            return str(explicit)
        return None

    @staticmethod
    def _write_sample_metadata_tsv(sample_metadata, path: Path) -> Optional[str]:
        """Materialise ``exp_ctx['sample_metadata']`` (``{sample: {col: v}}`` or a
        list of records carrying ``sample_id``) as the TSV the DA script reads.
        Returns ``None`` when no usable records exist (never fabricates rows)."""
        records: list[dict] = []
        if isinstance(sample_metadata, dict):
            for sid, meta in sample_metadata.items():
                if isinstance(meta, dict):
                    records.append({"sample_id": str(sid), **meta})
        elif isinstance(sample_metadata, list):
            for meta in sample_metadata:
                if isinstance(meta, dict) and meta.get("sample_id") is not None:
                    records.append({**meta, "sample_id": str(meta["sample_id"])})
        if not records:
            return None
        cols = ["sample_id"]
        for rec in records:
            for key in rec:
                if key not in cols:
                    cols.append(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        import csv
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=cols, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            for rec in records:
                writer.writerow({k: ("" if rec.get(k) is None else rec.get(k))
                                 for k in cols})
        return str(path)

    def _ingest_precomputed_counts_matrix(self, exp_ctx: dict,
                                          matrix_path: str) -> dict:
        """Validate a pre-computed peak-count matrix and record honest direct-
        matrix provenance. Comma-separated input is re-written as the TSV the DA
        script consumes; the numeric content is never altered."""
        import pandas as pd

        src = Path(str(matrix_path))
        base = {
            "analysis": "peak_count_matrix",
            "data_type": "bulk_ATAC",
            "source": "precomputed_matrix",
            "input_counts_matrix_path": str(src),
        }
        if not src.is_file():
            return {**base, "status": "skipped", "ran": False,
                    "validation_level": "production",
                    "reason": f"counts matrix not found: {src}"}
        name = src.name.lower()
        sep = "," if name.endswith((".csv", ".csv.gz")) else "\t"
        try:
            df = pd.read_csv(src, sep=sep)
        except Exception as exc:  # malformed table -> honest skip
            return {**base, "status": "skipped", "ran": False,
                    "validation_level": "production",
                    "reason": f"counts matrix could not be parsed: {exc}"}
        if df.shape[1] < 2 or df.empty:
            return {**base, "status": "skipped", "ran": False,
                    "validation_level": "production",
                    "reason": "counts matrix must contain peak_id plus sample "
                              "columns"}
        if str(df.columns[0]) not in self._PEAK_ID_COLUMNS:
            return {**base, "status": "skipped", "ran": False,
                    "validation_level": "production",
                    "reason": ("first counts matrix column must be one of "
                               f"{list(self._PEAK_ID_COLUMNS)}; got "
                               f"{df.columns[0]!r}")}

        out_dir = src.parent / "bulk_atac_counts"
        counts_path = str(src)
        if sep != "\t" or name.endswith(".gz"):
            out_dir.mkdir(parents=True, exist_ok=True)
            counts_path = str(out_dir / "bulk_atac_peak_counts.tsv")
            df.to_csv(counts_path, sep="\t", index=False)

        meta_path = exp_ctx.get("sample_metadata_path")
        meta_source = "exp_context.sample_metadata_path" if meta_path else None
        if not meta_path:
            meta_path = self._write_sample_metadata_tsv(
                exp_ctx.get("sample_metadata"),
                out_dir / "bulk_atac_sample_metadata.tsv")
            meta_source = ("synthesized_from_exp_context.sample_metadata"
                           if meta_path else None)
        return {
            **base,
            "status": "success",
            "ran": True,
            "validation_level": "production",
            "counting_method": "pre-computed peak-count matrix (supplied)",
            "counts_matrix_path": counts_path,
            "sample_metadata_path": str(meta_path) if meta_path else None,
            "sample_metadata_source": meta_source,
            "n_peaks": int(df.shape[0]),
            "n_samples": int(df.shape[1] - 1),
            "sample_ids": [str(c) for c in df.columns[1:]],
            "warnings": [],
        }

    @staticmethod
    def _merge_peak_frip_into_qc(findings: dict, peaks_result: dict) -> None:
        """Close the bulk ATAC QC gate with the post-peak-calling FRiP.

        BAM-level QC cannot measure FRiP before peaks exist; once MACS3 has run
        and ``chromatin_peaks`` measured FRiP, merge it into ``findings['qc']``.
        ``pass_qc`` is decided only when FRiP and TSS enrichment are both
        measured, against the thresholds the QC step recorded (ENCODE defaults
        FRiP >= 0.20 / TSSe >= 5 when none were recorded)."""
        qc = findings.get("qc")
        if not isinstance(qc, dict) or not isinstance(peaks_result, dict):
            return
        if peaks_result.get("frip") is None:
            return
        qc["frip"] = peaks_result["frip"]
        if peaks_result.get("frip_per_sample"):
            qc["frip_per_sample"] = peaks_result["frip_per_sample"]
        qc["frip_source"] = "chromatin_peaks"
        not_computed = [
            m for m in (qc.get("metrics_not_computed") or [])
            if not str(m).strip().lower().startswith("frip")
        ]
        qc["metrics_not_computed"] = not_computed
        if qc.get("frip") is not None and qc.get("tss_enrichment") is not None:
            th = qc.get("thresholds") if isinstance(qc.get("thresholds"),
                                                    dict) else {}
            min_frip = float(th.get("min_frip", 0.2))
            min_tss = float(th.get("min_tss", 5.0))
            qc["qc_complete"] = True
            qc["pass_qc"] = bool(float(qc["frip"]) >= min_frip
                                 and float(qc["tss_enrichment"]) >= min_tss)

    def _run_bulk_atac(self, experiment_id: str, exp_ctx: dict,
                       intent: dict, files: list) -> dict:
        """Bulk ATAC-seq lane (production spine, ADR-077).

        FASTQ -> BAM (bwa-mem2) -> QC -> MACS3 -> peak-count matrix ->
        replicate-gated DESeq2 DA -> motifs / annotation / ORA (+ beta TOBIAS).
        Raw FASTQ inputs are aligned first; a pre-computed peak x sample count
        matrix skips straight to DA. The DA branch fires whenever the confirmed
        design contract (or the parsed intent) carries a comparison, so a
        ``design.contrasts`` / ``exp_ctx['comparisons']`` run is no longer lost
        when the intent lacks a free-text comparison. Under-specified designs
        skip honestly inside the DA script. scATAC FASTQ is NOT handled here.
        """
        findings = {}

        # FASTQ -> BAM: raw bulk ATAC reads must be aligned before QC/peaks.
        if _is_fastq(files):
            bam_files, preprocessing = self._run_bulk_atac_preprocessing(
                experiment_id, files, exp_ctx)
            findings["preprocessing"] = preprocessing
            if not bam_files:
                return {"status": "failed", "findings": findings,
                        "reason": "bulk_atac_preprocessing_failed"}
            # ADR-078: the aligner names each BAM by its declared sample, so
            # carry that identity forward (counts / replicate design / footprint
            # batches) instead of re-deriving it from BAM filenames.
            if (isinstance(preprocessing, dict)
                    and preprocessing.get("status") == "success"
                    and not exp_ctx.get("sample_ids")):
                aligned_sids = [
                    str(rec.get("sample_id") or rec.get("name"))
                    for rec in (preprocessing.get("bam_files") or [])
                    if isinstance(rec, dict) and rec.get("status") == "success"
                    and rec.get("bam")
                    and (rec.get("sample_id") or rec.get("name"))
                ]
                if len(aligned_sids) == len(bam_files):
                    exp_ctx = {**exp_ctx, "sample_ids": aligned_sids}
            files = bam_files

        # E7: resolve the shared design contract BEFORE deciding whether the
        # count-matrix + DA branch runs.
        contract = self._resolve_da_design(exp_ctx)
        precomputed_matrix = self._detect_precomputed_counts_matrix(
            exp_ctx, files)

        if precomputed_matrix:
            self.publish_status(
                experiment_id, "bulk ATAC pre-computed peak-count matrix...",
                0.60)
            count_result = self._ingest_precomputed_counts_matrix(
                exp_ctx, precomputed_matrix)
            findings["peak_counts"] = count_result
            self._run_bulk_atac_da_branch(
                experiment_id, exp_ctx, intent, contract, findings,
                count_result, peaks_result=None,
                files=[count_result.get("counts_matrix_path")
                       or precomputed_matrix],
                run_footprinting=False)
            return {"status": "done", "findings": findings}

        reference_fasta = exp_ctx.get("reference_fasta")
        if not reference_fasta and any(
                str(f).lower().endswith(".cram") for f in files):
            reference_fasta = exp_ctx.get("genome_fasta")
        gtf_path = (exp_ctx.get("gtf") or exp_ctx.get("gtf_path")
                    or exp_ctx.get("gtf_file"))

        # QC
        self.publish_status(experiment_id, "bulk ATAC QC...", 0.35)
        qc_params = {
            "data_type": "bulk_ATAC",
            "files":     files,
            "genome":    exp_ctx.get("genome", "hg38"),
        }
        if gtf_path:
            qc_params["gtf"] = str(gtf_path)
        if reference_fasta:
            qc_params["reference_fasta"] = str(reference_fasta)
        qc_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_qc.py",
            params=qc_params,
        )
        if isinstance(qc_result, dict) and qc_result.get("status") == "success":
            qc_result["validation_level"] = "production"
            qc_result.setdefault("data_type", "bulk_ATAC")
        findings["qc"] = qc_result
        if qc_result.get("status") == "success":
            self._publish_qc_finding(experiment_id, qc_result, "bulk_ATAC")

        # Peak calling
        self.publish_status(experiment_id, "bulk ATAC peak calling...", 0.55)
        peaks_params = {
            "data_type":    "bulk_ATAC",
            "files":        files,
            "genome":       exp_ctx.get("genome", "hg38"),
            "macs3_params": self.MACS3_PARAMS["bulk_ATAC"],
            "output_dir":   str(Path(files[0]).parent / "peaks"),
        }
        if exp_ctx.get("blacklist_bed"):
            peaks_params["blacklist_bed"] = str(exp_ctx["blacklist_bed"])
        if reference_fasta:
            peaks_params["reference_fasta"] = str(reference_fasta)
        peaks_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_peaks.py",
            params=peaks_params,
        )
        if isinstance(peaks_result, dict) and peaks_result.get("status") == "success":
            peaks_result["validation_level"] = "production"
            peaks_result.setdefault("data_type", "bulk_ATAC")
            # FRiP is measured only after peaks exist: close the QC gate here.
            self._merge_peak_frip_into_qc(findings, peaks_result)
        findings["peaks"] = peaks_result
        if peaks_result.get("status") == "success":
            self._publish_peaks_finding(experiment_id, peaks_result, "bulk_ATAC")

        # Peak-count matrix + replicate-gated DA whenever a comparison exists in
        # the confirmed design contract or the parsed intent.
        run_da = bool(contract.get("comparisons")
                      or intent.get("comparisons")
                      or intent.get("comparison"))
        if run_da:
            peak_universe = (
                peaks_result.get("consensus_peaks_path")
                or peaks_result.get("peaks_path")
            ) if isinstance(peaks_result, dict) else None
            count_result = None
            if peak_universe:
                self.publish_status(
                    experiment_id, "bulk ATAC peak-count matrix...", 0.70)
                count_params = {
                    "data_type": "bulk_ATAC",
                    "files": files,
                    "peaks_path": peak_universe,
                    "sample_ids": exp_ctx.get("sample_ids"),
                    "sample_metadata": exp_ctx.get("sample_metadata"),
                    "output_dir": str(Path(files[0]).parent /
                                      "bulk_atac_counts"),
                }
                if reference_fasta:
                    count_params["reference_fasta"] = str(reference_fasta)
                count_result = self.env.run_in_stack(
                    stack="chromatin",
                    script_path="aria/scripts/chromatin_peak_counts.py",
                    params=count_params,
                )
                if (isinstance(count_result, dict)
                        and count_result.get("status") == "success"):
                    count_result["validation_level"] = "production"
                    count_result.setdefault("data_type", "bulk_ATAC")
                findings["peak_counts"] = count_result
            else:
                findings["peak_counts"] = {
                    "status": "skipped",
                    "reason": "no_peak_universe_for_bulk_atac_counts",
                    "analysis": "peak_count_matrix",
                    "data_type": "bulk_ATAC",
                    "validation_level": "production",
                }
            self._run_bulk_atac_da_branch(
                experiment_id, exp_ctx, intent, contract, findings,
                count_result, peaks_result=peaks_result, files=files,
                run_footprinting=True)

        return {"status": "done", "findings": findings}

    def _run_bulk_atac_da_branch(self, experiment_id: str, exp_ctx: dict,
                                 intent: dict, contract: dict, findings: dict,
                                 count_result: Optional[dict],
                                 peaks_result: Optional[dict], files: list,
                                 run_footprinting: bool = True) -> None:
        """Replicate-gated DESeq2 DA + downstream interpretation over a bulk
        ATAC peak-count matrix (BAM-derived or pre-computed). Mutates
        ``findings`` in place; every non-run stage records an honest reason."""
        if not (isinstance(count_result, dict)
                and count_result.get("status") == "success"):
            findings["differential_accessibility"] = {
                "status": "skipped",
                "reason": "bulk_atac_peak_count_matrix_not_available",
                "analysis": "differential_accessibility",
                "data_type": "bulk_ATAC",
                "validation_level": "production",
                "details": findings.get("peak_counts"),
            }
            return
        if not count_result.get("sample_metadata_path"):
            findings["differential_accessibility"] = {
                "status": "skipped",
                "reason": ("bulk_atac_sample_metadata_missing: provide "
                           "exp_context['sample_metadata'] or "
                           "'sample_metadata_path' with condition/replicate "
                           "columns"),
                "analysis": "differential_accessibility",
                "data_type": "bulk_ATAC",
                "validation_level": "production",
            }
            return

        self.publish_status(
            experiment_id, "bulk ATAC differential accessibility...", 0.82)
        # E7: condition/replicate/covariates come from the shared design
        # contract; the historical literals apply only when no design confirmed
        # a column at all.
        out_base = Path(str(files[0])).parent
        da_params = {
            "data_type": "bulk_ATAC",
            "counts_matrix_path": count_result.get("counts_matrix_path"),
            "sample_metadata_path": count_result.get("sample_metadata_path"),
            "condition_col": contract.get("condition_col") or "condition",
            "replicate_col": contract.get("replicate_col") or "replicate",
            "comparisons": (
                contract.get("comparisons")
                or intent.get("comparisons")
                or intent.get("comparison")
            ),
            "covariates": contract.get("covariates") or [],
            "exp_context": exp_ctx,
            "output_dir": str(out_base / "bulk_atac_da"),
        }
        # Honor an explicit replicate-floor override (e.g. n=2 ENCODE
        # isogenic-replicate designs run with a low_power_warning). Absent an
        # override, the script keeps its production floor.
        min_reps = (exp_ctx.get("min_replicates_per_condition")
                    or intent.get("min_replicates_per_condition"))
        if min_reps is not None:
            da_params["min_replicates_per_condition"] = int(min_reps)
        # DA runs in the rna stack: pydeseq2 lives in aria-rna-env, not the
        # chromatin env. The DA script operates on TSVs (no chromatin deps).
        da_result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/chromatin_bulk_diffacc.py",
            params=da_params,
        )
        if not (isinstance(da_result, dict)
                and da_result.get("status") == "success"
                and da_result.get("ran")):
            findings["differential_accessibility"] = {
                "status": "skipped",
                "reason": (
                    da_result.get("reason")
                    if isinstance(da_result, dict)
                    else "bulk_atac_da_not_run"
                ),
                "analysis": "differential_accessibility",
                "data_type": "bulk_ATAC",
                "validation_level": "production",
                "details": da_result,
            }
            return

        da_result["validation_level"] = "production"
        findings["differential_accessibility"] = da_result
        # TF motif enrichment over the bulk ATAC DA peak sets (offline,
        # self-gating on genome FASTA + motif collection).
        motif_result = self._run_bulk_atac_motifs(
            experiment_id, exp_ctx, da_result, files)
        if motif_result is not None:
            findings["motifs"] = motif_result
        # Genomic annotation of the DA peaks (B2). Honest-skip without a GTF.
        ann_result = self._run_bulk_atac_annotation(
            experiment_id, exp_ctx, da_result, files)
        if ann_result is not None:
            findings["peak_annotation"] = ann_result
        # Functional ORA of genes near the DA peaks (B3).
        ora_result = self._run_bulk_atac_ora(
            experiment_id, exp_ctx, da_result, files)
        if ora_result is not None:
            findings["peak_ora"] = ora_result
        # C2 / ADR-078: condition-level TOBIAS footprinting (production) needs
        # BAMs, so it only runs on the alignment-level lane; it self-gates on
        # assets/env, and ``enable_footprinting=False`` opts out.
        if (run_footprinting and isinstance(peaks_result, dict)
                and exp_ctx.get("enable_footprinting") is not False):
            footprint_result = self._run_bulk_atac_footprinting(
                experiment_id, exp_ctx, da_result, peaks_result, files,
                contract=contract)
            if footprint_result is not None:
                findings["footprinting"] = footprint_result

    def _footprint_common_params(self, exp_ctx: dict, output_dir: str) -> dict:
        """Shared TOBIAS params with auto-resolved reference assets (ADR-078).

        The genome FASTA (+ sibling ``.fai``) and the versioned MEME motif
        collection are resolved from the local ARIA stores when not supplied, so
        the user never has to pass paths. The resolved motif provenance travels as
        the private ``_motif_source`` key, which the dispatchers pop before IPC and
        stamp onto the result as ``motif_source``.
        """
        from aria.utils import genomes, motifs

        genome_fasta = exp_ctx.get("genome_fasta")
        if not genome_fasta:
            try:
                local, _src = genomes.resolve_local_genome_fasta(
                    exp_ctx.get("genome") or "hg38")
            except Exception as exc:  # resolution must never fail the run
                log.warning("footprint genome resolution failed: %s", exc)
                local = None
            if local:
                genome_fasta = local
        genome_fai = exp_ctx.get("genome_fai")
        if (not genome_fai and genome_fasta
                and Path(f"{genome_fasta}.fai").is_file()):
            genome_fai = f"{genome_fasta}.fai"
        motif_meme = (
            exp_ctx.get("motif_meme")
            or exp_ctx.get("motif_meme_path")
            or exp_ctx.get("motif_collection_meme")
        )
        motif_source: Optional[dict] = None
        if not motif_meme:
            try:
                loaded = motifs.load_local_motif_collection(
                    exp_ctx.get("motif_collection") or motifs.DEFAULT_COLLECTION)
            except Exception as exc:
                log.warning("footprint motif resolution failed: %s", exc)
                loaded = None
            if loaded is not None:
                motif_meme, motif_source = loaded
        return {
            "genome_fasta": genome_fasta,
            "genome_fai": genome_fai,
            "peaks_bed": (
                exp_ctx.get("peaks_bed")
                or exp_ctx.get("peaks_path")
                or exp_ctx.get("peak_universe_bed")
            ),
            "motif_meme": motif_meme,
            "_motif_source": motif_source,
            "group_a": (
                exp_ctx.get("footprint_group_a")
                or exp_ctx.get("tobias_group_a")
                or exp_ctx.get("group_a")
            ),
            "group_b": (
                exp_ctx.get("footprint_group_b")
                or exp_ctx.get("tobias_group_b")
                or exp_ctx.get("group_b")
            ),
            "output_dir": output_dir,
            "cores": _positive_int(
                exp_ctx.get("tobias_cores")
                or exp_ctx.get("footprint_cores")
                or 8,
                8,
            ),
            "min_replicates_per_condition": _positive_int(
                exp_ctx.get("footprint_min_replicates_per_condition")
                or exp_ctx.get("min_replicates_per_condition")
                or 3,
                3,
            ),
            "footprint_fdr": float(exp_ctx.get("footprint_fdr", 0.05)),
            "max_label_permutations": _positive_int(
                exp_ctx.get("footprint_max_label_permutations") or 100, 100),
            "permutation_seed": int(exp_ctx.get("footprint_permutation_seed", 0)),
            "make_aggregate_plots": bool(
                exp_ctx.get("footprint_make_aggregate_plots", True)),
            "make_volcano_plot": bool(
                exp_ctx.get("footprint_make_volcano_plot", True)),
            "resume": bool(exp_ctx.get("footprint_resume", True)),
        }

    def _run_scatac_footprinting(self, experiment_id: str, exp_ctx: dict,
                                 clustered: str, out_dir: str) -> Optional[dict]:
        """Inline scATAC TOBIAS dispatch (C3; production per ADR-078).

        The backend refuses uncorrected footprints and returns ran:false for
        missing TOBIAS/assets/groups, so the live flow can attempt footprinting
        without fabricating a result. When the peak BED or the barcode design is
        not supplied, both are derived from the clustered ``.h5ad`` itself (peak
        coordinates from ``var_names``; condition/replicate/batch from ``obs``
        using only the confirmed design columns — never transferred labels).
        """
        fragments_files = exp_ctx.get("fragments_files")
        if not (isinstance(fragments_files, list) and fragments_files):
            fragments_files = None
        fragments_file = (
            exp_ctx.get("fragments_file")
            or exp_ctx.get("fragments_path")
            or exp_ctx.get("fragment_file")
            or (fragments_files[0] if fragments_files else None)
        )
        barcode_groups = (
            exp_ctx.get("barcode_groups")
            or exp_ctx.get("barcode_groups_tsv")
            or exp_ctx.get("barcode_group_tsv")
        )
        footprint_dir = Path(out_dir) / "footprinting"
        params = {
            **self._footprint_common_params(exp_ctx, str(footprint_dir)),
            "mode": "scatac",
            "data_type": "scATAC",
            "data_path": clustered,
            "fragments_file": str(fragments_file) if fragments_file else None,
            "barcode_groups": str(barcode_groups) if barcode_groups else None,
        }
        motif_source = params.pop("_motif_source", None)
        if fragments_files:
            params["fragments_files"] = [str(f) for f in fragments_files]
        sample_ids = exp_ctx.get("sample_ids")
        n_frag_files = len(fragments_files) if fragments_files else 1
        if (isinstance(sample_ids, list) and sample_ids
                and len(sample_ids) == n_frag_files):
            params["sample_ids"] = [str(x) for x in sample_ids]

        clustered_ok = bool(clustered) and Path(str(clustered)).is_file()
        if params.get("peaks_bed") is None and clustered_ok:
            derived = self._derive_peaks_bed_from_h5ad(
                clustered, footprint_dir / "derived_peaks.bed")
            if derived:
                params["peaks_bed"] = derived
        if barcode_groups is None and clustered_ok:
            design = self._derive_scatac_barcode_design(
                {**exp_ctx, "group_a": params.get("group_a"),
                 "group_b": params.get("group_b")},
                clustered, footprint_dir / "derived_barcode_design.tsv")
            if design is not None:
                tsv_path, group_a, group_b = design
                params["barcode_groups"] = tsv_path
                if group_a and group_b:
                    params["group_a"], params["group_b"] = group_a, group_b

        self.publish_status(
            experiment_id, "scATAC TOBIAS footprinting...", 0.95)
        result = self.env.run_in_stack(
            stack="tobias",
            script_path="aria/scripts/chromatin_footprint_tobias.py",
            params=params,
        )
        if isinstance(result, dict):
            result.setdefault("data_type", "scATAC")
            result.setdefault("analysis", "differential_tf_footprinting")
            result["validation_level"] = "production"
            if motif_source and result.get("ran"):
                result["motif_source"] = motif_source
            self._attach_footprint_rna_cross_evidence(result, exp_ctx)
        return result

    # ── ADR-078 helpers: design / peak derivation + RNA cross-evidence ──────

    # Labels that must never become an inferential footprint grouping unless the
    # user names them explicitly via ``footprint_groupby``.
    _NON_INFERENTIAL_GROUP_COLS = frozenset({
        "predicted_celltype", "transferred_label", "predicted_label",
        "predicted_cell_type", "leiden", "louvain", "cluster", "clusters",
    })

    @staticmethod
    def _h5_strings(values) -> list:
        out = []
        for value in values:
            if isinstance(value, bytes):
                value = value.decode("utf-8", errors="replace")
            out.append(value)
        return out

    @classmethod
    def _h5ad_index(cls, group) -> list:
        """Read the index (obs_names / var_names) of an AnnData h5 group."""
        key = group.attrs.get("_index", "_index")
        if isinstance(key, bytes):
            key = key.decode()
        if key not in group:
            return []
        return [str(v) for v in cls._h5_strings(group[key][()])]

    @classmethod
    def _h5ad_obs_column(cls, obs, column: str) -> Optional[list]:
        """Read one ``obs`` column via h5py (modern + legacy categorical)."""
        if column not in obs:
            return None
        node = obs[column]
        import h5py

        if isinstance(node, h5py.Group):
            if "categories" in node and "codes" in node:
                cats = cls._h5_strings(node["categories"][()])
                return [None if int(c) < 0 else cats[int(c)]
                        for c in node["codes"][()]]
            if "values" in node:  # nullable string/int arrays
                values = cls._h5_strings(node["values"][()])
                mask = node["mask"][()] if "mask" in node else None
                return [None if (mask is not None and bool(mask[i])) else v
                        for i, v in enumerate(values)]
            return None
        values = cls._h5_strings(node[()])
        legacy = obs.get("__categories")
        if legacy is not None and column in legacy:
            cats = cls._h5_strings(legacy[column][()])
            return [None if int(c) < 0 else cats[int(c)] for c in values]
        return values

    @classmethod
    def _read_h5ad_obs(cls, h5ad_path) -> tuple[list, dict]:
        """Return ``(obs_names, {column: values})`` from a clustered ``.h5ad``.

        Uses anndata (backed, read-only) when importable, else a pure h5py
        reader so the derivation works even in a minimal environment."""
        try:
            import anndata

            adata = anndata.read_h5ad(str(h5ad_path), backed="r")
            try:
                obs = adata.obs
                columns = {col: [None if _is_missing(v) else v
                                 for v in obs[col].tolist()]
                           for col in obs.columns}
                return [str(x) for x in obs.index], columns
            finally:
                file = getattr(adata, "file", None)
                if file is not None:
                    file.close()
        except Exception:
            pass
        import h5py

        with h5py.File(str(h5ad_path), "r") as handle:
            obs = handle["obs"]
            names = cls._h5ad_index(obs)
            columns = {}
            order = obs.attrs.get("column-order", [])
            for col in cls._h5_strings(list(order)) or [
                    k for k in obs.keys() if k not in ("_index", "__categories")]:
                values = cls._h5ad_obs_column(obs, str(col))
                if values is not None and len(values) == len(names):
                    columns[str(col)] = [None if _is_missing(v) else v
                                         for v in values]
            return names, columns

    @classmethod
    def _read_h5ad_var_names(cls, h5ad_path) -> list:
        try:
            import h5py

            with h5py.File(str(h5ad_path), "r") as handle:
                return cls._h5ad_index(handle["var"])
        except Exception:
            pass
        try:
            import anndata

            adata = anndata.read_h5ad(str(h5ad_path), backed="r")
            try:
                return [str(x) for x in adata.var_names]
            finally:
                file = getattr(adata, "file", None)
                if file is not None:
                    file.close()
        except Exception:
            return []

    def _derive_peaks_bed_from_h5ad(self, clustered, dest) -> Optional[str]:
        """Write a clean 3-column BED from the matrix's peak ``var_names``.

        Accepts ``chr:start-end``, ``chr_start_end`` and ``chr-start-end`` with
        non-negative integer ``start < end``; anything else is skipped (never
        guessed). Returns the BED path when at least one interval was written."""
        import re

        pattern = re.compile(r"^(.+?)[:_-](\d+)[-_:](\d+)$")
        intervals = []
        for name in self._read_h5ad_var_names(clustered):
            match = pattern.match(str(name).strip())
            if not match:
                continue
            chrom, start, end = match.group(1), int(match.group(2)), int(match.group(3))
            if chrom and 0 <= start < end:
                intervals.append((chrom, start, end))
        if not intervals:
            return None
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "w", encoding="utf-8") as fh:
            for chrom, start, end in intervals:
                fh.write(f"{chrom}\t{start}\t{end}\n")
        return str(dest)

    @staticmethod
    def _first_contrast(comparisons) -> tuple[Optional[str], Optional[str]]:
        for comp in comparisons or []:
            if isinstance(comp, dict):
                test, ref = comp.get("test"), comp.get("reference")
            elif isinstance(comp, (list, tuple)) and len(comp) == 2:
                test, ref = comp
            else:
                continue
            if test is not None and ref is not None:
                return str(test), str(ref)
        return None, None

    def _derive_scatac_barcode_design(self, exp_ctx: dict, clustered,
                                      dest) -> Optional[tuple]:
        """Write ``barcode<TAB>group[<TAB>replicate][<TAB>batch]`` from ``obs``.

        The grouping column is the explicit ``footprint_groupby`` or the
        confirmed design condition (never transferred scRNA labels or Leiden
        clusters by default), and the contrast must be explicit — ARIA does not
        guess ``(group_a, group_b)`` from a two-level column. Returns
        ``(tsv_path, group_a, group_b)`` or ``None`` when no valid grouping
        column exists."""
        contract = self._resolve_da_design(exp_ctx)
        try:
            barcodes, columns = self._read_h5ad_obs(clustered)
        except Exception as exc:
            log.warning("footprint barcode-design derivation failed: %s", exc)
            return None
        if not barcodes or not columns:
            return None

        explicit = exp_ctx.get("footprint_groupby")
        group_col = None
        if explicit:
            group_col = str(explicit) if str(explicit) in columns else None
        else:
            for cand in (contract.get("condition_col"), exp_ctx.get("groupby"),
                         "condition"):
                if (cand and str(cand) in columns
                        and str(cand) not in self._NON_INFERENTIAL_GROUP_COLS):
                    group_col = str(cand)
                    break
        if group_col is None:
            return None

        rep_col = None
        for cand in (contract.get("replicate_col"), exp_ctx.get("replicate_col"),
                     "donor", "replicate", "sample_id", "sample"):
            if cand and str(cand) in columns and str(cand) != group_col:
                rep_col = str(cand)
                break

        batch_col = None
        batch_candidates = list(contract.get("covariates") or [])
        if exp_ctx.get("batch_key"):
            batch_candidates.append(exp_ctx["batch_key"])
        batch_candidates.append("batch")
        for cand in batch_candidates:
            if (cand and str(cand) in columns
                    and str(cand) not in (group_col, rep_col)):
                batch_col = str(cand)
                break

        group_a = exp_ctx.get("group_a")
        group_b = exp_ctx.get("group_b")
        if not (group_a and group_b):
            group_a, group_b = self._first_contrast(contract.get("comparisons"))
        if not (group_a and group_b):
            group_a, group_b = self._first_contrast(exp_ctx.get("comparisons"))

        def _cell(values, idx):
            value = values[idx] if values is not None else None
            return "" if value is None else str(value).strip()

        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        header = ["barcode", "group"]
        if rep_col:
            header.append("replicate")
        if batch_col:
            header.append("batch")
        n_written = 0
        with open(dest, "w", encoding="utf-8") as fh:
            fh.write("\t".join(header) + "\n")
            for idx, barcode in enumerate(barcodes):
                group = _cell(columns[group_col], idx)
                if not group or not str(barcode).strip():
                    continue
                row = [str(barcode).strip(), group]
                if rep_col:
                    row.append(_cell(columns[rep_col], idx))
                if batch_col:
                    row.append(_cell(columns[batch_col], idx))
                fh.write("\t".join(row) + "\n")
                n_written += 1
        if n_written == 0:
            return None
        return (str(dest),
                str(group_a) if group_a else None,
                str(group_b) if group_b else None)

    def _rna_group_means_from_path(self, exp_ctx: dict, genes: list,
                                   group_a: str, group_b: str) -> Optional[dict]:
        """Mean RNA expression of candidate TF genes per footprint group, read
        from the paired scRNA ``.h5ad`` / ``.h5mu`` (``rna`` modality). Returns
        ``None`` whenever the data or matching group labels are unavailable."""
        path = exp_ctx.get("rna_data_path") or exp_ctx.get("scrna_data_path")
        if not path or not Path(str(path)).is_file() or not genes:
            return None
        import numpy as np

        if str(path).lower().endswith(".h5mu"):
            import mudata

            adata = mudata.read_h5mu(str(path))["rna"]
        else:
            import anndata

            adata = anndata.read_h5ad(str(path))
        label_col = next((c for c in (exp_ctx.get("rna_label_col"),
                                      exp_ctx.get("condition_col"),
                                      "condition", "cell_type")
                          if c and c in adata.obs.columns), None)
        if label_col is None:
            return None
        labels = adata.obs[label_col].astype(str).to_numpy()
        mask_a, mask_b = labels == str(group_a), labels == str(group_b)
        if not mask_a.any() or not mask_b.any():
            return None
        lookup = {str(v).upper(): i for i, v in enumerate(adata.var_names)}
        means: dict = {}
        for gene in genes:
            idx = lookup.get(str(gene).upper())
            if idx is None:
                continue
            column = adata.X[:, idx]
            column = (np.asarray(column.toarray()).ravel()
                      if hasattr(column, "toarray") else np.asarray(column).ravel())
            means[str(gene).upper()] = {
                str(group_a): float(column[mask_a].mean()),
                str(group_b): float(column[mask_b].mean()),
            }
        return means or None

    def _attach_footprint_rna_cross_evidence(self, result: dict,
                                             exp_ctx: dict) -> None:
        """Annotate top differential footprints with paired-RNA concordance.

        Associative cross-modal support only; any failure (no RNA, mismatched
        labels, missing readers) leaves footprinting untouched."""
        try:
            if not (isinstance(result, dict) and result.get("ran")
                    and isinstance(result.get("differential_summary"), dict)
                    and result["differential_summary"].get("parsed")):
                return
            ga, gb = result.get("group_a"), result.get("group_b")
            if not ga or not gb:
                return
            summary = result["differential_summary"]
            rna_means = exp_ctx.get("rna_group_means")
            if not isinstance(rna_means, dict):
                from aria.agents.narrative.synthesis.footprint_rna import _tf_genes

                genes = []
                for favored in (ga, gb):
                    for row in summary.get(f"top_toward_{favored}") or []:
                        for gene in _tf_genes((row or {}).get("tf", "")):
                            if gene not in genes:
                                genes.append(gene)
                try:
                    rna_means = self._rna_group_means_from_path(
                        exp_ctx, genes, str(ga), str(gb))
                except Exception as exc:
                    log.info("footprint RNA cross-evidence unavailable: %s", exc)
                    rna_means = None
            if isinstance(rna_means, dict) and rna_means:
                from aria.agents.narrative.synthesis.footprint_rna import (
                    footprint_rna_concordance,
                )

                result["rna_cross_evidence"] = footprint_rna_concordance(
                    summary, rna_means, str(ga), str(gb))
        except Exception as exc:  # cross-evidence must never fail footprinting
            log.warning("footprint RNA cross-evidence failed: %s", exc)

    def _bulk_condition_bams(self, exp_ctx: dict, files: list,
                             condition_col: str = "condition") -> dict:
        condition_col = condition_col or "condition"
        provided = exp_ctx.get("condition_bams")
        if isinstance(provided, dict):
            return provided
        sample_metadata = exp_ctx.get("sample_metadata") or {}
        if not isinstance(sample_metadata, dict):
            return {}
        sample_ids = exp_ctx.get("sample_ids")
        if not sample_ids or len(sample_ids) != len(files):
            from aria.scripts.chromatin_peak_counts import _strip_known_suffixes
            sample_ids = [_strip_known_suffixes(f) for f in files]
        out: dict[str, list[str]] = {}
        for sample_id, bam in zip(sample_ids, files):
            meta = sample_metadata.get(sample_id) or {}
            if not isinstance(meta, dict):
                continue
            condition = meta.get(condition_col) or meta.get("condition")
            if condition:
                out.setdefault(str(condition), []).append(str(bam))
        return out

    def _bulk_replicate_bams(self, exp_ctx: dict, files: list,
                             condition_col: str = "condition",
                             replicate_col: str = "replicate") -> dict:
        """Build an explicit biological-replicate design without filename guesses.

        ``condition_col`` / ``replicate_col`` come from the shared design
        contract (E7) so a confirmed factor named otherwise groups correctly."""
        provided = exp_ctx.get("replicate_bams")
        if isinstance(provided, dict):
            return provided
        sample_metadata = exp_ctx.get("sample_metadata") or {}
        sample_ids = exp_ctx.get("sample_ids")
        if (not isinstance(sample_metadata, dict) or not sample_ids
                or len(sample_ids) != len(files)):
            return {}
        condition_col = condition_col or "condition"
        replicate_col = replicate_col or "replicate"
        out: dict[str, dict[str, list[str]]] = {}
        for sample_id, bam in zip(sample_ids, files):
            meta = sample_metadata.get(sample_id) or {}
            if not isinstance(meta, dict):
                continue
            condition = meta.get(condition_col) or meta.get("condition")
            replicate = meta.get(replicate_col)
            if replicate is None:
                replicate = meta.get("replicate")
            if condition is None or replicate is None:
                continue
            out.setdefault(str(condition), {}).setdefault(
                str(replicate), []).append(str(bam))
        return out

    def _bulk_replicate_batches(self, exp_ctx: dict, files: list,
                                condition_col: str = "condition",
                                replicate_col: str = "replicate",
                                batch_col: Optional[str] = None) -> dict:
        """``{condition: {replicate: batch}}`` from explicit sample metadata.

        Same strict identity rule as ``_bulk_replicate_bams`` (explicit
        ``sample_ids`` aligned 1:1 with ``files``; no filename guesses). A
        biological replicate is batched only when every one of its technical
        libraries agrees on one non-empty batch; otherwise it is left unadjusted
        and the backend keeps the unadjusted test."""
        provided = exp_ctx.get("replicate_batches")
        if isinstance(provided, dict):
            return provided
        if not batch_col:
            batch_col = next((str(c) for c in (
                *(exp_ctx.get("covariates") or []),
                exp_ctx.get("batch_col"), exp_ctx.get("batch_key"), "batch")
                if c and str(c).strip()), None)
        sample_metadata = exp_ctx.get("sample_metadata") or {}
        sample_ids = exp_ctx.get("sample_ids")
        if (not batch_col or not isinstance(sample_metadata, dict)
                or not sample_ids or len(sample_ids) != len(files)):
            return {}
        condition_col = condition_col or "condition"
        replicate_col = replicate_col or "replicate"
        seen: dict[tuple[str, str], set] = {}
        for sample_id in sample_ids:
            meta = sample_metadata.get(sample_id) or {}
            if not isinstance(meta, dict):
                continue
            condition = meta.get(condition_col) or meta.get("condition")
            replicate = meta.get(replicate_col)
            if replicate is None:
                replicate = meta.get("replicate")
            if condition is None or replicate is None:
                continue
            batch = meta.get(batch_col)
            batch = "" if batch is None else str(batch).strip()
            seen.setdefault((str(condition), str(replicate)), set()).add(batch)
        out: dict[str, dict[str, str]] = {}
        for (condition, replicate), batches in seen.items():
            if len(batches) == 1 and "" not in batches:
                out.setdefault(condition, {})[replicate] = next(iter(batches))
        return out

    def _run_bulk_atac_footprinting(self, experiment_id: str, exp_ctx: dict,
                                    da_result: dict, peaks_result: dict,
                                    files: list,
                                    contract: Optional[dict] = None,
                                    ) -> Optional[dict]:
        """Inline bulk ATAC condition-level TOBIAS dispatch (C2; production
        per ADR-078)."""
        comparisons = [
            c for c in (da_result.get("comparisons") or [])
            if isinstance(c, dict) and c.get("status") == "success"
        ]
        comp = comparisons[0] if comparisons else {}
        if contract is None:
            contract = self._resolve_da_design(exp_ctx)
        condition_col = (contract.get("condition_col")
                         or exp_ctx.get("condition_col") or "condition")
        replicate_col = (contract.get("replicate_col")
                         or exp_ctx.get("replicate_col") or "replicate")
        params = {
            **self._footprint_common_params(
                exp_ctx, str(Path(files[0]).parent / "bulk_atac_footprinting")),
            "mode": "bulk",
            "data_type": "bulk_ATAC",
            "condition_bams": self._bulk_condition_bams(
                exp_ctx, files, condition_col=condition_col),
            "replicate_bams": self._bulk_replicate_bams(
                exp_ctx, files, condition_col=condition_col,
                replicate_col=replicate_col),
            "group_a": (
                exp_ctx.get("footprint_group_a")
                or exp_ctx.get("tobias_group_a")
                or comp.get("test")
            ),
            "group_b": (
                exp_ctx.get("footprint_group_b")
                or exp_ctx.get("tobias_group_b")
                or comp.get("reference")
            ),
            "peaks_bed": (
                exp_ctx.get("peaks_bed")
                or exp_ctx.get("peaks_path")
                or peaks_result.get("consensus_peaks_path")
                or peaks_result.get("peaks_path")
            ),
            "replicate_batches": self._bulk_replicate_batches(
                exp_ctx, files, condition_col=condition_col,
                replicate_col=replicate_col,
                batch_col=(contract.get("batch_covariate")
                           or (contract.get("covariates") or [None])[0])),
        }
        if bool(exp_ctx.get("paired_replicates")
                or (exp_ctx.get("design") or {}).get("paired")):
            params["paired_replicates"] = True
        motif_source = params.pop("_motif_source", None)
        if any(str(f).lower().endswith(".cram") for f in files):
            reference_fasta = str(exp_ctx.get("reference_fasta")
                                  or params.get("genome_fasta") or "")
            if reference_fasta:
                params["reference_fasta"] = reference_fasta
        self.publish_status(
            experiment_id, "bulk ATAC TOBIAS footprinting...", 0.99)
        result = self.env.run_in_stack(
            stack="tobias",
            script_path="aria/scripts/chromatin_footprint_tobias.py",
            params=params,
        )
        if isinstance(result, dict):
            result.setdefault("data_type", "bulk_ATAC")
            result.setdefault("analysis", "differential_tf_footprinting")
            result["validation_level"] = "production"
            if motif_source and result.get("ran"):
                result["motif_source"] = motif_source
            self._attach_footprint_rna_cross_evidence(result, exp_ctx)
        return result

    def _run_bulk_atac_motifs(self, experiment_id: str, exp_ctx: dict,
                              da_result: dict, files: list) -> Optional[dict]:
        """TF motif enrichment over bulk ATAC DA peak sets.

        Builds direction-split DA region groups + the tested-peak background
        from the DESeq2 comparisons and reuses ``chromatin_motifs.py`` (snapatac2,
        chromatin stack) for the actual enrichment. The script self-gates on a
        local genome FASTA + a versioned motif collection (honest skip, never a
        fabricated enrichment). Returns ``None`` only when there are no DA peaks
        to interpret, so the caller leaves ``findings["motifs"]`` unset.
        """
        regions, background, warnings = _bulk_da_motif_regions(
            da_result.get("comparisons"))
        if not regions:
            return None

        assembly = exp_ctx.get("genome")
        genome_fasta = exp_ctx.get("genome_fasta")
        out_dir = str(Path(files[0]).parent / "bulk_atac_motifs")
        self.publish_status(
            experiment_id, "bulk ATAC TF motif enrichment...", 0.9)
        motif_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_motifs.py",
            params={
                "regions": regions,
                "background": background,
                "genome_fasta": genome_fasta,
                "assembly": assembly,
                "allow_genome_fetch": bool(
                    exp_ctx.get("allow_genome_fetch", False)),
                "motif_collection": exp_ctx.get("motif_collection"),
                "method": "hypergeometric",
                "output_dir": out_dir,
                "exp_context": exp_ctx,
                # F7: resolve the threshold policy by the REAL modality, not the
                # hardcoded scATAC default the engine assumes.
                "modality": "bulk_ATAC",
                "foreground_truncation_strategy": (
                    "rank_by_padj_then_abs_log2fc_before_cap"),
            },
        )
        if isinstance(motif_result, dict):
            if warnings:
                motif_result.setdefault("warnings", []).extend(warnings)
            motif_result.setdefault("data_type", "bulk_ATAC")
            if motif_result.get("ran"):
                motif_result["validation_level"] = "production"
        return motif_result

    def _run_bulk_atac_annotation(self, experiment_id: str, exp_ctx: dict,
                                  da_result: dict,
                                  files: list) -> Optional[dict]:
        """Genomic annotation of bulk ATAC DA peaks (B2).

        Annotates each comparison's significant peaks with the nearest gene,
        signed distance-to-TSS, and a feature class (Promoter/Exonic/Intronic/
        Distal Intergenic) over an auto-resolved GTF (``~/.aria/genomes/*``,
        ``exp_ctx['gtf']`` hint, or ``ARIA_GTF``). Runs in the ``rna`` stack
        (pure pandas, no chromatin deps, same place the DA CSVs were produced).
        Honest-skip when there are no DA comparisons; the script self-skips when
        no GTF is locatable. Returns ``None`` only when nothing was tested."""
        comparisons = da_result.get("comparisons") or []
        if not comparisons:
            return None
        self.publish_status(
            experiment_id, "bulk ATAC peak annotation...", 0.95)
        ann_result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/chromatin_peak_annotation.py",
            params={
                "data_type": "bulk_ATAC",
                "comparisons": comparisons,
                "genome": exp_ctx.get("genome"),
                "gtf": (exp_ctx.get("gtf") or exp_ctx.get("gtf_path")
                        or exp_ctx.get("gtf_file")),
                "promoter_upstream": exp_ctx.get("promoter_upstream"),
                "promoter_downstream": exp_ctx.get("promoter_downstream"),
                "output_dir": str(Path(files[0]).parent / "bulk_atac_annotation"),
            },
        )
        if isinstance(ann_result, dict):
            ann_result.setdefault("data_type", "bulk_ATAC")
            if ann_result.get("status") == "success":
                ann_result["validation_level"] = "production"
        return ann_result

    def _run_bulk_atac_ora(self, experiment_id: str, exp_ctx: dict,
                           da_result: dict, files: list) -> Optional[dict]:
        """Functional ORA of genes near bulk ATAC DA peaks (B3).

        Reuses the B2 peak→gene assignment + the bulk RNA ORA engine (local
        hypergeometric over versioned GMTs). Runs in the ``rna`` stack (where the
        ORA engine + its GMTs + seaborn live). Honest-skip when there are no DA
        comparisons; the script self-skips without a GTF or a provisioned GMT.
        Returns ``None`` only when nothing was tested."""
        comparisons = da_result.get("comparisons") or []
        if not comparisons:
            return None
        self.publish_status(
            experiment_id, "bulk ATAC pathway ORA...", 0.97)
        ora_result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/chromatin_peak_ora.py",
            params={
                "data_type": "bulk_ATAC",
                "comparisons": comparisons,
                "genome": exp_ctx.get("genome"),
                "organism": exp_ctx.get("organism"),
                "gtf": exp_ctx.get("gtf") or exp_ctx.get("gtf_path"),
                "output_dir": str(Path(files[0]).parent / "bulk_atac_ora"),
            },
        )
        if isinstance(ora_result, dict):
            ora_result.setdefault("data_type", "bulk_ATAC")
            if ora_result.get("status") == "success":
                ora_result["validation_level"] = "production"
        return ora_result

    # ── ChIP-seq pipeline ─────────────────────────────────────────────────

    def _run_chip(self, experiment_id: str, exp_ctx: dict,
                  intent: dict, files: list) -> dict:
        """
        ChIP-seq: QC, peak calling with input control, IDR for replicates.
        Separates histone mark (broad peaks) from TF binding (narrow peaks).
        """
        findings = {}

        # Classify: histone mark vs TF binding
        assay_class = self._classify_chip_target(files, intent)
        findings["assay_class"] = assay_class

        # QC
        qc_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_qc.py",
            params={
                "data_type":   "ChIP",
                "files":       files,
                "genome":      exp_ctx.get("genome", "hg38"),
                "assay_class": assay_class,
            },
        )
        findings["qc"] = qc_result

        # Peak calling — narrow for TF, broad for histones
        macs3_p = dict(self.MACS3_PARAMS["ChIP"])
        if assay_class == "histone":
            macs3_p["broad"] = True
        else:
            macs3_p["broad"] = False

        peaks_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_peaks.py",
            params={
                "data_type":    "ChIP",
                "files":        files,
                "genome":       exp_ctx.get("genome", "hg38"),
                "macs3_params": macs3_p,
                "output_dir":   str(Path(files[0]).parent / "peaks"),
            },
        )
        findings["peaks"] = peaks_result

        if peaks_result.get("status") == "success":
            self._publish_peaks_finding(experiment_id, peaks_result, "ChIP")

        return {"status": "done", "findings": findings}

    # ── CUT&RUN / CUT&TAG pipeline ────────────────────────────────────────

    def _run_cut_and_run(self, experiment_id: str, exp_ctx: dict,
                          intent: dict, files: list,
                          assay_type: str = "CUT_AND_RUN") -> dict:
        """
        CUT&RUN / CUT&TAG pipeline.
        Key differences from ChIP: very low background, --nolambda for MACS3,
        fragment sizes differ (nucleosome-protected ~200bp, TF ~100bp).
        """
        findings = {}

        qc_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_qc.py",
            params={
                "data_type": assay_type,
                "files":     files,
                "genome":    exp_ctx.get("genome", "hg38"),
            },
        )
        findings["qc"] = qc_result

        peaks_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_peaks.py",
            params={
                "data_type":    assay_type,
                "files":        files,
                "genome":       exp_ctx.get("genome", "hg38"),
                "macs3_params": self.MACS3_PARAMS[assay_type],
                "output_dir":   str(Path(files[0]).parent / "peaks"),
            },
        )
        findings["peaks"] = peaks_result

        return {"status": "done", "findings": findings}

    # ── ParameterAdvisor for LSI ──────────────────────────────────────────

    def _advise_lsi_params(self, experiment_id: str,
                            intent: dict,
                            qc_result: dict) -> dict:
        """
        Use ParameterAdvisor to decide LSI parameters.

        Key decision: how many SVD components to use.
        Component 1 is ALWAYS discarded (depth correlation).
        Typical range: use components 2-30 for most datasets.

        The ParameterAdvisor stores this decision in ARIAMemory so
        future experiments with similar cell types learn from it.
        """
        n_cells = qc_result.get("n_cells_after", 5000)

        # Biological context informs the component range
        # More heterogeneous tissue → more components needed
        complexity = intent.get("complexity", "moderate")
        if complexity == "complex" or n_cells > 20000:
            n_components_range = (30, 60)
        elif complexity == "simple" or n_cells < 2000:
            n_components_range = (15, 30)
        else:
            n_components_range = (20, 40)

        # Store decision with mandatory note about discarding component 1
        decision = {
            "method":            "LSI (TF-IDF + SVD)",
            "components_range":  n_components_range,
            "discard_component": 1,
            "rationale": (
                f"Component 1 discarded: correlates with sequencing depth "
                f"(technical artifact, not biology). "
                f"Range {n_components_range} based on dataset complexity "
                f"({n_cells} cells, {complexity} experiment)."
            ),
            "n_cells": n_cells,
        }

        # Store in memory for future reference
        room_id = f"{experiment_id}_scATAC"
        try:
            self.memory.create_room(
                room_id=f"{room_id}_lsi",
                hall_id=f"{experiment_id}_scATAC",
                analysis="lsi_dimensionality_reduction",
                params=decision,
                tool="TF-IDF + SVD (episcanpy/muon)",
            )
            self.memory.store_finding(
                finding_id=str(uuid.uuid4())[:8],
                room_id=f"{room_id}_lsi",
                content=json.dumps(decision),
                confidence="high",
            )
        except Exception as e:
            log.warning(f"Could not store LSI decision: {e}")

        return decision

    # ── TF motif enrichment with DebateCouncil ────────────────────────────

    def _run_motif_enrichment(self, experiment_id: str, exp_ctx: dict,
                               intent: dict, peaks_result: dict,
                               assay_type: str) -> dict:
        """
        TF motif enrichment analysis.
        DebateCouncil is MANDATORY here — Tn5 bias must be addressed
        before claiming any TF is active based on motif enrichment alone.
        """
        script_path = "aria/scripts/chromatin_motifs.py"
        if not self._script_exists(script_path):
            return self._planned_script_blocker(script_path, "motif_enrichment")

        motif_result = self.env.run_in_stack(
            stack="chromatin",
            script_path=script_path,
            params={
                "peaks_path": peaks_result.get("peaks_path", ""),
                "genome":     exp_ctx.get("genome", "hg38"),
                "organism":   exp_ctx.get("organism", "Homo sapiens"),
            },
        )

        if motif_result.get("status") != "success":
            return motif_result

        # DebateCouncil review — always for motif claims
        # Critic must check Tn5 bias before accepting TF activity claims
        try:
            from aria.agents.debate_council import DebateCouncil
            council = DebateCouncil(llm=self.llm, max_rounds=2)

            top_motifs = motif_result.get("top_motifs", [])[:5]
            claim = (
                f"The following transcription factors show significant motif "
                f"enrichment in accessible chromatin regions: "
                f"{', '.join(top_motifs)}"
            )

            debate_result = council.resolve(
                topic="tf_motif_enrichment",
                initial_claim=claim,
                evidence={
                    "top_motifs":     top_motifs,
                    "enrichment_scores": motif_result.get("scores", {}),
                    "assay_type":     assay_type,
                    "n_peaks":        motif_result.get("n_peaks", 0),
                    "tn5_bias_corrected": motif_result.get(
                        "tn5_bias_corrected", False
                    ),
                },
                biological_context=intent,
            )

            motif_result["debate_verdict"]    = debate_result.verdict.value
            motif_result["debate_consensus"]  = debate_result.consensus
            motif_result["debate_limitations"] = debate_result.limitations
            motif_result["confidence"]         = debate_result.confidence

            # Publish finding with debate-informed confidence
            conf_map = {
                "high":         Confidence.HIGH,
                "medium":       Confidence.MEDIUM,
                "low":          Confidence.LOW,
                "insufficient": Confidence.INSUFFICIENT,
            }
            self.publish_finding(
                experiment_id,
                {"summary": debate_result.consensus,
                 "top_motifs": top_motifs,
                 "verdict": debate_result.verdict.value},
                conf_map.get(debate_result.confidence, Confidence.MEDIUM),
            )

        except Exception as e:
            log.warning(f"DebateCouncil unavailable for motif review: {e}")
            # Publish with medium confidence if debate failed
            self.publish_finding(
                experiment_id,
                {"summary": f"Motif enrichment: {top_motifs}",
                 "note": "DebateCouncil review pending"},
                Confidence.MEDIUM,
            )

        return motif_result

    # ── Differential accessibility ────────────────────────────────────────

    def _run_differential_accessibility(
        self, experiment_id: str, exp_ctx: dict,
        intent: dict, peaks_result: dict
    ) -> dict:
        """Differential accessibility between conditions using DESeq2."""
        script_path = "aria/scripts/chromatin_differential.py"
        if not self._script_exists(script_path):
            return self._planned_script_blocker(
                script_path, "differential_accessibility"
            )

        return self.env.run_in_stack(
            stack="chromatin",
            script_path=script_path,
            params={
                "peaks_path":  peaks_result.get("consensus_peaks_path", ""),
                "comparison":  intent.get("comparison", ""),
                "genome":      exp_ctx.get("genome", "hg38"),
                "organism":    exp_ctx.get("organism", "Homo sapiens"),
            },
        )

    # ── Helpers ───────────────────────────────────────────────────────────

    def _classify_chip_target(self, files: list, intent: dict) -> str:
        """
        Classify ChIP target as 'histone' or 'tf' based on filenames
        and biological context. Determines peak calling mode.
        """
        HISTONE_MARKS = {
            "h3k4me3", "h3k4me1", "h3k27ac", "h3k27me3",
            "h3k9me3", "h3k36me3", "h3k9ac", "h4k20me1",
        }
        all_text = " ".join(files + [intent.get("user_question", "")]).lower()

        if any(mark in all_text for mark in HISTONE_MARKS):
            return "histone"

        # Default to TF if no histone marks detected
        return "tf"

    def _needs_tf_analysis(self, intent: dict) -> bool:
        """Check if the biological question requires TF analysis."""
        TF_KEYWORDS = [
            "transcription factor", "tf ", "motif", "footprint",
            "regulatory", "binding", "enhancer", "promoter",
        ]
        question = intent.get("user_question", "").lower()
        return any(kw in question for kw in TF_KEYWORDS)

    @staticmethod
    def _script_exists(script_path: str) -> bool:
        root = Path(__file__).resolve().parents[2]
        return (root / script_path).exists()

    @staticmethod
    def _planned_script_blocker(script_path: str, analysis: str) -> dict:
        return {
            "status": "skipped",
            "reason": "script_not_implemented",
            "analysis": analysis,
            "script_path": script_path,
            "validation_level": "scaffold",
            "details": (
                f"{script_path} is planned for the v4.6+ chromatin roadmap "
                "and is not available in this build."
            ),
        }

    def _publish_qc_finding(self, experiment_id: str,
                             qc_result: dict, modality: str):
        warnings = qc_result.get("warnings", [])
        conf = Confidence.HIGH if not warnings else Confidence.MEDIUM

        # FRiP/TSS are honestly None for a pre-called peak matrix (.h5mu) until a
        # reference TSS annotation / called peaks exist (ADR-002 / B9). Only flag
        # low quality when they were actually computed — never fabricate a value
        # or crash on the None contract.
        frip      = qc_result.get("frip")
        tss_score = qc_result.get("tss_enrichment")

        if (frip is not None and frip < 0.2) or \
           (tss_score is not None and tss_score < 4):
            conf = Confidence.LOW

        def _fmt(value, fmt: str) -> str:
            return format(value, fmt) if isinstance(value, (int, float)) \
                else "not computed"

        summary = f"{modality} QC: FRiP={_fmt(frip, '.3f')}, " \
                  f"TSS={_fmt(tss_score, '.2f')}"
        n_cells = qc_result.get("n_cells")
        n_peaks = qc_result.get("n_peaks")
        if isinstance(n_cells, int) and isinstance(n_peaks, int):
            summary = (f"{modality} QC: {n_cells:,} cells x {n_peaks:,} peaks; "
                       f"FRiP={_fmt(frip, '.3f')}, TSS={_fmt(tss_score, '.2f')}")

        self.publish_finding(
            experiment_id,
            {"summary":    summary,
             "modality":   modality,
             "frip":       frip,
             "tss_score":  tss_score,
             "n_cells":    n_cells,
             "n_peaks":    n_peaks,
             "warnings":   warnings},
            conf,
        )

    def _publish_peaks_finding(self, experiment_id: str,
                                peaks_result: dict, modality: str):
        n_peaks = peaks_result.get("n_peaks", 0)
        conf    = (Confidence.HIGH   if n_peaks > 50000 else
                   Confidence.MEDIUM if n_peaks > 10000 else
                   Confidence.LOW)
        self.publish_finding(
            experiment_id,
            {"summary":  f"{modality} peaks: {n_peaks:,} peaks called",
             "modality": modality,
             "n_peaks":  n_peaks,
             "peaks_path": peaks_result.get("peaks_path", "")},
            conf,
        )

    def _format_qc_checkpoint(self, qc_result: dict,
                                modality: str) -> str:
        lines = [f"{modality} QC completed with warnings:\n"]
        for w in qc_result.get("warnings", []):
            lines.append(f"  * {w}")
        frip = qc_result.get("frip", "N/A")
        tss  = qc_result.get("tss_enrichment", "N/A")
        lines.append(f"\n  FRiP score:      {frip}")
        lines.append(f"  TSS enrichment:  {tss}")
        lines.append(
            "\n  FRiP < 0.2 or TSS < 4 indicates poor library quality."
        )
        lines.append("How would you like to proceed?")
        return "\n".join(lines)

    def receive(self, message):
        pass
