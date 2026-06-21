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
        h5mu = next((f for f in files
                     if str(f).lower().endswith(".h5mu")), None)
        if h5mu:
            return self._run_scatac_matrix(experiment_id, exp_ctx, intent, h5mu)

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
            matrix_result = self._run_scatac_matrix(
                experiment_id, exp_ctx, intent,
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
        frag = next(
            (f for f in files
             if "fragment" in str(f).lower()
             and str(f).lower().endswith((".tsv", ".tsv.gz", ".bed", ".bed.gz"))),
            None,
        )
        if frag is None:
            # A single tsv/bed.gz with no "fragment" token: accept only when it
            # is the lone input (avoid grabbing peak/annotation BEDs).
            tsvs = [f for f in files
                    if str(f).lower().endswith((".tsv", ".tsv.gz"))]
            frag = tsvs[0] if len(tsvs) == 1 else None
        if frag is None:
            return {"status": "skipped", "ran": False,
                    "analysis": "fragments_to_peak_matrix",
                    "reason": "no_fragments_file"}

        peak_file = None
        if isinstance(peaks_result, dict) and peaks_result.get("status") == "success":
            peak_file = (peaks_result.get("consensus_peaks_path")
                         or peaks_result.get("peaks_path"))
        return self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_fragments_to_matrix.py",
            params={
                "fragments_file": frag,
                "genome": exp_ctx.get("genome", "hg38"),
                "peak_file": peak_file,
                "output_dir": str(Path(frag).parent / "aria_chromatin"),
                "min_fragments": int(exp_ctx.get("min_fragments", 1000)),
                "resolution": exp_ctx.get("leiden_resolution", 1.0),
            },
        )

    def _run_scatac_fastq_to_fragments(self, experiment_id: str, exp_ctx: dict,
                                       files: list) -> dict:
        """Dispatch the chromap scATAC FASTQ->fragments front-end (atacseq stack).

        The barcode read cannot be inferred from the genomic FASTQs, so the
        barcode read + whitelist must be explicit (exp_context); otherwise this
        honestly skips and the raw FASTQ goes no further. The genomic mates fall
        back to positional R1/R3 detection among the inputs.
        """
        fastqs = [f for f in files if _is_fastq([f])]
        r1 = exp_ctx.get("r1_fastq") or _pick_read(fastqs, ("_r1", "_R1"))
        r3 = (exp_ctx.get("r3_fastq") or exp_ctx.get("r2_fastq")
              or _pick_read(fastqs, ("_r3", "_R3")))
        barcode = (exp_ctx.get("barcode_fastq") or exp_ctx.get("cell_barcode_fastq"))
        whitelist = (exp_ctx.get("barcode_whitelist")
                     or exp_ctx.get("cell_barcode_whitelist"))

        if not barcode or not whitelist:
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
        return self.env.run_in_stack(
            stack="atacseq",
            script_path="aria/scripts/chromatin_scatac_align.py",
            params={
                "r1_fastq": r1, "r3_fastq": r3, "barcode_fastq": barcode,
                "genome_fasta": genome_fasta,
                "chromap_index": exp_ctx.get("chromap_index"),
                "barcode_whitelist": whitelist,
                "output_dir": str(Path((r1 or barcode)).parent
                                  / "aria_scatac_aligned"),
                "threads": int(exp_ctx.get("threads", 8)),
            },
        )

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
        design = exp_ctx.get("design") or {}
        pseudobulk = design.get("pseudobulk") or {}
        condition_col = (
            exp_ctx.get("condition_col")
            or pseudobulk.get("condition_col")
            or design.get("condition_col")
            or design.get("main_factor")
        )
        replicate_col = (
            exp_ctx.get("replicate_col")
            or pseudobulk.get("replicate_col")
            or design.get("replicate_col")
        )
        batch_col = (
            exp_ctx.get("batch_col")
            or exp_ctx.get("batch_covariate")
            or design.get("batch_col")
            or design.get("batch_factor")
            or design.get("batch_covariate")
        )
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
            # Forward confirmed cross-condition design only when present; absent
            # metadata makes the pseudobulk lane honestly skip (P0-5).
            for k in ("condition_col", "replicate_col", "comparisons"):
                if exp_ctx.get(k) is not None:
                    da_params[k] = exp_ctx[k]
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

            # 4. TF motif enrichment in the DA peak sets (offline, self-gating)
            da_csv = ((da_result.get("per_cluster") or {}).get("output_csv")
                      if da_result.get("status") == "success" else None)

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
                },
            )
            findings["regulatory"] = regulatory_result

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
                "analysis": "bulk_atac_alignment", "validation_level": "beta",
                "reason": "no_reference_genome_fasta",
                "message": (
                    f"Bulk ATAC FASTQ alignment needs a reference genome FASTA "
                    f"for assembly '{assembly}', which is not staged locally. "
                    f"Provide exp_context['genome_fasta'] or stage "
                    f"~/.aria/genomes/{assembly}/genome.fa."),
            }

        # 1. fastp trim/QC (reuse the validated bulk RNA trimmer).
        self.publish_status(experiment_id, "bulk ATAC FASTQ QC (fastp)...", 0.10)
        qc = self.env.run_in_stack(
            stack="rnaseq",
            script_path="aria/scripts/rna_fastq_qc.py",
            params={"fastq_dir": fastq_dir,
                    "output_dir": str(out_dir / "qc"), "threads": 8},
        )
        if qc.get("status") != "success" or not qc.get("samples"):
            return None, {"status": "failed", "ran": False,
                          "analysis": "bulk_atac_alignment",
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
                          "analysis": "bulk_atac_alignment", "stage": "align",
                          "details": align}

        bams = [b.get("bam") for b in align.get("bam_files", [])
                if b.get("status") == "success" and b.get("bam")]
        if not bams:
            return None, {"status": "failed", "ran": False,
                          "analysis": "bulk_atac_alignment", "stage": "align",
                          "reason": "no_bam_produced", "details": align}

        preprocessing = {
            "status": "success", "ran": True,
            "analysis": "bulk_atac_alignment", "validation_level": "beta",
            "data_type": "bulk_ATAC",
            "aligner": align.get("aligner", "bwa-mem2"),
            "genome_fasta": genome_fasta, "assembly": assembly,
            "filter": align.get("filter"),
            "n_samples_aligned": len(bams),
            "fastq_qc": {"n_samples": qc.get("n_samples"),
                         "multiqc_report": qc.get("multiqc_report")},
            "bam_files": align.get("bam_files"),
        }
        return bams, preprocessing

    def _run_bulk_atac(self, experiment_id: str, exp_ctx: dict,
                       intent: dict, files: list) -> dict:
        """Bulk ATAC-seq V47 slice: QC + MACS3 peak calling.

        Comparison requests build the technical peak-count matrix and then run
        replicate-gated DESeq2 DA only when explicit condition/replicate/
        comparison metadata are present; under-specified designs skip honestly.

        Raw FASTQ inputs are aligned to BAM first (bwa-mem2, ATAC filtering),
        mirroring the bulk RNA FASTQ->counts prepend; the rest of the lane is
        unchanged. scATAC FASTQ is intentionally NOT handled here.
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
            files = bam_files

        # QC
        self.publish_status(experiment_id, "bulk ATAC QC...", 0.35)
        qc_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_qc.py",
            params={
                "data_type": "bulk_ATAC",
                "files":     files,
                "genome":    exp_ctx.get("genome", "hg38"),
            },
        )
        if isinstance(qc_result, dict) and qc_result.get("status") == "success":
            qc_result.setdefault("validation_level", "beta")
        findings["qc"] = qc_result
        if qc_result.get("status") == "success":
            self._publish_qc_finding(experiment_id, qc_result, "bulk_ATAC")

        # Peak calling
        self.publish_status(experiment_id, "bulk ATAC peak calling...", 0.55)
        peaks_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_peaks.py",
            params={
                "data_type":    "bulk_ATAC",
                "files":        files,
                "genome":       exp_ctx.get("genome", "hg38"),
                "macs3_params": self.MACS3_PARAMS["bulk_ATAC"],
                "output_dir":   str(Path(files[0]).parent / "peaks"),
            },
        )
        if isinstance(peaks_result, dict) and peaks_result.get("status") == "success":
            peaks_result.setdefault("validation_level", "beta")
        findings["peaks"] = peaks_result
        if peaks_result.get("status") == "success":
            self._publish_peaks_finding(experiment_id, peaks_result, "bulk_ATAC")

        # Peak-count matrix + replicate-gated DA if comparison defined.
        if intent.get("comparison"):
            peak_universe = (
                peaks_result.get("consensus_peaks_path")
                or peaks_result.get("peaks_path")
            ) if isinstance(peaks_result, dict) else None
            count_result = None
            if peak_universe:
                self.publish_status(
                    experiment_id, "bulk ATAC peak-count matrix...", 0.70)
                count_result = self.env.run_in_stack(
                    stack="chromatin",
                    script_path="aria/scripts/chromatin_peak_counts.py",
                    params={
                        "data_type": "bulk_ATAC",
                        "files": files,
                        "peaks_path": peak_universe,
                        "sample_ids": exp_ctx.get("sample_ids"),
                        "sample_metadata": exp_ctx.get("sample_metadata"),
                        "output_dir": str(Path(files[0]).parent /
                                          "bulk_atac_counts"),
                    },
                )
                findings["peak_counts"] = count_result
            else:
                findings["peak_counts"] = {
                    "status": "skipped",
                    "reason": "no_peak_universe_for_bulk_atac_counts",
                    "analysis": "peak_count_matrix",
                    "validation_level": "scaffold",
                }
            if isinstance(count_result, dict) and count_result.get("status") == "success":
                self.publish_status(
                    experiment_id, "bulk ATAC differential accessibility...",
                    0.82)
                da_params = {
                    "data_type": "bulk_ATAC",
                    "counts_matrix_path": count_result.get("counts_matrix_path"),
                    "sample_metadata_path": count_result.get(
                        "sample_metadata_path"),
                    "condition_col": exp_ctx.get("condition_col", "condition"),
                    "replicate_col": exp_ctx.get(
                        "replicate_col",
                        exp_ctx.get("replicate_column", "replicate"),
                    ),
                    "comparisons": (
                        exp_ctx.get("comparisons")
                        or intent.get("comparisons")
                        or intent.get("comparison")
                    ),
                    "covariates": exp_ctx.get("covariates", []),
                    "exp_context": exp_ctx,
                    "output_dir": str(Path(files[0]).parent / "bulk_atac_da"),
                }
                # Honor an explicit replicate-floor override (e.g. n=2 ENCODE
                # isogenic-replicate designs run with a low_power_warning). Absent
                # an override, the script keeps its production floor.
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
                if (isinstance(da_result, dict)
                        and da_result.get("status") == "success"
                        and da_result.get("ran")):
                    findings["differential_accessibility"] = da_result
                    # TF motif enrichment over the bulk ATAC DA peak sets
                    # (offline, self-gating on genome FASTA + motif collection).
                    # Reuses the scATAC chromatin_motifs engine with explicit
                    # direction-split region groups; both directions reported.
                    motif_result = self._run_bulk_atac_motifs(
                        experiment_id, exp_ctx, da_result, files)
                    if motif_result is not None:
                        findings["motifs"] = motif_result
                    # Genomic annotation of the DA peaks (B2): nearest gene +
                    # distance-to-TSS + feature class over an auto-resolved GTF.
                    # Pure-python (pandas only) → runs in the rna stack alongside
                    # the DA CSVs. Honest-skip without a GTF.
                    ann_result = self._run_bulk_atac_annotation(
                        experiment_id, exp_ctx, da_result, files)
                    if ann_result is not None:
                        findings["peak_annotation"] = ann_result
                    # Functional ORA of genes near the DA peaks (B3): reuses the
                    # bulk RNA ORA engine over the B2 peak->gene assignment. Honest-
                    # skip without a GTF or a provisioned GMT.
                    ora_result = self._run_bulk_atac_ora(
                        experiment_id, exp_ctx, da_result, files)
                    if ora_result is not None:
                        findings["peak_ora"] = ora_result
                    # C2: condition-level TOBIAS footprinting is part of the
                    # live bulk ATAC flow now. The TOBIAS script self-gates on
                    # assets/env and returns ran:false instead of fabricating
                    # uncorrected footprints.
                    footprint_result = self._run_bulk_atac_footprinting(
                        experiment_id, exp_ctx, da_result, peaks_result, files)
                    if footprint_result is not None:
                        findings["footprinting"] = footprint_result
                else:
                    findings["differential_accessibility"] = {
                        "status": "skipped",
                        "reason": (
                            da_result.get("reason")
                            if isinstance(da_result, dict)
                            else "bulk_atac_da_not_run"
                        ),
                        "analysis": "differential_accessibility",
                        "validation_level": "beta",
                        "details": da_result,
                    }
            else:
                findings["differential_accessibility"] = {
                    "status": "skipped",
                    "reason": "bulk_atac_peak_count_matrix_not_available",
                    "analysis": "differential_accessibility",
                    "validation_level": "beta",
                    "details": findings.get("peak_counts"),
                }

        return {"status": "done", "findings": findings}

    def _footprint_common_params(self, exp_ctx: dict, output_dir: str) -> dict:
        return {
            "genome_fasta": exp_ctx.get("genome_fasta"),
            "genome_fai": exp_ctx.get("genome_fai"),
            "peaks_bed": (
                exp_ctx.get("peaks_bed")
                or exp_ctx.get("peaks_path")
                or exp_ctx.get("peak_universe_bed")
            ),
            "motif_meme": (
                exp_ctx.get("motif_meme")
                or exp_ctx.get("motif_meme_path")
                or exp_ctx.get("motif_collection_meme")
            ),
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
        }

    def _run_scatac_footprinting(self, experiment_id: str, exp_ctx: dict,
                                 clustered: str, out_dir: str) -> Optional[dict]:
        """Inline scATAC TOBIAS dispatch (C3).

        The backend refuses uncorrected footprints and returns ran:false for
        missing TOBIAS/assets/groups, so the live flow can attempt footprinting
        without fabricating a result.
        """
        fragments_file = (
            exp_ctx.get("fragments_file")
            or exp_ctx.get("fragments_path")
            or exp_ctx.get("fragment_file")
        )
        barcode_groups = (
            exp_ctx.get("barcode_groups")
            or exp_ctx.get("barcode_groups_tsv")
            or exp_ctx.get("barcode_group_tsv")
        )
        params = {
            **self._footprint_common_params(
                exp_ctx, str(Path(out_dir) / "footprinting")),
            "mode": "scatac",
            "data_type": "scATAC",
            "data_path": clustered,
            "fragments_file": fragments_file,
            "barcode_groups": barcode_groups,
        }
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
            result.setdefault("validation_level", "beta")
        return result

    def _bulk_condition_bams(self, exp_ctx: dict, files: list,
                             condition_col: str = "condition") -> dict:
        provided = exp_ctx.get("condition_bams")
        if isinstance(provided, dict):
            return provided
        sample_metadata = exp_ctx.get("sample_metadata") or {}
        if not isinstance(sample_metadata, dict):
            return {}
        sample_ids = exp_ctx.get("sample_ids")
        if not sample_ids or len(sample_ids) != len(files):
            sample_ids = [Path(f).stem.replace(".bam", "") for f in files]
        out: dict[str, list[str]] = {}
        for sample_id, bam in zip(sample_ids, files):
            meta = sample_metadata.get(sample_id) or {}
            if not isinstance(meta, dict):
                continue
            condition = meta.get(condition_col) or meta.get("condition")
            if condition:
                out.setdefault(str(condition), []).append(str(bam))
        return out

    def _run_bulk_atac_footprinting(self, experiment_id: str, exp_ctx: dict,
                                    da_result: dict, peaks_result: dict,
                                    files: list) -> Optional[dict]:
        """Inline bulk ATAC condition-level TOBIAS dispatch (C2)."""
        comparisons = [
            c for c in (da_result.get("comparisons") or [])
            if isinstance(c, dict) and c.get("status") == "success"
        ]
        comp = comparisons[0] if comparisons else {}
        condition_col = exp_ctx.get("condition_col", "condition")
        params = {
            **self._footprint_common_params(
                exp_ctx, str(Path(files[0]).parent / "bulk_atac_footprinting")),
            "mode": "bulk",
            "data_type": "bulk_ATAC",
            "condition_bams": self._bulk_condition_bams(
                exp_ctx, files, condition_col=condition_col),
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
        }
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
            result.setdefault("validation_level", "beta")
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
                motif_result.setdefault("validation_level", "beta")
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
                "gtf": exp_ctx.get("gtf") or exp_ctx.get("gtf_path"),
                "promoter_upstream": exp_ctx.get("promoter_upstream"),
                "promoter_downstream": exp_ctx.get("promoter_downstream"),
                "output_dir": str(Path(files[0]).parent / "bulk_atac_annotation"),
            },
        )
        if isinstance(ann_result, dict):
            ann_result.setdefault("data_type", "bulk_ATAC")
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
