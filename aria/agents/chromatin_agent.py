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
    validation_level = "scaffold"
    dispatch_enabled = False
    REQUIRED_SCRIPTS = (
        "aria/scripts/chromatin_qc.py",
        "aria/scripts/chromatin_peaks.py",
    )
    PLANNED_SCRIPTS = (
        "aria/scripts/chromatin_motifs.py",
        "aria/scripts/chromatin_differential.py",
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
        findings = {}

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

        # 2. LSI dimensionality reduction
        # ParameterAdvisor decides how many SVD components to use
        # and explicitly recommends discarding component 1
        self.publish_status(experiment_id, "scATAC LSI reduction...", 0.3)
        lsi_params = self._advise_lsi_params(
            experiment_id, intent, qc_result
        )
        findings["lsi_params"] = lsi_params

        # 3. Peak calling
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

        # 4. TF motif enrichment (if requested)
        if self._needs_tf_analysis(intent):
            self.publish_status(experiment_id,
                                "TF motif enrichment...", 0.7)
            motif_result = self._run_motif_enrichment(
                experiment_id, exp_ctx, intent, peaks_result, "scATAC"
            )
            findings["motifs"] = motif_result

        return {"status": "done", "findings": findings}

    # ── Bulk ATAC-seq pipeline ────────────────────────────────────────────

    def _run_bulk_atac(self, experiment_id: str, exp_ctx: dict,
                       intent: dict, files: list) -> dict:
        """Bulk ATAC-seq: QC, peak calling, differential accessibility."""
        findings = {}

        # QC
        qc_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/chromatin_qc.py",
            params={
                "data_type": "bulk_ATAC",
                "files":     files,
                "genome":    exp_ctx.get("genome", "hg38"),
            },
        )
        findings["qc"] = qc_result
        if qc_result.get("status") == "success":
            self._publish_qc_finding(experiment_id, qc_result, "bulk_ATAC")

        # Peak calling
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
        findings["peaks"] = peaks_result

        # Differential accessibility if comparison defined
        if intent.get("comparison") and qc_result.get("status") == "success":
            diff_result = self._run_differential_accessibility(
                experiment_id, exp_ctx, intent, peaks_result
            )
            findings["differential_accessibility"] = diff_result

        return {"status": "done", "findings": findings}

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

        # Flag critical QC failures
        frip      = qc_result.get("frip", 1.0)
        tss_score = qc_result.get("tss_enrichment", 10.0)

        if frip < 0.2 or tss_score < 4:
            conf = Confidence.LOW

        self.publish_finding(
            experiment_id,
            {"summary":    f"{modality} QC: "
                           f"FRiP={frip:.3f}, TSS={tss_score:.2f}",
             "modality":   modality,
             "frip":       frip,
             "tss_score":  tss_score,
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
