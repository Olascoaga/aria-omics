"""
ARIA GenomeArchAgent
--------------------
3D genome organization analysis: TADs, loops, compartments A/B.

Supported formats:
  .cool / .mcool  — Cooler format (multi-resolution)
  .hic            — Juicer/4D Nucleome format
  .pairs / .pairs.gz — Pre-alignment pairs (upstream of .cool)

Analysis hierarchy (from coarse to fine resolution):
  1. Quality control & ICE/KR balancing  (all resolutions)
  2. Compartments A/B                    (100kb - 1Mb resolution)
  3. TAD calling — Insulation Score      (10kb - 40kb resolution)
  4. Loop calling — dots/chromosight     (5kb - 10kb resolution)

Memory strategy (critical for HiC):
  - NEVER load the full genome-wide matrix into RAM
  - All analysis uses out-of-core operations (cooler chromosome slices)
  - Resolution selection at Checkpoint 2 determines RAM requirements
  - GenomeArchAgent estimates RAM before dispatching to EnvironmentManager

RAM estimates by resolution (diploid human genome):
  1Mb  → ~50MB   (compartments only)
  100kb → ~1GB   (compartments + coarse TADs)
  40kb  → ~8GB   (TADs)
  10kb  → ~100GB (fine TADs + loops — requires HPC)
  5kb   → ~400GB (loop calling — requires HPC cluster)

ParameterAdvisor role:
  - Advises on Insulation Score window_size
  - Runs calibration on chr1 only (fast proxy for full genome)
  - Presents boundary_strength vs window_size tradeoff to user
  - Stores approved window_size in ARIAMemory for reproducibility

DebateCouncil role:
  - Compartment A/B interpretation (B ≠ absolute silencing)
  - TAD boundary conservation claims across conditions
  - Enhancer-promoter loop predictions (distance and score thresholds)
"""

from __future__ import annotations

import json
import logging
import uuid
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence, MessageType
from aria.llm.provider import LLMProvider, TaskTier
from aria.llm.parameter_advisor import ParameterAdvisor
from aria.memory.memory import ARIAMemory

log = logging.getLogger("aria.hic")


GENOME_ARCH_SYSTEM = """
You are ARIA's GenomeArchAgent — a specialist in 3D genome organization.

Your expertise:
- Hi-C / Micro-C data processing: ICE/KR balancing, quality assessment
- Compartments A/B: eigenvector decomposition, biological interpretation
- TAD calling: Insulation Score (primary), boundary strength assessment
- Chromatin loops: enhancer-promoter interactions, CTCF-anchored loops
- Resolution selection: matching analysis to biological question

Critical knowledge:

COMPARTMENTS:
  - Compartment A = active chromatin (high gene density, H3K27ac, accessible)
  - Compartment B = repressed chromatin (heterochromatin, H3K27me3, AT-rich)
  - BUT: Compartment B is NOT always transcriptionally silent. Some B-compartment
    genes are expressed (especially at A/B boundaries). Never claim "B = silent".
  - Compartment identity can flip between cell types (A-to-B or B-to-A switches)
    — these are among the most biologically interesting findings.
  - PC1 sign is arbitrary: always validate with gene density or H3K27ac if available.

TADs:
  - Insulation Score is the standard (Crane 2015). It is robust to sparse data.
  - window_size controls sensitivity: large window = mega-domains,
    small window = sub-TADs. This is the key hyperparameter.
  - TAD boundaries are enriched for CTCF binding and convergent CTCF motifs.
  - ~30% of TAD boundaries are NOT conserved across cell types.
  - Do NOT claim TADs are "cell type specific" without cross-comparison.

LOOPS:
  - Enhancer-promoter loops typically span 10kb-2Mb.
  - Not all Hi-C loops are functional. Validation with CTCF ChIP or ATAC
    footprinting is required before claiming regulatory significance.
  - Resolution < 5kb is required for reliable loop calling.
    At 10kb, many loops are missed or merged.

RESOLUTION and MEMORY:
  - Always confirm available RAM before selecting analysis resolution.
  - Chromosome-by-chromosome analysis is mandatory for resolutions < 40kb.
  - .mcool files contain multiple resolutions: always inspect before loading.
""".strip()


# ── RAM estimation table ──────────────────────────────────────────────────────
# Approximate RAM requirements for human genome (hg38) full-genome analysis
# Scale down ~10x for mouse, ~50x for drosophila
RAM_ESTIMATES_GB = {
    1_000_000: 0.05,    # 1Mb
    500_000:   0.2,     # 500kb
    100_000:   1.0,     # 100kb
    40_000:    8.0,     # 40kb
    10_000:    100.0,   # 10kb
    5_000:     400.0,   # 5kb
}

GENOME_SCALE = {
    "Homo sapiens":            1.0,
    "Mus musculus":            0.9,
    "Drosophila melanogaster": 0.02,
    "C. elegans":              0.003,
    "Danio rerio":             0.5,
    "S. cerevisiae":           0.001,
}


class GenomeArchAgent(BaseAgent):

    name        = "genome_arch_agent"
    description = (
        "3D genome organization: Hi-C/Micro-C QC, balancing, "
        "compartments A/B, TAD calling, loop calling."
    )

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm, api_key)
        self.advisor = ParameterAdvisor(memory, llm)

        from aria.utils.environment_manager import env_manager
        self.env = env_manager

    # ── Main entry point ──────────────────────────────────────────────────

    def run(self, experiment_id: str, context: dict) -> dict:
        exp_ctx = context.get("exp_context", {})
        intent  = context.get("biological_intent", {})
        mods    = exp_ctx.get("modalities", {})
        hic_files = mods.get("HiC", [])

        if not hic_files:
            self.publish_finding(
                experiment_id,
                {"error": "No Hi-C/Micro-C files detected"},
                Confidence.INSUFFICIENT,
            )
            return {"status": "failed", "reason": "no_hic_data"}

        self.publish_status(experiment_id,
                            "GenomeArchAgent starting...", 0.0)

        # Phase 1: Inspect files and available resolutions
        file_info = self._inspect_hic_files(hic_files, exp_ctx)

        # Checkpoint 2: Resolution selection
        # This determines everything downstream (RAM, analysis depth)
        resolution_decision = self._advise_resolution(
            experiment_id, intent, file_info, exp_ctx
        )

        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=2,
            question=self._format_resolution_checkpoint(
                file_info, resolution_decision
            ),
            options=[
                f"Use recommended resolution ({resolution_decision['recommended_resolution']:,}bp)",
                "Select different resolution",
                "Skip 3D genome analysis",
            ],
            context={
                "file_info":           file_info,
                "resolution_decision": resolution_decision,
            },
        )

        # Phase 2: QC and balancing
        self.publish_status(experiment_id, "Hi-C QC and balancing...", 0.2)
        qc_result = self.env.run_in_stack(
            stack="hic",
            script_path="aria/scripts/hic_qc_and_balance.py",
            params={
                "files":      hic_files,
                "genome":     exp_ctx.get("genome", "hg38"),
                "resolution": resolution_decision["recommended_resolution"],
            },
        )

        if qc_result.get("status") == "error":
            log.warning(f"Hi-C QC failed: {qc_result.get('details','')}")

        # Checkpoint 3: QC review
        if qc_result.get("warnings"):
            self.publish_escalation(
                experiment_id=experiment_id,
                checkpoint=3,
                question=self._format_qc_checkpoint(qc_result),
                options=[
                    "Continue — warnings noted in report",
                    "Abort — data quality insufficient",
                ],
                context={"qc": qc_result},
            )

        # Phase 3: Topology
        self.publish_status(experiment_id, "Computing 3D topology...", 0.5)
        topology_result = self._run_topology(
            experiment_id, exp_ctx, intent,
            hic_files, resolution_decision, qc_result
        )

        findings = {
            "file_info":   file_info,
            "qc":          qc_result,
            "topology":    topology_result,
        }

        self.publish_status(experiment_id,
                            "GenomeArchAgent complete.", 1.0)
        return {"status": "done", "findings": findings}

    # ── Phase 1: File inspection ──────────────────────────────────────────

    def _inspect_hic_files(self, files: list, exp_ctx: dict) -> dict:
        """
        Inspect Hi-C files: format, available resolutions, estimated size.
        Uses out-of-core inspection (reads metadata only, not data).
        """
        result = self.env.run_in_stack(
            stack="hic",
            script_path="aria/scripts/hic_inspect.py",
            params={
                "files":  files,
                "genome": exp_ctx.get("genome", "hg38"),
            },
        )

        if result.get("status") == "error":
            # Fallback: infer from filenames
            log.warning("File inspection failed — using filename inference")
            result = self._infer_file_info(files)

        return result

    def _infer_file_info(self, files: list) -> dict:
        """Infer file info from extensions when cooler is unavailable."""
        formats = []
        for f in files:
            if f.endswith(".cool") or f.endswith(".mcool"):
                formats.append("cooler")
            elif f.endswith(".hic"):
                formats.append("hic")
            elif f.endswith(".pairs") or f.endswith(".pairs.gz"):
                formats.append("pairs")

        return {
            "status":               "partial",
            "files":                files,
            "formats":              list(set(formats)),
            "available_resolutions": [1_000_000, 500_000, 100_000, 40_000,
                                       10_000, 5_000],
            "note": "Full inspection requires aria-hic-env",
        }

    # ── Resolution ParameterAdvisor ───────────────────────────────────────

    def _advise_resolution(self, experiment_id: str, intent: dict,
                            file_info: dict, exp_ctx: dict) -> dict:
        """
        Advise on analysis resolution based on:
          1. Biological question (compartments → coarse, loops → fine)
          2. Available resolutions in the file
          3. Estimated RAM requirements
          4. Lab's historical decisions (ARIAMemory)

        This is the most critical parameter decision in Hi-C analysis.
        The ParameterAdvisor stores the approved resolution so future
        experiments with similar questions learn from this choice.
        """
        question  = intent.get("user_question", "").lower()
        organism  = exp_ctx.get("organism", "Homo sapiens")
        scale     = GENOME_SCALE.get(organism, 1.0)
        available = file_info.get("available_resolutions",
                                   [1_000_000, 100_000, 40_000, 10_000])

        # ── Layer 1: Intent-based resolution range ────────────────────────
        LOOP_KEYWORDS = ["loop", "enhancer", "promoter", "ctcf anchor",
                          "hi-c loop", "chromatin loop"]
        TAD_KEYWORDS  = ["tad", "domain", "topological", "boundary",
                          "insulation", "compartment switch"]
        COMP_KEYWORDS = ["compartment", "a/b", "heterochromatin",
                          "euchromatin", "a-b", "lamina"]

        needs_loops       = any(k in question for k in LOOP_KEYWORDS)
        needs_tads        = any(k in question for k in TAD_KEYWORDS)
        needs_compartments = any(k in question for k in COMP_KEYWORDS)

        if needs_loops:
            target_res = 5_000
            analysis   = "loop_calling"
        elif needs_tads:
            target_res = 40_000
            analysis   = "tad_calling"
        elif needs_compartments:
            target_res = 100_000
            analysis   = "compartments"
        else:
            # Default: TAD level (most common use case)
            target_res = 40_000
            analysis   = "tad_calling"

        # ── Layer 2: RAM feasibility check ───────────────────────────────
        base_ram = RAM_ESTIMATES_GB.get(target_res, 100.0)
        req_ram  = base_ram * scale

        # If required resolution not available, use nearest coarser
        if target_res not in available:
            coarser = [r for r in sorted(available, reverse=True)
                       if r >= target_res]
            target_res = coarser[0] if coarser else available[0]
            req_ram    = RAM_ESTIMATES_GB.get(target_res, 1.0) * scale

        # ── Layer 3: Memory-based candidates for ParameterAdvisor ────────
        # Offer 3 resolution options: recommended + one finer + one coarser
        sorted_avail = sorted(available)
        idx          = sorted_avail.index(target_res) \
                       if target_res in sorted_avail else 0

        candidates = []
        for i in [max(0, idx - 1), idx, min(len(sorted_avail) - 1, idx + 1)]:
            res      = sorted_avail[i]
            ram_need = RAM_ESTIMATES_GB.get(res, 100.0) * scale
            candidates.append({
                "resolution":    res,
                "ram_required_gb": round(ram_need, 1),
                "analysis_depth": self._resolution_to_depth(res),
                "recommended":   (res == target_res),
            })

        # ── Store decision in memory ──────────────────────────────────────
        decision = {
            "recommended_resolution": target_res,
            "analysis_type":          analysis,
            "ram_required_gb":        round(req_ram, 1),
            "candidates":             candidates,
            "rationale":              (
                f"Resolution {target_res:,}bp selected for {analysis}. "
                f"Requires ~{req_ram:.1f}GB RAM for {organism}. "
                f"Intent: '{question[:80]}'"
            ),
        }

        # Recall lab's historical decisions
        hist = self.advisor._recall_similar_decisions(
            experiment_id, "hic_resolution", intent
        )
        if hist:
            decision["historical_precedent"] = hist[:2]

        # Store for reproducibility
        try:
            self.memory.store_decision(
                decision_id=str(uuid.uuid4())[:8],
                wing_id=experiment_id,
                checkpoint=2,
                question=f"Hi-C resolution for {analysis}",
                decision=str(target_res),
                rationale=decision["rationale"],
                made_by="advisor",
            )
        except Exception as e:
            log.warning(f"Could not store resolution decision: {e}")

        return decision

    # ── Phase 3: Topology ─────────────────────────────────────────────────

    def _run_topology(self, experiment_id: str, exp_ctx: dict,
                       intent: dict, files: list,
                       resolution_decision: dict,
                       qc_result: dict) -> dict:
        """
        Run topology analysis: compartments + TADs + optional loops.
        Uses out-of-core chromosome-by-chromosome processing.
        """
        resolution  = resolution_decision["recommended_resolution"]
        analysis    = resolution_decision["analysis_type"]
        organism    = exp_ctx.get("organism", "Homo sapiens")

        # ── Compartments A/B (always run at coarsest available res) ──────
        comp_resolution = max(r for r in
                              resolution_decision.get("candidates", [{}])
                              if r.get("resolution", 0) >= 100_000
                              for r in [r["resolution"]]) \
                          if any(c.get("resolution", 0) >= 100_000
                                 for c in resolution_decision.get("candidates", [])
                                 ) else 100_000

        compartments = self.env.run_in_stack(
            stack="hic",
            script_path="aria/scripts/hic_topology.py",
            params={
                "files":           files,
                "genome":          exp_ctx.get("genome", "hg38"),
                "organism":        organism,
                "analysis":        "compartments",
                "resolution":      comp_resolution,
                "chromosomes":     "all",
                "out_of_core":     True,
            },
        )

        if compartments.get("status") == "success":
            self._publish_compartments_finding(
                experiment_id, compartments, intent
            )

        # ── TAD calling via Insulation Score ─────────────────────────────
        tads = {}
        if resolution <= 100_000:
            # ParameterAdvisor calibrates window_size on chr1 first
            window_decision = self._advise_insulation_window(
                experiment_id, intent, files,
                resolution, exp_ctx.get("genome", "hg38")
            )

            tads = self.env.run_in_stack(
                stack="hic",
                script_path="aria/scripts/hic_topology.py",
                params={
                    "files":       files,
                    "genome":      exp_ctx.get("genome", "hg38"),
                    "organism":    organism,
                    "analysis":    "tads",
                    "resolution":  resolution,
                    "window_size": window_decision["recommended_window"],
                    "chromosomes": "all",
                    "out_of_core": True,
                },
            )

            if tads.get("status") == "success":
                self._publish_tads_finding(experiment_id, tads)

        # ── Loop calling (only if high-res data available) ────────────────
        loops = {}
        if resolution <= 10_000 and analysis == "loop_calling":
            loops = self.env.run_in_stack(
                stack="hic",
                script_path="aria/scripts/hic_topology.py",
                params={
                    "files":       files,
                    "genome":      exp_ctx.get("genome", "hg38"),
                    "organism":    organism,
                    "analysis":    "loops",
                    "resolution":  resolution,
                    "chromosomes": "all",
                    "out_of_core": True,
                },
            )

            if loops.get("status") == "success":
                self._publish_loops_finding(
                    experiment_id, loops, intent
                )

        return {
            "compartments": compartments,
            "tads":         tads,
            "loops":        loops,
        }

    # ── Insulation Score window_size ParameterAdvisor ────────────────────

    def _advise_insulation_window(self, experiment_id: str,
                                   intent: dict, files: list,
                                   resolution: int,
                                   genome: str) -> dict:
        """
        Calibrate Insulation Score window_size using chr1 as proxy.

        Runs 3 window sizes on chr1 only (fast), computes mean
        boundary strength for each, presents to user at Checkpoint 3.

        This is the ParameterAdvisor pattern applied to Hi-C:
          Layer 1: Intent-constrained range
          Layer 2: Objective metric (boundary_strength)
          Layer 3: Lab memory recall
        """
        question = intent.get("user_question", "").lower()

        # Layer 1: Intent-based window range
        MEGADOMAIN_KEYWORDS = ["megadomain", "large domain", "broad",
                                "topological domain large"]
        SUBDOMAIN_KEYWORDS  = ["sub-tad", "subdomain", "fine", "small",
                                "local"]

        if any(k in question for k in MEGADOMAIN_KEYWORDS):
            windows = [300_000, 500_000, 1_000_000]
        elif any(k in question for k in SUBDOMAIN_KEYWORDS):
            windows = [50_000, 100_000, 200_000]
        else:
            # Standard: 3 windows at 3x, 5x, 10x the resolution
            windows = [resolution * 3, resolution * 5, resolution * 10]

        # Layer 2: Calibrate on chr1 (fast proxy)
        calibration = self.env.run_in_stack(
            stack="hic",
            script_path="aria/scripts/hic_topology.py",
            params={
                "files":          files,
                "genome":         genome,
                "analysis":       "insulation_calibration",
                "resolution":     resolution,
                "windows":        windows,
                "chromosomes":    ["chr1"],  # chr1 only for speed
                "out_of_core":    True,
            },
        )

        # Choose window with highest mean boundary strength
        best_window  = windows[1]  # default: middle option
        best_strength = 0.0

        if calibration.get("status") == "success":
            for w, strength in calibration.get("boundary_strengths", {}).items():
                if float(strength) > best_strength:
                    best_strength = float(strength)
                    best_window   = int(w)

        decision = {
            "recommended_window":  best_window,
            "calibration_results": calibration.get("boundary_strengths", {}),
            "windows_tested":      windows,
            "rationale": (
                f"window_size={best_window:,}bp selected based on highest "
                f"mean boundary strength ({best_strength:.3f}) across chr1. "
                f"Reflects {intent.get('analysis_type','TAD')} biological scale."
            ),
        }

        # Store in memory
        try:
            self.memory.store_decision(
                decision_id=str(uuid.uuid4())[:8],
                wing_id=experiment_id,
                checkpoint=3,
                question=f"Insulation Score window_size at {resolution:,}bp",
                decision=str(best_window),
                rationale=decision["rationale"],
                made_by="advisor",
            )
        except Exception as e:
            log.warning(f"Could not store window decision: {e}")

        # Escalate to user
        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=3,
            question=self._format_window_checkpoint(decision, resolution),
            options=[
                f"Use recommended window ({best_window:,}bp)",
                "Enter custom window size",
                "Skip TAD calling",
            ],
            context={"window_decision": decision},
        )

        return decision

    # ── Findings publishers ───────────────────────────────────────────────

    def _publish_compartments_finding(self, experiment_id: str,
                                       result: dict, intent: dict):
        """
        Publish compartment A/B finding.
        DebateCouncil review for any A-to-B or B-to-A switch claims.
        """
        n_switches = result.get("n_ab_switches", 0)
        summary    = (
            f"Compartment A/B: {result.get('pct_A', 0):.1f}% A-compartment, "
            f"{result.get('pct_B', 0):.1f}% B-compartment"
        )

        if n_switches > 0:
            summary += f", {n_switches} A/B switches detected"

            # DebateCouncil for switch claims (biologically significant)
            try:
                from aria.agents.debate_council import DebateCouncil
                council = DebateCouncil(llm=self.llm, max_rounds=2)
                claim = (
                    f"{n_switches} genomic regions show compartment switching "
                    f"(A-to-B or B-to-A) compared to the reference condition."
                )
                debate = council.resolve(
                    topic="compartment_switching",
                    initial_claim=claim,
                    evidence={
                        "n_switches":       n_switches,
                        "switch_regions":   result.get("switch_regions", [])[:5],
                        "pct_A":            result.get("pct_A", 0),
                        "pct_B":            result.get("pct_B", 0),
                        "pc1_validated":    result.get("pc1_validated", False),
                    },
                    biological_context=intent,
                )
                result["debate_verdict"]   = debate.verdict.value
                result["debate_consensus"] = debate.consensus
                result["limitations"]      = debate.limitations
            except Exception as e:
                log.warning(f"DebateCouncil unavailable: {e}")

        conf = (Confidence.HIGH   if result.get("pc1_validated") else
                Confidence.MEDIUM)
        self.publish_finding(
            experiment_id,
            {"summary": summary, "n_ab_switches": n_switches},
            conf,
        )

    def _publish_tads_finding(self, experiment_id: str, result: dict):
        n_tads = result.get("n_tads", 0)
        conf   = (Confidence.HIGH   if n_tads > 1000 else
                  Confidence.MEDIUM if n_tads > 100  else
                  Confidence.LOW)
        self.publish_finding(
            experiment_id,
            {"summary": f"TADs: {n_tads:,} domains called",
             "n_tads":  n_tads,
             "median_tad_size_kb": result.get("median_size_kb", 0)},
            conf,
        )

    def _publish_loops_finding(self, experiment_id: str,
                                result: dict, intent: dict):
        n_loops = result.get("n_loops", 0)

        if n_loops > 0:
            # DebateCouncil for loop regulatory claims
            try:
                from aria.agents.debate_council import DebateCouncil
                council = DebateCouncil(llm=self.llm, max_rounds=2)
                top_loops = result.get("top_loops", [])[:5]
                claim = (
                    f"{n_loops:,} chromatin loops identified. "
                    f"Top loops span genes/enhancers: "
                    f"{', '.join(str(l) for l in top_loops)}"
                )
                debate = council.resolve(
                    topic="chromatin_loop_regulatory_significance",
                    initial_claim=claim,
                    evidence={
                        "n_loops":          n_loops,
                        "top_loops":        top_loops,
                        "resolution":       result.get("resolution", 0),
                        "ctcf_validated":   result.get("ctcf_validated", False),
                        "atac_corroborated": result.get("atac_corroborated", False),
                    },
                    biological_context=intent,
                )
                result["debate_verdict"]   = debate.verdict.value
                result["debate_consensus"] = debate.consensus
            except Exception as e:
                log.warning(f"DebateCouncil unavailable: {e}")

        conf = (Confidence.HIGH   if n_loops > 5000 else
                Confidence.MEDIUM if n_loops > 500  else
                Confidence.LOW)
        self.publish_finding(
            experiment_id,
            {"summary": f"Loops: {n_loops:,} chromatin loops called",
             "n_loops":  n_loops},
            conf,
        )

    # ── Checkpoint formatters ─────────────────────────────────────────────

    def _format_resolution_checkpoint(self, file_info: dict,
                                       decision: dict) -> str:
        lines = ["Hi-C Resolution Selection\n"]
        lines.append(f"Files: {len(file_info.get('files', []))} Hi-C file(s)")
        lines.append(
            f"Available resolutions: "
            f"{[f'{r:,}bp' for r in file_info.get('available_resolutions', [])]}"
        )
        lines.append("")
        lines.append("Resolution options:")
        for c in decision.get("candidates", []):
            marker = " [RECOMMENDED]" if c.get("recommended") else ""
            lines.append(
                f"  {c['resolution']:>10,}bp  "
                f"RAM ~{c['ram_required_gb']:.1f}GB  "
                f"{c['analysis_depth']}{marker}"
            )
        lines.append("")
        lines.append(f"Rationale: {decision.get('rationale', '')}")
        lines.append("\nSelect resolution to proceed:")
        return "\n".join(lines)

    def _format_qc_checkpoint(self, qc_result: dict) -> str:
        lines = ["Hi-C QC Warnings:\n"]
        for w in qc_result.get("warnings", []):
            lines.append(f"  * {w}")
        lines.append(
            f"\n  Cis/trans ratio: {qc_result.get('cis_trans_ratio', 'N/A')}"
        )
        lines.append(
            f"  Valid pairs: {qc_result.get('n_valid_pairs', 'N/A'):,}"
            if isinstance(qc_result.get("n_valid_pairs"), int) else
            f"  Valid pairs: {qc_result.get('n_valid_pairs', 'N/A')}"
        )
        lines.append("\nHow would you like to proceed?")
        return "\n".join(lines)

    def _format_window_checkpoint(self, decision: dict,
                                   resolution: int) -> str:
        lines = [
            f"Insulation Score Window Size (resolution: {resolution:,}bp)\n"
        ]
        lines.append("Calibration results (chr1 boundary strength):")
        for w, strength in decision.get("calibration_results", {}).items():
            marker = " [RECOMMENDED]" if int(w) == decision[
                "recommended_window"] else ""
            lines.append(
                f"  window={int(w):>10,}bp  "
                f"boundary_strength={float(strength):.4f}{marker}"
            )
        lines.append(f"\n  {decision.get('rationale', '')}")
        lines.append("\nApprove window size to proceed with full genome TAD calling:")
        return "\n".join(lines)

    # ── Helpers ───────────────────────────────────────────────────────────

    def _resolution_to_depth(self, resolution: int) -> str:
        DEPTH = {
            1_000_000: "compartments only",
            500_000:   "compartments",
            100_000:   "compartments + coarse TADs",
            40_000:    "TADs",
            10_000:    "TADs + loops",
            5_000:     "fine loops (HPC recommended)",
        }
        return DEPTH.get(resolution, "unknown")

    def receive(self, message):
        pass
