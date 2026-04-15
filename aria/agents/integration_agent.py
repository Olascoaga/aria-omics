"""
ARIA IntegrationAgent
---------------------
Synthesizes findings across omics modalities into a coherent biological picture.

Three integration strategies (in order of priority):

  1. WNN (Weighted Nearest Neighbors)
     - For: 2 modalities from the SAME cells (multiome: RNA + ATAC)
     - What it does: builds a joint embedding weighting each modality
       by its information content per cell
     - Output: shared cell embedding, modality weights per cell,
       joint clusters that reflect both transcriptome and chromatin
     - ParameterAdvisor: advises on k (neighbors)
     - DebateCouncil: invoked when modality weights are imbalanced
       (e.g. RNA=0.91, ATAC=0.09 suggests ATAC QC failure)

  2. Peak-to-gene linking
     - For: any experiment with both RNA and ATAC data
     - What it does: correlates peak accessibility with nearby gene
       expression across cells, identifying regulatory relationships
     - Output: peak-gene links stored as Tunnels in ARIAMemory
       (this is the answer to "what regulates this gene?")
     - DebateCouncil: invoked for high-confidence links that lack
       supporting CTCF or footprinting evidence

  3. MOFA+ (Multi-Omics Factor Analysis)
     - For: 3+ modalities OR when user wants latent factor decomposition
     - What it does: decomposes variance into shared and modality-specific
       factors, identifies which biological programs drive each factor
     - Output: factor scores per cell, factor loadings per feature,
       variance explained per modality
     - ParameterAdvisor: advises on number of factors
     - DebateCouncil: always invoked for factor interpretation claims

Integration order:
  Checkpoint 2 reveals available modalities
  OrchestratorAgent selects strategy based on modalities present
  IntegrationAgent executes and stores cross-modal links as Tunnels

Cross-modal conflicts:
  When RNA and ATAC disagree (open chromatin but no expression,
  or expressed gene but closed chromatin), the DebateCouncil is
  called to reason about the biology before reporting.
  These conflicts are often the most biologically interesting findings.
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

log = logging.getLogger("aria.integration")


INTEGRATION_SYSTEM = """
You are ARIA's IntegrationAgent — a specialist in multi-omics data integration.

Your expertise:
- WNN (Weighted Nearest Neighbors): joint embedding of RNA + ATAC
- Peak-to-gene linking: correlating chromatin accessibility with expression
- MOFA+: latent factor decomposition across 3+ modalities
- Cross-modal conflict resolution: reasoning about discordant signals

Critical knowledge:

WNN:
  - Modality weights reflect information content per cell, not quality alone
  - RNA weight > 0.85 often means ATAC data is low quality OR the biological
    question is primarily transcriptional (both are valid — distinguish them)
  - Cells near cluster boundaries often switch modality dominance
  - WNN clusters can differ from RNA-only clusters — this is expected and
    biologically meaningful, not a bug

PEAK-TO-GENE LINKING:
  - Correlation-based linking has high false positive rate
  - Distance cutoff (typically 500kb) must be justified, not arbitrary
  - Positive correlation (open + expressed) is the easy case
  - Negative correlation (open + NOT expressed) is the interesting case:
    could be poised enhancer, silencer, or technical artifact
  - Always require at least 2 independent lines of evidence before
    claiming a regulatory relationship (correlation + CTCF/footprinting
    or correlation + HiC loop)

MOFA+:
  - Number of factors is the critical hyperparameter
  - Too few: biological programs are confounded
  - Too many: factors become noise-driven
  - Factor 1 often captures cell cycle or technical variation — check this
  - Variance explained per modality tells you which omics drives the biology
  - If one modality explains <10% variance across all factors, question its quality

CROSS-MODAL CONFLICTS:
  Gene expressed but chromatin closed:
    → Possible: distal regulatory element, trans-regulation, technical artifact
    → Always flag — do not silently discard
  Chromatin open but gene not expressed:
    → Possible: poised enhancer, TF binding without activation, cell state transition
    → Biologically important — report with confidence MEDIUM and explanation
  HiC loop present but no expression link:
    → Possible: structural loop (CTCF-anchored) without regulatory function
    → Requires additional evidence before claiming enhancer-promoter interaction
""".strip()


class IntegrationAgent(BaseAgent):

    name        = "integration_agent"
    description = (
        "Multi-omics integration: WNN (RNA+ATAC), "
        "peak-to-gene linking, MOFA+ for 3+ modalities."
    )

    # Minimum cells required per modality for reliable WNN
    MIN_CELLS_WNN = 200

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm, api_key)
        self.advisor = ParameterAdvisor(memory, llm)

        from aria.utils.environment_manager import env_manager
        self.env = env_manager

    # ── Main entry point ──────────────────────────────────────────────────

    def run(self, experiment_id: str, context: dict) -> dict:
        """
        Select and execute the appropriate integration strategy.

        context must contain:
          - exp_context:        from DataAuditAgent
          - biological_intent:  from OrchestratorAgent
          - rna_results:        from RNAAgent (optional)
          - chromatin_results:  from ChromatinAgent (optional)
          - hic_results:        from GenomeArchAgent (optional)
        """
        exp_ctx  = context.get("exp_context", {})
        intent   = context.get("biological_intent", {})
        mods     = exp_ctx.get("modalities", {})

        self.publish_status(experiment_id,
                            "IntegrationAgent: selecting strategy...", 0.0)

        # Determine which modalities are available
        has_scrna  = "scRNA"    in mods
        has_scatac = "scATAC"   in mods
        has_rna    = has_scrna or "bulk_RNA" in mods
        has_atac   = has_scatac or "bulk_ATAC" in mods
        has_hic    = "HiC"      in mods
        has_chip   = any(m in mods for m in ("ChIP", "CUT_AND_RUN", "CUT_AND_TAG"))

        n_modalities = sum([has_rna, has_atac, has_hic, has_chip])

        if n_modalities < 2:
            self.publish_finding(
                experiment_id,
                {"error": "Integration requires at least 2 modalities. "
                          f"Only found: {list(mods.keys())}"},
                Confidence.INSUFFICIENT,
            )
            return {"status": "failed", "reason": "insufficient_modalities"}

        strategy  = self._select_strategy(
            has_scrna, has_scatac, has_hic, n_modalities, intent
        )
        results   = {}

        self.publish_status(
            experiment_id,
            f"Integration strategy: {strategy}", 0.1
        )

        # ── Execute selected strategy ─────────────────────────────────────

        if strategy in ("wnn", "wnn_then_mofa"):
            if has_scrna and has_scatac:
                self.publish_status(experiment_id,
                                    "Running WNN (RNA + ATAC)...", 0.2)
                wnn_result = self._run_wnn(
                    experiment_id, exp_ctx, intent, mods
                )
                results["wnn"] = wnn_result

        if has_rna and has_atac:
            self.publish_status(experiment_id,
                                "Linking peaks to genes...", 0.5)
            p2g_result = self._run_peak_to_gene(
                experiment_id, exp_ctx, intent, mods,
                context.get("rna_results", {}),
                context.get("chromatin_results", {}),
            )
            results["peak_to_gene"] = p2g_result

        if strategy in ("mofa", "wnn_then_mofa") and n_modalities >= 2:
            self.publish_status(experiment_id,
                                "Running MOFA+ factor decomposition...", 0.7)
            mofa_result = self._run_mofa(
                experiment_id, exp_ctx, intent, mods
            )
            results["mofa"] = mofa_result

        # ── Resolve cross-modal conflicts ─────────────────────────────────
        if results:
            self.publish_status(experiment_id,
                                "Resolving cross-modal conflicts...", 0.9)
            conflicts = self._resolve_conflicts(
                experiment_id, exp_ctx, intent, results
            )
            results["conflicts"] = conflicts

        self.publish_status(experiment_id,
                            "IntegrationAgent complete.", 1.0)
        return {"status": "done", "strategy": strategy, "findings": results}

    # ── Strategy selection ────────────────────────────────────────────────

    def _select_strategy(self, has_scrna: bool, has_scatac: bool,
                          has_hic: bool, n_modalities: int,
                          intent: dict) -> str:
        """
        Select the optimal integration strategy.

        Decision tree:
          same-cell RNA + ATAC (multiome) → WNN (primary) + peak-to-gene
          3+ modalities OR factor request → MOFA+
          RNA + ATAC (bulk or different cells) → peak-to-gene only
        """
        question = intent.get("user_question", "").lower()

        FACTOR_KEYWORDS = [
            "factor", "latent", "mofa", "decomposition",
            "shared variation", "modality-specific",
        ]

        user_wants_mofa = any(k in question for k in FACTOR_KEYWORDS)

        if has_scrna and has_scatac:
            if n_modalities >= 3 or user_wants_mofa:
                return "wnn_then_mofa"
            return "wnn"
        elif n_modalities >= 3 or user_wants_mofa:
            return "mofa"
        else:
            return "peak_to_gene_only"

    # ── WNN integration ───────────────────────────────────────────────────

    def _run_wnn(self, experiment_id: str, exp_ctx: dict,
                  intent: dict, mods: dict) -> dict:
        """
        WNN integration of scRNA-seq + scATAC-seq from the same cells.
        Uses ParameterAdvisor to select k (number of neighbors).
        Invokes DebateCouncil when modality weights are imbalanced.
        """
        # ParameterAdvisor for WNN k
        k_decision = self.advisor.advise_wnn_k(
            multimodal_adata=None,  # real adata loaded in script
            experiment_id=experiment_id,
            biological_context=intent,
        )

        # Escalate to user (Checkpoint 3)
        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=3,
            question=self._format_wnn_checkpoint(k_decision),
            options=[
                f"Use recommended k={k_decision.chosen_value}",
                "Enter custom k value",
                "Skip WNN integration",
            ],
            context={"wnn_k_decision": k_decision.decision_id},
        )

        # Execute WNN in isolated environment
        rna_files  = mods.get("scRNA", [])
        atac_files = mods.get("scATAC", [])

        wnn_result = self.env.run_in_stack(
            stack="integration",
            script_path="aria/scripts/integration_wnn.py",
            params={
                "rna_files":   rna_files,
                "atac_files":  atac_files,
                "genome":      exp_ctx.get("genome", "hg38"),
                "organism":    exp_ctx.get("organism", "Homo sapiens"),
                "k":           k_decision.chosen_value,
                "output_dir":  self._get_output_dir(rna_files),
            },
        )

        if wnn_result.get("status") == "success":
            self._evaluate_wnn_weights(
                experiment_id, wnn_result, intent
            )

        return wnn_result

    def _evaluate_wnn_weights(self, experiment_id: str,
                               wnn_result: dict, intent: dict):
        """
        Evaluate WNN modality weights and invoke DebateCouncil
        when imbalanced weights suggest data quality issues.
        """
        rna_weight  = wnn_result.get("mean_rna_weight", 0.6)
        atac_weight = wnn_result.get("mean_atac_weight", 0.4)

        summary = (
            f"WNN integration: RNA weight={rna_weight:.2f}, "
            f"ATAC weight={atac_weight:.2f}"
        )

        if rna_weight > 0.85:
            # Imbalanced — debate required
            try:
                from aria.agents.debate_council import DebateCouncil
                council = DebateCouncil(llm=self.llm, max_rounds=2)

                claim = (
                    f"WNN modality weights are severely imbalanced: "
                    f"RNA={rna_weight:.2f}, ATAC={atac_weight:.2f}. "
                    f"The joint embedding is dominated by RNA signal."
                )
                debate = council.resolve(
                    topic="wnn_modality_weight_imbalance",
                    initial_claim=claim,
                    evidence={
                        "mean_rna_weight":  rna_weight,
                        "mean_atac_weight": atac_weight,
                        "n_cells":          wnn_result.get("n_cells", 0),
                        "atac_frip":        wnn_result.get("atac_frip", None),
                    },
                    biological_context=intent,
                )
                wnn_result["weight_debate"] = debate.consensus
                wnn_result["weight_verdict"] = debate.verdict.value

                conf = Confidence.LOW if rna_weight > 0.90 else Confidence.MEDIUM
            except Exception as e:
                log.warning(f"DebateCouncil unavailable: {e}")
                conf = Confidence.MEDIUM
        else:
            conf = Confidence.HIGH

        self.publish_finding(
            experiment_id,
            {"summary": summary,
             "rna_weight": rna_weight,
             "atac_weight": atac_weight},
            conf,
        )

    # ── Peak-to-gene linking ──────────────────────────────────────────────

    def _run_peak_to_gene(self, experiment_id: str, exp_ctx: dict,
                           intent: dict, mods: dict,
                           rna_results: dict,
                           chromatin_results: dict) -> dict:
        """
        Link ATAC peaks to nearby genes using correlation across cells.
        Stores high-confidence links as Tunnels in ARIAMemory.

        This answers: "What regulatory elements control this gene?"
        """
        rna_files  = mods.get("scRNA", mods.get("bulk_RNA", []))
        atac_files = mods.get("scATAC", mods.get("bulk_ATAC", []))
        peaks_path = (chromatin_results.get("scATAC", {})
                      .get("findings", {})
                      .get("peaks", {})
                      .get("peaks_path", ""))

        p2g_result = self.env.run_in_stack(
            stack="integration",
            script_path="aria/scripts/integration_peak2gene.py",
            params={
                "rna_files":    rna_files,
                "atac_files":   atac_files,
                "peaks_path":   peaks_path,
                "genome":       exp_ctx.get("genome", "hg38"),
                "organism":     exp_ctx.get("organism", "Homo sapiens"),
                "distance_kb":  500,   # 500kb window (standard)
                "min_corr":     0.3,   # minimum correlation threshold
                "output_dir":   self._get_output_dir(rna_files),
            },
        )

        if p2g_result.get("status") == "success":
            self._store_peak_gene_tunnels(
                experiment_id, p2g_result
            )
            self._debate_peak_gene_links(
                experiment_id, p2g_result, intent
            )

        return p2g_result

    def _store_peak_gene_tunnels(self, experiment_id: str,
                                  p2g_result: dict):
        """
        Store high-confidence peak-to-gene links as Tunnels in ARIAMemory.
        Tunnels represent cross-modal connections (ATAC <-> RNA).
        These persist across sessions and build the lab's regulatory map.
        """
        top_links = p2g_result.get("top_links", [])[:50]

        for link in top_links:
            if link.get("correlation", 0) >= 0.4:
                try:
                    self.memory.create_tunnel(
                        tunnel_id=str(uuid.uuid4())[:8],
                        wing_id=experiment_id,
                        from_hall=f"{experiment_id}_scATAC",
                        to_hall=f"{experiment_id}_scRNA",
                        entity=link.get("gene", "unknown"),
                        description=(
                            f"Peak {link.get('peak', '?')} -> "
                            f"Gene {link.get('gene', '?')}: "
                            f"r={link.get('correlation', 0):.3f}, "
                            f"dist={link.get('distance_kb', 0):.0f}kb"
                        ),
                        confidence=(
                            "high"   if link.get("correlation", 0) >= 0.6 else
                            "medium"
                        ),
                    )
                except Exception as e:
                    log.warning(f"Could not store tunnel: {e}")

        n_stored = min(len(top_links), 50)
        log.info(f"Stored {n_stored} peak-gene links as Tunnels in ARIAMemory")

    def _debate_peak_gene_links(self, experiment_id: str,
                                 p2g_result: dict, intent: dict):
        """
        Invoke DebateCouncil for high-confidence peak-gene links.
        Requires 2 independent lines of evidence before claiming
        a regulatory relationship.
        """
        n_links   = p2g_result.get("n_links", 0)
        top_links = p2g_result.get("top_links", [])[:5]

        if n_links == 0:
            self.publish_finding(
                experiment_id,
                {"summary": "No significant peak-gene links found",
                 "n_links": 0},
                Confidence.INSUFFICIENT,
            )
            return

        try:
            from aria.agents.debate_council import DebateCouncil
            council = DebateCouncil(llm=self.llm, max_rounds=2)

            top_genes = [l.get("gene", "?") for l in top_links]
            claim = (
                f"{n_links} peak-gene regulatory links identified. "
                f"Top linked genes: {', '.join(top_genes)}"
            )

            debate = council.resolve(
                topic="peak_to_gene_regulatory_links",
                initial_claim=claim,
                evidence={
                    "n_links":           n_links,
                    "top_links":         top_links,
                    "distance_cutoff_kb": 500,
                    "min_correlation":    0.3,
                    "ctcf_validated":     p2g_result.get("ctcf_validated", False),
                    "hic_corroborated":   p2g_result.get("hic_corroborated", False),
                },
                biological_context=intent,
            )

            conf_map = {
                "high":         Confidence.HIGH,
                "medium":       Confidence.MEDIUM,
                "low":          Confidence.LOW,
                "insufficient": Confidence.INSUFFICIENT,
            }
            self.publish_finding(
                experiment_id,
                {"summary":     debate.consensus,
                 "n_links":     n_links,
                 "top_genes":   top_genes,
                 "verdict":     debate.verdict.value,
                 "limitations": debate.limitations},
                conf_map.get(debate.confidence, Confidence.MEDIUM),
            )

        except Exception as e:
            log.warning(f"DebateCouncil unavailable for p2g: {e}")
            conf = (Confidence.HIGH   if n_links > 1000 else
                    Confidence.MEDIUM if n_links > 100  else
                    Confidence.LOW)
            self.publish_finding(
                experiment_id,
                {"summary": f"{n_links} peak-gene links identified",
                 "top_genes": [l.get("gene") for l in top_links]},
                conf,
            )

    # ── MOFA+ integration ─────────────────────────────────────────────────

    def _run_mofa(self, experiment_id: str, exp_ctx: dict,
                   intent: dict, mods: dict) -> dict:
        """
        MOFA+ factor decomposition for 2+ modalities.
        Identifies shared and modality-specific sources of variation.
        Always uses DebateCouncil for factor interpretation.
        """
        n_factors_decision = self._advise_mofa_factors(
            experiment_id, intent, mods
        )

        mofa_result = self.env.run_in_stack(
            stack="integration",
            script_path="aria/scripts/integration_mofa.py",
            params={
                "modalities":    {
                    k: v for k, v in mods.items()
                    if k in ("scRNA", "bulk_RNA", "scATAC", "bulk_ATAC", "HiC")
                },
                "genome":        exp_ctx.get("genome", "hg38"),
                "organism":      exp_ctx.get("organism", "Homo sapiens"),
                "n_factors":     n_factors_decision,
                "output_dir":    "/tmp/aria_mofa",
            },
        )

        if mofa_result.get("status") == "success":
            self._interpret_mofa_factors(
                experiment_id, mofa_result, intent
            )

        return mofa_result

    def _advise_mofa_factors(self, experiment_id: str,
                              intent: dict, mods: dict) -> int:
        """
        Advise on number of MOFA+ factors using biological context.
        Stores decision in ARIAMemory for reproducibility.
        """
        question    = intent.get("user_question", "").lower()
        n_modalities = len([m for m in mods if m != "unknown"])

        COMPLEX_KEYWORDS = ["atlas", "trajectory", "heterogeneous",
                             "many cell types", "complex"]
        complexity = "complex" if any(k in question for k in COMPLEX_KEYWORDS) \
                     else "moderate"

        # Layer 1: intent-based range
        if complexity == "complex":
            n_factors = 20
        elif n_modalities >= 3:
            n_factors = 15
        else:
            n_factors = 10

        # Store decision
        try:
            self.memory.store_decision(
                decision_id=str(uuid.uuid4())[:8],
                wing_id=experiment_id,
                checkpoint=3,
                question="MOFA+ number of factors",
                decision=str(n_factors),
                rationale=(
                    f"n_factors={n_factors} based on {n_modalities} modalities "
                    f"and {complexity} experimental complexity."
                ),
                made_by="advisor",
            )
        except Exception:
            pass

        return n_factors

    def _interpret_mofa_factors(self, experiment_id: str,
                                 mofa_result: dict, intent: dict):
        """
        Interpret MOFA+ factors using DebateCouncil.
        Always check if Factor 1 captures technical variation.
        """
        try:
            from aria.agents.debate_council import DebateCouncil
            council = DebateCouncil(llm=self.llm, max_rounds=2)

            factors        = mofa_result.get("top_factors", [])
            variance_exp   = mofa_result.get("variance_explained", {})
            factor1_top    = mofa_result.get("factor1_top_features", [])

            claim = (
                f"MOFA+ identified {len(factors)} factors. "
                f"Top factor explains: {list(variance_exp.values())[:3]}. "
                f"Factor 1 top features: {factor1_top[:5]}"
            )

            debate = council.resolve(
                topic="mofa_factor_interpretation",
                initial_claim=claim,
                evidence={
                    "n_factors":        len(factors),
                    "variance_explained": variance_exp,
                    "factor1_features": factor1_top,
                    "cell_cycle_check": mofa_result.get("cell_cycle_factor", False),
                },
                biological_context=intent,
            )

            conf_map = {
                "high": Confidence.HIGH, "medium": Confidence.MEDIUM,
                "low": Confidence.LOW,  "insufficient": Confidence.INSUFFICIENT,
            }
            self.publish_finding(
                experiment_id,
                {"summary":    debate.consensus,
                 "n_factors":  len(factors),
                 "verdict":    debate.verdict.value,
                 "limitations": debate.limitations},
                conf_map.get(debate.confidence, Confidence.MEDIUM),
            )

        except Exception as e:
            log.warning(f"MOFA factor debate failed: {e}")
            self.publish_finding(
                experiment_id,
                {"summary": f"MOFA+ factors: {mofa_result.get('n_factors', 0)}"},
                Confidence.MEDIUM,
            )

    # ── Cross-modal conflict resolution ───────────────────────────────────

    def _resolve_conflicts(self, experiment_id: str, exp_ctx: dict,
                            intent: dict, results: dict) -> dict:
        """
        Identify and reason about cross-modal discordances.

        Cases:
          - Open chromatin but no expression → poised enhancer?
          - Expressed gene but closed chromatin → distal regulation?
          - WNN cluster != RNA cluster → chromatin drives cell identity?
          - HiC loop without expression correlation → structural loop?

        Each conflict is reported with biological interpretation
        and stored as a finding with explicit confidence and limitations.
        """
        conflicts = []

        # Check WNN vs RNA cluster discordance
        wnn    = results.get("wnn", {})
        p2g    = results.get("peak_to_gene", {})
        mofa   = results.get("mofa", {})

        if wnn.get("status") == "success":
            n_discordant = wnn.get("n_discordant_clusters", 0)
            if n_discordant > 0:
                conflict = {
                    "type":        "wnn_rna_cluster_discordance",
                    "description": (
                        f"{n_discordant} cell clusters differ between "
                        f"RNA-only and WNN joint embedding. "
                        f"Chromatin accessibility provides additional "
                        f"cell identity information beyond transcriptomics."
                    ),
                    "confidence":  "medium",
                    "action":      "Review discordant clusters manually",
                }
                conflicts.append(conflict)

        # Check peak-gene negative correlations (open but silent)
        if p2g.get("status") == "success":
            n_negative = p2g.get("n_negative_correlations", 0)
            if n_negative > 0:
                conflict = {
                    "type":        "open_chromatin_silent_gene",
                    "description": (
                        f"{n_negative} peaks are accessible but their "
                        f"linked genes are not expressed. "
                        f"These may represent poised enhancers, silencers, "
                        f"or developmental regulatory elements."
                    ),
                    "n_events":    n_negative,
                    "confidence":  "medium",
                    "action":      "Prioritize for footprinting and HiC validation",
                }
                conflicts.append(conflict)
                # Store as findings for NarrativeAgent
                self.publish_finding(
                    experiment_id,
                    {"summary": conflict["description"],
                     "n_events": n_negative,
                     "type": "poised_regulatory_elements"},
                    Confidence.MEDIUM,
                )

        if conflicts:
            log.info(f"Identified {len(conflicts)} cross-modal conflicts")

        return {
            "n_conflicts": len(conflicts),
            "conflicts":   conflicts,
        }

    # ── Checkpoint formatters ─────────────────────────────────────────────

    def _format_wnn_checkpoint(self, k_decision) -> str:
        lines = ["WNN Integration — k Parameter\n"]
        lines.append("Candidates evaluated:\n")
        for c in k_decision.candidates:
            marker = " [RECOMMENDED]" if c.recommended else ""
            rna_w  = c.metrics.get("rna_weight", 0)
            atac_w = c.metrics.get("atac_weight", 0)
            lines.append(
                f"  k={c.value:2d}  RNA_weight={rna_w:.2f}  "
                f"ATAC_weight={atac_w:.2f}{marker}"
            )
        lines.append(f"\n  {k_decision.justification}")
        lines.append("\nApprove k to proceed with WNN integration:")
        return "\n".join(lines)

    # ── Helpers ───────────────────────────────────────────────────────────

    def _get_output_dir(self, files: list) -> str:
        if files:
            return str(Path(files[0]).parent / "aria_integration")
        return "/tmp/aria_integration"

    def receive(self, message):
        pass
