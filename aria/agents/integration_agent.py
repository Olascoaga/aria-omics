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
  - Positive correlation and negative correlation are both descriptive
    observations until independent evidence supports interpretation
  - Always require at least 2 independent lines of evidence before
    claiming a regulatory relationship

MOFA+:
  - Number of factors is the critical hyperparameter
  - Too few: biological programs are confounded
  - Too many: factors become noise-driven
  - Factor interpretation requires explicit evidence; do not infer biology
    from feature names alone
  - Variance explained per modality tells you which omics drives the biology
  - If one modality explains <10% variance across all factors, question its quality

CROSS-MODAL CONFLICTS:
  Gene expressed but chromatin closed:
    → Discordant evidence; flag it and avoid causal interpretation
  Chromatin open but gene not expressed:
    → Discordant evidence; flag it and avoid causal interpretation
  HiC loop present but no expression link:
    → Structural contact without expression support
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
        mods     = dict(exp_ctx.get("modalities", {}) or {})

        self.publish_status(experiment_id,
                            "IntegrationAgent: selecting strategy...", 0.0)

        from aria.utils.multiome_contracts import (
            enforce_same_cell_integration_gate,
            infer_multiome_contract,
        )

        multiome_contract = infer_multiome_contract(exp_ctx)
        is_paired_mudata = (
            multiome_contract.get("object_type") == "paired_mudata"
            and multiome_contract.get("same_cell") is True
        )

        # Determine which modalities are available (a paired .h5mu carries both
        # scRNA and scATAC in a single shared-cell container).
        has_scrna  = ("scRNA" in mods) or is_paired_mudata
        has_scatac = ("scATAC" in mods) or is_paired_mudata
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

        # V48 I1 hard gate: any same-cell integration strategy MUST satisfy
        # the paired multiome contract before touching downstream algorithms.
        if strategy in ("wnn", "wnn_then_mofa"):
            gate = enforce_same_cell_integration_gate(exp_ctx)
            results["multiome_gate"] = gate
            if not gate["allowed"]:
                blocker_msgs = "; ".join(
                    str(b.get("message") or b.get("check"))
                    for b in gate.get("blockers", [])
                )
                self.publish_finding(
                    experiment_id,
                    {
                        "error": (
                            "Same-cell multi-omics integration blocked by "
                            f"multiome contract gate: {blocker_msgs}"
                        ),
                        "blockers": gate.get("blockers", []),
                    },
                    Confidence.INSUFFICIENT,
                )
                return {
                    "status": "blocked",
                    "reason": "multiome_same_cell_gate_blocked",
                    "strategy": strategy,
                    "multiome_gate": gate,
                }

            if gate["contract"].get("object_type") == "paired_mudata":
                self.publish_status(
                    experiment_id,
                    "Running paired .h5mu contract & joint RNA+ATAC QC...",
                    0.05,
                )
                paired_qc_res = self._run_paired_qc(
                    experiment_id, exp_ctx, intent, gate
                )
                results["paired_qc"] = paired_qc_res
                if paired_qc_res.get("status") != "success":
                    self.publish_finding(
                        experiment_id,
                        {
                            "error": (
                                "Paired .h5mu joint QC failed: "
                                f"{paired_qc_res.get('details', paired_qc_res.get('error_type', 'unknown'))}"
                            ),
                            "paired_qc": paired_qc_res,
                        },
                        Confidence.INSUFFICIENT,
                    )
                    return {
                        "status": "failed",
                        "reason": "paired_qc_failed",
                        "strategy": strategy,
                        "findings": results,
                    }
                validated_h5mu = paired_qc_res.get("output_path")
                if validated_h5mu:
                    mods = {
                        **mods,
                        "scRNA": [validated_h5mu],
                        "scATAC": [validated_h5mu],
                    }

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

    # ── Paired .h5mu contract & joint QC (V48 I1) ─────────────────────────

    def _run_paired_qc(
        self,
        experiment_id: str,
        exp_ctx: dict,
        intent: dict,
        gate: dict,
    ) -> dict:
        """Validate paired `.h5mu` and run joint RNA+ATAC QC over shared cells."""
        evidence_paths = list(gate.get("evidence_paths") or [])
        if not evidence_paths:
            return {
                "status": "error",
                "error_type": "MissingPairedH5mu",
                "details": "No paired .h5mu path resolved from multiome contract.",
            }

        h5mu_path = evidence_paths[0]
        qc_overrides = dict(exp_ctx.get("paired_qc_params") or {})
        params = {
            "h5mu_path": h5mu_path,
            "output_dir": qc_overrides.pop(
                "output_dir",
                self._get_output_dir(evidence_paths),
            ),
            "organism": exp_ctx.get("organism", "Homo sapiens"),
            "genome": exp_ctx.get("genome"),
        }
        for opt_key in (
            "sample_id",
            "reference_h5ad",
            "reference_label_col",
            "target_label_col",
            "reference_sample_col",
        ):
            if opt_key in exp_ctx and exp_ctx[opt_key] is not None:
                params[opt_key] = exp_ctx[opt_key]
        params.update(qc_overrides)

        qc_res = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/integration_paired_qc.py",
            params=params,
        )
        if qc_res.get("status") == "success":
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        f"Paired multiome joint QC retained "
                        f"{qc_res.get('n_cells_after')}/{qc_res.get('n_cells_before')} "
                        f"shared cells ({qc_res.get('n_genes')} genes, "
                        f"{qc_res.get('n_peaks')} peaks)."
                    ),
                    "n_cells_before": qc_res.get("n_cells_before"),
                    "n_cells_after": qc_res.get("n_cells_after"),
                    "label_join": qc_res.get("label_join"),
                },
                Confidence.HIGH,
            )
        return qc_res

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

        # Escalate to user (Checkpoint 3) and block until the dispatch-loop UI
        # resolves it. Otherwise the advisor value would run before user input.
        _, user_decision = self.publish_blocking_escalation(
            experiment_id=experiment_id,
            checkpoint=3,
            question=self._format_wnn_checkpoint(k_decision),
            options=[
                f"Use recommended k={k_decision.chosen_value}",
                "Other: Enter custom k value",
                "Skip WNN integration",
            ],
            context={
                "wnn_k_decision": k_decision.decision_id,
                "analysis_type": "wnn_integration",
                "parameter_name": "k",
            },
        )
        choice = self.checkpoint_choice_text(user_decision)
        if not choice:
            return {
                "status": "error",
                "error_type": "CheckpointUnresolved",
                "details": "WNN k checkpoint was not resolved.",
            }
        if "skip" in choice.lower():
            k_decision.chosen_by = "user"
            k_decision.approved_by_user = True
            return {
                "status": "skipped",
                "reason": "user_skipped_wnn",
                "details": "User skipped WNN integration at checkpoint 3.",
            }
        override = self.numeric_checkpoint_override(user_decision, int)
        if hasattr(self.advisor, "approve_decision"):
            k_decision = self.advisor.approve_decision(
                k_decision,
                user_override=override,
            )
        else:
            if override is not None:
                k_decision.chosen_value = override
                k_decision.chosen_by = "user"
            k_decision.approved_by_user = True

        # Execute WNN in isolated environment (aria-chromatin-env provides the
        # locked muon + mudata + scanpy + leidenalg + umap-learn stack).
        rna_files  = list(mods.get("scRNA", []) or [])
        atac_files = list(mods.get("scATAC", []) or [])
        wnn_overrides = dict(exp_ctx.get("wnn_params") or {})

        h5mu_path = wnn_overrides.pop("h5mu_path", None)
        if not h5mu_path:
            for cand in [*rna_files, *atac_files]:
                if str(cand).lower().endswith(".h5mu"):
                    h5mu_path = str(cand)
                    break

        ref_paths = [h5mu_path] if h5mu_path else (rna_files or atac_files)
        wnn_params = {
            "rna_files":   rna_files,
            "atac_files":  atac_files,
            "genome":      exp_ctx.get("genome", "hg38"),
            "organism":    exp_ctx.get("organism", "Homo sapiens"),
            "k":           k_decision.chosen_value,
            "seed":        int(exp_ctx.get("seed", 0)),
            "output_dir":  wnn_overrides.pop(
                "output_dir",
                self._get_output_dir(ref_paths),
            ),
        }
        if h5mu_path:
            wnn_params["h5mu_path"] = h5mu_path
        wnn_params.update(wnn_overrides)

        wnn_result = self.env.run_in_stack(
            stack="chromatin",
            script_path="aria/scripts/integration_wnn.py",
            params=wnn_params,
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
        rna_weight  = wnn_result.get("mean_rna_weight")
        atac_weight = wnn_result.get("mean_atac_weight")

        if rna_weight is None or atac_weight is None:
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        "WNN integration completed without measured modality "
                        "weights (reported as null; not estimated)."
                    ),
                    "rna_weight": rna_weight,
                    "atac_weight": atac_weight,
                },
                Confidence.LOW,
            )
            return

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
        if n_factors_decision is None:
            return {
                "status": "skipped",
                "reason": "user_skipped_mofa",
                "details": "User skipped MOFA+ at checkpoint 3.",
            }

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
                              intent: dict, mods: dict) -> int | None:
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

        question = (
            "MOFA+ factor count\n\n"
            f"Recommended n_factors={n_factors} based on {n_modalities} "
            f"modalities and {complexity} experimental complexity.\n\n"
            "Approve the factor count to proceed with MOFA+."
        )
        _, user_decision = self.publish_blocking_escalation(
            experiment_id=experiment_id,
            checkpoint=3,
            question=question,
            options=[
                f"Use recommended n_factors={n_factors}",
                "Other: Enter custom n_factors",
                "Skip MOFA+",
            ],
            context={
                "analysis_type": "mofa_integration",
                "parameter_name": "n_factors",
                "recommended_value": n_factors,
            },
        )
        choice = self.checkpoint_choice_text(user_decision)
        if not choice:
            return None
        if "skip" in choice.lower():
            return None
        override = self.numeric_checkpoint_override(user_decision, int)
        chosen_factors = override if override is not None else n_factors

        self.memory.store_decision(
            decision_id=str(uuid.uuid4())[:8],
            wing_id=experiment_id,
            checkpoint=3,
            question="MOFA+ number of factors",
            decision=str(chosen_factors),
            rationale=(
                f"n_factors={chosen_factors} based on user approval of the "
                f"checkpoint recommendation for {n_modalities} modalities and "
                f"{complexity} experimental complexity."
            ),
            made_by="user" if override is not None else "advisor",
        )

        return chosen_factors

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
                    "technical_factor_flag": mofa_result.get(
                        "technical_factor_flag", False
                    ),
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

        Cases are descriptive discordances only; ARIA stores them with explicit
        confidence and limitations, without assigning a mechanism from pattern
        alone.
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

        # Check peak-gene negative correlations (discordant evidence)
        if p2g.get("status") == "success":
            n_negative = p2g.get("n_negative_correlations", 0)
            if n_negative > 0:
                conflict = {
                    "type":        "accessibility_expression_discordance",
                    "description": (
                        f"{n_negative} peak-gene pair(s) show negative "
                        f"accessibility-expression correlation. "
                        f"Treat these as discordant evidence until an "
                        f"independent validation assay supports mechanism."
                    ),
                    "n_events":    n_negative,
                    "confidence":  "medium",
                    "action":      "Prioritize for independent validation",
                }
                conflicts.append(conflict)
                # Store as findings for NarrativeAgent
                self.publish_finding(
                    experiment_id,
                    {"summary": conflict["description"],
                     "n_events": n_negative,
                     "type": "accessibility_expression_discordance"},
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
            if c.metrics.get("measured") is False:
                basis = c.metrics.get("basis", "prior")
                lines.append(
                    f"  k={c.value:2d}  basis={basis}  "
                    f"pre_run_metrics=not_computed{marker}"
                )
            else:
                rna_w = c.metrics.get("rna_weight")
                atac_w = c.metrics.get("atac_weight")
                if rna_w is not None and atac_w is not None:
                    lines.append(
                        f"  k={c.value:2d}  RNA_weight={rna_w:.2f}  "
                        f"ATAC_weight={atac_w:.2f}{marker}"
                    )
                else:
                    lines.append(f"  k={c.value:2d}{marker}")
        lines.append(f"\n  {k_decision.justification}")
        lines.append(
            "\n  Note: RNA/ATAC modality weights are reported only after the "
            "WNN script executes; this checkpoint does not display estimated "
            "weights as measured metrics."
        )
        lines.append("\nApprove k to proceed with WNN integration:")
        return "\n".join(lines)

    # ── Helpers ───────────────────────────────────────────────────────────

    def _get_output_dir(self, files: list) -> str:
        if files:
            return str(Path(files[0]).parent / "aria_integration")
        return "/tmp/aria_integration"

    def receive(self, message):
        pass
