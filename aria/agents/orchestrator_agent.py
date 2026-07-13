"""
ARIA OrchestratorAgent
----------------------
The central coordinator. Receives the user's biological question,
designs the analysis plan, and dispatches specialized agents.

Flow:
  1. Parse biological question → structured analysis plan
  2. Delegate to DataAuditAgent (always first)
  3. Wait for Checkpoint 1 confirmation
  4. INTERACTIVE EXPERIMENTAL DESIGN via DesignAgent (v4.0)
  5. Present analysis plan → Checkpoint 2
  6. Intercept modifications via Checkpoint 3 (Parameter Tuning)
  7. DISPATCH modality agents in dependency order
  8. Collect findings from all agents
  9. Trigger IntegrationAgent if multimodal
  10. Dispatch NarrativeAgent for final report
  11. Present Checkpoint 5 (final review)
"""

from __future__ import annotations

import logging
import os
import threading
import uuid
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Message, MessageType, Confidence
from aria.llm.provider import LLMProvider
from aria.llm.prompt_boundary import (
    PromptDataField,
    build_untrusted_prompt,
    system_with_untrusted_boundary,
)
from aria.memory.memory import ARIAMemory
from aria.runtime.experiment_session import ExperimentSession
from aria.utils.privacy import use_egress_policy
from aria.utils.provenance import collect_provenance

log = logging.getLogger("aria.orchestrator")


ORCHESTRATOR_SYSTEM = """
You are ARIA's Orchestrator — the strategic brain of a bioinformatics agent system.

Your job:
1. Understand the biological question behind a user's request
2. Design an analysis plan that uses available data modalities optimally
3. Identify dependencies between analyses (what must run before what)
4. Assess whether findings across modalities converge or conflict
5. Decide when integration analysis adds real biological value

You know these agents are available:
- data_audit_agent: always runs first
- scrna_agent: single-cell RNA-seq (QC, clustering, annotation, DE)
- bulk_rna_agent: bulk RNA-seq (DESeq2, pathway enrichment, plots)
- chromatin_agent: handles ATAC, ChIP, CUT&RUN, CUT&TAG
- genome_arch_agent: handles HiC, compartments, TADs, loops
- integration_agent: multimodal integration (needs 2+ modalities)
- narrative_agent: generates the final report

Always think about the biology first, then the methods.
""".strip()


# ── Agent registry ─────────────────────────────────────────────────────────────

AGENT_REGISTRY = {
    "setup_agent": {
        "module": "aria.agents.setup_agent",
        "class":  "SetupAgent",
    },
    "raw_ingestion_agent": {
        "module": "aria.agents.raw_ingestion_agent",
        "class":  "RawIngestionAgent",
    },
    "scrna_agent": {
        "module": "aria.agents.scrna_agent",
        "class":  "scRNAAgent",
    },
    "bulk_rna_agent": {
        "module": "aria.agents.bulk_rna_agent",
        "class":  "BulkRNAAgent",
    },
    "chromatin_agent": {
        "module": "aria.agents.chromatin_agent",
        "class":  "ChromatinAgent",
    },
    "genome_arch_agent": {
        "module": "aria.agents.genome_arch_agent",
        "class":  "GenomeArchAgent",
    },
    "integration_agent": {
        "module": "aria.agents.integration_agent",
        "class":  "IntegrationAgent",
    },
    "narrative_agent": {
        "module": "aria.agents.narrative_agent",
        "class":  "NarrativeAgent",
    },
    "rna_agent": {
        "module": "aria.agents.scrna_agent",
        "class":  "scRNAAgent",
    },
    "audit_agent": {
        "module": "aria.agents.audit_agent",
        "class":  "AuditAgent",
    },
}

MODALITY_TO_AGENT = {
    "scRNA":        "scrna_agent",
    "bulk_RNA":     "bulk_rna_agent",
    "bulk_RNA_raw": "bulk_rna_agent",
    "scATAC":      "chromatin_agent",
    "bulk_ATAC":   "chromatin_agent",
    "ChIP":        "chromatin_agent",
    "CUT_AND_RUN": "chromatin_agent",
    "CUT_AND_TAG": "chromatin_agent",
    "HiC":         "genome_arch_agent",
}

MODALITY_VALIDATION = {
    "scRNA": {
        "level": "production",
        "dispatch_enabled": True,
    },
    "bulk_RNA": {
        "level": "production",
        "dispatch_enabled": True,
    },
    "bulk_RNA_raw": {
        "level": "beta",
        "dispatch_enabled": True,
    },
    "scATAC": {
        "level": "beta",
        "dispatch_enabled": True,
        "reason": (
            "scATAC beta dispatch requires explicit CP3.5 acknowledgement. "
            "De-alpha closed 2026-06-15 (ADR-048): P0-P3 done; clustering + motif "
            "are externally concordant (SnapATAC2 on HC11), and the pseudobulk "
            "condition-DA lane is validated against edgeR/limma/DESeq2 on real "
            "biological replicates (5 young vs 5 old GSE278576 donors; ARIA a clean "
            "subset of R-DESeq2, LFC Spearman 0.61-0.76). The per-cluster wilcoxon "
            "marker path is exploratory/single-sample-fragile and stays caveated. "
            "Beta + requires_ack: dispatchable after ack, not fully autonomous."
        ),
    },
    "bulk_ATAC": {
        "level": "beta",
        "dispatch_enabled": True,
        "reason": (
            "Bulk ATAC V47 is open as a beta, acknowledgement-gated lane for "
            "measured QC, MACS3 peak calling, and replicate-aware DESeq2 "
            "differential accessibility over explicit peak-count matrices. "
            "DA requires explicit condition, replicate, and comparison metadata; "
            "under-specified designs skip honestly."
        ),
    },
    "ChIP": {
        "level": "scaffold",
        "dispatch_enabled": False,
        "reason": "ChIP validation is not closed.",
    },
    "CUT_AND_RUN": {
        "level": "scaffold",
        "dispatch_enabled": False,
        "reason": "CUT&RUN validation is not closed.",
    },
    "CUT_AND_TAG": {
        "level": "scaffold",
        "dispatch_enabled": False,
        "reason": "CUT&TAG validation is not closed.",
    },
    "HiC": {
        "level": "scaffold",
        # P0-3: Hi-C is dispatch-OFF by default. Its scripts run, but
        # `hic_topology` emits TADs/loops (unvalidated conclusions under a
        # scaffold level), so a default run is not publication-grade. Opt in
        # explicitly with ARIA_ALLOW_EXPERIMENTAL_HIC=1; such runs are stamped
        # experimental. This does not add or change any Hi-C science.
        "dispatch_enabled": False,
        "reason": (
            "Hi-C E2E (TAD/loop/compartment) validation is not closed; set "
            "ARIA_ALLOW_EXPERIMENTAL_HIC=1 to run it as experimental, "
            "not-publication-grade output."
        ),
        "experimental_env_flag": "ARIA_ALLOW_EXPERIMENTAL_HIC",
    },
}

MODALITY_VALIDATION_LEVELS = {
    modality: meta["level"]
    for modality, meta in MODALITY_VALIDATION.items()
}

INTEGRATION_VALIDATION = {
    "level": "scaffold",
    "dispatch_enabled": False,
    "reason": (
        "IntegrationAgent is scaffolded until WNN/MOFA+/peak-to-gene "
        "validation closes."
    ),
}


class OrchestratorAgent(BaseAgent):

    name        = "orchestrator"
    description = "Central coordinator — routes tasks to specialized agents."

    MODALITY_TO_AGENT = MODALITY_TO_AGENT
    MODALITY_VALIDATION = MODALITY_VALIDATION
    INTEGRATION_VALIDATION = INTEGRATION_VALIDATION

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider = None,
                 api_key: str = None):
        super().__init__(memory, llm=llm, api_key=api_key)
        self._sessions             = {}
        self._pending_checkpoints  = {}
        self._experiment_plans     = {}
        self._agent_results        = {}
        # v4.0: active DesignAgent instance (state machine)
        self._active_design_agent  = None
        # v4.1: pending dispatch held while audit CP 3.5 is open
        self._pending_dispatch     = {}
        # Per-experiment file log handlers (attached in run, detached at end)
        self._log_handlers         = {}

    def _get_session(self, experiment_id: str) -> ExperimentSession:
        sessions = getattr(self, "_sessions", None)
        if sessions is None:
            sessions = {}
            self._sessions = sessions
        if experiment_id not in sessions:
            sessions[experiment_id] = ExperimentSession(experiment_id=experiment_id)
        session = sessions[experiment_id]
        session.message_bus.register(self.name, self)
        if session.usage_log is None:
            session.usage_log = (
                Path.home() / ".aria" / "workspace" / experiment_id
                / "llm_usage.jsonl"
            )
        if session.llm_provider is None:
            provider = getattr(self, "llm", None)
            if provider is not None and hasattr(provider, "for_execution"):
                provider = provider.for_execution(
                    experiment_id=experiment_id,
                    usage_log=session.usage_log,
                    egress_policy=session.egress_policy,
                )
            session.llm_provider = provider
        return session

    def _register_session_agent(self, experiment_id: str, agent):
        """Attach a run-owned agent to the matching isolated broker."""
        self._get_session(experiment_id).message_bus.register(agent.name, agent)
        return agent

    def _sync_plan_record(self, session: ExperimentSession) -> None:
        plans = getattr(self, "_experiment_plans", None)
        if plans is None:
            plans = {}
            self._experiment_plans = plans
        plans[session.experiment_id] = session.as_plan_record()

    def _handle_design_checkpoint(
        self,
        experiment_id: str,
        checkpoint: int | float,
        user_decision: str,
    ) -> dict:
        session = self._get_session(experiment_id)
        design_agent = session.design_agent
        if design_agent is None:
            return {"status": "error", "message": "No active design session"}

        with use_egress_policy(session.egress_policy):
            result = design_agent.handle_user_response(
                experiment_id=experiment_id,
                checkpoint_num=checkpoint,
                choice=user_decision,
            )
        if result.get("status") == "cancelled":
            session.design_agent = None
            self._sync_plan_record(session)
            return {"status": "cancelled", "reason": "Design cancelled by user"}

        if result.get("status") != "done":
            self._sync_plan_record(session)
            return {"status": "design_in_progress", "next_step": result.get("step")}

        # Diseño completado → enriquecer contexto y pasar a CP2
        design = result["design"]
        exp_context = session.exp_context or {}
        return self._complete_design_and_publish_plan(
            experiment_id, exp_context, design)

    def _complete_design_and_publish_plan(
        self, experiment_id: str, exp_context: dict, design: dict,
    ) -> dict:
        """Enrich the context with the confirmed (or skipped) design, run design
        intelligence, build the analysis plan, and publish CP2.

        Shared by the normal DE design-completion path and the chromatin-only
        ``no_de_design`` skip path (DesignAgent returns ``status="skipped"``),
        so both reach plan/dispatch through identical wiring.
        """
        session = self._get_session(experiment_id)
        exp_context["organism"] = design["organism"]
        exp_context["genome"]   = design["genome"]
        exp_context["design"]   = design
        design_intelligence = self._run_design_intelligence(
            exp_context,
            session.intent,
        )
        exp_context["design_intelligence"] = design_intelligence
        self.memory.store_decision(
            decision_id=f"{experiment_id}_design_intelligence",
            wing_id=experiment_id,
            checkpoint="2.pre",
            question="Design Intelligence assessment",
            decision=design_intelligence.get("summary", "completed"),
            rationale=(
                "ARIA evaluated modality feasibility, recommended "
                "analyses, optional analyses, unsupported analyses, "
                "covariates, and design limitations before compute."
            ),
            made_by="design_intelligence",
        )
        session.exp_context = exp_context
        session.design_agent = None
        # Backward-compat mirror only; session.design_agent is authoritative.
        self._active_design_agent = None

        plan = self._design_analysis_plan(experiment_id, exp_context)
        plan["design_intelligence"] = design_intelligence
        # Annotate plan with resolved thresholds — single source of truth
        # so the preview and the actual run always agree.
        # F1 (ADR-055): the plan's |log2FC| cutoff is the fixed default, NOT
        # inferred from the question text. A TF/KO/OE study can opt into the
        # Exploratory/TF profile at CP3; the suggestion is surfaced but never
        # auto-applied.
        from aria.agents.bulk_rna_agent import (
            _default_lfc_threshold, suggest_lfc_profile)
        plan["lfc_threshold"]  = _default_lfc_threshold()
        plan["lfc_threshold_provenance"] = "default"
        plan["lfc_threshold_suggestion"] = suggest_lfc_profile(session.intent)
        plan["padj_threshold"] = 0.05
        session.plan = plan
        self._sync_plan_record(session)
        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=2,
            question=self._format_plan_summary(plan),
            options=[
                "Run recommended plan only",
                "Run recommended + optional supported analyses",
                "Modify plan",
                "Cancel",
            ],
            context={"plan": plan, "exp_context": exp_context},
        )
        return {"status": "plan_ready", "plan": plan}

    def run(self, experiment_id: str, context: dict) -> dict:
        session = self._get_session(experiment_id)
        # Open a per-experiment log file so background-thread agent failures
        # are recoverable without scrolling the TUI.
        try:
            from aria.utils.logging_setup import attach_experiment_log
            session.log_handler = attach_experiment_log(experiment_id)
            self._log_handlers[experiment_id] = session.log_handler
        except Exception as e:
            log.debug(f"Could not attach experiment log: {e}")

        # R6: append-only bus durability under a per-experiment path so a
        # mid-run crash keeps the findings/decisions of completed stages, and
        # concurrent runs never share a log file.
        try:
            session.message_bus.enable_persistence(
                str(Path.home() / ".aria" / "workspace" / experiment_id
                    / "bus_log.jsonl")
            )
        except Exception as e:
            log.debug(f"Could not enable bus persistence: {e}")

        with use_egress_policy(session.egress_policy):
            self.publish_status(experiment_id, "ARIA starting analysis...", 0.0)
            intent = self._parse_question(
                context["user_question"], llm=session.llm_provider
            )
            context["provenance"] = collect_provenance()
        session.context = context
        session.intent = intent
        session.status = "auditing"
        self._sync_plan_record(session)

        return {"status": "started", "intent": intent, "next": "data_audit"}

    def run_audit(self, experiment_id: str) -> dict:
        from aria.agents.data_audit_agent import DataAuditAgent
        session = self._get_session(experiment_id)
        plan = self._experiment_plans.get(experiment_id, {})
        agent = self._register_session_agent(
            experiment_id,
            DataAuditAgent(self.memory, llm=session.llm_provider),
        )
        with use_egress_policy(session.egress_policy):
            return agent.run(experiment_id, plan["context"])

    # ── Callback principal para checkpoints resueltos ───────────────────
    def on_checkpoint_resolved(self, message_id: str,
                                user_decision: str,
                                experiment_id: str,
                                corrections: dict | None = None) -> dict:
        session_bus = self._get_session(experiment_id).message_bus
        session_bus.resolve_checkpoint(
            message_id, {"choice": user_decision, "corrections": corrections})

        resolved_msg = next(
            (m for m in session_bus.get_log(experiment_id) if m.id == message_id),
            None,
        )
        if not resolved_msg:
            return {"status": "error", "message": "Checkpoint not found"}

        cp = resolved_msg.checkpoint
        context = resolved_msg.payload.get("context", {}) or {}

        if context.get("agent_parameter_checkpoint"):
            return {
                "status": "parameter_checkpoint_resolved",
                "checkpoint": cp,
                "sender": resolved_msg.sender,
            }

        # ── v4.0: Delegar checkpoints de diseño al DesignAgent ──────────
        if cp in (2.1, 2.2, 2.3, 2.4, 2.5, 2.6):
            session = self._get_session(experiment_id)
            if session.design_agent is not None:
                return self._handle_design_checkpoint(
                    experiment_id=experiment_id,
                    checkpoint=cp,
                    user_decision=user_decision,
                )

        # ── Checkpoints normales ────────────────────────────────────────
        if cp == 1:
            return self._after_checkpoint_1(
                experiment_id, user_decision, resolved_msg, corrections=corrections)
        elif cp == 2:
            return self._after_checkpoint_2(experiment_id, user_decision, resolved_msg)
        elif cp == 3:
            return self._after_checkpoint_3(experiment_id, user_decision, resolved_msg)
        elif cp == 3.5:
            return self._after_audit_checkpoint(experiment_id, user_decision, resolved_msg)

        return {"status": "ok"}

    # ── Checkpoint 1: Auditoría completada ──────────────────────────────
    def _after_checkpoint_1(self, experiment_id: str, decision: str, msg: Message,
                            corrections: dict | None = None) -> dict:
        if "cancel" in decision.lower():
            return {"status": "cancelled"}

        exp_context = msg.payload.get("context", {}).get("exp_context", {})

        # CHECKPOINT-1 "Correct metadata": apply the user's modality/organism/
        # genome corrections to exp_context BEFORE design/dispatch, so the right
        # pipeline runs. Previously this decision was a silent no-op and ARIA
        # proceeded with the (possibly wrong) auto-detected metadata.
        if "correct" in decision.lower() and corrections:
            from aria.agents.data_audit_agent import apply_metadata_corrections
            apply_metadata_corrections(exp_context, corrections)
            self.publish_finding(
                experiment_id,
                {"summary": (
                    "User corrected the data audit at checkpoint 1: "
                    f"modality(ies)={list((exp_context.get('modalities') or {}).keys())}, "
                    f"organism={exp_context.get('organism')}, "
                    f"genome={exp_context.get('genome')}."
                )},
                Confidence.HIGH,
            )
            self.memory.store_decision(
                decision_id=str(uuid.uuid4())[:8],
                wing_id=experiment_id,
                checkpoint=1,
                question="Correct data-audit metadata",
                decision=str(corrections),
                rationale="User overrode the auto-detected modality/organism/genome.",
                made_by="user",
            )

        # P1-8a / W-PRIV: honor the user's air-gapped choice from the sensitivity
        # checkpoint. Enabling it here (before design/dispatch) makes BOTH the
        # in-process LLM and every dispatched subprocess refuse network egress.
        from aria.utils.sensitivity import decision_enables_air_gapped
        session = self._get_session(experiment_id)
        if decision_enables_air_gapped(decision):
            from aria.utils.privacy import enable_air_gapped_runtime
            enable_air_gapped_runtime(
                reason="sensitivity_checkpoint",
                policy=session.egress_policy,
            )
            exp_context["air_gapped"] = True
            sensitivity = (msg.payload.get("context", {}) or {}).get("sensitivity", {})
            self.publish_finding(
                experiment_id,
                {"summary": (
                    "Air-gapped mode enabled at the sensitivity checkpoint: ALL "
                    "network egress (cloud LLM, Enrichr ORA, GEO/SRA) is blocked "
                    f"for this run (sensitivity={sensitivity.get('level', 'n/a')})."
                )},
                Confidence.HIGH,
            )
            self.memory.store_decision(
                decision_id=str(uuid.uuid4())[:8],
                wing_id=experiment_id,
                checkpoint=1,
                question="Data sensitivity / egress policy",
                decision=decision,
                rationale=(
                    "User enabled air-gapped mode for a sensitive input; ARIA "
                    "blocks all network egress for the run."
                ),
                made_by="user",
            )

        self.memory.store_decision(
            decision_id=str(uuid.uuid4())[:8],
            wing_id=experiment_id,
            checkpoint=1,
            question=msg.payload.get("question", ""),
            decision=decision,
            rationale="User confirmed data audit results",
            made_by="user",
        )

        # ── v4.0: Iniciar máquina de estados del DesignAgent ────────────
        from aria.agents.design_agent import DesignAgent
        design_agent = self._register_session_agent(
            experiment_id,
            DesignAgent(memory=self.memory, llm=session.llm_provider),
        )
        with use_egress_policy(session.egress_policy):
            result = design_agent.start_design(
                experiment_id=experiment_id,
                exp_context=exp_context,
                biological_intent=self._experiment_plans.get(
                    experiment_id, {}
                ).get("intent", {}),
            )

        if result.get("status") == "failed":
            return {"status": "cancelled", "reason": result.get("reason", "Design failed")}

        # Chromatin-only (non-DE) run: DesignAgent skipped the DE design phase
        # and synthesized a minimal design. Proceed straight to plan/CP2 so the
        # modality agent (e.g. ChromatinAgent for scATAC) reaches dispatch.
        if result.get("status") == "skipped":
            session.exp_context = exp_context
            self._experiment_plans[experiment_id]["exp_context"] = exp_context
            self._sync_plan_record(session)
            return self._complete_design_and_publish_plan(
                experiment_id, exp_context, result["design"])

        # Guardar el agente activo para seguir recibiendo sus checkpoints
        session.design_agent = design_agent
        session.exp_context = exp_context
        # Backward-compat mirror only; session.design_agent is authoritative.
        self._active_design_agent = design_agent
        self._experiment_plans[experiment_id]["exp_context"] = exp_context
        self._sync_plan_record(session)

        # El primer checkpoint de diseño (2.1) ya fue publicado.
        return {"status": "design_in_progress", "next_checkpoint": 2.1}

    # ── Checkpoint 2: Plan de análisis ─────────────────────────────────
    def _after_checkpoint_2(self, experiment_id: str, decision: str, msg: Message) -> dict:
        if "cancel" in decision.lower():
            return {"status": "cancelled"}

        plan        = msg.payload.get("context", {}).get("plan", {})
        exp_context = msg.payload.get("context", {}).get("exp_context", {})

        if "modify" in decision.lower():
            self.publish_escalation(
                experiment_id=experiment_id,
                checkpoint=3,
                question="Select the statistical thresholds for Differential Expression (padj & log2FC):\n\n"
                         "  1. Strict Profile (padj < 0.01, |log2FC| > 1.5)\n"
                         "  2. Standard Profile (padj < 0.05, |log2FC| > 1.0)\n"
                         "  3. Exploratory/TF Profile (padj < 0.05, |log2FC| > 0.58)\n"
                         "  4. High Sensitivity Profile (padj < 0.10, |log2FC| > 0.58)\n"
                         "  5. Keep Default (Agent decides)",
                options=["Strict", "Standard", "Exploratory/TF", "High Sensitivity", "Keep Default"],
                context={"plan": plan, "exp_context": exp_context}
            )
            return {"status": "modifying_plan"}

        if "optional" in decision.lower():
            exp_context["design_intelligence_selection"] = "recommended_plus_optional"
            exp_context["run_optional_supported"] = True
            rationale_suffix = " User requested optional supported analyses too."
        else:
            exp_context["design_intelligence_selection"] = "recommended_only"
            exp_context["run_optional_supported"] = False
            rationale_suffix = " Optional supported analyses were not added."

        self.memory.store_decision(
            decision_id=str(uuid.uuid4())[:8],
            wing_id=experiment_id,
            checkpoint=2,
            question="Analysis plan confirmed",
            decision=decision,
            rationale="User approved analysis plan directly." + rationale_suffix,
            made_by="user",
        )

        # Propagate resolved thresholds into exp_context so BulkRNAAgent
        # uses the same values shown in the plan preview.
        exp_context.setdefault("global_padj", plan.get("padj_threshold", 0.05))
        exp_context.setdefault("global_lfc",  plan.get("lfc_threshold",  1.0))
        exp_context.setdefault("lfc_threshold_provenance",
                               plan.get("lfc_threshold_provenance", "default"))

        # v4.0: inject plan contrasts into design for BulkRNAAgent
        plan_contrasts = plan.get("contrasts")
        if plan_contrasts and exp_context.get("design"):
            exp_context["design"]["plan_contrasts"] = plan_contrasts

        return self._gate_and_dispatch(experiment_id, plan, exp_context)

    # ── Checkpoint 3: Perfiles de umbral ───────────────────────────────
    def _after_checkpoint_3(self, experiment_id: str, decision: str, msg: Message) -> dict:
        if "cancel" in decision.lower():
            return {"status": "cancelled"}
            
        plan        = msg.payload.get("context", {}).get("plan", {})
        exp_context = msg.payload.get("context", {}).get("exp_context", {})
        
        rationale = "User maintained default/agent-inferred analytical thresholds."
        
        # F1 (ADR-055): a CP3 profile is the EXPLICIT, user-confirmed override
        # path. Each branch stamps its own threshold provenance so the report
        # records that the cutoff was a user choice, not a prompt-inferred value.
        if "strict" in decision.lower():
            exp_context["global_padj"] = 0.01
            exp_context["global_lfc"] = 1.5
            exp_context["lfc_threshold_provenance"] = "cp3_profile:strict"
            rationale = "User enforced Strict analytical thresholds (padj < 0.01, LFC > 1.5) to isolate high-confidence targets."
        elif "standard" in decision.lower():
            exp_context["global_padj"] = 0.05
            exp_context["global_lfc"] = 1.0
            exp_context["lfc_threshold_provenance"] = "cp3_profile:standard"
            rationale = "User enforced Standard analytical thresholds (padj < 0.05, LFC > 1.0) for conventional DE."
        elif "exploratory" in decision.lower() or "tf" in decision.lower():
            exp_context["global_padj"] = 0.05
            exp_context["global_lfc"] = 0.58
            exp_context["lfc_threshold_provenance"] = "cp3_profile:exploratory_tf"
            rationale = "User enforced Exploratory/TF analytical thresholds (padj < 0.05, LFC > 0.58) to capture subtle regulatory signals."
        elif "sensitivity" in decision.lower():
            exp_context["global_padj"] = 0.10
            exp_context["global_lfc"] = 0.58
            exp_context["lfc_threshold_provenance"] = "cp3_profile:sensitivity"
            rationale = "User enforced High Sensitivity analytical thresholds (padj < 0.10, LFC > 0.58) to maximize target discovery."

        self.memory.store_decision(
            decision_id=str(uuid.uuid4())[:8],
            wing_id=experiment_id,
            checkpoint=3,
            question="Statistical Threshold Tuning (padj & log2FC)",
            decision=decision,
            rationale=rationale,
            made_by="user",
        )
        
        # v4.0: inject plan contrasts into design (también en CP3)
        plan_contrasts = plan.get("contrasts")
        if plan_contrasts and exp_context.get("design"):
            exp_context["design"]["plan_contrasts"] = plan_contrasts

        return self._gate_and_dispatch(experiment_id, plan, exp_context)

    # ── Audit gate (v4.1) ────────────────────────────────────────────────

    def _gate_and_dispatch(self, experiment_id: str,
                           plan: dict, exp_context: dict) -> dict:
        """
        Run AuditAgent synchronously. If blocking issues found, publish CP 3.5
        and hold dispatch. Otherwise start the analysis thread immediately.
        """
        from aria.agents.audit_agent import AuditAgent

        session = self._get_session(experiment_id)
        audit = self._register_session_agent(
            experiment_id,
            AuditAgent(memory=self.memory, llm=session.llm_provider),
        )
        with use_egress_policy(session.egress_policy):
            audit_result = audit.run_audit(
                exp_context,
                experiment_id,
                modality_validation=MODALITY_VALIDATION,
            )
        exp_context = self._apply_capability_dispatch_policy(
            exp_context,
            audit_result,
        )
        findings = audit_result.get("findings", [])

        # Store audit findings in exp_context so NarrativeAgent can surface them
        exp_context["audit_findings"] = findings
        if audit_result.get("capability_matrix"):
            exp_context["capability_matrix"] = audit_result["capability_matrix"]
            exp_context["readiness_cards"] = audit_result.get("readiness_cards", {})

        if audit_result["status"] == "blocking":
            # Serialize the blocking issues for the user-facing question
            blocking = [f for f in findings if f["severity"] == "blocking"]
            lines = []
            for i, f in enumerate(blocking, 1):
                lines.append(f"  [{i}] {f['message']}")
                lines.append(f"      → {f['recommendation']}")
            if not blocking:
                matrix = audit_result.get("capability_matrix", {}) or {}
                cards = matrix.get("cards", {}) or {}
                for modality in matrix.get("dispatch", {}).get("requires_ack", []):
                    card = cards.get(modality, {})
                    lines.append(
                        f"  [{len(lines) + 1}] {modality} is "
                        f"{card.get('status', 'yellow')} and requires explicit "
                        "acknowledgement before dispatch."
                    )
                    for finding in card.get("findings", []):
                        lines.append(f"      -> {finding.get('message', '')}")
            issue_text = "\n".join(lines)

            warnings = [f for f in findings if f["severity"] == "warning"]
            warn_text = ""
            if warnings:
                warn_lines = [f"  • {w['message']}" for w in warnings]
                warn_text = "\n\nAdditional warnings:\n" + "\n".join(warn_lines)

            question = (
                "Quality Audit found dispatch readiness issue(s) that may "
                "compromise your results:\n\n"
                f"{issue_text}{warn_text}\n\n"
                f"How would you like to proceed?"
            )

            # Hold the dispatch payload
            session.pending_dispatch = (plan, exp_context)
            self._pending_dispatch[experiment_id] = session.pending_dispatch

            self.publish_escalation(
                experiment_id=experiment_id,
                checkpoint=3.5,
                question=question,
                options=(
                    ["Proceed anyway", "Cancel analysis"]
                    if exp_context.get("modalities")
                    else ["Cancel analysis"]
                ),
                context={"plan": plan, "exp_context": exp_context,
                         "audit_findings": findings},
            )
            return {"status": "audit_blocking", "audit": audit_result}

        # No blocking issues — start analysis immediately
        if findings:
            warn_msgs = "; ".join(f["message"][:80] for f in findings)
            self.publish_status(
                experiment_id,
                f"Audit warnings (non-blocking): {warn_msgs}",
                0.09,
            )

        threading.Thread(
            target=self._dispatch_agents,
            args=(experiment_id, plan, exp_context),
            daemon=True,
        ).start()
        return {"status": "analysis_running", "plan": plan, "exp_context": exp_context}

    @staticmethod
    def _apply_capability_dispatch_policy(
        exp_context: dict,
        audit_result: dict,
    ) -> dict:
        matrix = audit_result.get("capability_matrix") or {}
        blocked = set((matrix.get("dispatch") or {}).get("blocked") or [])
        if not blocked:
            return exp_context

        modalities = exp_context.get("modalities") or {}
        filtered = {
            modality: files
            for modality, files in modalities.items()
            if modality not in blocked
        }
        updated = {**exp_context, "modalities": filtered}
        updated["blocked_modalities_by_capability"] = sorted(
            modality for modality in blocked if modality in modalities
        )
        return updated

    def _after_audit_checkpoint(self, experiment_id: str,
                                decision: str, msg: Message) -> dict:
        """CP 3.5 — user decided whether to proceed despite blocking audit issues."""
        session = self._get_session(experiment_id)
        pending = session.pending_dispatch
        session.pending_dispatch = None
        self._pending_dispatch.pop(experiment_id, None)

        if "cancel" in decision.lower() or pending is None:
            return {"status": "cancelled"}

        plan, exp_context = pending
        if not (exp_context.get("modalities") or {}):
            return {"status": "cancelled", "reason": "no_dispatchable_modalities"}
        self.publish_status(
            experiment_id, "Proceeding with analysis despite audit warnings.", 0.10
        )
        threading.Thread(
            target=self._dispatch_agents,
            args=(experiment_id, plan, exp_context),
            daemon=True,
        ).start()
        return {"status": "analysis_running", "plan": plan, "exp_context": exp_context}

    # ── Despacho de agentes (hilo separado) ─────────────────────────────
    def _dispatch_agents(self, experiment_id: str, plan: dict, exp_context: dict):
        self.publish_status(experiment_id, "Dispatching agents...", 0.1)

        modalities    = exp_context.get("modalities", {})
        intent        = self._experiment_plans.get(experiment_id, {}).get("intent", {})
        agent_results = {}

        if "bulk_RNA_raw" in modalities:
            raw_files = modalities["bulk_RNA_raw"]
            existing  = modalities.get("bulk_RNA", [])
            modalities = {**modalities, "bulk_RNA": list(existing) + list(raw_files)}
            exp_context = {**exp_context, "modalities": modalities}

        blocked_modalities = self._blocked_modalities(modalities)
        if blocked_modalities:
            modalities = {
                m: files for m, files in modalities.items()
                if m not in blocked_modalities
            }
            exp_context = {**exp_context, "modalities": modalities}
            for modality, meta in blocked_modalities.items():
                agent_name = MODALITY_TO_AGENT.get(modality, "unknown")
                agent_results.setdefault(agent_name, {
                    "status": "skipped",
                    "reason": "modality_not_validated",
                    "validation_level": meta.get("level", "unknown"),
                    "modalities": [],
                    "details": [],
                })
                agent_results[agent_name]["modalities"].append(modality)
                agent_results[agent_name]["details"].append(
                    meta.get("reason", "Modality is not validated for dispatch.")
                )
                self.publish_finding(
                    experiment_id,
                    {
                        "summary": (
                            f"{modality} is {meta.get('level', 'unvalidated')} "
                            "and was not dispatched."
                        ),
                        "modality": modality,
                        "validation_level": meta.get("level", "unknown"),
                        "reason": meta.get(
                            "reason", "Modality is not validated for dispatch."
                        ),
                    },
                    Confidence.INSUFFICIENT,
                )

        # P0-3: modalities that run ONLY because their experimental opt-in flag
        # is set are dispatched, but the report is stamped as experimental /
        # not publication-grade so the unvalidated conclusions are never read as
        # final results.
        for modality, meta in self._experimental_modalities(modalities).items():
            exp_context = {
                **exp_context,
                "experimental_modalities": sorted(
                    set(exp_context.get("experimental_modalities", []))
                    | {modality}
                ),
            }
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        f"{modality} ran in EXPERIMENTAL mode "
                        f"(via {meta.get('experimental_env_flag')}); its output "
                        "is NOT publication-grade and is not validated."
                    ),
                    "modality": modality,
                    "validation_level": meta.get("level", "scaffold"),
                    "experimental": True,
                    "reason": meta.get(
                        "reason", "Modality is dispatch-gated; run is experimental."
                    ),
                },
                Confidence.INSUFFICIENT,
            )

        self.publish_status(experiment_id, "Checking computational environment...", 0.05)
        setup_result = self._run_agent(
            agent_name="setup_agent",
            experiment_id=experiment_id,
            context={"exp_context": exp_context, "biological_intent": intent},
        )
        agent_results["setup_agent"] = setup_result

        if setup_result.get("status") == "done":
            genome_config = setup_result.get("genome_config", {})
            exp_context   = {**exp_context, "genome_config": genome_config}

        raw_ingestion_result = self._run_agent(
            agent_name="raw_ingestion_agent",
            experiment_id=experiment_id,
            context={"exp_context": exp_context, "biological_intent": intent},
        )
        agent_results["raw_ingestion_agent"] = raw_ingestion_result
        if raw_ingestion_result.get("status") == "done":
            updates = raw_ingestion_result.get("exp_context_updates", {})
            exp_context = {**exp_context, **updates}
            modalities = exp_context.get("modalities", modalities)
        elif raw_ingestion_result.get("status") == "error":
            self.publish_status(
                experiment_id,
                "Raw ingestion failed before analysis dispatch.",
                1.0,
            )
            self._get_session(experiment_id).agent_results = agent_results
            self._agent_results[experiment_id] = agent_results
            self._present_final_summary(experiment_id, agent_results, [])
            return

        steps = plan.get("steps", []) or self._infer_steps(modalities)
        ordered = self._resolve_execution_order(steps)

        agents_needed = {MODALITY_TO_AGENT[m] for m in modalities if m in MODALITY_TO_AGENT}
        if "bulk_RNA_raw" in modalities:
            agents_needed.add("bulk_rna_agent")

        runnable = sorted(agents_needed - {"integration_agent", "narrative_agent"})

        def _step_hint(step: dict) -> str:
            text = (step.get("agent", "") + " " + str(step.get("analysis", "")) + " " + str(step.get("description", ""))).lower()
            for a in runnable:
                if a in text or a.replace("_agent","") in text:
                    return a
            return ""

        ordered_by_plan = []
        for step in ordered:
            hint = _step_hint(step)
            if hint and hint not in ordered_by_plan:
                ordered_by_plan.append(hint)
        for a in runnable:
            if a not in ordered_by_plan:
                ordered_by_plan.append(a)

        agents_run: set = set()
        n_agents = max(len(ordered_by_plan), 1)

        for i, agent_name in enumerate(ordered_by_plan):
            if agent_name in agents_run:
                continue

            self.publish_status(experiment_id, f"Running {agent_name}...", 0.15 + (i / n_agents) * 0.50)
            result = self._run_agent(
                agent_name=agent_name,
                experiment_id=experiment_id,
                context={"exp_context": exp_context, "biological_intent": intent},
            )
            agent_results[agent_name] = result
            agents_run.add(agent_name)
            log.info(f"{agent_name}: {result.get('status', '?')}")

        n_mods = len([m for m in modalities if m in MODALITY_TO_AGENT])
        integration_requested = n_mods >= 2 or plan.get("integration_needed")
        if integration_requested and not self._integration_dispatch_enabled():
            meta = INTEGRATION_VALIDATION
            agent_results["integration_agent"] = {
                "status": "skipped",
                "reason": "agent_not_validated",
                "validation_level": meta.get("level", "unknown"),
                "details": meta.get(
                    "reason", "IntegrationAgent is not validated for dispatch."
                ),
            }
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        "IntegrationAgent is scaffolded and was not dispatched."
                    ),
                    "agent": "integration_agent",
                    "validation_level": meta.get("level", "unknown"),
                    "reason": meta.get(
                        "reason", "IntegrationAgent is not validated for dispatch."
                    ),
                },
                Confidence.INSUFFICIENT,
            )
        elif integration_requested:
            self.publish_status(experiment_id, "Running IntegrationAgent...", 0.72)
            result = self._run_agent(
                agent_name="integration_agent",
                experiment_id=experiment_id,
                context={
                    "exp_context":       exp_context,
                    "biological_intent": intent,
                    "rna_results":       agent_results.get("scrna_agent", agent_results.get("bulk_rna_agent", agent_results.get("rna_agent", {}))),
                    "chromatin_results": agent_results.get("chromatin_agent", {}),
                    "hic_results":       agent_results.get("genome_arch_agent", {}),
                },
            )
            agent_results["integration_agent"] = result

        self.publish_status(experiment_id, "Generating report...", 0.88)
        findings = self._get_session(experiment_id).message_bus.get_findings(
            experiment_id
        )

        narrative_result = self._run_agent(
            agent_name="narrative_agent",
            experiment_id=experiment_id,
            context={
                "exp_context":       exp_context,
                "biological_intent": intent,
                "agent_results":     agent_results,
                "findings":          [m.payload for m in findings],
            },
        )
        agent_results["narrative_agent"] = narrative_result
        self._get_session(experiment_id).agent_results = agent_results
        self._agent_results[experiment_id] = agent_results

        self._present_final_summary(experiment_id, agent_results, findings)
        self.publish_status(experiment_id, "Analysis complete.", 1.0)

    def _run_agent(self, agent_name: str, experiment_id: str, context: dict) -> dict:
        if agent_name not in AGENT_REGISTRY:
            return {"status": "error", "error_type": "UnknownAgent", "details": f"'{agent_name}' not in AGENT_REGISTRY."}
        entry = AGENT_REGISTRY[agent_name]
        try:
            import importlib
            module = importlib.import_module(entry["module"])
            AgentClass = getattr(module, entry["class"])
            session = self._get_session(experiment_id)
            agent = self._register_session_agent(
                experiment_id,
                AgentClass(memory=self.memory, llm=session.llm_provider),
            )
            with use_egress_policy(session.egress_policy):
                return agent.run(experiment_id, context)
        except Exception as e:
            log.error(f"{agent_name} raised: {e}", exc_info=True)
            return {"status": "error", "error_type": type(e).__name__, "details": str(e), "agent": agent_name}

    @staticmethod
    def _experimental_override_active(meta: dict) -> bool:
        """A dispatch-gated modality may carry an explicit experimental opt-in
        env flag (P0-3). The override is active only when that flag is set to a
        truthy value; it never silently enables a modality."""
        flag = meta.get("experimental_env_flag")
        if not flag:
            return False
        return os.environ.get(flag, "").strip().lower() in {"1", "true", "yes", "on"}

    @classmethod
    def _blocked_modalities(cls, modalities: dict) -> dict:
        blocked = {}
        for modality in modalities:
            meta = MODALITY_VALIDATION.get(modality, {})
            if meta and not meta.get("dispatch_enabled", True):
                # An explicit experimental opt-in unblocks dispatch but the run
                # is stamped experimental (see _experimental_modalities).
                if cls._experimental_override_active(meta):
                    continue
                blocked[modality] = meta
        return blocked

    @classmethod
    def _experimental_modalities(cls, modalities: dict) -> dict:
        """Modalities that are dispatch-gated by default but allowed to run only
        because their experimental env flag is set. These produce
        not-publication-grade output and must be stamped as such."""
        experimental = {}
        for modality in modalities:
            meta = MODALITY_VALIDATION.get(modality, {})
            if (meta and not meta.get("dispatch_enabled", True)
                    and cls._experimental_override_active(meta)):
                experimental[modality] = meta
        return experimental

    @staticmethod
    def _integration_dispatch_enabled() -> bool:
        return bool(INTEGRATION_VALIDATION.get("dispatch_enabled", True))

    def _resolve_execution_order(self, steps: list) -> list:
        remaining, ordered, completed = list(steps), [], set()
        for _ in range(len(steps) + 1):
            if not remaining: break
            for step in remaining[:]:
                if set(step.get("depends_on", [])) <= completed:
                    ordered.append(step)
                    completed.add(step.get("agent", ""))
                    remaining.remove(step)
        ordered.extend(remaining) 
        return ordered

    def _infer_steps(self, modalities: dict) -> list:
        steps, seen, order = [], set(), 1
        for m in modalities:
            agent = MODALITY_TO_AGENT.get(m)
            if agent and agent not in seen:
                steps.append({"order": order, "agent": agent, "analysis": f"Analyze {m}", "depends_on": [], "can_parallel": False})
                seen.add(agent)
                order += 1
        return steps

    def _present_final_summary(self, experiment_id: str, agent_results: dict, findings: list):
        n_done   = sum(1 for r in agent_results.values() if r.get("status") == "done")
        n_errors = sum(1 for r in agent_results.values() if r.get("status") == "error")
        report   = agent_results.get("narrative_agent", {}).get("report_path")

        lines = [f"Analysis complete.\n", f"Agents completed:  {n_done}/{len(agent_results)}", f"Findings reported: {len(findings)}"]
        if n_errors: lines.append(f"Errors in:         {', '.join(k for k, r in agent_results.items() if r.get('status') == 'error')}")
        if report: lines.append(f"Report:            {report}")

        self.publish_escalation(
            experiment_id=experiment_id, checkpoint=5,
            question="\n".join(lines), options=["Accept report", "Request revision", "Export methods"],
            context={"n_findings": len(findings), "report_path": report},
        )

        # Detach the per-experiment log handler so the file is flushed and
        # the next experiment starts with a fresh handler list.
        session = self._get_session(experiment_id)
        handler = session.log_handler
        mirror = self._log_handlers.pop(experiment_id, None)
        if handler is None:
            handler = mirror
        session.log_handler = None
        if handler is not None:
            try:
                from aria.utils.logging_setup import detach_experiment_log
                detach_experiment_log(handler)
            except Exception:
                pass

    def _parse_question(self, user_question: str, llm=None) -> dict:
        prompt = build_untrusted_prompt(
            task="Extract the biological analysis intent from the submitted question.",
            fields=[PromptDataField(
                name="user_question",
                value=user_question,
                kind="user_text",
                source="user",
            )],
            response_contract='''Return JSON with this schema:
{
  "analysis_type": "differential|regulatory|structural|temporal|cell_type|integration",
  "biological_entities": ["genes, TFs, cell types, conditions mentioned"],
  "comparison": "what vs what (if applicable)",
  "key_modalities_needed": ["RNA", "ATAC", "HiC", "ChIP"],
  "complexity": "simple|moderate|complex",
  "summary": "one sentence biological question summary"
}''',
        )
        return self.think_structured(
            prompt,
            system=system_with_untrusted_boundary(ORCHESTRATOR_SYSTEM),
            schema_hint="Return analysis intent as JSON.",
            llm=llm,
        )

    def _design_analysis_plan(self, experiment_id: str, exp_context: dict) -> dict:
        prompt = build_untrusted_prompt(
            task="Design the analysis pipeline using the submitted context.",
            fields=[
                PromptDataField(
                    "available_modalities",
                    list(exp_context.get("modalities", {}).keys()),
                    "identifier_list",
                    "data_audit",
                ),
                PromptDataField("user_question", exp_context.get(
                    "user_question", ""
                ), "user_text", "user"),
                PromptDataField("organism", exp_context.get("organism"),
                                "metadata", "data_audit"),
                PromptDataField("genome", exp_context.get("genome"),
                                "metadata", "data_audit"),
                PromptDataField("is_multimodal", exp_context.get(
                    "is_multimodal"
                ), "metadata", "data_audit"),
            ],
            response_contract='''Return JSON with this schema:
{
  "steps": [{"order": 1, "agent": "rna_agent", "analysis": "description", "depends_on": [], "can_parallel": false}],
  "contrasts": [{"numerator": "groupA", "denominator": "groupB"}],
  "integration_needed": true,
  "integration_type": "WNN|MOFA|none",
  "estimated_complexity": "low|medium|high",
  "rationale": "brief explanation"
}''',
        )
        try:
            session = self._get_session(experiment_id)
            with use_egress_policy(session.egress_policy):
                plan = self.think_structured(
                    prompt,
                    system=system_with_untrusted_boundary(ORCHESTRATOR_SYSTEM),
                    schema_hint="Return analysis plan as JSON.",
                    llm=session.llm_provider,
                )
            if not plan or "steps" not in plan:
                raise ValueError("Invalid plan from LLM")
        except Exception as e:
            log.warning(f"LLM plan generation failed: {e}. Using fallback plan.")
            plan = self._fallback_plan(exp_context)
        
        self._experiment_plans[experiment_id]["plan"] = plan
        return plan

    def _run_design_intelligence(self, exp_context: dict, intent: dict) -> dict:
        try:
            from aria.agents.design_intelligence import DesignIntelligence
            return DesignIntelligence().evaluate(exp_context, intent)
        except Exception as e:
            log.warning(f"Design intelligence failed: {e}", exc_info=True)
            return {
                "status": "error",
                "recommended": [],
                "optional": [],
                "unsupported": [],
                "warnings": [f"Design intelligence unavailable: {e}"],
            }

    def _fallback_plan(self, exp_context: dict) -> dict:
        """Generate a minimal plan if LLM fails."""
        modalities = exp_context.get("modalities", {})
        steps = []
        order = 1
        for mod in modalities:
            if mod in MODALITY_TO_AGENT:
                steps.append({
                    "order": order,
                    "agent": MODALITY_TO_AGENT[mod],
                    "analysis": f"Analyze {mod}",
                    "depends_on": [],
                    "can_parallel": False
                })
                order += 1
        return {
            "steps": steps,
            "integration_needed": len(steps) > 1,
            "integration_type": "WNN" if len(steps) > 1 else "none",
            "estimated_complexity": "low",
            "rationale": "Fallback plan (LLM unavailable)",
            "contrasts": []
        }
    
    def _format_plan_summary(self, plan: dict) -> str:
        lines = ["ARIA Analysis Plan:\n"]
        for step in plan.get("steps", []):
            lines.append(
                f"  Step {step['order']}: [{step['agent']}] {step['analysis']}"
                f"{' (parallel)' if step.get('can_parallel') else ''}"
            )
        has_bulk = any(
            "bulk_rna" in s.get("agent", "") or "rna" in s.get("agent", "")
            for s in plan.get("steps", [])
        )
        if has_bulk:
            lfc  = plan.get("lfc_threshold",  1.0)
            padj = plan.get("padj_threshold", 0.05)
            lines.append(f"\n  DE thresholds: padj < {padj}, |log2FC| > {lfc}")
        if plan.get("integration_needed"):
            lines.append(f"\n  Integration: {plan.get('integration_type', 'TBD')}")
        try:
            from aria.agents.design_intelligence import format_design_intelligence
            di_block = format_design_intelligence(plan.get("design_intelligence", {}))
            if di_block:
                lines.append(di_block)
        except Exception:
            pass
        lines.append(
            f"\n  Complexity: {plan.get('estimated_complexity', '?')}"
            f"\n  {plan.get('rationale', '')}"
            "\n\nChoose option 1 to run only recommended analyses, or option 2 "
            "to add optional supported analyses. Unsupported/not-recommended "
            "items will not run unless the plan is modified explicitly."
            f"\nProceed with this plan?"
        )
        return "\n".join(lines)
