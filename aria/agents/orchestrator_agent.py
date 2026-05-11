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
import threading
import uuid
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Message, MessageType, Confidence, bus
from aria.llm.provider import LLMProvider
from aria.memory.memory import ARIAMemory

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


class OrchestratorAgent(BaseAgent):

    name        = "orchestrator"
    description = "Central coordinator — routes tasks to specialized agents."

    MODALITY_TO_AGENT = MODALITY_TO_AGENT

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider = None,
                 api_key: str = None):
        super().__init__(memory, llm=llm, api_key=api_key)
        self._pending_checkpoints  = {}
        self._experiment_plans     = {}
        self._agent_results        = {}
        # v4.0: active DesignAgent instance (state machine)
        self._active_design_agent  = None
        # v4.1: pending dispatch held while audit CP 3.5 is open
        self._pending_dispatch     = {}

    def run(self, experiment_id: str, context: dict) -> dict:
        self.publish_status(experiment_id, "ARIA starting analysis...", 0.0)
        intent = self._parse_question(context["user_question"])

        self._experiment_plans[experiment_id] = {
            "user_question": context["user_question"],
            "intent":        intent,
            "context":       context,
            "status":        "auditing",
            "findings":      [],
        }

        return {"status": "started", "intent": intent, "next": "data_audit"}

    def run_audit(self, experiment_id: str) -> dict:
        from aria.agents.data_audit_agent import DataAuditAgent
        plan = self._experiment_plans.get(experiment_id, {})
        return DataAuditAgent(self.memory).run(experiment_id, plan["context"])

    # ── Callback principal para checkpoints resueltos ───────────────────
    def on_checkpoint_resolved(self, message_id: str,
                                user_decision: str,
                                experiment_id: str) -> dict:
        bus.resolve_checkpoint(message_id, {"choice": user_decision})

        resolved_msg = next((m for m in bus.get_log() if m.id == message_id), None)
        if not resolved_msg:
            return {"status": "error", "message": "Checkpoint not found"}

        cp = resolved_msg.checkpoint

        # ── v4.0: Delegar checkpoints de diseño al DesignAgent ──────────
        if cp in (2.1, 2.2, 2.3, 2.4, 2.5, 2.6) and self._active_design_agent is not None:
            result = self._active_design_agent.handle_user_response(
                experiment_id=experiment_id,
                checkpoint_num=cp,
                choice=user_decision,
            )
            if result.get("status") == "cancelled":
                self._active_design_agent = None
                return {"status": "cancelled", "reason": "Design cancelled by user"}

            elif result.get("status") == "done":
                # Diseño completado → enriquecer contexto y pasar a CP2
                design = result["design"]
                exp_context = self._experiment_plans.get(experiment_id, {}).get("exp_context", {})
                exp_context["organism"] = design["organism"]
                exp_context["genome"]   = design["genome"]
                exp_context["design"]   = design
                self._experiment_plans[experiment_id]["exp_context"] = exp_context
                self._active_design_agent = None

                plan = self._design_analysis_plan(experiment_id, exp_context)
                # Annotate plan with resolved thresholds — single source of truth
                # so the preview and the actual run always agree.
                intent = self._experiment_plans.get(experiment_id, {}).get("intent", {})
                from aria.agents.bulk_rna_agent import _infer_lfc_threshold
                plan["lfc_threshold"]  = _infer_lfc_threshold(intent)
                plan["padj_threshold"] = 0.05
                self.publish_escalation(
                    experiment_id=experiment_id,
                    checkpoint=2,
                    question=self._format_plan_summary(plan),
                    options=["Confirm and run", "Modify plan", "Cancel"],
                    context={"plan": plan, "exp_context": exp_context},
                )
                return {"status": "plan_ready", "plan": plan}

            else:
                # El DesignAgent ya publicó el siguiente checkpoint
                return {"status": "design_in_progress", "next_step": result.get("step")}

        # ── Checkpoints normales ────────────────────────────────────────
        if cp == 1:
            return self._after_checkpoint_1(experiment_id, user_decision, resolved_msg)
        elif cp == 2:
            return self._after_checkpoint_2(experiment_id, user_decision, resolved_msg)
        elif cp == 3:
            return self._after_checkpoint_3(experiment_id, user_decision, resolved_msg)
        elif cp == 3.5:
            return self._after_audit_checkpoint(experiment_id, user_decision, resolved_msg)

        return {"status": "ok"}

    # ── Checkpoint 1: Auditoría completada ──────────────────────────────
    def _after_checkpoint_1(self, experiment_id: str, decision: str, msg: Message) -> dict:
        if "cancel" in decision.lower():
            return {"status": "cancelled"}

        exp_context = msg.payload.get("context", {}).get("exp_context", {})
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
        design_agent = DesignAgent(memory=self.memory, llm=self.llm)
        result = design_agent.start_design(
            experiment_id=experiment_id,
            exp_context=exp_context,
            biological_intent=self._experiment_plans.get(experiment_id, {}).get("intent", {}),
        )

        if result.get("status") == "failed":
            return {"status": "cancelled", "reason": result.get("reason", "Design failed")}

        # Guardar el agente activo para seguir recibiendo sus checkpoints
        self._active_design_agent = design_agent
        self._experiment_plans[experiment_id]["exp_context"] = exp_context

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

        self.memory.store_decision(
            decision_id=str(uuid.uuid4())[:8],
            wing_id=experiment_id,
            checkpoint=2,
            question="Analysis plan confirmed",
            decision=decision,
            rationale="User approved analysis plan directly.",
            made_by="user",
        )

        # Propagate resolved thresholds into exp_context so BulkRNAAgent
        # uses the same values shown in the plan preview.
        exp_context.setdefault("global_padj", plan.get("padj_threshold", 0.05))
        exp_context.setdefault("global_lfc",  plan.get("lfc_threshold",  1.0))

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
        
        if "strict" in decision.lower():
            exp_context["global_padj"] = 0.01
            exp_context["global_lfc"] = 1.5
            rationale = "User enforced Strict analytical thresholds (padj < 0.01, LFC > 1.5) to isolate high-confidence targets."
        elif "standard" in decision.lower():
            exp_context["global_padj"] = 0.05
            exp_context["global_lfc"] = 1.0
            rationale = "User enforced Standard analytical thresholds (padj < 0.05, LFC > 1.0) for conventional DE."
        elif "exploratory" in decision.lower() or "tf" in decision.lower():
            exp_context["global_padj"] = 0.05
            exp_context["global_lfc"] = 0.58
            rationale = "User enforced Exploratory/TF analytical thresholds (padj < 0.05, LFC > 0.58) to capture subtle regulatory signals."
        elif "sensitivity" in decision.lower():
            exp_context["global_padj"] = 0.10
            exp_context["global_lfc"] = 0.58
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

        audit = AuditAgent(memory=self.memory, llm=self.llm)
        audit_result = audit.run_audit(exp_context, experiment_id)
        findings = audit_result.get("findings", [])

        # Store audit findings in exp_context so NarrativeAgent can surface them
        exp_context["audit_findings"] = findings

        if audit_result["status"] == "blocking":
            # Serialize the blocking issues for the user-facing question
            blocking = [f for f in findings if f["severity"] == "blocking"]
            lines = []
            for i, f in enumerate(blocking, 1):
                lines.append(f"  [{i}] {f['message']}")
                lines.append(f"      → {f['recommendation']}")
            issue_text = "\n".join(lines)

            warnings = [f for f in findings if f["severity"] == "warning"]
            warn_text = ""
            if warnings:
                warn_lines = [f"  • {w['message']}" for w in warnings]
                warn_text = "\n\nAdditional warnings:\n" + "\n".join(warn_lines)

            question = (
                f"Quality Audit found {len(blocking)} blocking issue(s) "
                f"that may compromise your results:\n\n"
                f"{issue_text}{warn_text}\n\n"
                f"How would you like to proceed?"
            )

            # Hold the dispatch payload
            self._pending_dispatch[experiment_id] = (plan, exp_context)

            self.publish_escalation(
                experiment_id=experiment_id,
                checkpoint=3.5,
                question=question,
                options=["Proceed anyway", "Cancel analysis"],
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

    def _after_audit_checkpoint(self, experiment_id: str,
                                decision: str, msg: Message) -> dict:
        """CP 3.5 — user decided whether to proceed despite blocking audit issues."""
        pending = self._pending_dispatch.pop(experiment_id, None)

        if "cancel" in decision.lower() or pending is None:
            return {"status": "cancelled"}

        plan, exp_context = pending
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
        if n_mods >= 2 or plan.get("integration_needed"):
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
        findings = bus.get_findings(experiment_id)

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
            agent = AgentClass(memory=self.memory, llm=self.llm)
            return agent.run(experiment_id, context)
        except Exception as e:
            log.error(f"{agent_name} raised: {e}", exc_info=True)
            return {"status": "error", "error_type": type(e).__name__, "details": str(e), "agent": agent_name}

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

    def _parse_question(self, user_question: str) -> dict:
        prompt = f"""
User question: "{user_question}"
Extract the biological analysis intent. Return JSON:
{{
  "analysis_type": "differential|regulatory|structural|temporal|cell_type|integration",
  "biological_entities": ["genes, TFs, cell types, conditions mentioned"],
  "comparison": "what vs what (if applicable)",
  "key_modalities_needed": ["RNA", "ATAC", "HiC", "ChIP"],
  "complexity": "simple|moderate|complex",
  "summary": "one sentence biological question summary"
}}
"""
        return self.think_structured(prompt, system=ORCHESTRATOR_SYSTEM, schema_hint="Return analysis intent as JSON.")

    def _design_analysis_plan(self, experiment_id: str, exp_context: dict) -> dict:
        prompt = f"""
Available modalities: {list(exp_context.get("modalities", {}).keys())}
User question: "{exp_context.get("user_question", "")}"
Organism: {exp_context.get("organism")}
Genome: {exp_context.get("genome")}
Is multimodal: {exp_context.get("is_multimodal")}

Design the analysis pipeline. Return JSON:
{{
  "steps": [{{"order": 1, "agent": "rna_agent", "analysis": "description", "depends_on": [], "can_parallel": false}}],
  "contrasts": [{{"numerator": "groupA", "denominator": "groupB"}}],
  "integration_needed": true,
  "integration_type": "WNN|MOFA|none",
  "estimated_complexity": "low|medium|high",
  "rationale": "brief explanation"
}}
"""
        try:
            plan = self.think_structured(prompt, system=ORCHESTRATOR_SYSTEM,
                                         schema_hint="Return analysis plan as JSON.")
            if not plan or "steps" not in plan:
                raise ValueError("Invalid plan from LLM")
        except Exception as e:
            log.warning(f"LLM plan generation failed: {e}. Using fallback plan.")
            plan = self._fallback_plan(exp_context)
        
        self._experiment_plans[experiment_id]["plan"] = plan
        return plan

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
        lines.append(
            f"\n  Complexity: {plan.get('estimated_complexity', '?')}"
            f"\n  {plan.get('rationale', '')}"
            f"\nProceed with this plan?"
        )
        return "\n".join(lines)