"""
ARIA OrchestratorAgent
----------------------
The central coordinator. Receives the user's biological question,
designs the analysis plan, and dispatches specialized agents.

Flow:
  1. Parse biological question → structured analysis plan
  2. Delegate to DataAuditAgent (always first)
  3. Wait for Checkpoint 1 confirmation
  4. Present analysis plan → Checkpoint 2
  5. DISPATCH modality agents in dependency order   ← the missing piece
  6. Collect findings from all agents
  7. Trigger IntegrationAgent if multimodal
  8. Dispatch NarrativeAgent for final report
  9. Present Checkpoint 5 (final review)
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
- preprocessing_agent: raw FASTQs → fastp → STAR → featureCounts
- scrna_agent: single-cell RNA-seq (QC, clustering, annotation, DE)
- bulk_rna_agent: bulk RNA-seq (DESeq2, pathway enrichment, plots)
- chromatin_agent: handles ATAC, ChIP, CUT&RUN, CUT&TAG
- genome_arch_agent: handles HiC, compartments, TADs, loops
- integration_agent: multimodal integration (needs 2+ modalities)
- narrative_agent: generates the final report

Always think about the biology first, then the methods.
""".strip()


# ── Agent registry ─────────────────────────────────────────────────────────────
# Lazy imports prevent circular dependencies at module load time.

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
    "preprocessing_agent": {
        "module": "aria.agents.preprocessing_agent",
        "class":  "PreprocessingAgent",
    },
    # Legacy alias — kept for backward compatibility
    "rna_agent": {
        "module": "aria.agents.scrna_agent",
        "class":  "scRNAAgent",
    },
}

MODALITY_TO_AGENT = {
    "scRNA":        "scrna_agent",
    "bulk_RNA_raw": "preprocessing_agent",  # raw FASTQs → preprocessing first
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
        super().__init__(memory, api_key)
        self.llm                   = llm or LLMProvider.from_config()
        self._pending_checkpoints  = {}
        self._experiment_plans     = {}
        self._agent_results        = {}

    # ── Entry point ─────────────────────────────────────────────────────────

    def run(self, experiment_id: str, context: dict) -> dict:
        self.publish_status(experiment_id,
                            "ARIA starting analysis...", 0.0)

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

    # ── Checkpoint routing ──────────────────────────────────────────────────

    def on_checkpoint_resolved(self, message_id: str,
                                user_decision: str,
                                experiment_id: str) -> dict:
        bus.resolve_checkpoint(message_id, {"choice": user_decision})

        resolved_msg = next(
            (m for m in bus.get_log() if m.id == message_id), None
        )
        if not resolved_msg:
            return {"status": "error", "message": "Checkpoint not found"}

        cp = resolved_msg.checkpoint
        if cp == 1:
            return self._after_checkpoint_1(
                experiment_id, user_decision, resolved_msg)
        elif cp == 2:
            return self._after_checkpoint_2(
                experiment_id, user_decision, resolved_msg)
        return {"status": "ok"}

    def _after_checkpoint_1(self, experiment_id: str,
                             decision: str, msg: Message) -> dict:
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

        plan = self._design_analysis_plan(experiment_id, exp_context)

        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=2,
            question=self._format_plan_summary(plan),
            options=["Confirm and run", "Modify plan", "Cancel"],
            context={"plan": plan, "exp_context": exp_context},
        )

        return {"status": "plan_ready", "plan": plan}

    def _after_checkpoint_2(self, experiment_id: str,
                             decision: str, msg: Message) -> dict:
        """
        After plan confirmation — dispatch analysis agents.
        Runs in a background thread so the TUI remains responsive.
        """
        if "cancel" in decision.lower():
            return {"status": "cancelled"}

        plan        = msg.payload.get("context", {}).get("plan", {})
        exp_context = msg.payload.get("context", {}).get("exp_context", {})

        self.memory.store_decision(
            decision_id=str(uuid.uuid4())[:8],
            wing_id=experiment_id,
            checkpoint=2,
            question="Analysis plan confirmed",
            decision=decision,
            rationale="User approved analysis plan",
            made_by="user",
        )

        threading.Thread(
            target=self._dispatch_agents,
            args=(experiment_id, plan, exp_context),
            daemon=True,
        ).start()

        return {
            "status":        "analysis_running",
            "plan":          plan,
            "exp_context":   exp_context,
        }

    # ── Dispatcher ──────────────────────────────────────────────────────────

    def _dispatch_agents(self, experiment_id: str,
                          plan: dict, exp_context: dict):
        """
        Instantiate and run each agent in dependency order.

        This is the piece that was missing — the actual connection
        between the LLM-designed plan and the agents that execute it.

        Order:
          1. Modality agents (RNA, Chromatin, HiC)
          2. IntegrationAgent (if 2+ modalities)
          3. NarrativeAgent (always)
        """
        self.publish_status(experiment_id,
                            "Dispatching agents...", 0.1)

        modalities    = exp_context.get("modalities", {})
        intent        = self._experiment_plans.get(
                            experiment_id, {}
                        ).get("intent", {})
        agent_results = {}

        # ── Step 0: SetupAgent — provision environment before analysis ────
        # Installs conda envs, downloads genome, builds STAR index.
        # Transparent to the user — runs silently if already set up.
        self.publish_status(experiment_id,
                            "Checking computational environment...", 0.05)
        setup_result = self._run_agent(
            agent_name="setup_agent",
            experiment_id=experiment_id,
            context={
                "exp_context":       exp_context,
                "biological_intent": intent,
            },
        )
        agent_results["setup_agent"] = setup_result

        # Inject genome_config into exp_context for downstream agents
        if setup_result.get("status") == "done":
            genome_config = setup_result.get("genome_config", {})
            exp_context   = {**exp_context, "genome_config": genome_config}

        steps = plan.get("steps", []) or self._infer_steps(modalities)
        ordered = self._resolve_execution_order(steps)

        agents_needed = {
            MODALITY_TO_AGENT[m]
            for m in modalities
            if m in MODALITY_TO_AGENT
        }

        # If raw FASTQs detected, preprocessing must run before bulk DE
        has_raw = "bulk_RNA_raw" in modalities
        if has_raw:
            agents_needed.add("preprocessing_agent")
            agents_needed.add("bulk_rna_agent")

        # ── Modality agents ──────────────────────────────────────────────
        # Deduplicate: each agent runs at most once per experiment.
        # The LLM plan may list multiple steps under the same agent
        # (e.g. "FASTQ QC", "Alignment", "DE" all → bulk_rna_agent).
        # The agent itself handles its internal pipeline stages.
        agents_run: set = set()
        n_steps = max(len(ordered), 1)

        for i, step in enumerate(ordered):
            agent_name = step.get("agent", "")

            if agent_name not in agents_needed:
                continue
            if agent_name in ("integration_agent", "narrative_agent"):
                continue
            if agent_name in agents_run:
                # Agent already executed — skip duplicate steps
                log.debug(f"Skipping duplicate step for {agent_name}")
                continue

            self.publish_status(
                experiment_id,
                f"Running {agent_name}...",
                0.15 + (i / n_steps) * 0.50,
            )

            result = self._run_agent(
                agent_name=agent_name,
                experiment_id=experiment_id,
                context={
                    "exp_context":       exp_context,
                    "biological_intent": intent,
                },
            )
            agent_results[agent_name] = result
            agents_run.add(agent_name)
            log.info(f"{agent_name}: {result.get('status', '?')}")

        # ── IntegrationAgent ─────────────────────────────────────────────
        n_mods = len([m for m in modalities if m in MODALITY_TO_AGENT])
        if n_mods >= 2 or plan.get("integration_needed"):
            self.publish_status(experiment_id,
                                "Running IntegrationAgent...", 0.72)
            result = self._run_agent(
                agent_name="integration_agent",
                experiment_id=experiment_id,
                context={
                    "exp_context":       exp_context,
                    "biological_intent": intent,
                    "rna_results":       agent_results.get("scrna_agent",
                                        agent_results.get("bulk_rna_agent",
                                        agent_results.get("rna_agent", {}))),
                    "chromatin_results": agent_results.get("chromatin_agent", {}),
                    "hic_results":       agent_results.get("genome_arch_agent", {}),
                },
            )
            agent_results["integration_agent"] = result

        # ── NarrativeAgent ────────────────────────────────────────────────
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

        # ── Checkpoint 5: final summary ───────────────────────────────────
        self._present_final_summary(
            experiment_id, agent_results, findings
        )
        self.publish_status(experiment_id, "Analysis complete.", 1.0)

    def _run_agent(self, agent_name: str,
                   experiment_id: str,
                   context: dict) -> dict:
        """
        Instantiate and run one agent by name.
        Uses lazy import to avoid circular dependencies.
        Returns a structured error dict if the agent is unavailable.
        """
        if agent_name not in AGENT_REGISTRY:
            return {
                "status":     "error",
                "error_type": "UnknownAgent",
                "details":    f"'{agent_name}' not in AGENT_REGISTRY.",
            }

        entry = AGENT_REGISTRY[agent_name]

        try:
            import importlib
            module     = importlib.import_module(entry["module"])
            AgentClass = getattr(module, entry["class"])
            agent      = AgentClass(memory=self.memory, llm=self.llm)
            return agent.run(experiment_id, context)

        except ModuleNotFoundError as e:
            log.warning(f"Module not found for {agent_name}: {e}")
            return {"status": "skipped", "reason": str(e),
                    "agent": agent_name}

        except Exception as e:
            log.error(f"{agent_name} raised: {e}", exc_info=True)
            return {"status": "error", "error_type": type(e).__name__,
                    "details": str(e), "agent": agent_name}

    # ── Helpers ─────────────────────────────────────────────────────────────

    def _resolve_execution_order(self, steps: list) -> list:
        """Topological sort by depends_on."""
        remaining = list(steps)
        ordered   = []
        completed = set()

        for _ in range(len(steps) + 1):
            if not remaining:
                break
            for step in remaining[:]:
                if set(step.get("depends_on", [])) <= completed:
                    ordered.append(step)
                    completed.add(step.get("agent", ""))
                    remaining.remove(step)

        ordered.extend(remaining)  # add anything with circular deps
        return ordered

    def _infer_steps(self, modalities: dict) -> list:
        """Fallback step list when LLM plan fails."""
        steps = []
        seen  = set()
        order = 1
        for m in modalities:
            agent = MODALITY_TO_AGENT.get(m)
            if agent and agent not in seen:
                steps.append({
                    "order":       order,
                    "agent":       agent,
                    "analysis":    f"Analyze {m}",
                    "depends_on":  [],
                    "can_parallel": False,
                })
                seen.add(agent)
                order += 1
        return steps

    def _present_final_summary(self, experiment_id: str,
                                agent_results: dict,
                                findings: list):
        n_done   = sum(1 for r in agent_results.values()
                       if r.get("status") == "done")
        n_errors = sum(1 for r in agent_results.values()
                       if r.get("status") == "error")
        report   = agent_results.get("narrative_agent", {}).get(
                        "report_path")

        lines = [
            f"Analysis complete.",
            f"",
            f"Agents completed:  {n_done}/{len(agent_results)}",
            f"Findings reported: {len(findings)}",
        ]
        if n_errors:
            errors = [k for k, r in agent_results.items()
                      if r.get("status") == "error"]
            lines.append(f"Errors in:         {', '.join(errors)}")
        if report:
            lines.append(f"Report:            {report}")

        self.publish_escalation(
            experiment_id=experiment_id,
            checkpoint=5,
            question="\n".join(lines),
            options=["Accept report", "Request revision", "Export methods"],
            context={"n_findings": len(findings), "report_path": report},
        )

    # ── Analysis planning ────────────────────────────────────────────────────

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
        return self.think_structured(
            prompt,
            system=ORCHESTRATOR_SYSTEM,
            schema_hint="Return analysis intent as JSON.",
        )

    def _design_analysis_plan(self, experiment_id: str,
                               exp_context: dict) -> dict:
        modalities    = list(exp_context.get("modalities", {}).keys())
        user_question = exp_context.get("user_question", "")

        prompt = f"""
Available modalities: {modalities}
User question: "{user_question}"
Organism: {exp_context.get("organism")}
Genome: {exp_context.get("genome")}
Is multimodal: {exp_context.get("is_multimodal")}

Design the analysis pipeline. Return JSON:
{{
  "steps": [
    {{
      "order": 1,
      "agent": "rna_agent",
      "analysis": "description",
      "depends_on": [],
      "can_parallel": false
    }}
  ],
  "integration_needed": true,
  "integration_type": "WNN|MOFA|none",
  "estimated_complexity": "low|medium|high",
  "rationale": "brief explanation"
}}
"""
        plan = self.think_structured(
            prompt,
            system=ORCHESTRATOR_SYSTEM,
            schema_hint="Return analysis plan as JSON.",
        )
        self._experiment_plans[experiment_id]["plan"] = plan
        return plan

    def _format_plan_summary(self, plan: dict) -> str:
        steps = plan.get("steps", [])
        lines = ["ARIA Analysis Plan:\n"]
        for step in steps:
            p = " (parallel)" if step.get("can_parallel") else ""
            lines.append(f"  Step {step['order']}: "
                         f"[{step['agent']}] {step['analysis']}{p}")
        if plan.get("integration_needed"):
            lines.append(
                f"\n  Integration: {plan.get('integration_type', 'TBD')}")
        lines.append(
            f"\n  Complexity: {plan.get('estimated_complexity', '?')}")
        lines.append(f"\n  {plan.get('rationale', '')}")
        lines.append("\nProceed with this plan?")
        return "\n".join(lines)

    def receive(self, message: Message):
        if message.type == MessageType.FINDING:
            exp_id = message.experiment_id
            if exp_id in self._experiment_plans:
                self._experiment_plans[exp_id]["findings"].append(
                    message.payload)
        elif message.type == MessageType.ESCALATION:
            self._pending_checkpoints[message.id] = message
