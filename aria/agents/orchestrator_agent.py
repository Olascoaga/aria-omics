"""
ARIA OrchestratorAgent
----------------------
Central coordinator. Parses the biological question, designs the
analysis plan, and dispatches specialized agents in order.

Checkpoints:
  1 -> Data audit confirmation
  2 -> Analysis plan approval
  3 -> QC / parameter decisions (raised by sub-agents)
  4 -> Preliminary findings review
  5 -> Final report approval
"""

from __future__ import annotations
import uuid
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Message, MessageType, Confidence, CavemanMode, bus
from aria.memory.memory import ARIAMemory


ORCHESTRATOR_SYSTEM = """
You are ARIA's Orchestrator — the strategic coordinator of a bioinformatics agent system.

Your responsibilities:
1. Understand the biological question behind the user's request
2. Design an analysis plan that uses available data modalities optimally
3. Identify dependencies between analyses (what must run before what)
4. Assess whether findings across modalities converge or conflict
5. Decide when integration analysis adds real biological value

Available agents:
- data_audit_agent:  always runs first
- rna_agent:         bulk RNA-seq and scRNA-seq (including spatial)
- chromatin_agent:   ATAC, ChIP, CUT&RUN, CUT&TAG
- genome_arch_agent: HiC, TADs, loops, compartments
- integration_agent: multimodal integration (requires >= 2 modalities)
- narrative_agent:   final report and visualizations

Always think about the biology first, then the methods.
""".strip()


class OrchestratorAgent(BaseAgent):

    name        = "orchestrator"
    description = "Central coordinator — routes tasks to specialized agents."

    MODALITY_TO_AGENT = {
        "scRNA":       "rna_agent",
        "bulk_RNA":    "rna_agent",
        "spatial":     "rna_agent",
        "scATAC":      "chromatin_agent",
        "bulk_ATAC":   "chromatin_agent",
        "ChIP":        "chromatin_agent",
        "CUT_AND_RUN": "chromatin_agent",
        "CUT_AND_TAG": "chromatin_agent",
        "HiC":         "genome_arch_agent",
    }

    def __init__(self, memory: ARIAMemory, llm=None, api_key: str = None):
        super().__init__(memory, llm, api_key)
        self._pending_checkpoints: dict[str, Message] = {}
        self._experiment_plans:    dict[str, dict]    = {}

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
        plan        = self._experiment_plans.get(experiment_id, {})
        audit_agent = DataAuditAgent(self.memory, self.llm)
        return audit_agent.run(experiment_id, plan["context"])

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
            return self._after_checkpoint_1(experiment_id, user_decision, resolved_msg)
        elif cp == 2:
            return self._after_checkpoint_2(experiment_id, user_decision, resolved_msg)
        return {"status": "ok"}

    def _after_checkpoint_1(self, experiment_id, decision, msg) -> dict:
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

    def _after_checkpoint_2(self, experiment_id, decision, msg) -> dict:
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
        )
        return {"status": "analysis_running", "plan": plan,
                "exp_context": exp_context}

    def _parse_question(self, user_question: str) -> dict:
        prompt = f"""
User question: "{user_question}"

Extract biological analysis intent. Return JSON:
{{
  "analysis_type": "differential|regulatory|structural|temporal|cell_type|integration",
  "biological_entities": ["genes, TFs, cell types, conditions mentioned"],
  "comparison": "what vs what (if applicable)",
  "key_modalities_needed": ["RNA","ATAC","HiC","ChIP","CUT_AND_RUN","CUT_AND_TAG","spatial"],
  "complexity": "simple|moderate|complex",
  "summary": "one sentence biological question summary"
}}
"""
        try:
            return self.think_structured(prompt, system=ORCHESTRATOR_SYSTEM)
        except Exception:
            return {
                "analysis_type": "differential",
                "summary": user_question,
                "complexity": "moderate",
            }

    def _design_analysis_plan(self, experiment_id: str,
                               exp_context: dict) -> dict:
        modalities    = list(exp_context.get("modalities", {}).keys())
        user_question = exp_context.get("user_question", "")
        prompt = f"""
Available modalities: {modalities}
User question: "{user_question}"
Organism: {exp_context.get("organism")}
Genome:   {exp_context.get("genome")}
Multimodal: {exp_context.get("is_multimodal")}

Design the optimal analysis pipeline. Return JSON:
{{
  "steps": [
    {{"order":1,"agent":"rna_agent","analysis":"description",
      "depends_on":[],"can_parallel":false}}
  ],
  "integration_needed": true,
  "integration_type": "WNN|MOFA|ArchR|scglue|none",
  "estimated_complexity": "low|medium|high",
  "rationale": "brief explanation"
}}
"""
        try:
            plan = self.think_structured(prompt, system=ORCHESTRATOR_SYSTEM)
        except Exception:
            plan = {
                "steps": [{"order": 1, "agent": "rna_agent",
                            "analysis": "scRNA-seq analysis",
                            "depends_on": [], "can_parallel": False}],
                "integration_needed": False,
                "integration_type": "none",
                "estimated_complexity": "medium",
                "rationale": "Default single-modality plan",
            }
        self._experiment_plans[experiment_id]["plan"] = plan
        return plan

    def _format_plan_summary(self, plan: dict) -> str:
        steps = plan.get("steps", [])
        lines = ["ARIA Analysis Plan:\n"]
        for step in steps:
            parallel = " (parallel)" if step.get("can_parallel") else ""
            lines.append(
                f"  Step {step['order']}: [{step['agent']}] "
                f"{step['analysis']}{parallel}"
            )
        if plan.get("integration_needed"):
            lines.append(
                f"\n  Integration: {plan.get('integration_type','TBD')}"
            )
        lines.append(f"\n  Complexity: {plan.get('estimated_complexity','?')}")
        lines.append(f"  Rationale:  {plan.get('rationale','')}")
        lines.append("\nProceed with this plan?")
        return "\n".join(lines)

    def receive(self, message: Message):
        if message.type == MessageType.FINDING:
            exp_id = message.experiment_id
            if exp_id in self._experiment_plans:
                self._experiment_plans[exp_id]["findings"].append(
                    message.payload
                )
        elif message.type == MessageType.ESCALATION:
            self._pending_checkpoints[message.id] = message
