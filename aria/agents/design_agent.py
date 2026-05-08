"""
ARIA DesignAgent (v4.0)
------------------------
Interactive experimental design agent for ARIA.

Runs AFTER DataAuditAgent (Checkpoint 1) and BEFORE the analysis plan.
Converts raw file detection into a formal, user‑confirmed experimental design
suitable for differential expression (DESeq2).

State machine design:
  - STEP_GROUPS → STEP_ORGANISM → STEP_FACTOR → STEP_BATCH →
    STEP_PSEUDOREP → STEP_CONFIRM → DONE
  - Each step publishes exactly one escalation (checkpoint) and returns.
  - The orchestrator calls handle_user_response() when the user answers.
  - When DONE, the design is stored in memory.
"""

from __future__ import annotations

import json
import logging
import re
import uuid
from enum import Enum, auto
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence, bus
from aria.llm.provider import LLMProvider, TaskTier
from aria.memory.memory import ARIAMemory

log = logging.getLogger("aria.design")


DESIGN_SYSTEM = """
You are ARIA's DesignAgent — an expert in RNA‑seq experimental design.

Your job is to help the user formalise their experiment:
  - Identify biological groups from file names
  - Propose a short, descriptive group name
  - Suggest the most likely experimental factor (genotype, treatment, time…)

Rules:
  - Replicate identifiers (rep1, rep2, R1, R2, _1, _2) are NOT groups.
  - Group names should be simple: "WT", "KO", "Treated", "Control".
  - Never invent groups that are not supported by the file names.
  - If unsure, ask the user.
""".strip()


class DesignStep(Enum):
    START          = auto()
    GROUPS         = auto()
    ORGANISM       = auto()
    FACTOR         = auto()
    BATCH          = auto()
    PSEUDOREP      = auto()
    CONFIRM        = auto()
    DONE           = auto()


class DesignAgent(BaseAgent):
    name        = "design_agent"
    description = "Interactive experimental design — groups, factor, batch, organism."

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm, api_key)
        # State machine
        self._step: DesignStep = DesignStep.START
        # Data accumulated across steps
        self._parsed_samples: list[dict] = []
        self._proposed_groups: dict = {}
        self._confirmed_groups: Optional[dict] = None
        self._organism: Optional[str] = None
        self._genome: Optional[str] = None
        self._main_factor: Optional[str] = None
        self._batch_covariate: Optional[str] = None
        # IDs for current escalation (needed to resolve later)
        self._pending_escalation_id: Optional[str] = None
        self._experiment_id: Optional[str] = None

    # ── Public entry: called by orchestrator after CP1 ──────────────────

    def start_design(self, experiment_id: str,
                     exp_context: dict,
                     biological_intent: dict) -> dict:
        """
        Kick off the design process. Publishes the first checkpoint (groups)
        and returns immediately. The orchestrator must call
        handle_user_response() for each subsequent step.
        """
        self._experiment_id = experiment_id
        self._step = DesignStep.START

        # Collect sample filenames from the context
        raw_files = []
        for mod in ("bulk_RNA", "bulk_RNA_raw", "scRNA"):
            raw_files.extend(exp_context.get("modalities", {}).get(mod, []))
        if not raw_files:
            self.publish_finding(
                experiment_id,
                {"summary": "No sample files found for design phase."},
                Confidence.INSUFFICIENT,
            )
            self._step = DesignStep.DONE
            return {"status": "failed", "reason": "no_samples"}

        # Store parsed samples and intent for later steps
        self._parsed_samples = self._parse_samples(raw_files)
        self._biological_intent = biological_intent

        # Ask LLM to propose groups (not blocking)
        self._proposed_groups = self._propose_groups(self._parsed_samples,
                                                      biological_intent)

        # Publish first checkpoint: group confirmation
        self._step = DesignStep.GROUPS
        self._publish_groups_checkpoint()

        return {"status": "awaiting_user", "step": "groups"}

    def handle_user_response(self, experiment_id: str,
                              checkpoint_num: float,
                              choice: str) -> dict:
        """
        Process the user's answer for the current step and advance the
        state machine. Publishes the next checkpoint or finalises.
        Returns a dict with at least {"status": "awaiting_user" | "done" | "cancelled"}.
        """
        if self._step == DesignStep.GROUPS:
            return self._handle_groups_response(choice)
        elif self._step == DesignStep.ORGANISM:
            return self._handle_organism_response(choice)
        elif self._step == DesignStep.FACTOR:
            return self._handle_factor_response(choice)
        elif self._step == DesignStep.BATCH:
            return self._handle_batch_response(choice)
        elif self._step == DesignStep.PSEUDOREP:
            return self._handle_pseudorep_response(choice)
        elif self._step == DesignStep.CONFIRM:
            return self._handle_confirm_response(choice)
        else:
            return {"status": "error", "reason": f"Unknown step {self._step}"}

    # ── Step handlers ───────────────────────────────────────────────────

    def _handle_groups_response(self, choice: str) -> dict:
        if choice == "Yes — confirm groups":
            self._confirmed_groups = self._proposed_groups.get("groups", {})
        elif choice == "Re‑infer groups":
            self._proposed_groups = self._propose_groups(
                self._parsed_samples, self._biological_intent
            )
            self._publish_groups_checkpoint()
            return {"status": "awaiting_user", "step": "groups"}
        else:
            # Any other choice: keep original proposal for now
            self._confirmed_groups = self._proposed_groups.get("groups", {})

        # Move to organism step
        self._step = DesignStep.ORGANISM
        self._publish_organism_checkpoint()
        return {"status": "awaiting_user", "step": "organism"}

    def _handle_organism_response(self, choice: str) -> dict:
        import re as _re
        org_map = {
            "Homo sapiens (hg38)":      ("Homo sapiens",      "hg38"),
            "Mus musculus (mm39)":      ("Mus musculus",      "mm39"),
            "Rattus norvegicus (rn7)":  ("Rattus norvegicus", "rn7"),
            "Danio rerio (danRer11)":   ("Danio rerio",       "danRer11"),
        }
        if choice in org_map:
            self._organism, self._genome = org_map[choice]
        else:
            # Free-text: expect "Species name (assembly)" — e.g. "Gallus gallus (galGal6)"
            m = _re.match(r'^(.+?)\s*\(([^)]+)\)\s*$', choice.strip())
            if m:
                self._organism = m.group(1).strip()
                self._genome   = m.group(2).strip()
            elif choice.strip() and not choice.lower().startswith("other"):
                # Plain organism name without assembly
                self._organism = choice.strip()
                self._genome   = "unknown"
            else:
                # Empty input or unrecognised "Other..." string — warn and fall back
                log.warning(
                    "Organism free-text was empty or not parseable ('%s'). "
                    "Defaulting to Homo sapiens / hg38.", choice
                )
                self._organism, self._genome = "Homo sapiens", "hg38"

        self._step = DesignStep.FACTOR
        self._publish_factor_checkpoint()
        return {"status": "awaiting_user", "step": "factor"}

    def _handle_factor_response(self, choice: str) -> dict:
        factor_map = {
            "genotype": "genotype",
            "treatment": "treatment",
            "condition": "condition",
            "time": "time",
        }
        self._main_factor = factor_map.get(choice, "condition")
        self._step = DesignStep.BATCH
        self._publish_batch_checkpoint()
        return {"status": "awaiting_user", "step": "batch"}

    def _handle_batch_response(self, choice: str) -> dict:
        if choice == "Yes — include as covariate":
            # Pick the first detected batch keyword as name
            batch_keywords = {"lane", "batch", "plate", "run", "flowcell"}
            for sample in self._parsed_samples:
                for tok in sample["tokens"]:
                    if tok in batch_keywords:
                        self._batch_covariate = tok
                        break
                if self._batch_covariate:
                    break
            if not self._batch_covariate:
                self._batch_covariate = "batch"
        else:
            self._batch_covariate = None

        self._step = DesignStep.PSEUDOREP
        self._publish_pseudorep_checkpoint()
        return {"status": "awaiting_user", "step": "pseudorep"}

    def _handle_pseudorep_response(self, choice: str) -> dict:
        # Just log the warning; no change to design
        self._step = DesignStep.CONFIRM
        self._publish_confirm_checkpoint()
        return {"status": "awaiting_user", "step": "confirm"}

    def _handle_confirm_response(self, choice: str) -> dict:
        if choice != "Yes — proceed":
            return {"status": "cancelled", "reason": "User did not approve design"}

        # Build final design
        design = self._build_design()
        self._store_design(design)
        self._step = DesignStep.DONE
        return {"status": "done", "design": design}

    # ── Checkpoint publishers ────────────────────────────────────────────

    def _publish_groups_checkpoint(self):
        groups = self._proposed_groups.get("groups", {})
        summary_lines = []
        for g, members in groups.items():
            summary_lines.append(
                f"  {g}: {len(members)} sample(s) — {', '.join(members[:3])}{'...' if len(members)>3 else ''}"
            )
        question = (
            "ARIA proposes the following biological groups:\n\n" +
            "\n".join(summary_lines) +
            f"\n\nConfidence: {self._proposed_groups.get('confidence', '?')}\n"
            f"Reasoning: {self._proposed_groups.get('reasoning', 'none')}\n\n"
            "Is this correct?"
        )
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.1,
            question=question,
            options=[
                "Yes — confirm groups",
                "Re‑infer groups",
                "Edit group names (using defaults)",
                "Cancel experiment",
            ],
            context={"proposed_groups": groups},
        )
        self._pending_escalation_id = msg.id

    def _publish_organism_checkpoint(self):
        question = (
            "Please confirm the organism and genome assembly:\n\n"
            "[1] Homo sapiens (hg38)\n"
            "[2] Mus musculus (mm39)\n"
            "[3] Rattus norvegicus (rn7)\n"
            "[4] Danio rerio (danRer11)\n"
            "[5] Other — type name and assembly (e.g. Gallus gallus (galGal6))\n"
        )
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.2,
            question=question,
            options=[
                "Homo sapiens (hg38)",
                "Mus musculus (mm39)",
                "Rattus norvegicus (rn7)",
                "Danio rerio (danRer11)",
                "Other organism / assembly...",
            ],
        )
        self._pending_escalation_id = msg.id

    def _publish_factor_checkpoint(self):
        group_names = list(self._confirmed_groups.keys())
        # Quick LLM suggestion (non‑blocking, optional)
        factor_guess = "condition"
        try:
            prompt = (
                f"Groups: {', '.join(group_names)}\n"
                "Most likely experimental factor? One word."
            )
            raw = self.think(prompt, system=DESIGN_SYSTEM,
                             tier=TaskTier.LIGHT, max_tokens=15)
            factor_guess = raw.strip().lower()
        except Exception:
            pass

        question = (
            f"Main experimental factor: proposed '{factor_guess}'\n\n"
            "Choose the factor that is explicitly varied in this experiment:"
        )
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.3,
            question=question,
            options=[
                "genotype",
                "treatment",
                "condition",
                "time",
                "Cancel",
            ],
        )
        self._pending_escalation_id = msg.id

    def _publish_batch_checkpoint(self):
        # Detect batch keywords
        batch_keywords = {"lane", "batch", "plate", "run", "flowcell"}
        possible_batch = set()
        for p in self._parsed_samples:
            for tok in p["tokens"]:
                if tok in batch_keywords or re.match(r'^(lane|batch)\d+$', tok):
                    possible_batch.add(tok)

        if not possible_batch:
            # No batch detected → skip to pseudorep
            self._batch_covariate = None
            self._step = DesignStep.PSEUDOREP
            self._publish_pseudorep_checkpoint()
            return

        question = (
            f"Possible batch variable detected: {', '.join(sorted(possible_batch))}\n\n"
            "Do you want to include this as a covariate in the design formula?"
        )
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.4,
            question=question,
            options=[
                "Yes — include as covariate",
                "No — ignore batch",
                "Cancel",
            ],
        )
        self._pending_escalation_id = msg.id

    def _publish_pseudorep_checkpoint(self):
        # Heuristic pseudoreplication check
        def _root(stem: str) -> str:
            return re.sub(r'[_\-]?(rep\d*|r\d*|_\d+)$', '', stem, flags=re.IGNORECASE)

        warning_groups = []
        for group, members in (self._confirmed_groups or {}).items():
            roots = set(_root(m) for m in members)
            if len(roots) < len(members):
                warning_groups.append(group)

        if not warning_groups:
            self._step = DesignStep.CONFIRM
            self._publish_confirm_checkpoint()
            return

        question = (
            f"Group(s) {', '.join(warning_groups)} appear to have samples derived "
            f"from the same biological source (e.g., technical replicates).\n\n"
            "Are these truly independent biological replicates?"
        )
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.5,
            question=question,
            options=[
                "Yes — they are independent biological replicates",
                "No — they are technical replicates (consider pairing)",
                "Not sure — proceed anyway",
            ],
        )
        self._pending_escalation_id = msg.id

    def _publish_confirm_checkpoint(self):
        design = self._build_design()
        summary = [
            f"Organism: {design['organism']}",
            f"Genome:   {design['genome']}",
            f"Groups:   {', '.join(design['groups'].keys())}",
            f"Factor:   {design['main_factor']}",
            f"Design formula: {design['design_formula']}",
            f"Replicates: {design['replicates']}",
        ]
        if design.get("batch_covariate"):
            summary.append(f"Batch covariate: {design['batch_covariate']}")

        question = "Experimental design ready:\n\n" + "\n".join(summary) + "\n\nProceed with analysis?"
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.6,
            question=question,
            options=[
                "Yes — proceed",
                "Cancel",
            ],
        )
        self._pending_escalation_id = msg.id

    # ── Helpers ──────────────────────────────────────────────────────────

    @staticmethod
    def _parse_samples(file_paths: list[str]) -> list[dict]:
        """Parse sample names, merging paired-end read pairs (R1/R2)."""
        # Agrupar por nombre base (quitando _1/_2 o _R1/_R2)
        pairs = {}
        for f in file_paths:
            stem = Path(f).stem
            while stem.endswith(('.fastq', '.fq', '.gz', '.tsv', '.csv', '.txt')):
                stem = Path(stem).stem
            # Normalizar nombres paired-end: quitar sufijo _1/_2 o _R1/_R2
            base = re.sub(r'(_[12]|_R[12])$', '', stem, flags=re.IGNORECASE)
            pairs.setdefault(base, []).append(f)

        parsed = []
        for base, files in pairs.items():
            tokens = re.split(r'[_\-\.\s]+', base.lower())
            tokens = [t for t in tokens if t]
            parsed.append({
                "raw": files[0],         # usar el primer archivo del par
                "stem": base,
                "tokens": tokens,
                "paired_files": files    # todos los archivos del par
            })
        return parsed

    def _propose_groups(self, parsed: list[dict], intent: dict) -> dict:
        stems = [p["stem"] for p in parsed]
        prompt = (
            "You are analyzing RNA‑seq sample names.\n\n"
            "Samples:\n" +
            "\n".join(f"- {s}" for s in stems) +
            "\n\nTask:\n"
            "1. Infer biological groups (NOT technical replicates).\n"
            "2. Assign each sample to a group.\n"
            "3. Return ONLY valid JSON.\n\n"
            'Output format:\n'
            '{\n  "groups": {\n    "group_name": ["sample1", "sample2"]\n  },\n'
            '  "confidence": "high|medium|low",\n'
            '  "reasoning": "short explanation"\n}'
        )
        try:
            result = self.think_structured(prompt, system=DESIGN_SYSTEM,
                                           schema_hint="JSON with groups, confidence, reasoning.",
                                           tier=TaskTier.MEDIUM)
            if isinstance(result, dict) and "groups" in result:
                return result
        except Exception as e:
            log.warning(f"LLM group proposal failed: {e}")
        # Heuristic fallback
        groups = {}
        for p in parsed:
            first = p["tokens"][0] if p["tokens"] else "group"
            groups.setdefault(first, []).append(p["stem"])
        return {"groups": groups, "confidence": "low",
                "reasoning": "Heuristic fallback (LLM unavailable)"}

    def _build_design(self) -> dict:
        if self._batch_covariate:
            formula = f"~ {self._batch_covariate} + {self._main_factor}"
        else:
            formula = f"~ {self._main_factor}"

        replicates = {g: len(mems) for g, mems in (self._confirmed_groups or {}).items()}

        return {
            "organism":        self._organism,
            "genome":          self._genome,
            "groups":          self._confirmed_groups,
            "main_factor":     self._main_factor,
            "design_formula":  formula,
            "batch_covariate": self._batch_covariate,
            "replicates":      replicates,
            "n_total_samples": sum(replicates.values()),
        }

    def _store_design(self, design: dict):
        try:
            self.memory.store_decision(
                decision_id=str(uuid.uuid4())[:8],
                wing_id=self._experiment_id,
                checkpoint=2,
                question="Experimental design (organism, groups, factor, formula)",
                decision=json.dumps(design, indent=2),
                rationale="User confirmed design interactively",
                made_by="user",
            )
            log.info(f"Design stored for {self._experiment_id}")
        except Exception as e:
            log.warning(f"Failed to store design decision: {e}")

    # ── Satisfy BaseAgent's abstract method ──────────────────────────────
    def run(self, experiment_id: str, context: dict) -> dict:
        """
        Not used in state-machine mode; use start_design() and
        handle_user_response() instead.
        """
        raise NotImplementedError(
            "DesignAgent uses start_design() + handle_user_response(), not run()"
        )

    def receive(self, message):
        pass
