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

        geo_metadata    = exp_context.get("geo_metadata")
        inferred_design = exp_context.get("inferred_design", {})

        if geo_metadata and inferred_design.get("groups"):
            # GEO dataset: groups, organism, and genome are already known from
            # SOFT metadata — always prefer them over filename-based inference.
            geo_samples = [
                s
                for samples in inferred_design["groups"].values()
                for s in samples
            ]
            self._parsed_samples = [
                {"raw": s, "stem": s,
                 "tokens": s.lower().split(), "paired_files": [s]}
                for s in geo_samples
            ]
            self._proposed_groups = {
                "groups":     inferred_design["groups"],
                "confidence": "high",
                "reasoning":  "Groups inferred from GEO SOFT metadata",
            }
            # Pre-seed organism/genome so CP2.2 shows the right value
            self._organism = inferred_design.get("organism", "")
            self._genome   = inferred_design.get("genome", "")
        elif not raw_files:
            self.publish_finding(
                experiment_id,
                {"summary": "No sample files found for design phase."},
                Confidence.INSUFFICIENT,
            )
            self._step = DesignStep.DONE
            return {"status": "failed", "reason": "no_samples"}
        else:
            # Normal file-based path: infer groups from sample filenames via LLM
            self._parsed_samples = self._parse_samples(raw_files)
            self._proposed_groups = self._propose_groups(self._parsed_samples,
                                                          biological_intent)

        self._biological_intent = biological_intent

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
        elif choice == "Cancel experiment":
            return {"status": "cancelled", "reason": "user_cancelled_groups"}
        else:
            # Free-text manual assignment, fed back via the TUI's "Other —"
            # path. Accepted formats:
            #   "young=hc1153,hc1265; old=hc11,hc77"
            #   "young: hc1153 hc1265 | old: hc11 hc77"
            # Sample matching is substring-based against parsed stems so the
            # user does not have to type the full 10x suffix.
            parsed = self._parse_manual_groups(choice)
            if parsed:
                self._proposed_groups = {
                    "groups":     parsed,
                    "confidence": "high",
                    "reasoning":  "User-provided manual group assignment",
                }
                self._confirmed_groups = parsed
            else:
                log.warning(
                    "Could not parse manual group assignment '%s'. "
                    "Keeping the prior proposal.", choice
                )
                self._confirmed_groups = self._proposed_groups.get("groups", {})

        # Move to organism step
        self._step = DesignStep.ORGANISM
        self._publish_organism_checkpoint()
        return {"status": "awaiting_user", "step": "organism"}

    def _parse_manual_groups(self, text: str) -> dict:
        """
        Parse a free-text manual group assignment from the user. Supports
        common shapes:
          - "young=hc1153,hc1265; old=hc11,hc77"
          - "young: hc1153 hc1265 | old: hc11 hc77"
          - JSON: {"young": ["hc1153","hc1265"], "old": ["hc11","hc77"]}
        Returns {} if nothing usable could be parsed (caller decides what to
        do with that — we never silently invent groups).
        """
        if not text or not isinstance(text, str):
            return {}
        stripped = text.strip()
        # JSON shortcut
        if stripped.startswith("{") and stripped.endswith("}"):
            try:
                obj = json.loads(stripped)
                if isinstance(obj, dict) and obj:
                    return self._normalise_manual_groups(obj)
            except (json.JSONDecodeError, ValueError):
                pass
        # Generic "g1=s1,s2; g2=s3,s4" / "g1: s1 s2 | g2: s3 s4"
        pieces = re.split(r"[;|]+", stripped)
        proposed: dict = {}
        for chunk in pieces:
            if not chunk or ("=" not in chunk and ":" not in chunk):
                continue
            sep = "=" if "=" in chunk else ":"
            name, body = chunk.split(sep, 1)
            name = name.strip()
            if not name:
                continue
            members = [m for m in re.split(r"[,\s]+", body) if m]
            if members:
                proposed[name] = members
        if not proposed:
            return {}
        return self._normalise_manual_groups(proposed)

    def _normalise_manual_groups(self, proposed: dict) -> dict:
        """
        Map user-supplied sample tokens (e.g. 'hc1153') back to the canonical
        parsed stems (e.g. 'GSE278576_hc1153_raw_feature_bc_matrix') using
        substring match. Drops tokens that match no parsed sample.
        """
        if not self._parsed_samples:
            # No parsed samples to bind against — accept the user's strings
            # verbatim so manual workflows are not blocked.
            return {str(g): [str(s) for s in mems]
                    for g, mems in proposed.items() if mems}
        stems = [p["stem"] for p in self._parsed_samples]
        out: dict = {}
        for grp, tokens in proposed.items():
            if not isinstance(tokens, (list, tuple)):
                tokens = [tokens]
            resolved: list[str] = []
            for tok in tokens:
                tok_str = str(tok).strip()
                if not tok_str:
                    continue
                # Exact match wins
                if tok_str in stems:
                    resolved.append(tok_str)
                    continue
                # Substring match — first hit wins (preserves user order)
                hit = next((s for s in stems if tok_str in s), None)
                if hit:
                    resolved.append(hit)
                else:
                    log.warning(
                        "Manual group assignment: token '%s' did not match "
                        "any sample stem — dropping.", tok_str
                    )
            if resolved:
                out[str(grp)] = resolved
        return out

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
            stripped = choice.strip()
            # Free-text: expect "Species name (assembly)" — e.g. "Gallus gallus (galGal6)"
            m = _re.match(r'^(.+?)\s*\(([^)]+)\)\s*$', stripped)
            if m:
                self._organism = m.group(1).strip()
                self._genome   = m.group(2).strip()
            elif stripped and not stripped.lower().startswith("other"):
                # Plain organism name without assembly
                self._organism = stripped
                self._genome   = "unknown"
            else:
                # Empty input or literal "Other..." sentinel — don't silently
                # default. Re-publish the organism checkpoint so the user gets
                # another chance instead of an inferred hg38 they never picked.
                log.warning(
                    "Organism free-text was empty or unparseable ('%s'). "
                    "Re-asking the user.", choice
                )
                self._publish_organism_checkpoint()
                return {"status": "awaiting_user", "step": "organism"}

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
            "Is this correct?\n\n"
            "If filenames carry no comparison signal (donor barcodes etc.) "
            "pick the manual option below and type the mapping like:\n"
            "  young=hc1153,hc1265; old=hc11,hc77"
        )
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.1,
            question=question,
            options=[
                "Yes — confirm groups",
                "Re‑infer groups",
                "Other — manually assign groups (e.g. young=hc1153,hc1265; old=hc11,hc77)",
                "Cancel experiment",
            ],
            context={"proposed_groups": groups},
        )
        self._pending_escalation_id = msg.id

    def _publish_organism_checkpoint(self):
        std_options = [
            "Homo sapiens (hg38)",
            "Mus musculus (mm39)",
            "Rattus norvegicus (rn7)",
            "Danio rerio (danRer11)",
            "Other organism / assembly...",
        ]

        # If organism is already known (e.g. from GEO), surface it as first option
        if self._organism and self._genome:
            geo_opt = f"{self._organism} ({self._genome})"
            if geo_opt in std_options:
                options = std_options
                geo_note = ""
            else:
                options = [geo_opt] + [o for o in std_options
                                        if o != "Other organism / assembly..."] \
                          + ["Other organism / assembly..."]
                geo_note = f"\n\n  [GEO] Detected: {geo_opt}"
        else:
            options     = std_options
            geo_note    = ""

        lines = [f"[{i}] {o}" for i, o in enumerate(options, 1)]
        question = (
            "Please confirm the organism and genome assembly:"
            + geo_note + "\n\n"
            + "\n".join(lines) + "\n"
        )
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.2,
            question=question,
            options=options,
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
        # The user's biological question often names the comparison (e.g.
        # "young vs old", "WT vs KO"). Surface it to the LLM — without it the
        # model can only pattern-match on filenames, which fails when sample
        # IDs are opaque (donor barcodes hc11/hc77/...).
        intent_block = ""
        question = (intent or {}).get("question") or (intent or {}).get("text")
        comparison = (intent or {}).get("comparison")
        if question or comparison:
            intent_block = (
                "\n\nUser's biological question (use it to name the groups, "
                "but only assign a sample to a group when the filename or "
                "external metadata supports it — do NOT guess):\n"
                f"  {question or comparison}\n"
            )
        prompt = (
            "You are analyzing RNA‑seq sample names.\n\n"
            "Samples:\n" +
            "\n".join(f"- {s}" for s in stems) +
            intent_block +
            "\n\nTask:\n"
            "1. Infer biological groups (NOT technical replicates).\n"
            "2. Assign each sample to a group.\n"
            "3. If the filenames carry no signal that maps samples to the\n"
            "   comparison the user asked about (e.g. donor barcodes with no\n"
            "   age/condition encoding), return a SINGLE group named 'all'\n"
            "   containing every sample and set confidence='low' — the user\n"
            "   will be offered a manual-assignment option downstream.\n"
            "4. Return ONLY valid JSON.\n\n"
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
