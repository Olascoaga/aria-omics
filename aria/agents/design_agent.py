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
        self._inferred_design: dict = {}
        self._pseudobulk_design: dict = {}
        self._replicate_handling: Optional[dict] = None
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
        self._replicate_handling = None

        # Collect sample filenames from the context
        raw_files = []
        for mod in ("bulk_RNA", "bulk_RNA_raw", "scRNA"):
            raw_files.extend(exp_context.get("modalities", {}).get(mod, []))

        geo_metadata    = exp_context.get("geo_metadata")
        inferred_design = exp_context.get("inferred_design", {})
        self._inferred_design = inferred_design or {}

        if inferred_design.get("groups"):
            # External metadata or h5ad obs already carries the experimental
            # design; prefer that over filename-based inference.
            inferred_samples = [
                s
                for samples in inferred_design["groups"].values()
                for s in samples
            ]
            self._parsed_samples = [
                {"raw": s, "stem": s,
                 "tokens": s.lower().split(), "paired_files": [s]}
                for s in inferred_samples
            ]
            source = inferred_design.get(
                "source",
                "GEO SOFT metadata" if geo_metadata else "inferred metadata",
            )
            self._proposed_groups = {
                "groups":     inferred_design["groups"],
                "confidence": inferred_design.get("confidence", "high"),
                "reasoning":  inferred_design.get(
                    "reasoning", f"Groups inferred from {source}"
                ),
            }
            # Pre-seed known values so checkpoints surface metadata from the
            # dataset instead of asking the LLM to infer them from sample IDs.
            self._organism = inferred_design.get("organism", "") or exp_context.get("organism")
            self._genome   = inferred_design.get("genome", "") or exp_context.get("genome")
            self._main_factor = (
                inferred_design.get("main_factor")
                or inferred_design.get("condition_col")
            )
            self._batch_covariate = inferred_design.get("batch_covariate")
            self._pseudobulk_design = inferred_design.get("pseudobulk", {}) or {}
        elif not raw_files:
            # No RNA samples and no inferred DE groups. A chromatin-only run
            # (scATAC / bulk_ATAC / ChIP / CUT&RUN / CUT&TAG) needs no
            # differential-expression design phase: skip it honestly with a
            # minimal design so the orchestrator can dispatch the modality
            # agent directly. This keeps the DE design path byte-identical for
            # RNA runs while unblocking the live scATAC orchestrator dispatch.
            modalities = exp_context.get("modalities", {}) or {}
            rna_modalities = [
                m for m in modalities
                if m in ("bulk_RNA", "bulk_RNA_raw", "scRNA")
            ]
            if modalities and not rna_modalities:
                self._step = DesignStep.DONE
                minimal_design = {
                    "organism": exp_context.get("organism") or "unknown",
                    "genome":   exp_context.get("genome") or "unknown",
                    "groups": {},
                    "main_factor": None,
                    "design_formula": None,
                    "n_total_samples": 0,
                    "pseudobulk": {},
                    "design_type": "no_de_design",
                    "modalities_without_de": sorted(modalities.keys()),
                }
                self.publish_finding(
                    experiment_id,
                    {"summary": (
                        "No differential-expression design needed for "
                        f"modality(ies) {sorted(modalities.keys())}; "
                        "proceeding directly to modality dispatch."
                    )},
                    Confidence.HIGH,
                )
                return {"status": "skipped",
                        "reason": "no_de_design_needed",
                        "design": minimal_design}
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
            #   "treated=s1,s2,s3; control=s4,s5,s6"
            #   "treated: s1 s2 s3 | control: s4 s5 s6"
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
          - "treated=s1,s2,s3; control=s4,s5,s6"
          - "treated: s1 s2 s3 | control: s4 s5 s6"
          - JSON: {"treated": ["s1","s2","s3"], "control": ["s4","s5","s6"]}
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
        Map user-supplied sample tokens (e.g. 'sample01') back to the canonical
        parsed stems (e.g. 'sample01_raw_feature_bc_matrix') using
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
        if choice == "Cancel":
            return {"status": "cancelled", "reason": "user_cancelled_factor"}
        self._main_factor = factor_map.get(choice, choice.strip() or "condition")
        self._step = DesignStep.BATCH
        self._publish_batch_checkpoint()
        return {"status": "awaiting_user", "step": "batch"}

    def _handle_batch_response(self, choice: str) -> dict:
        if choice == "Yes — include as covariate":
            # Preserve metadata-seeded batch/covariate when present. Otherwise
            # pick the first detected batch keyword as name.
            batch_keywords = {"lane", "batch", "plate", "run", "flowcell"}
            if not self._batch_covariate:
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
        choice_l = str(choice or "").lower()
        if choice_l.startswith("yes") and "independent" in choice_l:
            self._replicate_handling = self._build_replicate_handling(
                technical=False
            )
        elif choice_l.startswith("no") and "technical" in choice_l:
            self._replicate_handling = self._build_replicate_handling(
                technical=True
            )
        else:
            # B1: uncertainty cannot authorize independent biological n.
            return {
                "status": "cancelled",
                "reason": "replicate_structure_unresolved",
            }
        self._step = DesignStep.CONFIRM
        self._publish_confirm_checkpoint()
        return {"status": "awaiting_user", "step": "confirm"}

    def _handle_confirm_response(self, choice: str) -> dict:
        design = self._build_design()
        if design.get("design_status") == "blocking":
            if design.get("confounded_covariates"):
                reason = "condition_covariate_confounding"
            elif design.get("replicate_handling", {}).get("conflicts"):
                reason = "replicate_mapping_conflict"
            else:
                reason = "insufficient_biological_replicates"
            return {
                "status": "cancelled",
                "reason": reason,
                "design": design,
            }
        if choice != "Yes — proceed":
            return {"status": "cancelled", "reason": "User did not approve design"}

        # Build final design
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
            "If filenames carry no comparison signal (donor barcodes, "
            "single-file h5ad without condition encoding, etc.) pick the "
            "manual option below and type the mapping as "
            "<group>=<sample1>,<sample2>;<group2>=<sample3>,<sample4>. "
            "For example, two conditions with three replicates each:\n"
            "  treated=s1,s2,s3; control=s4,s5,s6"
        )
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.1,
            question=question,
            options=[
                "Yes — confirm groups",
                "Re‑infer groups",
                "Other — manually assign groups (e.g. treated=s1,s2,s3; control=s4,s5,s6)",
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
        factor_guess = self._main_factor or "condition"
        if not self._main_factor:
            try:
                prompt = (
                    f"Groups: {', '.join(group_names)}\n"
                    "Most likely experimental factor? One word."
                )
                raw = self.think(prompt, system=DESIGN_SYSTEM,
                                 tier=TaskTier.LIGHT, max_tokens=15)
                factor_guess = self._sanitize_factor_guess(raw)
            except Exception:
                factor_guess = "condition"

        question = (
            f"Main experimental factor: proposed '{factor_guess}'\n\n"
            "Choose the factor that is explicitly varied in this experiment:"
        )
        options = [
            "genotype",
            "treatment",
            "condition",
            "stimulation",
            "time",
            "Cancel",
        ]
        if factor_guess and factor_guess not in options:
            options = [factor_guess] + options
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.3,
            question=question,
            options=options,
        )
        self._pending_escalation_id = msg.id

    @staticmethod
    def _sanitize_factor_guess(raw: str) -> str:
        """
        Reduce a free-text LLM reply to a single short identifier suitable for
        a factor name. Rejects anything that looks like a sentence, contains
        whitespace, punctuation, or runs over a sensible length cap. Returns
        'condition' as a neutral fallback when nothing safe can be salvaged.

        The bug this guards against: when obs inference fails and the LLM
        free-texts something like "need file names to identify groups" or
        "paste your matrix file", the previous code passed the whole reply
        as a CP 2.3 option. Now the menu only ever sees a column-like token.
        """
        if not raw or not isinstance(raw, str):
            return "condition"
        token = raw.strip()
        # Strip simple quoting and trailing punctuation
        token = token.strip("'\"`.;:,!?()[]{}")
        # Factor names are short identifiers. A multi-word reply (any
        # internal whitespace) is the LLM emitting a sentence, not a
        # column name, so abort instead of taking the first word.
        if not token or any(c.isspace() for c in token):
            return "condition"
        token = token.lower()
        # Must be a plausible identifier: letters/digits/underscore/dash,
        # 1-40 chars, starts with a letter.
        import re as _re
        if not _re.match(r"^[a-z][a-z0-9_\-]{0,39}$", token):
            return "condition"
        # Reject common imperative verbs that the LLM emits when confused
        if token in {"need", "paste", "please", "provide", "give", "type",
                     "enter", "specify", "choose", "select"}:
            return "condition"
        return token

    def _publish_batch_checkpoint(self):
        # Detect batch keywords
        batch_keywords = {"lane", "batch", "plate", "run", "flowcell"}
        possible_batch = set()
        for p in self._parsed_samples:
            for tok in p["tokens"]:
                if tok in batch_keywords or re.match(r'^(lane|batch)\d+$', tok):
                    possible_batch.add(tok)

        if not possible_batch:
            if self._batch_covariate:
                possible_batch.add(self._batch_covariate)
            else:
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
        warning_groups = []
        for group, members in (self._confirmed_groups or {}).items():
            roots = {self._technical_replicate_root(m) for m in members}
            if len(roots) < len(members):
                warning_groups.append(group)

        if not warning_groups:
            self._replicate_handling = self._build_replicate_handling(
                technical=False
            )
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
                "No — technical replicates; merge by biological unit",
                "Not sure — stop and provide a sample sheet",
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
            f"Experimental unit: {design['experimental_unit']}",
            f"Analysis design: {design['analysis_design_formula']}",
            f"Nominal residual df: {design['nominal_residual_df']}",
        ]
        if design.get("design_status") == "blocking":
            summary.append(
                "BLOCKING: " + ", ".join(design.get("blocking_reasons", []))
            )
            if design.get("confounded_covariates"):
                summary.append(
                    "Confounded with "
                    f"'{design['main_factor']}': "
                    f"{', '.join(design['confounded_covariates'])} "
                    "(non-identifiable; cannot be resolved by ignoring the term)."
                )
        if design.get("batch_covariate"):
            summary.append(f"Batch covariate: {design['batch_covariate']}")

        question = "Experimental design ready:\n\n" + "\n".join(summary) + "\n\nProceed with analysis?"
        options = (
            ["Yes — proceed", "Cancel"]
            if design.get("design_status") != "blocking"
            else ["Cancel — provide an explicit sample sheet / resolve confounding"]
        )
        msg = self.publish_escalation(
            experiment_id=self._experiment_id,
            checkpoint=2.6,
            question=question,
            options=options,
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
        # "treated vs control", "WT vs KO", "stim vs ctrl"). Surface it to
        # the LLM — without it the model can only pattern-match on
        # filenames, which fails when sample IDs are opaque (donor barcodes,
        # internal accession IDs, or single-file h5ad inputs).
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
        main_factor = self._main_factor or "condition"
        if self._batch_covariate:
            formula = f"~ {self._batch_covariate} + {main_factor}"
        else:
            formula = f"~ {main_factor}"

        handling = self._replicate_handling or self._build_replicate_handling(
            technical=False
        )
        replicates = {
            str(group): sum(
                1 for unit in handling["units"]
                if unit["group"] == str(group)
            )
            for group in (self._confirmed_groups or {})
        }
        n_input_libraries = sum(
            len(members) for members in (self._confirmed_groups or {}).values()
        )
        n_units = sum(replicates.values())
        n_groups = sum(1 for count in replicates.values() if count > 0)
        nominal_residual_df = max(0, n_units - n_groups)
        blocking_groups = (
            {
                group: count for group, count in replicates.items() if count < 2
            }
            if handling["mode"] == "technical_aggregate"
            else {}
        )
        if handling["mode"] == "technical_aggregate":
            analysis_design_formula = (
                f"sum technical libraries by biological_unit; {formula}"
            )
            experimental_unit = "biological_unit"
        else:
            analysis_design_formula = formula
            experimental_unit = "sample"
        analysis_groups = {
            str(group): [
                unit["id"] for unit in handling["units"]
                if unit["group"] == str(group)
            ]
            for group in (self._confirmed_groups or {})
        }
        pseudobulk_design = dict(self._pseudobulk_design or {})
        if handling["mode"] == "technical_aggregate":
            pseudobulk_design.setdefault("condition_col", main_factor)
            pseudobulk_design["replicate_col"] = "biological_unit"

        mapping_conflicts = handling.get("conflicts", []) or []
        confounded_covariates = self._confounded_covariates(main_factor)
        blocking_reasons = []
        if blocking_groups:
            blocking_reasons.append("insufficient_biological_replicates")
        if mapping_conflicts:
            blocking_reasons.append("replicate_mapping_conflict")
        if confounded_covariates:
            blocking_reasons.append("condition_covariate_confounding")

        return {
            "organism":        self._organism,
            "genome":          self._genome,
            "groups":          self._confirmed_groups,
            "analysis_groups": analysis_groups,
            "main_factor":     main_factor,
            "design_formula":  formula,
            "batch_covariate": self._batch_covariate,
            "pseudobulk":      pseudobulk_design,
            "source":          self._inferred_design.get("source"),
            "sample_aliases":  self._inferred_design.get("sample_aliases", {}),
            "replicate_handling": handling,
            "experimental_unit": experimental_unit,
            "analysis_design_formula": analysis_design_formula,
            "replicates":      replicates,
            "n_input_libraries": n_input_libraries,
            "n_total_samples": n_units,
            "nominal_residual_df": nominal_residual_df,
            "design_status": "blocking" if blocking_reasons else "ready",
            "blocking_groups": blocking_groups,
            "blocking_reasons": blocking_reasons,
            "confounded_covariates": confounded_covariates,
        }

    def _confounded_covariates(self, main_factor: str) -> list[str]:
        """Factors from the inferred sample sheet perfectly aliased with condition.

        Confounding is a property of the metadata, not of the CP2.4 batch choice:
        declining to model a confounded batch does not make biology identifiable,
        so the block cannot be bypassed by the headless "ignore batch" default (B2).
        The explicit sample sheet is authoritative; a precomputed flag from the
        connector is the fallback when no sheet is present.
        """
        from aria.utils.design_matrix import factors_confounded_with_condition

        inferred = self._inferred_design or {}
        condition_col = inferred.get("condition_col") or main_factor
        sheet = inferred.get("sample_sheet")
        covariates = inferred.get("covariates") or []
        if sheet and covariates:
            return factors_confounded_with_condition(
                sheet, condition_col, covariates
            )
        return list(inferred.get("confounded_covariates") or [])

    @staticmethod
    def _technical_replicate_root(sample: str) -> str:
        """Remove only an explicit technical suffix (``_rep1``, ``_r1``, ``_1``).

        Plain biological identifiers ending in a digit (for example ``donor1``)
        are preserved; the old permissive regex could truncate those names.
        """
        value = str(sample)
        root = re.sub(
            r"(?:[_-](?:rep|r)?\d+)$", "", value, flags=re.IGNORECASE
        )
        return root or value

    def _build_replicate_handling(self, *, technical: bool) -> dict:
        """Build the CP2.5 execution contract without changing raw group membership."""
        mode = "technical_aggregate" if technical else "independent_biological"
        units: list[dict] = []
        sample_to_unit: dict[str, str] = {}
        conflicts: list[str] = []
        for group, members in (self._confirmed_groups or {}).items():
            by_root: dict[str, list[str]] = {}
            for member in members or []:
                root = (
                    self._technical_replicate_root(member)
                    if technical else str(member)
                )
                by_root.setdefault(root, []).append(str(member))
            for root, unit_members in by_root.items():
                unit_id = f"{group}::{root}"
                units.append({
                    "id": unit_id,
                    "group": str(group),
                    "root": root,
                    "members": unit_members,
                })
                for member in unit_members:
                    prior = sample_to_unit.get(member)
                    if technical and prior is not None and prior != unit_id:
                        conflicts.append(member)
                    sample_to_unit[member] = unit_id
        return {
            "mode": mode,
            "source": "user_confirmed_cp2.5" if technical else "cp2.5_independent",
            "aggregation": "sum_raw_counts" if technical else "none",
            "units": units,
            "sample_to_unit": sample_to_unit,
            "conflicts": sorted(set(conflicts)),
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
