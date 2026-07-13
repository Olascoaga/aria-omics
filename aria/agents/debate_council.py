"""
ARIA DebateCouncil
------------------
Internal peer review system for high-stakes biological interpretations.

Prevents sycophancy by forcing the Critic to formulate the most
plausible alternative hypothesis BEFORE considering acceptance.

Design principles (from Gemini's review):
  1. Critic must state the alternative hypothesis explicitly
  2. Critic demands cross-validation, not just positive markers
  3. Critic is topologically aware — challenges algorithmic fit
  4. Debate uses CavemanMode.FULL internally (token efficiency)
  5. Only the final consensus is in normal prose (user-facing)

When to use:
  - Cell type annotation (confidence MEDIUM or LOW)
  - Integration decisions (RNA/ATAC weight imbalance)
  - Final narrative conclusions (overclaiming prevention)
  - Any finding where a wrong interpretation affects the manuscript

When NOT to use:
  - File detection (deterministic)
  - QC metrics (objective)
  - Parameter search space (heuristic)
  - Status messages

Max 3 debate rounds per finding (cost control).
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from aria.llm.provider import LLMProvider, TaskTier
from aria.llm.prompt_boundary import (
    PromptDataField,
    build_untrusted_prompt,
    system_with_untrusted_boundary,
)
from aria.bus.message_bus import CavemanMode

log = logging.getLogger("aria.debate")


class DebateVerdict(Enum):
    ACCEPT          = "accept"           # Critic accepts Proposer's claim
    ACCEPT_REVISED  = "accept_revised"   # Accept with mandatory revision
    REJECT          = "reject"           # Claim not supported by evidence
    INSUFFICIENT    = "insufficient"     # Not enough data to conclude either way


@dataclass
class DebateRound:
    round_number:   int
    proposer_claim: str
    critic_response: str
    alternative_hypothesis: str   # Critic's counter-hypothesis
    evidence_requested: list[str] # What the Critic demanded
    verdict:        Optional[DebateVerdict] = None


@dataclass
class DebateResult:
    topic:          str
    final_claim:    str
    verdict:        DebateVerdict
    rounds:         list[DebateRound]
    consensus:      str            # Final user-facing summary
    confidence:     str            # "high" / "medium" / "low" / "insufficient"
    limitations:    list[str]      # Mandatory limitations section
    was_revised:    bool = False   # Did Critic force a revision?


# ── System prompts ────────────────────────────────────────────────────────────

PROPOSER_SYSTEM = """
You are ARIA's Proposer — a bioinformatics expert presenting a biological interpretation.

Your role:
- State your claim clearly with supporting evidence
- Cite specific metrics, marker genes, or statistical values
- Be direct — do not hedge excessively
- Acknowledge known weaknesses proactively

Format your response as:
CLAIM: [your interpretation]
EVIDENCE: [specific supporting data]
KNOWN_LIMITATIONS: [what could challenge this]
""".strip()


CRITIC_SYSTEM = """
You are ARIA's Critic — a rigorous scientific reviewer whose ONLY job is to find flaws.

Your mandate (from Gemini's review framework):

1. ALTERNATIVE HYPOTHESIS FIRST
   Before evaluating the Proposer's claim, you MUST state the most plausible
   alternative biological or technical explanation. Examples:
   - "This could be a doublet artifact rather than a transitional cell type"
   - "The clustering topology may be driven by library size, not biology"
   - "This TF motif enrichment could reflect open chromatin noise, not regulation"

2. DEMAND CROSS-VALIDATION
   Never accept a positive claim without demanding negative evidence.
   For cell type annotation: require ABSENCE of conflicting lineage markers,
   not just PRESENCE of positive markers.
   For TF activity: require footprinting support, not just motif enrichment.

3. TOPOLOGICAL AWARENESS
   Challenge algorithmic fit for the experimental design:
   - Is LIONESS appropriate for single-sample GRN inference?
   - Is WNN valid when RNA quality is poor?
   - Does this HiC resolution support the TAD calls being made?

4. VERDICT OPTIONS:
   ACCEPT — evidence is sufficient and alternative hypotheses are ruled out
   ACCEPT_REVISED — claim is plausible but must be reworded to be more conservative
   REJECT — evidence is insufficient or alternative hypothesis is more likely
   INSUFFICIENT — not enough data to conclude; recommend additional experiments

Respond with ONLY valid JSON (no preamble, no markdown fences) matching:
{
  "alternative_hypothesis": "string — the most likely alternative explanation",
  "challenges":             ["string", ...]  — specific evidence gaps,
  "evidence_requested":     ["string", ...]  — what would change your verdict,
  "verdict":                "ACCEPT" | "ACCEPT_REVISED" | "REJECT" | "INSUFFICIENT",
  "revised_claim":          "string — more conservative wording, empty if not ACCEPT_REVISED"
}
""".strip()


CONSENSUS_SYSTEM = """
You are ARIA's NarrativeAgent summarizing a completed peer review debate.
Write a clear, concise scientific summary in normal prose (not caveman mode).
Include: the final accepted claim, key caveats, and mandatory limitations.
This will appear in the final report — be precise and honest about uncertainty.
""".strip()


class DebateCouncil:
    """
    Two-agent internal peer review for high-stakes biological interpretations.

    The Proposer presents a claim with evidence.
    The Critic challenges it by first stating the alternative hypothesis.
    They iterate up to max_rounds times.
    The final consensus is written in normal prose for the report.

    All inter-agent communication uses CavemanMode.FULL (token efficient).
    Only the final consensus uses normal language (user-facing).
    """

    def __init__(self, llm: LLMProvider, max_rounds: int = 3):
        self.llm        = llm
        self.max_rounds = max_rounds

    def resolve(
        self,
        topic:              str,
        initial_claim:      str,
        evidence:           dict,
        biological_context: dict = None,
    ) -> DebateResult:
        """
        Run a full debate on a biological interpretation.

        Args:
            topic:              What is being debated ("cell_type_annotation", etc.)
            initial_claim:      The Proposer's starting interpretation
            evidence:           Supporting data (markers, metrics, etc.)
            biological_context: From OrchestratorAgent (organism, condition, etc.)

        Returns:
            DebateResult with final consensus and mandatory limitations
        """
        rounds     = []
        current_claim = initial_claim
        context_str   = self._format_context(biological_context or {})
        evidence_str  = self._format_evidence(evidence)

        for round_num in range(1, self.max_rounds + 1):
            log.debug(f"Debate round {round_num}/{self.max_rounds}: {topic}")

            # ── Proposer presents / defends ──────────────────────────────
            proposer_response = self._proposer_turn(
                claim=current_claim,
                evidence=evidence_str,
                context=context_str,
                round_number=round_num,
                previous_rounds=rounds,
            )

            # ── Critic challenges ────────────────────────────────────────
            critic_response = self._critic_turn(
                proposer_claim=proposer_response,
                evidence=evidence_str,
                context=context_str,
                topic=topic,
            )

            # ── Parse critic verdict ─────────────────────────────────────
            verdict, alternative, challenges, revised = self._parse_critic(
                critic_response
            )

            round_result = DebateRound(
                round_number=round_num,
                proposer_claim=proposer_response,
                critic_response=critic_response,
                alternative_hypothesis=alternative,
                evidence_requested=challenges,
                verdict=verdict,
            )
            rounds.append(round_result)

            # ── Check if debate should end ───────────────────────────────
            if verdict in (DebateVerdict.ACCEPT, DebateVerdict.INSUFFICIENT):
                break

            if verdict == DebateVerdict.ACCEPT_REVISED and revised:
                current_claim = revised
                round_result.verdict = DebateVerdict.ACCEPT_REVISED
                break

            if verdict == DebateVerdict.REJECT:
                # Give Proposer one more round to address the rejection
                if round_num == self.max_rounds:
                    break  # No more rounds — reject stands
                current_claim = f"{current_claim}\n[Addressing critic: {challenges}]"

        # ── Final consensus ──────────────────────────────────────────────
        final_verdict = rounds[-1].verdict or DebateVerdict.INSUFFICIENT
        was_revised   = any(r.verdict == DebateVerdict.ACCEPT_REVISED for r in rounds)

        consensus, limitations = self._write_consensus(
            topic=topic,
            final_claim=current_claim,
            rounds=rounds,
            verdict=final_verdict,
        )

        confidence = self._verdict_to_confidence(final_verdict, len(rounds))

        return DebateResult(
            topic=topic,
            final_claim=current_claim,
            verdict=final_verdict,
            rounds=rounds,
            consensus=consensus,
            confidence=confidence,
            limitations=limitations,
            was_revised=was_revised,
        )

    # ── LLM calls ─────────────────────────────────────────────────────────

    def _proposer_turn(self, claim: str, evidence: str, context: str,
                       round_number: int, previous_rounds: list) -> str:
        history = None
        if previous_rounds:
            last = previous_rounds[-1]
            history = last.critic_response

        prompt = build_untrusted_prompt(
            task=(
                "Present or defend the biological interpretation using only the "
                "available evidence. Address the previous critic challenge when present."
            ),
            fields=[
                PromptDataField("topic", claim, "generated_text", "proposer"),
                PromptDataField("evidence", evidence, "structured_result",
                                "agent_results"),
                PromptDataField("biological_context", context, "metadata",
                                "experiment_context"),
                PromptDataField("round", round_number, "metadata", "debate"),
                PromptDataField("previous_critic_challenge", history,
                                "generated_text", "critic"),
            ],
            response_contract="Return concise scientific prose.",
        )
        return self.llm.complete(
            prompt=prompt,
            system=system_with_untrusted_boundary(PROPOSER_SYSTEM),
            tier=TaskTier.MEDIUM,
            max_tokens=600,
        )

    def _critic_turn(self, proposer_claim: str, evidence: str,
                     context: str, topic: str) -> str:
        prompt = build_untrusted_prompt(
            task=(
                "Challenge the interpretation. State the most plausible alternative "
                "hypothesis first and require evidence for acceptance."
            ),
            fields=[
                PromptDataField("topic", topic, "generated_text", "debate"),
                PromptDataField("proposer_claim", proposer_claim,
                                "generated_text", "proposer"),
                PromptDataField("evidence", evidence, "structured_result",
                                "agent_results"),
                PromptDataField("biological_context", context, "metadata",
                                "experiment_context"),
            ],
            response_contract="Return the challenge as scientific prose.",
        )
        return self.llm.complete(
            prompt=prompt,
            system=system_with_untrusted_boundary(CRITIC_SYSTEM),
            tier=TaskTier.MEDIUM,
            max_tokens=600,
        )

    def _write_consensus(self, topic: str, final_claim: str,
                          rounds: list, verdict: DebateVerdict
                          ) -> tuple[str, list[str]]:
        """Write the final user-facing consensus in normal prose."""
        rounds_summary = [{
            "round": r.round_number,
            "verdict": r.verdict.value if r.verdict else "ongoing",
            "alternative_hypothesis": r.alternative_hypothesis,
        } for r in rounds]

        prompt = build_untrusted_prompt(
            task=(
                "Write a concise scientific summary of the final accepted "
                "interpretation (2-3 sentences), followed by a Limitations list. "
                "If the code-owned verdict is reject or insufficient, say so clearly."
            ),
            fields=[
                PromptDataField("topic", topic, "generated_text", "debate"),
                PromptDataField("final_claim", final_claim, "generated_text",
                                "debate"),
                PromptDataField("verdict", verdict.value, "metadata", "code"),
                PromptDataField("rounds", rounds_summary, "structured_result",
                                "debate"),
            ],
            response_contract="Return prose followed by bullet-point limitations.",
        )
        response = self.llm.complete(
            prompt=prompt,
            system=system_with_untrusted_boundary(CONSENSUS_SYSTEM),
            tier=TaskTier.HEAVY,    # Normal prose for user — use best model
            max_tokens=400,
        )

        # Extract limitations if present
        limitations = []
        if "Limitations" in response or "limitation" in response.lower():
            lines = response.split("\n")
            in_limitations = False
            for line in lines:
                if "limitation" in line.lower():
                    in_limitations = True
                    continue
                if in_limitations and line.strip().startswith(("-", "*", "•")):
                    limitations.append(line.strip().lstrip("-*• "))
                elif in_limitations and line.strip() == "":
                    continue
                elif in_limitations:
                    in_limitations = False

        return response, limitations

    # ── Parsing ───────────────────────────────────────────────────────────

    _VERDICT_MAP = {
        "ACCEPT":         DebateVerdict.ACCEPT,
        "ACCEPT_REVISED": DebateVerdict.ACCEPT_REVISED,
        "REJECT":         DebateVerdict.REJECT,
        "INSUFFICIENT":   DebateVerdict.INSUFFICIENT,
    }

    def _parse_critic(self, critic_response: str
                       ) -> tuple[DebateVerdict, str, list[str], str]:
        """
        Extract structured fields from critic response.
        Primary path: JSON (as instructed in CRITIC_SYSTEM).
        Fallback:     legacy line-prefix parser (handles older formats).
        """
        # Primary: JSON
        json_blob = self._extract_json(critic_response)
        if json_blob is not None:
            try:
                data = json.loads(json_blob)
                verdict = self._VERDICT_MAP.get(
                    str(data.get("verdict", "")).upper().strip(),
                    DebateVerdict.INSUFFICIENT,
                )
                alternative = str(data.get("alternative_hypothesis", "") or "").strip()
                raw_challenges = data.get("challenges", []) or []
                raw_evidence   = data.get("evidence_requested", []) or []
                if isinstance(raw_challenges, str):
                    raw_challenges = [raw_challenges]
                if isinstance(raw_evidence, str):
                    raw_evidence = [raw_evidence]
                challenges = [str(c).strip() for c in raw_challenges if str(c).strip()]
                challenges.extend(str(e).strip() for e in raw_evidence if str(e).strip())
                revised_claim = str(data.get("revised_claim", "") or "").strip()
                return verdict, alternative, challenges, revised_claim
            except (json.JSONDecodeError, TypeError, AttributeError) as e:
                log.debug(f"Critic JSON parse failed: {e}; falling back to line parser")

        # Fallback: legacy prefix parser (case-insensitive, tolerates markdown)
        verdict       = DebateVerdict.INSUFFICIENT
        alternative   = ""
        challenges    = []
        revised_claim = ""

        for raw_line in critic_response.split("\n"):
            # Strip markdown inline marks (** # `) anywhere on the line.
            # Keep underscores: the section names contain them (e.g.
            # ALTERNATIVE_HYPOTHESIS) and stripping them silently broke the
            # legacy fallback path.
            line = re.sub(r"[*#`]+", "", raw_line).strip()
            upper = line.upper()
            if upper.startswith("VERDICT:"):
                v = line.split(":", 1)[1].strip().upper()
                verdict = self._VERDICT_MAP.get(v, DebateVerdict.INSUFFICIENT)
            elif upper.startswith("ALTERNATIVE_HYPOTHESIS:") or \
                 upper.startswith("ALTERNATIVE HYPOTHESIS:"):
                alternative = line.split(":", 1)[1].strip()
            elif upper.startswith("CHALLENGES:") or \
                 upper.startswith("EVIDENCE_REQUESTED:") or \
                 upper.startswith("EVIDENCE REQUESTED:"):
                challenges.append(line.split(":", 1)[1].strip())
            elif upper.startswith("REVISED_CLAIM:") or \
                 upper.startswith("REVISED CLAIM:"):
                revised_claim = line.split(":", 1)[1].strip()

        return verdict, alternative, challenges, revised_claim

    @staticmethod
    def _extract_json(text: str) -> Optional[str]:
        """Return the first balanced JSON object found in text, or None."""
        # Fast path: text is already pure JSON
        stripped = text.strip()
        if stripped.startswith("{") and stripped.endswith("}"):
            return stripped

        # Strip ```json ... ``` fences
        fence = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", text, re.DOTALL)
        if fence:
            return fence.group(1)

        # Last resort: greedy first {...} block
        m = re.search(r"\{.*\}", text, re.DOTALL)
        return m.group(0) if m else None

    # ── Helpers ───────────────────────────────────────────────────────────

    def _format_evidence(self, evidence: dict) -> str:
        if not evidence:
            return "No additional evidence provided."
        parts = []
        for k, v in evidence.items():
            if isinstance(v, list):
                parts.append(f"{k}: {', '.join(str(x) for x in v[:10])}")
            elif isinstance(v, dict):
                parts.append(f"{k}: {list(v.items())[:5]}")
            else:
                parts.append(f"{k}: {v}")
        return "\n".join(parts)

    def _format_context(self, context: dict) -> str:
        if not context:
            return "No biological context provided."
        return (
            f"Organism: {context.get('organism', 'unknown')}\n"
            f"Question: {context.get('user_question', 'unknown')}\n"
            f"Analysis type: {context.get('analysis_type', 'unknown')}"
        )

    def _verdict_to_confidence(self, verdict: DebateVerdict,
                                n_rounds: int) -> str:
        if verdict == DebateVerdict.ACCEPT and n_rounds == 1:
            return "high"
        elif verdict in (DebateVerdict.ACCEPT, DebateVerdict.ACCEPT_REVISED):
            return "medium"
        elif verdict == DebateVerdict.REJECT:
            return "low"
        else:
            return "insufficient"
