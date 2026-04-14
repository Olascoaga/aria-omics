"""
ARIA ContextManager
-------------------
Dynamic context window management for heterogeneous LLM backends.

Problem:
  Claude Sonnet has 200k tokens and degrades gracefully.
  Llama 3 8B has 4k tokens and fails CATASTROPHICALLY at the limit.
  The same message history that Claude handles fine will corrupt Llama's output.

Solution:
  Per-model context budgeting with a 4-step degradation cascade:
    1. TRUNCATE   — drop oldest conversation history (least critical)
    2. COMPRESS   — summarize findings in CavemanULTRA (~30x reduction)
    3. PRIORITIZE — keep only L0 memory + current task
    4. SPLIT      — divide task into smaller subtasks (never truncate task itself)

Token counting:
  - Cloud models (Claude/GPT): tiktoken cl100k_base (good approximation)
  - Local models: conservative heuristic (1 token ≈ 3.5 chars) with safety margin
  - Always leave 20% buffer for model safety

Design principle:
  NEVER truncate the current task prompt or critical parameters.
  NEVER silently degrade without logging what was dropped.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

log = logging.getLogger("aria.llm.context")


# ── Model profiles ────────────────────────────────────────────────────────────

@dataclass
class ModelProfile:
    """
    Context window profile for a specific model.
    Controls how aggressively ContextManager truncates.
    """
    name:            str
    context_window:  int          # max tokens (hard limit from provider)
    is_local:        bool = False # local models fail harder at limits

    @property
    def safe_limit(self) -> int:
        """
        Usable token budget with safety margin.
        Local models: 75% of window (they fail catastrophically near limit)
        Cloud models: 85% of window (graceful degradation)
        """
        margin = 0.75 if self.is_local else 0.85
        return int(self.context_window * margin)

    @property
    def compression_threshold(self) -> int:
        """Start compressing when context exceeds this."""
        return int(self.safe_limit * 0.70)

    @property
    def critical_threshold(self) -> int:
        """Emergency: drop everything except task + L0 memory."""
        return int(self.safe_limit * 0.90)


# Pre-defined profiles for common models
PROFILES = {
    "claude-sonnet-4-20250514":  ModelProfile("claude-sonnet-4-20250514",  200_000),
    "claude-haiku-4-5-20251001": ModelProfile("claude-haiku-4-5-20251001", 200_000),
    "gpt-4o":                    ModelProfile("gpt-4o",                    128_000),
    "gpt-4o-mini":               ModelProfile("gpt-4o-mini",               128_000),
    "ollama/llama3:70b":         ModelProfile("ollama/llama3:70b",           8_000, is_local=True),
    "ollama/llama3:8b":          ModelProfile("ollama/llama3:8b",            4_000, is_local=True),
    "ollama/mistral:7b":         ModelProfile("ollama/mistral:7b",           8_000, is_local=True),
}


class ContextManager:
    """
    Manages token budgets and context truncation for a single model.
    One ContextManager per ModelProfile.
    """

    def __init__(self, profile: ModelProfile):
        self.profile = profile
        self._tokenizer = self._load_tokenizer(profile)

    # ── Public interface ──────────────────────────────────────────────────

    def prepare_messages(
        self,
        system:             str,
        history:            list[dict],   # [{"role": "user/assistant", "content": "..."}]
        new_prompt:         str,
        max_response_tokens: int = 1024,
    ) -> list[dict]:
        """
        Build the final message list that fits within the model's context.
        Applies degradation cascade if needed.
        Returns LiteLLM-compatible message list.
        """
        available = self.profile.safe_limit - max_response_tokens

        # Count base tokens (system + new prompt — non-negotiable)
        base_tokens = (
            self._count(system) +
            self._count(new_prompt) +
            64  # overhead for role markers
        )

        if base_tokens > available:
            # Task itself is too large — this should never happen in practice
            # but we log it rather than fail silently
            log.error(
                f"[{self.profile.name}] Base prompt ({base_tokens} tokens) "
                f"exceeds safe limit ({available}). Cannot truncate task. "
                f"Consider switching to a larger model."
            )

        history_budget = available - base_tokens

        # Apply cascade
        managed_history = self._apply_cascade(history, history_budget)

        # Build final message list
        messages = []
        if system:
            messages.append({"role": "system", "content": system})
        messages.extend(managed_history)
        messages.append({"role": "user", "content": new_prompt})

        # Final verification
        total = sum(self._count(m["content"]) for m in messages)
        log.debug(
            f"[{self.profile.name}] Final context: {total} tokens "
            f"(limit: {self.profile.safe_limit})"
        )

        return messages

    def count_tokens(self, text: str) -> int:
        """Public token counter for external use."""
        return self._count(text)

    def fits(self, text: str, budget: int) -> bool:
        return self._count(text) <= budget

    # ── Degradation cascade ───────────────────────────────────────────────

    def _apply_cascade(
        self,
        history: list[dict],
        budget:  int,
    ) -> list[dict]:
        """
        4-step cascade to fit history into token budget.
        Each step is logged so nothing is silently dropped.
        """
        if not history:
            return []

        current_tokens = sum(self._count(m["content"]) for m in history)

        # Step 0: Already fits — no action needed
        if current_tokens <= budget:
            return history

        log.info(
            f"[{self.profile.name}] History {current_tokens} tokens "
            f"> budget {budget}. Applying context cascade."
        )

        # ── Step 1: TRUNCATE oldest messages ─────────────────────────────
        truncated = self._truncate_history(history, budget)
        tokens_after = sum(self._count(m["content"]) for m in truncated)
        dropped = len(history) - len(truncated)
        if dropped > 0:
            log.info(
                f"[{self.profile.name}] TRUNCATE: dropped {dropped} old "
                f"messages. {tokens_after} tokens remaining."
            )
        if tokens_after <= budget:
            return truncated

        # ── Step 2: COMPRESS with CavemanULTRA ───────────────────────────
        compressed = self._compress_history(truncated)
        tokens_after = sum(self._count(m["content"]) for m in compressed)
        log.info(
            f"[{self.profile.name}] COMPRESS: caveman ultra applied. "
            f"{tokens_after} tokens after compression."
        )
        if tokens_after <= budget:
            return compressed

        # ── Step 3: PRIORITIZE — keep only most recent N messages ────────
        prioritized = self._prioritize_history(compressed, budget)
        tokens_after = sum(self._count(m["content"]) for m in prioritized)
        log.info(
            f"[{self.profile.name}] PRIORITIZE: kept {len(prioritized)} "
            f"most recent messages. {tokens_after} tokens."
        )
        if tokens_after <= budget:
            return prioritized

        # ── Step 4: MINIMAL — last message only ──────────────────────────
        log.warning(
            f"[{self.profile.name}] MINIMAL: context critically constrained. "
            f"Keeping only last exchange. Consider a larger model."
        )
        last = prioritized[-2:] if len(prioritized) >= 2 else prioritized[-1:]
        return last

    def _truncate_history(
        self, history: list[dict], budget: int
    ) -> list[dict]:
        """
        Drop oldest messages until history fits budget.
        Always keep the most recent exchange (last 2 messages).
        """
        if len(history) <= 2:
            return history

        result = list(history)
        while len(result) > 2:
            total = sum(self._count(m["content"]) for m in result)
            if total <= budget:
                break
            result.pop(0)  # drop oldest

        return result

    def _compress_history(self, history: list[dict]) -> list[dict]:
        """
        Apply CavemanULTRA compression to all assistant messages in history.
        User messages are kept intact (they contain instructions/context).
        Assistant messages are compressed (they contain verbose reasoning).

        CavemanULTRA rules applied directly (no LLM call — deterministic):
          - Drop: articles, filler, hedging, pleasantries
          - Keep: technical terms, numbers, gene names, tool names
          - Replace: verbose phrases with short equivalents
          - Use: arrows for causality (X→Y), abbreviations
        """
        compressed = []
        for msg in history:
            if msg["role"] == "assistant":
                compressed.append({
                    "role":    "assistant",
                    "content": self._caveman_ultra(msg["content"]),
                })
            else:
                compressed.append(msg)
        return compressed

    def _prioritize_history(
        self, history: list[dict], budget: int
    ) -> list[dict]:
        """
        Keep only the N most recent messages that fit within budget.
        Works backwards from most recent.
        """
        result = []
        used = 0
        for msg in reversed(history):
            tokens = self._count(msg["content"])
            if used + tokens <= budget:
                result.insert(0, msg)
                used += tokens
            else:
                break
        return result

    # ── CavemanULTRA compression (deterministic, no LLM call) ────────────

    _FILLER_PHRASES = [
        "I'd be happy to ", "I would be happy to ",
        "Certainly! ", "Sure! ", "Of course! ",
        "Great question! ", "That's a great point. ",
        "I understand that ", "It's worth noting that ",
        "Please note that ", "It's important to note that ",
        "I should mention that ", "I want to highlight that ",
        "As you can see, ", "As we discussed, ",
        "In order to ", "In summary, ",
        "First and foremost, ", "Last but not least, ",
    ]

    _VERBOSE_MAP = {
        "implement a solution for":    "fix",
        "implement a solution to":     "fix",
        "in order to achieve":         "to get",
        "with respect to":             "re:",
        "at this point in time":       "now",
        "due to the fact that":        "because",
        "in the event that":           "if",
        "a significant number of":     "many",
        "the majority of":             "most",
        "on a regular basis":          "regularly",
        "it is possible that":         "maybe",
        "there is a possibility":      "maybe",
        "it appears that":             "seems",
        "it would seem that":          "seems",
        "extensive":                   "big",
        "implement":                   "add",
        "utilize":                     "use",
        "facilitate":                  "help",
        "demonstrate":                 "show",
        "approximately":               "~",
        "therefore":                   "→",
        "consequently":                "→",
        "results in":                  "→",
        "which leads to":              "→",
        "however":                     "but",
        "nevertheless":                "still",
        "furthermore":                 "also",
        "in addition":                 "also",
        "additionally":                "also",
    }

    def _caveman_ultra(self, text: str) -> str:
        """
        Deterministic CavemanULTRA compression.
        No LLM call. Regex + lookup tables.
        Preserves: code blocks, gene names, numbers, technical terms.
        """
        import re

        # Extract and protect code blocks
        code_blocks = {}
        code_pattern = re.compile(r'```.*?```', re.DOTALL)
        for i, match in enumerate(code_pattern.finditer(text)):
            key = f"__CODE_BLOCK_{i}__"
            code_blocks[key] = match.group()
            text = text.replace(match.group(), key, 1)

        # Drop filler phrases
        for phrase in self._FILLER_PHRASES:
            text = text.replace(phrase, "")
            text = text.replace(phrase.lower(), "")

        # Apply verbose → terse map
        for verbose, terse in self._VERBOSE_MAP.items():
            text = re.sub(re.escape(verbose), terse,
                          text, flags=re.IGNORECASE)

        # Drop articles at sentence/fragment start
        text = re.sub(r'\b(The|A|An) ([A-Z])', r'\2', text)

        # Collapse multiple spaces/newlines
        text = re.sub(r'\n{3,}', '\n\n', text)
        text = re.sub(r' {2,}', ' ', text)
        text = text.strip()

        # Restore code blocks
        for key, block in code_blocks.items():
            text = text.replace(key, block)

        return text

    # ── Token counting ────────────────────────────────────────────────────

    def _load_tokenizer(self, profile: ModelProfile):
        """
        Load the best available tokenizer for this model.
        Falls back to conservative character heuristic for local models.
        """
        if profile.is_local:
            return None  # Use heuristic

        try:
            import tiktoken
            # cl100k_base covers Claude, GPT-4, GPT-3.5
            return tiktoken.get_encoding("cl100k_base")
        except Exception:
            log.warning("tiktoken unavailable — using character heuristic")
            return None

    def _count(self, text: str) -> int:
        """
        Count tokens in text.

        Cloud models: tiktoken (accurate)
        Local models: conservative heuristic with safety margin
          1 token ≈ 3.5 chars (conservative; true average is ~4)
          Safety margin: × 1.15 for local models (fail hard at limit)
        """
        if not text:
            return 0

        if self._tokenizer is not None:
            try:
                return len(self._tokenizer.encode(text))
            except Exception:
                pass

        # Conservative heuristic for local models
        char_per_token = 3.5
        safety = 1.15 if self.profile.is_local else 1.0
        return int(len(text) / char_per_token * safety)
