"""
ARIA LLMProvider
----------------
Universal LLM abstraction layer built on LiteLLM.

Decouples every agent from any specific provider.
Supports: Anthropic, OpenAI, Google Gemini, Ollama, vLLM, and any
          LiteLLM-compatible endpoint — swapped via config, zero code change.

Routing strategy (LLMOps best practice):
  HEAVY  → frontier model  (complex reasoning, integration, narrative)
  MEDIUM → mid-size model  (analysis decisions, parameter advice)
  LIGHT  → small/local     (caveman compression, classification, status)

Each tier falls back to the next available model automatically.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
import logging
from datetime import datetime, timezone
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional

import litellm
from litellm import completion

from aria.llm.context_manager import ContextManager, ModelProfile
from aria.utils.env_loader import load_aria_env
from aria.utils.provenance import record_llm_usage

log = logging.getLogger("aria.llm")

# Silence LiteLLM's aggressive logging
litellm.suppress_debug_info = True
logging.getLogger("LiteLLM").setLevel(logging.WARNING)


# ── Routing tiers ────────────────────────────────────────────────────────────

class TaskTier(Enum):
    """
    Complexity tiers for LLM routing.

    HEAVY  — biological reasoning, cross-modal integration, report generation
    MEDIUM — parameter decisions, QC interpretation, structured analysis
    LIGHT  — caveman compression, file classification, status summaries
    """
    HEAVY  = "heavy"
    MEDIUM = "medium"
    LIGHT  = "light"


@dataclass
class ModelConfig:
    """Configuration for a single model endpoint."""
    provider:    str            # "anthropic", "openai", "ollama", "gemini"
    model:       str            # model identifier for LiteLLM
    context_window: int         # max tokens this model handles reliably
    is_local:    bool = False   # True for Ollama/vLLM (no API cost)
    api_base:    Optional[str] = None  # for local endpoints


# ── Default model profiles (overridden by config.yaml) ───────────────────────

DEFAULT_MODELS: dict[TaskTier, list[ModelConfig]] = {
    TaskTier.HEAVY: [
        ModelConfig("anthropic", "claude-opus-4-7",            200_000),
        ModelConfig("anthropic", "claude-sonnet-4-6",          200_000),
        ModelConfig("openai",    "gpt-4o",                     128_000),
        ModelConfig("gemini",    "gemini/gemini-1.5-pro",      1_000_000),
        ModelConfig("ollama",    "ollama/llama3:70b",           8_000,  is_local=True,
                    api_base="http://localhost:11434"),
    ],
    TaskTier.MEDIUM: [
        ModelConfig("anthropic", "claude-sonnet-4-6",          200_000),
        ModelConfig("anthropic", "claude-haiku-4-5-20251001",  200_000),
        ModelConfig("openai",    "gpt-4o-mini",                128_000),
        ModelConfig("ollama",    "ollama/llama3:8b",            8_000,  is_local=True,
                    api_base="http://localhost:11434"),
        ModelConfig("ollama",    "ollama/mistral:7b",           8_000,  is_local=True,
                    api_base="http://localhost:11434"),
    ],
    TaskTier.LIGHT: [
        ModelConfig("anthropic", "claude-haiku-4-5-20251001",  200_000),
        ModelConfig("ollama",    "ollama/llama3:8b",            8_000,  is_local=True,
                    api_base="http://localhost:11434"),
        ModelConfig("ollama",    "ollama/mistral:7b",           8_000,  is_local=True,
                    api_base="http://localhost:11434"),
    ],
}


class LLMProvider:
    """
    Universal LLM provider for all ARIA agents.

    Usage:
        provider = LLMProvider.from_config()
        response = provider.complete(
            prompt="...",
            system="...",
            tier=TaskTier.HEAVY,
            max_tokens=1024
        )
    """

    DETERMINISTIC_TEMPERATURE = 0.0
    DETERMINISTIC_SEED = 0

    def __init__(
        self,
        models:      dict[TaskTier, list[ModelConfig]] = None,
        api_keys:    dict[str, str] = None,
        cache_dir:   Optional[str] = None,
    ):
        load_aria_env()
        self.models   = models or DEFAULT_MODELS
        self.api_keys = api_keys or {}
        self._context_managers: dict[str, ContextManager] = {}
        self._inject_api_keys()
        # File-backed prompt cache. Disabled when ARIA_LLM_CACHE=0.
        # Cache key: sha256(model + system + prompt + max_tokens + deterministic
        # generation controls).
        self._cache_enabled = os.environ.get("ARIA_LLM_CACHE", "1") != "0"
        if self._cache_enabled:
            base = cache_dir or os.environ.get("ARIA_LLM_CACHE_DIR") \
                   or str(Path.home() / ".aria" / "llm_cache")
            self._cache_dir = Path(base)
            try:
                self._cache_dir.mkdir(parents=True, exist_ok=True)
            except OSError:
                self._cache_enabled = False
                self._cache_dir = None
        else:
            self._cache_dir = None

    # ── Public interface ─────────────────────────────────────────────────

    def complete(
        self,
        prompt:     str,
        system:     str       = "",
        tier:       TaskTier  = TaskTier.MEDIUM,
        max_tokens: int       = 1024,
        messages:   list      = None,   # full conversation history (optional)
    ) -> str:
        """
        Complete a prompt using the appropriate model for the task tier.
        Automatically falls back through the tier's model list on failure.
        Context is managed and truncated per model profile.

        Returns the response text.
        """
        candidates = self.models.get(tier, self.models[TaskTier.MEDIUM])

        last_error = None
        for model_cfg in candidates:
            try:
                return self._call(model_cfg, prompt, system,
                                  max_tokens, messages, tier=tier)
            except Exception as e:
                log.warning(
                    f"Model {model_cfg.model} failed: {e}. "
                    f"Trying next fallback..."
                )
                last_error = e
                time.sleep(0.5)

        raise RuntimeError(
            f"All models failed for tier {tier.value}. "
            f"Last error: {last_error}"
        )

    def complete_heavy(self, prompt: str, system: str = "",
                       max_tokens: int = 2048, messages: list = None) -> str:
        return self.complete(prompt, system, TaskTier.HEAVY,
                             max_tokens, messages)

    def complete_medium(self, prompt: str, system: str = "",
                        max_tokens: int = 1024, messages: list = None) -> str:
        return self.complete(prompt, system, TaskTier.MEDIUM,
                             max_tokens, messages)

    def complete_light(self, prompt: str, system: str = "",
                       max_tokens: int = 512, messages: list = None) -> str:
        return self.complete(prompt, system, TaskTier.LIGHT,
                             max_tokens, messages)

    def get_active_model(self, tier: TaskTier) -> Optional[ModelConfig]:
        """Return the first available model for a tier (for display)."""
        candidates = self.models.get(tier, [])
        return candidates[0] if candidates else None

    # ── Private methods ──────────────────────────────────────────────────

    def _call(
        self,
        cfg:        ModelConfig,
        prompt:     str,
        system:     str,
        max_tokens: int,
        messages:   list = None,
        tier:       TaskTier = TaskTier.MEDIUM,
    ) -> str:
        """Execute a single LiteLLM call with context management."""

        # Cache lookup (only for single-prompt calls — multi-turn histories
        # are skipped because the cache key would explode in size).
        cache_key = None
        if self._cache_enabled and messages is None:
            cache_key = self._cache_key(
                cfg.model,
                system,
                prompt,
                max_tokens,
                self.DETERMINISTIC_TEMPERATURE,
                self.DETERMINISTIC_SEED,
            )
            cached = self._cache_get(cache_key)
            if cached is not None:
                log.debug(f"LLM cache hit: {cfg.model} key={cache_key[:8]}")
                record_llm_usage({
                    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                    "provider": cfg.provider,
                    "model": cfg.model,
                    "tier": tier.value,
                    "temperature": self.DETERMINISTIC_TEMPERATURE,
                    "seed": self.DETERMINISTIC_SEED,
                    "deterministic": True,
                    "cache_hit": True,
                    "prompt_tokens": 0,
                    "completion_tokens": 0,
                    "total_tokens": 0,
                    "estimated_cost_usd": 0.0,
                })
                return cached

        # Get or build ContextManager for this model
        ctx_mgr = self._get_context_manager(cfg)

        # Build message list
        if messages:
            # Manage context window for this specific model
            msgs = ctx_mgr.prepare_messages(
                system=system,
                history=messages,
                new_prompt=prompt,
                max_response_tokens=max_tokens,
            )
        else:
            msgs = ctx_mgr.prepare_messages(
                system=system,
                history=[],
                new_prompt=prompt,
                max_response_tokens=max_tokens,
            )

        kwargs = dict(
            model=cfg.model,
            messages=msgs,
            max_tokens=max_tokens,
            temperature=self.DETERMINISTIC_TEMPERATURE,
            seed=self.DETERMINISTIC_SEED,
        )

        # Inject API base for local models
        if cfg.api_base:
            kwargs["api_base"] = cfg.api_base

        # Inject API key for cloud providers
        api_key = self.api_keys.get(cfg.provider) or \
                  os.environ.get(self._key_env(cfg.provider))
        if api_key and not cfg.is_local:
            kwargs["api_key"] = api_key

        response = completion(**kwargs)
        text = response.choices[0].message.content
        usage = getattr(response, "usage", None)
        prompt_tokens = int(getattr(usage, "prompt_tokens", 0) or 0)
        completion_tokens = int(getattr(usage, "completion_tokens", 0) or 0)
        total_tokens = int(
            getattr(usage, "total_tokens", 0)
            or (prompt_tokens + completion_tokens)
        )
        try:
            cost = float(litellm.completion_cost(completion_response=response) or 0.0)
        except Exception:
            cost = 0.0
        record_llm_usage({
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "provider": cfg.provider,
            "model": cfg.model,
            "tier": tier.value,
            "temperature": self.DETERMINISTIC_TEMPERATURE,
            "seed": self.DETERMINISTIC_SEED,
            "deterministic": True,
            "cache_hit": False,
            "prompt_tokens": prompt_tokens,
            "completion_tokens": completion_tokens,
            "total_tokens": total_tokens,
            "estimated_cost_usd": cost,
        })
        if cache_key is not None and text:
            self._cache_put(cache_key, text)
        return text

    # ── Prompt cache helpers ─────────────────────────────────────────────

    @staticmethod
    def _cache_key(
        model: str,
        system: str,
        prompt: str,
        max_tokens: int,
        temperature: float,
        seed: int,
    ) -> str:
        blob = json.dumps(
            {
                "m": model,
                "s": system,
                "p": prompt,
                "t": max_tokens,
                "temperature": temperature,
                "seed": seed,
            },
            sort_keys=True,
        )
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()

    def _cache_path(self, key: str) -> Path:
        # 2-char shard prefix to keep any single directory shallow.
        return self._cache_dir / key[:2] / f"{key}.json"

    def _cache_get(self, key: str) -> Optional[str]:
        try:
            path = self._cache_path(key)
            if not path.exists():
                return None
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f).get("text")
        except (OSError, json.JSONDecodeError):
            return None

    def _cache_put(self, key: str, text: str) -> None:
        try:
            path = self._cache_path(key)
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(".tmp")
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump({"text": text}, f)
            tmp.replace(path)
        except OSError:
            pass  # cache is best-effort

    def _get_context_manager(self, cfg: ModelConfig) -> ContextManager:
        """Get or create a ContextManager for a model."""
        if cfg.model not in self._context_managers:
            profile = ModelProfile(
                name=cfg.model,
                context_window=cfg.context_window,
                is_local=cfg.is_local,
            )
            self._context_managers[cfg.model] = ContextManager(profile)
        return self._context_managers[cfg.model]

    def _inject_api_keys(self):
        """Set API keys from environment if not explicitly provided."""
        env_map = {
            "anthropic": "ANTHROPIC_API_KEY",
            "openai":    "OPENAI_API_KEY",
            "gemini":    "GEMINI_API_KEY",
        }
        for provider, env_var in env_map.items():
            if provider not in self.api_keys:
                val = os.environ.get(env_var)
                if val:
                    self.api_keys[provider] = val

    @staticmethod
    def _key_env(provider: str) -> str:
        return {
            "anthropic": "ANTHROPIC_API_KEY",
            "openai":    "OPENAI_API_KEY",
            "gemini":    "GEMINI_API_KEY",
        }.get(provider, f"{provider.upper()}_API_KEY")

    # ── Factory ──────────────────────────────────────────────────────────

    @classmethod
    def from_config(cls, config_path: str = "~/.aria/config.yaml") -> "LLMProvider":
        """
        Build LLMProvider from ARIA config file.
        Falls back to defaults + environment variables if config not found.

        Example config.yaml:
        ---
        llm:
          heavy:
            provider: anthropic
            model: claude-sonnet-4-20250514
          medium:
            provider: ollama
            model: llama3:70b
            api_base: http://localhost:11434
          light:
            provider: ollama
            model: llama3:8b
            api_base: http://localhost:11434
        """
        import yaml
        from pathlib import Path

        config_file = Path(config_path).expanduser()

        if not config_file.exists():
            log.info("No config.yaml found — using defaults + env vars")
            return cls()

        with open(config_file) as f:
            config = yaml.safe_load(f)

        llm_config = config.get("llm", {})
        models = {}

        tier_map = {
            "heavy":  TaskTier.HEAVY,
            "medium": TaskTier.MEDIUM,
            "light":  TaskTier.LIGHT,
        }

        for tier_name, tier_enum in tier_map.items():
            if tier_name in llm_config:
                entry = llm_config[tier_name]
                provider = entry.get("provider", "anthropic")
                model_id = entry.get("model", "")
                # Prefix with provider for LiteLLM if local
                if provider == "ollama" and not model_id.startswith("ollama/"):
                    model_id = f"ollama/{model_id}"
                cfg = ModelConfig(
                    provider=provider,
                    model=model_id,
                    context_window=entry.get("context_window", 8_000),
                    is_local=(provider in ("ollama", "vllm")),
                    api_base=entry.get("api_base"),
                )
                # Primary config + defaults as fallback
                models[tier_enum] = [cfg] + DEFAULT_MODELS.get(tier_enum, [])[1:]

        return cls(models=models)
