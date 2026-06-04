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

from aria.version import __version__ as ARIA_VERSION
from aria.llm.context_manager import ContextManager, ModelProfile
from aria.utils.env_loader import load_aria_env
from aria.utils.privacy import air_gapped_enabled
from aria.utils.provenance import record_llm_usage

log = logging.getLogger("aria.llm")

# Silence LiteLLM's aggressive logging
litellm.suppress_debug_info = True
logging.getLogger("LiteLLM").setLevel(logging.WARNING)
# Defensive: never crash a call because a provider rejects an optional param
# (e.g. anthropic/gemini reject `seed`). LiteLLM drops the unsupported param for
# that provider instead of raising litellm.UnsupportedParamsError.
litellm.drop_params = True

# Providers that actually honor an integer `seed` (deterministic sampling). For
# everyone else we omit it — temperature=0.0 remains the primary determinism
# control — so the call does not depend on litellm.drop_params and the recorded
# provenance is honest about whether a seed was applied.
_SEED_SUPPORTED_PROVIDERS = {"openai", "ollama", "vllm"}


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
        ModelConfig("anthropic", "claude-opus-4-8",            200_000),
        ModelConfig("anthropic", "claude-sonnet-4-6",          200_000),
        ModelConfig("openai",    "gpt-4o",                     128_000),
        ModelConfig("gemini",    "gemini/gemini-2.5-pro",      1_000_000),
        ModelConfig("ollama",    "ollama/llama3:70b",           8_000,  is_local=True,
                    api_base="http://localhost:11434"),
    ],
    TaskTier.MEDIUM: [
        ModelConfig("anthropic", "claude-sonnet-4-6",          200_000),
        ModelConfig("anthropic", "claude-haiku-4-5-20251001",  200_000),
        ModelConfig("openai",    "gpt-4o-mini",                128_000),
        # Cloud fallback for a Gemini-only user — must precede the local Ollama
        # entries, which are off unless an Ollama server is actually running.
        ModelConfig("gemini",    "gemini/gemini-2.5-flash",  1_000_000),
        ModelConfig("ollama",    "ollama/llama3:8b",            8_000,  is_local=True,
                    api_base="http://localhost:11434"),
        ModelConfig("ollama",    "ollama/mistral:7b",           8_000,  is_local=True,
                    api_base="http://localhost:11434"),
    ],
    TaskTier.LIGHT: [
        ModelConfig("anthropic", "claude-haiku-4-5-20251001",  200_000),
        ModelConfig("gemini",    "gemini/gemini-2.5-flash",  1_000_000),
        ModelConfig("ollama",    "ollama/llama3:8b",            8_000,  is_local=True,
                    api_base="http://localhost:11434"),
        ModelConfig("ollama",    "ollama/mistral:7b",           8_000,  is_local=True,
                    api_base="http://localhost:11434"),
    ],
}


def diagnose_llm_failure(error: Exception | str) -> str:
    """Translate an LLM-provider failure into an actionable hint (no traceback).

    ARIA reaches no configured model when the configured provider has no API key
    and the local Ollama fallback is not running. This inspects which API keys are
    actually present and tells the user the smallest fix.
    """
    keys = {
        "anthropic": bool(os.environ.get("ANTHROPIC_API_KEY")),
        "openai": bool(os.environ.get("OPENAI_API_KEY")),
        "gemini": bool(os.environ.get("GEMINI_API_KEY")
                       or os.environ.get("GOOGLE_API_KEY")),
    }
    present = [p for p, v in keys.items() if v]
    lines = [
        "ARIA could not reach any configured LLM model.",
        f"  Detail: {str(error)[:200]}",
    ]
    if present:
        lines += [
            f"  Detected API key(s) in your environment: {', '.join(present)}.",
            "  Your ~/.aria/config.yaml likely points at a provider you have no "
            "key for.",
            f"  Fix: set each llm tier's provider/model to one you have a key for, "
            f"e.g. provider: {present[0]}.",
        ]
    else:
        lines += [
            "  No ANTHROPIC / OPENAI / GEMINI API key was found in your "
            "environment,",
            "  and the local Ollama fallback is not running (connection refused).",
            "  Fix: export an API key (e.g. GEMINI_API_KEY=...) or start an "
            "Ollama server.",
        ]
    lines.append("  Diagnose: run `aria doctor --llm`.")
    return "\n".join(lines)


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
    # R3 (audit 2026-05-29): every LLM call gets a wall-clock timeout so a hung
    # provider cannot block the dispatch thread forever. Overridable via
    # ARIA_LLM_TIMEOUT (seconds). A hit timeout raises and the tier falls back.
    DEFAULT_TIMEOUT_S = 120.0

    def __init__(
        self,
        models:      dict[TaskTier, list[ModelConfig]] = None,
        api_keys:    dict[str, str] = None,
        cache_dir:   Optional[str] = None,
    ):
        load_aria_env()
        self._air_gapped = air_gapped_enabled()
        configured_models = models or DEFAULT_MODELS
        if self._air_gapped:
            configured_models = {
                tier: [m for m in cfgs if m.is_local]
                for tier, cfgs in configured_models.items()
            }
        self.models   = configured_models
        self.api_keys = api_keys or {}
        # R3: per-call timeout (seconds). Env override, clamped to a sane floor.
        try:
            self._timeout_s = max(
                5.0, float(os.environ.get("ARIA_LLM_TIMEOUT",
                                          self.DEFAULT_TIMEOUT_S))
            )
        except (TypeError, ValueError):
            self._timeout_s = self.DEFAULT_TIMEOUT_S
        self._context_managers: dict[str, ContextManager] = {}
        self._inject_api_keys()
        # File-backed prompt cache. Disabled when ARIA_LLM_CACHE=0.
        # Cache key: sha256(model + system + prompt + max_tokens + deterministic
        # generation controls + ARIA/cache salt). C6/X10: cache entries expire
        # after ARIA_LLM_CACHE_TTL_DAYS (default 30); set <=0 to disable expiry.
        self._cache_enabled = os.environ.get("ARIA_LLM_CACHE", "1") != "0"
        self._cache_version_salt = os.environ.get(
            "ARIA_LLM_CACHE_SALT", f"aria-{ARIA_VERSION}"
        )
        try:
            self._cache_ttl_s = int(
                float(os.environ.get("ARIA_LLM_CACHE_TTL_DAYS", "30"))
                * 86400
            )
        except (TypeError, ValueError):
            self._cache_ttl_s = 30 * 86400
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
        # P0-2: resolve lazily. `self.models.get(tier, self.models[MEDIUM])`
        # evaluated the default eagerly, so a config without a MEDIUM tier raised
        # KeyError on EVERY call — even for a tier that IS configured. Fall back
        # to MEDIUM only if needed, then fail with an explicit, actionable error.
        candidates = self.models.get(tier) or self.models.get(TaskTier.MEDIUM)
        if not candidates:
            raise RuntimeError(
                f"No model is configured for tier {tier.value} and no MEDIUM "
                f"fallback tier is available. Configure at least the MEDIUM tier "
                f"(or this tier) in config.yaml / LLMProvider(models=...)."
            )
        if self._air_gapped:
            candidates = [c for c in candidates if c.is_local]
            if not candidates:
                raise RuntimeError(
                    f"ARIA_AIR_GAPPED is enabled and no local model is "
                    f"configured for tier {tier.value}."
                )

        last_error = None
        for attempt, model_cfg in enumerate(candidates):
            try:
                # R4: when a later candidate answers, the call is a degraded
                # fallback — record which model it replaced and why so the
                # report can disclose that the section did not run on the
                # primary model.
                return self._call(
                    model_cfg, prompt, system, max_tokens, messages, tier=tier,
                    attempt=attempt,
                    fallback_from=(candidates[attempt - 1].model
                                   if attempt else None),
                    fallback_reason=(str(last_error)[:200]
                                     if last_error else None),
                )
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
        attempt:    int = 0,
        fallback_from: Optional[str] = None,
        fallback_reason: Optional[str] = None,
    ) -> str:
        """Execute a single LiteLLM call with context management."""
        if self._air_gapped and not cfg.is_local:
            raise RuntimeError(
                "ARIA_AIR_GAPPED is enabled; refusing cloud LLM call to "
                f"{cfg.provider}/{cfg.model}."
            )

        # R4: degradation provenance — stamped on every usage event (cache hit
        # or live) so the report can show whether a section ran on the primary
        # model or a fallback, and why the primary was skipped.
        degradation = {
            "is_fallback": bool(attempt),
            "fallback_attempt": int(attempt),
            "fallback_from": fallback_from,
            "fallback_reason": fallback_reason,
        }

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
                self._cache_version_salt,
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
                    **degradation,
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

        # seed only where the provider honors it (anthropic/gemini reject it).
        seed_applied = cfg.provider in _SEED_SUPPORTED_PROVIDERS
        seed_value = self.DETERMINISTIC_SEED if seed_applied else None
        kwargs = dict(
            model=cfg.model,
            messages=msgs,
            max_tokens=max_tokens,
            temperature=self.DETERMINISTIC_TEMPERATURE,
            timeout=self._timeout_s,   # R3: bound a hung provider
        )
        if seed_applied:
            kwargs["seed"] = self.DETERMINISTIC_SEED

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
            "seed": seed_value,
            "seed_applied": seed_applied,
            "deterministic": True,
            "cache_hit": False,
            "prompt_tokens": prompt_tokens,
            "completion_tokens": completion_tokens,
            "total_tokens": total_tokens,
            "estimated_cost_usd": cost,
            **degradation,
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
        version_salt: str = "",
    ) -> str:
        blob = json.dumps(
            {
                "m": model,
                "s": system,
                "p": prompt,
                "t": max_tokens,
                "temperature": temperature,
                "seed": seed,
                "version_salt": version_salt,
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
                payload = json.load(f)
            if self._cache_ttl_s > 0:
                created = float(payload.get("created_unix", 0) or 0)
                if created <= 0 or time.time() - created > self._cache_ttl_s:
                    try:
                        path.unlink()
                    except OSError:
                        pass
                    return None
            if payload.get("version_salt") not in {None, self._cache_version_salt}:
                return None
            return payload.get("text")
        except (OSError, json.JSONDecodeError):
            return None

    def _cache_put(self, key: str, text: str) -> None:
        try:
            path = self._cache_path(key)
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(".tmp")
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump({
                    "text": text,
                    "created_unix": time.time(),
                    "version_salt": self._cache_version_salt,
                }, f)
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
