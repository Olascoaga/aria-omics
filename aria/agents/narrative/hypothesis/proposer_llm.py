"""LLM proposer for the SPECULATIVE tier (ADR-057 S5; rail #10 provenance).

This is the single place the LLM is given freedom: it reads ONLY the grounded
``EvidenceSignal`` list and proposes a competing set of hypotheses. It invents no
facts — the four publication gates (grounding, falsifiability, language,
devils_advocate) and the quarantine downstream filter anything it gets wrong, so
the model is free over the *connection* while the *facts* stay real.

The proposer is decoupled from the concrete LLM: it takes a ``complete`` callable
``(prompt, system) -> str`` so the agent core stays deterministic and the proposer
is trivially testable with a fake. ``LLMProposer.from_provider`` wires ARIA's
``LLMProvider`` (HEAVY tier, file-backed cache). Every hypothesis is stamped with
model provenance (model label, temperature, prompt + input-evidence hashes) so the
inputs are reproducible even though the generation is not.
"""

from __future__ import annotations

import hashlib
import json
import logging
from typing import Any, Callable

from .types import Hypothesis

log = logging.getLogger(__name__)

CompleteFn = Callable[[str, str], str]

_SYSTEM = (
    "You are a careful molecular-biology hypothesis generator. You are given ONLY "
    "a list of audited measurements (entities, directions, the analysis node each "
    "came from, and any confounds that analysis already flagged). Propose competing, "
    "falsifiable hypotheses that connect these measurements. RULES: (1) name ONLY "
    "entities present in the evidence — this applies to the mechanism prose too: "
    "do NOT name a gene/protein/TF in the mechanism that is not in the evidence "
    "list (a mechanism naming an un-measured entity is rejected); (2) cite ONLY "
    "the given audited_node_ref "
    "values in observation_refs; (3) the mechanism must be hedged speculation "
    "(may/could/suggests/we hypothesize), never an assertion of causation or finding; "
    "(4) every hypothesis needs a concrete discriminating experiment; (5) each "
    "hypothesis must offer a simpler/competing explanation and acknowledge every "
    "confound flagged on the evidence it uses. If the evidence does not support a "
    "defensible hypothesis, return an empty list."
)

_SCHEMA_HINT = (
    'Return a JSON array of objects with keys: '
    '"id" (short string), '
    '"mechanism" (hedged speculative sentence), '
    '"entities" (list of entity names taken from the evidence), '
    '"observation_refs" (list of audited_node_ref values from the evidence), '
    '"experiment" {"perturbation","readout","predicted_direction","refuting_outcome"}, '
    '"devils_advocate" {"simpler_explanation", "confounds" (list)}.'
)


def build_proposer_prompt(
    signals: list, exp_ctx: dict | None, n_hypotheses: int
) -> str:
    """Render the audited evidence into a proposer prompt."""
    question = str((exp_ctx or {}).get("biological_question") or "").strip()
    lines: list[str] = []
    if question:
        lines.append(f"Biological question: {question}")
    lines.append(
        f"Propose up to {n_hypotheses} competing, falsifiable hypotheses from "
        "the audited evidence below."
    )
    lines.append("Audited evidence:")
    for sig in signals:
        d = sig.to_dict() if hasattr(sig, "to_dict") else dict(sig)
        caveats = ", ".join(d.get("caveats_inherited") or []) or "none"
        value = d.get("value")
        value_str = f" value={value}" if value is not None else ""
        lines.append(
            f"- {d.get('entity')} [{d.get('entity_kind')}] "
            f"{d.get('measure')} {d.get('direction')}{value_str} "
            f"from {d.get('audited_node_ref')}; confounds: {caveats}"
        )
    lines.append("")
    lines.append(_SCHEMA_HINT)
    return "\n".join(lines)


def _strip_fences(raw: str) -> str:
    text = (raw or "").strip()
    if text.startswith("```"):
        parts = text.split("```")
        if len(parts) >= 2:
            text = parts[1]
        if text.startswith("json"):
            text = text[4:]
    return text.strip("` \n")


def parse_hypotheses(raw: str) -> list[Hypothesis]:
    """Lenient parse of an LLM response into Hypothesis objects.

    Skips malformed items rather than crashing — a bad generation yields fewer (or
    zero) candidates, which the gates then treat as honest-null. Never fabricates.
    """
    text = _strip_fences(raw)
    if not text:
        return []
    try:
        data = json.loads(text)
    except (ValueError, TypeError):
        return []
    if isinstance(data, dict):
        data = data.get("hypotheses") or data.get("results") or []
    if not isinstance(data, list):
        return []
    out: list[Hypothesis] = []
    for i, item in enumerate(data):
        if not isinstance(item, dict):
            continue
        item.setdefault("id", f"h{i + 1}")
        try:
            out.append(Hypothesis.from_dict(item))
        except (TypeError, ValueError, KeyError):
            continue
    return out


def classify_parse(raw: str, parsed: int) -> dict:
    """Classify why a proposer response did or did not yield usable hypotheses.

    ADR-057 rail #10 (C): an honest-null must never hide the fact that the model
    actually answered. A non-empty response that does not parse — almost always a
    JSON truncated by ``max_tokens`` — is a generation failure, NOT "no defensible
    hypothesis", and the agent/report must say so instead of staying mute.
    """
    text = _strip_fences(raw or "")
    diag = {"raw_chars": len(raw or ""), "parsed": int(parsed)}
    if not (raw or "").strip():
        diag["status"] = "empty_response"
        return diag
    try:
        data = json.loads(text)
    except (ValueError, TypeError):
        # Non-empty text that is not valid JSON: the dominant cause is a response
        # cut off mid-object by the token budget (an "Unterminated string").
        diag["status"] = "parse_error"
        diag["detail"] = (
            "LLM response was non-empty but not valid JSON "
            "(likely truncated by max_tokens)"
        )
        return diag
    if isinstance(data, dict):
        data = data.get("hypotheses") or data.get("results") or []
    if not isinstance(data, list) or not data:
        # Valid JSON, empty list: the model legitimately declined to speculate.
        diag["status"] = "no_candidates"
        return diag
    if parsed == 0:
        diag["status"] = "malformed_items"
        return diag
    diag["status"] = "ok"
    return diag


class LLMProposer:
    """A ``Proposer`` that turns audited evidence into candidate hypotheses via an LLM."""

    def __init__(
        self,
        complete: CompleteFn,
        *,
        n_hypotheses: int = 4,
        model_label: str = "aria-llm",
        temperature: float | None = None,
        max_tokens: int = 8192,
        max_tokens_cap: int = 32768,
        max_retries: int = 2,
        served_meta: Callable[[], dict | None] | None = None,
    ) -> None:
        self._complete = complete
        self._n = n_hypotheses
        self._model_label = model_label
        self._temperature = temperature
        # H6 (issue 7): adaptive token budget. Start here and escalate on a
        # truncation, up to the cap, so a large evidence set cannot silently
        # truncate the JSON into an honest-null.
        self._max_tokens = int(max_tokens)
        self._max_tokens_cap = int(max_tokens_cap)
        self._max_retries = int(max_retries)
        # H6 (issue 8): a callable returning the provider's REAL last-completion
        # metadata (served model, fallback, finish_reason). None for fakes.
        self._served_meta = served_meta
        # ADR-057 rail #10 (C): last parse diagnostic, so an honest-null can
        # distinguish "model declined" from "model answered, we failed to parse".
        self.last_diagnostics: dict | None = None

    @classmethod
    def from_provider(
        cls, llm: Any, *, n_hypotheses: int = 4, max_tokens: int = 8192, **kwargs
    ) -> "LLMProposer":
        """Wire ARIA's LLMProvider (HEAVY tier) into a proposer.

        ADR-057 rail #10 (A): the default ``max_tokens`` is 8192 and escalates on
        truncation (H6). A full set of competing hypotheses — each with a hedged
        mechanism AND a four-field discriminating experiment AND a devils_advocate
        — runs to several thousand tokens on realistic evidence; a smaller budget
        truncates the JSON mid-object, which parses to zero candidates and
        silently collapses to honest-null (the v4.7.0 bug). The model_label and
        served metadata are read from the provider so provenance reflects the REAL
        model (issue 8), not a static label.
        """
        from aria.llm.provider import TaskTier

        def _complete(prompt: str, system: str, max_tokens: int = max_tokens) -> str:
            return llm.complete(
                prompt=prompt,
                system=system,
                tier=TaskTier.HEAVY,
                max_tokens=max_tokens,
            )

        label = "aria-llm:heavy"
        try:
            active = llm.get_active_model(TaskTier.HEAVY)
            if active is not None:
                label = "/".join(
                    p for p in (
                        str(getattr(active, "provider", "") or ""),
                        str(getattr(active, "model", "") or ""),
                    ) if p
                ) or label
        except Exception:
            pass

        return cls(
            _complete,
            n_hypotheses=n_hypotheses,
            model_label=label,
            max_tokens=max_tokens,
            served_meta=lambda: getattr(llm, "_last_completion", None),
            **kwargs,
        )

    def _invoke(self, prompt: str, system: str, budget: int) -> str:
        """Call the completion fn with an explicit budget, tolerating fakes.

        Injected test fakes use the simple ``(prompt, system)`` signature; the
        real provider-backed callable accepts a ``max_tokens`` override so the
        budget can escalate on truncation.
        """
        try:
            return self._complete(prompt, system, max_tokens=budget)
        except TypeError:
            return self._complete(prompt, system)

    def _provenance(
        self, prompt: str, signals: list, diag: dict, served: dict | None
    ) -> dict:
        evidence = json.dumps(
            [s.to_dict() if hasattr(s, "to_dict") else dict(s) for s in signals],
            sort_keys=True,
            default=str,
        )
        prov = {
            "generator": "HypothesisAgent.LLMProposer",
            "model_label": self._model_label,
            "temperature": self._temperature,
            "prompt_sha256": hashlib.sha256(prompt.encode("utf-8")).hexdigest(),
            "input_evidence_sha256": hashlib.sha256(
                evidence.encode("utf-8")
            ).hexdigest(),
            # H6: real generation provenance.
            "max_tokens": diag.get("max_tokens", self._max_tokens),
            "retries": diag.get("retries", 0),
            "finish_reason": diag.get("finish_reason"),
        }
        if served:
            prov["served_provider"] = served.get("provider")
            prov["served_model"] = served.get("model")
            prov["is_fallback"] = served.get("is_fallback")
            prov["fallback_from"] = served.get("fallback_from")
        return prov

    def __call__(self, signals: list, exp_ctx: dict | None) -> list[Hypothesis]:
        signals = list(signals or [])
        if not signals:
            self.last_diagnostics = {"status": "no_evidence", "parsed": 0}
            return []
        prompt = build_proposer_prompt(signals, exp_ctx, self._n)

        # H6 (issue 7): generate, and retry with a doubled budget ONLY on a
        # truncation (parse_error), up to max_retries / the cap. A clean parse,
        # an honest empty list, or malformed items are not retried.
        budget = self._max_tokens
        attempts = 0
        candidates: list[Hypothesis] = []
        diag: dict = {}
        while True:
            raw = self._invoke(prompt, _SYSTEM, budget)
            candidates = parse_hypotheses(raw)
            diag = classify_parse(raw, len(candidates))
            if (
                diag.get("status") != "parse_error"
                or attempts >= self._max_retries
                or budget >= self._max_tokens_cap
            ):
                break
            attempts += 1
            budget = min(budget * 2, self._max_tokens_cap)

        diag["max_tokens"] = budget
        diag["retries"] = attempts
        # Inferred finish_reason; overridden by the provider's real value below.
        diag["finish_reason"] = (
            "length" if diag.get("status") == "parse_error" else "stop"
        )
        served = None
        if self._served_meta is not None:
            try:
                served = self._served_meta()
            except Exception:
                served = None
        if served:
            diag["served_model"] = served.get("model")
            diag["served_fallback"] = bool(served.get("is_fallback"))
            if served.get("finish_reason"):
                diag["finish_reason"] = served["finish_reason"]
        self.last_diagnostics = diag

        if diag.get("status") not in ("ok", "no_candidates"):
            # Visible, non-silent: the model answered but we got nothing usable.
            log.warning(
                "HypothesisAgent proposer yielded no usable hypotheses "
                "(status=%s, raw_chars=%s, retries=%s, max_tokens=%s). %s",
                diag.get("status"),
                diag.get("raw_chars"),
                diag.get("retries"),
                diag.get("max_tokens"),
                diag.get("detail", ""),
            )
        provenance = self._provenance(prompt, signals, diag, served)
        for hyp in candidates:
            # H9 (F1): provenance is ALWAYS code-stamped, never inherited from the
            # model's JSON (from_dict already drops a model-supplied provenance;
            # overwrite unconditionally so no future path can let one through).
            hyp.provenance = dict(provenance)
        return candidates
