"""Protocols for modality-specific narrative plugins."""

from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Any, Protocol

from aria.agents.narrative.types import Caveat, NarrativeBlock
from aria.utils.stats import finite_float

OK_STAGE_STATUSES = frozenset({"success", "ok", "done"})
_TAG_STRIP_RE = re.compile(r"<[^>]+>")


def clean_token(val: Any, default: str = "unspecified") -> str:
    """Sanitize dynamic identifiers for narrative blocks/methods (strip HTML & control chars)."""
    if val is None:
        return default
    txt = _TAG_STRIP_RE.sub("", str(val))
    txt = re.sub(r"[\r\n\t]+", " ", txt).strip()
    return txt or default


def finite_round(val: Any, digits: int = 4) -> float | None:
    """Coerce ``val`` to a finite float rounded to ``digits`` decimals, else ``None``."""
    f = finite_float(val)
    return round(f, digits) if f is not None else None


def as_int(val: Any) -> int | None:
    """Coerce ``val`` to ``int``, rejecting ``bool`` and non-finite floats (no OverflowError)."""
    if val is None or isinstance(val, bool):
        return None
    val_type = type(val)
    if val_type.__name__ in ("bool", "bool_") and val_type.__module__.startswith("numpy"):
        return None
    try:
        if isinstance(val, float) and not math.isfinite(val):
            return None
        return int(val)
    except (TypeError, ValueError, OverflowError):
        return None


def stage_state(stage: Any) -> str:
    """Fail-closed stage status classifier: ``'ok'``, ``'skipped'``, ``'error'``, or ``'absent'``."""
    if not isinstance(stage, dict) or not stage:
        return "absent"
    st = str(stage.get("status") or "").lower().strip()
    if st in {"error", "failed"}:
        return "error"
    if st in {"skipped", "not_run", "blocked"} or stage.get("ran") is False:
        return "skipped"
    if st in OK_STAGE_STATUSES or stage.get("ran") is True:
        return "ok"
    return "skipped"


def diagnostic_block(
    block_id: str,
    analysis: str,
    title: str,
    stage: dict[str, Any],
    state: str,
    *,
    modality: str,
    block_type: str = "limitation",
    skipped_status: str = "limitation",
    caveat_template: str | None = None,
    metadata: dict[str, Any] | None = None,
) -> NarrativeBlock:
    """Build a fail-closed diagnostic/limitation NarrativeBlock for skipped or errored stages.

    Unpacks nested ``details`` dicts cleanly and populates ``error`` only when
    ``status == "error"`` (never attaching a spurious error string to benign
    skipped/limitation blocks).
    """
    raw_details = stage.get("details")
    if isinstance(raw_details, dict):
        raw_details = (
            raw_details.get("details")
            or raw_details.get("reason")
            or raw_details.get("error_type")
            or raw_details.get("message")
            or stage.get("stage")
        )
    reason = clean_token(
        stage.get("reason")
        or stage.get("message")
        or raw_details
        or stage.get("error_type")
        or f"stage reported status={stage.get('status')}",
        default="stage did not execute",
    )
    status = "error" if state == "error" else skipped_status
    tmpl = caveat_template or (
        "{title} did not produce inferential output ({reason}); "
        "no unmeasured values were substituted."
    )
    return NarrativeBlock(
        id=block_id,
        modality=modality,
        analysis=analysis,
        block_type=block_type,
        title=title,
        status=status,
        confidence="insufficient",
        claim="",
        error=reason if status == "error" else None,
        caveats=[
            Caveat(
                text=tmpl.format(title=title, status=status, reason=reason),
                severity="warning" if status == "error" else "info",
            )
        ],
        metadata=metadata or {},
    )


def record_figure_error(
    container: Any,
    key: str,
    exc: BaseException | str,
) -> None:
    """Record a non-fatal figure or table staging error for W-WARN routing.

    Appends a formatted message to ``container['figures']['_errors']``, which
    ``warning_router.route_degradation_signals`` routes onto the modality's
    narrative blocks as a warning caveat without aborting report generation.
    """
    if not isinstance(container, dict):
        return
    figs = container.get("figures")
    if not isinstance(figs, dict):
        figs = {}
        container["figures"] = figs
    errs = figs.get("_errors")
    if not isinstance(errs, list):
        errs = []
        figs["_errors"] = errs
    if isinstance(exc, BaseException):
        detail = f"{type(exc).__name__}: {exc}"
    else:
        detail = str(exc)
    errs.append(f"{key}: {detail}")


class ModalityNarrator(Protocol):
    name: str

    def accepts(self, agent_name: str, agent_result: dict) -> bool:
        """Return True when this narrator can handle the agent result."""

    def collect(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[NarrativeBlock]:
        """Convert an agent result into structured narrative blocks."""

    def methods(self, agent_name: str, agent_result: dict,
                context: dict | None = None) -> list[str]:
        """Return copy-pasteable methods text for this modality."""

    def figures(self, agent_name: str, agent_result: dict,
                report_dir: Path | None = None) -> list[dict]:
        """Return figure descriptors, optionally staging files first."""

    def tables(self, agent_name: str, agent_result: dict,
               report_dir: Path | None = None) -> list[dict]:
        """Return table descriptors, optionally staging files first."""
