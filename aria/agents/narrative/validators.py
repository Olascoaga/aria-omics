"""Integrity validators for structured narrative blocks."""

from __future__ import annotations

import re
from pathlib import Path

from aria.agents.narrative.types import Caveat, NarrativeBlock


class NarrativeValidationError(ValueError):
    pass


CAUSAL_PATTERNS = (
    # original set
    "drives",
    "enforces",
    "acts as",
    "governs",
    "position as",
    "positions as",
    "directly regulates",
    "binds to",
    "hierarchical gatekeeper",
    "hierarchical gatekeepers",
    # broadened lexicon (audit 2026-05-28, F-SCI-CAUSAL): associative
    # transcriptomic/communication evidence does not license causal verbs.
    "causes",
    "caused by",
    "induces",
    "induced by",
    "triggers",
    "leads to",
    "results in",
    "gives rise to",
    "controls",
    "regulates",
    "drives the",
    "promotes",
    "inhibits",
    "suppresses",
    "activates",
    "represses",
    "responsible for",
    "master regulator",
    "master regulators",
    "upstream of",
    "downstream effector",
    "orchestrates",
    "dictates",
    "determines",
    "mediates",
)

_CAUSAL_RE = re.compile(
    r"\b(" + "|".join(re.escape(p) for p in CAUSAL_PATTERNS) + r")\b",
    flags=re.IGNORECASE,
)


def find_causal_language(text: str) -> str | None:
    """Return the first causal pattern found in ``text``, or None.

    Reusable so any layer (validators, renderer) can scan free text — not just
    the structured claim/evidence — for unlicensed causal phrasing.
    """
    if not text:
        return None
    match = _CAUSAL_RE.search(text)
    return match.group(0) if match else None

_CONFIDENCE_DOWN = {
    "high": "medium",
    "medium": "low",
    "low": "insufficient",
    "insufficient": "insufficient",
}


def validate_blocks(blocks: list[NarrativeBlock],
                    base_dir: str | Path | None = None,
                    check_files: bool = False) -> list[NarrativeBlock]:
    return [
        validate_block(block, base_dir=base_dir, check_files=check_files)
        for block in blocks
    ]


def validate_block(block: NarrativeBlock,
                   base_dir: str | Path | None = None,
                   check_files: bool = False) -> NarrativeBlock:
    _validate_required_content(block)
    _apply_low_confidence_warning(block)
    _apply_causal_guard(block)
    _apply_trajectory_guard(block)
    if check_files:
        _validate_referenced_files(block, base_dir=base_dir)
    return block


def _validate_required_content(block: NarrativeBlock) -> None:
    if block.status == "success":
        if not block.claim.strip():
            raise NarrativeValidationError(
                f"{block.id}: successful block requires claim"
            )
        if not block.evidence:
            raise NarrativeValidationError(
                f"{block.id}: successful block requires evidence"
            )
    elif not block.claim.strip() and not block.caveats and not block.error:
        raise NarrativeValidationError(
            f"{block.id}: non-success block requires claim, caveat, or error"
        )


def _apply_low_confidence_warning(block: NarrativeBlock) -> None:
    if block.confidence not in {"low", "insufficient"}:
        return
    msg = (
        "LOW/INSUFFICIENT confidence finding must remain visible in the "
        "limitations section."
    )
    if msg not in block.warnings:
        block.warnings.append(msg)


def _text_for_causal_scan(block: NarrativeBlock) -> str:
    parts = [block.claim]
    for ev in block.evidence:
        parts.append(str(ev.label))
        if ev.value is not None:
            parts.append(str(ev.value))
    return " ".join(parts)


def _apply_causal_guard(block: NarrativeBlock) -> None:
    if block.metadata.get("causal_evidence") is True:
        return
    pattern = find_causal_language(_text_for_causal_scan(block))
    if not pattern:
        return
    caveat_text = (
        f"Causal language pattern '{pattern}' was detected without explicit "
        "causal evidence; interpret this block as associative."
    )
    if not any(caveat_text == caveat.text for caveat in block.caveats):
        block.caveats.append(Caveat(caveat_text, "warning"))
    block.confidence = _CONFIDENCE_DOWN.get(block.confidence, "insufficient")


def _apply_trajectory_guard(block: NarrativeBlock) -> None:
    identity = f"{block.id} {block.analysis}".lower()
    if "trajectory" not in identity and "paga" not in identity and "dpt" not in identity:
        return
    if block.metadata.get("velocity_computed") or block.metadata.get("time_course"):
        return
    caveat_text = (
        "PAGA/DPT is exploratory without RNA velocity or time-course support; "
        "do not interpret ordering as active differentiation."
    )
    if not any("PAGA/DPT is exploratory" in c.text for c in block.caveats):
        block.caveats.append(Caveat(caveat_text, "warning"))


def _validate_referenced_files(block: NarrativeBlock,
                               base_dir: str | Path | None = None) -> None:
    base = Path(base_dir) if base_dir else None
    for kind, records in (("table", block.tables), ("figure", block.figures)):
        for rec in records or []:
            raw = rec.get("path") if isinstance(rec, dict) else None
            if not raw:
                continue
            path = Path(raw)
            if not path.is_absolute() and base is not None:
                path = base / path
            if not path.exists():
                raise NarrativeValidationError(
                    f"{block.id}: referenced {kind} does not exist: {raw}"
                )
