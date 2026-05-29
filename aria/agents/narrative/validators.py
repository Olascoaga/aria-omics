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


def find_causal_language(text: str,
                         exclude: list[str] | None = None) -> str | None:
    """Return the first causal pattern found in ``text``, or None.

    Reusable so any layer (validators, renderer) can scan free text — not just
    the structured claim/evidence — for unlicensed causal phrasing.

    ``exclude`` lists external named entities (database term names, gene
    symbols) that must be removed before scanning. Their names are data ARIA is
    reporting, not claims ARIA is asserting (audit 2026-05-29, B5): a GO/Reactome
    term such as "TP53 Regulates Transcription of Cell Cycle Genes" contains a
    causal verb but is not a causal claim.
    """
    if not text:
        return None
    if exclude:
        for name in exclude:
            token = str(name).strip()
            if len(token) < 2:
                continue
            text = re.sub(rf"\b{re.escape(token)}\b", " ", text,
                          flags=re.IGNORECASE)
    match = _CAUSAL_RE.search(text)
    return match.group(0) if match else None


# Evidence whose label/value carries an EXTERNAL named entity — a database
# term/gene-set identifier or a gene symbol — rather than an assertion authored
# by ARIA. The causal guard must read ARIA's claims, not the names of the
# entities it reports (audit 2026-05-29, B5).
_NAMED_ENTITY_SOURCES = {
    "bulk_pathways",
    "pseudobulk_pathways",
    "gsea_preranked",
}


def _evidence_is_named_entity(ev) -> bool:
    if getattr(ev, "source", None) in _NAMED_ENTITY_SOURCES:
        return True
    low = str(getattr(ev, "label", "") or "").lower()
    return (
        low.startswith("gsea term ")
        or low.startswith("top gene ")
        or low.endswith(" term")
        or low.endswith(" fdr")
    )


def collect_named_entities(block: NarrativeBlock) -> list[str]:
    """External identifiers (DB term names, gene symbols) carried by a block.

    These must not be read as ARIA's causal assertions by the causal guard. The
    precise token is extracted so it can be redacted from composed prose before
    scanning.
    """
    names: list[str] = []
    for ev in block.evidence:
        label = str(ev.label or "")
        low = label.lower()
        if low.startswith("gsea term "):
            names.append(label[len("GSEA term "):])
        elif low.startswith("top gene "):
            names.append(label[len("top gene "):])
        elif low.endswith(" term") and ev.value is not None:
            names.append(str(ev.value))
        elif ev.source in _NAMED_ENTITY_SOURCES and ev.value is not None:
            names.append(str(ev.value))
    for row in block.metrics.get("top_pathways") or []:
        if isinstance(row, dict) and row.get("term"):
            names.append(str(row["term"]))
    return [n for n in names if str(n).strip()]

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
    # Scan ARIA's asserted claim plus evidence that is NOT an external named
    # entity. DB term names and gene symbols are data ARIA reports, not claims
    # it makes, so they must not trip the causal guard (audit 2026-05-29, B5).
    parts = [block.claim]
    for ev in block.evidence:
        if _evidence_is_named_entity(ev):
            continue
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
