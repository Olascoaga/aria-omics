"""Presentation-agnostic read-model over the MessageBus (U0).

Both the interactive TUI (:mod:`aria.tui`) and the headless runner
(:mod:`aria.headless`) independently polled the process-global MessageBus and
re-derived the same run state: de-duplicating messages by id, classifying
STATUS/FINDING/ESCALATION, detecting completion (``narrative_agent`` progress
``>= 1.0``), and extracting the report path. A Textual cockpit (U1) would have
been a third copy of that logic.

This module turns the bus state (plus an optional :class:`ExperimentSession`
for the planned-vs-run ledger) into one immutable, serializable snapshot. It is
pure: it reads the bus, never writes, and imports no UI toolkit. Checkpoint
*resolution* stays with the orchestrator; this only *reads* pending checkpoints.

Normalization note (latent-bug fix): agents publish STATUS text under
``payload["message"]`` (``base_agent.publish_status``), while the legacy
consumers read ``payload["status"]``. The snapshot reads either key, so status
text and the ``"Report saved: <path>"`` completion line are surfaced correctly
regardless of which key carried them.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Optional

from aria.bus.message_bus import bus, Message, MessageType

# Canonical checkpoint titles. Single source of truth for any presentation
# layer (the legacy `aria.tui` keeps its own copy until U1 converges on this).
CHECKPOINT_TITLES: dict[int | float, str] = {
    1:   "Data Audit Results",
    2:   "Analysis Plan",
    2.1: "Experimental Groups",
    2.2: "Organism",
    2.3: "Experimental Factor",
    2.4: "Batch Effects",
    2.5: "Pseudoreplication Check",
    2.6: "Design Confirmation",
    3:   "Quality Control / Parameter Decision",
    4:   "Preliminary Findings",
    5:   "Final Report Ready",
}

# Pre-dispatch design checkpoints (CP1 audit, CP2.x design). Used for phase.
_AUDIT_CHECKPOINTS = {1}
_DESIGN_CHECKPOINTS = {2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6}

_REPORT_SAVED_RE = re.compile(r"Report saved:\s*(\S+)")

_CONFIDENCE_BUCKETS = ("HIGH", "MEDIUM", "LOW", "INSUFFICIENT")


def status_text(payload: Optional[dict]) -> str:
    """Normalize a STATUS payload's human-readable text.

    Agents write the text under ``message`` (``base_agent.publish_status``); the
    historical consumers read ``status``. Prefer an explicit ``status`` (legacy
    key) when present, otherwise fall back to ``message``.
    """
    if not payload:
        return ""
    return str(payload.get("status") or payload.get("message") or "")


# ── View dataclasses (immutable, serializable) ───────────────────────────────

@dataclass(frozen=True)
class ProgressEvent:
    ts: datetime
    sender: str
    text: str
    progress: float


@dataclass(frozen=True)
class FindingView:
    ts: datetime
    sender: str
    confidence: str          # HIGH | MEDIUM | LOW | INSUFFICIENT
    summary: str
    payload: dict = field(default_factory=dict)


@dataclass(frozen=True)
class CheckpointView:
    message_id: str
    number: int | float | str
    title: str
    question: str
    options: list[str]
    resolved: bool = False
    decision: Any = None
    # Structured escalation context (e.g. CP2.1 carries
    # {"proposed_groups": {group: [samples]}}) for richer editors like U2.
    context: dict = field(default_factory=dict)


@dataclass(frozen=True)
class LedgerNodeView:
    modality: str
    analysis: str
    label: str
    status: str              # ran | skipped | not_run | error
    planned: bool
    divergence: bool
    reason: Optional[str]
    node_id: str


@dataclass(frozen=True)
class ExperimentSnapshot:
    experiment_id: str
    phase: str               # audit | design | dispatch | report | done
    progress: float
    last_status: Optional[ProgressEvent]
    findings_by_confidence: dict[str, list[FindingView]]
    pending_checkpoint: Optional[CheckpointView]
    ledger: list[LedgerNodeView]
    report_path: Optional[str]
    done: bool
    elapsed_s: float
    silent_s: float


# ── Internal helpers ─────────────────────────────────────────────────────────

def _confidence_label(msg: Message) -> str:
    raw = msg.payload.get("confidence")
    if raw is None:
        raw = getattr(getattr(msg, "confidence", None), "value", None)
    label = str(raw or "medium").upper()
    return label if label in _CONFIDENCE_BUCKETS else "MEDIUM"


def _detect_completion(log: list[Message]) -> tuple[bool, Optional[str]]:
    """Return (done, report_path) from a run's STATUS history.

    Completion is the narrative agent reaching progress ``>= 1.0`` (the report
    is written), or an orchestrator "complete" at ``>= 1.0``. The report path is
    taken from a ``report_path`` payload key or the ``"Report saved: <path>"``
    status line — both now read correctly via :func:`status_text`.
    """
    done = False
    report_path: Optional[str] = None
    for m in log:
        payload = m.payload or {}
        if "report_path" in payload and payload["report_path"]:
            report_path = str(payload["report_path"])
        if m.type != MessageType.STATUS:
            continue
        progress = payload.get("progress") or 0
        try:
            progress = float(progress)
        except (TypeError, ValueError):
            progress = 0.0
        text = status_text(payload)
        if progress >= 1.0 and m.sender == "narrative_agent":
            done = True
            mm = _REPORT_SAVED_RE.search(text)
            if mm:
                report_path = mm.group(1)
        elif (progress >= 1.0 and m.sender == "orchestrator"
              and "complete" in text.lower()):
            done = True
    return done, report_path


def _phase(progress: float, done: bool,
           pending: Optional[CheckpointView]) -> str:
    if done:
        return "done"
    if pending is not None:
        n = pending.number
        if n in _AUDIT_CHECKPOINTS:
            return "audit"
        if n in _DESIGN_CHECKPOINTS:
            return "design"
        if n in (3, 4):
            return "dispatch"
        if n == 5:
            return "report"
    if progress < 0.1:
        return "audit"
    if progress < 0.7:
        return "dispatch"
    return "report"


def _ledger_views(session: Any) -> list[LedgerNodeView]:
    if session is None:
        return []
    exp_ctx = getattr(session, "exp_context", None) or {}
    agent_results = getattr(session, "agent_results", None) or {}
    if not exp_ctx and not agent_results:
        return []
    try:
        from aria.agents.narrative.run_ledger import build_run_ledger
        ledger = build_run_ledger(exp_ctx, agent_results)
    except Exception:
        return []
    out: list[LedgerNodeView] = []
    for e in ledger.get("entries", []) or []:
        if not isinstance(e, dict):
            continue
        out.append(LedgerNodeView(
            modality=str(e.get("modality", "")),
            analysis=str(e.get("analysis", "")),
            label=str(e.get("label", "")),
            status=str(e.get("status", "not_run")),
            planned=bool(e.get("planned")),
            divergence=bool(e.get("divergence")),
            reason=e.get("reason"),
            node_id=str(e.get("node_id", "")),
        ))
    return out


# ── Public API ───────────────────────────────────────────────────────────────

def build_snapshot(experiment_id: str, *,
                   session: Any = None,
                   start_time: Optional[datetime] = None,
                   now: Optional[datetime] = None) -> ExperimentSnapshot:
    """Build an immutable, UI-agnostic snapshot of one experiment's run state.

    Reads the process-global bus, scoped to ``experiment_id``. ``session`` (an
    :class:`ExperimentSession`) is optional and only used to reconcile the
    planned-vs-run ledger. ``start_time`` anchors ``elapsed_s``; ``now`` is
    injectable for deterministic tests.
    """
    now = now or datetime.now()
    log = bus.get_log(experiment_id)

    # Latest STATUS message -> progress + last_status banner.
    last_status: Optional[ProgressEvent] = None
    progress = 0.0
    last_msg_ts: Optional[datetime] = None
    for m in log:
        if last_msg_ts is None or m.timestamp > last_msg_ts:
            last_msg_ts = m.timestamp
        if m.type == MessageType.STATUS:
            payload = m.payload or {}
            p = payload.get("progress") or 0
            try:
                p = float(p)
            except (TypeError, ValueError):
                p = 0.0
            last_status = ProgressEvent(
                ts=m.timestamp, sender=m.sender,
                text=status_text(payload), progress=p,
            )
            progress = p

    # Findings grouped by confidence (served from the bus index).
    grouped: dict[str, list[FindingView]] = {b: [] for b in _CONFIDENCE_BUCKETS}
    for m in bus.get_findings(experiment_id):
        grouped[_confidence_label(m)].append(FindingView(
            ts=m.timestamp, sender=m.sender,
            confidence=_confidence_label(m),
            summary=str((m.payload or {}).get("summary", "")),
            payload=dict(m.payload or {}),
        ))

    # First unresolved checkpoint (both legacy loops act on pending[0]).
    pending: Optional[CheckpointView] = None
    pendings = [
        m for m in bus.get_pending_checkpoints(experiment_id=experiment_id)
        if not (m.payload or {}).get("resolved")
    ]
    if pendings:
        m = pendings[0]
        payload = m.payload or {}
        num = payload.get("checkpoint", m.checkpoint)
        pending = CheckpointView(
            message_id=m.id,
            number=num if num is not None else "?",
            title=CHECKPOINT_TITLES.get(num, f"Checkpoint {num}"),
            question=str(payload.get("question", "")),
            options=list(payload.get("options", ["Continue", "Cancel"])),
            resolved=bool(payload.get("resolved")),
            decision=payload.get("user_decision"),
            context=dict(payload.get("context") or {}),
        )

    done, report_path = _detect_completion(log)
    phase = _phase(progress, done, pending)

    elapsed_s = (now - start_time).total_seconds() if start_time else 0.0
    silent_s = (now - last_msg_ts).total_seconds() if last_msg_ts else 0.0

    return ExperimentSnapshot(
        experiment_id=experiment_id,
        phase=phase,
        progress=progress,
        last_status=last_status,
        findings_by_confidence=grouped,
        pending_checkpoint=pending,
        ledger=_ledger_views(session),
        report_path=report_path,
        done=done,
        elapsed_s=elapsed_s,
        silent_s=silent_s,
    )
