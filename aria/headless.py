"""
ARIA headless runner
--------------------

Non-interactive driver for the orchestrator + design checkpoint state machine.
It mirrors ``aria.tui.run_analysis`` but replaces the interactive
``print_checkpoint`` prompts with a deterministic, fully-logged answer policy.

Two uses:

  1. Batch / unattended analysis ("run ARIA over N datasets without a TTY").
  2. The integration-test surface for the checkpoint control path, which is
     otherwise only reachable through the interactive TUI.

Design principles honored:

  * No hardcoded biology. The default policy follows ARIA's OWN inference
    surfaced in each checkpoint payload (accept inferred groups/organism,
    pick the proposed factor parsed from the question, recommended+optional
    plan). It never injects gene/condition/dataset names (ADR-011).
  * Every checkpoint question, its options, and the chosen answer are recorded
    so a headless run is as auditable as an interactive one.
"""

from __future__ import annotations

import re
import time
import uuid
from dataclasses import dataclass, field
from typing import Callable, Optional

from aria.bus.message_bus import bus


# ── Answer policy ─────────────────────────────────────────────────────────

# A policy maps (checkpoint_number, question_text, options) -> chosen option.
AnswerPolicy = Callable[[object, str, list], str]


def default_answer_policy(cp_num, question: str, options: list) -> str:
    """Accept ARIA's high-confidence inference / recommended path.

    The choice for each checkpoint follows what the agent itself proposed,
    parsed from the payload — not any dataset-specific knowledge.
    """
    q = question or ""
    opts = list(options or [])

    def first_matching(substr: str) -> Optional[str]:
        for o in opts:
            if substr.lower() in o.lower():
                return o
        return None

    # CP1 — confirm detected data
    if cp_num == 1:
        return first_matching("confirm") or (opts[0] if opts else "Confirm and continue")

    # CP2.1 — confirm inferred groups
    if cp_num == 2.1:
        return first_matching("confirm groups") or (opts[0] if opts else opts[0])

    # CP2.2 — organism: opts[0] is ARIA's pre-seeded inferred organism
    if cp_num == 2.2:
        return opts[0]

    # CP2.3 — factor: pick the option equal to ARIA's proposed factor
    if cp_num == 2.3:
        m = re.search(r"proposed '([^']+)'", q)
        if m:
            proposed = m.group(1).strip()
            for o in opts:
                if o.strip().lower() == proposed.lower():
                    return o
        non_cancel = [o for o in opts if "cancel" not in o.lower()]
        return non_cancel[0] if non_cancel else opts[0]

    # CP2.4 — batch covariate: do not add unless ARIA pre-seeded one
    if cp_num == 2.4:
        return first_matching("ignore batch") or first_matching("no") or opts[0]

    # CP2.5 — pseudoreplication: treat as independent biological replicates
    if cp_num == 2.5:
        return first_matching("independent biological replicates") or opts[0]

    # CP2.6 — confirm assembled design
    if cp_num == 2.6:
        return first_matching("proceed") or opts[0]

    # CP3.5 — proceed despite non-blocking audit warnings
    if cp_num == 3.5:
        non_cancel = [o for o in opts if "cancel" not in o.lower()]
        return non_cancel[0] if non_cancel else opts[0]

    # CP2 — analysis plan: recommended + optional supported analyses
    if cp_num == 2:
        return first_matching("optional") or first_matching("recommended") or opts[0]

    # CP3 / CP4 / CP5 and anything else: take the recommended (first) option
    return opts[0] if opts else "Continue"


@dataclass
class HeadlessResult:
    experiment_id: str
    status: str                       # "completed" | "cancelled" | "timeout" | "failed"
    report_path: Optional[str] = None
    decisions: list = field(default_factory=list)   # [(cp_num, chosen), ...]


# ── Checkpoint draining (shared by runner and tests) ────────────────────────

def drain_pending_checkpoints(orchestrator, experiment_id: str,
                              policy: AnswerPolicy,
                              log: Callable[[str], None] = lambda _m: None,
                              max_rounds: int = 30) -> tuple[str, list]:
    """Resolve every currently-pending checkpoint for ``experiment_id``.

    Returns (status, decisions) where status is "ok" | "none" | "cancelled".
    """
    decisions: list = []
    resolved_any = False
    for _ in range(max_rounds):
        # R6: scope the pending-checkpoint read to this experiment at the bus
        # level so a concurrent run sharing the global bus is never observed.
        pending = [
            m for m in bus.get_pending_checkpoints(experiment_id=experiment_id)
            if not m.payload.get("resolved")
        ]
        if not pending:
            break
        msg = pending[0]
        cp_num = msg.payload.get("checkpoint", "?")
        question = msg.payload.get("question", "")
        options = msg.payload.get("options", ["Continue", "Cancel"])
        choice = policy(cp_num, question, options)

        log(f"CHECKPOINT {cp_num} options={options} -> {choice!r}")
        result = orchestrator.on_checkpoint_resolved(
            message_id=msg.id,
            user_decision=choice,
            experiment_id=experiment_id,
        )
        decisions.append((cp_num, choice))
        resolved_any = True
        if result.get("status") == "cancelled":
            return "cancelled", decisions
        time.sleep(0.05)
    return ("ok" if resolved_any else "none"), decisions


def run_headless(data_dir: str, question: str,
                 policy: AnswerPolicy = default_answer_policy,
                 reproducible_mode: bool = False,
                 timeout: float = 3600.0,
                 log: Callable[[str], None] = print) -> HeadlessResult:
    """Run a full ARIA analysis without a TTY.

    Drives the orchestrator through audit, the design checkpoint state machine,
    and the live dispatch loop, resolving checkpoints via ``policy``. Returns a
    :class:`HeadlessResult` with the report path on success.
    """
    from aria.memory.memory import ARIAMemory
    from aria.agents.orchestrator_agent import OrchestratorAgent

    experiment_id = str(uuid.uuid4())[:12]
    memory = ARIAMemory()
    orch = OrchestratorAgent(memory)

    ctx = {
        "data_dir": str(data_dir),
        "user_question": question,
        "reproducible_mode": reproducible_mode,
    }

    started = orch.run(experiment_id, ctx)
    if started.get("status") != "started":
        memory.close()
        return HeadlessResult(experiment_id, "failed")

    orch.run_audit(experiment_id)

    status, decisions = drain_pending_checkpoints(orch, experiment_id, policy, log)
    if status == "cancelled":
        memory.close()
        return HeadlessResult(experiment_id, "cancelled", decisions=decisions)

    # Live loop: keep draining late checkpoints until narrative signals done.
    # NOTE: agent STATUS messages do not always carry experiment_id, so we keep
    # reading exp-less messages (the narrative "Report saved" completion signal
    # can arrive without one). R6: but we now DROP messages explicitly tagged to
    # a DIFFERENT experiment, so a concurrent run sharing the global bus cannot
    # leak its events here. This mirrors aria.tui._live_analysis_loop.
    seen: set = set()
    start = time.time()
    report_path: Optional[str] = None
    while True:
        for m in bus.get_log():
            if m.experiment_id and m.experiment_id != experiment_id:
                continue
            if m.id in seen:
                continue
            seen.add(m.id)
            payload = m.payload or {}
            pr = payload.get("progress", "")
            st = str(payload.get("status", ""))
            if m.sender == "narrative_agent" and isinstance(pr, (int, float)) and pr >= 1.0:
                mm = re.search(r"Report saved:\s*(\S+)", st)
                if mm:
                    report_path = mm.group(1)
        _s, more = drain_pending_checkpoints(orch, experiment_id, policy, log)
        decisions.extend(more)
        if report_path is None:
            # Filesystem fallback: the report dir is named with the experiment
            # id suffix and contains report.html once narrative finishes.
            report_path = _find_report_on_disk(experiment_id)
        if report_path:
            memory.close()
            return HeadlessResult(experiment_id, "completed", report_path, decisions)
        if time.time() - start > timeout:
            memory.close()
            return HeadlessResult(experiment_id, "timeout", decisions=decisions)
        time.sleep(0.5)


def _find_report_on_disk(experiment_id: str) -> Optional[str]:
    """Return the path to a finished report.html for this experiment, if any.

    Report directories are named ``aria_<ts>_<slug>_<suffix>`` where suffix is
    derived from the experiment id; matching on the id's last segment is enough
    to disambiguate within a session.
    """
    from pathlib import Path

    reports = Path.home() / ".aria" / "reports"
    if not reports.is_dir():
        return None
    suffix = experiment_id.split("-")[-1]
    for d in reports.iterdir():
        if d.is_dir() and d.name.endswith(suffix):
            html = d / "report.html"
            if html.exists():
                return str(html)
    return None
