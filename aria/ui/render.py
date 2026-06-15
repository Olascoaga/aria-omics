"""Pure Rich renderers for the ARIA cockpit (U1 + U5).

These functions turn a :class:`~aria.runtime.experiment_view.ExperimentSnapshot`
into Rich renderables. They are deliberately free of any Textual import so they
can be unit-tested in the standard env and reused by both the Textual cockpit
and (if ever needed) a plain Rich view. No I/O, no bus reads, no mutation.
"""

from __future__ import annotations

from typing import Optional

from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from aria.runtime.experiment_view import (
    ExperimentSnapshot, ExperimentHistoryView,
)

# Ordered pipeline stages shown in the cockpit timeline.
_STAGES = ["audit", "design", "dispatch", "report"]

# Where each phase sits on the timeline (so "done" lights the whole bar).
_PHASE_INDEX = {
    "audit": 0, "design": 1, "dispatch": 2, "report": 3, "done": 4,
}

_CONF_STYLE = {
    "HIGH": "green", "MEDIUM": "yellow", "LOW": "red", "INSUFFICIENT": "dim",
}

_LEDGER_STATUS_STYLE = {
    "ran": "green", "skipped": "yellow", "not_run": "red", "error": "red bold",
}

_READINESS_STATUS_STYLE = {"green": "green", "yellow": "yellow", "red": "red"}
_READINESS_MARK = {"green": "●", "yellow": "▲", "red": "■"}
_RESOURCE_STATUS_STYLE = {
    "ready": "green",
    "missing": "yellow",
    "blocked": "red",
    "pending": "cyan",
    "info": "dim",
}
_RESOURCE_MARK = {
    "ready": "✓",
    "missing": "!",
    "blocked": "■",
    "pending": "…",
    "info": "·",
}
_ARTIFACT_STATUS_STYLE = {
    "present": "green",
    "ok": "green",
    "missing": "yellow",
    "warn": "yellow",
    "violation": "red bold",
}
_ARTIFACT_MARK = {
    "present": "✓",
    "ok": "✓",
    "missing": "!",
    "warn": "▲",
    "violation": "■",
}


def _fmt_duration(seconds: float) -> str:
    seconds = int(max(0.0, seconds))
    if seconds < 60:
        return f"{seconds}s"
    if seconds < 3600:
        return f"{seconds // 60}m{seconds % 60:02d}s"
    return f"{seconds // 3600}h{(seconds % 3600) // 60:02d}m"


def render_run_header(snap: ExperimentSnapshot, *,
                      version: str = "",
                      data_dir: Optional[str] = None,
                      modalities: Optional[list[str]] = None,
                      organism: Optional[str] = None,
                      air_gapped: bool = False) -> Panel:
    """Top header: compact run identity + environment context."""
    t = Text()
    t.append("Experiment ", style="cyan")
    t.append(str(snap.experiment_id), style="bold")
    if version:
        t.append("  |  ARIA ", style="dim")
        t.append(f"v{version}", style="cyan")
    if modalities:
        t.append("  |  ", style="dim")
        t.append(", ".join(modalities), style="cyan")
    if organism:
        t.append("  |  ", style="dim")
        t.append(str(organism), style="dim")
    t.append("  |  elapsed ", style="dim")
    t.append(_fmt_duration(snap.elapsed_s), style="dim")
    t.append("  |  egress ", style="dim")
    t.append("blocked" if air_gapped else "available",
             style="yellow" if air_gapped else "dim")
    if data_dir:
        t.append("\nData ", style="cyan")
        t.append(str(data_dir), style="dim")
    return Panel(t, title="[bold]Run[/]", border_style="cyan", padding=(0, 1))


def render_status_banner(message: str, version: str = "",
                         *, error: bool = False) -> Panel:
    """Left-panel banner shown while the run is starting (front-door transition).

    Used between the intake screen and the first snapshot so the cockpit shows a
    clear state instead of an empty panel; ``error=True`` styles a start failure.
    """
    style = "red" if error else "cyan"
    t = Text()
    if version:
        t.append("ARIA  ", style="cyan")
        t.append(f"v{version}\n\n", style="dim")
    t.append(message, style=style)
    return Panel(t, title="[bold]Run[/]", border_style=style, padding=(0, 1))


def _progress_value(value: float) -> float:
    return max(0.0, min(1.0, float(value or 0.0)))


def _progress_style(value: float) -> str:
    pct = _progress_value(value)
    if pct >= 1.0:
        return "bold green"
    if pct >= 0.75:
        return "bold cyan"
    if pct >= 0.35:
        return "cyan"
    return "blue"


def _progress_bar(value: float, *, width: int = 24) -> Text:
    """Terminal-safe segmented progress bar with a restrained sci-fi feel."""
    pct = _progress_value(value)
    filled = int(round(pct * width))
    style = _progress_style(pct)
    t = Text()
    t.append("▕", style="dim cyan")
    if filled:
        t.append("▰" * filled, style=style)
    if filled < width:
        t.append("▱" * (width - filled), style="dim")
    t.append("▏", style="dim cyan")
    return t


def _agent_progress_line(sender: str, progress: float, text: str = "",
                         *, width: int = 18) -> Text:
    t = Text()
    t.append(f"{sender:<22.22}", style="cyan")
    t.append(" ")
    t.append_text(_progress_bar(progress, width=width))
    t.append(f" {int(round(_progress_value(progress) * 100)):3d}%",
             style=_progress_style(progress))
    if text:
        t.append(f"  {text[:54]}", style="dim")
    return t


def _event_sort_key(event) -> float:
    ts = getattr(event, "ts", None)
    if hasattr(ts, "timestamp"):
        return float(ts.timestamp())
    try:
        return float(ts or 0.0)
    except (TypeError, ValueError):
        return 0.0


def render_agent_progress(snap: ExperimentSnapshot, *, limit: int = 10) -> Panel:
    """Left panel: latest progress/status per sender.

    This restores the cockpit's "agents are working" affordance: the global
    pipeline still shows the phase, while this panel shows which agents have
    emitted status and where each one is in its own work.
    """
    table = Table(expand=True, show_edge=False, pad_edge=False)
    table.add_column("agent progress")

    events = list(snap.agent_progress or [])
    events.sort(key=_event_sort_key)
    if not events:
        waiting = Text()
        waiting.append_text(_agent_progress_line(
            "orchestrator", 0.0, "Waiting for status."))
        table.add_row(waiting)
        return Panel(table, title="[bold]Agents[/]", border_style="dim",
                     padding=(0, 1))

    hidden = max(0, len(events) - limit)
    for event in events[-limit:]:
        status = (event.text or "").strip() or "working"
        table.add_row(_agent_progress_line(event.sender, event.progress, status))
    if hidden:
        table.add_row(Text(f"… {hidden} older agent status rows hidden",
                           style="dim"))
    border = "green" if snap.done else "cyan"
    return Panel(table, title="[bold]Agents[/]", border_style=border,
                 padding=(0, 1))


def render_overview(snap: ExperimentSnapshot) -> Panel:
    """Default cockpit view: clear run state for non-console users."""
    table = Table(expand=True, show_edge=False, pad_edge=False)
    table.add_column("area", style="cyan", no_wrap=True)
    table.add_column("status")

    if snap.done:
        state = Text("Analysis complete", style="bold green")
        if snap.report_path:
            state.append(f"\nReport: {snap.report_path}", style="dim")
    elif snap.pending_checkpoint is not None:
        state = Text("Decision required", style="bold yellow")
        state.append(f"\nCP{snap.pending_checkpoint.number}: "
                     f"{snap.pending_checkpoint.title}", style="dim")
    elif snap.last_status is not None:
        state = Text(snap.last_status.text or "Working", style="cyan")
        state.append(f"\n{snap.last_status.sender}", style="dim")
    else:
        state = Text("Preparing run", style="cyan")
    table.add_row("Now", state)

    active = list(snap.agent_progress or [])
    active.sort(key=_event_sort_key)
    if active:
        latest = active[-4:]
        agents = Text()
        for event in latest:
            agents.append_text(_agent_progress_line(
                event.sender, event.progress, event.text, width=14))
            agents.append("\n")
        table.add_row("Agents", agents)
    else:
        table.add_row("Agents", Text("Waiting for agent status.", style="dim"))

    counts = {k: len(v) for k, v in snap.findings_by_confidence.items()}
    findings = Text()
    for conf in ("HIGH", "MEDIUM", "LOW", "INSUFFICIENT"):
        findings.append(f"{conf[:4]}:{counts.get(conf, 0)}  ",
                        style=_CONF_STYLE[conf])
    table.add_row("Findings", findings)

    if snap.artifacts:
        arts = Text(f"{len(snap.artifacts)} artifact(s) available")
        if snap.done:
            arts.append(" · open [a] Artifacts", style="dim")
    else:
        arts = Text("Report artifacts appear when the run finishes.", style="dim")
    table.add_row("Outputs", arts)

    if snap.pending_checkpoint is not None:
        next_step = Text("Open [d] Decisions or press an option number.",
                         style="yellow")
    elif snap.done:
        next_step = Text("Review Artifacts or press q to exit.", style="green")
    else:
        next_step = Text("ARIA is running. Review [g] Agents for details.",
                         style="dim")
    table.add_row("Next", next_step)

    border = "green" if snap.done else (
        "yellow" if snap.pending_checkpoint is not None else "cyan")
    return Panel(table, title="[bold]Overview[/]", border_style=border,
                 padding=(0, 1))


def render_timeline(snap: ExperimentSnapshot) -> Panel:
    """Center-top panel: pipeline stage progress."""
    here = _PHASE_INDEX.get(snap.phase, 0)
    t = Text()
    for i, stage in enumerate(_STAGES):
        if i < here:
            t.append(f" ✓ {stage} ", style="green")
        elif i == here:
            t.append(f" ● {stage} ", style="bold cyan")
        else:
            t.append(f" ○ {stage} ", style="dim")
        if i < len(_STAGES) - 1:
            t.append("→", style="dim")
    pct = int(round(snap.progress * 100))
    sub = Text(f"\nprogress {pct}%", style="dim")
    if snap.last_status and snap.last_status.text:
        sub.append(f"  ·  {snap.last_status.text}", style="dim")
    t.append(sub)
    return Panel(t, title="[bold]Pipeline[/]", border_style="cyan",
                 padding=(0, 1))


# Cockpit tab views and their keys. Overview is the default view.
_CENTER_MODES = [
    ("overview", "o"),
    ("agents", "g"),
    ("decisions", "d"),
    ("findings", "f"),
    ("artifacts", "a"),
    ("resources", "u"),
    ("ledger", "l"),
    ("readiness", "r"),
]


def render_mode_bar(active_mode: str) -> Text:
    """Center-mode indicator: which view is showing and the key for each.

    Discoverability for the cockpit's swappable center. The active view is
    highlighted; the toggle key is shown in brackets (findings is the default).
    """
    t = Text()
    for i, (mode, key) in enumerate(_CENTER_MODES):
        if i:
            t.append("  ", style="dim")
        label = f"[{key}] {mode}" if key else mode
        if mode == active_mode:
            t.append(f" {label} ", style="bold black on cyan")
        else:
            t.append(label, style="dim")
    return t


def render_findings(snap: ExperimentSnapshot, *, limit: int = 12) -> Panel:
    """Center panel: findings stream with confidence badges."""
    counts = {k: len(v) for k, v in snap.findings_by_confidence.items()}
    header = Text()
    for conf in ("HIGH", "MEDIUM", "LOW", "INSUFFICIENT"):
        header.append(f"{conf[:4]}:{counts.get(conf, 0)}  ",
                      style=_CONF_STYLE[conf])

    flat = [f for bucket in snap.findings_by_confidence.values()
            for f in bucket]
    flat.sort(key=lambda f: f.ts)

    body = Text()
    body.append(header)
    body.append("\n")
    if not flat:
        body.append("No findings yet.", style="dim")
    else:
        hidden = max(0, len(flat) - limit)
        if hidden:
            # Disclose, never silently drop: older findings are off-screen.
            body.append(f"… +{hidden} older ({len(flat)} total)\n", style="dim")
        for f in flat[-limit:]:
            body.append(f"[{f.confidence[:4]}] ",
                        style=_CONF_STYLE.get(f.confidence, "white"))
            body.append(f"{f.summary[:90]}\n")
    return Panel(body, title="[bold]Findings[/]", border_style="cyan",
                 padding=(0, 1))


def render_checkpoint(snap: ExperimentSnapshot) -> Panel:
    """Right panel: the pending checkpoint (decision point), or idle state."""
    cp = snap.pending_checkpoint
    if cp is None:
        if snap.done:
            msg = Text("Run complete.\n", style="green")
            if snap.report_path:
                msg.append("\nReport\n", style="cyan")
                msg.append(snap.report_path, style="dim")
            msg.append("\n\nPress q to exit.", style="dim")
            border = "green"
        else:
            msg = Text("No pending checkpoint.", style="dim")
            border = "dim"
        return Panel(msg, title="[bold]Checkpoint[/]", border_style=border,
                     padding=(0, 1))
    t = Text()
    t.append(f"CP{cp.number} · {cp.title}\n", style="bold yellow")
    if cp.question:
        t.append(f"{cp.question[:240]}\n\n", style="white")
    for i, opt in enumerate(cp.options, start=1):
        t.append(f"  [{i}] ", style="cyan")
        t.append(f"{opt}\n")
    if cp.number == 2.1 and (cp.context or {}).get("proposed_groups"):
        t.append("\n  [e] ", style="bold green")
        t.append("open tabular design editor", style="green")
    return Panel(t, title="[bold]Decision required[/]", border_style="yellow",
                 padding=(0, 1))


def render_ledger(snap: ExperimentSnapshot) -> Panel:
    """U5 — run ledger: planned-vs-run governance made visible."""
    table = Table(expand=True, show_edge=False, pad_edge=False)
    table.add_column("modality", style="cyan", no_wrap=True)
    table.add_column("analysis", no_wrap=True)
    table.add_column("status", no_wrap=True)
    table.add_column("planned", justify="center", no_wrap=True)
    table.add_column("divergence", no_wrap=True)

    n_div = 0
    for node in snap.ledger:
        if node.divergence:
            n_div += 1
        status_txt = Text(node.status,
                          style=_LEDGER_STATUS_STYLE.get(node.status, "white"))
        div_txt = (Text(f"⚠ {node.reason or 'planned, did not run'}",
                        style="red")
                   if node.divergence else Text("—", style="dim"))
        table.add_row(
            node.modality,
            node.label or node.analysis,
            status_txt,
            "✓" if node.planned else "·",
            div_txt,
        )

    if not snap.ledger:
        table.add_row("—", "no ledger yet", "—", "—", "—")

    border = "red" if n_div else "green"
    title = (f"[bold]Run ledger[/]  "
             f"[{'red' if n_div else 'green'}]{n_div} divergence(s)[/]")
    return Panel(table, title=title, border_style=border, padding=(0, 1))


def render_readiness(snap: ExperimentSnapshot) -> Panel:
    """U3 — per-modality readiness cards (green/yellow/red), registry-driven."""
    t = Text()
    if not snap.readiness:
        t.append("No readiness cards yet "
                 "(computed at the quality-audit checkpoint).", style="dim")
        return Panel(t, title="[bold]Modality readiness[/]",
                     border_style="dim", padding=(0, 1))

    worst = "green"
    for card in snap.readiness:
        style = _READINESS_STATUS_STYLE.get(card.status, "white")
        mark = _READINESS_MARK.get(card.status, "?")
        if card.status == "red" or (card.status == "yellow" and worst != "red"):
            worst = card.status
        t.append(f"{mark} {card.modality}", style=f"bold {style}")
        t.append(f"  [{card.validation_level}] ", style="dim")
        t.append(f"{card.dispatch_policy}\n", style=style)
        if card.reason:
            t.append(f"    {card.reason}\n", style="dim")
        for msg in card.findings[:3]:
            t.append(f"    · {msg}\n", style=style)
    border = _READINESS_STATUS_STYLE.get(worst, "green")
    return Panel(t, title="[bold]Modality readiness[/]", border_style=border,
                 padding=(0, 1))


def render_resources(snap: ExperimentSnapshot) -> Panel:
    """U4 — local resource center: envs, references, caches, and egress policy."""
    table = Table(expand=True, show_edge=False, pad_edge=False)
    table.add_column("resource", style="cyan", no_wrap=True)
    table.add_column("status", no_wrap=True)
    table.add_column("detail")
    table.add_column("path/action")

    worst = "ready"
    rank = {"ready": 0, "info": 0, "pending": 1, "missing": 2, "blocked": 3}
    for res in snap.resources:
        status = str(res.status)
        if rank.get(status, 0) > rank.get(worst, 0):
            worst = status
        style = _RESOURCE_STATUS_STYLE.get(status, "white")
        mark = _RESOURCE_MARK.get(status, "?")
        tail = res.path or res.action or "—"
        table.add_row(
            f"{res.category}: {res.name}",
            Text(f"{mark} {status}", style=style),
            res.detail,
            Text(str(tail), style="dim"),
        )

    if not snap.resources:
        table.add_row("resources", "—", "No resource snapshot yet.", "—")

    border = _RESOURCE_STATUS_STYLE.get(worst, "cyan")
    return Panel(table, title="[bold]Resources[/]", border_style=border,
                 padding=(0, 1))


def render_artifacts(snap: ExperimentSnapshot) -> Panel:
    """U6 — report artifact browser: report/figs/tables/methodology/claims."""
    table = Table(expand=True, show_edge=False, pad_edge=False)
    table.add_column("artifact", style="cyan", no_wrap=True)
    table.add_column("status", no_wrap=True)
    table.add_column("detail")

    worst = "present"
    rank = {"present": 0, "ok": 0, "missing": 1, "warn": 1, "violation": 2}
    n_violation = 0
    for art in snap.artifacts:
        status = str(art.status)
        if status == "violation":
            n_violation += 1
        if rank.get(status, 0) > rank.get(worst, 0):
            worst = status
        style = _ARTIFACT_STATUS_STYLE.get(status, "white")
        mark = _ARTIFACT_MARK.get(status, "?")
        table.add_row(
            f"{art.category}: {art.name}",
            Text(f"{mark} {status}", style=style),
            art.detail,
        )

    if not snap.artifacts:
        table.add_row("artifacts", "—",
                      "No artifacts yet (the report is written at the end of the run).")

    border = _ARTIFACT_STATUS_STYLE.get(worst, "cyan")
    title = "[bold]Artifacts[/]"
    if n_violation:
        title += f"  [red]{n_violation} claim violation(s)[/]"
    return Panel(table, title=title, border_style=border, padding=(0, 1))


def render_history(history: list[ExperimentHistoryView]) -> Panel:
    """U7 — resume/history: prior experiments and their on-disk resume points.

    Rendered as compact text (not a wide table) so it stays readable in the
    narrow intake sidebar. Each entry shows identity, modalities, decision count,
    and whether a report bundle exists on disk (the honest resume point).
    """
    if not history:
        return Panel(Text("No prior experiments yet.", style="dim"),
                     title="[bold]Resume / history[/]", border_style="dim",
                     padding=(0, 1))

    t = Text()
    n_resumable = 0
    for h in history:
        if h.has_report:
            n_resumable += 1
        mark, mark_style = ("✓", "green") if h.has_report else ("·", "dim")
        t.append(f"{mark} ", style=mark_style)
        t.append(f"{h.name}", style="bold cyan")
        t.append(f"  [{(h.experiment_id or '?')[:8]}]\n", style="dim")
        meta = h.organism
        if h.modalities:
            meta += " · " + ", ".join(h.modalities)
        meta += f" · {h.n_decisions} decision(s)"
        t.append(f"    {meta}\n", style="dim")
        resume = ("✓ report on disk" if h.has_report else "— no report on disk")
        when = (h.updated_at or "")[:10]
        t.append(f"    {resume}", style=mark_style)
        if when:
            t.append(f" · {when}", style="dim")
        t.append("\n")

    title = "[bold]Resume / history[/]"
    if n_resumable:
        title += f"  [green]{n_resumable} with report[/]"
    return Panel(t, title=title, border_style="cyan", padding=(0, 1))
