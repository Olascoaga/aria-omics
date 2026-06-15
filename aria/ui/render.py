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

from aria.runtime.experiment_view import ExperimentSnapshot

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
    """Left panel: run identity + environment context."""
    t = Text()
    t.append("Experiment  ", style="cyan")
    t.append(f"{snap.experiment_id}\n")
    if version:
        t.append("ARIA        ", style="cyan")
        t.append(f"v{version}\n", style="dim")
    if modalities:
        t.append("Modality    ", style="cyan")
        t.append(f"{', '.join(modalities)}\n")
    if organism:
        t.append("Organism    ", style="cyan")
        t.append(f"{organism}\n", style="dim")
    if data_dir:
        t.append("Data        ", style="cyan")
        t.append(f"{data_dir}\n", style="dim")
    t.append("Air-gapped  ", style="cyan")
    t.append("yes\n" if air_gapped else "no\n",
             style="yellow" if air_gapped else "dim")
    t.append("Elapsed     ", style="cyan")
    t.append(_fmt_duration(snap.elapsed_s), style="dim")
    return Panel(t, title="[bold]Run[/]", border_style="cyan", padding=(0, 1))


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
        msg = (Text("Run complete.", style="green") if snap.done
               else Text("No pending checkpoint.", style="dim"))
        return Panel(msg, title="[bold]Checkpoint[/]", border_style="dim",
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
