"""
ARIA Terminal Interface
-----------------------
Professional TUI built with Rich.

Key fix: background checkpoint polling loop.
When the Dispatcher runs agents in a background thread, those agents
publish checkpoints (Leiden resolution, Hi-C resolution, final report)
that must be intercepted and shown to the user in real time.

The old TUI processed checkpoints once after audit, then stopped.
This version runs a polling loop that stays alive until the background
thread finishes, surfacing every checkpoint as it arrives.

Flow:
  1. User submits question + data dir
  2. Orchestrator.run() → intent parsed
  3. DataAuditAgent → Checkpoint 1 (confirm data detected)
  4. Checkpoint 2 (confirm analysis plan)
  5. Dispatcher thread launched in background
  6. [POLLING LOOP] — stays alive, showing agent progress +
     surfacing checkpoints 3-5 as agents publish them
  7. Thread completes → final summary shown
"""

from __future__ import annotations
import os
import sys
import uuid
import time
import threading
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.table import Table
from rich.prompt import Prompt, Confirm
from rich.rule import Rule
from rich.align import Align
from rich.live import Live
from rich.spinner import Spinner
from rich import box

from aria.memory.memory import ARIAMemory
from aria.agents.orchestrator_agent import OrchestratorAgent
from aria.bus.message_bus import bus, MessageType, Message

# ── Theme ─────────────────────────────────────────────────────────────────────

C = {
    "bg":      "#0a0f1a",
    "surface": "#111827",
    "border":  "#1e3a5f",
    "cyan":    "#00d4ff",
    "teal":    "#0891b2",
    "green":   "#10b981",
    "amber":   "#f59e0b",
    "red":     "#ef4444",
    "dim":     "#4b5563",
    "text":    "#e2e8f0",
    "muted":   "#94a3b8",
    "accent":  "#818cf8",
}

AGENT_COLORS = {
    "orchestrator":      C['cyan'],
    "data_audit_agent":  C['teal'],
    "setup_agent":       "#fbbf24",
    "scrna_agent":       C['green'],
    "bulk_rna_agent":    "#34d399",
    "rna_agent":         C['green'],
    "chromatin_agent":   C['accent'],
    "genome_arch_agent": "#a78bfa",
    "integration_agent": "#f472b6",
    "narrative_agent":   C['amber'],
    "debate_council":    C['red'],
}

AGENT_ICONS = {
    "orchestrator":      "(orch )",
    "data_audit_agent":  "(audit)",
    "setup_agent":       "(setup)",
    "scrna_agent":       "(scRNA)",
    "bulk_rna_agent":    "(bulk )",
    "rna_agent":         "(rna  )",
    "chromatin_agent":   "(chrom)",
    "genome_arch_agent": "( 3D  )",
    "integration_agent": "(integ)",
    "narrative_agent":   "(narr )",
    "debate_council":    "(deb  )",
}

# Human-readable names for display (paired with icons)
AGENT_NAMES = {
    "orchestrator":      "Orchestrator",
    "data_audit_agent":  "Data Audit",
    "setup_agent":       "Setup",
    "scrna_agent":       "scRNA",
    "bulk_rna_agent":    "Bulk RNA",
    "rna_agent":         "RNA",
    "chromatin_agent":   "Chromatin",
    "genome_arch_agent": "3D Genome",
    "integration_agent": "Integration",
    "narrative_agent":   "Report",
    "debate_council":    "Debate",
}

CHECKPOINT_TITLES = {
    1: "Data Audit Results",
    2: "Analysis Plan",
    3: "Quality Control / Parameter Decision",
    4: "Preliminary Findings",
    5: "Final Report Ready",
}

console = Console(highlight=False)

ARIA_BANNER = r"""
   ########    ########  ###    ########
  ##      ## ##      ## ###   ##      ##
  ##      ## ##      ## ###   ##      ##
  ##      ## ######### ###    ##      ##
 ########### ##      ## ###  ###########
  ##      ## ##      ## ###   ##      ##
  ##      ## ##      ## ###   ##      ##
  ##         ##      ## ##    ##
"""
TAGLINE = "Agentic Research Intelligence for -omics Analysis"
VERSION = "v0.2.0-alpha"


# ── Display helpers ───────────────────────────────────────────────────────────

def print_banner():
    console.print()
    console.print(Align.center(Text(ARIA_BANNER, style=f"bold {C['cyan']}")))
    console.print(Align.center(Text(TAGLINE, style=f"italic {C['muted']}")))
    console.print(Align.center(Text(f"  {VERSION}", style=C['dim'])))
    console.print()
    console.print(Rule(style=C['border']))
    console.print()


def print_section(title: str, color: str = None):
    console.print()
    console.print(Rule(
        f"[bold {color or C['teal']}]{title}[/]",
        style=C['border']
    ))
    console.print()


def print_agent_message(agent: str, message: str):
    color  = AGENT_COLORS.get(agent, C['muted'])
    icon   = AGENT_ICONS.get(agent, "[-]")
    prefix = Text(f" {icon} [{agent}] ", style=f"bold {color}")
    msg    = Text(message, style=C['text'])
    console.print(prefix + msg)


def print_checkpoint(number: int, title: str,
                     content: str, options: list[str]) -> str:
    """
    Display a checkpoint and block until user responds.
    Called from the main thread — safe for Rich prompts.
    """
    console.print()
    console.print(Align.center(
        Text(f"  CHECKPOINT {number}  ",
             style=f"bold white on {C['teal']}")
    ))
    console.print()
    console.print(Panel(
        content,
        border_style=C['teal'],
        padding=(1, 2),
        title=f"[bold {C['cyan']}]{title}[/]",
        title_align="left",
    ))
    console.print()
    console.print(f"  [{C['muted']}]Choose an option:[/]")
    console.print()
    for i, opt in enumerate(options, 1):
        style = C['green'] if i == 1 else C['amber'] if i == 2 else C['red']
        console.print(f"    [{style}][{i}][/] {opt}")
    console.print()

    while True:
        choice = Prompt.ask(
            f"  [bold {C['cyan']}]Enter choice[/]",
            choices=[str(i) for i in range(1, len(options) + 1)],
            console=console,
        )
        return options[int(choice) - 1]


def print_agent_progress(msg: Message):
    """Display a STATUS message from an agent with clear identification."""
    agent    = msg.sender
    progress = msg.payload.get("progress", 0)
    status   = msg.payload.get("status", "")
    pct      = f"{int(progress * 100):3d}%"
    color    = AGENT_COLORS.get(agent, C['muted'])
    icon     = AGENT_ICONS.get(agent, f"({agent[:5]:<5})")
    name     = AGENT_NAMES.get(agent, agent.replace("_agent", "").title())
    bar_len  = 16
    filled   = int(bar_len * progress)
    bar      = "█" * filled + "░" * (bar_len - filled)
    # Format: [icon] Name        [bar]  PCT  status text
    console.print(
        f"  [{color}]{icon:<8}[/] "
        f"[{color}]{name:<12}[/] "
        f"[{color}]{bar}[/] "
        f"[{C['dim']}]{pct}[/]  "
        f"[{C['muted']}]{status}[/]"
    )


def print_finding(msg: Message):
    """Display a FINDING from an agent with confidence badge."""
    conf    = str(msg.payload.get("confidence", "medium")).upper()
    summary = msg.payload.get("summary", "")[:120]
    agent   = msg.sender

    conf_colors = {
        "HIGH":         C['green'],
        "MEDIUM":       C['amber'],
        "LOW":          C['red'],
        "INSUFFICIENT": C['dim'],
    }
    badge_color = conf_colors.get(conf, C['muted'])
    icon        = AGENT_ICONS.get(agent, "[-]")

    console.print(
        f"  [{C['dim']}]{icon}[/] "
        f"[bold {badge_color}][{conf}][/] "
        f"[{C['text']}]{summary}[/]"
    )


# ── Input helpers ─────────────────────────────────────────────────────────────

def select_data_directory() -> Path:
    print_section("Data Directory", C['cyan'])
    console.print(
        f"  [{C['muted']}]Enter the path to your raw data directory.[/]\n"
        f"  [{C['dim']}]Example: /data/my_experiment  or  ~/lab/project1[/]\n"
    )
    while True:
        raw  = Prompt.ask(f"  [bold {C['cyan']}]Data path[/]", console=console)
        path = Path(raw).expanduser().resolve()
        if path.exists() and path.is_dir():
            n = sum(1 for f in path.rglob("*") if f.is_file())
            console.print()
            console.print(Panel(
                f"  [{C['green']}]✓[/] Directory found\n"
                f"  [{C['muted']}]Files detected: {n}[/]\n"
                f"  [{C['dim']}]{path}[/]",
                border_style=C['green'], padding=(0, 1),
            ))
            return path
        console.print(f"\n  [{C['red']}]Directory not found: {path}[/]\n")


def ask_biological_question() -> str:
    print_section("Biological Question", C['cyan'])
    console.print(
        f"  [{C['muted']}]What would you like to know about your data?[/]\n"
        f"  [{C['dim']}]Examples:[/]\n"
        f"  [{C['dim']}]  * What TFs are differentially active in condition A vs B?[/]\n"
        f"  [{C['dim']}]  * Which genes show coordinated RNA and chromatin changes?[/]\n"
        f"  [{C['dim']}]  * What cell types are present and how do they differ?[/]\n"
    )
    return Prompt.ask(
        f"\n  [bold {C['cyan']}]Your question[/]", console=console
    )


def show_existing_experiments(memory: ARIAMemory):
    wings = memory.list_wings()
    if not wings:
        return
    print_section("Previous Experiments", C['dim'])
    table = Table(box=box.SIMPLE, border_style=C['border'],
                  header_style=f"bold {C['teal']}", show_header=True)
    table.add_column("ID",       style=C['dim'],   width=10)
    table.add_column("Name",     style=C['text'],  width=25)
    table.add_column("Organism", style=C['muted'], width=20)
    table.add_column("Genome",   style=C['muted'], width=10)
    table.add_column("Updated",  style=C['dim'],   width=20)
    for w in wings[-5:]:
        table.add_row(
            w["id"][:8], w["name"],
            w.get("organism", "?"), w.get("genome", "?"),
            w.get("updated_at", "")[:16],
        )
    console.print(table)


# ── Core analysis loop ────────────────────────────────────────────────────────

def run_analysis(orchestrator: OrchestratorAgent,
                 experiment_id: str,
                 context: dict):
    """
    Run the full analysis with live checkpoint polling.

    This is the fix for the background thread checkpoint gap.
    After Checkpoint 2 launches the dispatcher thread, we enter
    a polling loop that:
      - Shows agent progress messages in real time
      - Intercepts pending checkpoints and blocks the main thread
        to collect user input (Checkpoints 3-5)
      - Exits when the dispatcher thread finishes

    The dispatcher thread and the main thread communicate exclusively
    through the MessageBus — no shared state, no locks needed.
    """
    # Start orchestrator
    result = orchestrator.run(experiment_id, context)
    if result.get("status") != "started":
        console.print(f"[{C['red']}]Failed to start analysis.[/]")
        return

    print_agent_message("orchestrator",
        f"Intent: {result['intent'].get('summary', '')}")

    # Phase 1: Data audit (synchronous)
    print_section("Scanning Data", C['teal'])
    with console.status(
        f"[{C['teal']}]DataAuditAgent scanning...[/]", spinner="dots"
    ):
        orchestrator.run_audit(experiment_id)

    # Phase 2: Process checkpoints 1 and 2 (pre-dispatch)
    _drain_checkpoints(orchestrator, experiment_id)

    # After CP2, dispatcher thread is running in background.
    # Phase 3: Live polling loop — stays alive until thread done.
    _live_analysis_loop(orchestrator, experiment_id)


def _drain_checkpoints(orchestrator: OrchestratorAgent,
                        experiment_id: str,
                        max_rounds: int = 10):
    """
    Process all currently pending checkpoints synchronously.
    Used for pre-dispatch checkpoints (1 and 2).
    """
    for _ in range(max_rounds):
        pending = [
            m for m in bus.get_pending_checkpoints()
            if m.experiment_id == experiment_id
            and not m.payload.get("resolved")
        ]
        if not pending:
            break

        msg    = pending[0]
        cp_num = msg.payload.get("checkpoint", "?")
        title  = CHECKPOINT_TITLES.get(cp_num, f"Checkpoint {cp_num}")

        choice = print_checkpoint(
            number=cp_num,
            title=title,
            content=msg.payload.get("question", ""),
            options=msg.payload.get("options", ["Continue", "Cancel"]),
        )

        result = orchestrator.on_checkpoint_resolved(
            message_id=msg.id,
            user_decision=choice,
            experiment_id=experiment_id,
        )

        if result.get("status") == "cancelled":
            console.print(f"\n  [{C['amber']}]Analysis cancelled.[/]\n")
            return "cancelled"

        print_agent_message(
            "orchestrator", f"Checkpoint {cp_num} resolved: {choice}"
        )
        time.sleep(0.2)

    return "ok"


def _live_analysis_loop(orchestrator: OrchestratorAgent,
                         experiment_id: str,
                         poll_interval: float = 0.5,
                         timeout: float = 43200.0,
                         idle_warning_after: float = 1800.0):
    """
    The polling loop that keeps the TUI alive during background dispatch.

    Runs until:
      a) A STATUS message from narrative_agent with progress=1.0 arrives
      b) Checkpoint 5 (final report) is processed
      c) Timeout reached (default 12 h)

    Idle handling:
      - If no new messages for `idle_warning_after` seconds (default 30 min),
        print a gentle reminder that the analysis is still working.
      - STAR alignment on a full human bulk RNA-seq dataset routinely takes
        2-3 hours per sample sequentially, so we expect long silences.

    On each tick:
      1. Check for new STATUS messages → show progress
      2. Check for new FINDING messages → show with confidence badge
      3. Check for pending CHECKPOINTS → block for user input
    """
    print_section("Analysis Running", C['green'])
    console.print(
        f"  [{C['muted']}]Agents are running. "
        f"Checkpoints will appear here as they arise.[/]\n"
    )

    seen_message_ids: set = set()
    start_time = time.time()
    last_message_time = start_time
    last_heartbeat_time = start_time
    idle_warning_shown = False
    analysis_done = False

    # Known long-running stages (published status before the silence)
    # that justify suppressing the "something is wrong" anxiety
    LONG_STAGES = ("aligning to genome", "star", "building", "star index",
                    "counting reads", "featurecounts", "differential expression")
    last_status_text = ""

    while not analysis_done:
        elapsed = time.time() - start_time
        silent  = time.time() - last_message_time

        if elapsed > timeout:
            console.print(
                f"\n  [{C['amber']}]TUI timeout reached ({timeout/3600:.0f}h). "
                f"The background pipeline may still be running.[/]\n"
                f"  [{C['muted']}]Check: ps aux | grep -E 'STAR|fastp|featureCounts'[/]\n"
                f"  [{C['muted']}]Report will be saved to ~/.aria/reports/ "
                f"when complete.[/]\n"
            )
            break

        # Idle reminder (once, when silence exceeds threshold)
        if silent > idle_warning_after and not idle_warning_shown:
            is_long_stage = any(s in last_status_text.lower()
                                 for s in LONG_STAGES)
            hint = (
                "This is expected during STAR alignment (silent stage). "
                if is_long_stage else ""
            )
            console.print(
                f"\n  [{C['muted']}]No updates for {silent/60:.0f} min. "
                f"{hint}Still running in background.[/]\n"
            )
            idle_warning_shown = True

        # Heartbeat every 5 min during extended silence (keeps you informed)
        if silent > 300 and (time.time() - last_heartbeat_time) > 300:
            stage_label = last_status_text[:60] if last_status_text else "working"
            console.print(
                f"  [{C['dim']}]● heartbeat: {stage_label} "
                f"({silent/60:.0f} min silent, {elapsed/60:.0f} min total)[/]"
            )
            last_heartbeat_time = time.time()

        # ── Collect new messages from MessageBus ──────────────────────────
        all_msgs = bus.get_log()
        new_msgs  = [
            m for m in all_msgs
            if m.experiment_id == experiment_id
            and m.id not in seen_message_ids
        ]

        for msg in new_msgs:
            seen_message_ids.add(msg.id)
            last_message_time = time.time()     # reset idle counter
            idle_warning_shown = False

            if msg.type == MessageType.STATUS:
                print_agent_progress(msg)
                last_status_text = str(msg.payload.get("status", ""))
                # CRITICAL: only the NarrativeAgent reaching 1.0 signals
                # the end of the pipeline. SetupAgent, BulkRNAAgent, etc.
                # also hit 1.0 when they finish their own work — but the
                # overall analysis isn't done until the report is written.
                if (msg.sender == "narrative_agent"
                        and msg.payload.get("progress", 0) >= 1.0):
                    analysis_done = True
                # Also accept an explicit "all done" from the orchestrator
                elif (msg.sender == "orchestrator"
                        and msg.payload.get("progress", 0) >= 1.0
                        and "complete" in str(msg.payload.get("status","")).lower()):
                    analysis_done = True

            elif msg.type == MessageType.FINDING:
                print_finding(msg)

        # ── Check for pending checkpoints ─────────────────────────────────
        pending = [
            m for m in bus.get_pending_checkpoints()
            if m.experiment_id == experiment_id
            and not m.payload.get("resolved")
            and m.id not in seen_message_ids
        ]

        for msg in pending:
            seen_message_ids.add(msg.id)
            cp_num = msg.payload.get("checkpoint", "?")
            title  = CHECKPOINT_TITLES.get(cp_num, f"Checkpoint {cp_num}")

            console.print()  # Spacing before checkpoint
            choice = print_checkpoint(
                number=cp_num,
                title=title,
                content=msg.payload.get("question", ""),
                options=msg.payload.get("options", ["Continue", "Cancel"]),
            )

            result = orchestrator.on_checkpoint_resolved(
                message_id=msg.id,
                user_decision=choice,
                experiment_id=experiment_id,
            )

            if result.get("status") == "cancelled":
                console.print(
                    f"\n  [{C['amber']}]Analysis cancelled.[/]\n"
                )
                return

            print_agent_message(
                "orchestrator", f"Checkpoint {cp_num} resolved: {choice}"
            )

            # Checkpoint 5 = final report — analysis complete
            if cp_num == 5:
                analysis_done = True

        time.sleep(poll_interval)

    # ── Final summary ──────────────────────────────────────────────────────
    _print_final_summary(experiment_id)


def _print_final_summary(experiment_id: str):
    """Show a summary of all findings after analysis completes."""
    print_section("Analysis Complete", C['green'])

    findings = bus.get_findings(experiment_id)
    if not findings:
        console.print(f"  [{C['muted']}]No findings recorded.[/]")
        return

    # Group by confidence
    by_conf: dict[str, list] = {
        "HIGH": [], "MEDIUM": [], "LOW": [], "INSUFFICIENT": []
    }
    for f in findings:
        conf = str(f.payload.get("confidence", "MEDIUM")).upper()
        by_conf.setdefault(conf, []).append(f)

    conf_styles = {
        "HIGH":         (C['green'],  "●"),
        "MEDIUM":       (C['amber'],  "◐"),
        "LOW":          (C['red'],    "○"),
        "INSUFFICIENT": (C['dim'],    "·"),
    }

    console.print(Panel(
        "\n".join([
            f"  [{color}]{symbol} {conf}: "
            f"{len(by_conf.get(conf, []))} findings[/]"
            for conf, (color, symbol) in conf_styles.items()
            if by_conf.get(conf)
        ]),
        title=f"[bold {C['cyan']}]Findings Summary[/]",
        border_style=C['teal'],
        padding=(0, 2),
    ))

    # Show HIGH confidence findings
    highs = by_conf.get("HIGH", [])
    if highs:
        console.print()
        console.print(f"  [bold {C['green']}]Key findings:[/]")
        for f in highs[:5]:
            summary = f.payload.get("summary", "")[:100]
            console.print(f"  [{C['dim']}]●[/] [{C['text']}]{summary}[/]")

    # Report path
    report_msgs = [
        m for m in bus.get_log()
        if m.experiment_id == experiment_id
        and "report_path" in m.payload
    ]
    if report_msgs:
        report_path = report_msgs[-1].payload["report_path"]
        console.print()
        console.print(Panel(
            f"  [{C['cyan']}]Report:[/] [{C['dim']}]{report_path}[/]\n"
            f"  [{C['muted']}]Open in browser to view the full analysis.[/]",
            border_style=C['green'],
            padding=(0, 1),
        ))


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    memory       = ARIAMemory()
    orchestrator = OrchestratorAgent(memory)

    os.system("clear" if os.name != "nt" else "cls")
    print_banner()

    startup_ctx = memory.startup_context()
    if "No experiments" not in startup_ctx:
        console.print(Panel(
            Text(startup_ctx, style=C['muted']),
            border_style=C['border'],
            title=f"[{C['dim']}]Memory[/]",
            title_align="left",
            padding=(0, 1),
        ))
        console.print()

    show_existing_experiments(memory)
    print_section("New Analysis", C['cyan'])

    action = Prompt.ask(
        f"  [{C['cyan']}]Action[/]",
        choices=["new", "exit"],
        default="new",
        console=console,
    )
    if action == "exit":
        console.print(f"\n  [{C['muted']}]Goodbye.[/]\n")
        return

    data_dir      = select_data_directory()
    question      = ask_biological_question()
    experiment_id = str(uuid.uuid4())[:12]

    console.print()
    console.print(Panel(
        f"  [{C['cyan']}]Experiment ID:[/] [{C['dim']}]{experiment_id}[/]\n"
        f"  [{C['cyan']}]Data:[/]          [{C['dim']}]{data_dir}[/]\n"
        f"  [{C['cyan']}]Question:[/]      [{C['text']}]{question}[/]",
        border_style=C['border'],
        padding=(0, 1),
    ))
    console.print()

    if not Confirm.ask(
        f"  [bold {C['cyan']}]Launch ARIA?[/]",
        default=True, console=console
    ):
        console.print(f"\n  [{C['muted']}]Cancelled.[/]\n")
        return

    try:
        run_analysis(
            orchestrator=orchestrator,
            experiment_id=experiment_id,
            context={"data_dir": str(data_dir), "user_question": question},
        )
    except KeyboardInterrupt:
        console.print(
            f"\n\n  [{C['amber']}]Interrupted. "
            f"Partial results may be in ~/.aria/reports/[/]\n"
        )
    finally:
        memory.close()


if __name__ == "__main__":
    main()
