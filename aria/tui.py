"""
ARIA Terminal Interface
-----------------------
Professional TUI built with Rich.
Scientific instrument aesthetic: dark, precise, data-forward.
"""

from __future__ import annotations
import os
import sys
import uuid
import time
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.table import Table
from rich.prompt import Prompt, Confirm
from rich.rule import Rule
from rich.align import Align
from rich import box

from aria.memory.memory import ARIAMemory
from aria.agents.orchestrator_agent import OrchestratorAgent
from aria.bus.message_bus import bus, MessageType

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
VERSION = "v0.1.0-alpha"


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
    console.print(Rule(f"[bold {color or C['teal']}]{title}[/]",
                       style=C['border']))
    console.print()


def print_agent_message(agent: str, message: str):
    COLORS = {
        "orchestrator":     C['cyan'],
        "data_audit_agent": C['teal'],
        "rna_agent":        C['green'],
        "chromatin_agent":  C['accent'],
        "genome_arch_agent":"#a78bfa",
        "integration_agent":"#f472b6",
        "narrative_agent":  C['amber'],
    }
    ICONS = {
        "orchestrator":     "[o]",
        "data_audit_agent": "[?]",
        "rna_agent":        "[r]",
        "chromatin_agent":  "[c]",
        "genome_arch_agent":"[g]",
        "integration_agent":"[i]",
        "narrative_agent":  "[n]",
    }
    color  = COLORS.get(agent, C['muted'])
    icon   = ICONS.get(agent, "[-]")
    prefix = Text(f" {icon} [{agent}] ", style=f"bold {color}")
    msg    = Text(message, style=C['text'])
    console.print(prefix + msg)


def print_checkpoint(number: int, title: str,
                     content: str, options: list[str]) -> str:
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

    # Active checkpoint: user must type the number (not just click)
    while True:
        choice = Prompt.ask(
            f"  [bold {C['cyan']}]Enter choice[/]",
            choices=[str(i) for i in range(1, len(options) + 1)],
            console=console,
        )
        return options[int(choice) - 1]


def select_data_directory() -> Path:
    print_section("Data Directory", C['cyan'])
    console.print(
        f"  [{C['muted']}]Enter the path to your raw data directory.[/]\n"
        f"  [{C['dim']}]Example: /data/my_experiment  or  ~/lab/project1[/]\n"
    )
    while True:
        raw_path = Prompt.ask(
            f"  [bold {C['cyan']}]Data path[/]", console=console
        )
        path = Path(raw_path).expanduser().resolve()
        if path.exists() and path.is_dir():
            file_count = sum(1 for f in path.rglob("*") if f.is_file())
            console.print()
            console.print(Panel(
                f"  [{C['green']}]v[/] Directory found\n"
                f"  [{C['muted']}]Files detected: {file_count}[/]\n"
                f"  [{C['dim']}]{path}[/]",
                border_style=C['green'],
                padding=(0, 1),
            ))
            return path
        else:
            console.print(f"\n  [{C['red']}]Directory not found: {path}[/]\n")


def ask_biological_question() -> str:
    print_section("Biological Question", C['cyan'])
    console.print(
        f"  [{C['muted']}]What would you like to know about your data?[/]\n"
        f"  [{C['dim']}]Examples:[/]\n"
        f"  [{C['dim']}]  * What transcription factors are differentially active in condition A vs B?[/]\n"
        f"  [{C['dim']}]  * Which genes show coordinated changes in expression and chromatin accessibility?[/]\n"
        f"  [{C['dim']}]  * What cell types are present and how do they differ between conditions?[/]\n"
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


def _process_pending_checkpoints(orchestrator: OrchestratorAgent,
                                  experiment_id: str):
    TITLES = {
        1: "Data Audit Results",
        2: "Analysis Plan",
        3: "Quality Control / Parameter Review",
        4: "Preliminary Findings",
        5: "Final Report Ready",
    }
    for _ in range(10):
        exp_pending = [
            m for m in bus.get_pending_checkpoints()
            if m.experiment_id == experiment_id
            and not m.payload.get("resolved")
        ]
        if not exp_pending:
            break
        msg     = exp_pending[0]
        cp_num  = msg.payload.get("checkpoint", "?")
        title   = TITLES.get(cp_num, f"Checkpoint {cp_num}")
        choice  = print_checkpoint(
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
            return
        print_agent_message("orchestrator",
                             f"Checkpoint {cp_num} resolved: {choice}")
        time.sleep(0.3)


def run_analysis_with_progress(orchestrator, experiment_id, context):
    result = orchestrator.run(experiment_id, context)
    if result.get("status") != "started":
        console.print(f"[{C['red']}]Failed to start analysis.[/]")
        return
    print_agent_message("orchestrator",
        f"Biological intent: {result['intent'].get('summary', '')}")
    with console.status(
        f"[{C['teal']}]DataAuditAgent scanning directory...[/]",
        spinner="dots"
    ):
        orchestrator.run_audit(experiment_id)
    _process_pending_checkpoints(orchestrator, experiment_id)


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

    print_section("Analysis Running", C['green'])
    run_analysis_with_progress(
        orchestrator=orchestrator,
        experiment_id=experiment_id,
        context={"data_dir": str(data_dir), "user_question": question},
    )
    memory.close()


if __name__ == "__main__":
    main()
