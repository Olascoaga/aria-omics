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
  4. DesignAgent → Checkpoints 2.1‑2.6 (experimental design)  [v4.0]
  5. Checkpoint 2 (confirm analysis plan)
  6. Dispatcher thread launched in background
  7. [POLLING LOOP] — stays alive, showing agent progress +
     surfacing checkpoints 3-5 as agents publish them
  8. Thread completes → final summary shown
"""

from __future__ import annotations
import os
import re
import sys
import uuid
import time
import threading
import select
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
from aria.runtime.experiment_view import status_text, build_snapshot
from aria.ui.brand import ARIA_BANNER, TAGLINE
from aria.version import __version__

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
    "raw_ingestion_agent": "#22d3ee",
    "setup_agent":       "#fbbf24",
    "scrna_agent":       C['green'],
    "bulk_rna_agent":    "#34d399",
    "rna_agent":         C['green'],
    "chromatin_agent":   C['accent'],
    "genome_arch_agent": "#a78bfa",
    "integration_agent": "#f472b6",
    "narrative_agent":   C['amber'],
    "debate_council":    C['red'],
    "design_agent":      C['teal'],   # v4.0
}

AGENT_ICONS = {
    "orchestrator":      "(orch )",
    "data_audit_agent":  "(audit)",
    "raw_ingestion_agent": "(ingest)",
    "setup_agent":       "(setup)",
    "scrna_agent":       "(scRNA)",
    "bulk_rna_agent":    "(bulk )",
    "rna_agent":         "(rna  )",
    "chromatin_agent":   "(chrom)",
    "genome_arch_agent": "( 3D  )",
    "integration_agent": "(integ)",
    "narrative_agent":   "(narr )",
    "debate_council":    "(deb  )",
    "design_agent":      "(design)",  # v4.0
}

# Human-readable names for display (paired with icons)
AGENT_NAMES = {
    "orchestrator":      "Orchestrator",
    "data_audit_agent":  "Data Audit",
    "raw_ingestion_agent": "Raw Ingest",
    "setup_agent":       "Setup",
    "scrna_agent":       "scRNA",
    "bulk_rna_agent":    "Bulk RNA",
    "rna_agent":         "RNA",
    "chromatin_agent":   "Chromatin",
    "genome_arch_agent": "3D Genome",
    "integration_agent": "Integration",
    "narrative_agent":   "Report",
    "debate_council":    "Debate",
    "design_agent":      "Design",   # v4.0
}

CHECKPOINT_TITLES = {
    1:   "Data Audit Results",
    2:   "Analysis Plan",
    2.1: "Experimental Groups",        # v4.0
    2.2: "Organism",                   # v4.0
    2.3: "Experimental Factor",        # v4.0
    2.4: "Batch Effects",              # v4.0
    2.5: "Pseudoreplication Check",    # v4.0
    2.6: "Design Confirmation",        # v4.0
    3:   "Quality Control / Parameter Decision",
    4:   "Preliminary Findings",
    5:   "Final Report Ready",
}

console = Console(highlight=False)

VERSION = f"v{__version__}"


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


def _prompt_action() -> str:
    """Read the top-level action, exiting cleanly when stdin is closed.

    `conda run aria` can be launched from non-interactive runners with stdin at
    EOF. Rich's Prompt.ask raises EOFError in that case; treat it as "exit" so
    the CLI does not print a traceback before any analysis has started.
    """
    try:
        return Prompt.ask(
            f"  [{C['cyan']}]Action[/]",
            choices=["new", "exit"],
            default="new",
            console=console,
        )
    except EOFError:
        console.print(
            f"\n  [{C['muted']}]No interactive terminal detected; exiting.[/]\n"
            f"  [{C['muted']}]Launch ARIA from a real shell with:[/]\n"
            f"  [{C['cyan']}]conda activate aria-env && aria[/]\n"
        )
        return "exit"


def print_agent_message(agent: str, message: str):
    color  = AGENT_COLORS.get(agent, C['muted'])
    icon   = AGENT_ICONS.get(agent, "[-]")
    prefix = Text(f" {icon} [{agent}] ", style=f"bold {color}")
    msg    = Text(message, style=C['text'])
    console.print(prefix + msg)


def _collect_metadata_corrections(exp_context: dict) -> dict:
    """CHECKPOINT-1 correction sub-flow: pick the real modality and species.

    Returns a corrections dict ``{modality, organism, genome}`` (only the keys the
    user actually changed). Empty dict if the user keeps everything.
    """
    from aria.agents.data_audit_agent import (
        SUPPORTED_MODALITIES, default_genome_for_organism,
    )

    corrections: dict = {}
    current_mods = list((exp_context.get("modalities") or {}).keys())
    console.print(
        f"\n  [{C['muted']}]Currently detected modality: "
        f"{', '.join(current_mods) or 'unknown'}[/]"
    )
    console.print(f"  [{C['muted']}]Select the correct modality:[/]")
    for i, m in enumerate(SUPPORTED_MODALITIES, 1):
        console.print(f"    [{C['cyan']}][{i}][/] {m}")
    keep_idx = len(SUPPORTED_MODALITIES) + 1
    console.print(f"    [{C['muted']}][{keep_idx}][/] keep current")
    m_choice = Prompt.ask(
        f"  [bold {C['cyan']}]Modality[/]",
        choices=[str(i) for i in range(1, keep_idx + 1)],
        default=str(keep_idx), console=console,
    )
    if int(m_choice) != keep_idx:
        corrections["modality"] = SUPPORTED_MODALITIES[int(m_choice) - 1]

    organisms = [
        "Homo sapiens", "Mus musculus", "Drosophila melanogaster",
        "C. elegans", "Danio rerio", "S. cerevisiae",
    ]
    console.print(
        f"\n  [{C['muted']}]Current organism: "
        f"{exp_context.get('organism', 'unknown')}[/]"
    )
    console.print(f"  [{C['muted']}]Select the correct species:[/]")
    for i, o in enumerate(organisms, 1):
        console.print(f"    [{C['cyan']}][{i}][/] {o}")
    other_idx = len(organisms) + 1
    keep_org_idx = other_idx + 1
    console.print(f"    [{C['cyan']}][{other_idx}][/] Other (type it)")
    console.print(f"    [{C['muted']}][{keep_org_idx}][/] keep current")
    o_choice = Prompt.ask(
        f"  [bold {C['cyan']}]Species[/]",
        choices=[str(i) for i in range(1, keep_org_idx + 1)],
        default=str(keep_org_idx), console=console,
    )
    organism = None
    if int(o_choice) == other_idx:
        organism = Prompt.ask(
            f"  [bold {C['cyan']}]Type the species "
            f"(e.g. \"Rattus norvegicus\")[/]", console=console,
        ).strip()
    elif int(o_choice) != keep_org_idx:
        organism = organisms[int(o_choice) - 1]
    if organism:
        corrections["organism"] = organism
        genome = default_genome_for_organism(organism)
        if genome:
            corrections["genome"] = genome
            console.print(
                f"  [{C['muted']}]Reference genome set to {genome} "
                f"for {organism}.[/]"
            )

    if corrections:
        console.print(
            f"\n  [{C['green']}]Corrections applied:[/] "
            f"{corrections}\n"
        )
    else:
        console.print(f"\n  [{C['muted']}]No changes made.[/]\n")
    return corrections


def print_checkpoint(number, title: str,
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
        _discard_queued_stdin_lines()
        choice = Prompt.ask(
            f"  [bold {C['cyan']}]Enter choice[/]",
            choices=[str(i) for i in range(1, len(options) + 1)],
            console=console,
        )
        selected = options[int(choice) - 1]
        # Free-text follow-up for "Other" options (e.g. custom organism)
        if selected.lower().startswith("other"):
            console.print(
                f"\n  [{C['muted']}]Type the value "
                f"(free text, e.g. species \"Genus species (assembly)\" or "
                f"a custom group mapping):[/]"
            )
            custom = Prompt.ask(
                f"  [bold {C['cyan']}]Value[/]",
                console=console,
            ).strip()
            return custom if custom else selected
        return selected


def print_agent_progress(msg: Message):
    """Display a STATUS message from an agent with clear identification."""
    agent    = msg.sender
    progress = msg.payload.get("progress", 0)
    # Agents publish the status text under "message"; read either key (U0).
    status   = status_text(msg.payload)
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

_GEO_ACCESSION_RE = re.compile(
    r"^(GSE\d+|SRP\d+|PRJNA\d+|ERP\d+|DRP\d+)$", re.IGNORECASE
)


def select_data_directory() -> tuple:
    """
    Returns (data_dir: Path, geo_meta: dict | None).
    Accepts local paths OR GEO/SRA accession numbers (e.g. GSE183948).
    """
    print_section("Data / Accession", C['cyan'])
    console.print(
        f"  [{C['muted']}]Enter a local data directory or a GEO/SRA accession.[/]\n"
        f"  [{C['dim']}]Local:   /data/my_experiment  or  ~/lab/project1[/]\n"
        f"  [{C['dim']}]Public:  GSE183948   SRP123456   PRJNA987654[/]\n"
    )
    while True:
        raw = Prompt.ask(f"  [bold {C['cyan']}]Data path or accession[/]",
                         console=console)
        raw = raw.strip()

        # ── GEO / SRA accession ──────────────────────────────────────────
        if _GEO_ACCESSION_RE.match(raw):
            geo_result = _resolve_geo_accession(raw.upper())
            if geo_result:
                return geo_result["local_dir"], geo_result
            continue

        # ── Local path ───────────────────────────────────────────────────
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
            return path, None
        console.print(f"\n  [{C['red']}]Directory not found: {path}[/]\n")


def _resolve_geo_accession(accession: str) -> dict | None:
    """Download a GEO/SRA accession and return GEOConnector result dict."""
    from aria.connectors.geo_connector import GEOConnector

    console.print()
    status_lines: list[str] = []

    def _status(msg: str):
        status_lines.append(msg)

    with console.status(
        f"  [{C['teal']}]Connecting to NCBI GEO...[/]", spinner="dots"
    ):
        try:
            gc     = GEOConnector()
            result = gc.fetch(accession, status_callback=_status)
        except Exception as e:
            console.print(f"\n  [{C['red']}]GEO fetch failed: {e}[/]\n")
            return None

    # Show what was found
    design = result.get("inferred_design", {})
    groups = design.get("groups", {})
    files  = result.get("files", {})
    n_files = sum(len(v) for v in files.values() if isinstance(v, list))

    group_str = "  ".join(
        f"[{C['green']}]{g}[/] ({len(s)})" for g, s in list(groups.items())[:5]
    )
    has_raw = files.get("fastq_pending")

    file_summary = []
    if files.get("counts"):
        file_summary.append(f"{len(files['counts'])} count matrix file(s)")
    if files.get("h5ad"):
        file_summary.append(f"{len(files['h5ad'])} h5ad file(s)")
    if files.get("h5"):
        file_summary.append(f"{len(files['h5'])} 10x h5 file(s)")
    if files.get("mtx"):
        file_summary.append(f"{len(files['mtx'])} MEX matrix file(s)")
    if has_raw:
        file_summary.append("FASTQs available via SRA (not auto-downloaded)")
    if not file_summary:
        file_summary.append("no processed files found")

    console.print()
    console.print(Panel(
        f"  [{C['cyan']}]Accession:[/]  {accession}\n"
        f"  [{C['cyan']}]Title:[/]      {result.get('title','')[:80]}\n"
        f"  [{C['cyan']}]Organism:[/]   {result.get('organism','')} "
        f"({result.get('genome','')})\n"
        f"  [{C['cyan']}]Platform:[/]   {result.get('data_type','unknown')}\n"
        f"  [{C['cyan']}]Samples:[/]    {result.get('n_samples', '?')}\n"
        f"  [{C['cyan']}]Groups:[/]     {group_str or 'not inferred'}\n"
        f"  [{C['cyan']}]Files:[/]      {', '.join(file_summary)}\n"
        f"  [{C['dim']}]Local cache: {result['local_dir']}[/]",
        title=f"[bold {C['teal']}]GEO Dataset Found[/]",
        border_style=C['teal'],
        padding=(0, 1),
    ))
    console.print()

    if n_files == 0 and not has_raw:
        console.print(
            f"  [{C['amber']}]Warning: no downloadable processed files found "
            f"for {accession}. You may need to provide local files.[/]\n"
        )

    return result


def _validate_textual_intake_data(raw: str) -> str | None:
    """Validate an intake data value without resolving GEO over the network.

    Returns an error message (kept visible in the intake) or ``None`` when the
    value is a GEO/SRA accession or an existing local directory. This mirrors the
    classic re-prompt loop so a bad path keeps the user in the front door instead
    of silently dropping to the terminal after the TUI closes.
    """
    raw = (raw or "").strip()
    if not raw:
        return "Enter a data directory or accession."
    if _GEO_ACCESSION_RE.match(raw):
        return None
    path = Path(raw).expanduser()
    if not path.exists():
        return f"Path not found: {path}"
    if not path.is_dir():
        return f"That is a file, not a folder. Enter its directory: {path.parent}"
    return None


def _resolve_textual_intake_data(raw: str) -> tuple[Path, dict | None] | None:
    """Resolve a Textual intake data value into the classic run inputs."""
    raw = raw.strip()
    if _GEO_ACCESSION_RE.match(raw):
        geo_result = _resolve_geo_accession(raw.upper())
        if geo_result:
            return Path(geo_result["local_dir"]), geo_result
        return None

    path = Path(raw).expanduser().resolve()
    if path.exists() and path.is_dir():
        return path, None

    console.print(f"\n  [{C['red']}]Directory not found: {path}[/]\n")
    return None


def ask_biological_question() -> str:
    print_section("Biological Question", C['cyan'])
    console.print(
        f"  [{C['muted']}]What would you like to know about your data?[/]\n"
        f"  [{C['dim']}]Examples:[/]\n"
        f"  [{C['dim']}]  * What TFs are differentially active in condition A vs B?[/]\n"
        f"  [{C['dim']}]  * Which genes show coordinated RNA and chromatin changes?[/]\n"
        f"  [{C['dim']}]  * What cell types are present and how do they differ?[/]\n"
        f"  [{C['dim']}]Paste long prompts freely. Finish with a line containing only END.[/]\n"
    )
    return _read_multiline_question()


def _read_multiline_question() -> str:
    lines: list[str] = []
    while True:
        prompt = (
            f"\n  [bold {C['cyan']}]Your question[/] "
            f"[{C['dim']}](END to finish)[/]: "
            if not lines else
            f"  [{C['dim']}]...[/] "
        )
        line = console.input(prompt)
        if line.strip() == "END":
            break
        lines.append(line.rstrip("\n"))
    return "\n".join(lines).strip()


def _read_pasted_stdin_lines(timeout: float = 0.08) -> list[str]:
    """
    Capture lines already queued by a terminal paste.

    Rich Prompt.ask reads a single line, which lets the rest of a pasted
    multi-line question leak into the following yes/no and checkpoint prompts.
    This helper drains only immediately available stdin lines; normal typed
    follow-up answers are left alone because they are not queued yet.
    """
    if not sys.stdin.isatty():
        return []

    lines: list[str] = []
    wait = timeout
    while True:
        try:
            ready, _, _ = select.select([sys.stdin], [], [], wait)
        except (OSError, ValueError):
            break
        if not ready:
            break
        line = sys.stdin.readline()
        if line == "":
            break
        lines.append(line.rstrip("\n"))
        wait = 0.01
    return lines


def _discard_queued_stdin_lines(timeout: float = 0.02) -> None:
    """Drop stale pasted input before yes/no or checkpoint prompts."""
    _read_pasted_stdin_lines(timeout=timeout)


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
    status = _drain_checkpoints(orchestrator, experiment_id)
    if status == "cancelled":
        return

    # After CP2, dispatcher thread is running in background.
    # Phase 3: Live polling loop — stays alive until thread done.
    _live_analysis_loop(orchestrator, experiment_id)


def _drain_checkpoints(orchestrator: OrchestratorAgent,
                        experiment_id: str,
                        max_rounds: int = 10):
    """
    Process all currently pending checkpoints synchronously.
    Used for pre-dispatch checkpoints (1 and 2, and now the design checkpoints).
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

        # CHECKPOINT-1 "Correct metadata": let the user pick the real modality and
        # species, and pass the corrections to the orchestrator.
        corrections = None
        if cp_num == 1 and "correct" in choice.lower():
            corrections = _collect_metadata_corrections(
                msg.payload.get("context", {}).get("exp_context", {})
            )

        result = orchestrator.on_checkpoint_resolved(
            message_id=msg.id,
            user_decision=choice,
            experiment_id=experiment_id,
            corrections=corrections,
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
            last_message_time = time.time()     # reset idle counter
            idle_warning_shown = False

            if msg.type == MessageType.STATUS:
                seen_message_ids.add(msg.id)
                print_agent_progress(msg)
                last_status_text = status_text(msg.payload)
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
                seen_message_ids.add(msg.id)
                print_finding(msg)
            elif msg.type != MessageType.ESCALATION:
                seen_message_ids.add(msg.id)

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

            # Normal checkpoints (1, 2, 3, 4, 5, and now the design sub-checkpoints)
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

    # Report path (resolved by the shared U0 read-model: payload key OR the
    # "Report saved: <path>" status line, which the old payload-only scan missed).
    report_path = build_snapshot(experiment_id).report_path
    if report_path:
        console.print()
        console.print(Panel(
            f"  [{C['cyan']}]Report:[/] [{C['dim']}]{report_path}[/]\n"
            f"  [{C['muted']}]Open in browser to view the full analysis.[/]",
            border_style=C['green'],
            padding=(0, 1),
        ))


# ── Main ──────────────────────────────────────────────────────────────────────

def _use_cockpit(reproducible_mode: bool) -> bool:
    """Use the Textual cockpit when it is available, on a TTY, and not opted out.

    Falls back to the classic Rich TUI otherwise. The headless runner remains the
    canonical reproducible path and is never affected by this choice.
    """
    if reproducible_mode:
        return False
    if os.environ.get("ARIA_NO_TUI"):
        return False
    if "--classic-tui" in sys.argv:
        return False
    if not (sys.stdin.isatty() and sys.stdout.isatty()):
        return False
    try:
        from aria.ui.cockpit import cockpit_available
        return cockpit_available()
    except Exception:
        return False


def _launch_context(data_dir: Path, question: str,
                    geo_meta: dict | None,
                    reproducible_mode: bool) -> dict:
    ctx: dict = {
        "data_dir": str(data_dir),
        "user_question": question,
        "reproducible_mode": reproducible_mode,
    }
    if geo_meta:
        ctx["geo_metadata"] = geo_meta
    return ctx


def _data_summary_lines(data_dir: Path, geo_meta: dict | None) -> str:
    if not geo_meta:
        return f"  [{C['cyan']}]Data:[/]          [{C['dim']}]{data_dir}[/]"

    acc = geo_meta.get("accession", "")
    title = geo_meta.get("title", "")[:70]
    org = geo_meta.get("organism", "")
    dtype = geo_meta.get("data_type", "")
    n_sam = geo_meta.get("n_samples", "")
    return (
        f"  [{C['cyan']}]Accession:[/]     [{C['dim']}]{acc}[/]\n"
        f"  [{C['cyan']}]Title:[/]         [{C['text']}]{title}[/]\n"
        + (f"  [{C['cyan']}]Organism:[/]     [{C['dim']}]{org}[/]\n" if org else "")
        + (f"  [{C['cyan']}]Type / Samples:[/] [{C['dim']}]{dtype}  ·  {n_sam} samples[/]\n" if dtype or n_sam else "")
        + f"  [{C['cyan']}]Local cache:[/]  [{C['dim']}]{data_dir}[/]"
    )


def _print_launch_summary(experiment_id: str, data_dir: Path,
                          geo_meta: dict | None, question: str) -> None:
    console.print(Panel(
        f"  [{C['cyan']}]Experiment ID:[/] [{C['dim']}]{experiment_id}[/]\n"
        + _data_summary_lines(data_dir, geo_meta) + "\n"
        + f"  [{C['cyan']}]Question:[/]      [{C['text']}]{question}[/]",
        border_style=C['border'],
        padding=(0, 1),
    ))


def _run_cockpit_front_door(memory: ARIAMemory,
                            orchestrator: OrchestratorAgent,
                            reproducible_mode: bool) -> bool:
    """Run the Textual intake first, then transition to the cockpit run view."""
    from aria.ui.cockpit import launch_cockpit
    from aria.ui.intake import launch_intake
    from aria.runtime.experiment_view import build_history

    intake = launch_intake(
        startup_context=memory.startup_context(),
        experiments=memory.list_wings(),
        history=build_history(memory),
        version=VERSION,
        data_validator=_validate_textual_intake_data,
    )
    if intake is None:
        console.print(f"\n  [{C['muted']}]Goodbye.[/]\n")
        return True

    resolved = _resolve_textual_intake_data(intake.data_input)
    if resolved is None:
        # Pre-validated for local paths, so this is a GEO/SRA accession that did
        # not resolve. Be explicit instead of dropping silently to the terminal.
        console.print(
            f"\n  [{C['red']}]Could not resolve "
            f"'{intake.data_input}'.[/] "
            f"[{C['muted']}]Check the accession or your network/air-gap "
            f"settings and run ARIA again.[/]\n"
        )
        return True
    data_dir, geo_meta = resolved
    experiment_id = str(uuid.uuid4())[:12]
    ctx = _launch_context(
        data_dir, intake.question, geo_meta, reproducible_mode
    )
    result = launch_cockpit(orchestrator, experiment_id, ctx)
    if result is None:
        console.print(
            f"\n  [{C['red']}]The analysis could not start.[/] "
            f"[{C['muted']}]See the log at "
            f"~/.aria/logs/exp_{experiment_id}.log[/]\n"
        )
        return True
    _print_final_summary(experiment_id)
    return True


def _print_llm_runtime_error(exc: RuntimeError) -> bool:
    """Return True if a RuntimeError was handled as an LLM config failure."""
    msg = str(exc)
    if "models failed for tier" not in msg and "No model is configured" not in msg:
        return False
    from aria.llm.provider import diagnose_llm_failure
    console.print(Panel(
        diagnose_llm_failure(exc),
        border_style=C['red'],
        title=f"[{C['red']}]LLM provider unavailable[/]",
        title_align="left",
        padding=(0, 1),
    ))
    return True


def main():
    args = [arg for arg in sys.argv[1:] if arg != "--reproducible"]
    if args and args[0] == "doctor":
        from aria.doctor import main as doctor_main

        code = doctor_main(args[1:])
        if code:
            raise SystemExit(code)
        return

    # W-LEDGER: `aria diff A B` and `aria export <reportDir>` over the run ledger.
    if args and args[0] in ("diff", "export"):
        from aria.agents.narrative.ledger_export import cli_main

        raise SystemExit(cli_main(args))

    reproducible_mode = "--reproducible" in sys.argv[1:]
    memory       = ARIAMemory()
    orchestrator = OrchestratorAgent(memory)

    try:
        if _use_cockpit(reproducible_mode):
            if _run_cockpit_front_door(
                memory, orchestrator, reproducible_mode
            ):
                return

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

        action = _prompt_action()
        if action == "exit":
            console.print(f"\n  [{C['muted']}]Goodbye.[/]\n")
            return

        data_dir, geo_meta = select_data_directory()
        question           = ask_biological_question()
        experiment_id      = str(uuid.uuid4())[:12]

        console.print()
        _print_launch_summary(experiment_id, data_dir, geo_meta, question)
        console.print()

        _discard_queued_stdin_lines()
        if not Confirm.ask(
            f"  [bold {C['cyan']}]Launch ARIA?[/]",
            default=True, console=console
        ):
            console.print(f"\n  [{C['muted']}]Cancelled.[/]\n")
            return

        ctx = _launch_context(
            data_dir, question, geo_meta, reproducible_mode
        )
        run_analysis(
            orchestrator=orchestrator,
            experiment_id=experiment_id,
            context=ctx,
        )
    except KeyboardInterrupt:
        console.print(
            f"\n\n  [{C['amber']}]Interrupted. "
            f"Partial results may be in ~/.aria/reports/[/]\n"
        )
    except RuntimeError as exc:
        # An unreachable LLM provider should give an actionable hint, not a raw
        # traceback (real-run bug 2026-06-04).
        if not _print_llm_runtime_error(exc):
            raise
    finally:
        memory.close()


if __name__ == "__main__":
    main()
