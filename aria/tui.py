"""
ARIA TUI — Mock Warning Additions
-----------------------------------
Merge these two functions into aria/tui.py before select_data_directory().

Fixes DeepSeek code review finding:
  "Users can receive mock results without knowing. The TUI must warn
   explicitly when a fallback mock is active."

Usage in agents:
    from aria.tui import print_mock_warning
    if result.get("note", "").startswith("Mock"):
        print_mock_warning("genome_arch_agent", "HiC", "aria-hic-env")

Usage at startup (add to main() after print_banner()):
    check_environment_status()
    console.print()
"""


def print_mock_warning(agent: str, modality: str, missing_env: str):
    """
    Unmissable warning when an agent returns simulated (mock) results.
    Called when a script result contains 'note: Mock ...' in its payload.
    """
    console.print()
    console.print(Panel(
        f"  [bold yellow]MOCK RESULTS -- NOT REAL ANALYSIS[/]\n\n"
        f"  {agent} returned SIMULATED data for: {modality}\n\n"
        f"  The '{missing_env}' environment is not installed.\n\n"
        f"  To enable real analysis:\n"
        f"  [bold cyan]conda env create -f envs/{missing_env}.yml[/]",
        border_style="yellow",
        title="[bold yellow]  !! SIMULATED RESULTS -- DO NOT PUBLISH !!  [/]",
        title_align="left",
        padding=(1, 2),
    ))
    console.print()


def check_environment_status() -> dict:
    """
    Check which analytical Conda environments are installed.
    Called at ARIA startup. Shows ready vs missing stacks clearly.
    Missing stacks will produce mock results -- user must know this.
    """
    try:
        from aria.utils.environment_manager import env_manager
        report  = env_manager.get_status_report()
        ready   = report.get("ready_stacks", [])
        missing = report.get("missing_stacks", [])

        if not report.get("conda_available", False):
            console.print(Panel(
                "Conda not found in PATH.\n"
                "Isolated environment execution disabled.\n"
                "Install Miniforge: https://github.com/conda-forge/miniforge",
                border_style="red",
                title="[bold red]Conda not available[/]",
                padding=(0, 1),
            ))
            return report

        if missing:
            console.print(Panel(
                f"Ready:   {', '.join(ready) if ready else 'none'}\n"
                f"Missing: {', '.join(missing)}\n\n"
                f"[yellow]Missing environments will produce MOCK results.[/]\n"
                f"[dim]Install: conda env create -f envs/<n>.yml[/]",
                border_style="yellow",
                title="[dim]Environment Status[/]",
                title_align="left",
                padding=(0, 1),
            ))
        else:
            console.print(Panel(
                "[green]All analytical environments installed and ready.[/]",
                border_style="green",
                title="[dim]Environment Status[/]",
                title_align="left",
                padding=(0, 1),
            ))
        return report
    except Exception:
        return {}
