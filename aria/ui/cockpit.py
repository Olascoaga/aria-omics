"""ARIA cockpit — Textual presentation shell.

The cockpit renders the U0 read-model on a timer and routes checkpoint decisions
back through the orchestrator's existing ``on_checkpoint_resolved`` — it is a skin
over the same governed decision path, never a new one. It owns no scientific
state. Requires the ``tui`` extra (``pip install aria-omics[tui]``); when Textual
is absent, ``aria.tui.main`` falls back to the classic Rich TUI and the headless
runner stays the canonical reproducible path.

Design for testability: :class:`AriaCockpit` takes a ``snapshot_provider`` and a
``checkpoint_resolver`` callable, so it can be pilot-driven with fakes. The thin
:func:`launch_cockpit` wires those to ``build_snapshot`` and the real orchestrator.
"""

from __future__ import annotations

from typing import Any, Callable, Optional

from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.widgets import Footer, Static

from aria.runtime.experiment_view import ExperimentSnapshot
from aria.ui import render

SnapshotProvider = Callable[[], ExperimentSnapshot]
# (message_id, user_decision, experiment_id) -> Any
CheckpointResolver = Callable[..., Any]
MetaProvider = Callable[[], dict]
# Front-door wiring (single-app intake -> run transition):
IntakeScreenFactory = Callable[[], Any]
# IntakeResult -> (experiment_id, provider, resolver, meta_provider) | None
RunStarter = Callable[[Any], Optional[tuple]]


class AriaCockpit(App):
    """Cockpit app. Pure presentation over the read-model."""

    CSS = """
    #body { height: 1fr; }
    #left { width: 34; }
    #right { width: 46; }
    #center { width: 1fr; }
    Static { height: auto; }
    #agents { height: 1fr; }
    #mode-bar { height: 1; margin-bottom: 1; }
    #findings { height: 1fr; }
    #ledger { height: 1fr; }
    #readiness { height: 1fr; }
    #resources { height: 1fr; }
    #artifacts { height: 1fr; }
    """

    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("l", "toggle_ledger", "Ledger"),
        Binding("r", "toggle_readiness", "Readiness"),
        Binding("u", "toggle_resources", "Resources"),
        Binding("a", "toggle_artifacts", "Artifacts"),
        Binding("e", "edit_design", "Edit groups"),
        Binding("1", "choose(1)", "Opt 1", show=False),
        Binding("2", "choose(2)", "Opt 2", show=False),
        Binding("3", "choose(3)", "Opt 3", show=False),
        Binding("4", "choose(4)", "Opt 4", show=False),
        Binding("5", "choose(5)", "Opt 5", show=False),
        Binding("6", "choose(6)", "Opt 6", show=False),
    ]

    def __init__(self, snapshot_provider: Optional[SnapshotProvider] = None,
                 checkpoint_resolver: Optional[CheckpointResolver] = None, *,
                 experiment_id: str = "",
                 version: str = "",
                 meta_provider: Optional[MetaProvider] = None,
                 poll_interval: float = 0.5,
                 exit_on_done: bool = False,
                 intake_screen_factory: Optional[IntakeScreenFactory] = None,
                 run_starter: Optional[RunStarter] = None):
        super().__init__()
        self._provider = snapshot_provider
        self._resolver = checkpoint_resolver
        self.experiment_id = experiment_id
        self._version = version
        self._meta_provider = meta_provider or (lambda: {})
        self._poll_interval = poll_interval
        self._exit_on_done = exit_on_done
        # Front-door mode: collect intake in this SAME app, then start the run,
        # so there is no flicker back to the console between two App.run() calls.
        self._intake_screen_factory = intake_screen_factory
        self._run_starter = run_starter
        self._snap: Optional[ExperimentSnapshot] = None
        self._done_seen = False
        # findings | ledger | readiness | resources | artifacts
        self._center_mode = "findings"

    def compose(self) -> ComposeResult:
        yield Horizontal(
            Vertical(Static(id="run-header"), Static(id="agents"), id="left"),
            Vertical(
                Static(id="timeline"),
                Static(id="mode-bar"),
                Static(id="findings"),
                Static(id="ledger"),
                Static(id="readiness"),
                Static(id="resources"),
                Static(id="artifacts"),
                id="center",
            ),
            Vertical(Static(id="checkpoint"), id="right"),
            id="body",
        )
        yield Footer()

    def on_mount(self) -> None:
        self._apply_center_mode()
        if self._intake_screen_factory is not None:
            # Front door: collect intake first (same app), then begin the run.
            self.push_screen(self._intake_screen_factory(), self._on_intake_result)
            return
        self._begin_polling()

    def _begin_polling(self) -> None:
        self.refresh_snapshot()
        self.set_interval(self._poll_interval, self.refresh_snapshot)

    # ── Front-door transition (intake screen -> run view, one app) ────────────
    def _on_intake_result(self, result: Any) -> None:
        if result is None:
            self.exit(None)
            return
        # The intake screen has popped; reveal the run view with a banner while
        # the run starts off the UI thread (data scan can take a moment).
        try:
            self.query_one("#run-header", Static).update(
                render.render_status_banner("Starting analysis…", self._version))
        except Exception:
            pass
        self.run_worker(lambda: self._start_run(result), thread=True,
                        exclusive=True, name="start_run")

    def _start_run(self, result: Any) -> None:
        handles: Optional[tuple] = None
        try:
            if self._run_starter is not None:
                handles = self._run_starter(result)
        except Exception:
            handles = None
        self.call_from_thread(self._begin_run, handles)

    def _begin_run(self, handles: Optional[tuple]) -> None:
        if not handles:
            # Honest, in-TUI failure instead of a silent drop to the console.
            try:
                self.query_one("#run-header", Static).update(
                    render.render_status_banner(
                        "Could not start the analysis. Press q to exit.",
                        self._version, error=True))
            except Exception:
                pass
            return
        experiment_id, provider, resolver, meta_provider = handles
        self.experiment_id = experiment_id
        self._provider = provider
        self._resolver = resolver
        self._meta_provider = meta_provider or (lambda: {})
        self._begin_polling()

    def _apply_center_mode(self) -> None:
        self.query_one("#findings", Static).display = \
            self._center_mode == "findings"
        self.query_one("#ledger", Static).display = \
            self._center_mode == "ledger"
        self.query_one("#readiness", Static).display = \
            self._center_mode == "readiness"
        self.query_one("#resources", Static).display = \
            self._center_mode == "resources"
        self.query_one("#artifacts", Static).display = \
            self._center_mode == "artifacts"
        self.query_one("#mode-bar", Static).update(
            render.render_mode_bar(self._center_mode))

    # ── Rendering ────────────────────────────────────────────────────────────
    def refresh_snapshot(self) -> None:
        if self._provider is None:
            return
        try:
            snap = self._provider()
        except Exception:
            return
        self._snap = snap
        meta = {}
        try:
            meta = self._meta_provider() or {}
        except Exception:
            meta = {}

        self.query_one("#run-header", Static).update(render.render_run_header(
            snap, version=meta.get("version", self._version),
            data_dir=meta.get("data_dir"),
            modalities=meta.get("modalities"),
            organism=meta.get("organism"),
            air_gapped=bool(meta.get("air_gapped")),
        ))
        self.query_one("#timeline", Static).update(render.render_timeline(snap))
        self.query_one("#agents", Static).update(
            render.render_agent_progress(snap))
        self.query_one("#findings", Static).update(render.render_findings(snap))
        self.query_one("#ledger", Static).update(render.render_ledger(snap))
        self.query_one("#readiness", Static).update(
            render.render_readiness(snap))
        self.query_one("#resources", Static).update(
            render.render_resources(snap))
        self.query_one("#artifacts", Static).update(
            render.render_artifacts(snap))
        self.query_one("#checkpoint", Static).update(
            render.render_checkpoint(snap))

        if snap.done and not self._done_seen:
            self._done_seen = True
            if snap.artifacts:
                self._center_mode = "artifacts"
                self._apply_center_mode()
        if snap.done and self._exit_on_done:
            self.exit(snap)

    # ── Actions ──────────────────────────────────────────────────────────────
    def action_toggle_ledger(self) -> None:
        self._center_mode = "findings" if self._center_mode == "ledger" \
            else "ledger"
        self._apply_center_mode()

    def action_toggle_readiness(self) -> None:
        self._center_mode = "findings" if self._center_mode == "readiness" \
            else "readiness"
        self._apply_center_mode()

    def action_toggle_resources(self) -> None:
        self._center_mode = "findings" if self._center_mode == "resources" \
            else "resources"
        self._apply_center_mode()

    def action_toggle_artifacts(self) -> None:
        self._center_mode = "findings" if self._center_mode == "artifacts" \
            else "artifacts"
        self._apply_center_mode()

    def action_edit_design(self) -> None:
        """Open the U2 tabular design editor for the groups checkpoint (CP2.1)."""
        snap = self._snap
        if snap is None or snap.pending_checkpoint is None:
            return
        cp = snap.pending_checkpoint
        proposed = (cp.context or {}).get("proposed_groups")
        if cp.number != 2.1 or not proposed:
            return

        def _on_submit(groups_json: str) -> None:
            try:
                self._resolver(
                    message_id=cp.message_id,
                    user_decision=groups_json,
                    experiment_id=self.experiment_id,
                )
            except Exception:
                return
            self.refresh_snapshot()

        from aria.ui.design_editor import DesignEditorScreen
        self.push_screen(DesignEditorScreen(proposed, _on_submit))

    def action_choose(self, idx: int) -> None:
        snap = self._snap
        if snap is None or snap.pending_checkpoint is None:
            return
        cp = snap.pending_checkpoint
        if not (1 <= idx <= len(cp.options)):
            return
        choice = cp.options[idx - 1]
        try:
            self._resolver(
                message_id=cp.message_id,
                user_decision=choice,
                experiment_id=self.experiment_id,
            )
        except Exception:
            # A failed resolution must not crash the cockpit; the next poll
            # re-surfaces the still-pending checkpoint.
            return
        self.refresh_snapshot()

    def action_quit(self) -> None:
        self.exit(self._snap)


# ── Launcher / availability ──────────────────────────────────────────────────

def cockpit_available() -> bool:
    """True if the Textual cockpit can be launched (the ``tui`` extra is present)."""
    try:
        import textual  # noqa: F401
        return True
    except Exception:
        return False


def _session_for(orchestrator: Any, experiment_id: str) -> Any:
    sessions = getattr(orchestrator, "_sessions", None)
    if isinstance(sessions, dict):
        return sessions.get(experiment_id)
    return None


def _meta_for(orchestrator: Any, experiment_id: str, *,
              version: str, data_dir: Optional[str]) -> dict:
    session = _session_for(orchestrator, experiment_id)
    exp_ctx = getattr(session, "exp_context", None) or {}
    modalities = list((exp_ctx.get("modalities") or {}).keys()) or None
    return {
        "version": version,
        "data_dir": data_dir,
        "modalities": modalities,
        "organism": exp_ctx.get("organism"),
        "air_gapped": bool(exp_ctx.get("air_gapped")),
    }


def _start_run_handles(orchestrator: Any, experiment_id: str, context: dict,
                       version: str) -> Optional[tuple]:
    """Start the orchestrator + audit and return the cockpit run handles.

    Returns ``(experiment_id, provider, resolver, meta_provider)`` or ``None``
    if the orchestrator failed to start. Blocking (``run_audit`` scans data); in
    the front door this runs on a worker thread so the UI stays responsive.
    """
    from datetime import datetime
    from aria.runtime.experiment_view import build_snapshot

    started = orchestrator.run(experiment_id, context)
    if started.get("status") != "started":
        return None
    orchestrator.run_audit(experiment_id)

    start = datetime.now()

    def provider() -> ExperimentSnapshot:
        return build_snapshot(
            experiment_id,
            session=_session_for(orchestrator, experiment_id),
            start_time=start,
        )

    def meta_provider() -> dict:
        return _meta_for(orchestrator, experiment_id,
                         version=version, data_dir=context.get("data_dir"))

    return (experiment_id, provider, orchestrator.on_checkpoint_resolved,
            meta_provider)


def launch_cockpit(orchestrator: Any, experiment_id: str,
                   context: dict) -> Optional[ExperimentSnapshot]:
    """Drive a full run through the Textual cockpit (immediate-run path).

    Mirrors ``aria.tui.run_analysis`` but renders via the cockpit and resolves
    checkpoints through the read-model. Returns the final snapshot (or ``None``
    if the orchestrator failed to start).
    """
    from aria.version import __version__

    handles = _start_run_handles(orchestrator, experiment_id, context,
                                 __version__)
    if handles is None:
        return None
    experiment_id, provider, resolver, meta_provider = handles
    app = AriaCockpit(
        provider, resolver,
        experiment_id=experiment_id,
        version=__version__,
        meta_provider=meta_provider,
    )
    return app.run()


def run_control_center(orchestrator: Any, *,
                       version: str,
                       intake_kwargs: dict,
                       resolve_context: Callable[[Any], Optional[tuple]],
                       ) -> Optional[ExperimentSnapshot]:
    """Run the intake and the cockpit as ONE Textual app (no console flicker).

    The intake screen is the app's first screen; on submit, the run starts on a
    worker thread and the cockpit run view takes over in the same app. The
    transition never leaves the alternate screen, so the user does not see the
    invoking console between the two stages.

    ``resolve_context`` maps an ``IntakeResult`` to ``(experiment_id, context)``
    (or ``None``); it runs on the worker thread. Returns the final snapshot, or
    ``None`` if the user exited the intake or the run failed to start.
    """
    from aria.ui.intake import AriaIntakeScreen

    def starter(result: Any) -> Optional[tuple]:
        resolved = resolve_context(result)
        if not resolved:
            return None
        experiment_id, context = resolved
        return _start_run_handles(orchestrator, experiment_id, context, version)

    app = AriaCockpit(
        version=version,
        intake_screen_factory=lambda: AriaIntakeScreen(**intake_kwargs),
        run_starter=starter,
    )
    return app.run()
