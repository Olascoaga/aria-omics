"""ARIA cockpit — Textual presentation shell (U1) with the run-ledger view (U5).

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


class AriaCockpit(App):
    """Cockpit app. Pure presentation over the read-model."""

    CSS = """
    #body { height: 1fr; }
    #left { width: 34; }
    #right { width: 46; }
    #center { width: 1fr; }
    Static { height: auto; }
    #findings { height: 1fr; }
    #ledger { height: 1fr; }
    """

    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("l", "toggle_ledger", "Ledger"),
        Binding("e", "edit_design", "Edit groups"),
        Binding("1", "choose(1)", "Opt 1", show=False),
        Binding("2", "choose(2)", "Opt 2", show=False),
        Binding("3", "choose(3)", "Opt 3", show=False),
        Binding("4", "choose(4)", "Opt 4", show=False),
        Binding("5", "choose(5)", "Opt 5", show=False),
        Binding("6", "choose(6)", "Opt 6", show=False),
    ]

    def __init__(self, snapshot_provider: SnapshotProvider,
                 checkpoint_resolver: CheckpointResolver, *,
                 experiment_id: str = "",
                 version: str = "",
                 meta_provider: Optional[MetaProvider] = None,
                 poll_interval: float = 0.5,
                 exit_on_done: bool = True):
        super().__init__()
        self._provider = snapshot_provider
        self._resolver = checkpoint_resolver
        self.experiment_id = experiment_id
        self._version = version
        self._meta_provider = meta_provider or (lambda: {})
        self._poll_interval = poll_interval
        self._exit_on_done = exit_on_done
        self._snap: Optional[ExperimentSnapshot] = None
        self._show_ledger = False

    def compose(self) -> ComposeResult:
        yield Horizontal(
            Vertical(Static(id="run-header"), id="left"),
            Vertical(
                Static(id="timeline"),
                Static(id="findings"),
                Static(id="ledger"),
                id="center",
            ),
            Vertical(Static(id="checkpoint"), id="right"),
            id="body",
        )
        yield Footer()

    def on_mount(self) -> None:
        self.query_one("#ledger", Static).display = False
        self.refresh_snapshot()
        self.set_interval(self._poll_interval, self.refresh_snapshot)

    # ── Rendering ────────────────────────────────────────────────────────────
    def refresh_snapshot(self) -> None:
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
        self.query_one("#findings", Static).update(render.render_findings(snap))
        self.query_one("#ledger", Static).update(render.render_ledger(snap))
        self.query_one("#checkpoint", Static).update(
            render.render_checkpoint(snap))

        if snap.done and self._exit_on_done:
            self.exit(snap)

    # ── Actions ──────────────────────────────────────────────────────────────
    def action_toggle_ledger(self) -> None:
        self._show_ledger = not self._show_ledger
        self.query_one("#ledger", Static).display = self._show_ledger
        self.query_one("#findings", Static).display = not self._show_ledger

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


def launch_cockpit(orchestrator: Any, experiment_id: str,
                   context: dict) -> Optional[ExperimentSnapshot]:
    """Drive a full run through the Textual cockpit.

    Mirrors ``aria.tui.run_analysis`` but renders via the cockpit and resolves
    checkpoints through the read-model. Returns the final snapshot (or ``None``
    if the orchestrator failed to start).
    """
    from datetime import datetime

    from aria.runtime.experiment_view import build_snapshot
    from aria.version import __version__

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

    app = AriaCockpit(
        provider,
        orchestrator.on_checkpoint_resolved,
        experiment_id=experiment_id,
        version=__version__,
        meta_provider=lambda: _meta_for(
            orchestrator, experiment_id,
            version=__version__, data_dir=context.get("data_dir"),
        ),
    )
    return app.run()
