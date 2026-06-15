"""U1 Textual cockpit shell — pilot tests (skipped if the `tui` extra is absent).

Drives the cockpit with a scripted snapshot provider + a fake checkpoint
resolver to assert that (a) a checkpoint decision routes back through the
resolver with the right option, and (b) the ledger view toggles. Uses
``asyncio.run`` so no pytest-asyncio dependency is required.
"""

from __future__ import annotations

import asyncio

import pytest

pytest.importorskip("textual")

from textual.app import App  # noqa: E402

from aria.runtime.experiment_view import (  # noqa: E402
    ArtifactView, ExperimentSnapshot, CheckpointView, ExperimentHistoryView,
)
from aria.ui.cockpit import AriaCockpit  # noqa: E402
from aria.ui.intake import (  # noqa: E402
    AriaIntakeApp, AriaIntakeScreen, IntakeResult,
)


def _snap(*, pending=None, done=False, ledger=None, artifacts=None,
          report_path=None) -> ExperimentSnapshot:
    return ExperimentSnapshot(
        experiment_id="e1",
        phase="audit" if pending else ("done" if done else "dispatch"),
        progress=0.0 if pending else 0.5,
        last_status=None,
        findings_by_confidence={"HIGH": [], "MEDIUM": [], "LOW": [],
                                "INSUFFICIENT": []},
        pending_checkpoint=pending,
        ledger=ledger or [],
        report_path=report_path,
        done=done,
        elapsed_s=1.0,
        silent_s=0.0,
        artifacts=artifacts or [],
    )


def test_checkpoint_choice_routes_to_resolver():
    cp = CheckpointView(message_id="m1", number=1, title="Data Audit Results",
                        question="Confirm modality?",
                        options=["Continue", "Correct"])
    calls: list[dict] = []

    state = {"snap": _snap(pending=cp)}

    def provider():
        return state["snap"]

    def resolver(*, message_id, user_decision, experiment_id):
        calls.append({"message_id": message_id,
                      "user_decision": user_decision,
                      "experiment_id": experiment_id})
        # After resolution, the next snapshot has no pending checkpoint.
        state["snap"] = _snap(pending=None)

    async def _run():
        app = AriaCockpit(provider, resolver, experiment_id="e1",
                          exit_on_done=False)
        async with app.run_test() as pilot:
            await pilot.pause()
            await pilot.press("2")        # choose option index 2 -> "Correct"
            await pilot.pause()

    asyncio.run(_run())

    assert len(calls) == 1
    assert calls[0]["user_decision"] == "Correct"
    assert calls[0]["message_id"] == "m1"
    assert calls[0]["experiment_id"] == "e1"


def test_choice_ignored_when_no_pending_checkpoint():
    calls: list[dict] = []

    def resolver(**kw):
        calls.append(kw)

    async def _run():
        app = AriaCockpit(lambda: _snap(pending=None), resolver,
                          experiment_id="e1", exit_on_done=False)
        async with app.run_test() as pilot:
            await pilot.pause()
            await pilot.press("1")
            await pilot.pause()

    asyncio.run(_run())
    assert calls == []


def test_design_editor_submits_json_to_callback():
    from textual.app import App

    from aria.ui.design_editor import DesignEditorScreen

    captured: list[str] = []

    class _Host(App):
        def __init__(self, screen):
            super().__init__()
            self._screen = screen

        def on_mount(self):
            self.push_screen(self._screen)

    screen = DesignEditorScreen({"treated": ["s1", "s2"], "control": ["s3"]},
                                lambda j: captured.append(j))

    async def _run():
        async with _Host(screen).run_test() as pilot:
            await pilot.pause()
            await pilot.press("c")     # row 0 = s1: treated -> control
            await pilot.press("s")     # submit
            await pilot.pause()

    asyncio.run(_run())

    assert len(captured) == 1
    import json
    groups = json.loads(captured[0])
    assert "s1" in groups["control"]      # s1 moved to control
    assert groups["treated"] == ["s2"]


def test_cockpit_e_opens_design_editor():
    from aria.ui.design_editor import DesignEditorScreen

    cp = CheckpointView(message_id="m1", number=2.1,
                        title="Experimental Groups",
                        question="Confirm groups?",
                        options=["Yes — confirm groups", "Cancel experiment"],
                        context={"proposed_groups": {"a": ["s1"], "b": ["s2"]}})

    async def _run():
        app = AriaCockpit(lambda: _snap(pending=cp), lambda **k: None,
                          experiment_id="e1", exit_on_done=False)
        async with app.run_test() as pilot:
            await pilot.pause()
            await pilot.press("e")
            await pilot.pause()
            assert isinstance(app.screen, DesignEditorScreen)

    asyncio.run(_run())


def test_ledger_toggle_shows_and_hides():
    async def _run():
        app = AriaCockpit(lambda: _snap(pending=None), lambda **k: None,
                          experiment_id="e1", exit_on_done=False)
        async with app.run_test() as pilot:
            await pilot.pause()
            from textual.widgets import Static
            ledger = app.query_one("#ledger", Static)
            findings = app.query_one("#findings", Static)
            assert ledger.display is False
            await pilot.press("l")
            await pilot.pause()
            assert ledger.display is True
            assert findings.display is False
            await pilot.press("l")
            await pilot.pause()
            assert ledger.display is False

    asyncio.run(_run())


def test_readiness_toggle_shows_and_hides():
    async def _run():
        app = AriaCockpit(lambda: _snap(pending=None), lambda **k: None,
                          experiment_id="e1", exit_on_done=False)
        async with app.run_test() as pilot:
            await pilot.pause()
            from textual.widgets import Static
            readiness = app.query_one("#readiness", Static)
            findings = app.query_one("#findings", Static)
            assert readiness.display is False
            await pilot.press("r")
            await pilot.pause()
            assert readiness.display is True
            assert findings.display is False
            await pilot.press("r")
            await pilot.pause()
            assert readiness.display is False
            assert findings.display is True

    asyncio.run(_run())


def test_resources_toggle_shows_and_hides():
    async def _run():
        app = AriaCockpit(lambda: _snap(pending=None), lambda **k: None,
                          experiment_id="e1", exit_on_done=False)
        async with app.run_test() as pilot:
            await pilot.pause()
            from textual.widgets import Static
            resources = app.query_one("#resources", Static)
            findings = app.query_one("#findings", Static)
            assert resources.display is False
            await pilot.press("u")
            await pilot.pause()
            assert resources.display is True
            assert findings.display is False
            await pilot.press("u")
            await pilot.pause()
            assert resources.display is False
            assert findings.display is True

    asyncio.run(_run())


def test_artifacts_toggle_shows_and_hides():
    async def _run():
        app = AriaCockpit(lambda: _snap(pending=None), lambda **k: None,
                          experiment_id="e1", exit_on_done=False)
        async with app.run_test() as pilot:
            await pilot.pause()
            from textual.widgets import Static
            artifacts = app.query_one("#artifacts", Static)
            findings = app.query_one("#findings", Static)
            assert artifacts.display is False
            await pilot.press("a")
            await pilot.pause()
            assert artifacts.display is True
            assert findings.display is False
            await pilot.press("a")
            await pilot.pause()
            assert artifacts.display is False
            assert findings.display is True

    asyncio.run(_run())


def test_cockpit_stays_open_on_done_by_default_and_shows_artifacts():
    artifacts = [
        ArtifactView(category="report", name="report.html", status="present",
                     detail="HTML report.", path="/tmp/report.html")
    ]

    async def _run():
        app = AriaCockpit(
            lambda: _snap(done=True, artifacts=artifacts,
                          report_path="/tmp/report.html"),
            lambda **k: None,
            experiment_id="e1",
        )
        async with app.run_test() as pilot:
            await pilot.pause()
            from textual.widgets import Static
            artifacts_widget = app.query_one("#artifacts", Static)
            findings = app.query_one("#findings", Static)
            agents = app.query_one("#agents", Static)
            assert app._exit_on_done is False
            assert app._done_seen is True
            assert app._center_mode == "artifacts"
            assert artifacts_widget.display is True
            assert findings.display is False
            assert agents.display is True

    asyncio.run(_run())


def test_cockpit_exit_on_done_remains_available_for_callers_that_request_it():
    async def _run():
        app = AriaCockpit(lambda: _snap(done=True), lambda **k: None,
                          experiment_id="e1", exit_on_done=True)
        async with app.run_test() as pilot:
            await pilot.pause()
            return app.return_value

    result = asyncio.run(_run())
    assert isinstance(result, ExperimentSnapshot)
    assert result.done is True


class _IntakeHost(App):
    """Hosts AriaIntakeScreen and captures its dismiss result WITHOUT exiting,
    so tests can inspect the screen and result freely."""

    def __init__(self, **kwargs):
        super().__init__()
        self._kwargs = kwargs
        self.result = "UNSET"
        self.intake = None

    def on_mount(self) -> None:
        self.intake = AriaIntakeScreen(**self._kwargs)
        self.push_screen(self.intake, self._captured)

    def _captured(self, result) -> None:
        self.result = result  # do not exit; let the test assert


async def _mounted_intake(app, pilot):
    """Wait until the intake screen is mounted and return it."""
    await pilot.pause()
    await pilot.pause()
    return app.intake


def test_intake_app_hosts_the_screen():
    """The standalone AriaIntakeApp (used by launch_intake) mounts the screen."""
    async def _run():
        app = AriaIntakeApp(version="4.6.1")
        async with app.run_test() as pilot:
            await pilot.pause()
            await pilot.pause()
            return type(app.screen).__name__

    assert asyncio.run(_run()) == "AriaIntakeScreen"


def test_intake_submits_data_and_question(tmp_path):
    async def _run():
        app = _IntakeHost(version="4.6.1")
        async with app.run_test() as pilot:
            scr = await _mounted_intake(app, pilot)
            from textual.widgets import Input, TextArea
            scr.query_one("#data-input", Input).value = str(tmp_path)
            scr.query_one("#question-input", TextArea).load_text(
                "Which cell types differ?"
            )
            scr.action_start()
            return scr.submitted

    assert asyncio.run(_run()) == IntakeResult(
        data_input=str(tmp_path),
        question="Which cell types differ?",
    )


def test_intake_requires_question(tmp_path):
    async def _run():
        app = _IntakeHost(version="4.6.1")
        async with app.run_test() as pilot:
            scr = await _mounted_intake(app, pilot)
            from textual.widgets import Input
            scr.query_one("#data-input", Input).value = str(tmp_path)
            scr.action_start()
            return scr.submitted, scr.status_message

    submitted, status = asyncio.run(_run())
    assert submitted is None
    assert "biological question" in status


def test_intake_pastes_path_into_data_input_once(tmp_path):
    """A terminal Paste into the focused path Input lands once (no double)."""
    from textual import events

    async def _run():
        app = _IntakeHost(version="4.6.1")
        async with app.run_test() as pilot:
            scr = await _mounted_intake(app, pilot)
            from textual.widgets import Input
            data_input = scr.query_one("#data-input", Input)
            app.set_focus(data_input)
            await pilot.pause()
            # Post the Paste the way the driver does: to the app, which forwards
            # it to the focused widget's native paste handler.
            app.post_message(events.Paste(text=str(tmp_path)))
            await pilot.pause()
            await pilot.pause()
            return data_input.value

    assert asyncio.run(_run()) == str(tmp_path)


def test_intake_pastes_multiline_question_once():
    """A multi-line Paste into the question TextArea is not double-inserted."""
    from textual import events

    async def _run():
        app = _IntakeHost(version="4.6.1")
        async with app.run_test() as pilot:
            scr = await _mounted_intake(app, pilot)
            from textual.widgets import TextArea
            question = scr.query_one("#question-input", TextArea)
            app.set_focus(question)
            await pilot.pause()
            app.post_message(events.Paste(text="line one\nline two"))
            await pilot.pause()
            await pilot.pause()
            return question.text

    # Exactly once: regression guard for the TextArea double-paste bug.
    assert asyncio.run(_run()) == "line one\nline two"


def test_intake_ctrl_v_pastes_os_clipboard_into_path(monkeypatch):
    """Ctrl+V reads the OS clipboard (not Textual's empty in-app clipboard).

    Regression guard for the 'nothing pastes on Ctrl+V' bug: terminals that pass
    Ctrl+V through trigger `action_paste`, which must read the OS clipboard.
    """
    from aria.ui import clipboard
    monkeypatch.setattr(clipboard, "read_clipboard",
                        lambda *a, **k: "/data/exp1\nignored second line")

    async def _run():
        app = _IntakeHost(version="4.6.1")
        async with app.run_test() as pilot:
            scr = await _mounted_intake(app, pilot)
            from textual.widgets import Input
            data_input = scr.query_one("#data-input", Input)
            app.set_focus(data_input)
            await pilot.pause()
            await pilot.press("ctrl+v")
            await pilot.pause()
            return data_input.value

    # Path field keeps the first line only.
    assert asyncio.run(_run()) == "/data/exp1"


def test_intake_invalid_data_keeps_app_open():
    """A failing data validation stays on the intake instead of dismissing."""
    def _validator(value: str):
        return None if value == "OK" else "Path not found: " + value

    async def _run():
        app = _IntakeHost(version="4.6.1", data_validator=_validator)
        async with app.run_test() as pilot:
            scr = await _mounted_intake(app, pilot)
            from textual.widgets import Input, TextArea
            scr.query_one("#data-input", Input).value = "/bad/path"
            scr.query_one("#question-input", TextArea).load_text("Q?")
            scr.action_start()
            await pilot.pause()
            return app.result, scr.submitted, scr.status_message

    result, submitted, status = asyncio.run(_run())
    assert result == "UNSET"                  # never dismissed
    assert submitted is None
    assert "Path not found" in status         # told the user why


def test_intake_valid_data_submits_with_validator():
    def _validator(value: str):
        return None if value == "OK" else "bad"

    async def _run():
        app = _IntakeHost(version="4.6.1", data_validator=_validator)
        async with app.run_test() as pilot:
            scr = await _mounted_intake(app, pilot)
            from textual.widgets import Input, TextArea
            scr.query_one("#data-input", Input).value = "OK"
            scr.query_one("#question-input", TextArea).load_text("Q?")
            scr.action_start()
            return scr.submitted

    assert asyncio.run(_run()) == IntakeResult(data_input="OK", question="Q?")


def test_intake_ctrl_v_pastes_os_clipboard_into_question(monkeypatch):
    from aria.ui import clipboard
    monkeypatch.setattr(clipboard, "read_clipboard",
                        lambda *a, **k: "q line 1\nq line 2")

    async def _run():
        app = _IntakeHost(version="4.6.1")
        async with app.run_test() as pilot:
            scr = await _mounted_intake(app, pilot)
            from textual.widgets import TextArea
            question = scr.query_one("#question-input", TextArea)
            app.set_focus(question)
            await pilot.pause()
            await pilot.press("ctrl+v")
            await pilot.pause()
            return question.text

    # Question field keeps the full multi-line text.
    assert asyncio.run(_run()) == "q line 1\nq line 2"


def test_intake_renders_resume_history():
    history = [
        ExperimentHistoryView(
            experiment_id="abcd1234efgh", name="H9 RNA timecourse",
            organism="Homo sapiens", genome="hg38",
            updated_at="2026-06-10T09:30:00", summary="DEGs found.",
            n_decisions=3, last_decision="CP2.6: confirm",
            modalities=["bulk_RNA"], report_path="/r/report.html",
            has_report=True),
    ]

    async def _run():
        app = _IntakeHost(version="4.6.1", history=history)
        async with app.run_test() as pilot:
            scr = await _mounted_intake(app, pilot)
            from rich.console import Console
            # history -> render_history -> Panel into #experiments.
            renderable = scr._experiments_renderable()
            console = Console(width=80, record=True)
            with console.capture() as cap:
                console.print(renderable)
            return cap.get()

    out = asyncio.run(_run())
    assert "H9 RNA timecourse" in out
    assert "report on disk" in out


def test_unified_front_door_transitions_intake_to_run(monkeypatch):
    """The intake screen and the run view live in ONE app (no console flash)."""
    from datetime import datetime
    from aria.ui.cockpit import AriaCockpit
    from aria.runtime.experiment_view import ExperimentSnapshot, ProgressEvent

    def _snap_provider():
        return ExperimentSnapshot(
            experiment_id="EXP1", phase="audit", progress=0.0,
            last_status=ProgressEvent(ts=datetime.now(), sender="orchestrator",
                                      text="scanning", progress=0.0),
            findings_by_confidence={"HIGH": [], "MEDIUM": [], "LOW": [],
                                    "INSUFFICIENT": []},
            pending_checkpoint=None, ledger=[], report_path=None, done=False,
            elapsed_s=1.0, silent_s=0.0)

    started = {}

    def starter(result):
        started["result"] = result
        return ("EXP1", _snap_provider, lambda **k: None, lambda: {})

    async def _run():
        app = AriaCockpit(
            version="4.6.1", exit_on_done=False,
            intake_screen_factory=lambda: AriaIntakeScreen(version="4.6.1"),
            run_starter=starter,
        )
        async with app.run_test() as pilot:
            await pilot.pause()
            await pilot.pause()
            from textual.widgets import Input, TextArea
            assert type(app.screen).__name__ == "AriaIntakeScreen"
            app.screen.query_one("#data-input", Input).value = "/tmp"
            app.screen.query_one("#question-input", TextArea).load_text("Q?")
            app.screen.action_start()
            await pilot.pause()
            await app.workers.wait_for_complete()
            await pilot.pause()
            return started.get("result"), app.experiment_id, app._provider is not None

    result, exp_id, provider_set = asyncio.run(_run())
    assert result == IntakeResult(data_input="/tmp", question="Q?")
    assert exp_id == "EXP1"           # transitioned to the run within one app
    assert provider_set is True
