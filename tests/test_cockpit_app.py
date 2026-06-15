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

from aria.runtime.experiment_view import (  # noqa: E402
    ExperimentSnapshot, CheckpointView,
)
from aria.ui.cockpit import AriaCockpit  # noqa: E402


def _snap(*, pending=None, done=False, ledger=None) -> ExperimentSnapshot:
    return ExperimentSnapshot(
        experiment_id="e1",
        phase="audit" if pending else ("done" if done else "dispatch"),
        progress=0.0 if pending else 0.5,
        last_status=None,
        findings_by_confidence={"HIGH": [], "MEDIUM": [], "LOW": [],
                                "INSUFFICIENT": []},
        pending_checkpoint=pending,
        ledger=ledger or [],
        report_path=None,
        done=done,
        elapsed_s=1.0,
        silent_s=0.0,
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
