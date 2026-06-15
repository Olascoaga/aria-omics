"""ARIA Control Center — Textual presentation layer (U1+).

This subpackage is pure presentation. It renders the U0 read-model
(:mod:`aria.runtime.experiment_view`) and never touches the scientific core,
the orchestrator's decision logic, or the governance walls. The headless runner
(:mod:`aria.headless`) remains the canonical, reproducible, non-interactive path.

``render`` holds Rich-only pure renderers (importable without Textual);
``cockpit`` holds the Textual app + launcher (requires the ``tui`` extra).
"""
