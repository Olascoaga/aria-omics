"""Preprint-readiness audit B3-multi (FASE 7 / Claim 6): the heavy, opt-in
end-to-end lane over ``run_headless``.

Fase A (``test_preprint_audit_b3_multi_provider.py``) proved the invariant at the
claim-compiler/diff level in-process. This module covers the two pieces that
only the real runner can prove:

1. Provider-injection plumbing: ``run_headless`` forwards an injected
   ``LLMProvider`` to the orchestrator, so a run can be driven under a chosen
   (or scripted, no-egress) provider without editing scientific computations.
   These guards are lightweight and run in ordinary CI.

2. A dataset-gated end-to-end reproducibility check: the SAME analysis, driven
   twice through the real runner via the injected-provider seam, must yield an
   equivalent reproducibility artifact — identical deterministic statistics and
   public accepted claims, with only LLM prose allowed to differ. This is the
   end-to-end counterpart of Fase A's in-process provider-variation proof: Fase
   A holds provider prose variable and shows the public boundary is invariant;
   this run holds the pipeline fixed and shows repeated real executions agree
   over the run ledger, claim set and reproducibility capsule.

   It is skipped unless a dataset dir is provided via ``ARIA_B3_E2E_DATA``. A
   full-pipeline no-egress provider double is deliberately NOT used: several
   agents parse structured JSON from the model, so a canned-prose double cannot
   faithfully drive the run. The heavy lane therefore uses the configured
   provider (``LLMProvider.from_config``) and is opt-in, not part of ordinary CI.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

# This module runs in the narrative-kernel CI lane, which installs the package
# (litellm included), so no litellm import stub is used: the plumbing guards
# monkeypatch the orchestrator away, and the heavy E2E needs the real backend.


# ── 1. Provider-injection plumbing (lightweight, CI) ──────────────────────────

def _capture_orch_llm(monkeypatch, **kwargs):
    """Drive ``run_headless`` just far enough to capture the ``llm`` argument it
    hands ``OrchestratorAgent``, short-circuiting before any heavy work."""
    import aria.headless as headless

    captured: dict = {}

    class _StubOrch:
        def __init__(self, *a, **k):
            captured["llm"] = k.get("llm", "MISSING")

        def run(self, experiment_id, ctx):
            return {"status": "not_started"}  # makes run_headless return early

    class _StubMemory:
        def close(self):
            pass

    monkeypatch.setattr(headless, "OrchestratorAgent", _StubOrch, raising=False)
    monkeypatch.setattr(
        "aria.agents.orchestrator_agent.OrchestratorAgent", _StubOrch,
        raising=False,
    )
    monkeypatch.setattr(
        "aria.memory.memory.ARIAMemory", _StubMemory, raising=False
    )
    headless.run_headless("/data", "why?", **kwargs)
    return captured


def test_run_headless_forwards_injected_provider(monkeypatch):
    sentinel = object()
    captured = _capture_orch_llm(monkeypatch, llm_provider=sentinel)
    assert captured["llm"] is sentinel


def test_run_headless_defaults_to_no_injected_provider(monkeypatch):
    captured = _capture_orch_llm(monkeypatch)
    assert captured["llm"] is None


# ── 2. Dataset-gated end-to-end repetition reproducibility ────────────────────

_E2E_DATA = os.environ.get("ARIA_B3_E2E_DATA")


@pytest.mark.skipif(
    not _E2E_DATA,
    reason="set ARIA_B3_E2E_DATA to a dataset dir to run the B3 heavy E2E",
)
def test_repeated_real_runs_yield_equivalent_public_output(tmp_path):
    from aria.headless import run_headless
    from aria.llm.provider import LLMProvider
    from aria.agents.narrative.ledger_export import (
        diff_methodologies,
        load_methodology,
        verify_reproducible_capsule,
        write_reproducible_capsule,
    )

    question = "Compare the conditions and report differential expression."

    def _run():
        # Exercise the injected-provider seam end to end with the configured
        # provider. Determinism comes from reproducible_mode + temperature 0 +
        # the deterministic prompt cache, not from a scripted double.
        return run_headless(
            _E2E_DATA, question, reproducible_mode=True,
            llm_provider=LLMProvider.from_config(),
        )

    run_a = _run()
    run_b = _run()
    assert run_a.status == "completed", run_a
    assert run_b.status == "completed", run_b

    meth_a = load_methodology(Path(run_a.report_path).parent)
    meth_b = load_methodology(Path(run_b.report_path).parent)

    diff = diff_methodologies(meth_a, meth_b)
    # Deterministic statistics + public accepted claims must be invariant across
    # repeated real executions; only prose (never compared by diff_methodologies)
    # may differ.
    assert diff["claims"]["added"] == []
    assert diff["claims"]["removed"] == []
    assert diff["claims"]["changed"] == []
    assert diff["ledger"]["status_changed"] == []

    # Each run's reproducibility capsule verifies independently.
    for run in (run_a, run_b):
        report_dir = Path(run.report_path).parent
        capsule = write_reproducible_capsule(
            report_dir, tmp_path / f"{report_dir.name}_capsule.zip"
        )
        assert verify_reproducible_capsule(capsule)["status"] == "pass"
