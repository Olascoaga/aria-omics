"""Preprint-readiness audit — Q4: reproduced-bug probe harness (TDD anchor).

Codex reproduced five concrete probes. This module encodes them as tests that
assert the DESIRED (fixed) behaviour and are marked ``xfail(strict=True)`` while
the bug is open:

  - open bug  -> the desired-behaviour assertion fails -> xfail (green suite);
  - once the fix lands -> the assertion passes -> strict xfail turns the xpass
    into a FAILURE, forcing whoever fixes it to drop the marker and promote the
    probe to a normal regression test.

Two probes are pure-Python and reproduced airtight here (B6, C2). Three depend on
runtime wiring (air-gap A1, MessageBus A3, front-door classification E2) and are
kept as documented ``skip`` placeholders so the harness lists all five; the FASE 1
/ FASE 4 slices that own those fixes convert them into live ratchets.

Tracker: memory/audit/ARIA_PLAN_AUDITORIA_preprint_journal_2026-07-09.md
"""
from __future__ import annotations

import os

import numpy as np
import pytest


# ── B6 · count classifier only inspects a 200-row sample ────────────────────
# aria/utils/count_classifier.py:33 (sample_row_indices) + rna_bulk_de.py:873.
# A fractional value OUTSIDE the sampled slice is never seen, so a non-count
# matrix is accepted as raw counts and later rounded into DESeq2.
@pytest.mark.xfail(strict=True, reason="B6 open: classify_matrix samples ~200 rows; "
                   "a fractional value outside the sample is missed (FASE 3)")
def test_probe_b6_fractional_row_outside_sample_is_caught():
    from aria.utils.count_classifier import classify_matrix, sample_row_indices

    n = 5000
    mat = np.full((n, 6), 10.0)                      # all-integer counts
    sampled = set(sample_row_indices(n, seed=0).tolist())
    victim = next(i for i in range(n) if i not in sampled)
    mat[victim, 0] = 3.5                             # fractional, unsampled

    res = classify_matrix(mat, seed=0)
    # DESIRED: any fractional value anywhere blocks the raw-counts verdict.
    assert res["is_raw_counts"] is False


# ── C2 · semantic facts + explicit negation (FIXED FASE 6) ───────────────────
# Title-case pure-alpha symbols are resolved from the typed DE predicate frame,
# and explicit null evidence contradicts an affirmed entity-level DE claim.
def test_probe_c2_gfap_claim_against_null_evidence_is_unsupported():
    from aria.agents.narrative.types import NarrativeBlock, EvidenceItem
    from aria.agents.narrative.evidence_verifier import verify_block_claim_support

    block = NarrativeBlock(
        id="probe.c2", modality="bulk_RNA", analysis="differential_expression",
        block_type="result", title="probe", status="success", confidence="medium",
        claim="Gfap was differentially expressed between conditions.",
        evidence=[EvidenceItem(label="Result",
                               value="No significant genes were detected",
                               source="pydeseq2")],
    )
    manifest = verify_block_claim_support(block, strict=False)
    # DESIRED: the fabricated positive claim is NOT supported by null evidence.
    assert manifest["status"] == "unsupported"


# ── A1 · air-gap resolved once at provider construction (FIXED FASE 1) ───────
# aria/llm/provider.py — air-gap is no longer snapshotted in __init__; the
# _air_gapped property re-evaluates live, so enabling air-gap AFTER the provider
# is built still refuses every cloud call with zero egress.
def test_probe_a1_airgap_after_construction_blocks_egress(monkeypatch):
    monkeypatch.delenv("ARIA_AIR_GAPPED", raising=False)
    from aria.llm import provider as prov_mod
    from aria.llm.provider import LLMProvider, ModelConfig, TaskTier
    from aria.utils import privacy

    # Restore global air-gap state after the test (enable_* mutates os.environ).
    prev = os.environ.get("ARIA_AIR_GAPPED")
    prev_reason = privacy._runtime_air_gapped_reason
    try:
        # Any real network call fails loudly -> proves zero egress.
        def _boom(*a, **k):
            raise AssertionError("EGRESS: litellm.completion was called")
        monkeypatch.setattr(prov_mod.litellm, "completion", _boom)

        cloud_only = {TaskTier.MEDIUM:
                      [ModelConfig("anthropic", "claude-cloud", 8000, is_local=False)]}
        p = LLMProvider(models=cloud_only)          # built NOT air-gapped
        assert p._air_gapped is False

        privacy.enable_air_gapped_runtime(reason="post_construction_optin")
        assert p._air_gapped is True                # live re-check (the A1 fix)

        # The cloud call is refused BEFORE any egress, not attempted.
        with pytest.raises(RuntimeError, match="AIR_GAPPED"):
            p.complete(prompt="hello", tier=TaskTier.MEDIUM)
    finally:
        if prev is None:
            os.environ.pop("ARIA_AIR_GAPPED", None)
        else:
            os.environ["ARIA_AIR_GAPPED"] = prev
        privacy._runtime_air_gapped_reason = prev_reason


# ── A3 · global MessageBus leaks across concurrent runs ──────────────────────
def test_probe_a3_concurrent_runs_do_not_cross_contaminate(tmp_path):
    from aria.bus.message_bus import Message, MessageBus, MessageType, bus
    from aria.runtime.experiment_session import ExperimentSession

    a = ExperimentSession(
        "probe-A", message_bus=MessageBus(persist_path=tmp_path / "A.jsonl")
    )
    b = ExperimentSession(
        "probe-B", message_bus=MessageBus(persist_path=tmp_path / "B.jsonl")
    )
    try:
        bus.publish(Message(
            sender="probe", receiver="orchestrator", type=MessageType.STATUS,
            payload={"run": "A"}, experiment_id="probe-A",
        ))
        bus.publish(Message(
            sender="probe", receiver="orchestrator", type=MessageType.STATUS,
            payload={"run": "B"}, experiment_id="probe-B",
        ))
        assert [m.payload["run"] for m in a.message_bus.get_log()] == ["A"]
        assert [m.payload["run"] for m in b.message_bus.get_log()] == ["B"]
    finally:
        bus.unbind_experiment("probe-A")
        bus.unbind_experiment("probe-B")


# ── E2 · filename R1/R2 rule routes ATAC FASTQ to bulk RNA ───────────────────
@pytest.mark.skip(reason="E2 probe owned by FASE 4: mandatory library-type manifest; "
                  "ATAC/scATAC FASTQ must not classify as bulk_RNA on filename order "
                  "(aria/agents/data_audit_agent.py:43)")
def test_probe_e2_atac_fastq_not_classified_as_bulk_rna():
    ...
