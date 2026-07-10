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


# ── C2 · evidence verification is lexical + negation-blind ───────────────────
# aria/agents/narrative/evidence_verifier.py:65 — Title-case pure-alpha symbols
# (Gfap) fall outside the gene regex, and "No significant genes" is not read as
# negating the claim, so a fabricated positive claim is marked supported.
@pytest.mark.xfail(strict=True, reason="C2 open: verifier is lexical/negation-blind; "
                   "Gfap claim vs a no-hits evidence card passes (FASE 6)")
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


# ── A1 · air-gap resolved once at provider construction ──────────────────────
@pytest.mark.skip(reason="A1 probe owned by FASE 1: needs a per-run egress context; "
                  "air-gap enabled after LLMProvider() must still block every call "
                  "(aria/llm/provider.py:180)")
def test_probe_a1_airgap_after_construction_blocks_egress():
    ...


# ── A3 · global MessageBus leaks across concurrent runs ──────────────────────
@pytest.mark.skip(reason="A3 probe owned by FASE 1: per-run bus/provider isolation; "
                  "messages from run A must not appear in run B's LLM log "
                  "(aria/bus/message_bus.py:91)")
def test_probe_a3_concurrent_runs_do_not_cross_contaminate():
    ...


# ── E2 · filename R1/R2 rule routes ATAC FASTQ to bulk RNA ───────────────────
@pytest.mark.skip(reason="E2 probe owned by FASE 4: mandatory library-type manifest; "
                  "ATAC/scATAC FASTQ must not classify as bulk_RNA on filename order "
                  "(aria/agents/data_audit_agent.py:43)")
def test_probe_e2_atac_fastq_not_classified_as_bulk_rna():
    ...
