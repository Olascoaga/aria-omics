"""P0-3 regression: Hi-C must not dispatch by default.

Hi-C was a `scaffold` modality with `dispatch_enabled=True` (ADR-012 tracked
exception), so `genome_arch_agent` could publish TADs/loops at MEDIUM/HIGH
confidence — unvalidated conclusions presented as publication-grade. Per the
master plan (P0-3) Hi-C is now dispatch-OFF by default; running it requires the
explicit `ARIA_ALLOW_EXPERIMENTAL_HIC=1` opt-in and is stamped experimental /
not publication-grade. This does not add or change any Hi-C science.
"""

import pytest

from aria.agents.orchestrator_agent import (
    MODALITY_VALIDATION,
    OrchestratorAgent,
)


def test_hic_is_a_dispatch_gated_scaffold_by_default():
    meta = MODALITY_VALIDATION["HiC"]
    assert meta["level"] == "scaffold"
    assert meta["dispatch_enabled"] is False
    assert meta.get("experimental_env_flag") == "ARIA_ALLOW_EXPERIMENTAL_HIC"


def test_hic_is_blocked_when_flag_unset(monkeypatch):
    monkeypatch.delenv("ARIA_ALLOW_EXPERIMENTAL_HIC", raising=False)
    blocked = OrchestratorAgent._blocked_modalities({"HiC": ["x.hic"]})
    assert "HiC" in blocked
    experimental = OrchestratorAgent._experimental_modalities({"HiC": ["x.hic"]})
    assert experimental == {}


def test_hic_experimental_flag_unblocks_but_marks_experimental(monkeypatch):
    monkeypatch.setenv("ARIA_ALLOW_EXPERIMENTAL_HIC", "1")
    # Under the explicit opt-in Hi-C is no longer hard-blocked...
    blocked = OrchestratorAgent._blocked_modalities({"HiC": ["x.hic"]})
    assert "HiC" not in blocked
    # ...but it is recorded as an experimental, not-publication-grade run.
    experimental = OrchestratorAgent._experimental_modalities({"HiC": ["x.hic"]})
    assert "HiC" in experimental


def test_other_scaffolds_ignore_the_hic_flag(monkeypatch):
    """The Hi-C opt-in must not unblock unrelated scaffold modalities.

    T10 (tri-audit 2026-06-14): scATAC is no longer a scaffold (it is alpha +
    dispatch-on-ack since ADR-033), so this regression must use a still-scaffolded
    modality. `bulk_ATAC`/`ChIP` are scaffold + dispatch_enabled=False and have no
    experimental opt-in, so the Hi-C flag must leave them blocked.
    """
    monkeypatch.setenv("ARIA_ALLOW_EXPERIMENTAL_HIC", "1")
    blocked = OrchestratorAgent._blocked_modalities(
        {"bulk_ATAC": ["a.bed"], "ChIP": ["b.bam"]}
    )
    assert "bulk_ATAC" in blocked
    assert "ChIP" in blocked
    # The Hi-C opt-in only ever affects Hi-C.
    experimental = OrchestratorAgent._experimental_modalities(
        {"bulk_ATAC": ["a.bed"], "ChIP": ["b.bam"]}
    )
    assert experimental == {}
