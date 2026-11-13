"""ADR-081 & ADR-082 promotion + OrchestratorAgent scaffold/experimental gate guard.

Hi-C (`ADR-082`) and ChIP / CUT&RUN / CUT&TAG (`ADR-081`) are now promoted to
`production` (`dispatch_enabled=True`). We verify both their production status
in `MODALITY_VALIDATION` and that `OrchestratorAgent._blocked_modalities` /
`_experimental_modalities` continue to fail closed on any hypothetical future
`scaffold` modality unless explicitly opted in via its own `experimental_env_flag`.
"""

import pytest

from aria.agents.orchestrator_agent import (
    MODALITY_VALIDATION,
    OrchestratorAgent,
)


def test_hic_and_epigenomics_are_production_and_unblocked_by_default(monkeypatch):
    monkeypatch.delenv("ARIA_ALLOW_EXPERIMENTAL_HIC", raising=False)
    meta = MODALITY_VALIDATION["HiC"]
    assert meta["level"] == "production"
    assert meta["dispatch_enabled"] is True
    assert meta["promotion_adr"] == "ADR-082"

    for mod in ("ChIP", "CUT_AND_RUN", "CUT_AND_TAG"):
        assert MODALITY_VALIDATION[mod]["level"] == "production"
        assert MODALITY_VALIDATION[mod]["dispatch_enabled"] is True
        assert MODALITY_VALIDATION[mod]["promotion_adr"] == "ADR-081"

    assert OrchestratorAgent._blocked_modalities(
        {"HiC": ["x.mcool"], "ChIP": ["b.bam"], "CUT_AND_RUN": ["c.bam"], "CUT_AND_TAG": ["d.bam"]}
    ) == {}
    assert OrchestratorAgent._experimental_modalities(
        {"HiC": ["x.mcool"], "ChIP": ["b.bam"], "CUT_AND_RUN": ["c.bam"], "CUT_AND_TAG": ["d.bam"]}
    ) == {}


def test_scaffold_gate_blocks_when_flag_unset_and_marks_experimental_when_set(monkeypatch):
    monkeypatch.setitem(
        MODALITY_VALIDATION,
        "SyntheticScaffold",
        {
            "level": "scaffold",
            "dispatch_enabled": False,
            "experimental_env_flag": "ARIA_ALLOW_EXPERIMENTAL_SYNTH",
        },
    )
    monkeypatch.setitem(
        MODALITY_VALIDATION,
        "UnflaggedScaffold",
        {"level": "scaffold", "dispatch_enabled": False},
    )

    monkeypatch.delenv("ARIA_ALLOW_EXPERIMENTAL_SYNTH", raising=False)
    assert "SyntheticScaffold" in OrchestratorAgent._blocked_modalities({"SyntheticScaffold": ["x.dat"]})
    assert OrchestratorAgent._experimental_modalities({"SyntheticScaffold": ["x.dat"]}) == {}

    monkeypatch.setenv("ARIA_ALLOW_EXPERIMENTAL_SYNTH", "1")
    assert "SyntheticScaffold" not in OrchestratorAgent._blocked_modalities({"SyntheticScaffold": ["x.dat"]})
    assert "SyntheticScaffold" in OrchestratorAgent._experimental_modalities({"SyntheticScaffold": ["x.dat"]})

    # Opt-in flag never leaks to unrelated scaffold entries
    assert "UnflaggedScaffold" in OrchestratorAgent._blocked_modalities({"UnflaggedScaffold": ["y.dat"]})
    assert OrchestratorAgent._experimental_modalities({"UnflaggedScaffold": ["y.dat"]}) == {}
