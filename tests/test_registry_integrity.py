from pathlib import Path


def test_agent_registry_imports_and_script_contracts_are_valid():
    from aria.utils.registry_integrity import check_registry_integrity

    issues = check_registry_integrity(Path(__file__).resolve().parents[1])
    errors = [issue for issue in issues if issue.severity == "error"]

    assert errors == []


def test_every_dispatch_modality_has_validation_metadata():
    from aria.agents.orchestrator_agent import (
        MODALITY_TO_AGENT,
        MODALITY_VALIDATION,
    )

    missing = set(MODALITY_TO_AGENT) - set(MODALITY_VALIDATION)

    assert missing == set()
    assert all(MODALITY_VALIDATION[m].get("level") for m in MODALITY_TO_AGENT)


def test_scaffold_chromatin_modalities_are_not_dispatched():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION

    for modality in ("scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"):
        assert MODALITY_VALIDATION[modality]["level"] == "scaffold"
        assert MODALITY_VALIDATION[modality]["dispatch_enabled"] is False


def test_scaffold_integration_agent_is_not_dispatched():
    from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION

    assert INTEGRATION_VALIDATION["level"] == "scaffold"
    assert INTEGRATION_VALIDATION["dispatch_enabled"] is False
