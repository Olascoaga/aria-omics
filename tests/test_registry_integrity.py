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


def test_scatac_production_is_dispatchable_with_adr076():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION

    # ADR-076: scATAC promoted to production.
    assert MODALITY_VALIDATION["scATAC"]["level"] == "production"
    assert MODALITY_VALIDATION["scATAC"]["dispatch_enabled"] is True
    assert MODALITY_VALIDATION["scATAC"]["promotion_adr"] == "ADR-076"


def test_bulk_atac_production_is_dispatchable_with_adr077():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION

    # ADR-077: bulk_ATAC promoted to production.
    assert MODALITY_VALIDATION["bulk_ATAC"]["level"] == "production"
    assert MODALITY_VALIDATION["bulk_ATAC"]["dispatch_enabled"] is True
    assert MODALITY_VALIDATION["bulk_ATAC"]["promotion_adr"] == "ADR-077"


def test_remaining_scaffold_chromatin_modalities_are_not_dispatched():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION

    for modality in ("ChIP", "CUT_AND_RUN", "CUT_AND_TAG"):
        assert MODALITY_VALIDATION[modality]["level"] == "scaffold"
        assert MODALITY_VALIDATION[modality]["dispatch_enabled"] is False


def test_paired_integration_agent_is_production_for_same_cell_scope():
    from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION

    assert INTEGRATION_VALIDATION["level"] == "production"
    assert INTEGRATION_VALIDATION["dispatch_enabled"] is True
    assert INTEGRATION_VALIDATION["requires_ack"] is False
    assert INTEGRATION_VALIDATION["production_promoted"] is True
    assert INTEGRATION_VALIDATION["promotion_adr"] == "ADR-079"
    assert INTEGRATION_VALIDATION["scope"] == "paired_same_cell_rna_atac"
