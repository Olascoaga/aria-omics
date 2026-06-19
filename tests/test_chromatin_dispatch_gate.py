def test_orchestrator_allows_beta_scatac_after_readiness_ack_gate():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION, OrchestratorAgent

    blocked = OrchestratorAgent._blocked_modalities({
        "scRNA": ["/tmp/rna.h5ad"],
        "scATAC": ["/tmp/fragments.tsv.gz"],
    })

    assert blocked == {}
    # De-alpha 2026-06-15 (ADR-048): scATAC is beta but still requires_ack.
    assert MODALITY_VALIDATION["scATAC"]["level"] == "beta"
    assert MODALITY_VALIDATION["scATAC"]["dispatch_enabled"] is True


def test_orchestrator_allows_beta_bulk_atac_after_readiness_ack_gate():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION, OrchestratorAgent

    blocked = OrchestratorAgent._blocked_modalities({
        "bulk_ATAC": ["/tmp/fragments.tsv.gz"],
    })

    assert blocked == {}
    assert MODALITY_VALIDATION["bulk_ATAC"]["level"] == "beta"
    assert MODALITY_VALIDATION["bulk_ATAC"]["dispatch_enabled"] is True


def test_chromatin_missing_planned_scripts_return_structured_blocker():
    from aria.agents.chromatin_agent import ChromatinAgent

    motif = ChromatinAgent._planned_script_blocker(
        "aria/scripts/chromatin_motifs.py",
        "motif_enrichment",
    )

    assert motif["status"] == "skipped"
    assert motif["reason"] == "script_not_implemented"
    assert motif["validation_level"] == "scaffold"
