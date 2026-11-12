def test_orchestrator_allows_production_scatac_auto_dispatch():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION, OrchestratorAgent

    blocked = OrchestratorAgent._blocked_modalities({
        "scRNA": ["/tmp/rna.h5ad"],
        "scATAC": ["/tmp/fragments.tsv.gz"],
    })

    assert blocked == {}
    # ADR-076: scATAC promoted to production with autonomous dispatch when contracts are green.
    assert MODALITY_VALIDATION["scATAC"]["level"] == "production"
    assert MODALITY_VALIDATION["scATAC"]["dispatch_enabled"] is True
    assert MODALITY_VALIDATION["scATAC"]["promotion_adr"] == "ADR-076"


def test_orchestrator_allows_production_bulk_atac_auto_dispatch():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION, OrchestratorAgent

    blocked = OrchestratorAgent._blocked_modalities({
        "bulk_ATAC": ["/tmp/sample1.bam"],
    })

    assert blocked == {}
    assert MODALITY_VALIDATION["bulk_ATAC"]["level"] == "production"
    assert MODALITY_VALIDATION["bulk_ATAC"]["dispatch_enabled"] is True
    assert MODALITY_VALIDATION["bulk_ATAC"]["promotion_adr"] == "ADR-077"


def test_chromatin_missing_planned_scripts_return_structured_blocker():
    from aria.agents.chromatin_agent import ChromatinAgent

    motif = ChromatinAgent._planned_script_blocker(
        "aria/scripts/chromatin_motifs.py",
        "motif_enrichment",
    )

    assert motif["status"] == "skipped"
    assert motif["reason"] == "script_not_implemented"
    assert motif["validation_level"] == "scaffold"
