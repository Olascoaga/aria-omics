def test_orchestrator_blocks_scaffold_chromatin_dispatch():
    from aria.agents.orchestrator_agent import OrchestratorAgent

    blocked = OrchestratorAgent._blocked_modalities({
        "scRNA": ["/tmp/rna.h5ad"],
        "scATAC": ["/tmp/fragments.tsv.gz"],
    })

    assert set(blocked) == {"scATAC"}
    assert blocked["scATAC"]["level"] == "scaffold"


def test_chromatin_missing_planned_scripts_return_structured_blocker():
    from aria.agents.chromatin_agent import ChromatinAgent

    motif = ChromatinAgent._planned_script_blocker(
        "aria/scripts/chromatin_motifs.py",
        "motif_enrichment",
    )

    assert motif["status"] == "skipped"
    assert motif["reason"] == "script_not_implemented"
    assert motif["validation_level"] == "scaffold"
