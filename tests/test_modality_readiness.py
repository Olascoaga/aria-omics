from __future__ import annotations

from aria.agents.modality_audit import (
    ScRNAAuditAgent,
    build_capability_matrix,
)


def _scrna_context(groups: dict, *, modality: str = "scRNA") -> dict:
    return {
        "modalities": {modality: ["/data/input.h5ad"]},
        "design": {
            "groups": groups,
            "pseudobulk": {
                "condition_col": "condition",
                "replicate_col": "donor",
                "groupby_col": "cell_type",
            },
        },
    }


def test_scrna_readiness_green_for_publication_ready_pseudobulk():
    card = ScRNAAuditAgent().audit(
        _scrna_context({
            "ctrl": ["c1", "c2", "c3"],
            "stim": ["s1", "s2", "s3"],
        })
    )

    assert card["status"] == "green"
    assert card["dispatch_policy"] == "auto"
    assert card["checks"]["pseudobulk_replicates"]["min_replicates"] == 3


def test_scrna_readiness_yellow_requires_ack_for_n2_pseudobulk():
    card = ScRNAAuditAgent().audit(
        _scrna_context({
            "ctrl": ["c1", "c2"],
            "stim": ["s1", "s2"],
        })
    )

    assert card["status"] == "yellow"
    assert card["dispatch_policy"] == "requires_ack"
    assert any(f["severity"] == "warning" for f in card["findings"])


def test_scrna_readiness_red_blocks_under_replicated_pseudobulk():
    card = ScRNAAuditAgent().audit(
        _scrna_context({
            "ctrl": ["c1"],
            "stim": ["s1", "s2"],
        })
    )

    assert card["status"] == "red"
    assert card["dispatch_policy"] == "blocked"
    assert any(f["severity"] == "blocking" for f in card["findings"])


def test_capability_matrix_marks_beta_as_yellow_and_scaffold_as_red():
    matrix = build_capability_matrix(
        {
            "modalities": {
                "bulk_RNA_raw": ["/data/a.fastq.gz"],
                "scATAC": ["/data/a.h5mu"],
            }
        },
        modality_validation={
            "bulk_RNA_raw": {"level": "beta", "dispatch_enabled": True},
            "scATAC": {
                "level": "scaffold",
                "dispatch_enabled": False,
                "reason": "scATAC validation is not closed.",
            },
        },
    )

    assert matrix["cards"]["bulk_RNA_raw"]["status"] == "yellow"
    assert matrix["cards"]["bulk_RNA_raw"]["dispatch_policy"] == "requires_ack"
    assert matrix["cards"]["scATAC"]["status"] == "red"
    assert matrix["cards"]["scATAC"]["dispatch_policy"] == "blocked"
    assert matrix["dispatch"]["requires_ack"] == ["bulk_RNA_raw"]
    assert matrix["dispatch"]["blocked"] == ["scATAC"]


def test_audit_agent_surfaces_capability_matrix_without_heavy_checks(monkeypatch):
    from aria.agents.audit_agent import AuditAgent

    agent = AuditAgent.__new__(AuditAgent)
    agent.publish_status = lambda *args, **kwargs: None
    monkeypatch.setattr(agent, "_check_replicate_correlation", lambda *args: [])
    monkeypatch.setattr(agent, "_check_pca_batch_dominance", lambda *args: [])
    monkeypatch.setattr(agent, "_check_design_matrix_sanity", lambda *args: [])
    monkeypatch.setattr(agent, "_check_star_alignment", lambda *args: [])

    result = agent.run_audit(
        _scrna_context({
            "ctrl": ["c1", "c2"],
            "stim": ["s1", "s2"],
        }),
        "exp-readiness",
    )

    assert result["status"] == "blocking"
    assert result["capability_matrix"]["cards"]["scRNA"]["status"] == "yellow"
    assert result["capability_matrix"]["dispatch"]["requires_ack"] == ["scRNA"]
    assert any(f["check"] == "scrna_pseudobulk_low_power"
               for f in result["findings"])


def test_orchestrator_filters_red_modalities_before_dispatch():
    from aria.agents.orchestrator_agent import OrchestratorAgent

    orch = OrchestratorAgent.__new__(OrchestratorAgent)
    exp_context = {
        "modalities": {
            "bulk_RNA": ["/data/counts.tsv"],
            "scATAC": ["/data/atac.h5mu"],
        }
    }
    audit_result = {
        "capability_matrix": {
            "dispatch": {
                "allowed": ["bulk_RNA"],
                "requires_ack": [],
                "blocked": ["scATAC"],
            }
        }
    }

    filtered = orch._apply_capability_dispatch_policy(exp_context, audit_result)

    assert filtered["modalities"] == {"bulk_RNA": ["/data/counts.tsv"]}
    assert filtered["blocked_modalities_by_capability"] == ["scATAC"]
