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


def test_capability_matrix_marks_beta_and_scatac_alpha_as_yellow():
    matrix = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {
                "bulk_RNA_raw": ["/data/a.fastq.gz"],
                "scATAC": ["/data/a.h5mu"],
            }
        },
        modality_validation={
            "bulk_RNA_raw": {"level": "beta", "dispatch_enabled": True},
            "scATAC": {
                "level": "alpha",
                "dispatch_enabled": True,
                "reason": "scATAC alpha requires acknowledgement.",
            },
        },
    )

    assert matrix["cards"]["bulk_RNA_raw"]["status"] == "yellow"
    assert matrix["cards"]["bulk_RNA_raw"]["dispatch_policy"] == "requires_ack"
    assert matrix["cards"]["scATAC"]["status"] == "yellow"
    assert matrix["cards"]["scATAC"]["dispatch_policy"] == "requires_ack"
    assert matrix["dispatch"]["requires_ack"] == ["bulk_RNA_raw", "scATAC"]
    assert matrix["dispatch"]["blocked"] == []


def test_bulk_rna_raw_default_registry_promoted_to_production_auto_dispatch():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION

    assert MODALITY_VALIDATION["bulk_RNA_raw"]["level"] == "production"
    matrix = build_capability_matrix(
        {
            "genome": "hg38",
            "modalities": {
                "bulk_RNA_raw": [
                    "/data/ctrl_1_R1.fastq.gz",
                    "/data/ctrl_2_R1.fastq.gz",
                    "/data/ctrl_3_R1.fastq.gz",
                    "/data/treat_1_R1.fastq.gz",
                    "/data/treat_2_R1.fastq.gz",
                    "/data/treat_3_R1.fastq.gz",
                ],
            },
            "design": {
                "groups": {
                    "ctrl": ["ctrl_1", "ctrl_2", "ctrl_3"],
                    "treat": ["treat_1", "treat_2", "treat_3"],
                },
            },
        },
        modality_validation=MODALITY_VALIDATION,
    )
    card = matrix["cards"]["bulk_RNA_raw"]
    assert card["validation_level"] == "production"
    assert card["status"] == "green"
    assert card["dispatch_policy"] == "auto"
    assert matrix["dispatch"]["allowed"] == ["bulk_RNA_raw"]
    assert matrix["dispatch"]["requires_ack"] == []
    assert matrix["dispatch"]["blocked"] == []


def test_scatac_default_registry_promoted_to_production_auto_dispatch():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION

    assert MODALITY_VALIDATION["scATAC"]["level"] == "production"
    assert MODALITY_VALIDATION["scATAC"]["promotion_adr"] == "ADR-076"

    # Exploratory / chromatin-characterization run (no cross-condition contrast)
    matrix = build_capability_matrix(
        {"genome": "hg38", "modalities": {"scATAC": ["/data/a.h5mu"]}},
        modality_validation=MODALITY_VALIDATION,
    )
    card = matrix["cards"]["scATAC"]
    assert card["validation_level"] == "production"
    assert card["status"] == "green"
    assert card["dispatch_policy"] == "auto"
    assert matrix["dispatch"]["allowed"] == ["scATAC"]
    assert matrix["dispatch"]["requires_ack"] == []
    assert matrix["dispatch"]["blocked"] == []
    readiness = card["checks"]["production_readiness"]
    assert readiness["production_promoted"] is True
    assert readiness["tier"] == "production"
    assert readiness["promotion_adr"] == "ADR-076"
    assert readiness["caveated_sub_layers"]["per_cluster_marker_peaks"] == "exploratory"
    assert readiness["caveated_sub_layers"]["differential_tf_footprinting"] == "beta"


def test_scatac_readiness_gates_pseudobulk_replicates_like_scrna():
    from aria.agents.modality_audit import ScATACAuditAgent

    # n>=3 donors per condition -> green / auto
    green_card = ScATACAuditAgent().audit({
        "genome": "hg38",
        "modalities": {"scATAC": ["/data/cohort.h5ad"]},
        "design": {
            "condition_col": "condition",
            "replicate_col": "donor_id",
            "groups": {"ctrl": ["d1", "d2", "d3"], "treat": ["d4", "d5", "d6"]},
            "comparisons": [{"test": "treat", "reference": "ctrl"}],
        },
    })
    assert green_card["status"] == "green"
    assert green_card["dispatch_policy"] == "auto"

    # n==2 donors -> yellow / requires_ack
    yellow_card = ScATACAuditAgent().audit({
        "genome": "hg38",
        "modalities": {"scATAC": ["/data/cohort.h5ad"]},
        "design": {
            "condition_col": "condition",
            "replicate_col": "donor_id",
            "groups": {"ctrl": ["d1", "d2"], "treat": ["d3", "d4"]},
            "comparisons": [{"test": "treat", "reference": "ctrl"}],
        },
    })
    assert yellow_card["status"] == "yellow"
    assert yellow_card["dispatch_policy"] == "requires_ack"
    assert any(f["check"] == "scatac_pseudobulk_low_power" for f in yellow_card["findings"])

    # n<2 donors on contrasted run -> red / blocked
    red_card = ScATACAuditAgent().audit({
        "genome": "hg38",
        "modalities": {"scATAC": ["/data/cohort.h5ad"]},
        "design": {
            "condition_col": "condition",
            "replicate_col": "donor_id",
            "groups": {"ctrl": ["d1"], "treat": ["d2", "d3"]},
            "comparisons": [{"test": "treat", "reference": "ctrl"}],
        },
    })
    assert red_card["status"] == "red"
    assert red_card["dispatch_policy"] == "blocked"
    assert any(f["check"] == "scatac_pseudobulk_under_replicated" for f in red_card["findings"])


def test_scatac_explicit_beta_override_still_requires_ack():
    # Explicit downgrade override in modality_validation still enforces beta + ack.
    matrix = build_capability_matrix(
        {"genome": "hg38", "modalities": {"scATAC": ["/data/a.h5mu"]}},
        modality_validation={"scATAC": {"level": "beta", "dispatch_enabled": True}},
    )
    card = matrix["cards"]["scATAC"]
    assert card["validation_level"] == "beta"
    assert card["status"] == "yellow"
    assert card["dispatch_policy"] == "requires_ack"
    codes = {f["check"] for f in card["findings"]}
    assert "chromatin_readiness_beta_ack_required" in codes
    assert card["checks"]["production_readiness"]["production_promoted"] is False


def test_capability_matrix_marks_bulk_atac_beta_requires_ack():
    matrix = build_capability_matrix(
        {"genome": "hg38", "modalities": {"bulk_ATAC": ["/data/fragments.tsv.gz"]}},
        modality_validation={
            "bulk_ATAC": {
                "level": "beta",
                "dispatch_enabled": True,
                "reason": "bulk ATAC QC/peak calling beta.",
            },
        },
    )

    card = matrix["cards"]["bulk_ATAC"]
    assert card["validation_level"] == "beta"
    assert card["status"] == "yellow"
    assert card["dispatch_policy"] == "requires_ack"
    assert matrix["dispatch"]["requires_ack"] == ["bulk_ATAC"]
    assert matrix["dispatch"]["blocked"] == []
    codes = {f["check"] for f in card["findings"]}
    assert "bulk_atac_beta_ack_required" in codes
    readiness = card["checks"]["production_readiness"]
    assert readiness["technical_blockers_cleared"] is True
    assert readiness["open_blockers"] == []
    assert readiness["production_promoted"] is False
    assert readiness["tier"] == "beta"


def test_bulk_atac_default_registry_stays_beta_until_tier_promotion():
    matrix = build_capability_matrix(
        {"genome": "hg38", "modalities": {"bulk_ATAC": ["/data/fragments.tsv.gz"]}},
        modality_validation={"bulk_ATAC": {"level": "beta",
                                           "dispatch_enabled": True}},
    )

    card = matrix["cards"]["bulk_ATAC"]
    assert card["validation_level"] == "beta"
    assert card["status"] == "yellow"
    assert card["dispatch_policy"] == "requires_ack"
    assert card["checks"]["production_readiness"]["production_promoted"] is False
    assert card["checks"]["production_readiness"]["technical_blockers_cleared"] is True


def test_adr056_bulk_atac_stays_beta_while_adr076_promotes_scatac_explicitly():
    # ADR-056: bulk_ATAC stays beta + requires_ack with governance blockers listed;
    # ADR-076 explicitly promotes scATAC to production while preserving ADR-056 policy.
    from aria.agents.modality_audit import (
        _ATAC_GOVERNANCE_BLOCKERS,
        _atac_production_readiness,
    )

    bulk = _atac_production_readiness("bulk_ATAC")
    assert bulk["production_promoted"] is False
    assert bulk["tier"] == "beta" and bulk["requires_ack"] is True
    assert bulk["promotion_policy_adr"] == "ADR-056"
    assert tuple(bulk["governance_blockers"]) == _ATAC_GOVERNANCE_BLOCKERS
    assert "explicit_production_promotion_adr" in bulk["governance_blockers"]

    sc = _atac_production_readiness("scATAC")
    assert sc["production_promoted"] is True
    assert sc["tier"] == "production" and sc["requires_ack"] is False
    assert sc["promotion_policy_adr"] == "ADR-056"
    assert sc["promotion_adr"] == "ADR-076"
    # No legacy ambiguous field that could be misread as "is production".
    assert "production_ready" not in sc and "production_ready" not in bulk


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
        "genome": "hg38",
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


def test_orchestrator_does_not_filter_requires_ack_scatac_before_dispatch():
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
                "allowed": ["bulk_RNA", "scATAC"],
                "requires_ack": ["scATAC"],
                "blocked": [],
            }
        }
    }

    filtered = orch._apply_capability_dispatch_policy(exp_context, audit_result)

    assert filtered["modalities"] == exp_context["modalities"]
    assert "blocked_modalities_by_capability" not in filtered
