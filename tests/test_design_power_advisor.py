from __future__ import annotations

from aria.agents.modality_audit import build_capability_matrix
from aria.utils.design_power import assess_design_power


def _ctx(groups: dict, **design_updates) -> dict:
    design = {
        "groups": groups,
        "comparisons": [["treated", "control"]],
        "pseudobulk": {
            "condition_col": "condition",
            "replicate_col": "donor",
            "groupby_col": "cell_type",
            "comparisons": [["treated", "control"]],
        },
    }
    design.update(design_updates)
    return {
        "genome": "hg38",
        "modalities": {"bulk_RNA": ["/data/counts.tsv"]},
        "design": design,
    }


def test_design_power_advisor_marks_well_replicated_design_green():
    result = assess_design_power(
        _ctx({
            "control": ["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"],
            "treated": ["t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8"],
        })
    )

    assert result["status"] == "green"
    contrast = result["contrasts"][0]
    assert contrast["assessment"] == "adequate_for_target_assumption"
    assert 0.0 <= contrast["power_estimate_at_target_log2fc"] <= 1.0
    assert contrast["minimum_detectable_log2fc_at_80_power"] is not None
    assert result["assumptions"]["advisory_only"] is True


def test_design_power_advisor_warns_for_n2_low_power():
    result = assess_design_power(
        _ctx({
            "control": ["c1", "c2"],
            "treated": ["t1", "t2"],
        })
    )

    assert result["status"] == "yellow"
    assert result["contrasts"][0]["assessment"] == "low_power_supported_with_caveat"
    assert any(
        finding["check"] == "design_power_n2_low_power"
        for finding in result["findings"]
    )


def test_design_power_advisor_blocks_n1_inferential_contrast():
    result = assess_design_power(
        _ctx({
            "control": ["c1"],
            "treated": ["t1", "t2"],
        })
    )

    assert result["status"] == "red"
    assert result["contrasts"][0]["assessment"] == "unsupported"
    assert any(
        finding["severity"] == "blocking"
        and finding["check"] == "design_power_unsupported_replicates"
        for finding in result["findings"]
    )


def test_design_power_advisor_blocks_batch_condition_confounding():
    result = assess_design_power(
        _ctx(
            {
                "control": ["c1", "c2", "c3"],
                "treated": ["t1", "t2", "t3"],
            },
            batch_map={
                "c1": "batch_a",
                "c2": "batch_a",
                "c3": "batch_a",
                "t1": "batch_b",
                "t2": "batch_b",
                "t3": "batch_b",
            },
        )
    )

    assert result["status"] == "red"
    assert result["checks"]["batch_condition_confounding"]["confounded"] is True
    assert any(
        finding["check"] == "design_power_batch_condition_confounding"
        for finding in result["findings"]
    )


def test_design_power_advisor_warns_for_imbalanced_contrast():
    result = assess_design_power(
        _ctx({
            "control": ["c1", "c2", "c3", "c4", "c5", "c6"],
            "treated": ["t1", "t2"],
        })
    )

    assert result["status"] == "yellow"
    assert result["contrasts"][0]["balance_ratio"] < 0.5
    assert any(
        finding["check"] == "design_power_imbalanced_contrast"
        for finding in result["findings"]
    )


def test_capability_matrix_surfaces_design_power_preflight():
    matrix = build_capability_matrix(
        _ctx({
            "control": ["c1", "c2"],
            "treated": ["t1", "t2"],
        }),
        modality_validation={"bulk_RNA": {"level": "production",
                                          "dispatch_enabled": True}},
    )

    assert "design_power" in matrix["preflight"]
    assert matrix["preflight"]["design_power"]["status"] == "yellow"
    assert matrix["dispatch"]["allowed"] == ["bulk_RNA"]
    assert matrix["dispatch"]["blocked"] == []
    assert any(
        finding["check"] == "design_power_n2_low_power"
        for finding in matrix["findings"]
    )
