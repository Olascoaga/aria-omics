import pandas as pd


def test_design_matrix_blocks_complete_batch_condition_confounding():
    from aria.utils.design_matrix import validate_design_matrix

    meta = pd.DataFrame({
        "condition": ["ctrl", "ctrl", "stim", "stim"],
        "batch": ["b1", "b1", "b2", "b2"],
    }, index=["c1", "c2", "s1", "s2"])

    result = validate_design_matrix(
        meta,
        condition_col="condition",
        covariates=["batch"],
        min_replicates_per_condition=2,
    )

    assert result["status"] == "blocking"
    checks = {issue["check"] for issue in result["issues"]}
    assert "condition_covariate_confounding" in checks
    assert "rank_deficient_design" in checks


def test_design_matrix_detects_continuous_covariate():
    from aria.utils.design_matrix import validate_design_matrix

    meta = pd.DataFrame({
        "condition": ["ctrl", "ctrl", "ctrl", "stim", "stim", "stim"],
        "composition": [0.10, 0.15, 0.20, 0.50, 0.55, 0.60],
    })

    result = validate_design_matrix(
        meta,
        condition_col="condition",
        covariates=["composition"],
        min_replicates_per_condition=3,
    )

    assert result["status"] == "clean"
    assert result["continuous_factors"] == ["composition"]
    assert result["categorical_factors"] == ["condition"]


def test_design_matrix_warns_on_n1_design_cells():
    from aria.utils.design_matrix import validate_design_matrix

    meta = pd.DataFrame({
        "condition": ["ctrl", "ctrl", "ctrl", "stim", "stim", "stim"],
        "batch": ["b1", "b1", "b2", "b1", "b2", "b2"],
    })

    result = validate_design_matrix(
        meta,
        condition_col="condition",
        covariates=["batch"],
        min_replicates_per_condition=3,
    )

    assert result["status"] == "warnings"
    assert any(issue["check"] == "n1_design_cells" for issue in result["issues"])


def test_audit_agent_surfaces_design_matrix_confounding(tmp_path):
    from aria.agents.audit_agent import AuditAgent

    counts = pd.DataFrame({
        "ctrl_1": [10, 20, 30],
        "ctrl_2": [11, 21, 31],
        "stim_1": [50, 60, 70],
        "stim_2": [51, 61, 71],
    }, index=["gene1", "gene2", "gene3"])
    counts_path = tmp_path / "counts.tsv"
    counts.to_csv(counts_path, sep="\t")

    design = {
        "main_factor": "condition",
        "groups": {"ctrl": ["ctrl_1", "ctrl_2"], "stim": ["stim_1", "stim_2"]},
        "batch_factor": "batch",
        "batch_map": {
            "ctrl_1": "b1",
            "ctrl_2": "b1",
            "stim_1": "b2",
            "stim_2": "b2",
        },
    }

    agent = AuditAgent.__new__(AuditAgent)
    findings = agent._check_design_matrix_sanity([str(counts_path)], design)

    assert any(
        finding["severity"] == "blocking"
        and "confounded" in finding["message"]
        for finding in findings
    )
