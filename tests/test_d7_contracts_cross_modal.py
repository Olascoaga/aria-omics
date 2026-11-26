"""Regression tests for D7 layer-contract fixes in cross-modal concordance & orchestration."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


def test_d7_003_orchestrator_forwards_sample_covariates_in_cross_modal_concordance(tmp_path):
    """D7-003: OrchestratorAgent._run_cross_modal_concordance must forward
    exp_context['sample_covariates'] to cross_modal_concordance.py so Layer D
    runs Freedman-Lane OLS covariate residualization instead of unadjusted."""
    pytest.importorskip("litellm")
    from aria.agents.orchestrator_agent import OrchestratorAgent
    from aria.memory.memory import ARIAMemory
    from aria.scripts.cross_modal_concordance import cross_modal_concordance

    rng = np.random.default_rng(11)
    samples = [f"S{i}" for i in range(12)]
    sex = {s: ("M" if i < 6 else "F") for i, s in enumerate(samples)}
    sex_num = np.array([1.0 if sex[s] == "M" else 0.0 for s in samples])

    rna_rows = {"KDM5D": 100.0 + 400.0 * sex_num + rng.normal(0, 5, 12)}
    atac_rows = {"chrY:1000-1500": 50.0 + 200.0 * sex_num + rng.normal(0, 3, 12)}
    for i in range(10):
        rna_rows[f"AUTO{i}"] = rng.uniform(50, 150, 12)
        atac_rows[f"chr1:{10000 + i * 2000}-{10500 + i * 2000}"] = rng.uniform(20, 80, 12)

    rna_csv = tmp_path / "rna_cpm.csv"
    atac_csv = tmp_path / "atac_cpm.csv"
    pd.DataFrame.from_dict(rna_rows, orient="index", columns=samples).to_csv(rna_csv, index_label="gene")
    pd.DataFrame.from_dict(atac_rows, orient="index", columns=samples).to_csv(atac_csv, index_label="peak")

    gtf = tmp_path / "genes.gtf"
    gtf_lines = [
        'chrY\tTEST\tgene\t1000\t2000\t.\t+\t.\tgene_id "KDM5D"; gene_name "KDM5D";',
    ] + [
        f'chr1\tTEST\tgene\t{10000 + i * 2000}\t{11000 + i * 2000}\t.\t+\t.\tgene_id "AUTO{i}"; gene_name "AUTO{i}";'
        for i in range(10)
    ]
    gtf.write_text("\n".join(gtf_lines) + "\n")

    captured = {}

    class _RunRealEnv:
        def run_in_stack(self, stack, script_path, params):
            captured["params"] = dict(params)
            return cross_modal_concordance(params)

    orch = OrchestratorAgent(ARIAMemory(":memory:"), llm=None)
    orch.env = _RunRealEnv()
    exp_context = {
        "modalities": {"bulk_RNA": [str(rna_csv)], "bulk_ATAC": [str(atac_csv)]},
        "output_dir": str(tmp_path / "out"),
        "rna_cpm_path": str(rna_csv),
        "atac_cpm_path": str(atac_csv),
        "gtf_path": str(gtf),
        "sample_covariates": {"sex": sex},
        "n_permutations": 100,
        "seed": 0,
    }
    agent_results = {
        "bulk_rna_agent": {"status": "done", "findings": {"contrasts": []}},
        "chromatin_agent": {"status": "done", "findings": {"bulk_ATAC": {"findings": {"data_type": "bulk_ATAC"}}}},
    }
    res = orch._run_cross_modal_concordance("exp1", exp_context, agent_results)
    assert "sample_covariates" in captured["params"], (
        "OrchestratorAgent silently dropped exp_context['sample_covariates']"
    )
    assert captured["params"]["sample_covariates"] == {"sex": sex}
    assert res["sample_matched_links"]["covariate_adjustment"]["adjusted"] is True
    assert res["sample_matched_links"]["covariate_adjustment"]["unadjusted_confounding_risk"] is False


def test_d7_003_orchestrator_forwards_every_cross_modal_contract_input_key(tmp_path):
    """Generic contract test over SCRIPT_CONTRACTS['aria/scripts/cross_modal_concordance.py']:
    every input field accepted by the script contract and present in exp_context
    must be forwarded by OrchestratorAgent._run_cross_modal_concordance."""
    pytest.importorskip("litellm")
    from aria.agents.orchestrator_agent import OrchestratorAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    contract = SCRIPT_CONTRACTS["aria/scripts/cross_modal_concordance.py"]
    agent_owned = {"rna_results", "chromatin_results"}

    def _typed_value(field):
        if field.name == "output_dir":
            return str(tmp_path / "study_root")
        if field.name == "padj_threshold":
            return 0.01
        if field.type == "str":
            return f"val_{field.name}"
        if field.type == "path":
            return str(tmp_path / f"{field.name}.dat")
        if field.type == "int":
            return 42
        if field.type == "float":
            return 0.125
        if field.type == "bool":
            return True
        if field.type == "dict":
            return {f"k_{field.name}": {"S1": "A"}}
        if field.type == "list":
            return [f"item_{field.name}"]
        return f"any_{field.name}"

    expected_exp_ctx = {
        field.name: _typed_value(field)
        for field in contract.inputs
        if field.name not in agent_owned
    }
    exp_context = {
        "modalities": {"bulk_RNA": [str(tmp_path / "r.tsv")], "bulk_ATAC": [str(tmp_path / "a.tsv")]},
        **expected_exp_ctx,
    }
    agent_results = {
        "bulk_rna_agent": {"status": "done", "findings": {"contrasts": []}},
        "chromatin_agent": {"status": "done", "findings": {"bulk_ATAC": {"findings": {"data_type": "bulk_ATAC"}}}},
    }

    captured = {}

    class _CapEnv:
        def run_in_stack(self, stack, script_path, params):
            captured["stack"] = stack
            captured["script_path"] = script_path
            captured["params"] = dict(params)
            return {"status": "success", "ran": False, "reason": "captured"}

    orch = OrchestratorAgent(ARIAMemory(":memory:"), llm=None)
    orch.env = _CapEnv()
    orch._run_cross_modal_concordance("exp_contract", exp_context, agent_results)

    forwarded = captured["params"]
    assert captured["script_path"] == "aria/scripts/cross_modal_concordance.py"
    for field in contract.inputs:
        assert field.name in forwarded, (
            f"Contract input {field.name!r} was not forwarded by _run_cross_modal_concordance"
        )
        if field.name in agent_owned:
            assert isinstance(forwarded[field.name], dict)
        elif field.name == "output_dir":
            assert str(forwarded["output_dir"]).startswith(str(expected_exp_ctx["output_dir"]))
        else:
            assert forwarded[field.name] == expected_exp_ctx[field.name], (
                f"Contract input {field.name!r} mismatch: expected {expected_exp_ctx[field.name]!r}, "
                f"got {forwarded[field.name]!r}"
            )
