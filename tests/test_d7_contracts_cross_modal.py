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


def test_d7_014_cross_modal_adapter_low_replication_on_low_n_cohort(tmp_path):
    """D7-014: sample_matched_layer must emit low_power_warning for low-N cohorts
    (n < 10) and cross_modal_evidence must inherit LOW_REPLICATION rather than
    silently masking an unwritten key."""
    from aria.agents.narrative.hypothesis.adapters.cross_modal import cross_modal_evidence
    from aria.agents.narrative.hypothesis.caveats import LOW_REPLICATION
    from aria.scripts.cross_modal_concordance import cross_modal_concordance

    rng = np.random.default_rng(31)
    samples = [f"S{i}" for i in range(6)]
    x = np.linspace(1.0, 6.0, 6)
    rna_rows = {"G0": 50.0 + 100.0 * x + rng.normal(0, 0.5, 6)}
    atac_rows = {"chr1:1000-1500": 20.0 + 50.0 * x + rng.normal(0, 0.5, 6)}
    for i in range(1, 10):
        rna_rows[f"G{i}"] = rng.uniform(20, 80, 6)
        atac_rows[f"chr1:{10000 + i * 2000}-{10500 + i * 2000}"] = rng.uniform(10, 40, 6)

    rna_csv = tmp_path / "rna.csv"
    atac_csv = tmp_path / "atac.csv"
    pd.DataFrame.from_dict(rna_rows, orient="index", columns=samples).to_csv(rna_csv, index_label="gene")
    pd.DataFrame.from_dict(atac_rows, orient="index", columns=samples).to_csv(atac_csv, index_label="peak")
    gtf = tmp_path / "g.gtf"
    gtf.write_text(
        'chr1\tTEST\tgene\t1000\t2000\t.\t+\t.\tgene_id "G0"; gene_name "G0";\n'
        + "\n".join(
            f'chr1\tTEST\tgene\t{10000 + i * 2000}\t{11000 + i * 2000}\t.\t+\t.\tgene_id "G{i}"; gene_name "G{i}";'
            for i in range(1, 10)
        )
        + "\n"
    )
    cm = cross_modal_concordance({
        "rna_results": {"findings": {"contrasts": []}},
        "chromatin_results": {"findings": {}},
        "rna_cpm_path": str(rna_csv),
        "atac_cpm_path": str(atac_csv),
        "gtf_path": str(gtf),
        "output_dir": str(tmp_path / "out"),
        "n_permutations": 100,
        "seed": 0,
    })
    links = cm["sample_matched_links"]
    assert links["ran"] is True
    assert links["n_shared_samples"] == 6
    assert links["low_power_warning"] is True
    assert links["top_links"]

    sigs = cross_modal_evidence({"cross_modal_concordance": cm}, run_ledger=None, exp_ctx={})
    assert sigs
    assert LOW_REPLICATION in sigs[0].caveats_inherited

    # Also verify adapter handles absent low_power_warning cleanly using n_shared_samples:
    links_no_flag = {k: v for k, v in links.items() if k != "low_power_warning"}
    sigs_absent = cross_modal_evidence(
        {"cross_modal_concordance": {**cm, "sample_matched_links": links_no_flag}},
        run_ledger=None,
        exp_ctx={},
    )
    assert sigs_absent and LOW_REPLICATION in sigs_absent[0].caveats_inherited


def test_d7_015_cross_modal_adapter_matches_pearson_direction_tested_by_padj(tmp_path):
    """D7-015: sample_matched_layer evaluates empirical_pvalue/padj and link polarity
    from pearson_r. cross_modal_evidence must use pearson_r rather than inverting
    direction via spearman_r on high-leverage links."""
    from aria.agents.narrative.hypothesis.adapters.cross_modal import cross_modal_evidence
    from aria.scripts.cross_modal_concordance import cross_modal_concordance

    samples = [f"S{i}" for i in range(12)]
    rna_sig = np.array([10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 1e6, 1.1e6])
    atac_sig = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 1e6, 1.1e6])
    rng = np.random.default_rng(7)
    rna_rows = {"LEVERAGE1": rna_sig}
    atac_rows = {"chr1:1000-1500": atac_sig}
    for i in range(10):
        rna_rows[f"BG{i}"] = rng.uniform(50, 150, 12)
        atac_rows[f"chr1:{20000 + i * 2000}-{20500 + i * 2000}"] = rng.uniform(20, 80, 12)

    rna_csv = tmp_path / "rna.csv"
    atac_csv = tmp_path / "atac.csv"
    pd.DataFrame.from_dict(rna_rows, orient="index", columns=samples).to_csv(rna_csv, index_label="gene")
    pd.DataFrame.from_dict(atac_rows, orient="index", columns=samples).to_csv(atac_csv, index_label="peak")
    gtf = tmp_path / "g.gtf"
    gtf.write_text(
        'chr1\tTEST\tgene\t1000\t2000\t.\t+\t.\tgene_id "LEVERAGE1"; gene_name "LEVERAGE1";\n'
        + "\n".join(
            f'chr1\tTEST\tgene\t{20000 + i * 2000}\t{21000 + i * 2000}\t.\t+\t.\tgene_id "BG{i}"; gene_name "BG{i}";'
            for i in range(10)
        )
        + "\n"
    )
    cm = cross_modal_concordance({
        "rna_results": {"findings": {"contrasts": []}},
        "chromatin_results": {"findings": {}},
        "rna_cpm_path": str(rna_csv),
        "atac_cpm_path": str(atac_csv),
        "gtf_path": str(gtf),
        "output_dir": str(tmp_path / "out"),
        "n_permutations": 200,
        "seed": 0,
    })
    links = cm["sample_matched_links"]
    assert links["ran"] is True
    assert links["n_positive_links"] >= 1 and links["n_negative_links"] == 0
    row = next(r for r in links["top_links"] if r["gene"] == "LEVERAGE1")
    assert row["pearson_r"] > 0.9 and row["spearman_r"] < 0.0 and row["padj"] < 0.05

    sigs = cross_modal_evidence({"cross_modal_concordance": cm}, run_ledger=None, exp_ctx={})
    sig = next(s for s in sigs if s.entity == "LEVERAGE1")
    assert sig.direction == "up"
    assert sig.value == pytest.approx(round(row["pearson_r"], 4))
