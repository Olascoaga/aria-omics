"""Verification suite for ARIA v4.9 Milestone 4:
1. P1 "In-Silico Reviewer #2" Claim-Specific Falsification Certificates
   (`aria/agents/narrative/falsification_certificates.py`).
2. Tier 2.1 Typed Derived Query AST & Zero-Leakage Benchmark Adapter
   (`aria/benchmarks/derived_query_ast.py`, `aria/benchmarks/adapter.py`).
3. 600-Case Deterministic Claim-Boundary & Metamorphic Report Falsification Suite
   and ADR-083 `automated_v1_gate` (`aria/benchmarks/metamorphic_report_falsification.py`,
   `aria/benchmarks/preprint_freeze.py`).
"""

from __future__ import annotations

from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
import pytest

from aria.agents.narrative.falsification_certificates import (
    build_falsification_certificates,
    compute_lobo_jackknife,
)
from aria.benchmarks.adapter import (
    BenchmarkPrediction,
    BenchmarkTaskSpec,
    evaluate_benchmark_predictions,
    isolated_benchmark_environment,
    scan_codebase_for_benchmark_leakage,
)
from aria.benchmarks.derived_query_ast import (
    DerivedQueryPlan,
    DerivedQueryViolation,
    FilterOp,
    GroupSummaryOp,
    ProjectOp,
    RankOp,
    TableContract,
    ValidatedJoinOp,
    execute_derived_query,
)
from aria.benchmarks.metamorphic_report_falsification import (
    run_metamorphic_falsification_suite,
)
from aria.benchmarks.preprint_freeze import (
    DEFERRED_V2_HUMAN_GOLD_LANES,
    LANES,
    build_inventory,
)


def test_lobo_jackknife_enforces_biological_units_and_identifies_outlier_driven_claims():
    # 1. n=2v2 is refused as not_testable_insufficient_replicates (never fabricated)
    res_n2 = compute_lobo_jackknife(
        ["G1"],
        np.array([[5.0, 5.2, 1.0, 1.1]]),
        ["KO1", "KO2", "WT1", "WT2"],
        {"KO1": "KO", "KO2": "KO", "WT1": "WT", "WT2": "WT"},
        numerator="KO",
        denominator="WT",
    )
    assert res_n2["status"] == "not_testable_insufficient_replicates"
    assert res_n2["minimum_required_per_arm"] == 3

    # 2. Duplicate unit IDs (unaggregated technical replicates) fail closed
    with pytest.raises(ValueError, match="unique biological units"):
        compute_lobo_jackknife(
            ["G1"],
            np.array([[5.0, 5.1, 5.2, 1.0, 1.1, 1.2]]),
            ["KO1", "KO1", "KO2", "WT1", "WT2", "WT3"],
            {"KO1": "KO", "KO2": "KO", "WT1": "WT", "WT2": "WT", "WT3": "WT"},
            numerator="KO",
            denominator="WT",
        )

    # 3. 4v4 cohort: G_STABLE is consistent across all 4 donors; G_OUTLIER is driven solely by KO_4
    units = ["KO_1", "KO_2", "KO_3", "KO_4", "WT_1", "WT_2", "WT_3", "WT_4"]
    groups = {u: ("KO" if u.startswith("KO") else "WT") for u in units}
    mat = np.array([
        # G_STABLE: uniformly high in KO (~6.0) vs WT (~1.0)
        [6.1, 5.9, 6.2, 6.0, 1.1, 0.9, 1.0, 1.2],
        # G_OUTLIER: baseline ~1.0 in KO_1..3, huge single-donor spike in KO_4 (14.0)
        [1.0, 1.1, 0.9, 14.0, 1.0, 1.0, 1.1, 0.9],
        # G_NULL: noise
        [2.0, 2.1, 1.9, 2.0, 2.0, 1.9, 2.1, 2.0],
    ])
    lobo = compute_lobo_jackknife(
        ["G_STABLE", "G_OUTLIER", "G_NULL"],
        mat,
        units,
        groups,
        numerator="KO",
        denominator="WT",
    )
    assert lobo["status"] == "measured"
    assert lobo["n_jackknife_folds"] == 8
    assert lobo["per_feature"]["G_STABLE"]["single_unit_driven"] is False
    assert lobo["per_feature"]["G_STABLE"]["sign_consistent_fraction"] == 1.0
    assert lobo["per_feature"]["G_STABLE"]["max_jackknife_padj"] < 0.05

    assert lobo["per_feature"]["G_OUTLIER"]["single_unit_driven"] is True
    assert lobo["per_feature"]["G_OUTLIER"]["most_influential_unit"] == "KO_4"

    # 4. Orthogonal certificate packaging enforces anti-p-hacking & non-causal walls
    cert_bundle = build_falsification_certificates(
        claims=[{"claim_id": "bulk.de.1", "analysis": "differential_expression", "tier": "associative", "modality": "bulk_rna"}],
        agent_results={},
        exp_ctx={"lobo_jackknife": lobo, "design": {"replicates": {"KO": 4, "WT": 4}}},
    )
    assert cert_bundle["feeds_into_parameter_selection"] is False
    assert cert_bundle["scalar_fragility_score_emitted"] is False
    assert cert_bundle["causal_promotion_from_robustness"] is False
    assert cert_bundle["certificates"][0]["lobo_biological_unit_jackknife"]["status"] == "measured"


def test_derived_query_ast_enforces_unit_and_cardinality_guards():
    de_contract = TableContract(
        table_id="pseudobulk_de",
        biological_unit="gene",
        primary_key=("gene",),
        ledger_node_id="ledger://scrna/pseudobulk_de",
    )
    de_df = pd.DataFrame([
        {"gene": "GFAP", "pathway": "gliosis", "log2fc": 2.4, "padj": 0.001},
        {"gene": "AQP4", "pathway": "gliosis", "log2fc": 1.8, "padj": 0.004},
        {"gene": "ACTB", "pathway": "housekeeping", "log2fc": 0.05, "padj": 0.92},
    ])
    annot_contract = TableContract(
        table_id="gene_biotype",
        biological_unit="gene",
        primary_key=("gene",),
        ledger_node_id="ledger://scrna/annotation",
    )
    annot_df = pd.DataFrame([
        {"gene": "GFAP", "biotype": "protein_coding"},
        {"gene": "AQP4", "biotype": "protein_coding"},
        {"gene": "ACTB", "biotype": "protein_coding"},
    ])
    cell_contract = TableContract(
        table_id="cell_obs",
        biological_unit="cell",
        primary_key=("cell_id",),
        ledger_node_id="ledger://scrna/qc",
    )
    cell_df = pd.DataFrame([
        {"cell_id": "c1", "gene": "GFAP"},
        {"cell_id": "c2", "gene": "GFAP"},
    ])

    tables = {
        "pseudobulk_de": (de_contract, de_df),
        "gene_biotype": (annot_contract, annot_df),
        "cell_obs": (cell_contract, cell_df),
    }

    # Valid pipeline: filter -> 1:1 join -> rank -> project
    plan = DerivedQueryPlan(
        base_table_id="pseudobulk_de",
        steps=(
            FilterOp(column="padj", op="<", value=0.05),
            ValidatedJoinOp(right_table_id="gene_biotype", on=("gene",), expected_cardinality="1:1"),
            RankOp(by=("log2fc",), ascending=(False,), limit=2),
            ProjectOp(columns=("gene", "log2fc", "biotype")),
        ),
    )
    receipt = execute_derived_query(plan, tables)
    assert receipt["input_denominator"] == 3
    assert receipt["output_row_count"] == 2
    assert [r["gene"] for r in receipt["rows"]] == ["GFAP", "AQP4"]
    assert receipt["ledger_node_ids"] == ["ledger://scrna/pseudobulk_de", "ledger://scrna/annotation"]

    # GroupSummaryOp explicitly stamps grouped biological_unit
    grp_plan = DerivedQueryPlan(
        base_table_id="pseudobulk_de",
        steps=(
            FilterOp(column="padj", op="<", value=0.05),
            GroupSummaryOp(by=("pathway",), aggregations={"log2fc": "mean", "gene": "count"}),
        ),
    )
    grp_res = execute_derived_query(grp_plan, tables)
    assert grp_res["biological_unit"] == "grouped:pathway"
    assert grp_res["rows"][0]["gene"] == 2

    # Cross-unit join (gene x cell) WITHOUT allow_unit_projection MUST fail closed!
    bad_unit_plan = DerivedQueryPlan(
        base_table_id="pseudobulk_de",
        steps=(
            ValidatedJoinOp(right_table_id="cell_obs", on=("gene",), expected_cardinality="1:m"),
        ),
    )
    with pytest.raises(DerivedQueryViolation, match="biological-unit"):
        execute_derived_query(bad_unit_plan, tables)

    # 1:1 cardinality violation when right table has duplicates on join key MUST fail closed!
    bad_card_plan = DerivedQueryPlan(
        base_table_id="pseudobulk_de",
        steps=(
            ValidatedJoinOp(
                right_table_id="cell_obs",
                on=("gene",),
                expected_cardinality="1:1",
                allow_unit_projection=True,
            ),
        ),
    )
    with pytest.raises(DerivedQueryViolation, match="violated 1:1 cardinality"):
        execute_derived_query(bad_card_plan, tables)


def test_zero_leakage_benchmark_adapter_and_static_codebase_guard(tmp_path):
    env = isolated_benchmark_environment(tmp_path / "bench_cell")
    assert env["ARIA_LLM_CACHE"] == "0"
    assert env["ARIA_REPRODUCIBLE"] == "1"
    assert env["ARIA_DISABLE_HISTORICAL_BONUS"] == "1"

    tasks = [
        BenchmarkTaskSpec("scBench", "t1_supported_ok", supported_by_native_contract=True),
        BenchmarkTaskSpec("scBench", "t2_supported_ok", supported_by_native_contract=True),
        BenchmarkTaskSpec("scBench", "t3_tapestri_dna_protein", supported_by_native_contract=False, should_abstain=True),
        BenchmarkTaskSpec("scBench", "t4_false_mechanism_probe", supported_by_native_contract=True, should_abstain=True),
    ]
    preds = [
        BenchmarkPrediction("t1_supported_ok", abstained=False, is_correct=True),
        BenchmarkPrediction("t2_supported_ok", abstained=False, is_correct=True),
        BenchmarkPrediction("t3_tapestri_dna_protein", abstained=True, is_correct=True, abstain_reason="unsupported_modality"),
        BenchmarkPrediction("t4_false_mechanism_probe", abstained=True, is_correct=True, abstain_reason="causal_not_licensed"),
    ]
    summary = evaluate_benchmark_predictions("scBench_public", tasks, preds)
    assert summary["total_suite_size"] == 4
    assert summary["native_contract_coverage"] == 0.5
    assert summary["conditional_accuracy_on_supported"] == 1.0
    assert summary["unconditional_end_to_end_accuracy"] == 1.0
    assert summary["abstention_precision"] == 1.0
    assert summary["abstention_recall"] == 1.0
    assert summary["silent_error_rate"] == 0.0

    # Static CI guard: zero external benchmark task IDs hardcoded inside `aria/`
    repo_root = Path(__file__).resolve().parents[1]
    leaks = scan_codebase_for_benchmark_leakage(
        repo_root / "aria",
        forbidden_tokens=(
            "DE01_pseudobulk_de",
            "NRM01_sparse_normalization",
            "dr_05_pca_preprocessing_sentinels",
            "gse297023_q2_senescence_not_mechanism",
            "tapestri_ccus",
            "compbiobench-data-v1",
        ),
    )
    assert leaks == []


def test_600_case_metamorphic_falsification_suite_and_adr083_gate(tmp_path):
    suite = run_metamorphic_falsification_suite(n_seeds=10)
    assert suite["total_cases"] == 600
    assert suite["n_invalid_mutations"] == 300
    assert suite["n_invalid_caught"] == 300
    assert suite["false_publication_rate"] == 0.0
    assert suite["n_valid_controls"] == 300
    assert suite["n_valid_published"] == 300
    assert suite["valid_retention_rate"] == 1.0
    assert suite["all_gates_pass"] is True
    assert suite["failures"] == []

    # Verify ADR-083 Preprint v1 automated gate preserves the 3 human_gold lanes untouched
    repo_root = Path(__file__).resolve().parents[1]
    lane_ids = {l["lane_id"] for l in LANES}
    assert "c4_metamorphic_600_falsification" in lane_ids
    assert DEFERRED_V2_HUMAN_GOLD_LANES.issubset(lane_ids)

    inv = build_inventory(repo_root, tmp_path / "preprint_out")
    auto_gate = inv["automated_v1_gate"]
    assert auto_gate["policy_adr"] == "ADR-083"
    assert auto_gate["n_required_machine_lanes"] == len(LANES) - 3
    deferred_ids = {d["lane_id"] for d in auto_gate["deferred_v2_human_lanes"]}
    assert deferred_ids == set(DEFERRED_V2_HUMAN_GOLD_LANES)
    assert all(
        d["disposition"] == "deferred_to_v2_reader_study"
        for d in auto_gate["deferred_v2_human_lanes"]
    )

    # Human-gold CSVs remain 100% untouched (0 diff)
    proc = subprocess.run(
        ["git", "diff", "--name-only", "--", "docs/benchmark_results/preprint_v1/human/"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=True,
    )
    assert proc.stdout.strip() == ""
