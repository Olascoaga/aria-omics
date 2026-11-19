"""Verification guards for external benchmark Track 3 (scBench adapter + CompBioBench v1)."""

from __future__ import annotations

from pathlib import Path

from aria.benchmarks.adapter import scan_codebase_for_benchmark_leakage
from scripts.external_benchmarks.live_codeact_agent import (
    assert_blind_workspace,
    compute_aria_package_sha256,
    extract_json_answer,
)
from scripts.external_benchmarks.run_compbiobench_omics import (
    load_compbiobench_tasks,
    run_compbiobench_omics,
)
from scripts.external_benchmarks.scbench_aria_agent import run_synthetic_smoke

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_blind_workspace_guard_blocks_solution_notebooks(tmp_path):
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "counts.csv").write_text("g,s1\nG1,10\n", encoding="utf-8")
    assert_blind_workspace(tmp_path)

    nb = tmp_path / "CapsuleNotebook-123" / "solution_executed.ipynb"
    nb.parent.mkdir()
    nb.write_text("{}", encoding="utf-8")
    with pytest.raises(RuntimeError, match="ANTI-CONTAMINATION"):
        assert_blind_workspace(tmp_path)


def test_extract_json_answer_parses_final_tags_and_fences():
    assert extract_json_answer('<FINAL_ANSWER>{"a": 1, "b": true}</FINAL_ANSWER>') == {"a": 1, "b": True}
    assert extract_json_answer('```json\n{"choice": "E", "abstained": true}\n```') == {
        "choice": "E",
        "abstained": True,
    }


def test_scbench_adapter_synthetic_smoke_runs_real_pseudobulk_and_abstains_non_rna():
    res = run_synthetic_smoke()
    assert res["status"] == "pass"
    assert res["pseudobulk_n_sig"] >= 6
    assert res["lobo_status"] == "measured"
    assert res["non_rna_abstain_reason"].startswith("unsupported_")
    assert len(res["aria_package_sha256"]) == 64


def test_compbiobench_loader_and_dry_run(tmp_path):
    tsv = tmp_path / "compbiobench.v1.tsv"
    tsv.write_text(
        "question_id\tcurator_name\tdomain\tquestion_style\tskills_tested\tquestion\tinternet_required\tgpu_preferred\tfile_paths\n"
        "q_sc_ok\tcur\tSingle-cell\tRoutine Analysis\tscanpy\tCompute mean\tFalse\tFalse\tpbmc3k.h5ad\n"
        "q_bam_oos\tcur\tSingle-cell\tMetadata Recovery\tsamtools\tDemux\tFalse\tFalse\tpooled.bam\n"
        "q_pop_oos\tcur\tPopulation Genetics\tRoutine Analysis\tplink\tGWAS\tFalse\tFalse\tgwas.bed\n",
        encoding="utf-8",
    )
    tasks = load_compbiobench_tasks(tsv)
    assert len(tasks) == 3
    by_id = {t["question_id"]: t for t in tasks}
    assert by_id["q_sc_ok"]["native_contract_eligible"] is True
    assert by_id["q_bam_oos"]["native_contract_eligible"] is False
    assert by_id["q_pop_oos"]["offline_omics_eligible"] is False

    out_json = tmp_path / "out.json"
    out_tsv = tmp_path / "out.tsv"
    summary = run_compbiobench_omics(
        data_root=tmp_path,
        output_json=out_json,
        submission_tsv=out_tsv,
        dry_run=True,
    )
    assert summary["full_compbiobench_v1_denominator"] == 3
    assert summary["n_offline_omics_tasks"] == 2
    assert summary["n_native_contract_eligible"] == 1
    assert out_json.is_file()


def test_aria_package_remains_zero_leakage():
    tokens = (
        "DE01_pseudobulk_de",
        "NRM01_sparse_normalization",
        "dr_05_pca_preprocessing_sentinels",
        "gse297023_q2_senescence_not_mechanism",
        "tapestri_ccus",
        "compbiobench-data-v1",
        "bix-1-q1",
        "CapsuleFolder-",
        "GSE182797",
    )
    assert scan_codebase_for_benchmark_leakage(REPO_ROOT / "aria", tokens) == []
    assert len(compute_aria_package_sha256()) == 64
