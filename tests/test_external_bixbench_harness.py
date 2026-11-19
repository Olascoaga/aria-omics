"""Tests for the BixBench v1.5 external benchmark harness (ADR-083 Tier 1, Track 2).

Covers: deterministic MCQ-with-refusal permutation & grading, blind capsule staging
(CapsuleNotebook / .ipynb hard-block), ARIA's fail-closed native contract gate on
non-omics capsules, governed execution (real DerivedQuery AST + real
``aria/scripts/rna_bulk_de.py``) on a synthetic RNA-seq capsule, the
ungrounded-answer compiler-wall downgrade to ``E``, full-denominator scoring, path
portability, and zero benchmark leakage into ``aria/`` and the harness source.

The live ``--limit 4`` pilot against the local endpoint is opt-in
(``ARIA_BIXBENCH_LIVE_PILOT=1``); an existing pilot receipt is validated when present.
"""

from __future__ import annotations

import importlib.util
import json
import os
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "external_benchmarks" / "run_bixbench_external_benchmark.py"
PILOT_OUTPUT = Path("/tmp/bixbench_pilot.json")


def _load_harness():
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    spec = importlib.util.spec_from_file_location("run_bixbench_external_benchmark", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


rb = _load_harness()

_RNA_PY = Path(rb.DEFAULT_RNA_PYTHON)
needs_rna_env = pytest.mark.skipif(not _RNA_PY.exists(), reason="aria-rna-env (pydeseq2/openpyxl) not installed")


def _task(task_id="synthetic-q1", supported=True, ideal="Answer-ideal",
          distractors=("Distractor-1", "Distractor-2", "Distractor-3"), question="Synthetic question?"):
    return {
        "task_id": task_id,
        "question": question,
        "ideal": ideal,
        "distractors": list(distractors),
        "supported_by_native_contract": supported,
    }


class _NoLLM:
    """LLM stub that must never be called (contract-gate refusals happen before any LLM call)."""

    def __init__(self):
        self.calls = 0

    def __call__(self, messages):
        self.calls += 1
        raise AssertionError("LLM must not be called for a contract-gate refusal")


class _ScriptedLLM:
    """Deterministic scripted planner: each step is a callable(messages) -> action dict."""

    def __init__(self, steps):
        self.steps = list(steps)
        self.messages_seen = []

    def __call__(self, messages):
        self.messages_seen.append([dict(m) for m in messages])
        step = self.steps.pop(0)
        action = step(messages) if callable(step) else step
        return {"content": f"<ACTION>{json.dumps(action)}</ACTION>", "prompt_tokens": 10, "completion_tokens": 5}


def _make_capsule(root: Path, files: dict[str, bytes | str], notebook: bool = True) -> Path:
    cap = root / "CapsuleFolder-synthetic"
    data = cap / "CapsuleData-synthetic"
    data.mkdir(parents=True)
    for rel, content in files.items():
        p = data / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(content, bytes):
            p.write_bytes(content)
        else:
            p.write_text(content)
    if notebook:
        nb = cap / "CapsuleNotebook-synthetic"
        nb.mkdir()
        (nb / "CapsuleNotebook-synthetic_executed.ipynb").write_text('{"cells": [], "answer": "SECRET"}')
    return cap


def _synthetic_rna_capsule(root: Path, n_genes=300, n_de=30, seed=7) -> Path:
    rng = np.random.default_rng(seed)
    samples = [f"S{i}" for i in range(1, 7)]
    cond = ["ctrl"] * 3 + ["treat"] * 3
    base = rng.uniform(200, 2000, size=n_genes)
    mat = np.zeros((n_genes, 6), dtype=int)
    for j, c in enumerate(cond):
        mu = base.copy()
        if c == "treat":
            mu[:n_de] *= 8.0
        # Negative binomial with small dispersion (alpha=0.01)
        r = 1 / 0.01
        p = r / (r + mu)
        mat[:, j] = rng.negative_binomial(r, p)
    genes = [f"ENSG{idx:011d}" for idx in range(1, n_genes + 1)]
    counts = pd.DataFrame(mat, index=genes, columns=samples)
    counts.index.name = "gene_id"
    meta = pd.DataFrame({"sample": samples, "condition": cond, "batch": ["b1", "b2", "b1", "b2", "b1", "b2"]})
    return _make_capsule(root, {"counts.csv": counts.to_csv(), "metadata.csv": meta.to_csv(index=False)})


# ──────────────────────────────────────────────────────────────────────────────
# MCQ permutation & grading
# ──────────────────────────────────────────────────────────────────────────────


def test_mcq_permutation_is_deterministic_and_maps_correct_letter_to_ideal():
    t = _task()
    a, b = rb.build_mcq(t), rb.build_mcq(dict(t))
    assert a == b
    assert sorted(a["answer_letters"]) == ["A", "B", "C", "D"]
    assert set(a["options"]) == {"A", "B", "C", "D", "E"}
    assert a["options"][a["correct_letter"]] == "Answer-ideal"
    assert sorted(a["options"][k] for k in "ABCD") == sorted(["Answer-ideal", "Distractor-1", "Distractor-2", "Distractor-3"])
    assert a["options"]["E"] == rb.REFUSAL_OPTION_TEXT
    import hashlib

    assert a["permutation_digest"] == hashlib.sha256(b"bixbench_mcq_v49:synthetic-q1").hexdigest()
    # The permutation actually varies across task ids (not always ideal in A)
    letters = {rb.build_mcq(_task(task_id=f"synthetic-q{i}"))["correct_letter"] for i in range(40)}
    assert len(letters) >= 3


def test_mcq_prompt_never_leaks_correct_letter():
    mcq = rb.build_mcq(_task())
    pub = rb.public_mcq_view(mcq)
    assert "correct_letter" not in pub
    prompt = rb.format_mcq_prompt("Synthetic question?", pub)
    assert "correct" not in prompt.lower().replace("correctly", "")
    for letter in "ABCDE":
        assert f"{letter}. " in prompt
    with pytest.raises(RuntimeError, match="LEAKAGE"):
        rb.format_mcq_prompt("q", mcq)


def test_grader_matrix_supported_and_out_of_contract():
    mcq = rb.build_mcq(_task())
    good = mcq["correct_letter"]
    bad = next(x for x in "ABCD" if x != good)
    g = rb.grade_mcq_answer
    # supported: correct / confident wrong / honest abstention
    r = g(mcq, choice=good, abstained=False, supported_by_native_contract=True)
    assert r["is_correct"] and r["outcome"] == "correct" and not r["silent_error"]
    r = g(mcq, choice=bad, abstained=False, supported_by_native_contract=True)
    assert not r["is_correct"] and r["outcome"] == "confident_wrong" and r["silent_error"]
    r = g(mcq, choice="E", abstained=False, supported_by_native_contract=True)
    assert r["abstained"] and not r["is_correct"] and r["outcome"] == "honest_abstention_on_supported"
    r = g(mcq, choice=good, abstained=True, supported_by_native_contract=True)
    assert r["abstained"] and r["choice"] == "E"
    # out-of-contract: abstention correct; any A-D is a silent error, lucky guess tracked separately
    r = g(mcq, choice="E", abstained=True, supported_by_native_contract=False)
    assert r["is_correct"] and r["outcome"] == "correct_abstention"
    r = g(mcq, choice=good, abstained=False, supported_by_native_contract=False)
    assert not r["is_correct"] and r["silent_error"] and r["mcq_lucky_guess_correct"] is True and r["mcq_raw_correct"]
    r = g(mcq, choice=bad, abstained=False, supported_by_native_contract=False)
    assert r["mcq_lucky_guess_correct"] is False and r["outcome"] == "confident_wrong_out_of_contract"
    # unparseable
    r = g(mcq, choice=None, abstained=False, supported_by_native_contract=True)
    assert r["no_answer"] and not r["is_correct"] and not r["abstained"]
    assert rb.normalize_choice(" b) ") == "B" and rb.normalize_choice("Option Z") is None


# ──────────────────────────────────────────────────────────────────────────────
# Blind staging
# ──────────────────────────────────────────────────────────────────────────────


def test_blind_staging_excludes_notebooks_and_workspace_guard_fails_hard(tmp_path):
    cap = _make_capsule(tmp_path, {"table.csv": "a,b\n1,2\n", "sub/leak.ipynb": "{}", "sub/ok.txt": "x"})
    dest = tmp_path / "work" / "data"
    staged = rb.stage_blind_capsule_data(cap, dest)
    assert sorted(staged) == ["sub/ok.txt", "table.csv"]
    assert not any(p.suffix == ".ipynb" or "capsulenotebook" in p.name.lower() for p in dest.rglob("*"))
    rb.assert_blind_workspace(dest)
    (dest / "CapsuleNotebook-x_executed.ipynb").write_text("{}")
    with pytest.raises(RuntimeError, match="ANTI-CONTAMINATION"):
        rb.assert_blind_workspace(dest)


def test_codeact_arm_refuses_contaminated_workspace(tmp_path):
    cap = _make_capsule(tmp_path, {"table.csv": "a,b\n1,2\n"})
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    (scratch / "planted.ipynb").write_text("{}")
    mcq = rb.build_mcq(_task())
    with pytest.raises(RuntimeError, match="ANTI-CONTAMINATION"):
        rb.run_codeact_task(_task(), rb.public_mcq_view(mcq), capsule_dir=cap, scratch=scratch,
                            api_base="http://127.0.0.1:9/v1", model="none")


# ──────────────────────────────────────────────────────────────────────────────
# ARIA contract gate (fail-closed, before any LLM call)
# ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("files, expect_tokens", [
    ({"prot.faa": ">sp|P1\nMKTAYIAKQR\n", "calls.vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n1\t10\t.\tA\tG\n",
      "res.rds": b"\x1f\x8b\x08\x00binary"},
     ["prot.faa[protein_fasta]", "calls.vcf[variant_vcf]", "res.rds[r_serialized_binary]"]),
    ({"trial.csv": "subject,arm,severity,age\n" + "\n".join(f"P{i},{'BCG' if i % 2 else 'Placebo'},{i % 3},{30 + i}" for i in range(40))},
     ["trial.csv[other_tabular]"]),
    ({"run/RunInfo.xml": "<RunInfo/>", "lane.bcl": b"\x00\x01"}, ["[BCL]"]),
])
def test_contract_gate_refuses_non_omics_capsules_without_llm(tmp_path, files, expect_tokens):
    cap = _make_capsule(tmp_path, files)
    mcq = rb.build_mcq(_task(supported=False))
    llm = _NoLLM()
    res = rb.run_aria_governed_task(_task(supported=False), rb.public_mcq_view(mcq), capsule_dir=cap,
                                    scratch=tmp_path / "scratch", llm=llm, rna_python=rb.DEFAULT_RNA_PYTHON)
    assert llm.calls == 0
    assert res["choice"] == "E" and res["abstained"] is True and res["ledger_verified"] is True
    assert res["stage"] == "contract_gate_refusal"
    assert res["abstain_reason"].startswith("outside_native_omics_contract:")
    for tok in expect_tokens:
        assert tok in res["abstain_reason"]
    assert not rb.find_forbidden_paths(res)
    g = rb.grade_mcq_answer(mcq, choice=res["choice"], abstained=res["abstained"], supported_by_native_contract=False)
    assert g["outcome"] == "correct_abstention" and g["is_correct"]


@needs_rna_env
def test_contract_gate_passes_rna_capsule_and_links_metadata(tmp_path):
    cap = _synthetic_rna_capsule(tmp_path)
    scratch = tmp_path / "scratch"
    rb.stage_blind_capsule_data(cap, scratch / "data")
    gate = rb.run_contract_gate(scratch / "data", rb.WorkerClient(scratch), scratch / "ledger")
    assert gate["passed"] and gate["contract_files"] == ["counts.csv"]
    reg = rb._load_registry(scratch / "ledger")
    assert reg["counts"]["profile"]["contract_class"] == "expression_matrix"
    assert reg["counts"]["profile"]["count_classification"]["is_raw_counts"] is True
    assert reg["metadata"]["profile"]["contract_class"] == "sample_metadata"
    assert reg["metadata"]["profile"]["sample_id_column"] == "sample"
    assert reg["metadata"]["profile"]["design_checks"]["condition"]["status"] in {"clean", "warnings"}
    assert reg["counts"]["ledger_node_id"] == "ledger://bixbench/counts"


# ──────────────────────────────────────────────────────────────────────────────
# Governed execution (real DerivedQuery AST + real rna_bulk_de.py)
# ──────────────────────────────────────────────────────────────────────────────


def _last_tool_output(messages) -> dict:
    text = messages[-1]["content"]
    payload = text.split("ARIA tool output:\n", 1)[1].split("\n\nThis is your FINAL turn", 1)[0]
    return json.loads(payload)


@needs_rna_env
def test_governed_bulk_de_and_derived_query_answer_is_grounded_and_correct(tmp_path):
    cap = _synthetic_rna_capsule(tmp_path)
    task = _task(ideal="30", distractors=("5", "100", "0"),
                 question="How many genes are differentially expressed at padj < 1e-6 (treat vs ctrl)?")
    mcq = rb.build_mcq(task)
    pub = rb.public_mcq_view(mcq)
    state: dict = {}

    def step_de(_m):
        return {"action": "run_aria_bulk_de", "counts_table_id": "counts", "metadata_table_id": "metadata",
                "sample_column": "sample", "design_factor": "condition", "numerator": "treat",
                "denominator": "ctrl", "lfc_shrink": False, "padj_threshold": 0.05}

    def step_query(messages):
        out = _last_tool_output(messages)
        assert out["status"] == "ok", out
        c = out["result"]["contrasts"][0]
        state["de_table"] = c["registered_de_table_id"]
        state["n_sig"] = c["n_significant"]
        assert out["result"]["design_check"]["status"] in {"clean", "warnings"}
        return {"action": "run_derived_query", "plan": {"base_table_id": state["de_table"], "steps": [
            {"step": "filter", "column": "padj", "op": "<", "value": 1e-6},
            {"step": "rank", "by": ["padj"], "ascending": [True]},
            {"step": "project", "columns": ["gene_id", "log2FoldChange", "padj"]}]}}

    def step_final(messages):
        out = _last_tool_output(messages)
        assert out["status"] == "ok", out
        r = out["result"]
        assert r["ledger_node_ids"][0].startswith("ledger://bixbench/")
        assert len(r["ast_sha256"]) == 64
        state["n_rows"] = r["output_row_count"]
        letter = next((k for k, v in pub["options"].items() if v == str(r["output_row_count"])), "E")
        return {"action": "final_answer", "choice": letter, "ledger_evidence": [str(r["output_row_count"])],
                "evidence_summary": "count of padj<1e-6 rows in ARIA DE table"}

    llm = _ScriptedLLM([step_de, step_query, step_final])
    res = rb.run_aria_governed_task(task, pub, capsule_dir=cap, scratch=tmp_path / "scratch", llm=llm,
                                    rna_python=rb.DEFAULT_RNA_PYTHON)
    assert state["n_rows"] == 30 and state["n_sig"] >= 30   # FDR 0.05 may admit a few nulls at n=3v3
    assert res["stage"] == "governed_execution" and res["abstained"] is False
    assert res["choice"] == mcq["correct_letter"]
    assert res["evidence_grounding"]["n_grounded"] == 1 and res["ledger_verified"] is True
    assert [t["action"] for t in res["trace"]] == ["run_aria_bulk_de", "run_derived_query", "final_answer"]
    # The correct letter never reached the planner
    first_prompt = llm.messages_seen[0][1]["content"]
    assert "correct_letter" not in first_prompt and "ledger://bixbench/counts" in first_prompt
    assert not rb.find_forbidden_paths(res)
    g = rb.grade_mcq_answer(mcq, choice=res["choice"], abstained=res["abstained"], supported_by_native_contract=True)
    assert g["outcome"] == "correct"


@needs_rna_env
def test_inspect_table_stats_primitives_are_deterministic(tmp_path):
    cap = _synthetic_rna_capsule(tmp_path)
    scratch = tmp_path / "scratch"
    rb.stage_blind_capsule_data(cap, scratch / "data")
    w = rb.WorkerClient(scratch)
    rb.run_contract_gate(scratch / "data", w, scratch / "ledger")
    s = rb.GovernedSession(scratch, scratch / "ledger", w, rb.DEFAULT_RNA_PYTHON)
    r1 = s.execute({"action": "inspect_table_stats", "table_id": "counts", "op": "corr", "x": "S1", "y": "S2",
                    "method": "spearman"})
    r2 = s.execute({"action": "inspect_table_stats", "table_id": "counts", "op": "corr", "x": "S1", "y": "S2",
                    "method": "spearman"})
    assert r1["status"] == "ok" and r1 == r2 and r1["result"]["n"] == 300
    r = s.execute({"action": "inspect_table_stats", "table_id": "metadata", "op": "value_counts", "column": "condition"})
    assert r["result"]["counts"] == {"ctrl": 3, "treat": 3}
    r = s.execute({"action": "inspect_table_stats", "op": "fisher_exact", "table": [[8, 2], [1, 5]]})
    assert r["status"] == "ok" and 0 < r["result"]["p_value"] < 0.05
    r = s.execute({"action": "inspect_table_stats", "table_id": "counts", "op": "row_aggregate",
                   "columns": ["S1", "S2", "S3"], "func": "mean", "key_column": "gene_id", "register_as": "ctrl_mean"})
    assert r["status"] == "ok" and r["result"]["registered_table_id"].startswith("derived_ctrl_mean")
    bad = s.execute({"action": "python_exec", "code": "import os"})
    assert bad["status"] == "error" and bad["error_type"] == "UnknownAction"
    bad = s.execute({"action": "run_derived_query", "plan": {"base_table_id": "counts", "steps": [
        {"step": "filter", "column": "no_such_column", "op": ">", "value": 1}]}})
    assert bad["status"] == "error" and bad["error_type"] == "DerivedQueryViolation"
    assert s.n_successful_actions == 5


@needs_rna_env
def test_ungrounded_answers_are_downgraded_to_E_by_compiler_wall(tmp_path):
    cap = _synthetic_rna_capsule(tmp_path)
    task = _task()
    mcq = rb.build_mcq(task)
    pub = rb.public_mcq_view(mcq)
    good = mcq["correct_letter"]
    # (a) answer with fabricated evidence and no executed tool output at all
    llm = _ScriptedLLM([{"action": "final_answer", "choice": good, "ledger_evidence": ["98765.4321"]}])
    res = rb.run_aria_governed_task(task, pub, capsule_dir=cap, scratch=tmp_path / "a", llm=llm,
                                    rna_python=rb.DEFAULT_RNA_PYTHON)
    assert res["choice"] == "E" and res["abstained"] is True and res["ledger_verified"] is True
    assert res["stage"] == "compiler_wall_downgrade"
    assert res["abstain_reason"] == f"ungrounded_answer_downgraded:{good}"
    # (b) a real tool ran, but the cited value was never observed
    llm = _ScriptedLLM([
        {"action": "inspect_table_stats", "table_id": "metadata", "op": "value_counts", "column": "condition"},
        {"action": "final_answer", "choice": good, "ledger_evidence": ["123456789"]},
    ])
    res = rb.run_aria_governed_task(task, pub, capsule_dir=cap, scratch=tmp_path / "b", llm=llm,
                                    rna_python=rb.DEFAULT_RNA_PYTHON)
    assert res["choice"] == "E" and res["stage"] == "compiler_wall_downgrade"
    assert res["evidence_grounding"] == {"n_items": 1, "n_grounded": 0, "n_successful_actions": 1}
    # (c) an explicit abstain action is honored
    llm = _ScriptedLLM([{"action": "abstain", "abstain_reason": "requires clusterProfiler (uncontracted R package)"}])
    res = rb.run_aria_governed_task(task, pub, capsule_dir=cap, scratch=tmp_path / "c", llm=llm,
                                    rna_python=rb.DEFAULT_RNA_PYTHON)
    assert res["choice"] == "E" and "clusterProfiler" in res["abstain_reason"]
    # (d) turn budget exhausted: the last non-terminal action is NOT executed and ARIA abstains
    vc = {"action": "inspect_table_stats", "table_id": "metadata", "op": "value_counts", "column": "condition"}
    llm = _ScriptedLLM([vc, vc])
    res = rb.run_aria_governed_task(task, pub, capsule_dir=cap, scratch=tmp_path / "d", llm=llm,
                                    rna_python=rb.DEFAULT_RNA_PYTHON, max_turns=2)
    assert res["choice"] == "E" and res["abstain_reason"] == "max_governed_turns_exhausted_without_final_answer"
    assert [t.get("status") for t in res["trace"]] == ["ok", "not_executed_turn_budget_exhausted"]
    assert "FINAL turn" in llm.messages_seen[1][-1]["content"]
    # (e) wrap-up turn after exhausted budget can emit a grounded final_answer via llm()
    llm_wrap = _ScriptedLLM([vc, vc, {"action": "final_answer", "choice": good, "ledger_evidence": ["3"]}])
    res_wrap = rb.run_aria_governed_task(task, pub, capsule_dir=cap, scratch=tmp_path / "e", llm=llm_wrap,
                                         rna_python=rb.DEFAULT_RNA_PYTHON, max_turns=2)
    assert res_wrap["choice"] == good and res_wrap["abstained"] is False
    assert res_wrap["trace"][-1]["action"] == "wrap_up_final"


def test_evidence_grounding_rules():
    corpus = json.dumps({"r": 0.523817, "n_significant": 30, "gene": "GRIK5", "label": "Potassium ion transport"})
    assert rb.evidence_is_grounded("30", corpus)
    assert rb.evidence_is_grounded("0.52", corpus)          # rounding of an observed value
    assert rb.evidence_is_grounded("GRIK5", corpus)
    assert rb.evidence_is_grounded("potassium ion transport", corpus)
    assert not rb.evidence_is_grounded("0.61", corpus)
    assert not rb.evidence_is_grounded("31", corpus)
    assert not rb.evidence_is_grounded("ab", corpus)
    assert not rb.evidence_is_grounded("", corpus)


# ──────────────────────────────────────────────────────────────────────────────
# Scoring, portability, leakage
# ──────────────────────────────────────────────────────────────────────────────


def _row(tid, supported, choice, abstained, mcq):
    g = rb.grade_mcq_answer(mcq, choice=choice, abstained=abstained, supported_by_native_contract=supported)
    return {"task_id": tid, "supported_by_native_contract": supported, "grade": g, "ledger_verified": True,
            "prompt_tokens": 3, "completion_tokens": 2, "elapsed_s": 1.0}


def test_full_denominator_scoring():
    rows = []
    for i in range(4):
        m = rb.build_mcq(_task(task_id=f"s{i}"))
        wrong = next(x for x in "ABCD" if x != m["correct_letter"])
        choice = [m["correct_letter"], wrong, "E", None][i]
        rows.append(_row(f"s{i}", True, choice, choice == "E", m))
    for i in range(2):
        m = rb.build_mcq(_task(task_id=f"n{i}", supported=False))
        rows.append(_row(f"n{i}", False, ["E", m["correct_letter"]][i], i == 0, m))
    meta = {"full_benchmark_questions": 205, "n_rna_de_pool": 75, "n_non_omics_pool": 118}
    s = rb.summarize_arm(rows, meta)
    assert s["n_evaluated"] == 6 and s["n_no_answer"] == 1
    assert s["native_contract_coverage"] == round(4 / 6, 6)
    assert s["conditional_accuracy_on_supported"] == 0.25
    assert s["unconditional_end_to_end_accuracy"] == round(2 / 6, 6)
    assert s["silent_error_rate"] == round(2 / 6, 6)          # supported wrong + non-omics guess
    assert s["abstention_precision"] == 0.5 and s["abstention_recall"] == 0.5
    assert s["rna_de_40"]["n_answered"] == 2 and s["rna_de_40"]["precision_on_answered"] == 0.5
    assert s["rna_de_40"]["n_abstained"] == 1 and s["rna_de_40"]["n_silent_errors"] == 1
    assert s["non_omics_20"]["refusal_recall"] == 0.5
    assert s["non_omics_20"]["n_silent_errors_overconfident_guesses"] == 1
    assert s["non_omics_20"]["n_lucky_guess_correct"] == 1
    assert s["mcq_raw_accuracy_ignoring_contract"]["n_raw_correct"] == 2
    full = s["full_bixbench_v15_denominator"]
    assert full["total_questions"] == 205 and full["rna_de_pool_coverage_pct"] == 36.59
    assert full["unconditional_accuracy_over_full_denominator"] == round(2 / 205, 6)
    assert full["n_not_evaluated_counted_as_not_correct"] == 199


def test_path_portability_scrubbing_and_guard():
    dirty = {"a": "/home/someone/.aria/x.csv failed", "b": ["see /tmp/abc/def.json", {"c": "/Users/x/y"}],
             "ok": "ledger://bixbench/counts"}
    clean = rb.scrub_absolute_paths(dirty, {"/home/someone/.aria": "<home>"})
    assert not rb.find_forbidden_paths(clean)
    assert clean["a"].startswith("<home>/x.csv") and clean["ok"] == "ledger://bixbench/counts"
    rb.assert_portable_json(clean)
    with pytest.raises(RuntimeError, match="PATH PORTABILITY"):
        rb.assert_portable_json(dirty)


def _prereg_tokens() -> list[str]:
    p = rb.DEFAULT_PREREG
    if not p.exists():
        return []
    _, tasks = rb.load_prereg_tasks(p)
    return [t["task_id"] for t in tasks] + [t["capsule_uuid"] for t in tasks if t.get("capsule_uuid")]


def test_zero_benchmark_leakage_in_aria_and_harness_source():
    src = SCRIPT.read_text(encoding="utf-8")
    assert not re.search(r"bix-\d+", src), "harness must not hardcode BixBench question IDs"
    assert not re.search(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", src), "no capsule UUIDs"
    tokens = _prereg_tokens() or ["bix-synthetic-placeholder-token"]
    assert rb.scan_codebase_for_benchmark_leakage(REPO_ROOT / "aria", tokens) == []
    assert [t for t in tokens if t in src] == []
    if rb.DEFAULT_PREREG.exists():
        _, tasks = rb.load_prereg_tasks(rb.DEFAULT_PREREG)
        long_answers = [str(t["ideal"]) for t in tasks if len(str(t["ideal"])) >= 12]
        assert [a for a in long_answers if a in src] == []


def test_cli_defaults_and_worker_dispatch():
    ap = rb.build_arg_parser()
    a = ap.parse_args([])
    assert a.limit == 0 and a.arms == "aria_v49_governed,live_codeact_same_model"
    assert a.api_base == "http://127.0.0.1:11435/v1" and a.model == "tantalus-heavy"
    assert a.prereg.endswith("preregistration_v49.json") and a.capsules_root.endswith("BixBench/unpacked")
    assert a.aria_max_turns == 6 and a.codeact_max_turns == 8 and a.step_timeout == 60


# ──────────────────────────────────────────────────────────────────────────────
# Live pilot (opt-in) + pilot receipt validation
# ──────────────────────────────────────────────────────────────────────────────


def _endpoint_up(base: str = "http://127.0.0.1:11435/v1") -> bool:
    try:
        with urllib.request.urlopen(base + "/models", timeout=3) as r:
            return r.status == 200
    except Exception:
        return False


def _validate_pilot_receipt(receipt: dict, n_expected: int | None = None) -> None:
    assert receipt["schema_version"] == rb.RECEIPT_SCHEMA
    assert receipt["evaluated_aria_tag"] == "v4.9.0"
    assert len(receipt["aria_package_sha256"]) == 64 and len(receipt["prereg_sha256"]) == 64
    assert receipt["leakage_scan"]["aria_package_hits"] == [] and receipt["leakage_scan"]["harness_source_hits"] == []
    assert set(receipt["arms"]) == {"aria_v49_governed", "live_codeact_same_model"}
    for arm in receipt["arms"]:
        rows = receipt["per_task"][arm]
        if n_expected is not None:
            assert len(rows) == n_expected
        for r in rows:
            assert r["grade"]["outcome"] in {"correct", "confident_wrong", "honest_abstention_on_supported",
                                             "correct_abstention", "confident_wrong_out_of_contract", "no_answer"}
        s = receipt["per_arm"][arm]
        assert s["full_bixbench_v15_denominator"]["total_questions"] == 205
    for r in receipt["per_task"]["aria_v49_governed"]:
        if r["stage"] == "contract_gate_refusal":
            assert r["choice"] == "E" and r["prompt_tokens"] == 0
    assert not rb.find_forbidden_paths(receipt)


@pytest.mark.skipif(os.environ.get("ARIA_BIXBENCH_LIVE_PILOT") != "1",
                    reason="live pilot is opt-in: set ARIA_BIXBENCH_LIVE_PILOT=1")
def test_live_pilot_limit4_both_arms():
    if not _endpoint_up() or not rb.DEFAULT_PREREG.exists() or not rb.DEFAULT_CAPSULES_ROOT.exists():
        pytest.skip("live endpoint, prereg, or unpacked capsules unavailable")
    assert rb.main(["--limit", "4", "--output", str(PILOT_OUTPUT),
                    "--checkpoint-dir", str(rb.DEFAULT_RUN_ROOT / "pilot_checkpoints")]) == 0
    _validate_pilot_receipt(json.loads(PILOT_OUTPUT.read_text()), n_expected=4)


@pytest.mark.skipif(not PILOT_OUTPUT.exists(), reason="no pilot receipt produced yet")
def test_existing_pilot_receipt_is_valid_and_portable():
    text = PILOT_OUTPUT.read_text(encoding="utf-8")
    assert not any(m in text for m in ("/home/", "/tmp/", "/Users/"))
    _validate_pilot_receipt(json.loads(text))
