"""Synthetic contract falsification and blind external-harness integration tests."""
from __future__ import annotations

import gzip
import json
from pathlib import Path
import subprocess

import pandas as pd
import pytest

from scripts.external_benchmarks import run_genotex_external_benchmark as h


@pytest.fixture
def matrix_factory(tmp_path):
    def make(name="fixture_series_matrix.txt.gz", *, groups=None, genders=None,
             features=60, all_nan=False, table=True, organism="Homo sapiens", title="Clinical cohort",
             ages=None, sample_prefix="sample", extra_rows=()):
        groups = groups if groups is not None else ["case"] * 4 + ["control"] * 4
        n = len(groups)
        genders = genders if genders is not None else ["Female", "Male"] * (n // 2)
        ages = ages if ages is not None else [str(20 + i % 4) for i in range(n)]
        samples = [f"{sample_prefix}{i}" for i in range(n)]
        rows = [["!Series_title", title], ["!Series_summary", "Patient gene expression"],
                ["!Series_overall_design", "Case-control tissue cohort"],
                ["!Sample_geo_accession", *samples], ["!Sample_title", *samples],
                ["!Sample_source_name_ch1", *(["blood"] * n)],
                ["!Sample_organism_ch1", *([organism] * n)],
                ["!Sample_characteristics_ch1", *[f"diagnosis: {x}" for x in groups]],
                ["!Sample_characteristics_ch1", *[f"age: {x}" for x in ages]],
                ["!Sample_characteristics_ch1", *[f"gender: {x}" for x in genders]], *extra_rows]
        path = tmp_path / name
        import csv
        with gzip.open(path, "wt") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerows(rows)
            if table:
                f.write("!series_matrix_table_begin\n")
                writer.writerow(["ID_REF", *samples])
                for i in range(features):
                    writer.writerow([f"gene{i}", *(["NaN"] * n if all_nan else [100 + i + j for j in range(n)])])
                f.write("!series_matrix_table_end\n")
        return path
    return make


def gate(path, **proposal):
    study = h.parse_series_matrix(path)
    proposed = {"is_human_clinical_cohort": True, "trait_column": "diagnosis",
                "age_column": "age", "gender_column": "gender", **proposal}
    return h.mechanical_gate(study, proposed)


def test_valid_clinical_design_uses_real_aria_contract(matrix_factory):
    path = matrix_factory()
    study = h.parse_series_matrix(path)
    result = h.mechanical_gate(study, h.deterministic_proposal(study, "case"))
    assert result["is_available"] and result["has_age"] and result["has_gender"]
    assert not result["abstained"]
    assert result["design_validation"]["status"] == "clean"
    assert study.expression["n_features"] == 60
    assert study.expression["n_samples_in_matrix"] == 8
    assert study.expression["numeric_non_null_fraction"] == 1
    assert study.expression["classification"]["kind"] in ("raw", "integer_small")


def test_governed_constant_gender_overclaim_is_vetoed(matrix_factory, monkeypatch):
    study = h.parse_series_matrix(matrix_factory(genders=["Female"] * 8))
    def mock_chat(messages, **kwargs):
        assert kwargs["temperature"] == 0
        payload = json.loads(messages[1]["content"])
        assert set(payload) == {"target_trait", "metadata"}
        assert "ground_truth" not in json.dumps(payload)
        return {"content": json.dumps({"is_human_clinical_cohort": True,
                "trait_column": "diagnosis", "trait_mapping_summary": "binary_case_control",
                "age_column": "age", "gender_column": "gender", "reasoning": "Sex is recorded."}),
                "prompt_tokens": 10, "completion_tokens": 20, "elapsed_s": 0.1}
    monkeypatch.setattr(h, "chat_completion", mock_chat)
    pred, usage = h.governed_prediction(study, "case", api_base="http://localhost/v1", model="mock")
    assert pred["is_available"]
    assert pred["has_gender"] is False
    assert pred["has_age"] is True
    assert "has_gender" in pred["constant_covariate_traps"]
    assert usage["completion_tokens"] == 20


@pytest.mark.parametrize("kwargs", [
    {"groups": ["case"] * 8}, {"groups": [""] * 8}, {"features": 0},
    {"all_nan": True}, {"table": False}, {"features": 49},
    {"title": "Human cell-line study"}, {"organism": "Mus musculus"},
    {"groups": ["case"] * 7 + ["control"]},
])
def test_unusable_abstains(matrix_factory, kwargs):
    pred = gate(matrix_factory(**kwargs))
    assert pred["is_available"] is False
    assert pred["abstained"] is True


def test_hallucinated_and_constant_columns_are_vetoed(matrix_factory):
    path = matrix_factory(ages=["40 years"] * 8, genders=["Female"] * 8)
    pred = gate(path, trait_column="invented", age_column="age", gender_column="gender")
    assert not pred["is_trait_available"] and not pred["is_available"]
    assert not pred["has_age"] and not pred["has_gender"]
    assert pred["design_validation"]["status"] == "blocking"
    assert not gate(path, trait_column="gender")["is_trait_available"]
    assert not gate(path)["has_gender"]


def test_missing_covariates_do_not_become_values(matrix_factory):
    pred = gate(matrix_factory(ages=["NA"] * 8, genders=["unknown"] * 8))
    assert pred["is_available"]
    assert not pred["has_age"] and not pred["has_gender"]


def test_parser_preserves_mixed_keys_and_raw_rows(matrix_factory):
    row = ["!Sample_characteristics_ch1", *(["other key: a: b; score: 4"] * 4 + ["different key: z"] * 4)]
    raw = ["!Sample_characteristics_ch1", *(["unkeyed value"] * 8)]
    study = h.parse_series_matrix(matrix_factory(extra_rows=[row, raw]))
    assert study.sample_meta.loc[0, "other_key"] == "a: b"
    assert study.sample_meta.loc[0, "score"] == "4"
    assert study.sample_meta.loc[7, "different_key"] == "z"
    assert study.sample_meta.loc[0, "characteristics_ch1_5"] == "unkeyed value"


def test_real_design_validator_blocks_confounded_covariate(matrix_factory):
    pred = gate(matrix_factory(genders=["Female"] * 4 + ["Male"] * 4))
    assert pred["has_gender"] and not pred["is_available"]
    assert any(i["check"] == "condition_covariate_confounding" for i in pred["design_validation"]["issues"])


def test_expression_metadata_alignment_prevents_phantom_controls(matrix_factory):
    path = matrix_factory()
    # Matrix samples are deliberately disjoint from metadata sample IDs.
    text = gzip.open(path, "rt").read()
    begin = text.index("!series_matrix_table_begin")
    with gzip.open(path, "wt") as f:
        f.write(text[:begin] + text[begin:].replace("sample", "other"))
    assert not gate(path)["is_available"]


def test_continuous_trait_is_recovered_but_frozen_de_contract_is_respected(matrix_factory):
    pred = gate(matrix_factory(groups=[str(i) for i in range(8)]))
    assert pred["is_trait_available"]
    assert not pred["is_available"]
    assert any(i["check"] == "numeric_condition_factor" for i in pred["design_validation"]["issues"])


def test_blind_staging_and_forbidden_entries(matrix_factory, tmp_path):
    raw = matrix_factory()
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    paths = h.stage_raw_files([raw], workspace)
    assert [p.name for p in paths] == [raw.name]
    assert not paths[0].stat().st_mode & 0o222
    h.assert_blind_workspace(workspace)
    for bad in ("cohort_info.json", "solution.ipynb", "output/preprocess/answer.txt"):
        p = workspace / bad
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("reference")
        with pytest.raises(RuntimeError, match="ANTI-CONTAMINATION"):
            h.assert_blind_workspace(workspace)
        p.unlink()
    with pytest.raises(RuntimeError):  # Even an empty preprocess directory fails.
        h.assert_blind_workspace(workspace)


def test_path_portability_including_free_text():
    obj = {"path": Path("/tmp/private/file"), "reason": "error at /home/someone/data and /Users/me/file",
           "nested": ["/tmp/a/b", "relative/file"], "endpoint": "http://localhost:11435/v1"}
    payload = json.dumps(h.portable(obj))
    assert all(root not in payload for root in ("/home/", "/tmp/", "/Users/"))
    assert h.portable(obj)["endpoint"] == obj["endpoint"]


def task(task_id, available, **gt):
    return {"task_id": task_id, "ground_truth": {"is_available": available,
            "is_trait_available": available, "has_age": None, "has_gender": None, **gt}}


def test_full_denominator_adapter_silent_errors_and_covariate_traps():
    tasks = [task("usable", True, has_age=True, has_gender=False), task("unusable", False)]
    records = [{"task_id": "usable", "prediction": {**h.unavailable("test"),
               "is_available": True, "has_age": True, "has_gender": False, "is_trait_available": True},
               "constant_covariate_traps": ["has_gender"]},
               {"task_id": "unusable", "prediction": {**h.unavailable("test"), "is_available": True}}]
    scores = h.score_arm(tasks, records)
    assert scores["dataset_filtering"]["TP"] == 1
    assert scores["dataset_filtering"]["FP"] == 1
    assert scores["aria_adapter"]["silent_error_rate"] == .5
    assert scores["covariate_recovery"]["constant_covariate_trap_rejection_rate"] == 1
    assert scores["covariate_recovery"]["age_n"] == 1
    # Missing tasks are not credited as correct abstentions, including in pilots.
    scores = h.score_arm(tasks, records[:1])
    assert scores["aria_adapter"]["total_suite_size"] == 2
    assert scores["aria_adapter"]["unconditional_end_to_end_accuracy"] == .5
    assert scores["dataset_filtering"]["accuracy"] == .5
    records[1]["prediction"] = h.unavailable("unusable")
    assert h.score_arm(tasks, records)["aria_adapter"]["unconditional_end_to_end_accuracy"] == 1


def test_frozen_engine_and_dynamic_preregistration_leakage_guard():
    h.assert_frozen_engine("HEAD")
    assert subprocess.check_output(["git", "diff", "HEAD", "--", "aria"], cwd=h.REPO_ROOT) == b""
    prereg = Path.home() / ".aria/benchmarks/external/preregistration_v49.json"
    tokens = ["synthetic_external_task_sentinel"]
    if prereg.exists():
        tasks = json.loads(prereg.read_text())["genotex_100"]["tasks"]
        tokens += [t[k] for t in tasks for k in ("task_id", "gse")]
        runner = Path(h.__file__).read_text()
        assert not any(t["gse"] in runner or t["trait"] in runner for t in tasks)
    assert h.scan_codebase_for_benchmark_leakage(h.REPO_ROOT / "aria", tokens) == []


def test_cli_checkpoint_resume_and_blind_multiplatform(matrix_factory, tmp_path, monkeypatch):
    root = tmp_path / "source"
    cohort = root / "input/GEO/synthetic/cohort"
    cohort.mkdir(parents=True)
    raw1 = matrix_factory(name="first_series_matrix.txt.gz", groups=["case"] * 8)
    raw2 = matrix_factory(name="second_series_matrix.txt.gz", groups=["control"] * 8)
    for raw in (raw1, raw2):
        (cohort / raw.name).write_bytes(raw.read_bytes())
    (root / "output/preprocess").mkdir(parents=True)
    (root / "output/preprocess/cohort_info.json").write_text("MUST NEVER STAGE")
    manifest = tmp_path / "prereg.json"
    t = {**task("external_fixture_cohort_q7z9", False), "trait": "case", "gse": "synthetic_accession",
         "matrix_files": [p.relative_to(root).as_posix() for p in cohort.iterdir()]}
    manifest.write_text(json.dumps({"genotex_100": {"tasks": [t]}}))
    out = tmp_path / "receipt.json"
    argv = ["--prereg", str(manifest), "--data-root", str(root), "--output", str(out),
            "--arms", h.ARMS[0]]
    assert h.main(argv) == 0
    receipt = json.loads(out.read_text())
    assert not receipt["results"][0]["prediction"]["is_available"]
    assert len(receipt["results"][0]["platform_predictions"]) == 2
    assert receipt["usage"]["llm_calls"] == 0
    def fail(*a, **kw):
        pytest.fail("resume executed an already checkpointed task")
    monkeypatch.setattr(h, "evaluate_task", fail)
    assert h.main([*argv, "--resume"]) == 0
    assert all(s not in out.read_text() for s in ("/home/", "/tmp/", "/Users/"))


def test_worker_isolation_and_same_endpoint_across_all_arms(matrix_factory, tmp_path, monkeypatch):
    from scripts.external_benchmarks.live_codeact_agent import CodeActResult
    root = tmp_path / "root"
    cohort = root / "input/GEO/fixture/cohort"
    cohort.mkdir(parents=True)
    raw = matrix_factory()
    target = cohort / raw.name
    target.write_bytes(raw.read_bytes())
    blind_task = {"task_id": "isolated_fixture_4x8", "trait": "case",
                  "matrix_files": [target.relative_to(root).as_posix()]}
    monkeypatch.setenv("ARIA_LLM_CACHE", "original")
    monkeypatch.setenv("ARIA_MEMORY_DB", "original-db")
    seen_memory = []
    def check_env():
        import os
        assert os.environ["ARIA_LLM_CACHE"] == "0"
        assert os.environ["ARIA_REPRODUCIBLE"] == "1"
        memory = os.environ["ARIA_MEMORY_DB"]
        assert not Path(memory).exists()
        seen_memory.append(memory)
    def chat(messages, **kwargs):
        check_env()
        assert kwargs["api_base"] == "http://same.endpoint/v1"
        assert kwargs["model"] == "same-model"
        return {"content": json.dumps({"is_human_clinical_cohort": True,
                "trait_column": "diagnosis", "age_column": "age", "gender_column": "gender"}),
                "prompt_tokens": 11, "completion_tokens": 22, "elapsed_s": .1}
    def live(prompt, **kwargs):
        check_env()
        assert kwargs["api_base"] == "http://same.endpoint/v1"
        assert kwargs["model"] == "same-model"
        assert kwargs["max_turns"] == 6 and kwargs["step_timeout"] == 45
        workspace = kwargs["work_dir"]
        h.assert_blind_workspace(workspace)
        assert [p.name for p in (workspace / "data").iterdir()] == [raw.name]
        assert not (workspace / "output").exists()
        assert "ground_truth" not in prompt
        pred = {k: True for k in h.FIELDS}
        pred.update(is_biased=False, abstained=False, abstain_reason=None)
        return CodeActResult(pred, "", 2, 1, 0, 33, 44, .2)
    monkeypatch.setattr(h, "chat_completion", chat)
    monkeypatch.setattr(h, "run_live_codeact_task", live)
    import os
    for arm in h.ARMS:
        record = h.evaluate_task(blind_task, arm, str(root), "http://same.endpoint/v1", "same-model")
        assert not record.get("error")
        assert record["prediction"]["is_available"]
        assert os.environ["ARIA_LLM_CACHE"] == "original"
        assert os.environ["ARIA_MEMORY_DB"] == "original-db"
        if arm == h.ARMS[0]:
            assert record["usage"]["llm_calls"] == 0
    assert len(seen_memory) == len(set(seen_memory)) == 2
    assert all(not Path(p).parent.exists() for p in seen_memory)


def test_severe_imbalance_and_nonfinite_expression_columns(matrix_factory):
    pred = gate(matrix_factory(groups=["case"] * 98 + ["control"] * 2))
    assert pred["is_trait_available"] and pred["is_biased"]
    assert not pred["is_available"]
    path = matrix_factory(name="missing_controls_series_matrix.txt.gz", groups=["case"] * 8 + ["control"] * 4)
    text = gzip.open(path, "rt").read()
    lines = []
    for line in text.splitlines():
        if line.startswith("gene"):
            parts = line.split("\t")
            parts[-4:] = ["NaN"] * 4
            line = "\t".join(parts)
        lines.append(line)
    with gzip.open(path, "wt") as f:
        f.write("\n".join(lines) + "\n")
    study = h.parse_series_matrix(path)
    assert study.expression["is_gene_available"]
    assert len(study.expression["usable_samples"]) == 8
    assert not gate(path)["is_trait_available"]


def test_colons_inside_realistic_geo_attribute_names(matrix_factory):
    row = ["!Sample_characteristics_ch1", *(["model (1: EXP, 2: poor-EXP): 2; note: a:b"] * 8)]
    study = h.parse_series_matrix(matrix_factory(extra_rows=[row]))
    assert study.sample_meta.loc[0, "model_1_exp_2_poor_exp"] == "2"
    assert study.sample_meta.loc[0, "note"] == "a:b"
    assert h.portable("./data/raw_series_matrix.txt.gz") == "./data/raw_series_matrix.txt.gz"


def test_nonexpression_assay_cannot_claim_gene_expression(matrix_factory):
    study = h.parse_series_matrix(matrix_factory(extra_rows=[["!Series_type", "Genome variation profiling by SNP array"]]))
    assert study.expression["n_numeric_features"] == 60
    assert not study.expression["is_gene_available"]
    assert not h.mechanical_gate(study, h.deterministic_proposal(study, "case"))["is_available"]


def test_llm_failure_keeps_trap_denominator_and_earns_no_abstention_credit(matrix_factory, tmp_path, monkeypatch):
    cohort = tmp_path / "input/GEO/fixture/cohort"
    cohort.mkdir(parents=True)
    raw = matrix_factory(genders=["Female"] * 8)
    target = cohort / raw.name
    target.write_bytes(raw.read_bytes())
    def fail(*args, **kwargs):
        raise ConnectionError("endpoint unavailable at /tmp/private/socket")
    monkeypatch.setattr(h, "chat_completion", fail)
    record = h.evaluate_task({"task_id": "failure_fixture_7w3", "trait": "case",
        "matrix_files": [target.relative_to(tmp_path).as_posix()]}, h.ARMS[1], str(tmp_path), "http://mock/v1", "mock")
    assert record["error"] and "has_gender" in record["constant_covariate_traps"]
    assert "/tmp/" not in json.dumps(record)
    scores = h.score_arm([task("failure_fixture_7w3", False)], [record])
    assert scores["aria_adapter"]["unconditional_end_to_end_accuracy"] == 0
    assert scores["covariate_recovery"]["constant_covariate_trap_n"] == 1
    assert scores["covariate_recovery"]["constant_covariate_trap_rejection_rate"] == 0


def test_codeact_turn_counter_preserves_full_budget_and_restores_request(monkeypatch):
    observed = []
    def fake(messages, **kwargs):
        observed.append(messages)
        assert kwargs == {"model": "same-model", "temperature": 0}
        return {"content": "mock"}
    monkeypatch.setattr(h.codeact_module, "chat_completion", fake)
    original_messages = [{"role": "system", "content": "agent"}, {"role": "user", "content": "task"}]
    with h.codeact_turn_budget():
        for _ in range(6):
            h.codeact_module.chat_completion(original_messages, model="same-model", temperature=0)
    assert h.codeact_module.chat_completion is fake
    assert "5 possible Python executions remaining" in observed[0][0]["content"]
    assert "FINAL RESPONSE ONLY" in observed[-1][0]["content"]
    assert observed[-1][-1]["role"] == "user"
    assert original_messages[0]["content"] == "agent"


def test_codeact_final_reporting_uses_only_its_own_observations(monkeypatch):
    seen = []
    answer = {k: False for k in h.FIELDS}
    answer.update(abstained=True, abstain_reason="insufficient evidence")
    def chat(messages, **kwargs):
        seen.append(messages)
        return {"content": json.dumps(answer)}
    monkeypatch.setattr(h.codeact_module, "chat_completion", chat)
    messages = [{"role": "system", "content": "coding instructions"},
                {"role": "user", "content": "original raw-file task"},
                {"role": "assistant", "content": "```python\nprint('observation')\n```"},
                {"role": "user", "content": "STDOUT: observation"}]
    with h.codeact_turn_budget():
        for _ in range(6):
            response = h.codeact_module.chat_completion(messages)
    final = seen[-1]
    assert [m["role"] for m in final] == ["system", "user"]
    payload = json.loads(final[1]["content"])
    assert payload["observed_python_results"] == ["STDOUT: observation"]
    assert payload["original_task"] == "original raw-file task"
    assert "print(" not in final[1]["content"]
    assert h.extract_json_answer(response["content"]) == answer


def test_source_symlink_cannot_expose_preprocessed_reference(matrix_factory, tmp_path):
    raw = matrix_factory()
    root = tmp_path / "root"
    refs = root / "output/preprocess"
    cohort = root / "input/GEO/fixture/cohort"
    refs.mkdir(parents=True)
    cohort.mkdir(parents=True)
    reference = refs / "hidden_series_matrix.txt.gz"
    reference.write_bytes(raw.read_bytes())
    link = cohort / "disguised_series_matrix.txt.gz"
    link.symlink_to(reference)
    with pytest.raises(ValueError, match="invalid_or_missing"):
        h.resolve_raw_files({"matrix_files": [link.relative_to(root).as_posix()]}, root)


def test_official_f1_uses_standard_confusion_formula():
    metrics = h.binary_metrics([(True, True), (True, False), (False, True), (False, None)])
    assert metrics["f1"] == 2 * metrics["TP"] / (2 * metrics["TP"] + metrics["FP"] + metrics["FN"])
    assert metrics["accuracy"] == .25
