"""Regression guards for public 11-modality atlas inferential stage contracts.

Covers four multi-assay edge cases surfaced during end-to-end public dataset
benchmarking:
1. Short-read bulk chromatin FASTQ preprocessing (`ChromatinAgent`) passes
   `min_length=20` by default so `2x25 bp` / `28+90 bp` CUT&Tag and CUT&RUN
   libraries are not discarded by RNA-seq's 36 bp fastp floor.
2. `chromatin_diffacc._pseudobulk_da` accepts both `{"test", "reference"}`
   mapping contrasts and `[test, reference]` 2-tuple/list contrasts.
3. `chromatin_qc` routes `scATAC` `.h5ad` peak matrices (e.g. merged multi-donor
   cohort union matrices) through `_h5ad_atac_qc` without fabricating TSSe/FRiP.
4. `multiome_contracts._paired_h5mu_records` inspects on-disk un-audited `.h5mu`
   files so ATAC-only `.h5mu` cohorts are not misclassified as paired Multiome.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


def test_bulk_chromatin_fastq_qc_defaults_to_20bp_min_length(tmp_path):
    pytest.importorskip("litellm")
    from aria.agents.chromatin_agent import ChromatinAgent

    captured: dict = {}

    class _FakeEnv:
        def run_in_stack(self, stack, script_path, params, **_kwargs):
            if script_path.endswith("rna_fastq_qc.py"):
                captured.update(params)
                return {"status": "failed", "samples": []}
            raise AssertionError(f"unexpected script {script_path}")

    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = _FakeEnv()
    agent.publish_status = lambda *_a, **_k: None

    fa = tmp_path / "genome.fa"
    fa.write_text(">chr1\nACGTACGTACGTACGTACGTACGT\n", encoding="utf-8")
    fq = tmp_path / "sample_R1.fastq.gz"
    fq.write_bytes(b"stub")

    agent._run_bulk_atac_preprocessing(
        experiment_id="exp1",
        fastq_files=[str(fq)],
        exp_ctx={"genome": "hg38", "genome_fasta": str(fa)},
        data_type="CUT_AND_TAG",
    )
    assert captured["min_length"] == 20


def test_pseudobulk_da_accepts_two_element_list_contrasts():
    ad = pytest.importorskip("anndata")
    pytest.importorskip("pydeseq2")
    from aria.scripts.chromatin_diffacc import _pseudobulk_da

    rng = np.random.default_rng(7)
    obs_rows = []
    for cond, donors in (("young", ("d1", "d2", "d3")), ("aged", ("d4", "d5", "d6"))):
        for d in donors:
            for i in range(15):
                obs_rows.append({"cell": f"{d}_{i}", "condition": cond, "donor_id": d})
    obs = pd.DataFrame(obs_rows).set_index("cell")
    X = rng.poisson(lam=3.0, size=(len(obs), 25)).astype(np.float32)
    # Inject strong differential signal on the first 5 peaks in "aged"
    aged_mask = (obs["condition"] == "aged").to_numpy()
    X[aged_mask, :5] += 25.0
    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=[f"chr1:{1000 * (i + 1)}-{1000 * (i + 1) + 500}" for i in range(25)]),
    )

    warnings: list[str] = []
    res = _pseudobulk_da(
        adata=adata,
        condition_col="condition",
        replicate_col="donor_id",
        comparisons=[["aged", "young"]],
        padj_max=0.05,
        lfc_min=0.5,
        top_n=10,
        min_cells=5,
        min_reps=3,
        allow_mock=False,
        warnings=warnings,
        n_cpus=1,
    )
    assert res["ran"] is True
    assert res["comparisons"][0]["test"] == "aged"
    assert res["comparisons"][0]["reference"] == "young"
    assert res["comparisons"][0]["status"] == "success"


def test_chromatin_qc_scatac_h5ad_peak_matrix_is_honest(tmp_path):
    ad = pytest.importorskip("anndata")
    from aria.scripts.chromatin_qc import chromatin_qc

    p = tmp_path / "merged_scatac_cohort.h5ad"
    obs = pd.DataFrame(
        {"n_fragments": [1200, 1800, 2400], "pct_mito": [0.02, 0.03, 0.01]},
        index=["c1", "c2", "c3"],
    )
    adata = ad.AnnData(
        X=np.ones((3, 12), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=[f"chr1:{i * 500}-{i * 500 + 300}" for i in range(12)]),
    )
    adata.write_h5ad(p)

    res = chromatin_qc({"data_type": "scATAC", "files": [str(p)], "genome": "hg38"})
    assert res["status"] == "success"
    assert res["input_kind"] == "h5ad"
    assert res["n_cells"] == 3
    assert res["n_peaks"] == 12
    assert res["median_fragments_per_cell"] == 1800.0
    assert res["frip"] is None
    assert res["tss_enrichment"] is None


def test_paired_h5mu_records_excludes_on_disk_atac_only_h5mu(tmp_path, monkeypatch):
    from aria.utils import mudata_io
    from aria.utils.multiome_contracts import infer_multiome_contract

    atac_only = tmp_path / "donor1_scatac.h5mu"
    atac_only.write_bytes(b"nonempty-h5mu-placeholder")

    monkeypatch.setattr(
        mudata_io,
        "read_h5mu_summary",
        lambda path, backed=True: {
            "status": "success",
            "path": str(path),
            "modalities": {"atac": {"n_obs": 100, "n_vars": 500, "obs_columns": []}},
            "atac_modality": "atac",
            "rna_modality": None,
            "is_paired_rna_atac": False,
        },
    )

    contract = infer_multiome_contract(
        {
            "genome": "hg38",
            "modalities": {"scATAC": [str(atac_only)]},
            # Simulate explicit library_types bypassing AssayDetector
            "assay_detections": [],
        }
    )
    assert contract["object_type"] != "paired_mudata"
    assert contract["same_cell"] is False


def test_llm_provider_retries_without_temperature_and_records_honest_provenance(
    tmp_path, monkeypatch
):
    pytest.importorskip("litellm")
    import json
    import aria.llm.provider as prov_mod
    from aria.llm.provider import LLMProvider, ModelConfig, TaskTier

    calls: list[dict] = []

    def _fake_completion(**kwargs):
        calls.append(dict(kwargs))
        if "temperature" in kwargs:
            raise ValueError("Unsupported parameter: 'temperature' is deprecated for this model.")
        return SimpleNamespace(
            choices=[
                SimpleNamespace(
                    message=SimpleNamespace(content="ok-frontier"),
                    finish_reason="stop",
                )
            ],
            usage=SimpleNamespace(prompt_tokens=12, completion_tokens=4, total_tokens=16),
        )

    monkeypatch.setattr(prov_mod, "completion", _fake_completion)
    monkeypatch.setattr(prov_mod, "_MODELS_REJECTING_TEMPERATURE", set())
    monkeypatch.setenv("ARIA_LLM_CACHE", "0")

    usage_log = tmp_path / "llm_usage.jsonl"
    cfg = ModelConfig(
        provider="openai",
        model="openai/wsl-test-frontier",
        context_window=128_000,
        is_local=True,
        api_base="http://127.0.0.1:11436/v1",
        extra_params={"reasoning_effort": "low"},
    )
    llm = LLMProvider(
        models={TaskTier.MEDIUM: [cfg]},
        api_keys={},
        experiment_id="exp-temp-test",
        usage_log=usage_log,
    )
    out = llm.complete("Hello", system="Sys", tier=TaskTier.MEDIUM, max_tokens=64)
    assert out == "ok-frontier"
    assert len(calls) == 2
    assert "temperature" in calls[0]
    assert "temperature" not in calls[1]
    assert calls[1]["api_key"] == "local-loopback-shim"
    assert calls[1]["reasoning_effort"] == "low"

    row = json.loads(usage_log.read_text(encoding="utf-8").strip())
    assert row["temperature"] is None
    assert row["temperature_controlled"] is False
    assert row["deterministic"] is False
    assert row["estimated_cost_usd"] == 0.0


def test_govern_executive_summary_replaces_empty_llm_output_with_fallback():
    pytest.importorskip("litellm")
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent.__new__(NarrativeAgent)
    grouped = {"high": [], "medium": [{"summary": "ok"}], "low": [], "insufficient": []}
    summary, warning, block = agent._build_executive_summary_block(
        *(
            *agent._govern_executive_summary(
                executive_summary="   ",
                grouped_findings=grouped,
                intent={},
                agent_results={},
                narrative_blocks=[],
            ),
            grouped,
            {},
            {},
            {},
            [],
        )
    )
    assert summary.strip() != ""
    assert "empty LLM summary" in warning
    assert block.claim.strip() == summary.strip()


def test_genome_arch_escalations_mark_agent_parameter_checkpoint():
    import ast
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "aria" / "agents" / "genome_arch_agent.py"
    tree = ast.parse(src.read_text(encoding="utf-8"))
    escalation_contexts = []
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "publish_escalation"
        ):
            for kw in node.keywords:
                if kw.arg == "context" and isinstance(kw.value, ast.Dict):
                    keys = [
                        k.value for k in kw.value.keys if isinstance(k, ast.Constant)
                    ]
                    escalation_contexts.append(keys)
    assert len(escalation_contexts) >= 3
    for keys in escalation_contexts:
        assert "agent_parameter_checkpoint" in keys
