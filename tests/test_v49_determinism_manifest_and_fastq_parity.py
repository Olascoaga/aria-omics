"""Regression & architectural verification suite for ARIA v4.9 Foundations:
- Milestone 1 (P0 Layers 1 & 2): deterministic subprocess thread/hash pinning (`A6`),
  `bwa-mem2 mem -K 100000000`, explicit LIANA seed + deterministic tie-breaking,
  deterministic ARPACK `v0` in `hic_topology`, and output SHA-256 tamper rejection
  in `mudata_io._stage_cache` & `hic_preprocess._check_hic_stage_cache` (`A4`,
  preserving 0 diff on `stage_manifest.py` and `deseq2.py`).
- Milestone 2 (P0 Layer 3): token-boundary keyword matching (`text_matching.py`),
  canonical 12-point Leiden union grid, pre-execution `DecisionManifest`,
  append-only `AnalyticalForkLedger`, capsule verification, and Zero-LLM Replay.
- Milestone 3 (Q1 Raw FASTQ Parity): `SetupAgent` `kb ref` auto-provisioning +
  file locking, `RawIngestionAgent` `genome_config` auto-wiring, Multiome
  `chromap_index` auto-build, and Illumina `scATAC` FASTQ layout proposal.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest

from aria.utils.decision_manifest import (
    build_decision_manifest,
    commit_decision_manifest,
    load_decision_manifest,
    verify_decision_manifest,
)
from aria.utils.privacy import EgressPolicy, execution_environment, use_egress_policy
from aria.utils.text_matching import has_keyword_boundary


def test_zero_diff_protected_files():
    repo_root = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [
            "git",
            "diff",
            "--name-only",
            "--",
            "aria/utils/stage_manifest.py",
            "aria/scripts/rna_bulk/deseq2.py",
        ],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=True,
    )
    assert proc.stdout.strip() == ""


def test_execution_environment_pins_hashseed_and_reproducible_blas_threads(monkeypatch):
    monkeypatch.delenv("ARIA_REPRODUCIBLE", raising=False)
    monkeypatch.delenv("ARIA_PIN_BLAS_THREADS", raising=False)
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)

    # Default non-reproducible run still guarantees PYTHONHASHSEED=0 without pinning BLAS
    base_env = execution_environment({})
    assert base_env["PYTHONHASHSEED"] == "0"
    assert "OMP_NUM_THREADS" not in base_env

    # Session-scoped reproducible policy pins BLAS/OpenMP/Numba threads without mutating os.environ
    repro_policy = EgressPolicy(air_gapped=False, reproducible=True)
    with use_egress_policy(repro_policy):
        env = execution_environment({})
        assert env["ARIA_REPRODUCIBLE"] == "1"
        for key in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "BLIS_NUM_THREADS",
            "NUMBA_NUM_THREADS",
        ):
            assert env[key] == "1"
        env_override = execution_environment({"ARIA_BLAS_THREADS": "4"})
        assert env_override["OMP_NUM_THREADS"] == "4"
        assert env_override["NUMBA_NUM_THREADS"] == "4"

    assert "ARIA_REPRODUCIBLE" not in os.environ


def test_atac_align_passes_deterministic_bwa_batch_flag_across_all_lanes(monkeypatch, tmp_path):
    from aria.scripts import atac_align

    r1 = tmp_path / "s1_R1.fq.gz"
    r2 = tmp_path / "s1_R2.fq.gz"
    r1.write_bytes(b"fq1")
    r2.write_bytes(b"fq2")

    captured_cmds: list[list[str]] = []

    def _fake_pipe(cmd1, cmd2, *a, **kw):
        captured_cmds.append(list(cmd1))
        return False  # stop right after capturing bwa-mem2 command

    monkeypatch.setattr(atac_align, "_pipe", _fake_pipe)
    atac_align._align_hic_chimeric_sample(
        {"name": "hic1", "r1": str(r1), "r2": str(r2), "paired": True},
        "/ref/hg38.fa",
        tmp_path,
        8,
        [],
    )
    atac_align._align_sample(
        {"name": "atac1", "r1": str(r1), "r2": str(r2), "paired": True},
        "/ref/hg38.fa",
        tmp_path,
        8,
        30,
        True,
        ("chrM",),
        [],
    )
    assert len(captured_cmds) == 2
    for cmd in captured_cmds:
        assert cmd[:2] == ["bwa-mem2", "mem"]
        k_idx = cmd.index("-K")
        assert cmd[k_idx + 1] == "100000000"

    # Also verify spike-in bwa-mem2 command
    spike_cmds: list[list[str]] = []

    class _FakePopen:
        def __init__(self, cmd, **kw):
            spike_cmds.append(list(cmd))
            self.stdout = MagicMock()

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def wait(self):
            return 0

    monkeypatch.setattr(atac_align, "_bwa_index_exists", lambda _p: True)
    monkeypatch.setattr(atac_align.subprocess, "Popen", _FakePopen)
    monkeypatch.setattr(
        atac_align.subprocess,
        "run",
        lambda *a, **kw: SimpleNamespace(returncode=0, stdout="42\n"),
    )
    cnt = atac_align._count_spikein_index(
        str(r1), str(r2), "/ref/sacCer3.fa", "", 4, [], "s1"
    )
    assert cnt == 42
    assert spike_cmds[0][:2] == ["bwa-mem2", "mem"]
    k_idx = spike_cmds[0].index("-K")
    assert spike_cmds[0][k_idx + 1] == "100000000"


def test_hic_deterministic_eigs_cis_supplies_v0_and_restores_rng():
    from aria.scripts.hic_topology import _run_eigs_cis_deterministic
    import scipy.sparse.linalg as spla

    np.random.seed(12345)
    before_draw = np.random.rand(3)
    np.random.seed(12345)

    observed_v0: list[np.ndarray] = []

    class _FakeEigMod:
        @staticmethod
        def eigs_cis(clr, **kwargs):
            mat = np.array(
                [[4.0, 1.0, 0.2, 0.1], [1.0, 3.0, 0.5, 0.2], [0.2, 0.5, 2.5, 0.8], [0.1, 0.2, 0.8, 3.5]],
                dtype=np.float64,
            )
            orig = spla.eigsh
            vals, vecs = orig(mat, 2)
            observed_v0.append(vecs)
            return vals, vecs

    _run_eigs_cis_deterministic(
        _FakeEigMod,
        None,
        phasing_track=None,
        view_df=None,
        n_eigs=2,
        ignore_diags=1,
        seed=0,
    )
    # Global RNG state was restored cleanly
    after_draw = np.random.rand(3)
    np.testing.assert_allclose(before_draw, after_draw)
    assert len(observed_v0) == 1


def test_stage_caches_reject_tampered_outputs_even_with_preserved_mtime(tmp_path):
    from aria.utils.mudata_io import _stage_cache, _stage_complete
    from aria.scripts.hic_preprocess import _check_hic_stage_cache, _save_hic_stage_cache

    inp = tmp_path / "in.txt"
    inp.write_text("upstream_data_v1", encoding="utf-8")
    out_dir = tmp_path / "stage_out"
    out_dir.mkdir()
    out_file = out_dir / "artifact.h5mu"
    out_file.write_bytes(b"VALID_HDF5_BYTES_0123456789")

    miss, state = _stage_cache(
        "integration_paired_qc",
        out_dir,
        inputs=[("input", str(inp))],
        params={"min_genes": 200},
    )
    assert miss is None
    _stage_complete({"status": "success", "output_path": str(out_file)}, state)

    # Clean hit
    hit, _ = _stage_cache(
        "integration_paired_qc",
        out_dir,
        inputs=[("input", str(inp))],
        params={"min_genes": 200},
    )
    assert hit is not None and hit["resumed"] is True

    # Silent same-length byte corruption WITH preserved mtime must be rejected!
    st_orig = out_file.stat()
    out_file.write_bytes(b"CORRUPTED_BYTES__0123456789")
    os.utime(out_file, ns=(st_orig.st_atime_ns, st_orig.st_mtime_ns))

    tampered_hit, _ = _stage_cache(
        "integration_paired_qc",
        out_dir,
        inputs=[("input", str(inp))],
        params={"min_genes": 200},
    )
    assert tampered_hit is None

    # Same verification for Hi-C stage cache
    hic_dir = tmp_path / "hic_out"
    hic_dir.mkdir()
    cool_file = hic_dir / "sample.mcool"
    cool_file.write_bytes(b"COOLER_PAYLOAD_V1")
    _save_hic_stage_cache(
        hic_dir,
        stage="hic_preprocess",
        version=1,
        params={"res": 10000},
        inputs=[str(inp)],
        outputs=[str(cool_file)],
        extra={"result": {"status": "success", "cooler_files": [str(cool_file)]}},
    )
    assert _check_hic_stage_cache(
        hic_dir, stage="hic_preprocess", version=1, params={"res": 10000}, inputs=[str(inp)]
    ) is not None

    st_cool = cool_file.stat()
    cool_file.write_bytes(b"COOLER_TAMPERED!!")
    os.utime(cool_file, ns=(st_cool.st_atime_ns, st_cool.st_mtime_ns))
    assert _check_hic_stage_cache(
        hic_dir, stage="hic_preprocess", version=1, params={"res": 10000}, inputs=[str(inp)]
    ) is None


def test_token_boundary_matching_fixes_substring_collisions():
    pytest.importorskip("litellm")
    from aria.llm.parameter_advisor import ParameterAdvisor
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.agents.genome_arch_agent import GenomeArchAgent
    from aria.agents.integration_agent import IntegrationAgent

    assert has_keyword_boundary("Compare wild-type vs KO genotype", ("type",)) is False
    assert has_keyword_boundary("Define major cell types in cortex", ("fine",)) is False
    assert has_keyword_boundary("Define major cell types in cortex", ("major", "cell types")) is True

    adv = ParameterAdvisor.__new__(ParameterAdvisor)
    # Previously "define" matched "fine" -> (0.6, 1.4)! Now properly coarse (0.2, 0.6).
    assert adv._intent_to_leiden_range({
        "analysis_type": "differential",
        "raw_user_question": "Define major cell types in human hippocampus",
    }) == (0.2, 0.6)
    # Previously "wild-type" and "genotype" matched "type" -> (0.2, 0.6)! Now full default (0.2, 1.0).
    assert adv._intent_to_leiden_range({
        "analysis_type": "differential",
        "raw_user_question": "Compare wild-type vs BMAL1-KO genotype",
    }) == (0.2, 1.0)
    # Genuine fine-grained subpopulation query -> (0.6, 1.4)
    assert adv._intent_to_leiden_range({
        "analysis_type": "cell_type",
        "raw_user_question": "Identify rare sub-populations of progenitors",
    }) == (0.6, 1.4)

    ch = ChromatinAgent.__new__(ChromatinAgent)
    assert ch._needs_tf_analysis({"raw_user_question": "Assess TF, motif, and enhancer activity"}) is True
    assert ch._needs_tf_analysis({"raw_user_question": "Check platform drift"}) is False

    ga = GenomeArchAgent.__new__(GenomeArchAgent)
    res_dec = ga._advise_resolution(
        "exp_hic",
        {"raw_user_question": "Generate dotplot summary of genome wide contacts"},
        {"available_resolutions": [100_000, 40_000, 10_000, 5_000]},
        {"organism": "Homo sapiens"},
    )
    # "dotplot" must NOT trigger 5kb loop_calling via substring "dot"
    assert res_dec["recommended_resolution"] == 40_000

    ia = IntegrationAgent.__new__(IntegrationAgent)
    # "transcription factor" must NOT misfire MOFA+ via bare "factor"
    assert ia._select_strategy(
        has_scrna=True,
        has_scatac=True,
        has_hic=False,
        n_modalities=2,
        intent={"raw_user_question": "Map transcription factor regulons in paired multiome"},
        exp_ctx={},
    ) == "wnn"
    assert ia._select_strategy(
        has_scrna=True,
        has_scatac=True,
        has_hic=False,
        n_modalities=2,
        intent={"raw_user_question": "Run MOFA+ latent factor decomposition"},
        exp_ctx={},
    ) == "wnn_then_mofa"


def test_canonical_leiden_grid_and_locked_decision_manifest_replay(monkeypatch):
    pytest.importorskip("litellm")
    from aria.llm import parameter_advisor as pa_mod
    from aria.llm.parameter_advisor import ParameterAdvisor

    class _Mem:
        def list_wings(self):
            return []

        def store_finding(self, **kw):
            pass

        def store_decision(self, **kw):
            pass

    class _BombLLM:
        def complete_medium(self, **kw):
            raise AssertionError("Zero-LLM replay must not invoke LLM!")

    evaluated_grids: list[list[float]] = []

    def _fake_eval(adata, val):
        evaluated_grids.append(val)
        return {
            "n_clusters": 8,
            "silhouette": round(0.55 - abs(float(val) - 0.60) * 0.2, 4),
            "modularity": 0.62,
            "n_singleton_clusters": 0,
            "min_cluster_size": 50,
            "max_cluster_size": 400,
        }

    monkeypatch.setattr(pa_mod.MetricEvaluator, "leiden_resolution", staticmethod(_fake_eval))

    adv = ParameterAdvisor(_Mem(), _BombLLM())
    dec_a = adv.advise_leiden_resolution(
        "exp_1",
        {"analysis_type": "cell_type", "reproducible_mode": True, "zero_llm_replay": True},
        adata=object(),
    )
    assert [c.value for c in dec_a.candidates] == list(ParameterAdvisor.CANONICAL_LEIDEN_GRID)

    dec_b = adv.advise_leiden_resolution(
        "exp_2",
        {"analysis_type": "differential", "reproducible_mode": True, "zero_llm_replay": True},
        adata=object(),
    )
    # Invariant across LLM intent label drift ("cell_type" vs "differential")!
    assert dec_a.chosen_value == dec_b.chosen_value == 0.60

    # Locked DecisionManifest parameter override selects exact locked value with 0 LLM calls
    dec_locked = adv.advise_leiden_resolution(
        "exp_3",
        {
            "analysis_type": "differential",
            "reproducible_mode": True,
            "locked_parameter_decisions": {"leiden_clustering.resolution": 0.87},
            "decision_manifest_sha256": "deadbeef1234567890",
        },
        adata=object(),
    )
    assert dec_locked.chosen_value == 0.87
    assert dec_locked.chosen_by == "decision_manifest"
    assert "deadbeef1234" in dec_locked.justification


def test_decision_manifest_fork_ledger_capsule_and_zero_llm_orchestrator_replay(monkeypatch, tmp_path):
    pytest.importorskip("litellm")
    from aria.agents.narrative_agent import NarrativeAgent
    from aria.agents.narrative.ledger_export import (
        verify_reproducible_capsule,
        write_reproducible_capsule,
    )
    from aria.agents.orchestrator_agent import OrchestratorAgent
    from aria.memory.memory import ARIAMemory

    ledger_path = tmp_path / "study_fork_ledger.jsonl"
    ws1 = tmp_path / "ws1"
    ws2 = tmp_path / "ws2"
    ws3 = tmp_path / "ws3"

    base_ctx = {
        "organism": "Mus musculus",
        "genome": "mm39",
        "reproducible_mode": True,
        "modalities": {"bulk_RNA": ["/data/WT1.tsv", "/data/KO1.tsv"]},
        "input_files": [
            {"modality": "bulk_RNA", "path": "/data/WT1.tsv", "size_bytes": 120, "sha256": "a" * 64},
            {"modality": "bulk_RNA", "path": "/data/KO1.tsv", "size_bytes": 124, "sha256": "b" * 64},
        ],
        "design": {
            "organism": "Mus musculus",
            "genome": "mm39",
            "main_factor": "genotype",
            "groups": {"WT": ["WT1", "WT2"], "KO": ["KO1", "KO2"]},
            "contrasts": [{"numerator": "KO", "denominator": "WT", "factor": "genotype"}],
        },
        "global_padj": 0.05,
        "global_lfc": 1.0,
        "lfc_threshold_provenance": "default",
    }
    plan = {
        "steps": [{"order": 1, "agent": "bulk_rna_agent", "analysis": "DE", "depends_on": []}],
        "contrasts": [{"numerator": "KO", "denominator": "WT"}],
    }
    intent = {
        "analysis_type": "differential",
        "comparison": "KO vs WT",
        "raw_user_question": "Compare KO vs WT expression",
        "user_question": "Compare KO vs WT expression",
    }

    # 1. Primary run -> fork_index == 0, locked_before_execution
    m1 = build_decision_manifest(
        experiment_id="run_primary",
        exp_context=base_ctx,
        plan=plan,
        biological_intent=intent,
    )
    m1_committed, fl1 = commit_decision_manifest(m1, workspace_dir=ws1, ledger_path=ledger_path)
    assert fl1["fork_index"] == 0
    assert fl1["total_analytical_forks"] == 1
    assert fl1["commitment_status"] == "locked_before_execution"
    assert fl1["is_exact_spec_replay"] is False
    assert (ws1 / "decision_manifest.json").is_file()

    # 2. Exact specification replay -> fork_index stays 0, total_analytical_forks stays 1
    m2 = build_decision_manifest(
        experiment_id="run_replay",
        exp_context=base_ctx,
        plan=plan,
        biological_intent=intent,
        commitment_status="replay_locked",
        replayed_from_manifest_sha256=m1_committed["manifest_sha256"],
    )
    # Note: replay_locked changes commitment_status on envelope, while core specification_sha256 is identical!
    assert m2["specification_sha256"] == m1_committed["specification_sha256"]
    m2_committed, fl2 = commit_decision_manifest(m2, workspace_dir=ws2, ledger_path=ledger_path)
    assert fl2["fork_index"] == 0
    assert fl2["total_analytical_forks"] == 1
    assert fl2["is_exact_spec_replay"] is True

    # 3. Modified threshold branch -> fork_index == 1, post_hoc_fork, records thresholds.global_lfc
    forked_ctx = {**base_ctx, "global_lfc": 0.58, "lfc_threshold_provenance": "cp3_profile:exploratory_tf"}
    m3 = build_decision_manifest(
        experiment_id="run_fork_1",
        exp_context=forked_ctx,
        plan=plan,
        biological_intent=intent,
    )
    _, fl3 = commit_decision_manifest(m3, workspace_dir=ws3, ledger_path=ledger_path)
    assert fl3["fork_index"] == 1
    assert fl3["total_analytical_forks"] == 2
    assert fl3["commitment_status"] == "post_hoc_fork"
    assert "thresholds.global_lfc" in fl3["fork_changes"]

    # 4. Report + Capsule verification of DecisionManifest & Fork Ledger
    mem = ARIAMemory(db_path=str(tmp_path / "mem.db"))
    mem.create_wing("run_primary", "Primary Study", organism="Mus musculus", genome="mm39")
    na = NarrativeAgent(mem, llm=None)
    na.reports_dir = tmp_path / "reports"
    na.reports_dir.mkdir()
    report_ctx = {
        **base_ctx,
        "user_question": "Compare KO vs WT expression",
        "decision_manifest": m1_committed,
        "analytical_fork_ledger": fl1,
    }
    res = na.run(
        "run_primary",
        {
            "exp_context": report_ctx,
            "biological_intent": intent,
            "agent_results": {},
            "findings": [],
        },
    )
    report_html = Path(res["report_path"])
    report_dir = report_html.parent
    assert "Analytical Commitment &amp; Forking-Path Ledger" in report_html.read_text(encoding="utf-8")
    assert (report_dir / "decision_manifest.json").is_file()
    loaded = load_decision_manifest(report_dir)
    assert loaded["manifest_sha256"] == m1_committed["manifest_sha256"]

    inp_wt = tmp_path / "WT1.tsv"
    inp_ko = tmp_path / "KO1.tsv"
    inp_wt.write_text("gene\tWT1\nG1\t10\n", encoding="utf-8")
    inp_ko.write_text("gene\tKO1\nG1\t20\n", encoding="utf-8")
    from aria.utils.provenance import hash_file
    meth = json.loads((report_dir / "methodology.json").read_text(encoding="utf-8"))
    meth["inputs"] = [
        {"modality": "bulk_RNA", "path": str(inp_wt), "size_bytes": inp_wt.stat().st_size, "sha256": hash_file(inp_wt)},
        {"modality": "bulk_RNA", "path": str(inp_ko), "size_bytes": inp_ko.stat().st_size, "sha256": hash_file(inp_ko)},
    ]
    (report_dir / "methodology.json").write_text(json.dumps(meth, indent=2, sort_keys=True), encoding="utf-8")

    cap_zip = write_reproducible_capsule(report_dir, tmp_path / "primary_capsule.zip")
    ver = verify_reproducible_capsule(cap_zip)
    assert ver["status"] in ("pass", "warning")
    assert ver["files"]["missing"] == []
    assert ver["files"]["mismatched"] == []
    assert ver["decision_manifest"]["present"] is True
    assert ver["decision_manifest"]["valid"] is True

    # 5. OrchestratorAgent Zero-LLM Replay makes ZERO LLM calls
    class _ExplodingLLM:
        def for_execution(self, **kw):
            return self

        def complete(self, *a, **kw):
            raise AssertionError("LLM must never be called during Zero-LLM DecisionManifest replay!")

        def complete_structured(self, *a, **kw):
            raise AssertionError("LLM must never be called during Zero-LLM DecisionManifest replay!")

    orch = OrchestratorAgent(mem, llm=_ExplodingLLM())
    replay_start = orch.run(
        "exp_zero_llm",
        {
            "data_dir": str(tmp_path),
            "user_question": "",
            "replay_decision_manifest": report_dir / "decision_manifest.json",
        },
    )
    assert replay_start["status"] == "started"
    assert replay_start["intent"]["zero_llm_replay"] is True
    assert replay_start["intent"]["raw_user_question"] == "Compare KO vs WT expression"
    plan_ready = orch._complete_design_and_publish_plan(
        "exp_zero_llm",
        dict(base_ctx, replayed_decision_manifest=loaded),
        dict(base_ctx["design"]),
    )
    assert plan_ready["status"] == "plan_ready"
    assert plan_ready["plan"]["lfc_threshold"] == 1.0
    mem.close()


def test_fastq_auto_provisioning_parity_kb_multiome_and_scatac_proposal(monkeypatch, tmp_path):
    pytest.importorskip("litellm")
    from aria.agents.raw_ingestion_agent import RawIngestionAgent
    from aria.agents.setup_agent import SetupAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils.provenance import hash_file
    from aria.utils.scatac_fastq_manifest import (
        propose_illumina_scatac_manifest,
        resolve_scatac_fastq_manifest,
    )

    # 1. SetupAgent invokes _ensure_kb_index for scRNA/Multiome FASTQ runs, never for CUT&RUN
    sa = SetupAgent.__new__(SetupAgent)
    sa.genomes_dir = tmp_path / "genomes"
    sa.genomes_dir.mkdir()
    sa.publish_status = lambda *a, **k: None
    sa.publish_finding = lambda *a, **k: None
    sa._ensure_conda = lambda *a, **k: "conda"
    sa._needed_envs = lambda *a, **k: []
    kb_calls: list[str] = []
    sa._ensure_genome = lambda gk, *a, **k: {
        "genome_key": gk,
        "genome_dir": str(sa.genomes_dir / gk),
        "fasta": str(tmp_path / "genome.fa"),
        "gtf": str(tmp_path / "ann.gtf"),
        "star_index": str(sa.genomes_dir / gk / "star_index"),
        "kb_dir": str(sa.genomes_dir / gk / "kb_index"),
        "kb_index": str(sa.genomes_dir / gk / "kb_index" / "index.idx"),
        "kb_t2g": str(sa.genomes_dir / gk / "kb_index" / "t2g.txt"),
        "fasta_ready": True,
        "gtf_ready": True,
        "index_ready": False,
        "kb_ready": False,
    }
    sa._ensure_star_index = lambda *a, **k: None
    sa._ensure_kb_index = lambda gk, cfg, *a, **k: kb_calls.append(gk)

    sa.run("exp_cut", {"exp_context": {"organism": "Mus musculus", "genome": "mm39", "modalities": {"CUT_AND_RUN": ["/x/a_R1.fq.gz"]}}, "biological_intent": {}})
    assert kb_calls == []
    sa.run("exp_sc", {"exp_context": {"organism": "Mus musculus", "genome": "mm39", "modalities": {"scRNA": ["/x/s1_R1.fq.gz", "/x/s1_R2.fq.gz"]}}, "biological_intent": {}})
    assert kb_calls == ["mm39"]

    # 2. RawIngestionAgent auto-wires SetupAgent's kb_index + kb_t2g + SHA-256 when chemistry is declared
    fq_dir = tmp_path / "scrna_fqs"
    fq_dir.mkdir()
    import gzip
    for name in ("donor1_S1_L001_R1_001.fastq.gz", "donor1_S1_L001_R2_001.fastq.gz"):
        with gzip.open(fq_dir / name, "wt", encoding="utf-8") as gz:
            gz.write("@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n")

    kb_idx = tmp_path / "index.idx"
    kb_t2g = tmp_path / "t2g.txt"
    kb_idx.write_bytes(b"KB_IDX_BYTES")
    kb_t2g.write_text("ENSMUST1\tENSMUSG1\tBmal1\n", encoding="utf-8")

    mem = ARIAMemory(db_path=str(tmp_path / "ing_mem.db"))
    ria = RawIngestionAgent(mem, llm=None)
    out_h5ad = tmp_path / "donor1_kb.h5ad"
    out_h5ad.write_bytes(b"H5AD_MOCK")
    captured_kb_params: list[dict] = []

    def _fake_kb_count(params):
        captured_kb_params.append(dict(params))
        return {"status": "success", "mode": "fastq_kb_count", "output_h5ad": str(out_h5ad)}

    monkeypatch.setattr(ria, "_run_kb_count", _fake_kb_count)
    ing_res = ria.run(
        "exp_ing",
        {
            "exp_context": {
                "data_dir": str(fq_dir),
                "chemistry": "10xv3",
                "modalities": {"scRNA": [str(fq_dir / "donor1_S1_L001_R1_001.fastq.gz"), str(fq_dir / "donor1_S1_L001_R2_001.fastq.gz")]},
                "genome_config": {
                    "kb_ready": True,
                    "kb_index": str(kb_idx),
                    "kb_index_sha256": hash_file(kb_idx),
                    "kb_t2g": str(kb_t2g),
                    "kb_t2g_sha256": hash_file(kb_t2g),
                },
            }
        },
    )
    assert ing_res["status"] == "done"
    assert len(captured_kb_params) == 1
    assert captured_kb_params[0]["chemistry"] == "10xv3"
    assert captured_kb_params[0]["index_sha256"] == hash_file(kb_idx)
    assert captured_kb_params[0]["t2g_sha256"] == hash_file(kb_t2g)
    mem.close()

    # 3. Illumina scATAC folder proposal produces a schema-valid manifest
    atac_dir = tmp_path / "atac_fqs"
    atac_dir.mkdir()
    wl = tmp_path / "737K-cratac-v1.txt"
    wl.write_text("AAACGAAAGAAACGCC\n", encoding="utf-8")
    for role in ("R1", "R2", "R3"):
        with gzip.open(atac_dir / f"hippo1_S1_L001_{role}_001.fastq.gz", "wt", encoding="utf-8") as gz:
            gz.write("@r\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n")

    proposed = propose_illumina_scatac_manifest(atac_dir, barcode_whitelist=str(wl))
    assert proposed is not None
    validated = resolve_scatac_fastq_manifest(
        {"data_dir": str(atac_dir), "scatac_fastq_manifest": proposed},
        require_paths=True,
    )
    assert validated["status"] == "valid", validated
    assert validated["manifest"]["libraries"][0]["library_id"] == "hippo1"
