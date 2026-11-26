"""Regression tests for Fix Round F4 slice D9 (hygiene, warnings, contracts, config, B2 determinism).

Built directly from D9_2026-11-26_verified.json repro_tests (D9-001..D9-009)
plus PYTHONHASHSEED byte-reproducibility of governance_b2.score_b2().
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def test_d9_001_pseudobulk_de_surfaces_real_pydeseq2_warnings(tmp_path: Path):
    ad = pytest.importorskip("anndata")
    pytest.importorskip("pydeseq2")
    from aria.scripts.rna_bulk.deseq2 import _run_deseq2
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de

    n_genes, n_cells_per_donor = 20, 12
    donors = ["d1", "d2", "d3", "d4", "d5", "d6"]
    conds = ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"]
    rows, obs_rows, pb_cols = [], [], {}
    for donor, cond in zip(donors, conds, strict=True):
        mat = np.zeros((n_cells_per_donor, n_genes), dtype=np.int32)
        mat[:, :10] = 5
        mat[:, 10:] = 200 if cond == "treat" else 0
        mat[0, 0] += 1 if cond == "treat" else 0
        rows.append(mat)
        pb_cols[donor] = mat.sum(axis=0)
        for _ in range(n_cells_per_donor):
            obs_rows.append({"cell_type": "Tcell", "condition": cond, "donor": donor})
    genes = [f"GENE{i}" for i in range(n_genes)]
    X = np.vstack(rows)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(obs_rows, index=[f"c{i}" for i in range(X.shape[0])]),
        var=pd.DataFrame(index=genes),
    )
    h5ad_path = tmp_path / "pb_degenerate.h5ad"
    adata.write_h5ad(h5ad_path)
    pb_df = pd.DataFrame(pb_cols, index=genes)
    pb_meta = pd.DataFrame({"condition": conds}, index=donors)
    _, bulk_warns = _run_deseq2(
        pb_df, pb_meta, "condition", "treat", "ctrl", 0.05, 0.5, lfc_shrink=True
    )
    pydeseq_warns = [w for w in bulk_warns if "pydeseq2 fit warning" in w]
    assert len(pydeseq_warns) > 0

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res = rna_pseudobulk_de(
            {
                "data_path": str(h5ad_path),
                "groupby": "cell_type",
                "condition_col": "condition",
                "replicate_col": "donor",
                "comparisons": [["treat", "ctrl"]],
                "use_raw": False,
                "output_dir": str(tmp_path / "out"),
            }
        )
        filters_after = list(warnings.filters)

    assert res["status"] == "success"
    comp = res["per_group"]["Tcell"]["per_comparison"]["treat_vs_ctrl"]
    blanket_ignore = any(
        f[0] == "ignore" and f[1] is None and f[2] is Warning for f in filters_after
    )
    surfaced = (
        list(res.get("warnings") or [])
        + list(comp.get("warnings") or [])
        + list(comp.get("fit_warnings") or [])
    )
    assert not blanket_ignore and (len(surfaced) > 0 or len(caught) > 0), (
        f"Control _run_deseq2 surfaced {len(pydeseq_warns)} warning(s), "
        f"but rna_pseudobulk_de installed blanket ignore={blanket_ignore}, "
        f"caught={len(caught)}, surfaced={surfaced!r}"
    )


def test_d9_002_base_run_script_does_not_silence_log2_runtime_warnings(
    tmp_path: Path, monkeypatch
):
    from aria.scripts._base import run_script

    in_json, out_json = tmp_path / "in.json", tmp_path / "out.json"
    in_json.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["dummy_script.py", str(in_json), str(out_json)])

    def analysis_fn(_params: dict) -> dict:
        vals = np.log2(np.array([-2.0, 0.0, 4.0], dtype=float))
        return {"status": "success", "n_nan": int(np.isnan(vals).sum())}

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        run_script(analysis_fn)
    log2_warns = [
        w
        for w in caught
        if issubclass(w.category, RuntimeWarning) and "log2" in str(w.message)
    ]
    assert len(log2_warns) >= 2, (
        f"_base.run_script silently suppressed log2 RuntimeWarnings: {log2_warns}"
    )


def test_d9_003_rna_advise_resolution_has_script_contract():
    pytest.importorskip("pydantic")
    from aria.utils.script_contracts import contract_for_script

    contract = contract_for_script("aria/scripts/rna_advise_resolution.py")
    assert contract is not None


def test_d9_004_aria_config_file_honored_by_llm_provider(tmp_path: Path, monkeypatch):
    pytest.importorskip("yaml")
    pytest.importorskip("litellm")
    from aria import doctor
    from aria.llm.provider import LLMProvider, TaskTier

    cfg_path = tmp_path / "custom_config.yaml"
    cfg_path.write_text(
        "llm:\n  heavy:\n    provider: ollama\n    model: local-audited-70b\n"
        "    api_base: http://127.0.0.1:11434\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("ARIA_CONFIG_FILE", str(cfg_path))
    monkeypatch.setenv("HOME", str(tmp_path / "empty_home"))
    assert doctor._read_llm_config_models().get("heavy") == [
        "ollama/local-audited-70b (local)"
    ]
    prov = LLMProvider.from_config()
    primary_heavy = prov.models[TaskTier.HEAVY][0]
    assert primary_heavy.model == "ollama/local-audited-70b"


def test_d9_005_gtf_env_vars_consistent_across_modalities(tmp_path: Path, monkeypatch):
    from aria.scripts.chromatin_peak_annotation import _resolve_gtf as chrom_resolve_gtf
    from aria.scripts.rna_bulk.gtf_io import _locate_gtf as bulk_locate_gtf
    from aria.utils.peak2gene_core import resolve_gtf_path as p2g_resolve_gtf

    gtf = tmp_path / "staged.gtf"
    gtf.write_text(
        'chr1\tARIA\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "TP53";\n',
        encoding="utf-8",
    )
    counts = tmp_path / "study" / "counts.tsv"
    counts.parent.mkdir(parents=True)
    counts.write_text("gene\ts1\nG1\t10\n", encoding="utf-8")
    monkeypatch.setenv("HOME", str(tmp_path / "empty_home"))
    monkeypatch.setenv("ARIA_GTF", str(gtf))
    monkeypatch.delenv("ARIA_GTF_PATH", raising=False)
    monkeypatch.delenv("ARIA_GENOME_DIR", raising=False)

    assert chrom_resolve_gtf({"genome": "hg38"}, tmp_path) == gtf
    assert p2g_resolve_gtf(explicit_gtf=None, genome="hg38") == str(gtf)
    assert bulk_locate_gtf(counts, hint=None, genome="hg38", organism="Homo sapiens") == gtf


def test_d9_006_data_audit_proposes_illumina_scatac_manifest_in_cp1(tmp_path: Path):
    pytest.importorskip("litellm")
    from aria.agents.data_audit_agent import DataAuditAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils.scatac_fastq_manifest import propose_illumina_scatac_manifest

    fq_dir = tmp_path / "scatac_10x"
    fq_dir.mkdir()
    for role in ("R1", "R2", "R3"):
        with gzip.open(fq_dir / f"donor1_S1_L001_{role}_001.fastq.gz", "wt") as gz:
            gz.write("@r1\n" + ("ACGT" * 12) + "\n+\n" + ("I" * 48) + "\n")
    proposed = propose_illumina_scatac_manifest(fq_dir)
    assert proposed is not None and len(proposed["libraries"]) == 1

    mem = ARIAMemory(db_path=str(tmp_path / "mem.db"))
    res = DataAuditAgent(memory=mem, llm=None).run(
        "exp_d9_scatac",
        {"data_dir": str(fq_dir), "user_question": "Analyze scATAC-seq accessibility"},
    )
    exp_ctx = res.get("exp_context") or {}
    assert exp_ctx.get("scatac_fastq_manifest") or exp_ctx.get(
        "proposed_scatac_fastq_manifest"
    )


def test_d9_007_undocumented_env_does_not_bypass_bulk_metadata_gate(monkeypatch):
    pytest.importorskip("litellm")
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.scripts.rna_bulk.counts import (
        _load_or_infer_metadata,
        _metadata_inference_allowed,
    )

    monkeypatch.delenv("ARIA_ALLOW_FILENAME_FALLBACK", raising=False)
    monkeypatch.setenv("ARIA_ALLOW_BULK_METADATA_INFERENCE", "1")
    assert not BulkRNAAgent._filename_fallback_allowed()
    counts = pd.DataFrame(
        {
            "ctrl_1": [100, 200],
            "ctrl_2": [110, 190],
            "treat_1": [500, 50],
            "treat_2": [520, 40],
        },
        index=["G1", "G2"],
    )
    meta, _warns = _load_or_infer_metadata(
        counts, "", "condition", allow_inference=_metadata_inference_allowed({})
    )
    assert meta is None


def test_d9_008_cross_modal_memory_tf_surfaces_ambiguous_direction():
    from aria.agents.narrative.synthesis.cross_modal_synthesis import _memory_tf

    tf_dirs = {
        "SOX10": {"tf": "SOX10", "directions": {1, -1}, "paralog": False},
        "OLIG2": {"tf": "OLIG2", "directions": {1}, "paralog": False},
    }
    rna_genes = {"SOX10": ("SOX10", 2.0), "OLIG2": ("OLIG2", 1.5)}
    res = _memory_tf(tf_dirs, rna_genes, flip=1)
    assert res is not None
    assert res.get("n_ambiguous_direction") == 1


def test_d9_009_llm_provider_from_config_preserves_extra_params(tmp_path: Path):
    pytest.importorskip("yaml")
    pytest.importorskip("litellm")
    from aria.llm.provider import LLMProvider, TaskTier

    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        "llm:\n  heavy:\n    provider: openai\n    model: o3-mini\n"
        "    extra_params:\n      reasoning_effort: low\n",
        encoding="utf-8",
    )
    prov = LLMProvider.from_config(str(cfg_path))
    assert prov.models[TaskTier.HEAVY][0].extra_params == {"reasoning_effort": "low"}


def test_b2_score_json_byte_identical_across_pythonhashseed():
    pytest.importorskip("pydantic")
    repo_root = str(Path(__file__).resolve().parents[1])
    code = (
        "import json; "
        "from aria.benchmarks.governance_b2 import score_b2; "
        "print(json.dumps(score_b2(), sort_keys=True, default=str))"
    )
    outputs = []
    for seed in ("0", "1"):
        env = dict(os.environ, PYTHONHASHSEED=seed, PYTHONPATH=repo_root)
        proc = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            check=True,
            env=env,
        )
        outputs.append(proc.stdout)
    assert outputs[0] == outputs[1], "score_b2() JSON differs between PYTHONHASHSEED=0 and 1"
