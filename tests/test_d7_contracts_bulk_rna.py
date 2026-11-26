"""Regression tests for D7 layer-contract fixes in bulk RNA-seq."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


def test_d7_007_rna_bulk_de_emits_fitted_formula_in_design_used_worker_env(tmp_path):
    """D7-007 (worker-env pure script check): rna_bulk_de.py must emit the actual
    fitted design formula (including modeled covariates) in design_used rather than
    just the main factor name."""
    pytest.importorskip("pydeseq2")
    from aria.scripts.rna_bulk_de import bulk_rna_de

    rng = np.random.default_rng(13)
    samples = ["c1", "c2", "c3", "c4", "t1", "t2", "t3", "t4"]
    counts = rng.poisson(150, size=(40, 8))
    counts[:6, 4:] = rng.poisson(900, size=(6, 4))
    counts_tsv = tmp_path / "counts_worker.tsv"
    pd.DataFrame(counts, index=[f"G{i}" for i in range(40)], columns=samples).to_csv(
        counts_tsv, sep="\t", index_label="gene_id"
    )
    meta_tsv = tmp_path / "meta_worker.tsv"
    pd.DataFrame({
        "sample": samples,
        "condition": ["ctrl"] * 4 + ["treat"] * 4,
        "batch": ["B1", "B1", "B2", "B2", "B1", "B1", "B2", "B2"],
    }).to_csv(meta_tsv, sep="\t", index=False)

    res = bulk_rna_de({
        "counts_path": str(counts_tsv),
        "metadata_path": str(meta_tsv),
        "design_factor": "condition",
        "covariates": ["batch"],
        "contrasts": [{"name": "treat_vs_ctrl", "numerator": "treat", "denominator": "ctrl"}],
        "output_dir": str(tmp_path / "de_worker_out"),
        "padj_threshold": 0.05,
        "lfc_threshold": 0.5,
        "min_counts": 5,
        "lfc_shrinkage": False,
        "run_ora": False,
        "run_plots": False,
    })
    assert res["status"] == "success", res
    assert res["design_used"] == "~ batch + condition", res.get("design_used")


def test_d7_001_bulk_rna_apply_design_must_not_fallback_positionally(tmp_path):
    """D7-001: When design sample IDs (GSM*) do not match count-matrix column names
    (SRR*) and no explicit GSM<->SRR map is supplied, BulkRNAAgent._apply_design
    must fail closed naming unmatched IDs on both sides rather than mapping by
    column position and inverting DESeq2 log2FC signs."""
    pytest.importorskip("litellm")
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.memory.memory import ARIAMemory

    rng = np.random.default_rng(7)
    n_genes = 60
    ctrl_counts = rng.poisson(200, size=(n_genes, 3))
    treat_counts = rng.poisson(200, size=(n_genes, 3))
    treat_counts[:10, :] = rng.poisson(1200, size=(10, 3))
    cols = ["SRR9001", "SRR9002", "SRR9003", "SRR9004", "SRR9005", "SRR9006"]
    df = pd.DataFrame(
        np.hstack([treat_counts, ctrl_counts]),
        index=[f"GENE{i:03d}" for i in range(n_genes)],
        columns=cols,
    )
    tsv = tmp_path / "counts.tsv"
    df.to_csv(tsv, sep="\t", index_label="gene_id")

    design = {
        "groups": {
            "ctrl": ["GSM1001", "GSM1002", "GSM1003"],
            "treat": ["GSM1004", "GSM1005", "GSM1006"],
        },
        "main_factor": "condition",
        "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
    }
    agent = BulkRNAAgent(ARIAMemory(":memory:"), llm=None)
    with pytest.raises(ValueError) as exc_info:
        agent._apply_design(design, [str(tsv)], "exp_pos")
    msg = str(exc_info.value)
    assert "GSM1001" in msg and "GSM1006" in msg
    assert "SRR9001" in msg and "SRR9006" in msg


def test_d7_001_bulk_rna_apply_design_maps_explicit_gsm_to_srr_and_geo_soft(tmp_path):
    """D7-001 positive control: explicit GSM<->SRR mapping (via design['gsm_to_srr']
    or GEO SOFT !Sample_relation per-sample SRA accessions) resolves reordered
    SRR matrix columns without positional guessing."""
    pytest.importorskip("litellm")
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.connectors.geo_connector import _infer_design, _parse_soft_text
    from aria.memory.memory import ARIAMemory

    # Matrix columns intentionally ordered [treat SRRs, ctrl SRRs] (reversed vs design).
    cols = ["SRR9004", "SRR9005", "SRR9006", "SRR9001", "SRR9002", "SRR9003"]
    tsv = tmp_path / "counts.tsv"
    pd.DataFrame(
        np.ones((5, len(cols)), dtype=int) * 50,
        index=[f"GENE{i}" for i in range(5)],
        columns=cols,
    ).to_csv(tsv, sep="\t", index_label="gene_id")

    agent = BulkRNAAgent(ARIAMemory(":memory:"), llm=None)

    # 1. Explicit gsm_to_srr map on design resolves exact labels despite reversed columns.
    design_explicit = {
        "groups": {
            "ctrl": ["GSM1001", "GSM1002", "GSM1003"],
            "treat": ["GSM1004", "GSM1005", "GSM1006"],
        },
        "main_factor": "condition",
        "gsm_to_srr": {
            "GSM1001": "SRR9001",
            "GSM1002": "SRR9002",
            "GSM1003": "SRR9003",
            "GSM1004": "SRR9004",
            "GSM1005": "SRR9005",
            "GSM1006": "SRR9006",
        },
        "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
    }
    names, labels, factor, contrasts = agent._apply_design(design_explicit, [str(tsv)], "exp_map")
    assert names == cols
    assert labels == {
        "SRR9004": "treat",
        "SRR9005": "treat",
        "SRR9006": "treat",
        "SRR9001": "ctrl",
        "SRR9002": "ctrl",
        "SRR9003": "ctrl",
    }
    assert factor == "condition"
    assert contrasts and contrasts[0]["numerator"] == "treat"

    # 2. GEO SOFT parser preserves per-sample SRA relations into sample_aliases.
    soft = "\n".join([
        "!Series_title = Demo GEO",
        "^SAMPLE = GSM1001",
        "!Sample_title = ctrl_a",
        "!Sample_characteristics_ch1 = condition: ctrl",
        "!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRR9001",
        "^SAMPLE = GSM1002",
        "!Sample_title = ctrl_b",
        "!Sample_characteristics_ch1 = condition: ctrl",
        "!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRR9002",
        "^SAMPLE = GSM1003",
        "!Sample_title = ctrl_c",
        "!Sample_characteristics_ch1 = condition: ctrl",
        "!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRR9003",
        "^SAMPLE = GSM1004",
        "!Sample_title = treat_a",
        "!Sample_characteristics_ch1 = condition: treat",
        "!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRR9004",
        "^SAMPLE = GSM1005",
        "!Sample_title = treat_b",
        "!Sample_characteristics_ch1 = condition: treat",
        "!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRR9005",
        "^SAMPLE = GSM1006",
        "!Sample_title = treat_c",
        "!Sample_characteristics_ch1 = condition: treat",
        "!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRR9006",
    ])
    meta = _parse_soft_text(soft)
    inferred = _infer_design(meta)
    assert "SRR9001" in inferred["sample_aliases"]["GSM1001"]
    _, geo_labels, _, _ = agent._apply_design(inferred, [str(tsv)], "exp_geo_sra")
    assert geo_labels == labels


def test_d7_007_bulk_rna_agent_records_actual_fitted_design_decision(tmp_path):
    """D7-007: rna_bulk_de must emit the actual fitted design formula (including
    covariates) in top-level design_used, and BulkRNAAgent._record_methodology_decisions
    must record that fitted formula and covariate rationale into ARIAMemory."""
    pytest.importorskip("litellm")
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.memory.memory import ARIAMemory
    from aria.utils.environment_manager import env_manager

    rng = np.random.default_rng(13)
    samples = ["c1_b1", "c2_b1", "c3_b2", "c4_b2", "t1_b1", "t2_b1", "t3_b2", "t4_b2"]
    n_genes = 50
    mat = rng.poisson(200, size=(n_genes, len(samples)))
    mat[:5, 4:] = rng.poisson(800, size=(5, 4))
    counts_tsv = tmp_path / "counts.tsv"
    pd.DataFrame(mat, index=[f"G{i}" for i in range(n_genes)], columns=samples).to_csv(
        counts_tsv, sep="\t", index_label="gene_id"
    )
    meta_tsv = tmp_path / "meta.tsv"
    pd.DataFrame({
        "sample": samples,
        "condition": ["ctrl"] * 4 + ["treat"] * 4,
        "batch": ["b1", "b1", "b2", "b2", "b1", "b1", "b2", "b2"],
    }).to_csv(meta_tsv, sep="\t", index=False)

    de_res = env_manager.run_in_stack(
        "rna",
        "aria/scripts/rna_bulk_de.py",
        {
            "files": [str(counts_tsv)],
            "metadata_file": str(meta_tsv),
            "design_factor": "condition",
            "covariates": ["batch"],
            "contrasts": [{"name": "treat_vs_ctrl", "numerator": "treat", "denominator": "ctrl"}],
            "padj_threshold": 0.05,
            "lfc_threshold": 0.5,
            "output_dir": str(tmp_path / "de_out"),
            "run_ora": False,
            "run_gsea": False,
            "run_qc_plots": False,
        },
    )
    assert de_res["status"] == "success", de_res
    assert "batch" in de_res["contrasts"][0]["fitted_design_formula"]
    assert "batch" in de_res["design_used"]

    mem = ARIAMemory(":memory:")
    agent = BulkRNAAgent(mem, llm=None)
    agent._record_methodology_decisions("exp_bulk_1", de_res)
    decisions = mem.get_decisions("exp_bulk_1")
    design_dec = next(d for d in decisions if d["question"] == "DESeq2 design formula")
    assert "batch" in design_dec["decision"], (
        f"Decisions table recorded {design_dec['decision']!r} with rationale "
        f"{design_dec['rationale']!r} despite DESeq2 fitting "
        f"{de_res['contrasts'][0]['fitted_design_formula']!r}"
    )
    assert "No batch or covariate adjustment applied" not in design_dec["rationale"]
    assert "batch" in design_dec["rationale"].lower()


def test_d7_012_bulk_rna_adapter_surfaces_unmodeled_or_confounded_batch_from_design_agent():
    """D7-012: bulk_rna_evidence must derive the BATCH confound caveat from real
    DesignAgent fields (confounded_covariates, inferred_design.covariates vs
    resolve_design_contract covariates) rather than reading nonexistent design['has_batch']."""
    pytest.importorskip("litellm")
    from aria.agents.design_agent import DesignAgent
    from aria.agents.narrative.hypothesis.adapters.bulk_rna import bulk_rna_evidence
    from aria.agents.narrative.hypothesis.caveats import BATCH
    from aria.memory.memory import ARIAMemory

    da = DesignAgent(ARIAMemory(":memory:"), llm=None)
    da._organism = "Homo sapiens"
    da._genome = "hg38"
    da._confirmed_groups = {"ctrl": ["c1", "c2", "c3"], "treat": ["t1", "t2", "t3"]}
    da._main_factor = "condition"
    da._batch_covariate = None  # unmodeled
    da._inferred_design = {
        "source": "sample_sheet",
        "covariates": ["batch"],
        "confounded_covariates": ["batch"],
    }
    design = da._build_design()
    assert design["confounded_covariates"] == ["batch"]
    assert "has_batch" not in design

    agent_results = {
        "bulk_rna_agent": {
            "status": "done",
            "findings": {
                "contrasts": [
                    {
                        "name": "treat_vs_ctrl",
                        "status": "success",
                        "top_genes": [{"gene": "TP53", "log2fc": 2.4, "padj": 1e-6}],
                    }
                ]
            },
        }
    }
    sigs = bulk_rna_evidence(
        agent_results,
        run_ledger=None,
        exp_ctx={"design": design, "inferred_design": da._inferred_design},
    )
    assert sigs and BATCH in sigs[0].caveats_inherited

    # Positive control: when batch is modeled and NOT confounded, BATCH caveat is absent.
    da_clean = DesignAgent(ARIAMemory(":memory:"), llm=None)
    da_clean._organism = "Homo sapiens"
    da_clean._genome = "hg38"
    da_clean._confirmed_groups = {"ctrl": ["c1", "c2", "c3"], "treat": ["t1", "t2", "t3"]}
    da_clean._main_factor = "condition"
    da_clean._batch_covariate = "batch"
    da_clean._inferred_design = {
        "source": "sample_sheet",
        "covariates": ["batch"],
        "confounded_covariates": [],
    }
    design_clean = da_clean._build_design()
    sigs_clean = bulk_rna_evidence(
        agent_results,
        run_ledger=None,
        exp_ctx={"design": design_clean, "inferred_design": da_clean._inferred_design},
    )
    assert sigs_clean and BATCH not in sigs_clean[0].caveats_inherited
