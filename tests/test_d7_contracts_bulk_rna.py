"""Regression tests for D7 layer-contract fixes in bulk RNA-seq."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


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
