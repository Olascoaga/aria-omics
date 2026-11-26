"""Regression tests for D4 silent data misalignment fixes in cross-modal concordance."""
from __future__ import annotations

from pathlib import Path

import pandas as pd


def test_d4_011_load_rna_table_does_not_splice_or_fabricate_duplicate_symbols():
    """D4-011: load_rna_table must never splice columns across distinct loci via
    groupby().first() or OR rna_known_sig onto an opposing/non-significant locus."""
    from aria.utils.cross_modal_concordance import load_rna_table

    # 1. In-memory up_gene_ids ('TBCE') + opposite-sign non-sig locus in top_genes
    # must NOT fabricate a chimeric significant downregulated TBCE.
    warns: list[str] = []
    comp = {
        "up_gene_ids": ["TBCE"],
        "top_genes": [
            {"gene_id": "ENSG00000284770", "gene_symbol": "TBCE", "log2FoldChange": -2.5, "padj": 0.85},
            {"gene_id": "ENSG00000141510", "gene_symbol": "TP53", "log2FoldChange": 1.8, "padj": 0.001},
        ],
    }
    df, _ = load_rna_table(comp, aliases={}, warnings=warns)
    tbce_rows = df.loc[df["gene_key"] == "TBCE"]
    if not tbce_rows.empty:
        tbce = tbce_rows.iloc[0]
        assert not (tbce["rna_known_sig"] and tbce["rna_direction"] < 0), (
            f"Chimeric significant downregulated TBCE fabricated: {tbce.to_dict()}"
        )
    assert any("ambiguous" in w.lower() or "collision" in w.lower() or "duplicate" in w.lower() for w in warns)

    # 2. Full table with two distinct Ensembl loci sharing symbol TBCE where the winning
    # locus has NaN log2fc must NOT splice log2fc from the second locus.
    warns2: list[str] = []
    comp_splice = {
        "top_genes": [
            {"gene_id": "ENSG00000285053", "gene_symbol": "TBCE", "log2FoldChange": None, "padj": 0.001},
            {"gene_id": "ENSG00000284770", "gene_symbol": "TBCE", "log2FoldChange": -3.0, "padj": 0.90},
            {"gene_id": "ENSG00000141510", "gene_symbol": "TP53", "log2FoldChange": 1.8, "padj": 0.001},
        ],
    }
    df2, _ = load_rna_table(comp_splice, aliases={}, warnings=warns2)
    tbce2 = df2.loc[df2["gene_key"] == "TBCE"]
    if not tbce2.empty:
        assert pd.isna(tbce2.iloc[0]["rna_log2fc"]), (
            f"Spliced log2fc across loci: {tbce2.iloc[0].to_dict()}"
        )
