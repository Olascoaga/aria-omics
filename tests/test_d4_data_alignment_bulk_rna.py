"""Regression tests for D4 silent data misalignment fixes in bulk RNA-seq."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def test_d4_001_version_suffix_stripping_eliminates_false_positive_de(tmp_path: Path):
    """D4-001: Mixed versioned/unversioned Ensembl IDs across per-sample files
    must normalize by stripping Ensembl versions before joining, producing 0 FP DE."""
    pytest.importorskip("pydeseq2")
    from aria.scripts.rna_bulk.counts import _load_counts
    from aria.scripts.rna_bulk.deseq2 import _run_deseq2

    rng = np.random.default_rng(123)
    n_shared, n_diff_ver = 60, 15
    base_expr = rng.poisson(250, size=(n_shared + n_diff_ver, 6))
    samples = ["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"]
    files = []
    for col_idx, sname in enumerate(samples):
        p = tmp_path / f"{sname}.tsv"
        rows = [f"ENSG{g:011d}\t{base_expr[g, col_idx]}" for g in range(n_shared)]
        for g in range(n_shared, n_shared + n_diff_ver):
            gid = f"ENSG{g:011d}.1" if sname.startswith("ctrl") else f"ENSG{g:011d}"
            rows.append(f"{gid}\t{base_expr[g, col_idx]}")
        p.write_text("\n".join(rows) + "\n", encoding="utf-8")
        files.append(str(p))

    counts, warns, meta = _load_counts(files)
    assert counts is not None, f"Expected clean merge after version stripping: {warns} {meta}"
    assert len(counts) == n_shared + n_diff_ver
    cu = meta.get("count_union")
    assert isinstance(cu, dict)
    assert cu["policy"] == "error"
    assert cu["n_genes_union"] == n_shared + n_diff_ver
    assert cu["n_genes_intersection"] == n_shared + n_diff_ver
    assert all(v == 0 for v in cu["per_file_missing"].values())

    meta_df = pd.DataFrame({"condition": ["ctrl"] * 3 + ["treat"] * 3}, index=samples)
    de_res, _ = _run_deseq2(counts, meta_df, "condition", "treat", "ctrl", 0.05, 0.5)
    assert de_res["n_sig"] == 0, (
        f"Fabricated {de_res['n_sig']} false-positive DE genes "
        f"(matrix={len(counts)} genes; warns={warns})"
    )


def test_d4_001_mismatched_gene_sets_fail_closed_and_support_explicit_policies(tmp_path: Path):
    """D4-001: Genuinely mismatched gene sets fail closed by default ('error'),
    naming count_union_policy + options + symmetric-difference examples; explicit
    'intersect' and 'union_zero_fill' record structured count_union metadata."""
    from aria.scripts.rna_bulk.counts import _load_counts

    f1 = tmp_path / "ctrl_1.tsv"
    f2 = tmp_path / "treat_1.tsv"
    f1.write_text(
        "ENSG00000000001\t100\nENSG00000000002\t200\nENSG00000000003\t300\n",
        encoding="utf-8",
    )
    f2.write_text(
        "ENSG00000000001\t110\nENSG00000000002\t190\nENSG00000000004\t400\n",
        encoding="utf-8",
    )
    files = [str(f1), str(f2)]

    # Default policy='error' must fail closed with actionable error message.
    counts_err, warns_err, meta_err = _load_counts(files)
    assert counts_err is None
    assert meta_err.get("refused") is True
    assert meta_err.get("error_type") == "CountGeneSetMismatch"
    details = meta_err.get("details", "")
    assert "count_union_policy" in details
    assert "intersect" in details and "union_zero_fill" in details and "error" in details
    assert "ENSG00000000003" in details or "ENSG00000000004" in details
    assert meta_err["count_union"]["n_genes_union"] == 4
    assert meta_err["count_union"]["n_genes_intersection"] == 2
    assert meta_err["count_union"]["per_file_missing"] == {"ctrl_1": 1, "treat_1": 1}

    # Explicit policy='intersect': keeps only the 2 shared genes, zero fabrication.
    counts_int, warns_int, meta_int = _load_counts(files, count_union_policy="intersect")
    assert counts_int is not None
    assert sorted(counts_int.index.tolist()) == ["ENSG00000000001", "ENSG00000000002"]
    assert meta_int["count_union"]["policy"] == "intersect"
    assert meta_int["count_union"]["n_genes_union"] == 4
    assert meta_int["count_union"]["n_genes_intersection"] == 2
    assert meta_int["count_union"]["per_file_missing"] == {"ctrl_1": 1, "treat_1": 1}

    # Explicit policy='union_zero_fill': explicit zero-fill with structured record.
    counts_uzf, warns_uzf, meta_uzf = _load_counts(files, count_union_policy="union_zero_fill")
    assert counts_uzf is not None
    assert len(counts_uzf) == 4
    assert int(counts_uzf.loc["ENSG00000000003", "treat_1"]) == 0
    assert int(counts_uzf.loc["ENSG00000000004", "ctrl_1"]) == 0
    assert meta_uzf["count_union"]["policy"] == "union_zero_fill"
    assert meta_uzf["count_union"]["per_file_missing"] == {"ctrl_1": 1, "treat_1": 1}

    # Invalid policy must fail closed.
    counts_bad, _, meta_bad = _load_counts(files, count_union_policy="guess")
    assert counts_bad is None and meta_bad.get("error_type") == "InvalidCountUnionPolicy"


def test_d4_001_bulk_narrator_surfaces_count_union_caveat_when_not_error():
    """D4-001: BulkRnaNarrator emits a warning caveat whenever count_union.policy != 'error'."""
    pytest.importorskip("litellm")
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    base_findings = {
        "sample_qc": {"n_samples": 4, "outliers": [], "size_ratio": 1.2},
        "contrasts": [{
            "name": "treat_vs_ctrl",
            "status": "success",
            "n_significant": 5,
            "n_upregulated": 3,
            "n_downregulated": 2,
            "top_genes": [{"gene": "G1", "symbol": "G1", "log2fc": 2.0, "padj": 0.001}],
        }],
    }
    narrator = BulkRnaNarrator()

    # policy == 'error': no count_union caveat emitted.
    blocks_err = narrator.collect(
        "bulk_rna_agent",
        {"findings": {**base_findings, "count_union": {
            "policy": "error", "n_genes_union": 100, "n_genes_intersection": 100,
            "per_file_missing": {"s1": 0, "s2": 0},
        }}},
    )
    assert not any("count_union_policy" in c.text for b in blocks_err for c in b.caveats)

    # policy == 'union_zero_fill' and 'intersect': warning caveat surfaced.
    for pol in ("intersect", "union_zero_fill"):
        blocks = narrator.collect(
            "bulk_rna_agent",
            {"findings": {**base_findings, "count_union": {
                "policy": pol, "n_genes_union": 100, "n_genes_intersection": 90,
                "per_file_missing": {"ctrl_1": 4, "treat_1": 6},
            }}},
        )
        matching = [
            c for b in blocks for c in b.caveats
            if pol in c.text and c.severity == "warning"
        ]
        assert matching, f"Expected warning caveat mentioning {pol!r} in BulkRnaNarrator blocks"
        assert any("100" in c.text and "90" in c.text for c in matching)


def test_d4_004_duplicate_metadata_and_count_ids_fail_closed(tmp_path: Path):
    """D4-004: Duplicate sample IDs in metadata TSV or count columns must fail
    closed rather than silently expanding rows and bypassing min_replicates=3."""
    from aria.scripts.rna_bulk.counts import (
        _enforce_metadata_correspondence,
        _load_or_infer_metadata,
    )

    rng = np.random.default_rng(7)
    mat = rng.poisson(200, size=(60, 5))
    mat[:10, 3:] = rng.poisson(800, size=(10, 2))
    counts = pd.DataFrame(
        mat,
        index=[f"ENSG{i:011d}" for i in range(60)],
        columns=["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2"],
    )
    meta_tsv = tmp_path / "metadata.tsv"
    dup_meta = pd.DataFrame(
        {"condition": ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"]},
        index=["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_2"],
    )
    dup_meta.to_csv(meta_tsv, sep="\t")

    # _load_or_infer_metadata must reject duplicate sample IDs naming 'treat_2'.
    with pytest.raises(ValueError, match="treat_2"):
        _load_or_infer_metadata(counts, str(meta_tsv), "condition")

    # _enforce_metadata_correspondence must reject duplicate metadata index AND duplicate count cols.
    with pytest.raises(ValueError, match="treat_2"):
        _enforce_metadata_correspondence(counts, dup_meta)

    dup_counts = counts[["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_2"]]
    clean_meta = pd.DataFrame(
        {"condition": ["ctrl", "ctrl", "ctrl", "treat", "treat"]},
        index=["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2"],
    )
    with pytest.raises(ValueError, match="treat_2"):
        _enforce_metadata_correspondence(dup_counts, clean_meta)

    # Direct _run_deseq2 call with duplicate metadata must also fail closed.
    from aria.scripts.rna_bulk.deseq2 import _run_deseq2

    res, _ = _run_deseq2(
        counts, dup_meta, "condition", "treat", "ctrl", 0.05, 0.5,
        min_replicates_per_condition=3,
    )
    assert res.get("status") == "error" and res.get("error_type") == "DuplicateSampleIds"

