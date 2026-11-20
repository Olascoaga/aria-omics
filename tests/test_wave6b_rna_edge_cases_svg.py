"""Wave 6B regression tests for Domain 1: Bulk RNA + scRNA edge cases & SVG parity.

Covers:
  - D1-H6: Numeric condition levels (0/24/48) across bulk DE, contrasts, sensitivity,
    pseudobulk DE, and differential abundance
  - D1-H7: gtf_io skips self-mapped / unresolved Ensembl IDs in _load_symbol_map & _to_symbols
  - D1-H8: bulk QC handles non-finite mean Spearman r (all-zero sample), respects min_reps
    in outlier pruning, and only reports outlier sensitivity status='success' when >=1 succeeds
  - D1-H9: rna_umi_dedup supports dual UMIs (ACGTACGT_TGCATGCA) from fastp per_read/per_index
  - D1-H10 & D1-M6: clamped PCA / n_pcs in advise_resolution & trajectory, missing X_pca
    handling in trajectory, and lognorm / negative-scaled .X guards before normalize_total+log1p
  - D1-H11 & D1-H12: lazy output_dir in rna_align & rna_quantify; per-chromosome exon cap in
    _load_strand_intervals
  - D1-H13 & D1-M9: unified _normalize_species in cellcomm/celltypist/qc, categorical autocrine
    comparison, rna_qc initial_min_genes & JSON-safe upper bounds, multi-sample scrna/qc metadata
    + <100-cell guard, lowercase barcodes in pooled plate barcode files
  - D1-M1 & D1-M8: sibling .svg emission in UMAP/PAGA/scRNA figures, PAGA matrix sanitization
    & category validation, error-row rendering in pseudobulk HTML table, semicolon-joined pathway
    genes in TSV, non-finite log2fc guard in _top_directional_genes
  - D1-M2..M5, D1-M7, D1-M9: BulkRnaNarrator TPM/GSEA/tables/pathways guards, log2FoldChange_raw
    export, _strip_gene_version preservation of dotted symbols, scRNA Methods gating, Wilcoxon
    method check, positive cluster markers for cluster-identity ORA, DA fisher fallback method,
    and pseudobulk batch_cov + numerator/denominator normalization.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

try:
    import anndata as ad  # noqa: F401
    _HAS_ANNDATA = True
except Exception:
    _HAS_ANNDATA = False

try:
    import scanpy as sc  # noqa: F401
    _HAS_SCANPY = True
except Exception:
    _HAS_SCANPY = False

try:
    import matplotlib  # noqa: F401
    _HAS_MPL = True
except Exception:
    _HAS_MPL = False

try:
    import litellm  # noqa: F401
    _HAS_LITELLM = True
except Exception:
    _HAS_LITELLM = False

needs_anndata = pytest.mark.skipif(not _HAS_ANNDATA, reason="anndata required")
needs_scanpy = pytest.mark.skipif(not _HAS_SCANPY, reason="scanpy required")
needs_mpl = pytest.mark.skipif(not _HAS_MPL, reason="matplotlib required")
needs_agent = pytest.mark.skipif(not _HAS_LITELLM, reason="litellm required for agents")


# ── D1-H6: Numeric condition levels (0 / 24 / 48) ────────────────────────────

def test_d1_h6_bulk_numeric_condition_levels_and_contrasts(tmp_path, monkeypatch):
    from aria.scripts import rna_bulk_de as rbd
    from aria.scripts.rna_bulk import contrasts, qc

    samples = ["s0_1", "s0_2", "s0_3", "s24_1", "s24_2", "s24_3", "s48_1", "s48_2", "s48_3"]
    meta = pd.DataFrame({"timepoint": [0, 0, 0, 24, 24, 24, 48, 48, 48]}, index=samples)

    # _suggest_contrasts must not crash with AttributeError on int levels
    sugg = contrasts._suggest_contrasts(meta, "timepoint")
    assert len(sugg) == 2
    assert all(isinstance(c["numerator"], str) and isinstance(c["denominator"], str) for c in sugg)

    # Full bulk_rna_de with numeric TSV metadata column
    counts_path = tmp_path / "counts.tsv"
    meta_path = tmp_path / "metadata.tsv"
    genes = [f"G{i}" for i in range(12)]
    pd.DataFrame(50, index=genes, columns=samples).to_csv(counts_path, sep="\t")
    meta.to_csv(meta_path, sep="\t")

    seen_calls = []

    def fake_deseq2(counts, metadata, design_factor, num, den, *a, **kw):
        seen_calls.append((num, den, list(metadata[design_factor].unique())))
        df = pd.DataFrame(
            {"baseMean": [50.0] * len(counts), "log2FoldChange": [1.5] * len(counts),
             "log2FoldChange_raw": [1.8] * len(counts), "pvalue": [0.001] * len(counts),
             "padj": [0.01] * len(counts)},
            index=counts.index,
        )
        return {
            "status": "success", "results": df,
            "sig_genes": list(counts.index[:2]), "up_genes": list(counts.index[:2]),
            "down_genes": [], "n_sig": 2, "n_up": 2, "n_down": 0,
        }, []

    monkeypatch.setattr(rbd, "_run_deseq2", fake_deseq2)
    monkeypatch.setattr(qc, "_run_deseq2", fake_deseq2)
    monkeypatch.setattr(rbd, "_sample_qc", lambda *a, **kw: {"outliers": [], "n_samples": 9})
    monkeypatch.setattr(rbd, "_generate_plots", lambda **kw: {})
    monkeypatch.setattr(rbd, "_load_symbol_map", lambda *a, **kw: {})
    monkeypatch.setattr(rbd, "_load_gene_annotation", lambda *a, **kw: {})
    monkeypatch.setattr(rbd, "compute_counts_lobo_jackknife", lambda *a, **kw: {"status": "skipped"})

    res = rbd.bulk_rna_de({
        "files": [str(counts_path)],
        "metadata_file": str(meta_path),
        "design_factor": "timepoint",
        "contrasts": [{"numerator": 24, "denominator": 0, "name": "24 vs 0"}],
        "run_pathways": False,
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "success", res
    assert seen_calls and seen_calls[0][:2] == ("24", "0")
    assert set(seen_calls[0][2]) == {"0", "24", "48"}

    # _run_outlier_sensitivity with numeric metadata column directly
    counts_df = pd.DataFrame(50, index=genes, columns=samples)
    sens = qc._run_outlier_sensitivity(
        counts=counts_df, metadata=meta, design_factor="timepoint",
        contrasts_in=[{"numerator": 24, "denominator": 0, "name": "24 vs 0"}],
        flagged_outliers=["s0_1"], removable_outliers=["s0_1"],
        primary_contrasts=res["contrasts"],
        padj_thr=0.05, lfc_thr=1.0, allow_mock=False, min_reps=2,
        covariates=[], lfc_shrink=False, fdr_family={"fdr_family": "per_contrast"},
        warnings=[],
    )
    assert sens["status"] == "success"
    assert sens["contrasts"][0]["status"] == "success"


@needs_anndata
def test_d1_h6_pseudobulk_and_diff_abundance_numeric_comparisons(tmp_path, monkeypatch):
    from aria.scripts import rna_diff_abundance as da_mod
    from aria.scripts import rna_pseudobulk_de as pb_mod

    n_cells = 24
    obs = pd.DataFrame({
        "cell_type": ["A", "B"] * (n_cells // 2),
        "timepoint": [0] * 12 + [24] * 12,
        "donor": [f"d{i // 4}" for i in range(n_cells)],
    }, index=[f"c{i}" for i in range(n_cells)])
    rng = np.random.default_rng(11)
    adata = ad.AnnData(
        X=rng.integers(5, 50, size=(n_cells, 4)),
        obs=obs,
        var=pd.DataFrame(index=["G0", "G1", "G2", "G3"]),
    )
    h5 = tmp_path / "numeric_cond.h5ad"
    adata.write_h5ad(h5)

    # Differential abundance with numeric comparison [[24, 0]]
    da_res = da_mod.rna_diff_abundance({
        "data_path": str(h5),
        "groupby": "cell_type",
        "condition_col": "timepoint",
        "replicate_col": "donor",
        "comparisons": [[24, 0]],
        "output_dir": str(tmp_path / "da"),
    })
    assert da_res["status"] == "success"
    assert da_res["per_comparison"]["24_vs_0"]["status"] == "success"

    # Pseudobulk DE with numeric comparison [[24, 0]]
    class _FakeDDS:
        def __init__(self, counts, metadata, **kw):
            self.counts = counts

        def deseq2(self):
            pass

    class _FakeStats:
        def __init__(self, dds, **kw):
            self.results_df = pd.DataFrame(
                {"pvalue": [0.001, 0.5, 0.5, 0.5], "padj": [0.004, 0.8, 0.8, 0.8],
                 "log2FoldChange": [1.5, 0.0, 0.0, 0.0], "baseMean": [50.0] * 4},
                index=dds.counts.columns,
            )

        def summary(self):
            pass

    monkeypatch.setitem(sys.modules, "pydeseq2.dds", SimpleNamespace(DeseqDataSet=_FakeDDS))
    monkeypatch.setitem(sys.modules, "pydeseq2.ds", SimpleNamespace(DeseqStats=_FakeStats))
    monkeypatch.setattr(pb_mod, "compute_counts_lobo_jackknife", lambda *a, **kw: {"status": "skipped"})

    pb_res = pb_mod.rna_pseudobulk_de({
        "data_path": str(h5),
        "groupby": "cell_type",
        "condition_col": "timepoint",
        "replicate_col": "donor",
        "comparisons": [[24, 0]],
        "min_cells_per_pseudosample": 2,
        "min_replicates_per_condition": 3,
        "use_raw": False,
        "lfc_shrink": False,
        "output_dir": str(tmp_path / "pb"),
    })
    assert pb_res["status"] == "success"
    assert pb_res["per_group"]["A"]["per_comparison"]["24_vs_0"]["status"] == "success"


# ── D1-H7: gtf_io skips self-mapped & unresolved Ensembl IDs ─────────────────

def test_d1_h7_gtf_io_skips_self_mapped_and_unresolved_ensembl(tmp_path):
    from aria.scripts.rna_bulk import gtf_io

    counts = tmp_path / "counts.tsv"
    counts.write_text("gene\ts1\nENSMUSG000001.7\t10\nENSMUSG000002.1\t20\n", encoding="utf-8")
    sym_file = tmp_path / "counts_with_symbols.tsv"
    pd.DataFrame({
        "gene_symbol": ["Sox2", "ENSMUSG000002.1", "ENSG000003", "TP53"],
        "s1": [10, 20, 30, 40],
    }, index=["ENSMUSG000001.7", "ENSMUSG000002.1", "ENSG000003.4", "TP53"]).to_csv(
        sym_file, sep="\t"
    )

    loaded = gtf_io._load_symbol_map([str(counts)], [])
    assert loaded == {"ENSMUSG000001": "Sox2"}

    # Even if a stale map contains self-mapped Ensembl entries, _to_symbols drops them
    stale_map = {"ENSMUSG000001": "Sox2", "ENSMUSG000002": "ENSMUSG000002.1"}
    assert gtf_io._to_symbols(
        ["ENSMUSG000001.7", "ENSMUSG000002.1", "ENSG000099.2", "HLA.DRA"], stale_map
    ) == ["Sox2", "HLA.DRA"]
    # And with empty symbol_map, raw Ensembl IDs are still excluded
    assert gtf_io._to_symbols(["ENSG000001.1", "ENSMUSG000002", "Sox2"], {}) == ["Sox2"]


# ── D1-H8: bulk QC non-finite mean_r, min_reps pruning, sensitivity status ───

def test_d1_h8_all_zero_sample_min_reps_and_sensitivity_status(tmp_path, monkeypatch):
    from aria.scripts.rna_bulk import qc

    rng = np.random.default_rng(3)
    samples = ["c1", "c2", "c3", "t1", "t2", "t3"]
    genes = [f"G{i}" for i in range(40)]
    mat = pd.DataFrame(rng.integers(10, 200, size=(40, 6)), index=genes, columns=samples)
    mat["c1"] = 0  # broken all-zero sample -> NaN Spearman correlation
    meta = pd.DataFrame({"condition": ["ctrl"] * 3 + ["treat"] * 3}, index=samples)

    monkeypatch.setattr(
        qc, "_plot_pca_mds",
        lambda vst_var, md, od, w, **kw: (None, None, np.zeros((6, 2)), [0.0, 0.0]),
    )
    warns = []
    sqc = qc._sample_qc(mat, meta, str(tmp_path), warns, design_factor="condition")
    assert sqc["replicate_correlations"]["c1"] is None
    assert "c1" in sqc["replicate_outliers"]
    assert "c1" in sqc["outliers"]
    # Strict JSON serialization (no NaN)
    json.dumps({k: v for k, v in sqc.items() if k != "vst_matrix"}, allow_nan=False)

    # _prune_outliers_for_design respects min_reps
    w3 = []
    assert qc._prune_outliers_for_design(["c1"], meta, "condition", w3, min_reps=3) == []
    assert any("<3 replicates" in w for w in w3)
    w2 = []
    assert qc._prune_outliers_for_design(["c1"], meta, "condition", w2, min_reps=2) == ["c1"]

    # _run_outlier_sensitivity must NOT report status='success' when all contrasts error
    monkeypatch.setattr(
        qc, "_run_deseq2",
        lambda *a, **kw: ({"status": "error", "details": "fit failed"}, ["err"]),
    )
    sens = qc._run_outlier_sensitivity(
        counts=mat, metadata=meta, design_factor="condition",
        contrasts_in=[{"numerator": "treat", "denominator": "ctrl", "name": "treat vs ctrl"}],
        flagged_outliers=["c1"], removable_outliers=["c1"],
        primary_contrasts=[{"name": "treat vs ctrl", "status": "success", "n_significant": 1}],
        padj_thr=0.05, lfc_thr=1.0, allow_mock=False, min_reps=2,
        covariates=[], lfc_shrink=False, fdr_family={"fdr_family": "per_contrast"},
        warnings=[],
    )
    assert sens["status"] == "error"
    assert sens["conclusion_robust"] is None


# ── D1-H9: Dual UMI support in rna_umi_dedup ─────────────────────────────────

def test_d1_h9_dual_umi_extraction_and_directional_collapse():
    from collections import Counter
    from aria.scripts import rna_umi_dedup as ud

    assert ud._extract_umi("read1:ACGTACGT") == "ACGTACGT"
    assert ud._extract_umi("read1:ACGTACGT_TGCATGCA") == "ACGTACGT_TGCATGCA"
    assert ud._extract_umi("read1:acgtacgt_tgcatgca/1") == "ACGTACGT_TGCATGCA"
    assert ud._extract_umi("read1:ACGT_TGCA_AAAA") is None
    assert ud._extract_umi("read1:ACGT_1234") is None

    winners = ud._directional_winning_umis(Counter({
        "ACGTACGT_TGCATGCA": 10,
        "ACGTACGA_TGCATGCA": 2,   # 1-MM error of first -> collapses
        "TTTTTTTT_CCCCCCCC": 5,   # distinct molecule -> survives
    }))
    assert winners == {"ACGTACGT_TGCATGCA", "TTTTTTTT_CCCCCCCC"}


# ── D1-H10 & D1-M6: PCA clamping, missing X_pca, lognorm/scaled .X guards ───

@needs_scanpy
@needs_anndata
def test_d1_h10_m6_advise_resolution_and_trajectory_guards(tmp_path, monkeypatch):
    import scanpy as sc
    from aria.scripts import rna_advise_resolution as adv
    from aria.scripts import rna_trajectory as traj
    from aria.scripts import rna_clustering as clu
    from aria.scripts import rna_integration as integ
    from aria.scripts import rna_figure_umap as fig_umap

    rng = np.random.default_rng(17)
    # 1. Already log-normalized matrix (non-negative, non-integer, max <= 50)
    log_x = np.log1p(rng.uniform(0.5, 10.0, size=(12, 8)))
    a_log = ad.AnnData(
        X=log_x.copy(),
        obs=pd.DataFrame({"batch": ["b1"] * 6 + ["b2"] * 6, "leiden": ["0"] * 6 + ["1"] * 6},
                         index=[f"c{i}" for i in range(12)]),
        var=pd.DataFrame(index=[f"G{i}" for i in range(8)]),
    )
    log_path = tmp_path / "lognorm.h5ad"
    a_log.write_h5ad(log_path)

    # 2. Scaled/negative matrix (no X_pca)
    scaled_x = rng.normal(0.0, 1.0, size=(12, 8))
    a_scaled = ad.AnnData(
        X=scaled_x.copy(),
        obs=a_log.obs.copy(),
        var=a_log.var.copy(),
    )
    scaled_path = tmp_path / "scaled.h5ad"
    a_scaled.write_h5ad(scaled_path)

    norm_calls = []
    orig_norm = sc.pp.normalize_total
    orig_log1p = sc.pp.log1p

    def spy_norm(*a, **kw):
        norm_calls.append("norm")
        return orig_norm(*a, **kw)

    def spy_log1p(*a, **kw):
        norm_calls.append("log1p")
        return orig_log1p(*a, **kw)

    monkeypatch.setattr(sc.pp, "normalize_total", spy_norm)
    monkeypatch.setattr(sc.pp, "log1p", spy_log1p)

    # rna_advise_resolution on small lognorm and scaled matrices: skips norm/log1p & clamps PCA
    for path in (log_path, scaled_path):
        norm_calls.clear()
        res = adv.rna_advise_resolution({
            "data_path": str(path),
            "resolutions": [0.5],
            "n_pcs": 40,
            "resume": False,
        })
        assert res["status"] == "success", res
        assert norm_calls == []

    # rna_trajectory without neighbors OR X_pca computes clamped PCA cleanly
    norm_calls.clear()
    t_res = traj.rna_trajectory({
        "data_path": str(scaled_path),
        "output_dir": str(tmp_path / "traj"),
    })
    assert t_res["status"] == "success", t_res
    assert norm_calls == []

    # rna_trajectory on degenerate 1-gene matrix returns clean skipped (not KeyError)
    tiny = ad.AnnData(
        X=np.ones((4, 1)),
        obs=pd.DataFrame({"leiden": ["0", "0", "1", "1"]}, index=[f"c{i}" for i in range(4)]),
        var=pd.DataFrame(index=["G0"]),
    )
    tiny_path = tmp_path / "tiny.h5ad"
    tiny.write_h5ad(tiny_path)
    t_tiny = traj.rna_trajectory({"data_path": str(tiny_path), "output_dir": str(tmp_path)})
    assert t_tiny["status"] == "skipped"

    # rna_clustering skips norm/log1p on scaled & lognorm .X
    for idx, path in enumerate((log_path, scaled_path)):
        norm_calls.clear()
        c_res = clu.rna_clustering({
            "data_path": str(path),
            "resolution": 0.5,
            "n_pcs": 40,
            "output_dir": str(tmp_path / f"clu_{idx}"),
        })
        assert c_res["status"] == "success", c_res
        assert norm_calls == []
        assert c_res["n_pcs"] == 7

    # rna_integration skips norm/log1p on scaled .X
    monkeypatch.setitem(
        sys.modules, "harmonypy",
        SimpleNamespace(run_harmony=lambda emb, *a, **kw: SimpleNamespace(Z_corr=emb)),
    )
    norm_calls.clear()
    i_res = integ.rna_integration({
        "data_path": str(scaled_path),
        "batch_col": "batch",
        "output_dir": str(tmp_path / "integ"),
    })
    assert i_res["status"] == "success", i_res
    assert norm_calls == []

    # rna_figure_umap._compute_umap skips norm/log1p on scaled .X
    norm_calls.clear()
    a_copy = a_scaled.copy()
    assert fig_umap._compute_umap(a_copy) == "X_umap"
    assert norm_calls == []
    assert np.isfinite(a_copy.obsm["X_umap"]).all()


# ── D1-H11 & D1-H12: Lazy output_dir & per-chromosome strand intervals ───────

def test_d1_h11_h12_lazy_output_dir_and_per_chrom_strand_intervals(tmp_path):
    from aria.scripts import rna_align, rna_quantify

    # Empty samples / failed upstream sample without r1_trimmed must not raise IndexError/KeyError
    r_empty = rna_align.rna_align({"samples": [], "genome_dir": str(tmp_path)})
    assert r_empty["status"] == "error" and r_empty["error_type"] == "NoSamples"

    r_failed = rna_align.rna_align({
        "samples": [{"name": "s1", "status": "failed"}],
        "genome_dir": str(tmp_path),
    })
    assert r_failed["status"] == "error" and r_failed["error_type"] == "UpstreamFastqQCFailure"

    # Empty bam_files / failed BAM without 'bam' key must not raise IndexError/KeyError
    q_empty = rna_quantify.rna_quantify({"bam_files": [], "gtf_file": str(tmp_path / "a.gtf")})
    assert q_empty["status"] == "error" and q_empty["error_type"] == "NoBamsFound"

    q_failed = rna_quantify.rna_quantify({
        "bam_files": [{"name": "s1", "status": "failed"}],
        "gtf_file": str(tmp_path / "a.gtf"),
    })
    assert q_failed["status"] == "error" and q_failed["error_type"] == "PartialAlignmentInput"

    # Per-chromosome exon cap: >5,000 exons on chr1 must NOT block chr2/chrM!
    gtf = tmp_path / "multi_chrom.gtf"
    lines = [
        f'chr1\tt\texon\t{i * 200 + 1}\t{i * 200 + 100}\t.\t+\t.\tgene_id "G{i}";\n'
        for i in range(6000)
    ]
    lines.append('chr2\tt\texon\t100\t200\t.\t-\t.\tgene_id "G_CHR2";\n')
    lines.append('chrM\tt\texon\t50\t150\t.\t+\t.\tgene_id "G_CHRM";\n')
    gtf.write_text("".join(lines), encoding="utf-8")

    compiled, err = rna_quantify._load_strand_intervals(str(gtf))
    assert err is None
    assert {"chr1", "chr2", "chrM"} <= set(compiled.keys())
    assert len(compiled["chr1"]["plus"]) == 5000
    assert compiled["chr2"]["minus"] == [(99, 200)]
    assert compiled["chrM"]["plus"] == [(49, 150)]


# ── D1-H13 & D1-M9: Species aliases, categorical autocrine, QC & barcodes ───

@pytest.mark.parametrize("alias", ["GRCm38", "GRCm39", "mm9", "m. musculus"])
def test_d1_h13_celltypist_mouse_aliases(alias):
    from aria.scripts.rna_celltypist import _resolve_celltypist_model

    model, source, _ = _resolve_celltypist_model(None, alias, "brain")
    assert (model, source) == ("Mouse_Whole_Brain.pkl", "tissue_hint")


@needs_anndata
def test_d1_h13_cellcomm_mouse_alias_and_categorical_autocrine(tmp_path, monkeypatch):
    from aria.scripts.rna_cellcomm import rna_cellcomm

    a = ad.AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"]}, index=["c0", "c1", "c2", "c3"]),
        var=pd.DataFrame(index=["G0", "G1"]),
    )
    h5 = tmp_path / "cc.h5ad"
    a.write_h5ad(h5)

    # Mismatched Categorical levels on source vs target used to raise TypeError on ==
    frame = pd.DataFrame({
        "source": pd.Categorical(["A", "A"], categories=["A", "C"]),
        "target": pd.Categorical(["A", "B"], categories=["A", "B"]),
        "ligand_complex": ["L1", "L2"],
        "receptor_complex": ["R1", "R2"],
        "specificity_rank": [0.01, 0.02],
    })
    seen = {}

    def fake_rank(matrix, **kw):
        seen.update(kw)
        matrix.uns["liana_res"] = frame

    monkeypatch.setitem(sys.modules, "liana", SimpleNamespace(mt=SimpleNamespace(rank_aggregate=fake_rank)))
    res = rna_cellcomm({"data_path": str(h5), "organism": "GRCm39", "output_dir": str(tmp_path)})
    assert res["status"] == "success"
    assert seen["resource_name"] == "mouseconsensus"
    assert res["n_autocrine_dropped"] == 1
    assert res["n_interactions"] == 1


@needs_scanpy
@needs_anndata
def test_d1_h13_rna_qc_min_genes_mouse_alias_and_json_bounds(tmp_path):
    from aria.scripts.rna_qc import _cache_params, rna_qc

    # Default initial_min_genes follows min_genes when min_genes < 200
    assert _cache_params({"min_genes": 50})["initial_min_genes"] == 50
    assert _cache_params({"min_genes": 500})["initial_min_genes"] == 200

    # Matrix with 60 genes expressed per cell (below 200, above min_genes=50) + mouse mt- gene
    rng = np.random.default_rng(5)
    x = rng.integers(1, 5, size=(10, 60)).astype(float)
    genes = ["mt-Nd1"] + [f"Gene{i}" for i in range(59)]
    a = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(10)]),
        var=pd.DataFrame(index=genes),
    )
    p = tmp_path / "low_gene_mouse.h5ad"
    a.write_h5ad(p)

    res = rna_qc({
        "data_path": str(p),
        "organism": "GRCm38",
        "min_genes": 50,
        "min_cells": 1,
        "initial_min_cells": 1,
        "run_scrublet": False,
        "output_dir": str(tmp_path / "qc_out"),
    })
    assert res["status"] == "success", res
    assert res["n_cells_after"] == 10
    assert not any("case-insensitive" in w for w in res["warnings"])
    json.dumps(res, allow_nan=False)

    # Non-finite max_counts/max_genes serialize as None, never Infinity
    res_inf = rna_qc({
        "data_path": str(p),
        "organism": "GRCm38",
        "min_genes": 50,
        "min_cells": 1,
        "initial_min_cells": 1,
        "max_counts": float("inf"),
        "max_genes": float("inf"),
        "run_scrublet": False,
        "output_dir": str(tmp_path / "qc_inf"),
    })
    assert res_inf["status"] == "success"
    assert res_inf["max_counts_used"] is None
    assert res_inf["max_genes_used"] is None
    json.dumps(res_inf, allow_nan=False)


@needs_agent
def test_d1_h13_multisample_scrna_qc_preserves_metadata_and_enforces_100_cell_guard(tmp_path):
    from aria.agents.scrna_agent import scRNAAgent

    def _run_multi(cells_per_sample: int):
        class _Env:
            def run_in_stack(self, *, stack, script_path, params):
                if script_path.endswith("rna_qc.py"):
                    sid = params["sample_id"]
                    return {
                        "status": "success",
                        "output_path": str(tmp_path / f"qc_{sid}.h5ad"),
                        "n_cells_before": cells_per_sample + 10,
                        "n_cells_after": cells_per_sample,
                        "pct_removed": 10.0,
                        "mt_threshold_used": 14.5,
                        "ambient_correction": {"ran": False, "method": None, "reason": "not_requested"},
                        "scrublet": {"ran": True, "n_doublets": 3},
                    }
                return {
                    "status": "success",
                    "output_path": str(tmp_path / "concat.h5ad"),
                    "n_cells_total": cells_per_sample * 2,
                    "n_genes_shared": 500,
                    "batch_col": "batch",
                }

        agent = scRNAAgent.__new__(scRNAAgent)
        agent.env = _Env()
        agent._workspace = lambda exp_id, *parts: tmp_path
        agent.publish_finding = lambda *a, **kw: None
        return agent._run_qc("exp", ["/a/s1.h5ad", "/b/s2.h5ad"], {"organism": "Homo sapiens"}, {})

    ok = _run_multi(80)
    assert ok["status"] == "success"
    assert ok["mt_threshold_used"] == 14.5
    assert ok["ambient_correction"] == {"ran": False, "method": None, "reason": "not_requested"}
    assert ok["scrublet"] == {"ran": True, "n_doublets": 6}

    too_few = _run_multi(30)  # 60 total < 100
    assert too_few["status"] == "error"
    assert too_few["error_type"] == "InsufficientCells"
    assert too_few["mt_threshold_used"] == 14.5


def test_d1_m9_lowercase_barcodes_in_pooled_plate_file(tmp_path):
    from aria.utils.rna_fastq_protocol import validate_pooled_plate_manifest

    bc_file = tmp_path / "barcodes.tsv"
    bc_file.write_text("barcode\tsample_id\nacgtacgtac\tS1\ntgcatgcatg\tS2\n", encoding="utf-8")
    manifest, errors = validate_pooled_plate_manifest({
        "preset": "drugseq",
        "barcode_to_sample": str(bc_file),
    })
    assert errors == []
    assert manifest["barcode_to_sample"] == {"ACGTACGTAC": "S1", "TGCATGCATG": "S2"}


# ── D1-M1 & D1-M8: SVG parity, PAGA sanitization, scRNA tables & _common ─────

@needs_mpl
def test_d1_m1_m8_script_svg_parity(tmp_path):
    from aria.scripts import rna_figure_umap as fig_umap
    from aria.scripts import rna_figure_paga as fig_paga

    # UMAP _plot_one emits sibling SVG
    coords = np.array([[0.0, 0.0], [1.0, 1.0], [0.5, 0.2], [1.5, 0.8]])
    umap_png = Path(fig_umap._plot_one(
        coords, np.array(["A", "A", "B", "B"]), "leiden",
        tmp_path / "umap_leiden.png", 10.0, "UMAP — leiden", "UMAP",
    ))
    assert umap_png.is_file() and umap_png.with_suffix(".svg").is_file()

    # PAGA _draw_paga emits sibling SVG and handles non-finite conn
    conn = np.array([[1.0, np.nan, 0.2], [np.nan, 1.0, np.inf], [0.2, np.inf, 1.0]])
    paga_png = Path(fig_paga._draw_paga(
        conn, ["0", "1", "2"], [10, 20, 30], tmp_path / "paga_graph.png",
    ))
    assert paga_png.is_file() and paga_png.with_suffix(".svg").is_file()


@needs_mpl
@needs_agent
def test_d1_m1_scrna_narrative_svg_parity(tmp_path):
    from aria.agents.narrative.scrna import figures as scrna_figs

    # scRNA narrative figures emit sibling SVGs
    findings = {
        "cell_communication": {
            "method": "liana_rank_aggregate (specificity_rank)",
            "top_interactions": [
                {"source": "T", "target": "B", "ligand": "L1", "receptor": "R1",
                 "score": 0.01, "rank_metric": "specificity_rank"},
                {"source": "B", "target": "T", "ligand": "L2", "receptor": "R2",
                 "score": 0.02, "rank_metric": "specificity_rank"},
            ],
        },
        "pseudobulk_de": {
            "per_group": {"T": {"per_comparison": {
                "trt_vs_ctrl": {"status": "success", "n_up": 5, "n_down": 3}
            }}},
        },
    }
    for fn, name in (
        (scrna_figs.render_cellcomm_heatmap, "cc_heat.png"),
        (scrna_figs.render_cellcomm_top_pairs_bar, "cc_bar.png"),
        (scrna_figs.render_per_celltype_de_bar, "de_bar.png"),
    ):
        png_path = Path(fn(findings, tmp_path / name))
        assert png_path.is_file() and png_path.with_suffix(".svg").is_file()


@needs_anndata
@needs_mpl
def test_d1_m8_paga_make_figures_defaults_and_shape_guard(tmp_path):
    from aria.scripts import rna_figure_paga as fig_paga

    a = ad.AnnData(
        X=np.ones((6, 2)),
        obs=pd.DataFrame({
            "cell_type": pd.Categorical(["T", "T", "B", "B", "M", "M"]),
            "leiden": pd.Categorical(["0", "0", "0", "1", "1", "1"]),
        }, index=[f"c{i}" for i in range(6)]),
        var=pd.DataFrame(index=["G0", "G1"]),
    )
    a.uns["paga"] = {
        "groups": "cell_type",
        "connectivities": np.array([[0.0, 0.02, np.nan], [0.02, 0.0, np.inf], [np.nan, np.inf, 0.0]]),
    }
    p = tmp_path / "traj.h5ad"
    a.write_h5ad(p)

    # Defaults groupby to uns['paga']['groups'] and sanitizes NaN/inf
    ok = fig_paga.make_paga_figures({"h5ad_path": str(p), "output_dir": str(tmp_path / "ok")})
    assert ok["status"] == "success"
    assert ok["groupby"] == "cell_type"
    assert ok["max_connectivity"] == pytest.approx(0.02)
    assert Path(ok["figures"]["paga_graph"]).with_suffix(".svg").is_file()

    # Passing mismatched groupby ('leiden' has 2 cats vs 3x3 conn) returns structured error
    bad = fig_paga.make_paga_figures({
        "h5ad_path": str(p), "output_dir": str(tmp_path / "bad"), "groupby": "leiden",
    })
    assert bad["status"] == "error" and bad["error_type"] == "PAGAShapeMismatch"


@needs_agent
def test_d1_m8_scrna_tables_error_rows_semicolon_genes_and_nonfinite_lfc(tmp_path):
    from aria.agents.narrative.scrna import _common, tables

    findings = {
        "pseudobulk_de": {
            "per_group": {
                "Tcell": {
                    "n_pseudosamples": 6,
                    "per_comparison": {
                        "trt_vs_ctrl": {"status": "error", "reason": "singular design"},
                    },
                }
            }
        },
        "pseudobulk_pathways": {
            "per_cluster": {
                "Tcell::trt_vs_ctrl": {
                    "results": {
                        "GO_BP": [{"term": "T cell activation", "padj": 0.01, "genes": ["CD3D", "CD3E"]}],
                    }
                }
            }
        },
    }
    html_rows = tables.extract_pseudobulk_de_table(findings)
    assert "error:" in html_rows and "singular design" in html_rows

    exported = tables.export_supplementary_tables(findings, tmp_path)
    tsv_text = Path(exported["pathway_enrichment"]).read_text(encoding="utf-8")
    assert "CD3D;CD3E" in tsv_text
    assert "['CD3D'" not in tsv_text

    # _top_directional_genes must not crash on None / NaN / inf log2fc
    comp = {
        "top_genes": [
            {"gene": "BadNone", "log2fc": None, "padj": 0.01},
            {"gene": "BadNan", "log2fc": float("nan"), "padj": 0.01},
            {"gene": "BadInf", "log2fc": float("inf"), "padj": 0.01},
            {"gene": "GoodUp", "log2fc": 2.5, "padj": 0.001},
            {"gene": "GoodDn", "log2fc": -1.8, "padj": 0.002},
        ]
    }
    up = _common._top_directional_genes(comp, "up")
    dn = _common._top_directional_genes(comp, "down")
    assert len(up) == 1 and "GoodUp" in up[0]
    assert len(dn) == 1 and "GoodDn" in dn[0]


# ── D1-M2..M5, D1-M7, D1-M9: Narrators, exports, methods, markers, DA, PB ───

@needs_agent
def test_d1_m2_bulk_narrator_tpm_gsea_tables_and_pathways(tmp_path):
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    narrator = BulkRnaNarrator()
    # Full-length run where TPM was not computed (e.g., no GTF) must NOT claim tag/UMI bias
    findings_fl = {
        "tpm_status": {"computed": False, "tpm_path": None, "library_protocol": "full_length"},
        "contrasts": [{"name": "c1", "status": "success"}],
    }
    methods_fl = " ".join(narrator.methods("bulk_rna_agent", {"findings": findings_fl}))
    assert "terminal tags or UMI" not in methods_fl

    # GSEA table with 0 pathways at FDR < 0.25 must be reflected honestly in claim
    gsea_csv = tmp_path / "gsea.csv"
    gsea_csv.write_text("Term,NES,FDR q-val\nPathA,1.2,0.45\n", encoding="utf-8")
    contrast = {
        "name": "trt_vs_ctrl",
        "status": "success",
        "pathways": {"GO_BP": [{"term": "immune", "padj": 0.01}, "not-a-dict"], "meta": "not-a-list"},
        "plots": {
            "tables": {"de_genes": None, "pathways": "/tmp/pw.tsv"},
            "gsea_table": str(gsea_csv),
            "ora_dotplots": {"GO_BP": None},
            "gsea_running_sums": [None],
        },
    }
    blocks = narrator.collect("bulk_rna_agent", {"findings": {"contrasts": [contrast]}})
    gsea_block = next(b for b in blocks if b.analysis == "gsea_preranked")
    assert "0 significant pathways" in gsea_block.claim
    assert "ranked pathway support" not in gsea_block.claim

    tbls = narrator.tables("bulk_rna_agent", {"findings": {"contrasts": [contrast]}})
    assert all(isinstance(t["path"], str) and t["path"] for t in tbls)
    c_block = next(b for b in blocks if b.analysis == "differential_expression")
    assert all(isinstance(t["path"], str) and t["path"] for t in c_block.tables)


def test_d1_m3_m4_lfc_raw_export_finite_heatmap_and_dotted_symbols(tmp_path):
    from aria.scripts import rna_pathway_viz as pviz
    from aria.scripts.rna_bulk import contrasts, plots, transforms

    # Dotted symbols preserved (not split at '.') while Ensembl versions ARE stripped
    df = pd.DataFrame({
        "baseMean": [100.0, 200.0, 150.0],
        "log2FoldChange": [2.0, -1.5, float("inf")],
        "log2FoldChange_raw": [2.4, -1.9, float("inf")],
        "padj": [0.001, 0.002, 0.0001],
    }, index=["ENSG000001.12", "HLA.DRA", "ENSG000002.1"])
    smap = {"ENSG000001": "TP53", "HLA.DRA": "HLA.DRA"}

    # export_de_table preserves log2FoldChange_raw and HLA.DRA
    out_tsv = tmp_path / "de_genes.tsv"
    res_path = pviz.export_de_table(df, smap, "trt vs ctrl", str(out_tsv), padj_thr=0.05)
    assert res_path is not None
    exported = pd.read_csv(res_path, sep="\t")
    assert "log2FoldChange_raw" in exported.columns
    assert set(exported["gene_symbol"]) >= {"TP53", "HLA.DRA"}

    # _ranked_signature_frame preserves HLA.DRA
    sig = pviz._ranked_signature_frame(df, smap)
    assert set(sig["gene"]) == {"TP53", "HLA.DRA"}

    # _format_top_genes preserves HLA.DRA
    top = contrasts._format_top_genes({"results": df, "sig_genes": ["ENSG000001.12", "HLA.DRA"]}, smap)
    assert [t["symbol"] for t in top] == ["TP53", "HLA.DRA"]

    # transforms preserve dotted gene symbols in length/biotype lookups
    counts = pd.DataFrame({"s1": [100, 100], "s2": [120, 80]}, index=["ENSG000001.12", "HLA.DRA"])
    tpm = transforms._compute_tpm(counts, {"ENSG000001": 1000, "HLA.DRA": 2000}, [])
    assert tpm is not None and list(tpm.index) == ["ENSG000001.12", "HLA.DRA"]

    # _generate_plots filters non-finite log2FC in heatmap_lfc
    if _HAS_MPL:
        hm_counts = pd.DataFrame(
            {"s1": [10, 20, 30], "s2": [40, 10, 25]},
            index=["ENSG000001.12", "HLA.DRA", "ENSG000002.1"],
        )
        p_out = plots._generate_plots(
            de_result={"results": df, "sig_genes": ["ENSG000001.12", "HLA.DRA"], "n_up": 1, "n_down": 1},
            sample_qc={"_symbol_map": smap},
            counts_filt=hm_counts,
            metadata=pd.DataFrame({"condition": ["a", "b"]}, index=["s1", "s2"]),
            design_factor="condition",
            output_dir=str(tmp_path / "figs"),
            padj_thr=0.05,
            lfc_thr=1.0,
        )
        assert "heatmap_lfc" in p_out


@needs_agent
def test_d1_m5_scrna_methods_gating():
    from aria.agents.narrative.scrna.text import build_scrna_methods

    # 1. Scrublet ran=False, predef_clusters=True, lognorm_recovered=False
    f1 = {
        "qc": {"n_cells_after": 500, "scrublet": {"ran": False}},
        "clustering": {"n_clusters": 4, "predef_clusters": True, "groupby": "cell_type"},
        "pseudobulk_de": {"groupby": "cell_type", "lognorm_recovered": False, "per_group": {}},
    }
    m1 = build_scrna_methods(f1)
    assert "Doublets were flagged with Scrublet" not in m1
    assert "PCA (50 components)" not in m1
    assert "k-NN graph construction" not in m1
    assert "For Seurat-derived h5ads" not in m1

    # 2. Scrublet ran=True, custom n_pcs=30, lognorm_recovered=True
    f2 = {
        "qc": {"n_cells_after": 500, "scrublet": {"ran": True}},
        "clustering": {"n_clusters": 5, "n_pcs": 30, "resolution": 0.6},
        "pseudobulk_de": {"groupby": "cell_type", "lognorm_recovered": True, "per_group": {}},
    }
    m2 = build_scrna_methods(f2)
    assert "Doublets were flagged with Scrublet." in m2
    assert "PCA (30 components)" in m2
    assert "PCA (50 components)" not in m2
    assert "For Seurat-derived h5ads" in m2


@needs_anndata
def test_d1_m7_de_per_cluster_rejects_non_wilcoxon_and_ora_prefers_positive_markers(tmp_path, monkeypatch):
    import types
    import aria.utils.safe_h5ad as safe_h5ad
    from aria.scripts import rna_de_per_cluster as de_mod
    from aria.scripts import rna_pathway_per_cluster as pw_mod
    from aria.utils import ora as local_ora

    # 1. Cached rank_genes_groups with method='t-test' must be rejected and recomputed
    dtype = [("0", "O"), ("1", "O")]
    rgg_ttest = {
        "params": {"groupby": "leiden", "use_raw": False, "method": "t-test"},
        "names": np.array([("G0", "G1")], dtype=dtype),
        "logfoldchanges": np.array([(2.0, 2.0)], dtype=dtype),
        "pvals_adj": np.array([(0.001, 0.001)], dtype=dtype),
    }
    rgg_wilcox = dict(rgg_ttest, params={"groupby": "leiden", "use_raw": False, "method": "wilcoxon"})
    a = ad.AnnData(
        X=np.ones((4, 2), dtype=float),
        obs=pd.DataFrame({"leiden": ["0", "0", "1", "1"]}, index=["c0", "c1", "c2", "c3"]),
        var=pd.DataFrame(index=["G0", "G1"]),
    )
    a.uns["rank_genes_groups"] = rgg_ttest
    recomputed = []

    def fake_rgg(matrix, **kw):
        recomputed.append(kw["method"])
        matrix.uns[kw["key_added"]] = rgg_wilcox

    monkeypatch.setitem(sys.modules, "scanpy", types.SimpleNamespace(tl=types.SimpleNamespace(rank_genes_groups=fake_rgg)))
    monkeypatch.setattr(safe_h5ad, "read_h5ad", lambda _p: a.copy())
    p = tmp_path / "cl.h5ad"
    p.write_text("stub", encoding="utf-8")

    out = de_mod.rna_de_per_cluster({"data_path": str(p), "output_dir": str(tmp_path)})
    assert recomputed == ["wilcoxon"]
    assert out["marker_source"] == "computed_aria_de"

    # 2. rna_pathway_per_cluster prioritizes positive markers for cluster-identity ORA
    submitted = {}
    monkeypatch.setattr(local_ora, "load_local_library", lambda name: ({"T_pathway": ["CD3D"]}, {"library": name}))

    def capture_ora(genes, gene_sets, bg, **kw):
        submitted.setdefault(tuple(genes), 0)
        submitted[tuple(genes)] += 1
        return []

    monkeypatch.setattr(local_ora, "run_ora", capture_ora)
    pw_mod.rna_pathway_per_cluster({
        "organism": "Homo sapiens",
        "output_dir": str(tmp_path),
        "de_genes_by_cluster": {
            "0": [
                {"gene": "MS4A1", "padj": 1e-50, "log2fc": -5.0},
                {"gene": "CD3D", "padj": 1e-4, "log2fc": 2.5},
            ],
            "Tcell::trt_vs_ctrl": [
                {"gene": "DOWN1", "padj": 1e-20, "log2fc": -3.0},
                {"gene": "UP1", "padj": 1e-10, "log2fc": 2.0},
            ],
        },
    })
    assert ("CD3D",) in submitted
    assert ("DOWN1", "UP1") in submitted


@needs_anndata
def test_d1_m9_da_fisher_fallback_method_and_pseudobulk_covariates(tmp_path, monkeypatch):
    import statsmodels.api as sm
    from aria.scripts.rna_diff_abundance import rna_diff_abundance

    a = ad.AnnData(
        X=np.ones((8, 2)),
        obs=pd.DataFrame({
            "cell": ["A", "B"] * 4,
            "condition": ["test"] * 4 + ["ref"] * 4,
            "donor": ["t1", "t1", "t2", "t2", "r1", "r1", "r2", "r2"],
        }, index=[f"c{i}" for i in range(8)]),
        var=pd.DataFrame(index=["G0", "G1"]),
    )
    h5 = tmp_path / "da_fail.h5ad"
    a.write_h5ad(h5)

    def broken_ols(**kw):
        raise ValueError("singular design")

    monkeypatch.setattr(sm, "OLS", broken_ols)
    res = rna_diff_abundance({
        "data_path": str(h5),
        "groupby": "cell",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons": [["test", "ref"]],
    })
    assert res["status"] == "success"
    assert res["method"] == "fisher_exact_fallback"
    assert res["model"]["primary"] == "fisher_exact_fallback"
    assert res["per_comparison"]["test_vs_ref"]["method"] == "fisher_exact_fallback"
    assert res["per_comparison"]["test_vs_ref"]["model"] == "fisher_exact_fallback"


@needs_agent
def test_d1_m9_pseudobulk_batch_cov_and_numerator_denominator_normalization(tmp_path):
    from aria.agents.scrna_agent import scRNAAgent

    norm = scRNAAgent._normalise_pseudobulk_comparisons
    assert norm([{"numerator": "treat", "denominator": "ctrl"}]) == [["treat", "ctrl"]]
    assert norm([{"numerator": 24, "denominator": 0}]) == [["24", "0"]]
    assert norm([[48, 0]]) == [["48", "0"]]

    captured = {}

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            captured[Path(script_path).name] = params
            if script_path.endswith("rna_diff_abundance.py"):
                return {"status": "success", "any_significant": False, "per_comparison": {}}
            return {"status": "success", "n_groups": 1, "per_group": {}}

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.env = _Env()
    agent._workspace = lambda exp_id, *parts: tmp_path
    agent._log_decision = lambda *a, **kw: None
    agent.publish_finding = lambda *a, **kw: None

    res = agent._run_pseudobulk(
        experiment_id="exp",
        current_h5ad=str(tmp_path / "in.h5ad"),
        exp_ctx={
            "organism": "Homo sapiens",
            "design": {
                "main_factor": "condition",
                "batch_covariate": "batch",
                "pseudobulk": {
                    "from_obs": True,
                    "condition_col": "condition",
                    "replicate_col": "donor",
                    "groupby_col": "cell_type",
                    "covariates": ["sex", "donor"],  # donor == replicate_col must be excluded
                    "comparisons": [{"numerator": "trt", "denominator": "ctrl"}],
                },
            },
        },
        intent={},
        annotation={"annotation_source": "input_obs", "label_col": "cell_type"},
    )
    assert res["status"] == "success"
    assert captured["rna_pseudobulk_de.py"]["covariates"] == ["sex", "batch"]
    assert captured["rna_diff_abundance.py"]["covariates"] == ["sex", "batch"]
