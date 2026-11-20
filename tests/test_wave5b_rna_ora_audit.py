"""Wave 5B regression & audit tests for Domain 1: Bulk RNA + scRNA + Local ORA.

Covers:
  - BI-01: pydeseq2 0.5.4 ref_level/lfc_shrink sign-flip & formulaic reference
  - CDX-01: local ORA BH denominator includes all non-empty in-background sets
  - BI-02 / D1-H3 / D1-H5 / D1-H14: BulkRNAAgent covariates, species/GTF,
    featureCounts headers/exact-before-substring matching, legacy list groups
  - D1-C2: scRNA generic sample-ID parent derivation, QC sid dedup & cache key
  - D1-C4: pseudobulk counts_adata.raw.X + lognorm recovery + error persistence
  - D1-C5: CellTypist never mutates original_adata.X or wipes original_adata.raw
  - BI-03: trajectory pseudotime_by_group finite filtering across all 3 sites
  - D1-H2: integration silhouette None/NaN guard in scrna/qc.py
  - BI-04 / BI-05: pseudobulk/adapter provenance forwarding & honest L-R wording
"""
from __future__ import annotations

import json
import math
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from aria.scripts.rna_bulk import deseq2 as bulk_deseq2
from aria.scripts import rna_pseudobulk_de as pb_mod
from aria.scripts import rna_narrative_adapter as adapter
from aria.utils import ora


try:
    import pydeseq2  # noqa: F401
    _HAS_PYDESEQ2 = True
except Exception:
    _HAS_PYDESEQ2 = False

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
    import litellm  # noqa: F401
    _HAS_LITELLM = True
except Exception:
    _HAS_LITELLM = False

needs_pydeseq2 = pytest.mark.skipif(not _HAS_PYDESEQ2, reason="pydeseq2 required")
needs_anndata = pytest.mark.skipif(not _HAS_ANNDATA, reason="anndata required")
needs_scanpy = pytest.mark.skipif(not _HAS_SCANPY, reason="scanpy required")
needs_agent = pytest.mark.skipif(not _HAS_LITELLM, reason="litellm required for agents")


# ── BI-01: formulaic reference level, _shrink_coeff, sign-flip guard ─────────

def test_bi01_build_design_formula_reference_and_quoting():
    assert bulk_deseq2._build_design_formula("condition", []) == "~ condition"
    assert bulk_deseq2._build_design_formula("condition", ["batch"]) == "~ batch + condition"
    assert (
        bulk_deseq2._build_design_formula("condition", [], ref_level="KO_ctrl")
        == "~ C(condition, Treatment(reference='KO_ctrl'))"
    )
    assert (
        bulk_deseq2._build_design_formula(
            "cell type", ["batch-id", "age score"], ref_level="Control"
        )
        == "~ `batch-id` + `age score` + C(`cell type`, Treatment(reference='Control'))"
    )
    assert (
        pb_mod._build_design_formula("time-point", ["donor"], ref_level="AA")
        == "~ donor + C(`time-point`, Treatment(reference='AA'))"
    )


@pytest.mark.parametrize("fn", [bulk_deseq2._shrink_coeff, pb_mod._shrink_coeff])
def test_bi01_shrink_coeff_exact_no_substring_collision(fn):
    def _dds(cols):
        return SimpleNamespace(varm={"LFC": pd.DataFrame(columns=cols)})

    # Substring traps MUST return None (never match the inverted ref coefficient!)
    assert fn(_dds(["Intercept", "condition[T.KO_ctrl]"]), "condition", "KO") is None
    assert fn(_dds(["Intercept", "condition[T.AA]"]), "condition", "A") is None
    assert fn(_dds(["Intercept", "condition[T.Control]"]), "condition", "C") is None
    assert fn(_dds(["Intercept", "condition[T.WT]"]), "condition", "T") is None
    assert fn(_dds(["Intercept", "batch[T.A]"]), "condition", "A") is None
    assert fn(_dds(["Intercept", "condition_extra[T.KO]"]), "condition", "KO") is None
    assert fn(
        _dds(["Intercept", "C(condition, Treatment(reference='KO'))[T.KO_ctrl]"]),
        "condition",
        "KO",
    ) is None

    # Exact plain & formulaic C(...) coefficients MUST match
    assert fn(_dds(["Intercept", "condition[T.KO]"]), "condition", "KO") == "condition[T.KO]"
    assert (
        fn(
            _dds(["Intercept", "C(condition, Treatment(reference='KO_ctrl'))[T.KO]"]),
            "condition",
            "KO",
        )
        == "C(condition, Treatment(reference='KO_ctrl'))[T.KO]"
    )
    assert (
        fn(
            _dds(["Intercept", "batch[T.A]", "C(condition, Treatment(reference='AA'))[T.A]"]),
            "condition",
            "A",
        )
        == "C(condition, Treatment(reference='AA'))[T.A]"
    )
    assert (
        fn(
            _dds(["Intercept", "C(`cell type`, Treatment(reference='Control'))[T.Case]"]),
            "cell type",
            "Case",
        )
        == "C(`cell type`, Treatment(reference='Control'))[T.Case]"
    )


def test_bi01_sign_flip_guard_and_nonfinite_lfc_excluded(monkeypatch):
    genes = [f"G{i}" for i in range(6)]
    samples = [f"s{i}" for i in range(6)]
    counts = pd.DataFrame(np.full((6, 6), 50), index=genes, columns=samples)
    metadata = pd.DataFrame(
        {"condition": ["KO_ctrl"] * 3 + ["KO"] * 3}, index=samples
    )

    class _FlippingDDS:
        def __init__(self, **kw):
            self.var = pd.DataFrame(index=genes)
            self.varm = {
                "LFC": pd.DataFrame(
                    columns=["Intercept", "C(condition, Treatment(reference='KO_ctrl'))[T.KO]"]
                )
            }

        def deseq2(self):
            pass

    class _FlippingStats:
        def __init__(self, dds, **kw):
            self.results_df = pd.DataFrame(
                {
                    "log2FoldChange": [2.5, -2.5, 3.0, -3.0, np.nan, np.inf],
                    "lfcSE": [0.2] * 6,
                    "padj": [0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
                    "pvalue": [1e-4] * 6,
                    "baseMean": [50.0] * 6,
                },
                index=genes,
            )

        def summary(self):
            pass

        def lfc_shrink(self, coeff):
            # Simulate an inverted coefficient flipping every sign
            self.results_df["log2FoldChange"] = -self.results_df["log2FoldChange"] * 0.9

    dds_mod = types.ModuleType("pydeseq2.dds")
    dds_mod.DeseqDataSet = _FlippingDDS
    ds_mod = types.ModuleType("pydeseq2.ds")
    ds_mod.DeseqStats = _FlippingStats
    monkeypatch.setitem(sys.modules, "pydeseq2", types.ModuleType("pydeseq2"))
    monkeypatch.setitem(sys.modules, "pydeseq2.dds", dds_mod)
    monkeypatch.setitem(sys.modules, "pydeseq2.ds", ds_mod)

    res, warns = bulk_deseq2._run_deseq2(
        counts, metadata, "condition", "KO", "KO_ctrl",
        padj_thr=0.05, lfc_thr=0.5, lfc_shrink=True, min_replicates_per_condition=3,
    )
    assert res["status"] == "success"
    assert res["lfc_shrinkage"]["applied"] is False
    assert "sign-flip" in (res["lfc_shrinkage"]["reason"] or "")
    assert any("sign" in w.lower() for w in warns)
    # Reverted to raw MLE signs AND excluded NaN/inf genes G4, G5 from sig/up/down!
    assert res["up_genes"] == ["G0", "G2"]
    assert res["down_genes"] == ["G1", "G3"]
    assert res["n_sig"] == 4


_BI01_LIVE_CASES = [
    ("KO_vs_KO_ctrl", "condition", "KO", "KO_ctrl", {}),
    ("A_vs_AA", "condition", "A", "AA", {}),
    ("Case_vs_Control", "condition", "Case", "Control", {}),
    ("WT_vs_KO", "condition", "WT", "KO", {}),
    ("space_col", "cell type", "Case", "Control", {}),
    ("hyphen_col", "time-point", "A", "AA", {}),
    (
        "covs_cat_cont",
        "treatment_group",
        "KO",
        "KO_ctrl",
        {
            "batch-id": ["b1", "b2", "b1", "b2", "b1", "b2", "b1", "b2"],
            "age score": [1.0, 2.0, 3.0, 4.0, 1.5, 2.5, 3.5, 4.5],
        },
    ),
]


def _make_bi01_counts(seed: int = 101):
    rng = np.random.default_rng(seed)
    n_genes = 60
    base = rng.negative_binomial(60, 0.4, (n_genes, 8))
    # First 10 genes strongly UP in test (cols 4..7) vs ref (cols 0..3)
    base[:10, 4:] = base[:10, 4:] * 6
    # Next 10 genes strongly DOWN in test vs ref
    base[10:20, :4] = base[10:20, :4] * 6
    genes = [f"G{i:02d}" for i in range(n_genes)]
    samples = [f"s{i}" for i in range(8)]
    return pd.DataFrame(base, index=genes, columns=samples)


@needs_pydeseq2
@pytest.mark.parametrize("label,col,test_lvl,ref_lvl,extra_covs", _BI01_LIVE_CASES)
def test_bi01_live_bulk_deseq2_no_sign_flip(label, col, test_lvl, ref_lvl, extra_covs):
    counts = _make_bi01_counts()
    meta_dict = {col: [ref_lvl] * 4 + [test_lvl] * 4, **extra_covs}
    metadata = pd.DataFrame(meta_dict, index=counts.columns)
    res, _ = bulk_deseq2._run_deseq2(
        counts,
        metadata,
        design_factor=col,
        numerator=test_lvl,
        denominator=ref_lvl,
        padj_thr=0.05,
        lfc_thr=0.5,
        covariates=list(extra_covs.keys()),
        lfc_shrink=True,
        min_replicates_per_condition=3,
    )
    assert res["status"] == "success", res
    assert res["lfc_shrinkage"]["applied"] is True, res["lfc_shrinkage"]
    df = res["results"]
    # All 10 planted UP genes must have positive shrunken & raw LFC
    assert (df.loc[[f"G{i:02d}" for i in range(10)], "log2FoldChange"] > 1.0).all()
    assert (df.loc[[f"G{i:02d}" for i in range(10)], "log2FoldChange_raw"] > 1.0).all()
    # All 10 planted DOWN genes must have negative shrunken & raw LFC
    assert (df.loc[[f"G{i:02d}" for i in range(10, 20)], "log2FoldChange"] < -1.0).all()
    assert (df.loc[[f"G{i:02d}" for i in range(10, 20)], "log2FoldChange_raw"] < -1.0).all()
    assert set(f"G{i:02d}" for i in range(10)).issubset(set(res["up_genes"]))
    assert set(f"G{i:02d}" for i in range(10, 20)).issubset(set(res["down_genes"]))


@needs_pydeseq2
@needs_anndata
@pytest.mark.parametrize("label,col,test_lvl,ref_lvl,extra_covs", _BI01_LIVE_CASES)
def test_bi01_live_pseudobulk_deseq2_no_sign_flip(
    tmp_path, label, col, test_lvl, ref_lvl, extra_covs
):
    from scipy import sparse

    counts = _make_bi01_counts(seed=202)
    rows, obs_rows = [], []
    for j, sid in enumerate(counts.columns):
        vec = counts[sid].to_numpy()
        half = vec // 2
        for k, part in enumerate((half, vec - half)):
            rows.append(part)
            row = {
                "cell_type": "TypeA",
                col: ref_lvl if j < 4 else test_lvl,
                "donor": f"d{j}",
            }
            for ck, cv in extra_covs.items():
                row[ck] = cv[j]
            obs_rows.append(row)
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.asarray(rows)),
        obs=pd.DataFrame(obs_rows, index=[f"c{i}" for i in range(len(rows))]),
        var=pd.DataFrame(index=counts.index),
    )
    h5 = tmp_path / f"{label}.h5ad"
    adata.write_h5ad(h5)

    res = pb_mod.rna_pseudobulk_de({
        "data_path": str(h5),
        "groupby": "cell_type",
        "condition_col": col,
        "replicate_col": "donor",
        "comparisons": [[test_lvl, ref_lvl]],
        "covariates": list(extra_covs.keys()),
        "min_cells_per_pseudosample": 2,
        "min_replicates_per_condition": 3,
        "use_raw": False,
        "lfc_shrink": True,
        "lfc_min": 0.5,
        "padj_max": 0.05,
        "output_dir": str(tmp_path / f"out_{label}"),
    })
    assert res["status"] == "success", res
    comp = res["per_group"]["TypeA"]["per_comparison"][f"{test_lvl}_vs_{ref_lvl}"]
    assert comp["status"] == "success", comp
    assert comp["lfc_shrinkage"]["applied"] is True, comp["lfc_shrinkage"]
    by_gene = {r["gene"]: r for r in comp["all_sig"]}
    for i in range(10):
        g = f"G{i:02d}"
        assert g in by_gene and by_gene[g]["log2fc"] > 1.0 and by_gene[g]["log2fc_raw"] > 1.0
    for i in range(10, 20):
        g = f"G{i:02d}"
        assert g in by_gene and by_gene[g]["log2fc"] < -1.0 and by_gene[g]["log2fc_raw"] < -1.0


# ── CDX-01: Local ORA BH denominator includes zero-overlap sets ──────────────

def test_cdx01_ora_bh_includes_zero_overlap_in_background_sets():
    bg = [f"G{i}" for i in range(2000)]
    query = [f"G{i}" for i in range(5)]
    # 1 hit set + 99 zero-overlap sets IN background + 50 sets OUTSIDE background (n=0)
    sets = {"HIT": [f"G{i}" for i in range(10)]}
    for j in range(99):
        sets[f"ZERO_{j}"] = [f"G{100 + j * 5}", f"G{101 + j * 5}"]
    for j in range(50):
        sets[f"OUTSIDE_{j}"] = [f"NOT_IN_BG_{j}"]

    res = ora.run_ora(query, sets, bg, padj_max=1.01, min_overlap=1, top=20)
    assert [r["term"] for r in res] == ["HIT"]
    hit = res[0]
    # Denominator m MUST be 100 (1 hit + 99 zero-overlap in-bg), excluding the 50 n==0 sets
    assert hit["padj"] == pytest.approx(min(1.0, hit["pvalue"] * 100.0), rel=1e-9)

    # Filtering by min_overlap=2 must keep the exact same m=100 denominator
    res_min2 = ora.run_ora(query, sets, bg, padj_max=1.01, min_overlap=2, top=20)
    assert res_min2[0]["padj"] == pytest.approx(hit["padj"], rel=1e-12)


# ── BI-02, D1-H3, D1-H5, D1-H14: BulkRNAAgent ────────────────────────────────

@needs_agent
def test_bi02_write_design_metadata_resolves_all_covariate_sources(tmp_path):
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    from aria.scripts.rna_bulk_de import _aggregate_technical_replicates

    group_labels = {
        "ctrl_u1_r1.bam": "ctrl",
        "ctrl_u1_r2.bam": "ctrl",
        "ctrl_u2_r1.bam": "ctrl",
        "ctrl_u2_r2.bam": "ctrl",
        "trt_u1_r1.bam": "trt",
        "trt_u1_r2.bam": "trt",
        "trt_u2_r1.bam": "trt",
        "trt_u2_r2.bam": "trt",
    }
    replicate_units = {
        "ctrl_u1_r1.bam": "ctrl_u1",
        "ctrl_u1_r2.bam": "ctrl_u1",
        "ctrl_u2_r1.bam": "ctrl_u2",
        "ctrl_u2_r2.bam": "ctrl_u2",
        "trt_u1_r1.bam": "trt_u1",
        "trt_u1_r2.bam": "trt_u1",
        "trt_u2_r1.bam": "trt_u2",
        "trt_u2_r2.bam": "trt_u2",
    }

    # Source 1: metadata_file TSV (without .bam suffix)
    meta_tsv = tmp_path / "user_meta.tsv"
    pd.DataFrame({
        "sample": [s.replace(".bam", "") for s in group_labels],
        "condition": list(group_labels.values()),
        "batch": ["b1", "b1", "b2", "b2", "b1", "b1", "b2", "b2"],
    }).to_csv(meta_tsv, sep="\t", index=False)

    # Source 2: inferred_design.sample_sheet (list of dicts) providing sex
    sample_sheet = [
        {"sample": s.replace(".bam", ""), "sex": "F" if "u1" in s else "M"}
        for s in group_labels
    ]
    # Source 3: sample_metadata (dict) providing age
    sample_metadata = {
        s.replace(".bam", ""): {"age": 30 if "u1" in s else 45}
        for s in group_labels
    }

    out_path = BulkRNAAgent._write_design_metadata(
        group_labels,
        "condition",
        tmp_path / "out",
        replicate_units=replicate_units,
        covariates=["batch", "sex", "age"],
        exp_ctx={
            "metadata_file": str(meta_tsv),
            "sample_metadata": sample_metadata,
            "inferred_design": {"sample_sheet": sample_sheet},
        },
        design={"batch_covariate": "batch", "covariates": ["sex", "age"]},
    )
    written = pd.read_csv(out_path, sep="\t", index_col=0)
    assert list(written.columns) == ["condition", "biological_unit", "batch", "sex", "age"]
    assert written["batch"].tolist() == ["b1", "b1", "b2", "b2", "b1", "b1", "b2", "b2"]
    assert written["sex"].tolist() == ["F", "F", "M", "M", "F", "F", "M", "M"]

    # Must not crash _aggregate_technical_replicates when covariates are present!
    counts = pd.DataFrame(
        {s: [10, 20] for s in group_labels}, index=["G1", "G2"]
    )
    agg_counts, agg_meta, _ = _aggregate_technical_replicates(
        counts, written, design_factor="condition",
        unit_col="biological_unit", covariates=["batch", "sex"],
    )
    assert agg_counts.shape == (2, 4)
    assert list(agg_meta.columns) == ["condition", "batch", "sex"]


@needs_agent
def test_d1_h3_species_genome_no_conflicting_defaults_and_gtf_forward(tmp_path):
    from aria.agents.bulk_rna_agent import BulkRNAAgent

    counts_path = tmp_path / "counts.tsv"
    counts_path.write_text(
        "gene\tctrl1\tctrl2\tctrl3\ttrt1\ttrt2\ttrt3\nG1\t10\t12\t11\t40\t42\t41\n"
    )
    design = {
        "groups": {"ctrl": ["ctrl1", "ctrl2", "ctrl3"], "trt": ["trt1", "trt2", "trt3"]},
        "main_factor": "condition",
        "contrasts": [{"numerator": "trt", "denominator": "ctrl"}],
    }

    def _run_with_ctx(exp_ctx):
        captured = {}
        agent = BulkRNAAgent.__new__(BulkRNAAgent)
        agent.env = SimpleNamespace(
            run_in_stack=lambda **kw: (captured.update(kw["params"]) or {"status": "success", "contrasts": []})
        )
        agent.memory = SimpleNamespace(store_decision=lambda **kw: None)
        agent.publish_status = lambda *a, **kw: None
        agent.publish_finding = lambda *a, **kw: None
        agent._publish_findings = lambda *a, **kw: None
        agent._record_methodology_decisions = lambda *a, **kw: None
        res = agent.run("exp-h3", {"exp_context": exp_ctx, "biological_intent": {}})
        assert res["status"] == "done"
        return captured

    # Only organism="Mus musculus" -> MUST NOT fabricate genome="hg38"
    p1 = _run_with_ctx({
        "modalities": {"bulk_RNA": [str(counts_path)]},
        "organism": "Mus musculus",
        "design": design,
    })
    assert p1.get("organism") == "Mus musculus"
    assert p1.get("genome") in (None, "")

    # Only genome="mm39" + genome_config.gtf -> MUST NOT fabricate organism="Homo sapiens"
    p2 = _run_with_ctx({
        "modalities": {"bulk_RNA": [str(counts_path)]},
        "genome": "mm39",
        "genome_config": {"gtf": "/ref/mm39.gtf"},
        "design": design,
    })
    assert p2.get("genome") == "mm39"
    assert p2.get("organism") in (None, "")
    assert p2.get("gtf") == "/ref/mm39.gtf"

    # Fallback to design["organism"] / design["genome"]
    p3 = _run_with_ctx({
        "modalities": {"bulk_RNA": [str(counts_path)]},
        "design": {**design, "organism": "Mus musculus", "genome": "mm10"},
    })
    assert p3.get("organism") == "Mus musculus"
    assert p3.get("genome") == "mm10"


@needs_agent
def test_d1_h5_featurecounts_comments_exact_match_and_same_group_aliases(tmp_path):
    from aria.agents.bulk_rna_agent import BulkRNAAgent

    fc_path = tmp_path / "featureCounts.txt"
    fc_path.write_text(
        "# Program:featureCounts v2.0.6; Command:featureCounts ...\n"
        "Geneid\tChr\tStart\tEnd\tStrand\tLength\tS1\tS2\tS10\tS11\n"
        "G1\tchr1\t1\t100\t+\t100\t10\t11\t50\t52\n",
        encoding="utf-8",
    )
    assert BulkRNAAgent._read_sample_names(str(fc_path)) == ["S1", "S2", "S10", "S11"]

    agent = BulkRNAAgent.__new__(BulkRNAAgent)
    design = {
        "groups": {"ctrl": ["S1", "S2"], "treat": ["S10", "S11"]},
        "main_factor": "condition",
        "sample_aliases": {
            "S1": ["S1", "s1"],
            "S10": ["S10", "s10"],
        },
        "contrasts": [{"numerator": "treat", "denominator": "ctrl"}],
    }
    samples, group_labels, factor, contrasts = agent._apply_design(
        design, [str(fc_path)], "exp-h5"
    )
    assert samples == ["S1", "S2", "S10", "S11"]
    assert group_labels == {"S1": "ctrl", "S2": "ctrl", "S10": "treat", "S11": "treat"}
    assert factor == "condition"
    assert len(contrasts) == 1


@needs_agent
def test_d1_h14_legacy_fallback_list_valued_group_labels_no_typeerror(tmp_path, monkeypatch):
    from aria.agents.bulk_rna_agent import BulkRNAAgent

    monkeypatch.setenv("ARIA_ALLOW_FILENAME_FALLBACK", "1")
    counts_path = tmp_path / "counts.tsv"
    counts_path.write_text(
        "gene\tctrl_1\tctrl_2\ttreat_1\ttreat_2\n"
        "G1\t10\t12\t40\t44\n",
        encoding="utf-8",
    )
    statuses = []
    agent = BulkRNAAgent.__new__(BulkRNAAgent)
    agent.publish_status = lambda exp_id, msg, prog: statuses.append(msg)
    agent.publish_finding = lambda *a, **kw: None
    agent.publish_escalation = lambda *a, **kw: None

    res = agent.run(
        "exp-h14",
        {
            "exp_context": {"modalities": {"bulk_RNA": [str(counts_path)]}},
            "biological_intent": {},
        },
    )
    assert res["status"] == "failed"
    assert res["reason"] == "explicit_contrast_required"
    assert any("Detected groups: ['ctrl', 'treat']" in s for s in statuses)


# ── D1-C2: scRNA generic sample_id, QC sid dedup, rna_qc cache key ──────────

@needs_agent
def test_d1_c2_sample_id_from_path_and_qc_dedup(tmp_path):
    from aria.agents.scrna_agent import scRNAAgent

    assert scRNAAgent._sample_id_from_path("/data/donor_A/outs/filtered_feature_bc_matrix.h5") == "donor_A"
    assert scRNAAgent._sample_id_from_path("/data/donor_B/outs/filtered_feature_bc_matrix") == "donor_B"
    assert scRNAAgent._sample_id_from_path("/data/sample_X/raw_feature_bc_matrix.h5") == "sample_X"
    assert scRNAAgent._sample_id_from_path("/data/sample_Y/matrix.h5ad") == "sample_Y"
    assert scRNAAgent._sample_id_from_path("/data/sample_Z/counts.h5ad") == "sample_Z"
    assert scRNAAgent._sample_id_from_path("/data/sample_W/adata.h5ad") == "sample_W"
    assert scRNAAgent._sample_id_from_path("/data/sample01_raw_feature_bc_matrix.h5") == "sample01"

    # Two different paths sharing the same derived parent name must be deduplicated in _run_qc
    seen_sids = []

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            if script_path.endswith("rna_qc.py"):
                sid = params["sample_id"]
                seen_sids.append(sid)
                return {
                    "status": "success",
                    "output_path": str(tmp_path / f"qc_filtered_{sid}.h5ad"),
                    "n_cells_before": 100,
                    "n_cells_after": 90,
                    "pct_removed": 10.0,
                }
            return {
                "status": "success",
                "output_path": str(tmp_path / "concatenated.h5ad"),
                "n_cells_total": 180,
                "n_genes_shared": 500,
                "batch_col": "batch",
            }

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.env = _Env()
    agent._workspace = lambda exp_id, *parts: tmp_path
    agent.publish_finding = lambda *a, **kw: None
    qc_res = agent._run_qc(
        "exp-c2",
        ["/run1/dup/matrix.h5ad", "/run2/dup/matrix.h5ad", "/run3/dup/matrix.h5ad"],
        {"organism": "Homo sapiens"},
        {},
    )
    assert qc_res["status"] == "success"
    assert seen_sids == ["dup", "dup_2", "dup_3"]


def test_d1_c2_rna_qc_cache_params_distinguishes_inputs(tmp_path):
    from aria.scripts.rna_qc import _cache_params

    f1 = tmp_path / "s1" / "filtered_feature_bc_matrix.h5ad"
    f2 = tmp_path / "s2" / "filtered_feature_bc_matrix.h5ad"
    f1.parent.mkdir()
    f2.parent.mkdir()
    f1.write_bytes(b"sample-1")
    f2.write_bytes(b"sample-2-longer")

    c1 = _cache_params({"data_path": str(f1), "sample_id": "same"})
    c2 = _cache_params({"data_path": str(f2), "sample_id": "same"})
    assert c1["data_path"] != c2["data_path"]
    assert c1 != c2


@needs_scanpy
@needs_anndata
def test_d1_c2_rna_qc_live_no_cross_sample_cache_hit(tmp_path):
    from aria.scripts.rna_qc import rna_qc

    p1 = tmp_path / "s1.h5ad"
    p2 = tmp_path / "s2.h5ad"
    out_dir = tmp_path / "qc_out"

    a1 = ad.AnnData(
        X=np.ones((6, 3), dtype=float),
        obs=pd.DataFrame(
            {"total_counts": [1000.0] * 6, "n_genes_by_counts": [300.0] * 6, "pct_counts_mt": [2.0] * 6},
            index=[f"s1_c{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    a2 = ad.AnnData(
        X=np.ones((10, 3), dtype=float),
        obs=pd.DataFrame(
            {"total_counts": [1000.0] * 10, "n_genes_by_counts": [300.0] * 10, "pct_counts_mt": [2.0] * 10},
            index=[f"s2_c{i}" for i in range(10)],
        ),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    a1.write_h5ad(p1)
    a2.write_h5ad(p2)

    r1 = rna_qc({"data_path": str(p1), "sample_id": "shared", "run_scrublet": False, "output_dir": str(out_dir)})
    r2 = rna_qc({"data_path": str(p2), "sample_id": "shared", "run_scrublet": False, "output_dir": str(out_dir)})
    assert r1["status"] == "success" and r1["n_cells_after"] == 6
    assert r2["status"] == "success" and r2["n_cells_after"] == 10
    assert not r2.get("resumed")


# ── D1-C4: counts_adata.raw.X fallback, lognorm recovery, agent error persist ─

@needs_anndata
def test_d1_c4_pseudobulk_counts_adata_raw_and_lognorm_recovery(monkeypatch, tmp_path):
    class _FakeDDS:
        def __init__(self, counts, metadata, **kw):
            self.counts = counts

        def deseq2(self):
            pass

    class _FakeStats:
        def __init__(self, dds, **kw):
            self.results_df = pd.DataFrame(
                {
                    "pvalue": [0.001, 0.5, 0.5, 0.5],
                    "padj": [0.004, 0.8, 0.8, 0.8],
                    "log2FoldChange": [2.0, 0.1, -0.1, 0.0],
                    "baseMean": [50.0] * 4,
                },
                index=dds.counts.columns,
            )

        def summary(self):
            pass

    monkeypatch.setitem(sys.modules, "pydeseq2.dds", SimpleNamespace(DeseqDataSet=_FakeDDS))
    monkeypatch.setitem(sys.modules, "pydeseq2.ds", SimpleNamespace(DeseqStats=_FakeStats))
    monkeypatch.setattr(pb_mod, "compute_counts_lobo_jackknife", lambda *a, **kw: {"status": "skipped"})

    rng = np.random.default_rng(7)
    n_cells = 12
    raw_counts = rng.integers(5, 80, size=(n_cells, 4)).astype(float)
    lib_sizes = raw_counts.sum(axis=1)
    lognorm = np.log1p(raw_counts / lib_sizes[:, None] * 10000.0)
    obs = pd.DataFrame(
        {
            "cell_type": ["A"] * n_cells,
            "condition": ["test"] * 6 + ["ref"] * 6,
            "donor": [f"d{i // 2}" for i in range(n_cells)],
            "nCount_RNA": lib_sizes,
        },
        index=[f"c{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=["G0", "G1", "G2", "G3"])

    # Labelled data_path (scaled X)
    labelled = ad.AnnData(X=np.full((n_cells, 4), -0.5), obs=obs, var=var)
    labelled_path = tmp_path / "labelled.h5ad"
    labelled.write_h5ad(labelled_path)

    # Case 1: counts_data_path has scaled X but raw integer counts in .raw.X
    c_raw = ad.AnnData(X=np.full((n_cells, 4), -0.2), obs=obs, var=var)
    c_raw.raw = ad.AnnData(X=raw_counts, obs=obs, var=var)
    c_raw_path = tmp_path / "counts_with_raw.h5ad"
    c_raw.write_h5ad(c_raw_path)

    res1 = pb_mod.rna_pseudobulk_de({
        "data_path": str(labelled_path),
        "counts_data_path": str(c_raw_path),
        "groupby": "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons": [["test", "ref"]],
        "min_cells_per_pseudosample": 2,
        "min_replicates_per_condition": 3,
        "lfc_shrink": False,
        "output_dir": str(tmp_path / "out1"),
    })
    assert res1["status"] == "success", res1
    assert res1["count_source"] == "raw_counts"
    assert res1["lognorm_recovered"] is False

    # Case 2: counts_data_path has log-normalized X (no integer counts), allow_lognorm_recovery=True
    c_log = ad.AnnData(X=lognorm, obs=obs, var=var)
    c_log_path = tmp_path / "counts_lognorm.h5ad"
    c_log.write_h5ad(c_log_path)

    res2 = pb_mod.rna_pseudobulk_de({
        "data_path": str(labelled_path),
        "counts_data_path": str(c_log_path),
        "groupby": "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons": [["test", "ref"]],
        "min_cells_per_pseudosample": 2,
        "min_replicates_per_condition": 3,
        "allow_lognorm_recovery": True,
        "lfc_shrink": False,
        "output_dir": str(tmp_path / "out2"),
    })
    assert res2["status"] == "success", res2
    assert res2["count_source"] == "recovered_from_lognorm"
    assert res2["lognorm_recovered"] is True


@needs_agent
@pytest.mark.parametrize("pb_out", [
    {"status": "error", "reason": "pseudobulk_de_failed: NonIntegerCounts", "details": "bad"},
    {"status": "skipped", "reason": "explicit_comparison_required"},
])
def test_d1_c4_scrna_agent_persists_error_and_skipped_pseudobulk(tmp_path, pb_out):
    from aria.agents.scrna_agent import scRNAAgent

    h5 = tmp_path / "qc.h5ad"
    h5.write_text("stub")
    agent = scRNAAgent.__new__(scRNAAgent)
    agent.publish_status = lambda *a, **kw: None
    agent.publish_finding = lambda *a, **kw: None
    agent._log_decision = lambda *a, **kw: None
    agent._prepare_focused_h5ads = lambda *a, **kw: {"status": "skipped"}
    agent._run_qc = lambda *a, **kw: {"status": "success", "output_path": str(h5)}
    agent._resolve_batch_column = lambda *a, **kw: None
    agent._predefined_celltype_col = lambda *a, **kw: None
    agent._run_clustering = lambda *a, **kw: (
        {"status": "success", "output_path": str(h5), "n_clusters": 2, "groupby": "leiden", "top_markers": {}, "cluster_sizes": {}},
        SimpleNamespace(chosen_value=0.5, justification="ok"),
    )
    agent._annotate_cell_types = lambda *a, **kw: {"cell_types": {}, "label_col": "leiden"}
    agent._differential_expression = lambda *a, **kw: {"status": "success", "n_significant_genes": 0}
    agent._needs_pseudobulk = lambda *a, **kw: True
    agent._run_pseudobulk = lambda *a, **kw: pb_out
    agent._needs_trajectory = lambda *a, **kw: False
    agent._needs_cell_communication = lambda *a, **kw: False

    out = agent.run(
        "exp-c4",
        {"exp_context": {"modalities": {"scRNA": [str(h5)]}}, "biological_intent": {}},
    )
    assert out["status"] == "done"
    assert out["findings"]["pseudobulk_de"] == pb_out


# ── D1-C5: CellTypist preserves original_adata.X and .raw ────────────────────

@needs_scanpy
@needs_anndata
def test_d1_c5_celltypist_never_mutates_original_raw_or_x(monkeypatch, tmp_path):
    from aria.scripts import rna_celltypist as ct_mod

    raw_x = np.array([[10, 20, 30, 40], [5, 15, 25, 35], [12, 8, 16, 24], [7, 9, 11, 13]], dtype=float)
    a = ad.AnnData(
        X=raw_x.copy(),
        obs=pd.DataFrame({"leiden": ["0", "0", "1", "1"]}, index=[f"c{i}" for i in range(4)]),
        var=pd.DataFrame(index=["G0", "G1", "G2", "G3"]),
    )
    # .raw stores raw integer counts (NOT lognorm_10k), so adata stayed original_adata!
    a.raw = a.copy()
    in_path = tmp_path / "in.h5ad"
    a.write_h5ad(in_path)

    def _annotate(matrix, **kw):
        assert matrix.raw is None
        return SimpleNamespace(
            predicted_labels=pd.DataFrame({"predicted_labels": ["T"] * 4}, index=matrix.obs_names),
            probability_matrix=pd.DataFrame({"T": [0.9] * 4}, index=matrix.obs_names),
        )

    monkeypatch.setitem(sys.modules, "celltypist", SimpleNamespace(annotate=_annotate, models=SimpleNamespace()))
    monkeypatch.setattr(ct_mod, "_celltypist_model_cached", lambda m: True)

    res = ct_mod.rna_celltypist({
        "data_path": str(in_path),
        "model": "explicit.pkl",
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "success", res
    saved = ad.read_h5ad(res["output_path"])
    assert saved.raw is not None, "original_adata.raw was wiped in place!"
    np.testing.assert_array_equal(saved.X, raw_x)
    np.testing.assert_array_equal(saved.raw.X, raw_x)
    assert list(saved.obs["cell_type_celltypist"]) == ["T"] * 4


# ── BI-03: Trajectory pseudotime_by_group finite filtering ─────────────────────

@needs_agent
def test_bi03_trajectory_pseudotime_none_nan_inf_bool_safe():
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator
    from aria.agents.narrative.scrna.text import (
        _describe_trajectory_context,
        summarize_scrna_text,
    )

    findings = {
        "trajectory": {
            "status": "success",
            "paga": {"n_connections": 3, "n_strong": 1, "max_connectivity": 0.2},
            "pseudotime": {
                "computed": True,
                "root_used": "stem",
                "pseudotime_by_group": {
                    "mature": 0.8,
                    "stem": 0.1,
                    "disconnected_A": None,
                    "disconnected_B": float("nan"),
                    "inf_group": float("inf"),
                    "bool_group": True,
                },
            },
            "velocity": {"computed": False},
        }
    }
    blocks = ScrnaNarrator()._trajectory_blocks(findings)
    assert len(blocks) == 1
    ev_map = {e.label: e.value for e in blocks[0].evidence}
    assert ev_map["DPT ordered groups"] == "stem -> mature"

    desc = _describe_trajectory_context(findings)
    assert any("from stem to mature" in line for line in desc)

    summary = summarize_scrna_text(findings)
    assert "stem → mature" in summary

    # All-None pseudotime_by_group must also not crash
    all_none = {
        "trajectory": {
            "status": "success",
            "paga": {"n_connections": 1, "n_strong": 0, "max_connectivity": 0.01},
            "pseudotime": {
                "computed": True,
                "root_used": "auto",
                "pseudotime_by_group": {"A": None, "B": float("nan")},
            },
            "velocity": {"computed": False},
        }
    }
    assert len(ScrnaNarrator()._trajectory_blocks(all_none)) == 1
    assert isinstance(_describe_trajectory_context(all_none), list)
    assert isinstance(summarize_scrna_text(all_none), str)


# ── D1-H2: Integration silhouette None/NaN guard in scrna/qc.py ──────────────

@needs_agent
@pytest.mark.parametrize("bad", [None, float("nan"), float("inf"), True])
def test_d1_h2_integration_nonfinite_silhouette_no_crash(tmp_path, bad):
    from aria.agents.scrna_agent import scRNAAgent
    from aria.bus.message_bus import Confidence

    published = []
    agent = scRNAAgent.__new__(scRNAAgent)
    agent._workspace = lambda exp_id, *parts: tmp_path
    agent.env = SimpleNamespace(
        run_in_stack=lambda **kw: {
            "status": "success",
            "n_batches": 2,
            "silhouette_before": bad,
            "silhouette_after": bad,
            "batch_correction_delta": bad,
            "output_path": str(tmp_path / "integrated.h5ad"),
        }
    )
    agent.publish_finding = lambda exp_id, payload, conf: published.append((payload, conf))

    res = agent._run_integration("exp-h2", str(tmp_path / "in.h5ad"), "batch")
    assert res["status"] == "success"
    assert len(published) == 1
    assert published[0][1] == Confidence.MEDIUM
    assert "unavailable" in published[0][0]["summary"]


# ── BI-04 & BI-05: Pseudobulk/adapter provenance forwarding & L-R wording ────

@needs_agent
def test_bi04_scrna_pseudobulk_forwards_all_provenance_keys(tmp_path):
    from aria.agents.scrna_agent import scRNAAgent

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            if script_path.endswith("rna_diff_abundance.py"):
                return {
                    "status": "success",
                    "method": "donor_clr_ols_hc3",
                    "model": {"primary": "donor_clr_ols_hc3", "paired_donor_fixed_effects": True},
                    "groupby": "cell_type",
                    "condition_col": "condition",
                    "replicate_col": "donor",
                    "covariates": [],
                    "paired_design": True,
                    "significance_alpha": 0.10,
                    "any_significant": False,
                    "n_replicates_per_group": {"test": 3, "ref": 3},
                    "per_comparison": {},
                    "output_path": str(tmp_path / "da.tsv"),
                    "warnings": [],
                }
            if script_path.endswith("rna_pseudobulk_de.py"):
                return {
                    "status": "success",
                    "groupby": "cell_type",
                    "condition_col": "condition",
                    "replicate_col": "donor",
                    "covariates": [],
                    "thresholds": {"padj_max": 0.05, "lfc_min": 0.5},
                    "multiple_testing": {"fdr_strategy": "per_cluster"},
                    "lfc_shrinkage": {"requested": True, "method": "apeGLM"},
                    "power": {"alpha_nominal": 0.05},
                    "robustness_multiverse": {"method": "per-block"},
                    "background_genes": ["G1", "G2"],
                    "background_size": 2,
                    "background_source": "dataset_expressed_genes",
                    "count_source": "raw_counts",
                    "lognorm_recovered": False,
                    "paired_design": True,
                    "auto_paired_donor_covariate": True,
                    "n_groups": 1,
                    "per_group": {
                        "T": {
                            "background_genes": ["G1", "G2"],
                            "per_comparison": {
                                "treat_vs_ctrl": {
                                    "status": "success",
                                    "n_significant": 1,
                                    "all_sig": [{"gene": "G1", "log2fc": 2.0, "padj": 0.001}],
                                }
                            },
                        }
                    },
                    "output_csv": str(tmp_path / "pb.csv"),
                    "full_results_csv": str(tmp_path / "pb_all.csv"),
                }
            if script_path.endswith("rna_pathway_per_cluster.py"):
                return {
                    "status": "success",
                    "organism": "Homo sapiens",
                    "databases": {"GO_BP": "GO_Biological_Process_2021"},
                    "ora_method": "local_hypergeometric",
                    "gene_set_versions": {"GO_BP": {"library": "GO_Biological_Process_2021", "release": "2021"}},
                    "databases_skipped": ["KEGG"],
                    "background_size": None,
                    "background_source": "per_cluster_expressed_genes",
                    "per_cluster": {"T::treat_vs_ctrl": {"n_significant": 1, "results": {}}},
                    "output_csv": str(tmp_path / "pw.csv"),
                }
            raise AssertionError(script_path)

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.env = _Env()
    agent._workspace = lambda exp_id, *parts: tmp_path
    agent._log_decision = lambda *a, **kw: None
    agent.publish_finding = lambda *a, **kw: None
    res = agent._run_pseudobulk(
        "exp-bi04",
        str(tmp_path / "in.h5ad"),
        {
            "organism": "Homo sapiens",
            "design": {
                "groups": {"ctrl": ["s1", "s2", "s3"], "treat": ["s4", "s5", "s6"]},
                "main_factor": "condition",
                "pseudobulk": {
                    "from_obs": True,
                    "condition_col": "condition",
                    "replicate_col": "donor",
                    "groupby_col": "cell_type",
                    "comparisons": [["treat", "ctrl"]],
                },
            },
        },
        {},
        {"annotation_source": "input_obs", "label_col": "cell_type"},
    )
    assert res["status"] == "success"
    pwp = res["pseudobulk_pathways"]
    assert pwp["ora_method"] == "local_hypergeometric"
    assert pwp["gene_set_versions"]["GO_BP"]["release"] == "2021"
    assert pwp["databases_skipped"] == ["KEGG"]
    assert pwp["output_csv"] == str(tmp_path / "pw.csv")

    pb = res["pseudobulk_de"]
    assert pb["lfc_shrinkage"] == {"requested": True, "method": "apeGLM"}
    assert pb["power"] == {"alpha_nominal": 0.05}
    assert pb["robustness_multiverse"] == {"method": "per-block"}
    assert pb["output_csv"] == str(tmp_path / "pb.csv")
    assert pb["full_results_csv"] == str(tmp_path / "pb_all.csv")

    da = res["differential_abundance"]
    assert da["model"]["paired_donor_fixed_effects"] is True


def test_bi04_narrative_adapter_forwards_provenance_and_mt_zero():
    # _first_mt_threshold reads mt_threshold_used first and preserves 0.0
    assert adapter._first_mt_threshold([{"mt_threshold_used": 0.0, "mt_threshold": 15.0}]) == 0.0
    assert adapter._first_mt_threshold([{"mt_threshold_used": 12.5}]) == 12.5
    assert adapter._first_mt_threshold([{"mt_threshold": 0.0, "mt_high": 20.0}]) == 0.0

    report = {
        "mode": "pseudobulk",
        "organism": "Homo sapiens",
        "stages": {
            "qc": {
                "status": "success",
                "n_samples": 1,
                "per_sample": [{"n_cells_before": 100, "n_cells_after": 90, "mt_threshold_used": 0.0}],
            },
            "celltypist": {
                "status": "success",
                "model_used": "Immune_All_Low.pkl",
                "model_source": "default_immune_fallback",
                "model_warning": "fallback warning",
                "confidence_available": True,
                "confidence_extraction_failed": False,
                "per_cluster": {"0": {"label": "T cell", "mean_confidence": 0.35, "n_cells": 50}},
                "warnings": ["ct warn"],
            },
            "pathways": {
                "status": "success",
                "organism": "Homo sapiens",
                "databases": {"GO_BP": "GO_Biological_Process_2021"},
                "ora_method": "local_hypergeometric",
                "gene_set_versions": {"GO_BP": {"release": "2021"}},
                "per_cluster": {"T::a_vs_b": {"n_significant": 1, "results": {}}},
            },
            "pseudobulk": {
                "status": "success",
                "groupby": "cell_type",
                "condition_col": "condition",
                "replicate_col": "donor",
                "multiple_testing": {"fdr_strategy": "per_cluster", "n_tests_global": 500},
                "lfc_shrinkage": {"requested": True, "method": "apeGLM"},
                "count_source": "recovered_from_lognorm",
                "lognorm_recovered": True,
                "n_groups": 1,
                "per_group": {},
            },
            "cell_communication": {
                "status": "success",
                "method": "LIANA (specificity_rank)",
                "n_cell_types": 3,
                "n_interactions": 80,
                "n_significant": 12,
                "n_autocrine_dropped": 5,
                "top_pairs": ["A→B"],
            },
        },
    }
    bundle = adapter.adapt(report)
    f = bundle["agent_results"]["scrna_agent"]["findings"]["scRNA"]["findings"]
    assert f["qc"]["mt_threshold"] == 0.0
    assert f["qc"]["mt_threshold_used"] == 0.0
    assert f["cell_types"]["celltypist"]["per_cluster"]["0"]["mean_confidence"] == 0.35
    assert f["cell_types"]["celltypist"]["model_warning"] == "fallback warning"
    assert f["cell_types"]["celltypist"]["warnings"] == ["ct warn"]
    assert f["pseudobulk_pathways"]["ora_method"] == "local_hypergeometric"
    assert f["pseudobulk_pathways"]["gene_set_versions"]["GO_BP"]["release"] == "2021"
    assert f["pseudobulk_de"]["multiple_testing"]["fdr_strategy"] == "per_cluster"
    assert f["pseudobulk_de"]["lfc_shrinkage"]["method"] == "apeGLM"
    assert f["pseudobulk_de"]["count_source"] == "recovered_from_lognorm"
    assert f["pseudobulk_de"]["lognorm_recovered"] is True
    assert f["cell_communication"]["n_significant"] == 12


@needs_agent
def test_bi05_cellcomm_wording_honest_significant_vs_tested():
    from aria.agents.narrative.scrna.text import summarize_scrna_text

    ccc_base = {
        "status": "done",
        "method": "LIANA (specificity_rank)",
        "n_cell_types": 3,
        "n_interactions": 80,
        "n_significant": 12,
        "n_autocrine_dropped": 5,
        "top_pairs": ["A→B"],
    }
    text_with_sig = summarize_scrna_text({"cell_communication": ccc_base})
    assert "12 significant L-R interactions (80 tested non-autocrine)" in text_with_sig

    text_zero_sig = summarize_scrna_text({
        "cell_communication": {**ccc_base, "n_significant": 0}
    })
    assert "0 significant L-R interactions (80 tested non-autocrine)" in text_zero_sig

    ccc_no_sig = {k: v for k, v in ccc_base.items() if k != "n_significant"}
    text_no_sig = summarize_scrna_text({"cell_communication": ccc_no_sig})
    assert "80 tested non-autocrine L-R interactions" in text_no_sig
    assert "80 significant L-R interactions" not in text_no_sig
