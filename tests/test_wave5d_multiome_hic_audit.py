"""Wave 5D audit remediation tests: Domain 3 (Multiome) + Domain 4 (Hi-C/Micro-C) + CI Docker."""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


# ═══════════════════════════════════════════════════════════════════════════
# PART 1: DOMAIN 3 (MULTIOME)
# ═══════════════════════════════════════════════════════════════════════════


def test_d3_c1_paired_qc_preserves_uns_and_multilevel_sample_id(tmp_path: Path):
    ad = pytest.importorskip("anndata")
    md = pytest.importorskip("mudata")
    from aria.scripts.integration_paired_qc import integration_paired_qc

    n = 12
    barcodes = [f"D{i % 2}#CELL_{i:02d}" for i in range(n)]
    sample_ids = ["donorA" if i % 2 == 0 else "donorB" for i in range(n)]
    rna = ad.AnnData(
        X=sp.csr_matrix(np.full((n, 6), 50, dtype=np.int64)),
        obs=pd.DataFrame({"sample_id": sample_ids, "donor_id": sample_ids}, index=barcodes),
        var=pd.DataFrame(index=[f"G{j}" for j in range(5)] + ["MT-CO1"]),
    )
    atac = ad.AnnData(
        X=sp.csr_matrix(np.full((n, 6), 50, dtype=np.int64)),
        obs=pd.DataFrame({"sample_id": sample_ids, "donor_id": sample_ids}, index=barcodes),
        var=pd.DataFrame(index=[f"chr1:{j * 500 + 100}-{j * 500 + 400}" for j in range(6)]),
    )
    mdata = md.MuData({"rna": rna, "atac": atac})
    mdata.obs["sample_id"] = sample_ids
    mdata.uns["aria_cohort_merge"] = {
        "n_samples": 2,
        "samples": ["donorA", "donorB"],
        "nullable_ignored": None,
    }
    in_path = tmp_path / "cohort_merged.h5mu"
    mdata.write_h5mu(in_path)

    res = integration_paired_qc(
        {
            "h5mu_path": str(in_path),
            "output_dir": str(tmp_path / "out_qc"),
            "sample_id": "cohort_scalar_should_not_overwrite",
            "min_genes": 2,
            "min_rna_counts": 10,
            "min_peaks": 2,
            "min_atac_counts": 10,
            "max_pct_mt": 50.0,
            "resume": False,
        }
    )
    assert res["status"] == "success", res
    out_md = md.read_h5mu(res["output_path"])
    assert "aria_cohort_merge" in out_md.uns
    assert list(out_md.uns["aria_cohort_merge"]["samples"]) == ["donorA", "donorB"]
    assert set(out_md.obs["sample_id"].astype(str)) == {"donorA", "donorB"}
    assert set(out_md.mod["rna"].obs["sample_id"].astype(str)) == {"donorA", "donorB"}
    assert set(out_md.mod["atac"].obs["sample_id"].astype(str)) == {"donorA", "donorB"}


def test_d3_h1_wnn_diagnostics_nan_weights_axis1_and_axis5():
    from aria.utils.wnn_diagnostics import compute_wnn_anti_degradation_diagnostics

    rng = np.random.default_rng(42)
    n = 40
    # Construct informative RNA manifold (2 tight clusters in 6D) and pure noise ATAC manifold.
    rna_emb = np.vstack([
        rng.normal(loc=4.0, scale=0.3, size=(20, 6)),
        rng.normal(loc=-4.0, scale=0.3, size=(20, 6)),
    ])
    atac_emb = rng.normal(loc=0.0, scale=1.0, size=(n, 6))

    # Build a joint kNN graph mostly preserving RNA neighborhoods.
    rows, cols, data = [], [], []
    for i in range(n):
        base = 0 if i < 20 else 20
        for off in range(1, 7):
            j = base + ((i - base + off) % 20)
            rows.append(i)
            cols.append(j)
            data.append(0.8)
    g = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    labels = ["c0"] * 20 + ["c1"] * 20

    # Inject a single NaN into weights while ATAC mean weight is low (0.10).
    rna_w = np.full(n, 0.90)
    atac_w = np.full(n, 0.10)
    rna_w[0] = np.nan
    atac_w[0] = np.nan

    diag = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=rna_emb,
        atac_embedding=atac_emb,
        wnn_distances=g,
        graph_is_connectivity=True,
        wnn_clusters=labels,
        rna_clusters=labels,
        atac_clusters=["a0"] * 20 + ["a1"] * 20,
        rna_weights=rna_w,
        atac_weights=atac_w,
        k=5,
        seed=0,
    )
    split_dim = diag["held_out_check"]["split_dimension_oof"]
    mae = split_dim["oof_vs_insample_rna_weight_mae"]
    assert mae is None or math.isfinite(mae)
    # NaN weight must NOT disable uninformative-modality alert!
    assert "atac_uninformative_downweighted" in diag["issue_codes"]

    # All-NaN weights must return None for oof_vs_insample_rna_weight_mae (never NaN).
    diag_all_nan = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=rna_emb,
        atac_embedding=atac_emb,
        wnn_distances=g,
        graph_is_connectivity=True,
        wnn_clusters=labels,
        rna_clusters=labels,
        atac_clusters=labels,
        rna_weights=np.full(n, np.nan),
        atac_weights=np.full(n, np.nan),
        k=5,
        seed=0,
    )
    assert diag_all_nan["held_out_check"]["split_dimension_oof"]["oof_vs_insample_rna_weight_mae"] is None


def test_d3_h2_h11_harmony_qa_collapse_atac_status_and_p2g_csv_and_dup_barcodes(tmp_path: Path):
    from aria.scripts.integration_wnn import _evaluate_harmony_qa
    from aria.scripts.integration_eregulons import _load_precomputed_p2g_csv

    # 1. Harmony 1-cluster biological collapse (wnn_cl_sil is None) must still trigger QA warning.
    qa_collapse = _evaluate_harmony_qa(
        should_harmonize=True,
        n_joint_clusters=1,
        rna_sil_before=0.35,
        rna_sil_after=0.01,
        atac_sil_before=0.30,
        atac_sil_after=0.02,
        wnn_cl_sil=None,
    )
    assert qa_collapse is not None
    assert qa_collapse["status"] == "warnings"
    assert qa_collapse["possible_overcorrection"] is True
    assert "collapsed to 1 cluster" in (qa_collapse.get("warning") or "")

    # 2. Clean RNA + ATAC-only residual_batch_effect must set qa["status"] = "warnings".
    qa_atac_only = _evaluate_harmony_qa(
        should_harmonize=True,
        n_joint_clusters=4,
        rna_sil_before=0.25,
        rna_sil_after=0.02,
        atac_sil_before=0.25,
        atac_sil_after=0.35,  # worsened + residual_batch_effect on ATAC only
        wnn_cl_sil=0.45,
    )
    assert qa_atac_only is not None
    assert qa_atac_only["status"] == "warnings"
    assert any(it.get("modality") == "atac" for it in qa_atac_only["issues"])
    assert qa_atac_only.get("warning")

    # 3. _load_precomputed_p2g_csv rejects padj=None when max_padj is set, and counts n_fdr_significant strictly <= 0.05.
    csv_p = tmp_path / "p2g.csv"
    pd.DataFrame(
        [
            {"gene": "G1", "peak": "chr1:100-500", "correlation": 0.4, "padj": 0.01, "distance_bp": 1000},
            {"gene": "G2", "peak": "chr1:600-900", "correlation": 0.5, "padj": "", "distance_bp": 500},
            {"gene": "G3", "peak": "chr1:1000-1500", "correlation": 0.6, "padj": 0.08, "distance_bp": 200},
        ]
    ).to_csv(csv_p, index=False)

    res_strict = _load_precomputed_p2g_csv(csv_p, min_corr=0.05, max_padj=0.05)
    assert [lk["gene"] for lk in res_strict["links"]] == ["G1"]
    assert res_strict["n_fdr_significant"] == 1

    res_loose = _load_precomputed_p2g_csv(csv_p, min_corr=0.05, max_padj=0.10)
    assert [lk["gene"] for lk in res_loose["links"]] == ["G1", "G3"]
    assert res_loose["n_fdr_significant"] == 1  # G3 (padj=0.08) NOT counted as <= 0.05!

    res_unfiltered = _load_precomputed_p2g_csv(csv_p, min_corr=0.05, max_padj=None)
    assert len(res_unfiltered["links"]) == 3
    assert res_unfiltered["n_fdr_significant"] == 1


def test_d3_h11_split_file_duplicate_barcodes_rejected(tmp_path: Path):
    ad = pytest.importorskip("anndata")
    pytest.importorskip("scanpy")
    from aria.scripts.integration_eregulons import integration_eregulons

    rna = ad.AnnData(
        X=sp.csr_matrix(np.ones((4, 3), dtype=np.int64)),
        obs=pd.DataFrame(index=["BC1", "BC1", "BC2", "BC3"]),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    atac = ad.AnnData(
        X=sp.csr_matrix(np.ones((4, 3), dtype=np.int64)),
        obs=pd.DataFrame(index=["BC1", "BC2", "BC3", "BC4"]),
        var=pd.DataFrame(index=["chr1:100-200", "chr1:300-400", "chr1:500-600"]),
    )
    rna_p = tmp_path / "rna_dup.h5ad"
    atac_p = tmp_path / "atac_ok.h5ad"
    rna.write_h5ad(rna_p)
    atac.write_h5ad(atac_p)

    # eRegulons split-file branch (needs minimal GTF + FASTA + motif file to reach split ingestion)
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('chr1\tsrc\tgene\t1000\t2000\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n')
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "ACGT" * 250 + "\n")
    (tmp_path / "ref.fa.fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    meme = tmp_path / "motifs.meme"
    meme.write_text(
        "MEME version 4\nALPHABET= ACGT\nMOTIF M1 G1\nletter-probability matrix: alength= 4 w= 4\n"
        "0.7 0.1 0.1 0.1\n0.1 0.7 0.1 0.1\n0.1 0.1 0.7 0.1\n0.1 0.1 0.1 0.7\n"
    )
    res_ereg = integration_eregulons(
        {
            "rna_files": [str(rna_p)],
            "atac_files": [str(atac_p)],
            "gtf_path": str(gtf),
            "genome_fasta": str(fa),
            "motif_file": str(meme),
            "output_dir": str(tmp_path / "ereg_out"),
            "resume": False,
        }
    )
    assert res_ereg["status"] == "error"
    assert res_ereg["error_type"] == "DuplicateBarcodes"


def test_d3_h11_split_file_duplicate_barcodes_rejected_wnn(tmp_path: Path):
    ad = pytest.importorskip("anndata")
    pytest.importorskip("muon")
    pytest.importorskip("scanpy")
    from aria.scripts.integration_wnn import integration_wnn

    rna = ad.AnnData(
        X=sp.csr_matrix(np.ones((4, 3), dtype=np.int64)),
        obs=pd.DataFrame(index=["BC1", "BC1", "BC2", "BC3"]),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    atac = ad.AnnData(
        X=sp.csr_matrix(np.ones((4, 3), dtype=np.int64)),
        obs=pd.DataFrame(index=["BC1", "BC2", "BC3", "BC4"]),
        var=pd.DataFrame(index=["chr1:100-200", "chr1:300-400", "chr1:500-600"]),
    )
    rna_p = tmp_path / "rna_dup.h5ad"
    atac_p = tmp_path / "atac_ok.h5ad"
    rna.write_h5ad(rna_p)
    atac.write_h5ad(atac_p)

    res_wnn = integration_wnn(
        {
            "rna_files": [str(rna_p)],
            "atac_files": [str(atac_p)],
            "output_dir": str(tmp_path / "wnn_out"),
            "resume": False,
        }
    )
    assert res_wnn["status"] == "error"
    assert res_wnn["error_type"] == "DuplicateBarcodes"


def test_d3_h6_h10_assemble_suffix_reconciliation_unmapped_fragments_and_contigs(tmp_path: Path):
    ad = pytest.importorskip("anndata")
    from aria.scripts.integration_assemble_multiome import (
        _quantify_fragments,
        _revcomp,
        _translate_barcodes,
        integration_assemble_multiome,
    )
    from aria.utils.scatac_cohort_union import _peak_intervals, resolve_sample_ids

    # 1. resolve_sample_ids raises ValueError on mismatched explicit sample_ids length.
    with pytest.raises(ValueError, match="length"):
        resolve_sample_ids(["a.h5ad", "b.h5ad"], sample_ids=["only_one"])

    # 2. _peak_intervals and _quantify_fragments normalize Ensembl vs UCSC contigs.
    adata_ens = ad.AnnData(
        X=sp.csr_matrix(np.ones((2, 2), dtype=np.int64)),
        obs=pd.DataFrame(index=["C1", "C2"]),
        var=pd.DataFrame(index=["1:100-200", "MT:10-50"]),
    )
    assert _peak_intervals(adata_ens) == [("chr1", 100, 200), ("chrM", 10, 50)]

    peaks_bed = tmp_path / "peaks_ens.bed"
    peaks_bed.write_text("1\t100\t300\n", encoding="utf-8")
    frags_ucsc = tmp_path / "frags_ucsc.tsv"
    frags_ucsc.write_text("chr1\t150\t200\tBC1\t7\n", encoding="utf-8")
    q_ad = _quantify_fragments(str(frags_ucsc), str(peaks_bed))
    assert list(q_ad.var_names) == ["chr1:100-300"]
    assert int(q_ad.X.toarray()[0, 0]) == 7

    # 3. 10x_arc_whitelist_pair reconciles -1 suffix asymmetry + reverse_complement orientation.
    gex_wl_entries = ["AAAA", "CCCC", "ACAC"]
    atac_wl_entries = ["AAAC", "CCCA", "ACCG"]
    gex_wl = tmp_path / "gex_wl.txt"
    atac_wl = tmp_path / "atac_wl.txt"
    gex_wl.write_text("\n".join(gex_wl_entries) + "\n")
    atac_wl.write_text("\n".join(atac_wl_entries) + "\n")

    rna_suf = ad.AnnData(
        X=np.ones((3, 2)),
        obs=pd.DataFrame(index=[b + "-1" for b in gex_wl_entries]),
    )
    # ATAC is bare (no -1) AND reverse-complemented (chromap-style).
    atac_bare_rc = ad.AnnData(
        X=np.ones((3, 2)),
        obs=pd.DataFrame(index=[_revcomp(b) for b in atac_wl_entries]),
    )
    tr, keep, names = _translate_barcodes(
        rna_suf,
        atac_bare_rc,
        {
            "barcode_namespace": "10x_arc_whitelist_pair",
            "gex_whitelist_path": str(gex_wl),
            "atac_whitelist_path": str(atac_wl),
        },
    )
    assert tr["orientation"] == "reverse_complement"
    assert keep == [0, 1, 2]
    assert names == ["AAAA-1", "CCCC-1", "ACAC-1"]

    # 4. Unmapped fragment barcodes skipped when translation is active (never contaminate RNA barcode).
    mudata = pytest.importorskip("mudata")
    rna_bcs = ["RNA_1", "RNA_2"]
    atac_bcs = ["ATAC_1", "ATAC_2"]
    rna_h5 = tmp_path / "rna.h5ad"
    ad.AnnData(
        X=sp.csr_matrix(np.ones((2, 3), dtype=np.int64)),
        obs=pd.DataFrame(index=rna_bcs),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    ).write_h5ad(rna_h5)
    bmap = tmp_path / "bmap.tsv"
    bmap.write_text("rna_barcode\tatac_barcode\nRNA_1\tATAC_1\nRNA_2\tATAC_2\n")
    frags = tmp_path / "frags_collide.tsv"
    # Unmapped ATAC fragment whose raw string happens to equal "RNA_1" must be IGNORED!
    frags.write_text(
        "chr1\t150\t200\tATAC_1\t4\n"
        "chr1\t150\t200\tATAC_2\t6\n"
        "chr1\t150\t200\tRNA_1\t999\n"
    )
    res = integration_assemble_multiome(
        {
            "mode": "paired_matrices",
            "rna_path": str(rna_h5),
            "fragments_path": str(frags),
            "peaks_path": str(peaks_bed),
            "barcode_map_path": str(bmap),
            "min_shared_cells": 2,
            "output_dir": str(tmp_path / "assembled"),
            "resume": False,
        }
    )
    assert res["status"] == "success", res
    m_out = mudata.read_h5mu(res["output_path"])
    assert float(m_out.mod["atac"].obs.loc["RNA_1", "atac_fragments"]) == 4.0
    assert float(m_out.mod["atac"].obs.loc["RNA_2", "atac_fragments"]) == 6.0


def test_d3_h12_multiome_contracts_fail_closed_arc_raw_and_unrecognized_type(tmp_path: Path, monkeypatch):
    from aria.utils.multiome_contracts import (
        _drop_arc_raw_siblings,
        _paired_h5mu_records,
        _unpaired_container_files,
        enforce_same_cell_integration_gate,
        validate_multiome_contract,
    )

    # 1. Empty 0-byte .h5mu on disk fails closed (NOT paired_rna=True).
    empty_h5mu = tmp_path / "empty.h5mu"
    empty_h5mu.write_bytes(b"")
    recs_empty = _paired_h5mu_records({"modalities": {"scATAC": [str(empty_h5mu)]}})
    assert recs_empty == []

    # 2. Corrupt non-HDF5 .h5mu on disk fails closed.
    corrupt_h5mu = tmp_path / "corrupt.h5mu"
    corrupt_h5mu.write_bytes(b"NOT_AN_HDF5_FILE_CONTENT")
    recs_corrupt = _paired_h5mu_records({"modalities": {"scATAC": [str(corrupt_h5mu)]}})
    assert recs_corrupt == []

    # 3. raw_feature_bc_matrix/matrix.mtx.gz sibling is NOT falsely flagged as unpaired container.
    outs = tmp_path / "sample1" / "outs"
    filt_mex = outs / "filtered_feature_bc_matrix"
    raw_mex = outs / "raw_feature_bc_matrix"
    filt_mex.mkdir(parents=True)
    raw_mex.mkdir(parents=True)
    (filt_mex / "matrix.mtx.gz").write_bytes(b"x")
    (raw_mex / "matrix.mtx.gz").write_bytes(b"x")
    other_raw = tmp_path / "sample2" / "outs" / "raw_feature_bc_matrix" / "matrix.mtx.gz"
    other_raw.parent.mkdir(parents=True)
    other_raw.write_bytes(b"x")

    exp_ctx = {
        "modalities": {
            "scRNA": [str(raw_mex / "matrix.mtx.gz")],
            "scATAC": [str(filt_mex / "matrix.mtx.gz")],
        }
    }
    assert _unpaired_container_files(exp_ctx, [str(filt_mex)]) == []
    assert _unpaired_container_files(exp_ctx, [str(filt_mex / "matrix.mtx.gz")]) == []

    # Separate donor's raw matrix IS still flagged as unpaired.
    exp_ctx_other = {"modalities": {"scRNA": [str(other_raw)], "scATAC": [str(filt_mex)]}}
    assert _unpaired_container_files(exp_ctx_other, [str(filt_mex)]) == [str(other_raw)]

    # 4. _drop_arc_raw_siblings & inspected_h5mu_paths normalize relative vs absolute paths.
    monkeypatch.chdir(tmp_path)
    dropped = _drop_arc_raw_siblings(
        [
            {"path": str(filt_mex.resolve())},
            {"path": "sample1/outs/raw_feature_bc_matrix/matrix.mtx.gz"},
        ]
    )
    assert len(dropped) == 1

    # Inspected unpaired .h5mu recorded as absolute path blocks relative path duplicate in modalities.
    unpaired_p = tmp_path / "inspected_unpaired.h5mu"
    unpaired_p.write_bytes(b"x")
    recs_insp = _paired_h5mu_records(
        {
            "assay_detections": [
                {
                    "path": str(unpaired_p.resolve()),
                    "evidence": {"format": "h5mu", "paired_rna": False},
                }
            ],
            "modalities": {"scATAC": ["inspected_unpaired.h5mu"]},
        }
    )
    assert recs_insp == []

    # 5. Unrecognized object_type emits blocking finding and fails gate.
    val = validate_multiome_contract(
        {
            "genome": "hg38",
            "multiome_contract": {
                "object_type": "bogus_unknown_type",
                "same_cell": True,
                "cell_namespace": "shared_barcodes",
            },
        }
    )
    assert val["status"] == "red"
    assert any(
        f["check"] == "multiome_contract_unrecognized_object_type" and f["severity"] == "blocking"
        for f in val["findings"]
    )
    gate = enforce_same_cell_integration_gate(
        {
            "genome": "hg38",
            "multiome_contract": {
                "object_type": "bogus_unknown_type",
                "same_cell": True,
                "cell_namespace": "shared_barcodes",
            },
        }
    )
    assert gate["allowed"] is False


def test_d3_h8_h9_noncausal_labels_batch_hatch_and_mofa_bio_count_parity(tmp_path: Path, monkeypatch):
    import aria.agents._narrative_integration as ni
    from aria.agents.narrative.narrators.integration import IntegrationNarrator
    from aria.agents.narrative.validators import find_causal_language
    from aria.agents.narrative_agent import NarrativeAgent

    # 1. _SIGN_CLASS_STYLE has zero causal terms ("activator", "repressive").
    for cls, (_col, _tag, label) in ni._SIGN_CLASS_STYLE.items():
        low = label.lower()
        assert "activator" not in low
        assert "repressive" not in low
        assert find_causal_language(label) is None

    # 2. render_mofa_variance_explained hatches batch_flag/batch_factor_ids even when exclude_batch_factors=False.
    captured = []
    real_save = ni._save_dual

    def _cap(fig, p, dpi=220):
        captured.append(fig)
        return real_save(fig, p, dpi=60)

    monkeypatch.setattr(ni, "_save_dual", _cap)
    mofa_unexcluded_batch = {
        "status": "success",
        "n_cells": 1000,
        "n_factors": 4,
        "n_active_factors": 3,
        "views": {"scRNA": {"kind": "rna"}, "scATAC": {"kind": "atac"}},
        "variance_explained": {"scRNA": 0.30, "scATAC": 0.20},
        "variance_explained_additive": {"scRNA": 0.32, "scATAC": 0.21},
        "top_factors": [
            {"factor_id": 1, "active": True, "category": "shared", "r2_by_view": {"scRNA": 0.15, "scATAC": 0.10}},
            {
                "factor_id": 2,
                "active": True,
                "category": "shared",
                "batch_flag": True,
                "batch_association": {"batch_key": "donor_id"},
                "excluded_from_biological_interpretation": False,
                "exclusion_reasons": [],
                "r2_by_view": {"scRNA": 0.10, "scATAC": 0.08},
            },
            {
                "factor_id": 3,
                "active": True,
                "category": "modality_specific",
                "excluded_from_biological_interpretation": True,
                "exclusion_reasons": ["batch_or_donor_confounded:donor_id(eta2=0.52)"],
                "r2_by_view": {"scRNA": 0.05, "scATAC": 0.02},
            },
        ],
        "shared_factors": [{"factor_id": 1}, {"factor_id": 2}],
        "modality_specific_factors": {"scRNA": [{"factor_id": 3}]},
        "excluded_factors": [
            {
                "factor_id": 3,
                "r2_by_view": {"scRNA": 0.05, "scATAC": 0.02},
                "exclusion_reasons": ["batch_or_donor_confounded:donor_id(eta2=0.52)"],
            }
        ],
        "batch_factor_ids": [2, 3],
    }
    assert ni.render_mofa_variance_explained(mofa_unexcluded_batch, tmp_path)
    fig = captured[-1]
    ax_hm, ax_bar = fig.axes[0], fig.axes[1]
    yt = [t.get_text() for t in ax_hm.get_yticklabels()]
    assert any("Factor 2" in s and "[EXCLUDED]" in s for s in yt)
    bar_hatches = [p.get_hatch() for p in ax_bar.patches]
    assert "///" in bar_hatches

    # 3. IntegrationNarrator: n_bio_active == n_shared + n_spec, and caveat/methods reflect actual batch/donor reasons (not depth/cell-cycle!).
    narr = IntegrationNarrator()
    block = narr._mofa_block(mofa_unexcluded_batch, {})
    assert block.metrics["n_active_biological_factors"] == 1
    assert block.metrics["n_shared_factors"] == 1
    assert block.metrics["n_modality_specific_factors"] == 0
    assert block.metrics["n_shared_factors"] + block.metrics["n_modality_specific_factors"] == block.metrics["n_active_biological_factors"]
    excl_cav = [c.text for c in block.caveats if "quarantined" in c.text][0]
    assert "batch/donor" in excl_cav
    assert "sequencing depth" not in excl_cav
    assert "cell-cycle" not in excl_cav

    methods_txt = " ".join(
        narr.methods("integration_agent", {"status": "done", "findings": {"mofa": mofa_unexcluded_batch}})
    )
    assert "batch/donor covariates" in methods_txt
    assert "|r| >= 0.70" not in methods_txt

    # 4. NarrativeAgent uses _mofa_n_bio_active (1), NOT raw n_active_factors (3).
    agent_results = {"integration_agent": {"status": "done", "findings": {"mofa": mofa_unexcluded_batch}}}
    conc = NarrativeAgent.__new__(NarrativeAgent)._summarize_agent_results_for_llm(agent_results)
    exec_sum = NarrativeAgent._deterministic_integration_executive_summary(agent_results)
    assert "MOFA+ retained 1 active biological of 4" in conc
    assert "MOFA+ retained 1 active biological of 4" in exec_sum


# ═══════════════════════════════════════════════════════════════════════════
# PART 2: DOMAIN 4 (HI-C / MICRO-C)
# ═══════════════════════════════════════════════════════════════════════════


def test_d4_c1_m2_hic_figure_unphased_none_and_cis_only_window(tmp_path: Path):
    from aria.scripts.hic_figure import (
        _pick_mappable_window,
        _render_compartment_insulation_track,
    )

    # 1. pc1_validated=None must fail closed to unphased labels (never "A (+E1)").
    comp_tsv = tmp_path / "comp.tsv"
    pd.DataFrame(
        {
            "chrom": ["chr1"] * 10,
            "start": np.arange(10) * 100_000,
            "end": (np.arange(10) + 1) * 100_000,
            "E1": [0.5, 0.4, -0.3, -0.5, 0.2, 0.6, -0.4, -0.2, 0.3, -0.1],
        }
    ).to_csv(comp_tsv, sep="\t", index=False)
    bundle = _render_compartment_insulation_track(str(comp_tsv), None, tmp_path, pc1_validated=None)
    assert bundle is not None
    svg = Path(bundle["svg"]).read_text(encoding="utf-8")
    assert "+E1 (unphased)" in svg
    assert "-E1 (unphased)" in svg
    assert "A (+E1)" not in svg

    # 2. _pick_mappable_window must use cis sparse fetch and NEVER call clr.pixels().
    class _CisOnlyMockClr:
        chromsizes = {"chr1": 100_000}

        def bins(self):
            df = pd.DataFrame(
                {
                    "chrom": ["chr1"] * 10,
                    "start": np.arange(10) * 10_000,
                    "end": (np.arange(10) + 1) * 10_000,
                    "weight": [np.nan] * 10,
                }
            )
            return SimpleNamespace(fetch=lambda ch: df)

        def pixels(self):
            raise AssertionError("clr.pixels() must NEVER be called for genome-wide load!")

        def matrix(self, balance=False, sparse=True):
            assert balance is False and sparse is True
            m = sp.coo_matrix(([5.0], ([4], [6])), shape=(10, 10)).tocsr()
            return SimpleNamespace(fetch=lambda ch: m)

    start, end = _pick_mappable_window(_CisOnlyMockClr(), "chr1", 10_000, 4, has_w=True)
    assert start == 40_000 and end == 80_000


def test_d4_c5_box_smooth_symmetry_and_constant_one_preservation():
    from aria.scripts.hic_qc_and_balance import _box_smooth_2d

    rng = np.random.default_rng(123)
    for n, h in [(15, 1), (25, 2), (40, 3)]:
        raw = rng.uniform(0.0, 10.0, size=(n, n))
        m_sym = sp.csr_matrix(0.5 * (raw + raw.T))
        sm = _box_smooth_2d(m_sym, n, h).toarray()
        max_asym = float(np.max(np.abs(sm - sm.T)))
        assert max_asym < 1e-12, f"Asymmetry {max_asym} for n={n}, h={h}"

        ones = sp.csr_matrix(np.ones((n, n), dtype=float))
        sm_ones = _box_smooth_2d(ones, n, h).toarray()
        assert np.allclose(sm_ones, 1.0, atol=1e-12), f"Constant-1 violated for n={n}, h={h}"


def test_d4_c4_list_mcool_resolutions_modern_and_legacy(monkeypatch, tmp_path: Path):
    import sys
    from aria.utils.hic_constants import list_mcool_resolutions

    # Pure unit test with mocked cooler module for legacy /0, /1 paths.
    fake_binsizes = {"/0": 100_000, "/1": 40_000, "/2": 10_000}

    class _FakeCooler:
        def __init__(self, uri: str):
            grp = uri.split("::", 1)[1]
            self.binsize = fake_binsizes[grp]

    fake_cooler = SimpleNamespace(
        fileops=SimpleNamespace(list_coolers=lambda p: ["/0", "/1", "/2"]),
        Cooler=_FakeCooler,
    )
    monkeypatch.setitem(sys.modules, "cooler", fake_cooler)
    res_map = list_mcool_resolutions(str(tmp_path / "legacy.mcool"))
    assert res_map == {10_000: "/2", 40_000: "/1", 100_000: "/0"}


def test_d4_c3_h4_genome_arch_agent_coarsest_fallback_and_balanced_sample_ids(tmp_path: Path):
    from aria.agents.genome_arch_agent import GenomeArchAgent
    from aria.memory.memory import ARIAMemory

    mem = ARIAMemory(db_path=str(tmp_path / "ga.sqlite"))
    agent = GenomeArchAgent(mem, MagicMock())
    agent.env = MagicMock()

    # 1. When target_res (40kb) > max(available), pick coarsest (10kb), NOT finest (1kb)!
    dec = agent._advise_resolution(
        "exp_coarsest",
        {"user_question": "Analyze TADs"},
        {"available_resolutions": [1_000, 2_000, 5_000, 10_000]},
        {"organism": "Homo sapiens"},
    )
    assert dec["recommended_resolution"] == 10_000

    # 2. _run_topology passes balanced_sample_ids aligned 1-to-1 with balanced_files to differential_3d.
    captured_diff: list[dict] = []

    def _fake_run(*, stack, script_path, params):
        if params.get("analysis") == "differential_3d":
            captured_diff.append(dict(params))
            return {"status": "success", "ran": True, "n_ab_switches": 0}
        if params.get("analysis") == "compartments":
            return {"status": "success", "ran": True, "pct_A": 50.0, "pct_B": 50.0, "pc1_validated": True}
        if params.get("analysis") == "insulation_calibration":
            return {"status": "success", "boundary_strengths": {"200000": 0.5}}
        return {"status": "success", "ran": True}

    agent.env.run_in_stack.side_effect = _fake_run
    qc_result = {
        "balanced_files": [str(tmp_path / "s1.bal.mcool"), str(tmp_path / "s3.bal.mcool")],
        "balanced_sample_ids": ["s1", "s3"],
    }
    agent._run_topology(
        "exp_sids",
        {
            "genome": "hg38",
            "organism": "Homo sapiens",
            "sample_ids": ["s1", "s2_failed", "s3"],
            "sample_metadata": [
                {"sample_id": "s1", "condition": "WT", "replicate": "r1"},
                {"sample_id": "s2_failed", "condition": "WT", "replicate": "r2"},
                {"sample_id": "s3", "condition": "KO", "replicate": "r1"},
            ],
            "comparisons": [("KO", "WT")],
        },
        {},
        qc_result["balanced_files"],
        {"recommended_resolution": 40_000, "analysis_type": "tad_calling", "available_resolutions": [40_000, 100_000]},
        qc_result,
        work_root=tmp_path,
    )
    assert len(captured_diff) == 1
    assert captured_diff[0]["sample_ids"] == ["s1", "s3"]


def test_d4_h1_h2_h7_hic_narrator_unphased_diff3d_and_skipped_contrasts():
    from aria.agents.narrative.narrators.hic import HiCNarrator
    from aria.agents.narrative.validators import find_causal_language

    narr = HiCNarrator()
    findings_unphased = {
        "topology": {
            "compartments": {
                "status": "success",
                "pct_A": 50.0,
                "pct_B": 50.0,
                "resolution": 100_000,
                "pc1_validated": False,
            },
            "differential_3d": {
                "status": "success",
                "ran": True,
                "descriptive_only": False,
                "pc1_validated": False,
                "n_ab_switches": 0,
                "n_unphased_sign_flips": 14,
                "n_differential_insulation_bins": 5,
                "comparisons": [
                    {
                        "contrast": "KO_vs_WT",
                        "status": "success",
                        "ran": True,
                        "descriptive_only": False,
                        "n_differential_compartment_bins": 20,
                        "n_ab_switches": 0,
                        "n_unphased_sign_flips": 14,
                        "n_differential_insulation_bins": 5,
                    },
                    {
                        "contrast": "MUT_vs_WT",
                        "status": "skipped",
                        "ran": False,
                        "descriptive_only": True,
                        "reason": "insufficient_biological_replicates (1 vs 3)",
                    },
                ],
            },
        }
    }
    blocks = {b.id: b for b in narr.collect("genome_arch_agent", {"findings": findings_unphased})}
    diff_blk = blocks["hic.differential_3d"]
    assert diff_blk.confidence == "low"
    assert diff_blk.metrics["n_comparisons"] == 1
    assert diff_blk.metrics["n_unphased_sign_flips"] == 14
    assert "14 unphased E1 sign flip(s)" in diff_blk.claim
    assert "0 compartment switch(es)" not in diff_blk.claim
    assert any("unphased sign flips" in c.text.lower() for c in diff_blk.caveats)
    assert any("MUT_vs_WT" in c.text and "skipped" in c.text for c in diff_blk.caveats)
    assert find_causal_language(diff_blk.claim) is None

    methods_txt = " ".join(narr.methods("genome_arch_agent", {"findings": findings_unphased}))
    assert "per-replicate phased E1" not in methods_txt
    assert "unphased E1" in methods_txt

    # When pc1_validated=True, methods() reports phased E1 and confidence stays medium.
    findings_phased = {
        "topology": {
            "differential_3d": {
                "status": "success",
                "ran": True,
                "descriptive_only": False,
                "pc1_validated": True,
                "n_ab_switches": 8,
                "n_differential_insulation_bins": 2,
                "comparisons": [
                    {
                        "contrast": "KO_vs_WT",
                        "status": "success",
                        "ran": True,
                        "n_differential_compartment_bins": 10,
                        "n_ab_switches": 8,
                        "n_differential_insulation_bins": 2,
                    }
                ],
            }
        }
    }
    blk_phased = narr._differential_3d_block(
        findings_phased["topology"]["differential_3d"], findings_phased
    )
    assert blk_phased.confidence == "medium"
    assert "8 compartment switch(es)" in blk_phased.claim
    methods_phased = " ".join(narr.methods("genome_arch_agent", {"findings": findings_phased}))
    assert "per-replicate phased E1" in methods_phased


def test_d4_h3_h8_h9_coords_match_ice_phasing_and_loops_early_skip(tmp_path: Path, monkeypatch):
    from aria.scripts import hic_topology as ht

    # 1. _coords_match ignores Categorical category-set mismatch.
    df_a = pd.DataFrame(
        {
            "chrom": pd.Categorical(["chr1", "chr1"], categories=["chr1", "chrUn_A"]),
            "start": [0, 10000],
            "end": [10000, 20000],
        }
    )
    df_b = pd.DataFrame(
        {
            "chrom": pd.Categorical(["chr1", "chr1"], categories=["chr1", "chrUn_B"]),
            "start": [0, 10000],
            "end": [10000, 20000],
        }
    )
    assert not df_a.equals(df_b)  # pandas Categorical .equals() falsely returns False!
    assert ht._coords_match(df_a, df_b) is True

    # 2. _phasing_valid_chroms intersects with valid ICE-balanced bins.
    phasing_df = pd.DataFrame(
        {
            "chrom": ["chr1"] * 10 + ["chr2"] * 10,
            "start": list(np.arange(10) * 10_000) + list(np.arange(10) * 10_000),
            "end": list((np.arange(10) + 1) * 10_000) + list((np.arange(10) + 1) * 10_000),
            # chr1 has GC variance on ICE-valid bins; chr2 has constant GC on ICE-valid bins and variance only on NaN-weight bins!
            "GC": [0.3, 0.4, 0.5, 0.6, 0.7, 0.35, 0.45, 0.55, 0.65, 0.5]
            + [0.5, 0.5, 0.5, 0.5, 0.5, 0.2, 0.8, 0.3, 0.7, 0.9],
        }
    )
    bins_with_w = phasing_df[["chrom", "start", "end"]].copy()
    bins_with_w["weight"] = [0.1] * 10 + [0.1] * 5 + [np.nan] * 5

    class _FakeBins:
        columns = ["chrom", "start", "end", "weight"]

        def __getitem__(self, sl):
            return bins_with_w

    fake_clr = SimpleNamespace(bins=lambda: _FakeBins())
    assert ht._phasing_valid_chroms(phasing_df, ["chr1", "chr2"], clr=fake_clr) == {"chr1"}


def test_d4_live_cooler_mcool_legacy_loops_skip_scc_and_partial_phasing(tmp_path: Path, monkeypatch):
    cooler = pytest.importorskip("cooler")
    pytest.importorskip("cooltools")
    import h5py
    from aria.scripts import hic_topology as ht
    from aria.scripts.hic_figure import _pick_mappable_window
    from aria.scripts.hic_qc_and_balance import _resolve_scc_uri, hic_qc_and_balance
    from aria.utils.hic_constants import list_mcool_resolutions

    # 1. Live legacy zoom-level .mcool (/0 -> 100kb, /1 -> 40kb) via list_mcool_resolutions.
    legacy_mcool = tmp_path / "legacy_zoom.mcool"
    for idx_g, (grp, bs) in enumerate([("/0", 100_000), ("/1", 40_000)]):
        nb = 2_000_000 // bs
        bins = pd.DataFrame(
            {
                "chrom": ["chr1"] * nb,
                "start": np.arange(nb) * bs,
                "end": (np.arange(nb) + 1) * bs,
            }
        )
        pix = pd.DataFrame(
            [(i, j, max(1, 20 - (j - i))) for i in range(nb) for j in range(i, min(nb, i + 10))],
            columns=["bin1_id", "bin2_id", "count"],
        )
        cooler.create_cooler(
            f"{legacy_mcool}::{grp}",
            bins,
            pix,
            mode="w" if idx_g == 0 else "a",
        )

    res_map = list_mcool_resolutions(str(legacy_mcool))
    assert res_map == {40_000: "/1", 100_000: "/0"}
    scc_uri, scc_res = _resolve_scc_uri(str(legacy_mcool), 10_000, None)
    assert scc_res == 40_000 and scc_uri.endswith("::/1")
    scc_uri_exp, scc_res_exp = _resolve_scc_uri(str(legacy_mcool), 40_000, 100_000)
    assert scc_res_exp == 100_000 and scc_uri_exp.endswith("::/0")

    # Verify _pick_mappable_window works live on a real cooler without calling clr.pixels().
    clr_live = cooler.Cooler(f"{legacy_mcool}::/1")
    s_win, e_win = _pick_mappable_window(clr_live, "chr1", 40_000, 20, has_w=False)
    assert s_win == 0 and e_win == 800_000

    # 2. _compute_loops skips coarse 100kb cooler BEFORE calling _ensure_balanced_cooler!
    coarse_cool = tmp_path / "coarse_100k.cool"
    bins_c = pd.DataFrame({"chrom": ["chr1"] * 10, "start": np.arange(10) * 100_000, "end": (np.arange(10) + 1) * 100_000})
    pix_c = pd.DataFrame({"bin1_id": [0, 1], "bin2_id": [1, 2], "count": [5, 5]})
    cooler.create_cooler(str(coarse_cool), bins_c, pix_c)

    def _boom(*args, **kwargs):
        raise AssertionError("_ensure_balanced_cooler must NOT be called before resolution gate!")

    monkeypatch.setattr(ht, "_ensure_balanced_cooler", _boom)
    loop_res = ht._compute_loops(str(coarse_cool), "hg38", 100_000, "all", str(tmp_path / "loops_out"))
    assert loop_res["status"] == "skipped"
    assert loop_res["reason"] == "resolution_too_coarse_for_dot_calling"
    monkeypatch.undo()

    # 3. hic_qc_and_balance works on legacy .mcool and emits balanced_sample_ids aligned 1-to-1 with balanced_files.
    qc_res = hic_qc_and_balance(
        {
            "files": [str(legacy_mcool)],
            "sample_ids": ["legacy_s1"],
            "resolution": 40_000,
            "output_dir": str(tmp_path / "qc_legacy"),
            "resume": False,
        }
    )
    assert qc_res["status"] == "success", qc_res
    assert qc_res["balanced_sample_ids"] == ["legacy_s1"]
    assert len(qc_res["balanced_files"]) == 1

    # 4. Partial phasing (chr1 phased in FASTA, chr2 missing from FASTA) in compartments + differential_3d!
    bs = 100_000
    n_per_ch = 24
    rng_fa = np.random.default_rng(77)
    fa_path = tmp_path / "only_chr1.fa"
    seq_parts = []
    for b in range(n_per_ch):
        p_gc = 0.65 if (b < 12) else 0.35
        probs = [(1 - p_gc) / 2, p_gc / 2, p_gc / 2, (1 - p_gc) / 2]
        seq_parts.append("".join(rng_fa.choice(["A", "C", "G", "T"], size=bs, p=probs)))
    full_seq = "".join(seq_parts)
    fa_path.write_text(">chr1\n" + full_seq + "\n", encoding="utf-8")

    def _make_2chrom_cooler(path: Path, flip_chr1: bool, flip_chr2: bool, seed: int, extra_contig: str | None = None):
        rng = np.random.default_rng(seed)
        bin_rows = []
        for ch in ["chr1", "chr2"] + ([extra_contig] if extra_contig else []):
            nb_c = n_per_ch if ch in {"chr1", "chr2"} else 2
            for b in range(nb_c):
                bin_rows.append((ch, b * bs, (b + 1) * bs))
        bins_df = pd.DataFrame(bin_rows, columns=["chrom", "start", "end"])
        pix_rows = []
        for ch_idx, (ch, flip) in enumerate([("chr1", flip_chr1), ("chr2", flip_chr2)]):
            off = ch_idx * n_per_ch
            # Checkerboard compartment pattern; KO flips bins [6..12) from +1 to -1.
            sgn = np.array([1 if b < 12 else -1 for b in range(n_per_ch)], dtype=int)
            if flip:
                sgn[6:12] = -1
            for i in range(n_per_ch):
                for j in range(i, n_per_ch):
                    mu = 30.0 * ((j - i + 1.0) ** (-0.7)) * (2.6 if sgn[i] == sgn[j] else 0.4)
                    cnt = int(rng.poisson(mu))
                    if cnt > 0:
                        pix_rows.append((off + i, off + j, cnt))
        pix_df = pd.DataFrame(pix_rows, columns=["bin1_id", "bin2_id", "count"])
        cooler.create_cooler(str(path), bins_df, pix_df)
        return str(path)

    c_files = []
    sids = ["WT_1", "WT_2", "KO_1", "KO_2"]
    meta = []
    for idx, sid in enumerate(sids):
        is_ko = sid.startswith("KO")
        # Give WT_2 a different extra unanalyzed contig to verify _coords_match ignores Categorical mismatch!
        extra = "chrUn_testA" if idx == 0 else ("chrUn_testB" if idx == 1 else None)
        c_files.append(
            _make_2chrom_cooler(tmp_path / f"{sid}.cool", is_ko, is_ko, 500 + idx, extra_contig=extra)
        )
        meta.append({"sample_id": sid, "condition": "KO" if is_ko else "WT", "replicate": f"r{(idx % 2) + 1}"})

    comp_res = ht._compute_compartments(
        c_files[0], "hg38", "Homo sapiens", bs, "all", str(tmp_path / "comp_partial"), genome_fasta=str(fa_path)
    )
    assert comp_res["status"] == "success", comp_res
    assert comp_res["pc1_validated"] is False
    assert comp_res["phased_chroms"] == ["chr1"]
    assert comp_res["unphased_chroms"] == ["chr2"]

    diff_res = ht._compute_differential_3d(
        {
            "files": c_files,
            "sample_ids": sids,
            "sample_metadata": meta,
            "comparisons": [("KO", "WT")],
            "comp_resolution": bs,
            "tad_resolution": bs,
            "window_size": bs * 3,
            "min_replicates_per_condition": 2,
            "genome_fasta": str(fa_path),
        },
        str(tmp_path / "diff_partial"),
    )
    assert diff_res["status"] == "success", diff_res
    assert diff_res["pc1_validated"] is False
    assert diff_res["phased_chroms"] == ["chr1"]
    assert diff_res["unphased_chroms"] == ["chr2"]
    comp_csv = pd.read_csv(diff_res["differential_compartments_csv"])
    chr2_classes = set(comp_csv[comp_csv["chrom"] == "chr2"]["switch_class"].unique())
    assert "A_to_B" not in chr2_classes and "B_to_A" not in chr2_classes
    assert "n_unphased_sign_flips" in diff_res
    if diff_res["n_unphased_sign_flips"] > 0:
        assert any(
            row["switch_class"] == "unphased_sign_flip"
            for row in diff_res["comparisons"][0]["top_switches"]
        )


# ═══════════════════════════════════════════════════════════════════════════
# PART 3: CI DOCKER FIX
# ═══════════════════════════════════════════════════════════════════════════


def test_ci01_dockerfile_rna_preinstalls_python_annoy_from_conda_forge():
    root = Path(__file__).resolve().parents[1]
    text = (root / "docker" / "Dockerfile.rna").read_text(encoding="utf-8")
    idx_annoy = text.find("python-annoy")
    idx_yml = text.find("/opt/aria/envs/aria-rna-env.yml")
    assert idx_annoy != -1 and idx_yml != -1
    assert idx_annoy < idx_yml
    assert "-c conda-forge python-annoy" in text
