"""Wave 6D remediation tests: Domain 3 (Multiome) + Domain 4 (Hi-C/Micro-C) polish."""

from __future__ import annotations

import json
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


def test_d3_h4_m1d_integration_adapter_padj_and_max_regulon_padj():
    from aria.agents.narrative.hypothesis.adapters.integration import integration_evidence

    ledger = {
        "entries": [
            {"node_id": "ledger://integration/peak_to_gene", "status": "ran"},
            {"node_id": "ledger://integration/eregulons", "status": "ran"},
        ]
    }
    # 1. Active FDR cutoff rejects missing padj rows; reads max_regulon_padj (0.01) over default 0.05.
    res_strict = {
        "integration_agent": {
            "status": "done",
            "findings": {
                "peak_to_gene": {
                    "status": "success",
                    "max_padj": 0.05,
                    "top_links": [
                        {"gene": "MISSING_PADJ", "peak": "chr1:100-200", "correlation": 0.9},
                        {"gene": "SIG_GENE", "peak": "chr1:300-400", "correlation": 0.5, "padj": 0.02},
                    ],
                },
                "eregulons": {
                    "status": "success",
                    "max_regulon_padj": 0.01,
                    "top_eregulons": [
                        {
                            "regulon_id": "NO_PADJ_pos",
                            "tf": "NOPADJ",
                            "sign_class": "pos_tf_pos_p2g",
                            "mean_r_tf_gene": 0.8,
                        },
                        {
                            "regulon_id": "ABOVE_CUT_pos",
                            "tf": "ABOVECUT",
                            "sign_class": "pos_tf_pos_p2g",
                            "mean_r_tf_gene": 0.7,
                            "regulon_padj": 0.03,  # > max_regulon_padj=0.01!
                        },
                        {
                            "regulon_id": "PASS_TF_pos",
                            "tf": "PASSTF",
                            "sign_class": "pos_tf_pos_p2g",
                            "mean_r_tf_gene": 0.6,
                            "regulon_padj": 0.005,
                        },
                    ],
                },
            },
        }
    }
    sigs = integration_evidence(res_strict, ledger)
    entities = {s.entity for s in sigs}
    assert entities == {"SIG_GENE", "PASSTF"}

    # 2. When FDR is disabled (max_padj=None / background_method="none"), missing padj
    #    sorts AFTER finite padj (as 1.0, NEVER 0.0!).
    res_nofdr = {
        "integration_agent": {
            "status": "done",
            "findings": {
                "peak_to_gene": {
                    "status": "success",
                    "background_method": "none",
                    "max_padj": None,
                    "top_links": [
                        {"gene": "UNSCORED", "peak": "chr1:100-200", "correlation": 0.95},
                        {"gene": "SCORED", "peak": "chr1:300-400", "correlation": 0.40, "padj": 0.04},
                    ],
                },
                "eregulons": {
                    "status": "success",
                    "max_regulon_padj": None,
                    "top_eregulons": [
                        {
                            "regulon_id": "TF_NOPADJ",
                            "tf": "TFNOPADJ",
                            "sign_class": "pos_tf_pos_p2g",
                            "mean_r_tf_gene": 0.9,
                        },
                        {
                            "regulon_id": "TF_SCORED",
                            "tf": "TFSCORED",
                            "sign_class": "pos_tf_pos_p2g",
                            "mean_r_tf_gene": 0.4,
                            "regulon_padj": 0.02,
                        },
                    ],
                },
            },
        }
    }
    sigs_nf = integration_evidence(res_nofdr, ledger)
    p2g_order = [s.entity for s in sigs_nf if s.audited_node_ref == "ledger://integration/peak_to_gene"]
    ereg_order = [s.entity for s in sigs_nf if s.audited_node_ref == "ledger://integration/eregulons"]
    assert p2g_order == ["SCORED", "UNSCORED"]
    assert ereg_order == ["TFSCORED", "TFNOPADJ"]


def test_d3_h5_h14_m4_m7_integration_agent_wiring(tmp_path: Path):
    from aria.agents.integration_agent import IntegrationAgent
    from aria.bus.message_bus import Confidence
    from aria.memory.memory import ARIAMemory

    mem = ARIAMemory(db_path=str(tmp_path / "integ.sqlite"))
    agent = IntegrationAgent(mem, MagicMock())
    agent.env = MagicMock()

    # 1. D3-H5: raw_multiome_fastq reads exp_ctx["multiome_contract"] and top-level rna_fastqs/atac_fastqs.
    captured_asm: list[dict] = []
    out_h5mu = tmp_path / "asm.h5mu"
    out_h5mu.write_bytes(b"ok")

    def _fake_run(*, stack, script_path, params):
        captured_asm.append(dict(params))
        return {"status": "success", "output_path": str(out_h5mu)}

    agent.env.run_in_stack.side_effect = _fake_run
    exp_ctx = {
        "output_dir": str(tmp_path / "out"),
        "multiome_contract": {
            "object_type": "raw_multiome_fastq",
            "platform": "share_seq",
            "rna_fastqs": [str(tmp_path / "r1.fq.gz"), str(tmp_path / "r2.fq.gz")],
            "atac_fastqs": [str(tmp_path / "a1.fq.gz"), str(tmp_path / "a2.fq.gz"), str(tmp_path / "ab.fq.gz")],
        },
    }
    gate = {"allowed": True, "evidence_paths": [], "contract": {"object_type": "raw_multiome_fastq"}}
    res = agent._run_multiome_preprocessing("exp1", exp_ctx, {}, gate, {})
    assert res is not None and res["status"] == "success"
    assert captured_asm[-1]["rna_fastqs"] == exp_ctx["multiome_contract"]["rna_fastqs"]
    assert captured_asm[-1]["atac_fastqs"] == exp_ctx["multiome_contract"]["atac_fastqs"]
    assert captured_asm[-1]["platform"] == "share_seq"

    # 2. D3-H14: paired-matrix assembly prefers explicit peak matrix even when fragments.tsv.gz sidecar present.
    captured_asm.clear()
    rna_p = tmp_path / "rna.h5ad"
    atac_mat = tmp_path / "atac_peaks.h5ad"
    frags_p = tmp_path / "atac_fragments.tsv.gz"
    for p in (rna_p, atac_mat, frags_p):
        p.write_bytes(b"x")
    gate_mat = {
        "allowed": True,
        "evidence_paths": [],
        "contract": {"object_type": "explicit_same_cell", "cell_namespace": "shared_barcodes"},
    }
    res2 = agent._run_multiome_preprocessing(
        "exp2",
        {"output_dir": str(tmp_path / "out2")},
        {},
        gate_mat,
        {"scRNA": [str(rna_p)], "scATAC": [str(frags_p), str(atac_mat)]},
    )
    assert res2 is not None and res2["status"] == "success"
    assert captured_asm[-1]["atac_path"] == str(atac_mat)
    assert captured_asm[-1]["fragments_path"] == str(frags_p)

    # 3. D3-M4: _evaluate_wnn_weights checks ATAC dominance symmetrically and downgrades >0.90 to LOW.
    published: list[tuple[dict, Confidence]] = []
    agent.publish_finding = lambda exp_id, payload, conf: published.append((payload, conf))
    # Force DebateCouncil failure path to also verify fallback downgrades >0.90 to LOW!
    agent.llm = None
    wnn_atac_dom = {"mean_rna_weight": 0.06, "mean_atac_weight": 0.94, "n_cells": 500}
    agent._evaluate_wnn_weights("exp3", wnn_atac_dom, {})
    assert published[-1][1] == Confidence.LOW

    wnn_atac_mod = {"mean_rna_weight": 0.12, "mean_atac_weight": 0.88, "n_cells": 500}
    agent._evaluate_wnn_weights("exp3", wnn_atac_mod, {})
    assert published[-1][1] == Confidence.MEDIUM

    # 4. D3-M7: _interpret_mofa_factors forwards batch fields; _resolve_conflicts excludes batch-only factors from depth/phase block.
    from aria.agents import debate_council as dc_mod

    debate_ev_captured: list[dict] = []

    class _FakeCouncil:
        def __init__(self, llm, max_rounds=2):
            pass

        def resolve(self, *, topic, initial_claim, evidence, biological_context):
            debate_ev_captured.append(dict(evidence))
            return SimpleNamespace(
                consensus="ok",
                verdict=SimpleNamespace(value="supported"),
                confidence="medium",
                limitations=[],
            )

    orig_dc = dc_mod.DebateCouncil
    dc_mod.DebateCouncil = _FakeCouncil
    try:
        mofa_res = {
            "status": "success",
            "top_factors": [{"factor_id": 1}, {"factor_id": 2}],
            "shared_factors": [{"factor_id": 1}],
            "excluded_factors": [
                {
                    "factor_id": 2,
                    "batch_flag": True,
                    "technical_depth_flag": False,
                    "cell_cycle_flag": False,
                    "exclusion_reasons": ["batch_or_donor_confounded:donor_id(eta2=0.810)"],
                }
            ],
            "batch_factor_flag": True,
            "batch_factor_ids": [2],
            "batch_key": "donor_id",
            "batch_condition_confounded": True,
            "variance_explained": {"rna": 0.25, "atac": 0.15},
            "factor1_top_features": ["G1"],
        }
        agent._interpret_mofa_factors("exp4", mofa_res, {})
        assert len(debate_ev_captured) == 1
        ev = debate_ev_captured[0]
        assert ev["batch_factor_flag"] is True
        assert ev["batch_factor_ids"] == [2]
        assert ev["batch_key"] == "donor_id"
        assert ev["batch_condition_confounded"] is True

        conflict_out = agent._resolve_conflicts("exp4", {}, {}, {"mofa": mofa_res})
        ctypes = [c["type"] for c in conflict_out["conflicts"]]
        assert "mofa_batch_confounded_factors" in ctypes
        assert "mofa_confounder_factor_excluded" not in ctypes
    finally:
        dc_mod.DebateCouncil = orig_dc


def test_d3_h11_h13_m6_mudata_cache_dup_barcodes_and_mofa_sparse_atac(tmp_path: Path):
    ad = pytest.importorskip("anndata")
    from aria.scripts.integration_mofa import _prepare_view_matrices, integration_mofa
    from aria.scripts.integration_peak2gene import integration_peak2gene
    from aria.utils.mudata_io import _stage_cache, _stage_complete

    # 1. mudata_io._stage_cache fails closed when receipt_sha256 is empty/missing.
    inp = tmp_path / "in.txt"
    inp.write_text("hello", encoding="utf-8")
    out_file = tmp_path / "out.csv"
    out_file.write_text("a,b\n1,2\n", encoding="utf-8")
    stage_dir = tmp_path / "stage_cache"
    cached0, st = _stage_cache(
        "test_stage_w6d",
        stage_dir,
        [("in", str(inp))],
        {"p": 1},
        resume=True,
        required=("output_path",),
    )
    assert cached0 is None
    _stage_complete({"status": "success", "output_path": str(out_file)}, st)
    cached1, _ = _stage_cache(
        "test_stage_w6d",
        stage_dir,
        [("in", str(inp))],
        {"p": 1},
        resume=True,
        required=("output_path",),
    )
    assert cached1 is not None and cached1["resumed"] is True

    # Strip receipt_sha256 from manifest -> MUST fail closed (return None)!
    manifest_p = stage_dir / "test_stage_w6d.stage_manifest.json"
    mpayload = json.loads(manifest_p.read_text(encoding="utf-8"))
    mpayload["extra"].pop("receipt_sha256", None)
    manifest_p.write_text(json.dumps(mpayload), encoding="utf-8")
    cached_missing_sha, _ = _stage_cache(
        "test_stage_w6d",
        stage_dir,
        [("in", str(inp))],
        {"p": 1},
        resume=True,
        required=("output_path",),
    )
    assert cached_missing_sha is None

    # 2. Split-file duplicate barcodes rejected in integration_peak2gene and integration_mofa.
    rna_dup = ad.AnnData(
        X=sp.csr_matrix(np.ones((4, 3), dtype=np.int64)),
        obs=pd.DataFrame(index=["BC1", "BC1", "BC2", "BC3"]),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    atac_ok = ad.AnnData(
        X=sp.csr_matrix(np.ones((4, 3), dtype=np.int64)),
        obs=pd.DataFrame(index=["BC1", "BC2", "BC3", "BC4"]),
        var=pd.DataFrame(index=["chr1:100-200", "chr1:300-400", "chr1:500-600"]),
    )
    rna_p = tmp_path / "rna_dup.h5ad"
    atac_p = tmp_path / "atac_ok.h5ad"
    rna_dup.write_h5ad(rna_p)
    atac_ok.write_h5ad(atac_p)

    gtf_p = tmp_path / "genes.gtf"
    gtf_p.write_text('chr1\tsrc\tgene\t1000\t2000\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n')

    try:
        import scanpy  # noqa: F401
        res_p2g = integration_peak2gene(
            {
                "rna_files": [str(rna_p)],
                "atac_files": [str(atac_p)],
                "gtf_path": str(gtf_p),
                "output_dir": str(tmp_path / "p2g_dup"),
                "resume": False,
            }
        )
        assert res_p2g["status"] == "error"
        assert res_p2g["error_type"] == "DuplicateBarcodes"
    except ImportError:
        pass

    try:
        import mofapy2  # noqa: F401
        import scanpy  # noqa: F401
        res_mofa = integration_mofa(
            {
                "modalities": {"scRNA": [str(rna_p)], "scATAC": [str(atac_p)]},
                "output_dir": str(tmp_path / "mofa_dup"),
                "resume": False,
            }
        )
        assert res_mofa["status"] == "error"
        assert res_mofa["error_type"] == "DuplicateBarcodes"
    except ImportError:
        pass

    # 3. D3-M6: _prepare_view_matrices on pre-normalized sparse ATAC selects n_top_peaks sparsely.
    try:
        import scanpy  # noqa: F401
        rng = np.random.default_rng(7)
        dense = rng.uniform(0.0, 2.0, size=(20, 12)).astype(np.float64)
        # Give columns 0..7 strictly higher variance than 8..11
        dense[:, :8] *= 5.0
        sp_mat = sp.csr_matrix(dense)
        prep = _prepare_view_matrices(
            mod_name="scATAC",
            mod_kind="atac",
            matrix_all=sp_mat,
            var_names=[f"chr1:{i*100}-{i*100+50}" for i in range(12)],
            train_indices=np.arange(20),
            is_raw_counts=False,
            n_top_genes=3,
            n_top_peaks=7,
            scale_views=True,
        )
        assert prep["n_features_selected"] == 7  # used n_top_peaks=7, NOT n_top_genes=3!
    except ImportError:
        pass


def test_d3_m1_eregulon_core_zero_var_depth_sign_flip_and_small_n_clamp(tmp_path: Path):
    from aria.utils.eregulon_core import (
        compute_candidate_eregulons,
        compute_gc_stratified_motif_enrichment,
    )

    # 1. CMH zero-variance motif gets pvalue == 1.0 (NOT 0.5).
    hits = np.array(
        [
            [False, False, False, False, False, False],  # 0 hits -> var == 0 -> p == 1.0
            [True, True, True, True, True, True],        # all hits -> var == 0 -> p == 1.0
        ],
        dtype=bool,
    )
    fg = np.array([True, True, True, False, False, False], dtype=bool)
    gc = np.array([0.4, 0.5, 0.6, 0.4, 0.5, 0.6], dtype=float)
    cmh = compute_gc_stratified_motif_enrichment(hits, fg, gc, n_gc_bins=2)
    assert np.allclose(cmh["pvalue"], [1.0, 1.0])

    # 2. Small cohort (n_cells=30, min_cells_detected=10) clamps eff_min_cells=6 matching peak2gene_core,
    #    and depth-driven sign reversal flags depth_confounded=True.
    n_cells = 30
    rng = np.random.default_rng(99)
    # Build synthetic FASTA + .fai + MEME motif matching an embedded sequence in chr1 peaks.
    fa_path = tmp_path / "mini.fa"
    motif_seq = "ACGTACGTACGT"
    chrom_seq = ("TTTT" * 500) + motif_seq + ("CCCC" * 500)
    fa_path.write_text(f">chr1\n{chrom_seq}\n", encoding="utf-8")
    (tmp_path / "mini.fa.fai").write_text(f"chr1\t{len(chrom_seq)}\t6\t{len(chrom_seq)}\t{len(chrom_seq)+1}\n")

    meme_path = tmp_path / "mini.meme"
    mat_lines = []
    base_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    for ch in motif_seq:
        probs = [0.01, 0.01, 0.01, 0.01]
        probs[base_idx[ch]] = 0.97
        mat_lines.append(" ".join(f"{p:.2f}" for p in probs))
    meme_path.write_text(
        "MEME version 4\nALPHABET= ACGT\n"
        f"MOTIF M1 TF1\nletter-probability matrix: alength= 4 w= {len(motif_seq)}\n"
        + "\n".join(mat_lines)
        + "\n",
        encoding="utf-8",
    )

    # TF1 detected in 8 of 30 cells (< 10 raw min_cells_detected, >= 6 clamped eff_min_cells!).
    n_genes = 30
    gene_names = ["TF1"] + [f"G{i}" for i in range(1, n_genes)]
    rna = np.zeros((n_cells, n_genes), dtype=np.int64)
    # High library depth in cells 0..7 so TF1 + targets are positively correlated raw,
    # while residualizing depth flips correlation sign!
    for c in range(8):
        rna[c, 0] = 20 + c * 2
        for g in range(1, 7):
            rna[c, g] = 30 - c  # inversely ordered within depth group -> negative partial r!
    # Background genes expressed across cells so library size tracks cells 0..7 vs 8..29.
    for c in range(n_cells):
        base_lib = 200 if c < 8 else 10
        rna[c, 7:] = rng.poisson(base_lib, size=n_genes - 7) + 1

    peak_names = [f"chr1:{1990}-{2020}"] + [f"chr1:{100 + i*20}-{120 + i*20}" for i in range(1, 10)]
    atac = np.zeros((n_cells, len(peak_names)), dtype=np.int64)
    atac[:8, 0] = 5
    atac[:, 1:] = 2

    tss_records = {g: ("chr1", 2000, "+") for g in gene_names}
    precomputed_p2g = {
        "links": [
            {
                "gene": f"G{i}",
                "peak": peak_names[0],
                "correlation": 0.7,
                "padj": 0.001,
                "distance_bp": 100,
                "region_class": "promoter_proximal",
            }
            for i in range(1, 7)
        ],
        "n_tested": 20,
        "n_fdr_significant": 6,
        "background_match": "precomputed_csv",
    }
    res = compute_candidate_eregulons(
        rna_matrix=sp.csr_matrix(rna),
        rna_gene_names=gene_names,
        atac_matrix=sp.csr_matrix(atac),
        atac_peak_names=peak_names,
        gene_tss_records=tss_records,
        meme_path=meme_path,
        fasta_path=fa_path,
        min_cells_detected=10,
        min_target_genes=3,
        max_regulon_padj=1.0,
        n_null_permutations=20,
        motif_pvalue=1e-2,
        scanner_backend="numpy",
        precomputed_p2g=precomputed_p2g,
    )
    # TF1 was retained because eff_min_cells clamped 10 -> 6 for n_cells=30!
    assert res["n_motif_tf_mappings"] >= 1
    assert res["n_eregulons"] >= 1
    reg = res["eregulons"][0]
    assert reg["mean_r_tf_gene"] > 0 and reg["mean_r_tf_gene_partial_depth"] < 0
    assert reg["depth_confounded"] is True


def test_d3_m2_m3b_m8_d2_m6_narrator_plots_and_cohort_union_provenance(tmp_path: Path, monkeypatch):
    import h5py
    import aria.agents._narrative_integration as ni
    from aria.agents.narrative.narrators.integration import IntegrationNarrator
    from aria.utils.scatac_cohort_union import merge_scatac_matrices

    # 1. IntegrationNarrator: float distance_kb=0.5, gated BH-FDR & rna_normalization, n_other_sign_classes.
    narr = IntegrationNarrator()
    p2g_nofdr = {
        "status": "success",
        "n_cells": 200,
        "distance_kb": 0.5,
        "min_corr": 0.25,
        "max_padj": None,
        "background_method": "none",
        "rna_normalization": "pre_normalized_input",
        "atac_normalization": "raw_counts",
        "n_tested": 100,
        "n_links": 12,
        "n_positive_correlations": 10,
        "n_negative_correlations": 2,
        "n_promoter_proximal": 5,
        "n_distal_cis": 7,
    }
    p2g_blk = narr._peak_to_gene_block(p2g_nofdr, {})
    assert p2g_blk.metrics["distance_kb"] == 0.5
    assert "0.5 kb" in p2g_blk.claim
    assert "BH-FDR" not in p2g_blk.claim

    methods_txt = " ".join(
        narr.methods("integration_agent", {"status": "done", "findings": {"peak_to_gene": p2g_nofdr}})
    )
    assert "0.5 kb" in methods_txt
    assert "FDR <= None" not in methods_txt
    assert "pre-normalized gene expression" in methods_txt
    assert "log1p-CP10k" not in methods_txt

    ereg_other = {
        "status": "success",
        "n_shared_cells": 600,
        "n_candidate_modules_tested": 50,
        "n_eregulons": 5,
        "n_triplets": 30,
        "n_unique_tfs": 4,
        "n_pos_tf_pos_p2g": 2,
        "n_neg_tf_pos_p2g": 1,
        "n_other_sign_classes": 2,
        "n_paralog_ambiguous": 0,
        "n_depth_confounded": 0,
    }
    ereg_blk = narr._eregulons_block(ereg_other, {})
    assert ereg_blk.metrics["n_other_sign_classes"] == 2
    assert "2 other sign-class" in ereg_blk.claim
    assert any(ev.label == "other sign-class modules" and ev.value == 2 for ev in ereg_blk.evidence)

    # 2. _narrative_integration: _read_wnn_h5mu matches "GEX"/"Peaks" via _match_modality + IQR _label_str(cl_key).
    h5_p = tmp_path / "custom_keys.h5mu"
    str_dt = h5py.string_dtype(encoding="utf-8")
    bcs = [f"c{i}" for i in range(10)]
    with h5py.File(h5_p, "w") as f:
        obs = f.create_group("obs")
        obs.attrs["_index"] = "_index"
        obs.create_dataset("_index", data=np.asarray(bcs, dtype=object), dtype=str_dt)
        obs.create_dataset("leiden_wnn", data=np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1], dtype=np.int64))
        obs.create_dataset("rna_weight", data=np.linspace(0.3, 0.8, 10))
        obsm = f.create_group("obsm")
        obsm.create_dataset("X_wnn_umap", data=np.zeros((10, 2), dtype=np.float32))
        mod = f.create_group("mod")
        for mk in ("GEX", "Peaks"):
            mg = mod.create_group(mk)
            mobs = mg.create_group("obs")
            mobs.attrs["_index"] = "_index"
            mobs.create_dataset("_index", data=np.asarray(bcs, dtype=object), dtype=str_dt)
            mobsm = mg.create_group("obsm")
            mobsm.create_dataset("X_umap", data=np.ones((10, 2), dtype=np.float32))
    bundle = ni._read_wnn_h5mu(h5_p)
    assert bundle is not None
    assert bundle["rna_umap"] is not None and bundle["atac_umap"] is not None

    # Render weights with float/bytes cluster keys (0.0, 1.0) -> IQR whiskers must match via _label_str!
    captured_figs = []
    real_save = ni._save_dual

    def _cap(fig, p, dpi=220):
        captured_figs.append(fig)
        return real_save(fig, p, dpi=60)

    monkeypatch.setattr(ni, "_save_dual", _cap)
    w_out = ni.render_wnn_modality_weights(
        {
            "status": "success",
            "output_path": str(h5_p),
            "mean_rna_weight": 0.55,
            "mean_atac_weight": 0.45,
            "per_cluster_modality_weights": {
                0.0: {"n_cells": 5, "mean_rna_weight": 0.4, "mean_atac_weight": 0.6},
                1.0: {"n_cells": 5, "mean_rna_weight": 0.7, "mean_atac_weight": 0.3},
            },
        },
        tmp_path,
    )
    assert w_out is not None
    # Verify IQR legend entry was added (meaning iqr_by_cluster matched!)
    leg_texts = [t.get_text() for t in captured_figs[-1].axes[0].get_legend().get_texts()]
    assert any("IQR" in txt for txt in leg_texts)

    # 3. render_mofa_variance_explained renders 3 views dynamically.
    mofa_3v = {
        "status": "success",
        "views": {"rna": {"kind": "rna"}, "atac": {"kind": "atac"}, "ADT": {"kind": "protein"}},
        "variance_explained": {"rna": 0.20, "atac": 0.15, "ADT": 0.30},
        "variance_explained_additive": {"rna": 0.22, "atac": 0.16, "ADT": 0.31},
        "top_factors": [
            {"factor_id": 1, "category": "shared", "r2_by_view": {"rna": 0.12, "atac": 0.08, "ADT": 0.18}},
            {"factor_id": 2, "category": "shared", "r2_by_view": {"rna": 0.08, "atac": 0.07, "ADT": 0.12}},
        ],
    }
    assert ni.render_mofa_variance_explained(mofa_3v, tmp_path) is not None
    ax_hm_3v = captured_figs[-1].axes[0]
    xt_3v = [t.get_text() for t in ax_hm_3v.get_xticklabels()]
    assert len(xt_3v) == 3
    assert any("ADT" in t for t in xt_3v)

    # 4. render_regulatory_triangulation avoids "total BH-FDR=" when background_method == "none".
    p2g_tri = {
        "status": "success",
        "n_links": 2,
        "background_method": "none",
        "max_padj": None,
        "top_links": [
            {"gene": "G1", "correlation": 0.5, "distance_bp": 500, "region_class": "promoter_proximal"},
            {"gene": "G2", "correlation": 0.4, "distance_bp": 15000, "region_class": "distal_cis"},
        ],
    }
    assert ni.render_regulatory_triangulation(p2g_tri, None, tmp_path) is not None
    title_txt = captured_figs[-1].axes[0].get_title()
    assert "BH-FDR" not in title_txt
    assert "total=2" in title_txt

    # 5. D2-M6: scatac_cohort_union writes structured dict to merged.uns["peak_provenance"].
    ad = pytest.importorskip("anndata")
    p1 = tmp_path / "s1.h5ad"
    p2 = tmp_path / "s2.h5ad"
    ad.AnnData(
        X=sp.csr_matrix([[3, 2]], dtype=np.int64),
        obs=pd.DataFrame(index=["C1"]),
        var=pd.DataFrame(index=["chr1:100-200", "chr1:300-400"]),
    ).write_h5ad(p1)
    ad.AnnData(
        X=sp.csr_matrix([[4]], dtype=np.int64),
        obs=pd.DataFrame(index=["C2"]),
        var=pd.DataFrame(index=["chr1:150-250"]),
    ).write_h5ad(p2)
    out_union = tmp_path / "union.h5ad"
    u_res = merge_scatac_matrices([str(p1), str(p2)], str(out_union))
    assert u_res["status"] == "success", u_res
    assert u_res["peak_provenance"] == "cohort_interval_union_consensus"
    merged_ad = ad.read_h5ad(out_union)
    assert merged_ad.uns["peak_provenance"] == {
        "method": "overlap_unified",
        "source": "cohort_interval_union_consensus",
        "n_samples": 2,
    }


# ═══════════════════════════════════════════════════════════════════════════
# PART 2: DOMAIN 4 (HI-C / MICRO-C)
# ═══════════════════════════════════════════════════════════════════════════


def test_d4_h5_h6_hic_adapter_low_replication_qval_and_disk_fallthrough(tmp_path: Path):
    from aria.agents.narrative.hypothesis.adapters.hic import _low_replication, hic_evidence
    from aria.agents.narrative.hypothesis.caveats import LOW_REPLICATION

    # 1. D4-H5: descriptive_only=True overrides n_samples >= 2; loops without n_samples >= 2 keeps LOW_REPLICATION.
    assert _low_replication(
        {
            "qc": {"n_samples": 4},
            "topology": {
                "loops": {"status": "success", "n_samples": 4},
                "differential_3d": {"descriptive_only": True},
            },
        }
    ) is True
    assert _low_replication(
        {
            "qc": {"n_samples": 4},
            "topology": {"loops": {"status": "success"}},  # default n_samples=1
        }
    ) is True
    assert _low_replication(
        {
            "qc": {"n_samples": 4},
            "topology": {
                "loops": {"status": "success", "n_samples": 4},
                "differential_3d": {"descriptive_only": False},
            },
        }
    ) is False

    # 2. D4-H6: NaN/missing qval rejected under active fdr_threshold; never substitutes count into value;
    #    falls through to annotated_loops_path when top_loops yields < _MAX_LOOP_GENES.
    ann_tsv = tmp_path / "annotated_loops.tsv"
    pd.DataFrame(
        [
            {
                "chrom1": "chr1", "start1": 10000, "end1": 20000,
                "chrom2": "chr1", "start2": 90000, "end2": 100000,
                "la_exp.donut.qval": 0.002, "count": 50,
                "anchor1_genes": "DISK_GENE1", "anchor2_genes": "DISK_GENE2",
            },
            {
                "chrom1": "chr1", "start1": 30000, "end1": 40000,
                "chrom2": "chr1", "start2": 110000, "end2": 120000,
                "la_exp.donut.qval": float("nan"), "count": 999,
                "anchor1_genes": "NAN_QVAL_GENE", "anchor2_genes": "",
            },
        ]
    ).to_csv(ann_tsv, sep="\t", index=False)

    ledger = {"entries": [{"node_id": "ledger://hic/loops", "status": "ran"}]}
    res = {
        "genome_arch_agent": {
            "status": "done",
            "findings": {
                "qc": {"n_samples": 4},
                "topology": {
                    "loops": {
                        "status": "success",
                        "ran": True,
                        "fdr_threshold": 0.05,
                        "annotated_loops_path": str(ann_tsv),
                        "top_loops": [
                            {
                                "chrom1": "chr1", "start1": 0, "end1": 10000,
                                "chrom2": "chr1", "start2": 50000, "end2": 60000,
                                "qval": float("nan"), "count": 123,
                                "anchor1_genes": "BAD_INLINE",
                            },
                            {
                                "chrom1": "chr1", "start1": 0, "end1": 10000,
                                "chrom2": "chr1", "start2": 60000, "end2": 70000,
                                "qval": 0.001, "count": 88,
                                "anchor1_genes": "TOP_GENE",
                            },
                        ],
                    }
                },
            },
        }
    }
    sigs = hic_evidence(res, run_ledger=ledger)
    by_ent = {s.entity: s for s in sigs}
    assert set(by_ent) == {"TOP_GENE", "DISK_GENE1", "DISK_GENE2"}
    assert "BAD_INLINE" not in by_ent
    assert "NAN_QVAL_GENE" not in by_ent
    assert by_ent["TOP_GENE"].value == pytest.approx(0.001)
    assert by_ent["DISK_GENE1"].value == pytest.approx(0.002)
    assert LOW_REPLICATION in by_ent["TOP_GENE"].caveats_inherited

    # When fdr_threshold is None and qval is NaN, value MUST be None (never count!).
    res_nofdr = {
        "genome_arch_agent": {
            "status": "done",
            "findings": {
                "topology": {
                    "loops": {
                        "status": "success",
                        "ran": True,
                        "top_loops": [
                            {
                                "chrom1": "chr1", "start1": 0, "end1": 10000,
                                "chrom2": "chr1", "start2": 50000, "end2": 60000,
                                "qval": float("nan"), "count": 777,
                                "anchor1_genes": "MYC",
                            }
                        ],
                    }
                }
            },
        }
    }
    sigs_nf = hic_evidence(res_nofdr, run_ledger=ledger)
    assert len(sigs_nf) == 1
    assert sigs_nf[0].entity == "MYC"
    assert sigs_nf[0].value is None


def test_d4_h10_m9_hic_narrator_diagnostics_tables_tads_fig_and_captions(tmp_path: Path):
    from aria.agents.narrative.narrators.hic import HiCNarrator

    narr = HiCNarrator()
    findings = {
        "fastq_alignment": {
            "status": "failed",
            "ran": False,
            "stage": "align",
            "details": {"error_type": "BwaFailed", "details": "index missing"},
        },
        "preprocessing": {
            "status": "skipped",
            "ran": False,
            "reason": "already_mcool",
        },
        "qc": {
            "status": "success",
            "n_samples": 2,
            "n_valid_pairs": 50000,
            "cis_trans_ratio": 0.72,
            "balanced_files": [],  # raw fallback!
            "ps_curve_tsv": str(tmp_path / "s1.expected_cis.tsv"),
            "per_file": [
                {"sample_id": "s1", "ps_curve_tsv": str(tmp_path / "s1.expected_cis.tsv"), "status": "success"},
                {"sample_id": "s2", "ps_curve_tsv": str(tmp_path / "s2.expected_cis.tsv"), "status": "success"},
            ],
        },
        "figures": {
            "hic_contact_map": str(tmp_path / "map.png"),
            "hic_compartment_insulation_track": str(tmp_path / "track.png"),
        },
        "topology": {
            "compartments": {
                "status": "skipped",
                "ran": False,
                "reason": "too_coarse",
                "compartments_bedgraph": str(tmp_path / "comp.bedGraph"),
            },
            "tads": {
                "status": "success",
                "n_tads": 42,
                "median_size_kb": 180.0,
            },
            "loops": {
                "status": "success",
                "n_loops": 10,
                "resolution": 10000,
                "loop_class_counts": {"promoter_enhancer": 4, "promoter_promoter": 2, "distal_distal": 4},
            },
            "differential_3d": {
                "status": "skipped",
                "ran": False,
                "descriptive_only": True,
                "differential_compartments_csv": str(tmp_path / "top_comp.csv"),
                "differential_insulation_csv": str(tmp_path / "top_ins.csv"),
                "comparisons": [
                    {
                        "contrast": "KO_vs_WT",
                        "status": "skipped",
                        "ran": False,
                        "descriptive_only": True,
                        "descriptive_csv": str(tmp_path / "desc_KO_vs_WT.csv"),
                    }
                ],
            },
        },
    }
    blocks = {b.id: b for b in narr.collect("genome_arch_agent", {"findings": findings})}
    assert "hic.fastq_alignment" in blocks
    assert blocks["hic.fastq_alignment"].status == "error"
    assert "index missing" in (blocks["hic.fastq_alignment"].error or "")
    assert "hic.preprocessing" in blocks
    assert blocks["hic.preprocessing"].status == "limitation"

    # Contact map caption must NOT say "ICE-balanced" when balanced_files is empty!
    qc_map_fig = [f for f in blocks["hic.qc"].figures if f["id"] == "hic.hic_contact_map"][0]
    assert "ICE-balanced" not in qc_map_fig["caption"]
    assert "unbalanced" in qc_map_fig["caption"].lower()

    # TADs block attaches hic_compartment_insulation_track when compartments != ok!
    assert any(f["id"] == "hic.hic_compartment_insulation_track" for f in blocks["hic.tads"].figures)

    # Loops block labels promoter_enhancer accurately!
    loop_ev_labels = [ev.label for ev in blocks["hic.loops"].evidence]
    assert "Promoter-enhancer loop candidates" in loop_ev_labels
    assert "Promoter-anchored loop candidates" not in loop_ev_labels

    # tables() registers compartments_bedgraph, per_file[1:] ps_curve_tsv, descriptive_csv, and top-level diff CSVs.
    tbl_paths = {t["path"] for t in narr.tables("genome_arch_agent", {"findings": findings})}
    assert str(tmp_path / "comp.bedGraph") in tbl_paths
    assert str(tmp_path / "s2.expected_cis.tsv") in tbl_paths
    assert str(tmp_path / "desc_KO_vs_WT.csv") in tbl_paths
    assert str(tmp_path / "top_comp.csv") in tbl_paths
    assert str(tmp_path / "top_ins.csv") in tbl_paths


def test_d4_m1_m4_m6_m7_m8_m10a_hic_unit_checks(tmp_path: Path, monkeypatch):
    from aria.agents.genome_arch_agent import GenomeArchAgent
    from aria.memory.memory import ARIAMemory
    from aria.scripts import hic_figure as hf
    from aria.scripts import hic_preprocess as hp
    from aria.scripts import hic_qc_and_balance as hq
    from aria.scripts import hic_topology as ht

    # 1. D4-M4: _measure_cooler_pixels accumulates fractional float counts in float64 without truncating to 0.
    bins_df = pd.DataFrame(
        {"chrom": ["chr1", "chr1", "chr2"], "start": [0, 30_000, 0]}
    )
    pix_df = pd.DataFrame(
        {"bin1_id": [0, 0], "bin2_id": [1, 2], "count": [0.4, 0.2]}
    )

    class _Sub:
        def __init__(self, fn, length):
            self._fn = fn
            self._len = length

        def __getitem__(self, sl):
            return self._fn(sl)

        def __len__(self):
            return self._len

    class _FloatClr:
        info = {"nnz": 2, "nbins": 3}

        def bins(self):
            return _Sub(lambda sl: bins_df, len(bins_df))

        def pixels(self):
            return _Sub(lambda sl: pix_df.iloc[sl], len(pix_df))

    stats = hq._measure_cooler_pixels(_FloatClr())
    # Exact float ratio 0.4 / 0.6 == 0.666667 (previously int64 truncated 0.4 -> 0 and gave 0.0!)
    assert stats["cis_trans_ratio"] == pytest.approx(2.0 / 3.0, abs=1e-4)
    assert stats["long_range_cis"] == pytest.approx(2.0 / 3.0, abs=1e-4)

    # 2. D4-M1: _compute_ps_decay masks NaN sum_col before grouping so denominator isn't inflated.
    import sys

    fake_exp_df = pd.DataFrame(
        {
            "region1": ["chr1", "chr2", "chr1", "chr2", "chr1", "chr2"],
            "dist": [2, 2, 3, 3, 4, 4],
            "dist_bp": [20_000, 20_000, 30_000, 30_000, 40_000, 40_000],
            "n_valid": [10, 100, 10, 100, 10, 100],
            "count.sum": [5.0, float("nan"), 2.0, float("nan"), 1.0, float("nan")],
            "count.avg": [0.5, float("nan"), 0.2, float("nan"), 0.1, float("nan")],
        }
    )

    class _PsClr:
        binsize = 10_000
        chromnames = ["chr1", "chr2"]
        chromsizes = {"chr1": 100_000, "chr2": 100_000}
        info = {"nbins": 20}

        def bins(self):
            return SimpleNamespace(columns=["chrom", "start", "end"])

    fake_eig = SimpleNamespace()
    fake_ins = SimpleNamespace()
    fake_cooltools_api = SimpleNamespace(
        expected=SimpleNamespace(expected_cis=lambda *a, **k: fake_exp_df),
        eigdecomp=fake_eig,
        insulation=fake_ins,
    )
    monkeypatch.setitem(sys.modules, "cooltools", SimpleNamespace(api=fake_cooltools_api))
    monkeypatch.setitem(sys.modules, "cooltools.api", fake_cooltools_api)
    monkeypatch.setitem(sys.modules, "cooltools.api.eigdecomp", fake_eig)
    monkeypatch.setitem(sys.modules, "cooltools.api.insulation", fake_ins)
    ps_out = hq._compute_ps_decay(_PsClr(), tmp_path / "ps.tsv")
    ps_df = pd.read_csv(ps_out["ps_curve_tsv"], sep="\t")
    # At dist=2: valid sum=5.0, valid denom=10 (NOT 110!) -> contact_prob == 0.5!
    assert float(ps_df.loc[ps_df["dist"] == 2, "contact_prob"].iloc[0]) == pytest.approx(0.5)
    assert float(ps_df.loc[ps_df["dist"] == 2, "n_valid"].iloc[0]) == pytest.approx(10.0)

    # 3. D4-M6: _moderated_linear_contrast rejects NaN covariate and supports continuous numeric covariates (1 df).
    rng = np.random.default_rng(11)
    y = rng.normal(size=(10, 6))
    y[:3, :3] += 3.0
    conds = ["KO", "KO", "KO", "WT", "WT", "WT"]
    with pytest.raises(ValueError, match="missing/NaN"):
        ht._moderated_linear_contrast(
            y, conds, "KO", "WT",
            covariate_rows=[{"age": 20.0}, {"age": None}, {"age": 35.0}, {"age": 22.0}, {"age": 30.0}, {"age": 40.0}],
        )
    # 6 unique continuous numeric values must take 1 df (df_resid = 6 - 3 = 3), NOT saturate df!
    fit_cont = ht._moderated_linear_contrast(
        y, conds, "KO", "WT",
        covariate_rows=[{"age": 21.5}, {"age": 34.2}, {"age": 48.1}, {"age": 23.0}, {"age": 36.8}, {"age": 51.4}],
    )
    assert np.all(np.isfinite(fit_cont["padj"]))

    # Empty metadata CSV returns skipped cleanly (no IndexError).
    empty_meta = tmp_path / "empty_meta.csv"
    empty_meta.write_text("", encoding="utf-8")
    f1 = tmp_path / "a.cool"
    f2 = tmp_path / "b.cool"
    f1.write_bytes(b"x")
    f2.write_bytes(b"x")
    monkeypatch.setitem(sys.modules, "cooler", SimpleNamespace())
    skip_res = ht._compute_differential_3d(
        {"files": [str(f1), str(f2)], "sample_metadata": str(empty_meta), "comparisons": [("KO", "WT")]},
        str(tmp_path / "diff_empty"),
    )
    assert skip_res["status"] == "skipped"
    assert skip_res["reason"] == "missing_sample_metadata_or_comparisons"

    # D4-H11: hic_preprocess stage cache handles chromsizes / reference_fasta Path aliases and invalidates on edit.
    cs_file = tmp_path / "hg38.chrom.sizes"
    fa_file = tmp_path / "hg38.fa"
    cs_file.write_text("chr1\t1000000\n", encoding="utf-8")
    fa_file.write_text(">chr1\nACGT\n", encoding="utf-8")
    hp_out = tmp_path / "hp_cache"
    hp_out.mkdir(parents=True, exist_ok=True)
    cached_mcool = hp_out / "s1.mcool"
    cached_mcool.write_bytes(b"fake_mcool")
    hp_params = {
        "files": [str(f1)],
        "sample_ids": ["s1"],
        "chromsizes": cs_file,
        "reference_fasta": fa_file,
        "base_resolution": 10_000,
        "zoomify": True,
        "output_dir": str(hp_out),
        "resume": True,
    }
    hp._save_hic_stage_cache(
        hp_out,
        stage="hic_preprocess",
        version=hp.STAGE_VERSION,
        params={
            "genome": "hg38",
            "base_res": 10_000,
            "zoom_res": sorted(hp.DEFAULT_ZOOM_RESOLUTIONS),
            "zoomify": True,
            "mapq": 30,
            "sample_ids": ["s1"],
            "chrom_sizes": str(cs_file),
            "genome_fasta": str(fa_file),
            "genome_fai": "",
        },
        inputs=[str(f1), str(cs_file), str(fa_file)],
        outputs=[str(cached_mcool)],
        extra={"result": {"status": "success", "cooler_files": [str(cached_mcool)]}},
    )
    hp_resumed = hp.hic_preprocess(hp_params)
    assert hp_resumed["status"] == "success" and hp_resumed.get("resumed") is True
    # Mutating chromsizes invalidates the cache (re-runs and fails on fake .cool input instead of stale resume!)
    cs_file.write_text("chr1\t2000000\nchr2\t500000\n", encoding="utf-8")
    hp_invalid = hp.hic_preprocess(hp_params)
    assert hp_invalid.get("resumed") is not True
    assert hp_invalid["status"] == "error"

    # 4. D4-M7: _pick_display_chrom honors explicit chrom/chromosomes and excludes chrMT/chrEBV/EBV.
    assert hf._pick_display_chrom(["chrMT", "chrEBV", "EBV", "chr2", "chr3"]) == "chr2"
    assert hf._pick_display_chrom(["chr1", "chr2", "chr3"], requested=["chr3"]) == "chr3"
    assert hf._pick_display_chrom(["1", "2", "3"], requested="chr2") == "2"

    # 5. D4-M7 & D4-M10a: GenomeArchAgent passes chromosomes to calibration & figures and populates switch_regions!
    mem = ARIAMemory(db_path=str(tmp_path / "ga.sqlite"))
    agent = GenomeArchAgent(mem, MagicMock())
    agent.env = MagicMock()
    calls: list[tuple[str, dict]] = []

    def _fake_ga_run(*, stack, script_path, params):
        calls.append((script_path, dict(params)))
        if script_path.endswith("hic_inspect.py"):
            return {"status": "success", "files": [str(f1), str(f2)], "available_resolutions": [40_000, 100_000]}
        if script_path.endswith("hic_qc_and_balance.py"):
            return {
                "status": "success",
                "balanced_files": [str(f1), str(f2)],
                "balanced_sample_ids": ["s1", "s2"],
                "ps_curve_tsv": str(tmp_path / "ps.tsv"),
                "per_file": [],
            }
        if script_path.endswith("hic_topology.py"):
            ana = params.get("analysis")
            if ana == "compartments":
                return {"status": "success", "ran": True, "pct_A": 50.0, "pct_B": 50.0, "pc1_validated": True}
            if ana == "insulation_calibration":
                return {"status": "success", "boundary_strengths": {"200000": 0.4}, "calibration_chrom": "chr2"}
            if ana == "tads":
                return {"status": "success", "ran": True, "n_tads": 10}
            if ana == "differential_3d":
                return {
                    "status": "success",
                    "ran": True,
                    "n_ab_switches": 2,
                    "comparisons": [
                        {
                            "contrast": "KO_vs_WT",
                            "ran": True,
                            "top_switches": [
                                {"chrom": "chr2", "start": 0, "end": 100000, "delta_E1": -0.8, "padj": 0.001, "switch_class": "A_to_B"},
                                {"chrom": "chr2", "start": 200000, "end": 300000, "delta_E1": 0.7, "padj": 0.01, "switch_class": "B_to_A"},
                            ],
                        }
                    ],
                }
        if script_path.endswith("hic_figure.py"):
            return {"status": "success", "figures": {}, "figures_svg": {}, "figure_bundles": {}}
        return {"status": "success"}

    agent.env.run_in_stack.side_effect = _fake_ga_run
    debate_sw_captured: list[dict] = []
    from aria.agents import debate_council as dc_mod

    class _FakeCouncilHiC:
        def __init__(self, llm, max_rounds=2):
            pass

        def resolve(self, *, topic, initial_claim, evidence, biological_context):
            debate_sw_captured.append(dict(evidence))
            return SimpleNamespace(consensus="ok", verdict=SimpleNamespace(value="supported"), limitations=[])

    monkeypatch.setattr(dc_mod, "DebateCouncil", _FakeCouncilHiC)
    out = agent.run(
        "exp_ga_w6d",
        {
            "exp_context": {
                "genome": "hg38",
                "chromosomes": ["chr2"],
                "modalities": {"HiC": [str(f1), str(f2)]},
                "sample_ids": ["s1", "s2"],
                "sample_metadata": [
                    {"sample_id": "s1", "condition": "KO", "replicate": "r1"},
                    {"sample_id": "s2", "condition": "WT", "replicate": "r1"},
                ],
                "comparisons": [("KO", "WT")],
                "output_dir": str(tmp_path / "ga_out"),
            },
            "biological_intent": {"user_question": "Analyze TADs and compartments"},
        },
    )
    assert out["status"] == "done"
    cal_calls = [p for s, p in calls if s.endswith("hic_topology.py") and p.get("analysis") == "insulation_calibration"]
    assert cal_calls[0]["chromosomes"] == ["chr2"]
    fig_calls = [p for s, p in calls if s.endswith("hic_figure.py")]
    assert fig_calls[0]["chromosomes"] == ["chr2"]
    assert len(debate_sw_captured) == 1
    assert len(debate_sw_captured[0]["switch_regions"]) == 2
    assert debate_sw_captured[0]["switch_regions"][0]["chrom"] == "chr2"

    # 6. D4-M8: _convert_bam_to_pairs removes stale stats_path before pairtools dedup.
    stale_stats = tmp_path / "s1.pairtools.stats.tsv"
    stale_stats.write_text("stale\t999\n", encoding="utf-8")
    seen_stats_existed: list[bool] = []

    def _fake_subproc(cmd, **kw):
        if "dedup" in cmd:
            seen_stats_existed.append(stale_stats.exists())
            stale_stats.write_text("total\t10\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(hp, "_stack_bin", lambda name: "/usr/bin/pairtools")
    monkeypatch.setattr(hp.subprocess, "run", _fake_subproc)
    _dp, pt_stats = hp._convert_bam_to_pairs(
        tmp_path / "in.bam", tmp_path / "chrom.sizes", tmp_path, "s1", mapq=30, threads=1
    )
    assert seen_stats_existed == [False]
    assert pt_stats == {"total": 10}


def test_d4_live_hic_env_multi_res_convergence_cache_and_chrom_figure(tmp_path: Path, monkeypatch):
    cooler = pytest.importorskip("cooler")
    pytest.importorskip("cooltools")
    from aria.scripts import hic_figure as hf
    from aria.scripts import hic_qc_and_balance as hq
    from aria.scripts import hic_topology as ht

    # Build a 2-chrom .cool with chr1 + chr2
    bs = 40_000
    n_per = 15
    bins = pd.DataFrame(
        [(ch, b * bs, (b + 1) * bs) for ch in ("chr1", "chr2") for b in range(n_per)],
        columns=["chrom", "start", "end"],
    )
    pix_rows = []
    for ch_i in range(2):
        off = ch_i * n_per
        for i in range(n_per):
            for j in range(i, min(n_per, i + 6)):
                pix_rows.append((off + i, off + j, max(1, 20 - 2 * (j - i))))
    cool_p = tmp_path / "two_chrom.cool"
    cooler.create_cooler(str(cool_p), bins, pd.DataFrame(pix_rows, columns=["bin1_id", "bin2_id", "count"]))

    # 1. D4-M3: _ensure_balanced_cooler fails when cooler.balance_cooler reports converged=False.
    orig_bal = cooler.balance_cooler

    def _nonconv_bal(clr, *a, **kw):
        w, st = orig_bal(clr, *a, **kw)
        st = dict(st or {})
        st["converged"] = False
        return w, st

    monkeypatch.setattr(cooler, "balance_cooler", _nonconv_bal)
    with pytest.raises(RuntimeError, match="did not converge"):
        ht._ensure_balanced_cooler(str(cool_p), bs, tmp_path / "bal_fail")
    monkeypatch.setattr(cooler, "balance_cooler", orig_bal)

    # 2. D4-H11: hic_qc_and_balance stage cache invalidates when min_valid_pairs changes.
    qc_dir = tmp_path / "qc_cache"
    r1 = hq.hic_qc_and_balance(
        {
            "files": [str(cool_p)],
            "sample_ids": ["s1"],
            "resolution": bs,
            "min_valid_pairs": 100,
            "output_dir": str(qc_dir),
            "resume": True,
        }
    )
    assert r1["status"] == "success" and r1["pass_qc"] is True
    assert r1.get("resumed") is not True
    assert r1["per_file"][0]["balance_info"]["converged_by_resolution"].get(str(bs)) is True
    r2 = hq.hic_qc_and_balance(
        {
            "files": [str(cool_p)],
            "sample_ids": ["s1"],
            "resolution": bs,
            "min_valid_pairs": 10_000_000,  # threshold changed -> MUST NOT resume stale pass_qc=True!
            "output_dir": str(qc_dir),
            "resume": True,
        }
    )
    assert r2["status"] == "success"
    assert r2.get("resumed") is not True
    assert r2["pass_qc"] is False

    # 3. D4-M7: _render_contact_map honors chrom="chr2" instead of hardcoded "chr1".
    figs_dir = tmp_path / "figs"
    figs_dir.mkdir(parents=True, exist_ok=True)
    fig_bundle = hf._render_contact_map(str(cool_p), bs, None, figs_dir, chrom="chr2")
    assert fig_bundle is not None
    svg_txt = Path(fig_bundle["svg"]).read_text(encoding="utf-8")
    assert "chr2" in svg_txt
    assert "(chr1," not in svg_txt


def test_d3_live_chromatin_env_wnn_hvg_fallback_and_cache_digest(tmp_path: Path, monkeypatch):
    ad = pytest.importorskip("anndata")
    md = pytest.importorskip("mudata")
    pytest.importorskip("muon")
    sc = pytest.importorskip("scanpy")
    from aria.scripts.integration_wnn import integration_wnn

    rng = np.random.default_rng(1234)
    n_cells = 30
    bcs = [f"CELL_{i:03d}" for i in range(n_cells)]
    rna_x = rng.poisson(3, size=(n_cells, 25)).astype(np.int64)
    rna_x[:15, :8] += 6
    atac_x = rng.poisson(2, size=(n_cells, 25)).astype(np.int64)
    atac_x[:15, :8] += 5
    mdata = md.MuData(
        {
            "rna": ad.AnnData(
                X=sp.csr_matrix(rna_x),
                obs=pd.DataFrame(index=bcs),
                var=pd.DataFrame(index=[f"G{i}" for i in range(25)]),
            ),
            "atac": ad.AnnData(
                X=sp.csr_matrix(atac_x),
                obs=pd.DataFrame(index=bcs),
                var=pd.DataFrame(index=[f"chr1:{i*500}-{i*500+300}" for i in range(25)]),
            ),
        }
    )
    h5mu_in = tmp_path / "in.h5mu"
    mdata.write_h5mu(h5mu_in)

    # 1. D3-M5: HVG returning < 3 genes triggers variance fallback + warning and succeeds!
    def _deg_hvg(adata, **kw):
        flags = np.zeros(adata.n_vars, dtype=bool)
        flags[0] = True  # only 1 HVG (< 3!)
        adata.var["highly_variable"] = flags

    monkeypatch.setattr(sc.pp, "highly_variable_genes", _deg_hvg)
    out_dir = tmp_path / "wnn_out"
    res = integration_wnn(
        {
            "h5mu_path": str(h5mu_in),
            "output_dir": str(out_dir),
            "min_cells": 20,
            "k": 5,
            "n_top_genes": 10,
            "n_pcs": 5,
            "n_lsi_comps": 5,
            "seed": 0,
            "resume": True,
        }
    )
    assert res["status"] == "success", res
    assert res["n_hvg"] == 10
    assert any("fell back to top 10" in w for w in res.get("warnings", []))

    # 2. D3-H13: stage manifest stores receipt_sha256 & output_sha256; tampering invalidates resume!
    manifest_p = out_dir / "integration_wnn.stage_manifest.json"
    mpay = json.loads(manifest_p.read_text(encoding="utf-8"))
    assert mpay["extra"].get("receipt_sha256")
    assert mpay["extra"].get("output_sha256")

    res_resumed = integration_wnn(
        {
            "h5mu_path": str(h5mu_in),
            "output_dir": str(out_dir),
            "min_cells": 20,
            "k": 5,
            "n_top_genes": 10,
            "n_pcs": 5,
            "n_lsi_comps": 5,
            "seed": 0,
            "resume": True,
        }
    )
    assert res_resumed["resumed"] is True

    # Tamper with cached receipt -> MUST recompute (resumed is False)!
    receipt_p = out_dir / ".integration_wnn.result.json"
    tampered = json.loads(receipt_p.read_text(encoding="utf-8"))
    tampered["n_cells"] = 999999
    receipt_p.write_text(json.dumps(tampered), encoding="utf-8")
    res_recomputed = integration_wnn(
        {
            "h5mu_path": str(h5mu_in),
            "output_dir": str(out_dir),
            "min_cells": 20,
            "k": 5,
            "n_top_genes": 10,
            "n_pcs": 5,
            "n_lsi_comps": 5,
            "seed": 0,
            "resume": True,
        }
    )
    assert res_recomputed["status"] == "success"
    assert res_recomputed["resumed"] is False
    assert res_recomputed["n_cells"] == n_cells
