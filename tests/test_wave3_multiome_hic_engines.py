"""Wave 3C regression tests: Multiome & Hi-C scientific engine fixes."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from aria.scripts.hic_inspect import _inspect_pairs
from aria.scripts.hic_topology import (
    _annotate_loops_with_gtf,
    _get_chromosomes,
    _parse_gtf_tss,
    _snap_window,
)
from aria.scripts.integration_eregulons import _load_precomputed_p2g_csv
from aria.scripts.integration_mofa import _compute_categorical_eta2
from aria.scripts.integration_wnn import _compute_graph_ilisi
from aria.utils.eregulon_core import (
    MemeMotif,
    compute_gc_stratified_motif_enrichment,
    fetch_peak_center_sequences_from_indexed_fasta,
    resolve_motifs_to_expressed_tfs,
)
from aria.utils.hic_constants import estimate_ram_gb
from aria.utils.peak2gene_core import parse_gtf_tss_detailed
from aria.utils.wnn_diagnostics import (
    _evaluate_rare_groups,
    compute_wnn_anti_degradation_diagnostics,
)


def test_eregulon_fasta_lookup_handles_ensembl_vs_ucsc_contigs(tmp_path: Path):
    seq = "ACGTACGTAC" * 100  # 1000 bp
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{seq}\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text(f"chr1\t{len(seq)}\t6\t{len(seq)}\t{len(seq) + 1}\n", encoding="utf-8")

    bundle = fetch_peak_center_sequences_from_indexed_fasta(
        ["1:200-400", "chr1:300-500"], fa, scan_window_bp=100
    )
    assert bundle["n_missing_contig"] == 0
    assert int(bundle["valid_mask"].sum()) == 2


def test_eregulon_cmh_excludes_skipped_strata_from_obs_hits():
    # 10 peaks: 5 low-GC (mixed fg/bg), 5 high-GC (all fg -> skipped stratum with n0_b=0).
    # Motif 0 hits ONLY the skipped high-GC fg peaks -> obs_hits and exp_hits in valid strata are both 0 -> log2_fe == 0.0!
    hits = np.array([[0, 0, 0, 0, 0, 1, 1, 1, 1, 1]], dtype=bool)
    fg = np.array([1, 0, 0, 0, 0, 1, 1, 1, 1, 1], dtype=bool)
    gc = np.array([0.30, 0.31, 0.32, 0.33, 0.34, 0.70, 0.71, 0.72, 0.73, 0.74])
    res = compute_gc_stratified_motif_enrichment(hits, fg, gc, n_gc_bins=2)
    assert res["log2_fold_enrichment"][0] == pytest.approx(0.0)


def test_eregulon_duplicate_identical_gene_symbols_not_ambiguous():
    motif = MemeMotif(
        motif_id="MA0143.1",
        raw_name="SOX2",
        pfm=np.full((4, 4), 0.25),
        bg_probs=(0.25, 0.25, 0.25, 0.25),
    )
    resolved, unmapped = resolve_motifs_to_expressed_tfs(
        [motif],
        ["Sox2", "Sox2", "Gapdh"],
        np.array([50, 40, 100]),
        min_cells_detected=10,
    )
    assert len(resolved) == 1
    assert resolved[0].tf_gene == "Sox2"
    assert not unmapped


def test_peak2gene_gtf_matches_unversioned_ensembl_and_casefolded_symbols(tmp_path: Path):
    gtf = tmp_path / "test.gtf"
    gtf.write_text(
        'chr1\tENSEMBL\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000186092.7"; gene_name "OR4F5"; gene_type "protein_coding";\n'
        'chr3\tENSEMBL\tgene\t2000\t6000\t.\t+\t.\tgene_id "ENSMUSG00000074637.3"; gene_name "Sox2"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    recs, _ = parse_gtf_tss_detailed(str(gtf), wanted={"ENSG00000186092", "SOX2"})
    assert "ENSG00000186092" in recs
    assert "SOX2" in recs


def test_wnn_diagnostics_singleton_nan_weights_and_falsy_zero():
    n = 30
    rng = np.random.default_rng(7)
    rna_emb = rng.normal(size=(n, 5))
    atac_emb = rng.normal(size=(n, 5))
    # Build simple kNN graph with self-loop + finite edges
    rows, cols, data = [], [], []
    for i in range(n):
        rows.append(i); cols.append(i); data.append(1.0)
        for off in (1, 2, 3, 4):
            j = (i + off) % n
            rows.append(i); cols.append(j); data.append(0.5)
    g = sp.csr_matrix((data, (rows, cols)), shape=(n, n))

    # Singleton rare group (size=1) must not raise ZeroDivisionError even if min_rare_cells=1
    labels = ["singleton"] + ["c0"] * 14 + ["c1"] * 15
    rare_res = _evaluate_rare_groups(
        labels,
        labels,
        [np.array([(i + 1) % n, (i + 2) % n]) for i in range(n)],
        np.array([[(i + 1) % n, (i + 2) % n] for i in range(n)]),
        None,
        source_name="test",
        max_rare_count=5,
        min_rare_cells=1,
        min_evaluable_rare_cells=2,
    )
    assert rare_res["n_rare_groups_detected"] == 0

    # Single NaN weight must not hide global_modality_imbalance or leak NaN into JSON
    rna_w = np.full(n, 0.95)
    atac_w = np.full(n, 0.05)
    rna_w[0] = np.nan
    atac_w[0] = np.nan
    diag = compute_wnn_anti_degradation_diagnostics(
        rna_embedding=rna_emb,
        atac_embedding=atac_emb,
        wnn_distances=g,
        graph_is_connectivity=True,
        wnn_clusters=labels,
        rna_clusters=labels,
        atac_clusters=labels,
        rna_weights=rna_w,
        atac_weights=atac_w,
        k=4,
    )
    wa = diag["per_cluster_weight_audit"]
    assert wa["evaluable"] is True
    assert math.isfinite(wa["global_mean_rna_weight"])
    assert "global_modality_imbalance" in diag["issue_codes"]


def test_mofa_anova_and_wnn_ilisi_never_fabricate_1_on_nan():
    z = np.array([1.0, 2.0, np.nan, 4.0, 5.0, 6.0])
    w = np.ones(6)
    labels = ["A", "A", "A", "B", "B", "B"]
    eta2, fstat = _compute_categorical_eta2(z, labels, weights=w)
    assert eta2 is not None and math.isfinite(eta2) and 0.0 <= eta2 <= 1.0
    assert fstat is not None and math.isfinite(fstat)
    # All-NaN returns None (not fabricated 1.0!)
    eta2_nan, _ = _compute_categorical_eta2(np.full(6, np.nan), labels, weights=w)
    assert eta2_nan is None

    # Graph iLISI must ignore diagonal self-loops and +inf edge weights
    n = 6
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        mat[i, i] = 100.0  # self-loop must be ignored
        mat[i, (i + 3) % n] = 1.0  # cross-batch neighbor
    mat[0, 1] = float("inf")
    med, norm = _compute_graph_ilisi(sp.csr_matrix(mat), ["b0", "b0", "b0", "b1", "b1", "b1"])
    assert med is not None and math.isfinite(med)
    assert norm is not None and math.isfinite(norm)


def test_precomputed_p2g_csv_loader_filters_nan_and_float_distances(tmp_path: Path):
    csv_p = tmp_path / "p2g.csv"
    pd.DataFrame(
        [
            {"gene": "G1", "peak": "chr1:100-500", "correlation": 0.4, "padj": 0.01, "distance_bp": "1000.0"},
            {"gene": "G2", "peak": "chr1:600-900", "correlation": float("nan"), "padj": 0.01, "distance_bp": "500"},
            {"gene": "G3", "peak": "chr1:1000-1500", "correlation": 0.01, "padj": 0.01, "distance_bp": "200"},
            {"gene": "G4", "peak": "chr1:2000-2500", "correlation": 0.6, "padj": 0.5, "distance_bp": "300"},
        ]
    ).to_csv(csv_p, index=False)
    res = _load_precomputed_p2g_csv(csv_p, min_corr=0.05, max_padj=0.05)
    assert len(res["links"]) == 1
    assert res["links"][0]["gene"] == "G1"
    assert res["links"][0]["distance_bp"] == 1000


def test_hic_constants_monotonic_and_pairs_header_no_500_cap(tmp_path: Path):
    resolutions = [1_000, 2_000, 5_000, 10_000, 25_000, 40_000, 100_000, 500_000, 1_000_000]
    rams = [estimate_ram_gb(r) for r in resolutions]
    for a, b in zip(rams, rams[1:]):
        assert a > b

    # 600 #chromsize lines must ALL be read without 500-line truncation.
    pairs = tmp_path / "big_header.pairs"
    lines = ["## pairs format v1.0"] + [f"#chromsize: chr{i} 1000000" for i in range(1, 601)] + ["r1\tchr1\t10\tchr1\t20\t+\t+"]
    pairs.write_text("\n".join(lines) + "\n", encoding="utf-8")
    info = _inspect_pairs(str(pairs))
    assert len(info["chromosomes"]) == 600


def test_hic_topology_snap_window_gff3_and_chrom_order(tmp_path: Path):
    assert _snap_window(300_000, 40_000) % 40_000 == 0
    assert _snap_window(500_000, 40_000) % 40_000 == 0
    assert _snap_window(10_000, 40_000) == 120_000  # >= 3 * res

    gff3 = tmp_path / "genes.gff3"
    gff3.write_text(
        "1\tENSEMBL\tgene\t10001\t15000\t.\t+\t.\tID=gene:ENSG1;Name=MYC\n",
        encoding="utf-8",
    )
    tss = _parse_gtf_tss(str(gff3))
    assert len(tss) == 1 and tss[0] == ("1", 10000, "MYC")

    # Loop annotation matches UCSC chr1 loop against Ensembl 1 GFF3 gene!
    loops = pd.DataFrame(
        [{"chrom1": "chr1", "start1": 9000, "end1": 12000, "chrom2": "chr1", "start2": 90000, "end2": 95000}]
    )
    ann, counts = _annotate_loops_with_gtf(loops, str(gff3))
    assert counts["promoter_enhancer"] == 1
    assert ann["anchor1_genes"].iloc[0] == "MYC"

    class _MockClr:
        chromnames = ["chr1", "chr2", "chr3"]

    assert _get_chromosomes(_MockClr(), "chr2") == ["chr2"]
    assert _get_chromosomes(_MockClr(), ["chr3", "chr1"]) == ["chr1", "chr3"]
