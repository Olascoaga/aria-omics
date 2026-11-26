"""Analytical oracle tests for aria/utils/peak2gene_core.py.

Verifies:
  - Contig normalization (UCSC/Ensembl/MT aliases), peak coordinate parsing, GTF TSS
    deduplication priority (PAR_Y drop, MANE_Select > protein_coding > basic > span > gene_id)
  - Indexed FASTA GC calculation (and off-by-one 0-based BED start check)
  - Normalization transforms: _log1p_cp10k_csc & _signac_tfidf_csc dense vs sparse parity
    and exact closed-form values
  - scalar_pearson_corr & vectorized_pearson_corr against scipy.stats.pearsonr
  - compute_peak_to_gene_links:
      * Exact Pearson r, trans-chromosomal background z-score, two-sided erfc p-value
        against scipy.stats.norm.sf, one-sided Signac p-value, and pre-filter BH padj
        against statsmodels.multipletests / scipy.stats.false_discovery_control
      * Dense vs sparse CSC/CSR exact parity
      * Cell-order, gene-order, and peak-order permutation invariance
      * Strand-aware signed tss_offset_bp (+ vs - strand) and distance modes
      * Single-chromosome same_chrom_distal fallback + warning
      * Edge cases: <3 cells, mismatched cell counts, all-zero genes/peaks, constant
        vectors, contig mismatch error
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse
from scipy.stats import false_discovery_control, norm, pearsonr

from aria.utils.peak2gene_core import (
    _log1p_cp10k_csc,
    _quantile_bin_labels,
    _signac_tfidf_csc,
    canonical_contig,
    compute_peak_gc_from_indexed_fasta,
    compute_peak_to_gene_links,
    genes_to_peak_indices,
    is_standard_contig,
    parse_gtf_tss_detailed,
    parse_peak_coord,
    scalar_pearson_corr,
    vectorized_pearson_corr,
)


class TestPeak2GeneHelpersOracle:
    def test_contig_normalization_and_peak_coord_parsing(self):
        assert canonical_contig("1") == "chr1"
        assert canonical_contig("chr1") == "chr1"
        assert canonical_contig("x") == "chrX"
        assert canonical_contig("MT") == "chrM"
        assert canonical_contig("chrMT") == "chrM"
        assert canonical_contig("GL000195.1") == "GL000195.1"
        assert canonical_contig("") is None
        assert canonical_contig(None) is None

        assert is_standard_contig("chr1") is True
        assert is_standard_contig("22") is True
        assert is_standard_contig("X") is True
        assert is_standard_contig("chrM") is False
        assert is_standard_contig("GL000195.1") is False

        assert parse_peak_coord("chr1:100-500") == ("chr1", 100, 500)
        assert parse_peak_coord("chr1-500-100") == ("chr1", 100, 500)
        assert parse_peak_coord("chr2\t200\t400\textra") == ("chr2", 200, 400)
        assert parse_peak_coord("not_a_peak") is None
        assert parse_peak_coord("") is None

    def test_parse_gtf_tss_detailed_priority_and_par_y_drop(self, tmp_path):
        gtf = tmp_path / "ann.gtf"
        gtf.write_text(
            "# header\n"
            # PAR_Y copy must be dropped
            'chrY\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG001_PAR_Y"; gene_name "SHOX"; gene_type "protein_coding";\n'
            'chrX\tHAVANA\tgene\t5000\t8000\t.\t+\t.\tgene_id "ENSG001"; gene_name "SHOX"; gene_type "protein_coding";\n'
            # Duplicate symbol: MANE_Select beats longer non-MANE
            'chr1\tHAVANA\tgene\t100\t90000\t.\t+\t.\tgene_id "ENSG002"; gene_name "DUP1"; gene_type "protein_coding";\n'
            'chr1\tHAVANA\tgene\t2000\t4000\t.\t-\t.\tgene_id "ENSG003"; gene_name "DUP1"; gene_type "protein_coding"; tag "MANE_Select";\n',
            encoding="utf-8",
        )
        recs, stats = parse_gtf_tss_detailed(str(gtf))
        assert stats["n_par_y_dropped"] == 1
        assert stats["n_duplicate_symbols_resolved"] == 1
        assert recs["SHOX"]["canon_chrom"] == "chrX"
        assert recs["SHOX"]["tss"] == 5000  # + strand -> start
        assert recs["DUP1"]["gene_id"] == "ENSG003"
        assert recs["DUP1"]["tss"] == 4000  # - strand -> end

    def test_log1p_cp10k_and_signac_tfidf_exact_dense_sparse_parity(self):
        raw = np.array([
            [10.0, 30.0, 0.0],
            [0.0, 0.0, 0.0],    # zero-count cell must yield all zeros, no NaN
            [5.0, 15.0, 20.0],
        ])
        # Closed-form log1p(1e4 * x / row_sum)
        exp_cp10k = np.array([
            [np.log1p(1e4 * 10 / 40), np.log1p(1e4 * 30 / 40), 0.0],
            [0.0, 0.0, 0.0],
            [np.log1p(1e4 * 5 / 40), np.log1p(1e4 * 15 / 40), np.log1p(1e4 * 20 / 40)],
        ])
        d_cp = _log1p_cp10k_csc(raw)
        s_cp = _log1p_cp10k_csc(sparse.csr_matrix(raw)).toarray()
        np.testing.assert_allclose(d_cp, exp_cp10k, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(s_cp, exp_cp10k, rtol=1e-12, atol=1e-12)

        # Closed-form Signac TF-IDF: log1p(1e4 * (x / cell_sum) * (n_cells / peak_sum))
        cell_sums = np.array([[40.0], [1.0], [40.0]])
        peak_sums = np.array([[15.0, 45.0, 20.0]])
        tf = np.where(raw.sum(axis=1, keepdims=True) > 0, raw / cell_sums, 0.0)
        idf = 3.0 / peak_sums
        exp_tfidf = np.log1p(tf * idf * 1e4)
        d_tf = _signac_tfidf_csc(raw)
        s_tf = _signac_tfidf_csc(sparse.csr_matrix(raw)).toarray()
        np.testing.assert_allclose(d_tf, exp_tfidf, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(s_tf, exp_tfidf, rtol=1e-12, atol=1e-12)

    def test_scalar_and_vectorized_pearson_match_scipy_pearsonr(self):
        rng = np.random.default_rng(14)
        x = rng.normal(size=40)
        cols = np.column_stack([
            2.0 * x + 1.0,                 # r = +1
            -3.0 * x + 5.0,                # r = -1
            rng.normal(size=40),           # arbitrary
            np.full(40, 7.0),              # constant -> NaN
        ])
        r_vec = vectorized_pearson_corr(x, cols)
        assert r_vec[0] == pytest.approx(1.0, abs=1e-12)
        assert r_vec[1] == pytest.approx(-1.0, abs=1e-12)
        assert r_vec[2] == pytest.approx(pearsonr(x, cols[:, 2]).statistic, abs=1e-12)
        assert np.isnan(r_vec[3])

        assert scalar_pearson_corr(x, cols[:, 2]) == pytest.approx(
            round(float(pearsonr(x, cols[:, 2]).statistic), 6), abs=1e-12
        )
        assert scalar_pearson_corr(x, cols[:, 3]) is None
        assert scalar_pearson_corr([1.0, 2.0], [1.0, 2.0]) is None
        assert np.all(np.isnan(vectorized_pearson_corr([1.0, 2.0], np.ones((2, 2)))))

    def test_quantile_bin_labels_and_genes_to_peak_indices(self):
        bins = _quantile_bin_labels(np.array([10.0, 30.0, 20.0, 40.0]), n_bins=2)
        assert bins.tolist() == [0, 1, 0, 1]
        assert _quantile_bin_labels(np.array([1.0, 2.0]), n_bins=1).tolist() == [0, 0]

        mapping = genes_to_peak_indices(
            ["chr1:1000-1200", "1:50000-50200", "chr2:1000-1200"],
            {"G1": ("1", 1100)},
            window_bp=5000,
        )
        assert mapping == {"G1": [0]}

    # D2-BUG-7 fixed by D4-017 on fix/d1-silent-failures; xfail marker removed at merge.
    def test_indexed_fasta_gc_exact_on_0based_bed_coordinates(self, tmp_path):
        # 0-based sequence: indices 0..9 = 'A'*10, 10..19 = 'G'*10, 20..29 = 'A'*10
        fa = tmp_path / "toy.fa"
        seq = "AAAAAAAAAA" + "GGGGGGGGGG" + "AAAAAAAAAA"
        fa.write_text(f">chr1\n{seq}\n", encoding="utf-8")
        fai = tmp_path / "toy.fa.fai"
        # Header ">chr1\n" is 6 bytes -> offset=6, line_bases=30, line_bytes=31
        fai.write_text("chr1\t30\t6\t30\t31\n", encoding="utf-8")

        gc, info = compute_peak_gc_from_indexed_fasta(
            str(fa), [("chr1", 0, 10), ("chr1", 10, 20), ("chr1", 20, 30)]
        )
        assert info["available"] is True
        # Peak [10, 20) is strictly 10 G's -> exact GC fraction = 1.0 (currently returns 10/11)
        np.testing.assert_allclose(gc, [0.0, 1.0, 0.0], atol=1e-12)


class TestComputePeakToGeneLinksOracle:
    def _make_multi_chrom_fixture(self, n_cells: int = 80, seed: int = 501):
        rng = np.random.default_rng(seed)
        # Latent cell state driver
        latent = rng.normal(0.0, 1.0, size=n_cells)

        # 2 genes on chr1: G_pos (+ strand, TSS=10000), G_neg (- strand, TSS=50000)
        g_pos = np.maximum(0.0, 5.0 + 3.0 * latent + rng.normal(0.0, 0.3, size=n_cells))
        g_neg = np.maximum(0.0, 5.0 - 2.5 * latent + rng.normal(0.0, 0.3, size=n_cells))
        g_zero = np.zeros(n_cells)
        rna = np.column_stack([g_pos, g_neg, g_zero])
        rna_names = ["G_pos", "G_neg", "G_zero"]

        gene_records = {
            "G_pos": {"gene": "G_pos", "gene_id": "E1", "chrom": "chr1", "canon_chrom": "chr1", "tss": 10_000, "strand": "+"},
            "G_neg": {"gene": "G_neg", "gene_id": "E2", "chrom": "chr1", "canon_chrom": "chr1", "tss": 50_000, "strand": "-"},
            "G_zero": {"gene": "G_zero", "gene_id": "E3", "chrom": "chr1", "canon_chrom": "chr1", "tss": 80_000, "strand": "+"},
            # Annotated gene on chr2 (absent from rna_var_names) so GTF contig universe includes chr2
            "G_chr2_ref": {"gene": "G_chr2_ref", "gene_id": "E4", "chrom": "chr2", "canon_chrom": "chr2", "tss": 500_000, "strand": "+"},
        }

        # 4 cis peaks on chr1 + 40 independent trans background peaks on chr2
        p_cis_pos = np.maximum(0.0, 4.0 + 2.8 * latent + rng.normal(0.0, 0.3, size=n_cells))   # strong + with G_pos, - with G_neg
        p_cis_neg = np.maximum(0.0, 4.0 - 2.5 * latent + rng.normal(0.0, 0.3, size=n_cells))   # strong - with G_pos, + with G_neg
        p_cis_prom = np.maximum(0.0, 4.0 + 2.6 * latent + rng.normal(0.0, 0.3, size=n_cells))  # overlapping G_pos TSS (dist=0)
        p_cis_null = rng.uniform(1.0, 5.0, size=n_cells)
        bg_trans = rng.uniform(1.0, 5.0, size=(n_cells, 40))

        atac = np.column_stack([p_cis_pos, p_cis_neg, p_cis_prom, p_cis_null, bg_trans])
        atac_names = [
            "chr1:15000-15400",    # +5000 bp downstream of G_pos (+ strand), +34800 bp upstream of G_neg (- strand)
            "chr1:45000-45400",    # +35000 bp from G_pos, +4800 bp downstream of G_neg (- strand: -(45200-50000)=+4800)
            "chr1:9900-10200",     # overlaps G_pos TSS=10000 -> distance_bp=0 (promoter_proximal)
            "chr1:20000-20400",    # null cis peak
        ] + [f"chr2:{i*2000}-{i*2000 + 250 + i * 11}" for i in range(40)]

        return rna, rna_names, atac, atac_names, gene_records

    def test_exact_pearson_zscore_pvalue_bh_and_dense_sparse_parity(self):
        sm_multitest = pytest.importorskip("statsmodels.stats.multitest")
        rna, rna_names, atac, atac_names, gene_records = self._make_multi_chrom_fixture(seed=501)

        res_dense = compute_peak_to_gene_links(
            rna_matrix=rna,
            rna_var_names=rna_names,
            atac_matrix=atac,
            atac_var_names=atac_names,
            gene_records=gene_records,
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="as_is",
            background_null=True,
            n_bg_bins_per_axis=1,   # single stratum -> all 40 trans peaks form exact known bg pool!
            n_bg_per_bin=200,
            min_bg_peaks=10,
            promoter_bp=2000,
            seed=0,
        )
        res_sparse = compute_peak_to_gene_links(
            rna_matrix=sparse.csr_matrix(rna),
            rna_var_names=rna_names,
            atac_matrix=sparse.csc_matrix(atac),
            atac_var_names=atac_names,
            gene_records=gene_records,
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="as_is",
            background_null=True,
            n_bg_bins_per_axis=1,
            n_bg_per_bin=200,
            min_bg_peaks=10,
            promoter_bp=2000,
            seed=0,
        )
        assert res_dense["status"] == "success"
        assert res_dense["n_genes_tested"] == 2  # G_zero excluded!
        assert res_dense["n_tested"] == 8        # 2 genes x 4 cis peaks on chr1
        assert res_dense["links"] == res_sparse["links"]

        # Independently compute EVERY statistic for all 8 pairs from first principles!
        links_by_key = {(r["gene"], r["peak"]): r for r in res_dense["links"]}
        p_two_ref = []
        keys_order = []
        for gi, gname in enumerate(["G_pos", "G_neg"]):
            bg_r = np.array([pearsonr(rna[:, gi], atac[:, 4 + k]).statistic for k in range(40)])
            mu_bg = float(np.mean(bg_r))
            sd_bg = float(np.std(bg_r, ddof=1))
            for pj in range(4):
                pname = atac_names[pj]
                rec = links_by_key[(gname, pname)]
                r_exact = float(pearsonr(rna[:, gi], atac[:, pj]).statistic)
                z_exact = (r_exact - mu_bg) / sd_bg
                p_two_exact = float(2.0 * norm.sf(abs(z_exact)))
                p_one_exact = 0.5 * p_two_exact

                assert rec["correlation"] == pytest.approx(round(r_exact, 6), abs=1e-6)
                assert rec["zscore"] == pytest.approx(round(z_exact, 6), abs=1e-5)
                assert rec["p_value"] == pytest.approx(p_two_exact, rel=1e-6, abs=1e-12)
                assert rec["p_value_signac"] == pytest.approx(p_one_exact, rel=1e-6, abs=1e-12)
                assert rec["background_scope"] == "trans"
                p_two_ref.append(p_two_exact)
                keys_order.append((gname, pname))

        # Verify BH padj matches BOTH statsmodels and scipy across all 8 tested pairs
        _, sm_padj, _, _ = sm_multitest.multipletests(p_two_ref, method="fdr_bh")
        sp_padj = false_discovery_control(p_two_ref, method="bh")
        got_padj = [links_by_key[k]["padj"] for k in keys_order]
        np.testing.assert_allclose(got_padj, sm_padj, rtol=1e-6, atol=1e-12)
        np.testing.assert_allclose(got_padj, sp_padj, rtol=1e-6, atol=1e-12)

        # Verify exact distance_bp, promoter_proximal vs distal_cis, and strand-aware tss_offset_bp
        prom_rec = links_by_key[("G_pos", "chr1:9900-10200")]
        assert prom_rec["distance_bp"] == 0
        assert prom_rec["region_class"] == "promoter_proximal"
        assert prom_rec["direction"] == "positive"

        neg_strand_rec = links_by_key[("G_neg", "chr1:45000-45400")]
        # Midpoint = 45200, TSS = 50000, strand = '-' -> offset = (45200 - 50000) * (-1) = +4800
        assert neg_strand_rec["tss_offset_bp"] == 4800
        assert neg_strand_rec["distance_bp"] == 4600
        assert neg_strand_rec["region_class"] == "distal_cis"
        assert neg_strand_rec["direction"] == "positive"

        # Anti-selection-bias invariant: filtering by min_corr=0.7 MUST NOT change padj of survivors!
        res_filt = compute_peak_to_gene_links(
            rna_matrix=rna,
            rna_var_names=rna_names,
            atac_matrix=atac,
            atac_var_names=atac_names,
            gene_records=gene_records,
            distance_bp=50_000,
            min_corr=0.7,
            preprocess="as_is",
            background_null=True,
            n_bg_bins_per_axis=1,
            n_bg_per_bin=200,
            min_bg_peaks=10,
            seed=0,
        )
        assert res_filt["n_tested"] == 8
        assert res_filt["n_links"] < 8
        for rec in res_filt["links"]:
            assert rec["padj"] == pytest.approx(
                links_by_key[(rec["gene"], rec["peak"])]["padj"], rel=1e-12, abs=1e-12
            )

    def test_cell_gene_and_peak_permutation_invariance(self):
        rna, rna_names, atac, atac_names, gene_records = self._make_multi_chrom_fixture(seed=502)
        base = compute_peak_to_gene_links(
            rna_matrix=rna,
            rna_var_names=rna_names,
            atac_matrix=atac,
            atac_var_names=atac_names,
            gene_records=gene_records,
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="paired",
            rna_is_raw_counts=True,
            atac_is_raw_counts=True,
            atac_transform="tfidf",
            background_null=True,
            n_bg_bins_per_axis=2,
            n_bg_per_bin=200,
            min_bg_peaks=10,
            seed=42,
        )
        rng = np.random.default_rng(999)
        c_perm = rng.permutation(rna.shape[0])
        g_perm = rng.permutation(len(rna_names))
        p_perm = rng.permutation(len(atac_names))

        perm_res = compute_peak_to_gene_links(
            rna_matrix=rna[np.ix_(c_perm, g_perm)],
            rna_var_names=[rna_names[i] for i in g_perm],
            atac_matrix=atac[np.ix_(c_perm, p_perm)],
            atac_var_names=[atac_names[i] for i in p_perm],
            gene_records=gene_records,
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="paired",
            rna_is_raw_counts=True,
            atac_is_raw_counts=True,
            atac_transform="tfidf",
            background_null=True,
            n_bg_bins_per_axis=2,
            n_bg_per_bin=200,
            min_bg_peaks=10,
            seed=42,
        )
        base_map = {(r["gene"], r["peak"]): r for r in base["links"]}
        perm_map = {(r["gene"], r["peak"]): r for r in perm_res["links"]}
        assert set(perm_map) == set(base_map)
        for k, b in base_map.items():
            p = perm_map[k]
            assert p["correlation"] == pytest.approx(b["correlation"], abs=1e-6)
            assert p["zscore"] == pytest.approx(b["zscore"], abs=1e-5)
            assert p["padj"] == pytest.approx(b["padj"], rel=1e-6, abs=1e-12)

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "D2-BUG-10: peak2gene_core._quantile_bin_labels breaks ties via mergesort input "
            "index, so fixed-width (e.g. 500bp) peaks are split across length strata by "
            "arbitrary column order, making background z-scores and padj depend on peak order"
        ),
    )
    def test_fixed_width_peaks_zscores_must_be_invariant_to_peak_order(self):
        rna, rna_names, atac, atac_names, gene_records = self._make_multi_chrom_fixture(seed=502)
        # All 44 peaks are standard fixed-width 500bp intervals (e.g. ArchR/Signac consensus peaks)
        fixed_names = [
            "chr1:15000-15500",
            "chr1:45000-45500",
            "chr1:9900-10400",
            "chr1:20000-20500",
        ] + [f"chr2:{i*2000}-{i*2000+500}" for i in range(40)]

        base = compute_peak_to_gene_links(
            rna_matrix=rna,
            rna_var_names=rna_names,
            atac_matrix=atac,
            atac_var_names=fixed_names,
            gene_records=gene_records,
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="as_is",
            background_null=True,
            n_bg_bins_per_axis=2,
            n_bg_per_bin=200,
            min_bg_peaks=10,
            seed=42,
        )
        p_perm = np.random.default_rng(123).permutation(len(fixed_names))
        perm = compute_peak_to_gene_links(
            rna_matrix=rna,
            rna_var_names=rna_names,
            atac_matrix=atac[:, p_perm],
            atac_var_names=[fixed_names[i] for i in p_perm],
            gene_records=gene_records,
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="as_is",
            background_null=True,
            n_bg_bins_per_axis=2,
            n_bg_per_bin=200,
            min_bg_peaks=10,
            seed=42,
        )
        base_map = {(r["gene"], r["peak"]): r for r in base["links"]}
        perm_map = {(r["gene"], r["peak"]): r for r in perm["links"]}
        for k, b in base_map.items():
            assert perm_map[k]["zscore"] == pytest.approx(b["zscore"], abs=1e-5)

    def test_edge_cases_insufficient_cells_contig_mismatch_and_distal_fallback(self):
        # < 3 cells -> InsufficientSharedCells
        res_2c = compute_peak_to_gene_links(
            rna_matrix=np.ones((2, 1)),
            rna_var_names=["G1"],
            atac_matrix=np.ones((2, 1)),
            atac_var_names=["chr1:100-200"],
            gene_records={"G1": {"gene": "G1", "chrom": "chr1", "canon_chrom": "chr1", "tss": 150}},
        )
        assert res_2c["status"] == "error"
        assert res_2c["error_type"] == "InsufficientSharedCells"

        # Mismatched cell counts -> ValueError
        with pytest.raises(ValueError, match="cell counts differ"):
            compute_peak_to_gene_links(
                rna_matrix=np.ones((5, 1)),
                rna_var_names=["G1"],
                atac_matrix=np.ones((4, 1)),
                atac_var_names=["chr1:100-200"],
                gene_records={},
            )

        # Contig mismatch (>= 20 valid peaks, < 50% overlap with GTF contigs) -> ContigMismatch
        rng = np.random.default_rng(1)
        res_cm = compute_peak_to_gene_links(
            rna_matrix=rng.uniform(1, 5, size=(10, 1)),
            rna_var_names=["G1"],
            atac_matrix=rng.uniform(1, 5, size=(10, 25)),
            atac_var_names=[f"chr9:{i*1000}-{i*1000+200}" for i in range(25)],
            gene_records={"G1": {"gene": "G1", "chrom": "chr1", "canon_chrom": "chr1", "tss": 500}},
        )
        assert res_cm["status"] == "error"
        assert res_cm["error_type"] == "ContigMismatch"

        # Single-chromosome distal fallback when no trans peaks exist
        rna_1c = rng.uniform(1, 5, size=(20, 1))
        atac_1c = rng.uniform(1, 5, size=(20, 10))
        atac_1c[:, 0] = rna_1c[:, 0] * 2.0 + 0.1
        peaks_1c = ["chr1:1000-1200"] + [f"chr1:{1_000_000 + i*5000}-{1_000_200 + i*5000}" for i in range(9)]
        res_dist = compute_peak_to_gene_links(
            rna_matrix=rna_1c,
            rna_var_names=["G1"],
            atac_matrix=atac_1c,
            atac_var_names=peaks_1c,
            gene_records={"G1": {"gene": "G1", "chrom": "chr1", "canon_chrom": "chr1", "tss": 1100}},
            distance_bp=50_000,
            min_corr=0.0,
            background_null=True,
            seed=0,
        )
        assert res_dist["status"] == "success"
        assert res_dist["links"][0]["background_scope"] == "same_chrom_distal"
        assert any("same_chrom_distal" in w for w in res_dist["warnings"])

    def test_peak2gene_helper_and_parameter_mutant_killers(self, tmp_path):
        import math
        import scipy.sparse as sp
        from aria.utils.peak2gene_core import (
            _log1p_cp10k_csc,
            _quantile_bin_labels,
            _signac_tfidf_csc,
            canonical_contig,
            is_standard_contig,
            parse_peak_coord,
            scalar_pearson_corr,
            vectorized_pearson_corr,
        )

        # 1. Contig canonicalization, standard check, and peak parsing edges
        assert canonical_contig("") is None
        assert canonical_contig("MT") == "chrM"
        assert canonical_contig("M") == "chrM"
        assert canonical_contig("chrMT") == "chrM"
        assert canonical_contig("1") == "chr1"
        assert canonical_contig("X") == "chrX"
        assert canonical_contig("chr22") == "chr22"
        assert canonical_contig("GL000195.1") == "GL000195.1"
        assert is_standard_contig("chr1") is True
        assert is_standard_contig("chrM") is False
        assert is_standard_contig("GL000195.1") is False
        assert parse_peak_coord("chr1:100-200") == ("chr1", 100, 200)
        assert parse_peak_coord("chr1-100-200") == ("chr1", 100, 200)
        assert parse_peak_coord("chr1:200-100") == ("chr1", 100, 200)
        assert parse_peak_coord("bad_peak") is None

        # 2. scalar & vectorized pearson edge cases (n<3, mismatched length, constant vector, 1D cols)
        assert scalar_pearson_corr([1.0, 2.0], [1.0, 2.0]) is None
        assert scalar_pearson_corr([1.0, 2.0, 3.0], [1.0, 2.0]) is None
        assert scalar_pearson_corr([2.0, 2.0, 2.0], [1.0, 2.0, 3.0]) is None
        assert scalar_pearson_corr([1.0, 2.0, 3.0], [2.0, 2.0, 2.0]) is None
        assert scalar_pearson_corr([1.0, 2.0, 3.0], [3.0, 2.0, 1.0]) == -1.0
        assert np.all(np.isnan(vectorized_pearson_corr([1.0, 2.0], np.ones((2, 3)))))
        r_1d = vectorized_pearson_corr([1.0, 2.0, 3.0, 4.0], np.array([2.0, 4.0, 6.0, 8.0]))
        np.testing.assert_allclose(r_1d, [1.0], atol=1e-12)
        r_const = vectorized_pearson_corr([1.0, 2.0, 3.0], np.column_stack([[5, 5, 5], [1, 2, 3]]))
        assert np.isnan(r_const[0]) and r_const[1] == pytest.approx(1.0, abs=1e-12)

        # 3. _quantile_bin_labels empty & n_bins<=1
        assert _quantile_bin_labels(np.array([]), 5).size == 0
        np.testing.assert_array_equal(_quantile_bin_labels(np.array([1.0, 2.0, 3.0]), 1), [0, 0, 0])

        # 4. Dense & sparse parity (including all-zero row/col!) for _log1p_cp10k_csc & _signac_tfidf_csc
        mat = np.array([[0.0, 0.0, 0.0], [2.0, 8.0, 0.0], [5.0, 5.0, 0.0]], dtype=float)
        d_cp = _log1p_cp10k_csc(mat, target_sum=1e4)
        s_cp = _log1p_cp10k_csc(sp.csr_matrix(mat), target_sum=1e4).toarray()
        np.testing.assert_allclose(d_cp, s_cp, rtol=1e-12, atol=1e-12)
        assert d_cp[0].sum() == 0.0
        assert d_cp[1, 0] == pytest.approx(math.log1p(2000.0), rel=1e-12)

        d_tf = _signac_tfidf_csc(mat, scale_factor=1e4)
        s_tf = _signac_tfidf_csc(sp.csr_matrix(mat), scale_factor=1e4).toarray()
        np.testing.assert_allclose(d_tf, s_tf, rtol=1e-12, atol=1e-12)
        assert d_tf[0].sum() == 0.0 and d_tf[:, 2].sum() == 0.0
        # Row 1 col 0: tf=2/10=0.2, idf=3/7, scaled=0.2*(3/7)*1e4
        assert d_tf[1, 0] == pytest.approx(math.log1p(0.2 * (3.0 / 7.0) * 1e4), rel=1e-12)

        # 5. compute_peak_to_gene_links options: standard_chroms_only, unparseable peaks,
        #    min_cells_detected on small (<50) and large (>=50) cohorts, min_distance_bp,
        #    distance_mode="endpoint_min", max_padj, max_pvalue_signac, gc_and_accessibility FASTA!
        fa = tmp_path / "mini.fa"
        # chr1: 30kb alternating GC regions; chr2: 40kb
        seq1 = ("GCGC" * 2500) + ("ATAT" * 2500) + ("ACGT" * 2500)
        seq2 = ("GCAT" * 10000)
        with open(fa, "wb") as fh:
            off1 = len(b">chr1\n")
            fh.write(b">chr1\n" + seq1.encode() + b"\n")
            off2 = fh.tell() + len(b">chr2\n")
            fh.write(b">chr2\n" + seq2.encode() + b"\n")
        (tmp_path / "mini.fa.fai").write_text(
            f"chr1\t{len(seq1)}\t{off1}\t{len(seq1)}\t{len(seq1)+1}\n"
            f"chr2\t{len(seq2)}\t{off2}\t{len(seq2)}\t{len(seq2)+1}\n",
            encoding="utf-8",
        )
        rng = np.random.default_rng(707)
        n_cells = 60
        z = rng.normal(0.0, 1.0, size=n_cells)
        rna_m = np.column_stack([
            np.maximum(0, np.round(20 + 10 * z + rng.normal(0, 1, n_cells))),
            np.zeros(n_cells),  # unexpressed gene skipped
        ])
        atac_cols = [
            np.maximum(0, np.round(15 + 8 * z + rng.normal(0, 1, n_cells))),   # overlaps TSS!
            np.maximum(0, np.round(15 - 8 * z + rng.normal(0, 1, n_cells))),   # negative link at 5kb!
            np.ones(n_cells),                                                  # unparseable peak
            np.ones(n_cells),                                                  # nonstandard contig
        ] + [rng.poisson(5, size=n_cells).astype(float) for _ in range(20)]
        atac_m = np.column_stack(atac_cols)
        atac_names = [
            "chr1:9800-10200",
            "chr1:15000-15400",
            "BAD_PEAK_NAME",
            "chrUn_gl000220:100-500",
        ] + [f"chr2:{1000 + i*1000}-{1400 + i*1000}" for i in range(20)]
        g_recs = {
            "G1": {"gene": "G1", "chrom": "chr1", "canon_chrom": "chr1", "tss": 10000, "strand": "-"},
            "G_ZERO": {"gene": "G_ZERO", "chrom": "chr1", "canon_chrom": "chr1", "tss": 10000, "strand": "+"},
            "G_CHR2": {"gene": "G_CHR2", "chrom": "chr2", "canon_chrom": "chr2", "tss": 1000, "strand": "+"},
        }
        res_gc = compute_peak_to_gene_links(
            rna_matrix=sp.csr_matrix(rna_m),
            rna_var_names=["G1", "G_ZERO"],
            atac_matrix=sp.csr_matrix(atac_m),
            atac_var_names=atac_names,
            gene_records=g_recs,
            distance_bp=20_000,
            min_distance_bp=0,
            min_corr=0.2,
            max_padj=0.05,
            max_pvalue_signac=0.05,
            preprocess="paired",
            rna_is_raw_counts=False,
            atac_is_raw_counts=False,
            min_cells_detected=5,
            standard_chroms_only=True,
            background_null=True,
            genome_fasta=str(fa),
            n_bg_bins_per_axis=2,
            n_bg_per_bin=5,
            min_bg_peaks=5,
            distance_mode="endpoint_min",
            seed=0,
        )
        assert res_gc["status"] == "success", res_gc
        assert res_gc["background_match"] == "gc_and_accessibility"
        assert res_gc["background_method"] == "binned_matched_background"
        assert res_gc["n_unparseable_peaks"] == 1
        assert res_gc["n_nonstandard_peaks_dropped"] == 1
        assert res_gc["min_cells_detected_effective"] == 5
        assert res_gc["rna_normalization"] == "pre_normalized_input"
        assert res_gc["atac_normalization"] == "pre_normalized_input"
        assert res_gc["n_positive"] == 1 and res_gc["n_negative"] == 1
        assert res_gc["n_promoter_proximal"] == 1 and res_gc["n_distal_cis"] == 1
        by_pk = {lk["peak"]: lk for lk in res_gc["links"]}
        # endpoint_min on chr1:9800-10200 vs tss=10000 gives dist=200 (not 0!)
        assert by_pk["chr1:9800-10200"]["distance_bp"] == 200
        assert by_pk["chr1:9800-10200"]["correlation"] > 0.7
        assert by_pk["chr1:15000-15400"]["correlation"] < -0.7
        assert by_pk["chr1:15000-15400"]["direction"] == "negative"
        # Minus strand: midpoint 15200 > tss 10000 -> tss_offset_bp == -5200
        assert by_pk["chr1:15000-15400"]["tss_offset_bp"] == -5200
        assert by_pk["chr1:15000-15400"]["distance_kb"] == pytest.approx(5.0, abs=1e-3)
        assert by_pk["chr1:15000-15400"]["p_value_signac"] == pytest.approx(
            0.5 * by_pk["chr1:15000-15400"]["p_value"], rel=1e-6, abs=1e-12
        )

    def test_peak2gene_dense_sorting_thresholds_and_contig_mutant_killers(self):
        import math
        import scipy.sparse as sp
        from aria.utils.peak2gene_core import (
            _log1p_cp10k_csc,
            _quantile_bin_labels,
            _signac_tfidf_csc,
            canonical_contig,
            is_standard_contig,
            parse_peak_coord,
            vectorized_pearson_corr,
        )

        # 1. Contig canonicalization & standard checks (chr0, chr99, chr100, chrY, chrM, non-chr)
        assert canonical_contig("chrM") == "chrM"
        assert canonical_contig("chrY") == "chrY"
        assert canonical_contig("Y") == "chrY"
        assert canonical_contig("chr1") == "chr1"
        assert canonical_contig("99") == "chr99"
        assert canonical_contig("0") == "0"
        assert canonical_contig("100") == "100"
        assert is_standard_contig("chrY") is True
        assert is_standard_contig("chrX") is True
        assert is_standard_contig("") is False
        assert is_standard_contig("unplaced_scaffold") is False

        # 2. parse_peak_coord tab-3, tab-4, underscore, extra colon
        assert parse_peak_coord("chr1\t10\t20") == ("chr1", 10, 20)
        assert parse_peak_coord("chr1\t10\t20\textra_name") == ("chr1", 10, 20)
        assert parse_peak_coord("chr1:10_20") == ("chr1", 10, 20)
        assert parse_peak_coord("chr1:10:20") is None

        # 3. _quantile_bin_labels non-trivial exact bins + n_bins==0
        np.testing.assert_array_equal(_quantile_bin_labels(np.array([1.0, 2.0]), 0), [0, 0])
        np.testing.assert_array_equal(
            _quantile_bin_labels(np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0]), 3),
            [0, 0, 1, 1, 2, 2],
        )

        # 4. Non-square (4 cells x 6 peaks) _signac_tfidf_csc & _log1p_cp10k_csc (kills shape[0] vs shape[1]!)
        rect = np.array([
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [2.0, 8.0, 0.0, 1.0, 0.0, 3.0],
            [0.0, 4.0, 6.0, 0.0, 0.0, 1.0],
            [5.0, 0.0, 5.0, 0.0, 0.0, 2.0],
        ], dtype=float)
        tf_d = _signac_tfidf_csc(rect, scale_factor=1e4)
        tf_s = _signac_tfidf_csc(sp.csr_matrix(rect), scale_factor=1e4).toarray()
        np.testing.assert_allclose(tf_d, tf_s, rtol=1e-12, atol=1e-12)
        # Cell 1 total=14, peak 0 total=7, n_cells=4 -> tf=2/14, idf=4/7
        assert tf_d[1, 0] == pytest.approx(math.log1p((2.0 / 14.0) * (4.0 / 7.0) * 1e4), rel=1e-12)
        assert tf_d[:, 4].sum() == 0.0 and tf_d[0, :].sum() == 0.0

        # 5. vectorized_pearson_corr constant x -> all NaN of exact length k
        r_cx = vectorized_pearson_corr(np.ones(5), np.arange(15, dtype=float).reshape(5, 3))
        assert r_cx.shape == (3,) and np.all(np.isnan(r_cx))
        r_n2_1d = vectorized_pearson_corr([1.0, 2.0], np.array([1.0, 2.0]))
        assert r_n2_1d.shape == (1,) and np.isnan(r_n2_1d[0])

        # 6. Dense ndarray inputs, exact n_cells==49 vs 50 eff_min_cells boundary + > vs >= peak_nnz,
        #    raw_counts label, gene_id preservation, zero-count peak drop, background_null=False,
        #    and descending |zscore| sort!
        rng = np.random.default_rng(808)
        for n_c, exp_eff, exp_tested in [(49, 9, 3), (50, 10, 2)]:
            z = rng.normal(0.0, 1.0, size=n_c)
            rna_d = np.column_stack([
                np.maximum(0.0, np.round(20 + 10 * z + rng.normal(0, 0.5, n_c))),
            ])
            # Peak with exact nnz == exp_eff: kept when n_cells < 50 (>=), dropped when n_cells >= 50 (>)!
            edge_pk = np.zeros(n_c)
            edge_pk[:exp_eff] = 5.0 + np.arange(exp_eff)
            atac_d = np.column_stack([
                np.maximum(0.0, np.round(20 + 10 * z + rng.normal(0, 0.5, n_c))),  # strong +
                np.maximum(0.0, np.round(20 + 5 * z + rng.normal(0, 2.0, n_c))),   # moderate +
                edge_pk,                                                           # exact boundary peak
                np.zeros(n_c),                                                     # all-zero peak -> dropped
            ] + [rng.poisson(5, size=n_c).astype(float) for _ in range(12)])
            pk_names = ["chr1:12000-12400", "chr1:18000-18400", "chr1:16000-16400", "chr1:14000-14400"] + [
                f"chr2:{1000 + i*2000}-{1400 + i*2000}" for i in range(12)
            ]
            g_recs = {
                "G_PLUS": {
                    "gene": "G_PLUS",
                    "gene_id": "ENSG000001",
                    "chrom": "chr1",
                    "canon_chrom": "chr1",
                    "tss": 10000,
                    "strand": "+",
                }
            }
            res = compute_peak_to_gene_links(
                rna_matrix=rna_d,
                rna_var_names=["G_PLUS"],
                atac_matrix=atac_d,
                atac_var_names=pk_names,
                gene_records=g_recs,
                distance_bp=50_000,
                min_corr=0.0,
                preprocess="paired",
                rna_is_raw_counts=True,
                atac_is_raw_counts=True,
                atac_transform="tfidf",
                min_cells_detected=10,
                background_null=True,
                n_bg_bins_per_axis=2,
                n_bg_per_bin=10,
                min_bg_peaks=5,
                seed=0,
            )
            assert res["status"] == "success"
            assert res["min_cells_detected_effective"] == exp_eff
            assert res["rna_normalization"] == "log1p_cp10k"
            assert res["atac_normalization"] == "signac_tfidf"
            assert res["background_match"] == "accessibility_and_length"
            assert res["gc_matching"] == {"available": False, "reason": "no_genome_fasta"}
            assert res["n_tested"] == exp_tested
            assert res["links"][0]["gene_id"] == "ENSG000001"
            assert res["links"][0]["tss_offset_bp"] > 0
            # Sorted by padj asc then |zscore| desc
            z_abs = [abs(lk["zscore"]) for lk in res["links"]]
            p_vals = [lk["padj"] for lk in res["links"]]
            if p_vals[0] == p_vals[1]:
                assert z_abs[0] >= z_abs[1]

        # preprocess="paired" with min_cells_detected=1 drops all-zero peak via peak_raw_sums <= 0.0
        res_m1 = compute_peak_to_gene_links(
            rna_matrix=rna_d,
            rna_var_names=["G_PLUS"],
            atac_matrix=atac_d[:, [0, 3]],
            atac_var_names=["chr1:12000-12400", "chr1:14000-14400"],
            gene_records=g_recs,
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="paired",
            rna_is_raw_counts=True,
            atac_is_raw_counts=True,
            min_cells_detected=1,
            background_null=False,
        )
        assert res_m1["status"] == "success"
        assert res_m1["n_tested"] == 1
        assert res_m1["atac_normalization"] == "raw_counts"
        assert res_m1["background_method"] == "none"
        assert res_m1["gc_matching"] == {"available": False, "reason": "not_requested"}
        assert all(lk["zscore"] is None for lk in res_m1["links"])

    def test_peak2gene_loop_continues_as_is_sorting_and_bg_boundaries(self):
        import scipy.sparse as sp
        from aria.utils.peak2gene_core import (
            canonical_contig,
            vectorized_pearson_corr,
        )

        # 1. vectorized_pearson_corr n<3 with 2D cols (shape (2,4)) returns shape (4,) of NaNs
        r_short = vectorized_pearson_corr([1.0, 2.0], np.ones((2, 4)))
        assert r_short.shape == (4,) and np.all(np.isnan(r_short))
        # Mixed-case & zero-padded contigs require stripping "chr" prefix to canonicalize!
        assert canonical_contig("CHR22") == "chr22"
        assert canonical_contig("chr01") == "chr1"
        assert canonical_contig("CHRX") == "chrX"
        assert canonical_contig("chrX") == "chrX" and canonical_contig("chr2") == "chr2"

        # 1b. Dense pre-normalized inputs in preprocess="paired" (kills issparse/asarray branch mutants)
        #     AND exact padj tie at 0.0 (n_c=80, r~0.99 vs r~0.95) kills -(abs(zscore)) sort mutant!
        rng_tie = np.random.default_rng(910)
        n_tie = 80
        zt = rng_tie.normal(0.0, 1.0, size=n_tie)
        rna_tie = np.column_stack([2.0 + zt + rng_tie.normal(0, 0.01, n_tie)])
        # Peak 0 has tiny r (~0.02) placed BEFORE Peak 1 (slightly larger r ~0.12); BOTH have
        # p_value > 0.5 so BH padj clips BOTH to 1.0 (exact tie!), and secondary sort key
        # -(abs(zscore)) MUST put Peak 1 first!
        orth = rng_tie.normal(0.0, 1.0, size=n_tie)
        orth = orth - (np.dot(orth, zt) / np.dot(zt, zt)) * zt
        orth /= np.linalg.norm(orth)
        zt_u = zt / np.linalg.norm(zt)
        pk_t0 = 2.0 + 0.02 * zt_u + np.sqrt(1.0 - 0.02**2) * orth
        pk_t1 = 2.0 + 0.10 * zt_u + np.sqrt(1.0 - 0.10**2) * orth
        atac_tie = np.column_stack([pk_t0, pk_t1] + [rng_tie.normal(2.0, 1.0, n_tie) for _ in range(15)])
        atac_tie_names = ["chr1:12000-12400", "chr1:15000-15400"] + [
            f"chr2:{1000 + i*1000}-{1400 + i*1000}" for i in range(15)
        ]
        res_tie = compute_peak_to_gene_links(
            rna_matrix=rna_tie,
            rna_var_names=["GT"],
            atac_matrix=atac_tie,
            atac_var_names=atac_tie_names,
            gene_records={
                "GT": {"gene": "GT", "chrom": "chr1", "canon_chrom": "chr1", "tss": 10000, "strand": "+"},
                "G2": {"gene": "G2", "chrom": "chr2", "canon_chrom": "chr2", "tss": 1000, "strand": "+"},
            },
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="paired",
            rna_is_raw_counts=False,
            atac_is_raw_counts=False,
            background_null=True,
            n_bg_bins_per_axis=2,
            n_bg_per_bin=10,
            min_bg_peaks=5,
            seed=0,
        )
        assert res_tie["status"] == "success"
        assert res_tie["rna_normalization"] == "pre_normalized_input"
        assert res_tie["atac_normalization"] == "pre_normalized_input"
        assert len(res_tie["links"]) == 2
        assert res_tie["links"][0]["padj"] == res_tie["links"][1]["padj"]
        assert abs(res_tie["links"][0]["zscore"]) > abs(res_tie["links"][1]["zscore"])

        # 1c. Integer dense inputs in preprocess="as_is" & background_null=True with zero-variance
        #     peak + single-peak trans pool (idx_pool.size == 1 -> sd_bg fallback 0.0)
        rna_int = np.array([[i % 5 + 1] for i in range(20)], dtype=np.int32)
        atac_int = np.column_stack([
            np.array([i % 5 + 1 for i in range(20)], dtype=np.int32),
            np.zeros(20, dtype=np.int32),
            np.zeros(20, dtype=np.int32),
            np.zeros(20, dtype=np.int32),
            np.array([(i * 3) % 7 + 1 for i in range(20)], dtype=np.int32),
        ])
        res_int = compute_peak_to_gene_links(
            rna_matrix=rna_int,
            rna_var_names=["GI"],
            atac_matrix=atac_int,
            atac_var_names=[
                "chr1:12000-12400", "chr1:13000-13400", "chr1:14000-14400",
                "chr1:15000-15400", "chr2:1000-1400",
            ],
            gene_records={
                "GI": {"gene": "GI", "chrom": "chr1", "canon_chrom": "chr1", "tss": 10000, "strand": "+"},
                "G2": {"gene": "G2", "chrom": "chr2", "canon_chrom": "chr2", "tss": 1000, "strand": "+"},
            },
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="as_is",
            background_null=True,
            min_bg_peaks=5,
            seed=0,
        )
        assert res_int["status"] == "success"
        assert res_int["links"][0]["correlation"] == pytest.approx(1.0, abs=1e-12)

        # 2. Leading skipped entries in BOTH peak and gene loops kill `continue -> break` mutants,
        #    empty canon_chrom falls back to chrom, exact signed tss_offset_bp kills strand_sign=2,
        #    odd peak start+end kills // 2 vs / 2, and weaker-z peak placed BEFORE stronger-z peak
        #    with equal padj kills -(abs(zscore)) -> +(abs(zscore)) sort mutant!
        rng = np.random.default_rng(909)
        n_c = 50
        z = rng.normal(0.0, 1.0, size=n_c)
        # Gene order: [0] unannotated, [1] all-zero, [2] exact nnz==10 (dropped when n_cells>=50), [3] valid
        g_edge = np.zeros(n_c)
        g_edge[:10] = 5.0 + np.arange(10)
        g_valid = np.maximum(0.0, np.round(30.0 + 12.0 * z + rng.normal(0, 0.4, n_c)))
        rna_mat = np.column_stack([g_valid, np.zeros(n_c), g_edge, g_valid])
        rna_names = ["G_MISSING", "G_ZERO", "G_EDGE_10", "G_TARGET"]

        # Peak order: [0] unparseable, [1] nonstandard, [2] nnz==10 edge, [3] weaker + link, [4] stronger + link,
        #             [5..24] 20 trans background peaks on chr2 (odd end coord for // 2 check!)
        pk_edge = np.zeros(n_c)
        pk_edge[:10] = 5.0 + np.arange(10)
        pk_weak = np.maximum(0.0, np.round(30.0 + 4.0 * z + rng.normal(0, 3.0, n_c)))
        pk_strong = np.maximum(0.0, np.round(30.0 + 12.0 * z + rng.normal(0, 0.2, n_c)))
        bg_cols = [rng.poisson(8, size=n_c).astype(float) for _ in range(20)]
        atac_mat = np.column_stack([pk_strong, pk_strong, pk_edge, pk_weak, pk_strong] + bg_cols)
        atac_names = [
            "NOT_A_COORD",
            "chrUn_gl000220:10-200",
            "chr1:11000-11200",
            "chr1:12000-12401",  # midpoint = (12000+12401)//2 = 12200 (odd sum!), weaker z
            "chr1:15000-15400",  # midpoint = 15200, stronger z -> MUST sort first if padj ties
        ] + [f"chr2:{1000 + i*1000}-{1401 + i*1000}" for i in range(20)]

        gene_records = {
            "G_ZERO": {"gene": "G_ZERO", "chrom": "chr1", "canon_chrom": "", "tss": 10000, "strand": "+"},
            "G_EDGE_10": {"gene": "G_EDGE_10", "chrom": "chr1", "canon_chrom": "", "tss": 10000, "strand": "+"},
            "G_TARGET": {"gene": "G_TARGET", "chrom": "chr1", "canon_chrom": "", "tss": 10000, "strand": "+"},
            "G_CHR2": {"gene": "G_CHR2", "chrom": "chr2", "canon_chrom": "chr2", "tss": 1000, "strand": "+"},
        }
        res50 = compute_peak_to_gene_links(
            rna_matrix=rna_mat,
            rna_var_names=rna_names,
            atac_matrix=atac_mat,
            atac_var_names=atac_names,
            gene_records=gene_records,
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="as_is",
            min_cells_detected=10,
            standard_chroms_only=True,
            background_null=True,
            n_bg_bins_per_axis=2,
            n_bg_per_bin=10,
            min_bg_peaks=20,
            seed=0,
        )
        assert res50["status"] == "success"
        assert res50["rna_normalization"] == "as_is"
        assert res50["atac_normalization"] == "as_is"
        assert res50["n_unparseable_peaks"] == 1
        assert res50["n_nonstandard_peaks_dropped"] == 1
        # Only G_TARGET passes (G_MISSING, G_ZERO, G_EDGE_10 all skipped) and only 2 cis peaks pass
        assert res50["n_tested"] == 2
        assert {lk["gene"] for lk in res50["links"]} == {"G_TARGET"}
        by_pk = {lk["peak"]: lk for lk in res50["links"]}
        # Exact integer midpoint (12000+12401)//2 = 12200 -> exact tss_offset_bp == +2200
        assert by_pk["chr1:12000-12401"]["tss_offset_bp"] == 2200
        assert isinstance(by_pk["chr1:12000-12401"]["tss_offset_bp"], int)
        assert by_pk["chr1:15000-15400"]["tss_offset_bp"] == 5200
        # Stronger peak placed AFTER weaker peak in input MUST sort first!
        assert res50["links"][0]["peak"] == "chr1:15000-15400"

        # Same matrix truncated to n_cells=49: G_EDGE_10 (nnz=10 >= 9) and pk_edge (nnz=10 >= 9) ARE KEPT!
        res49 = compute_peak_to_gene_links(
            rna_matrix=sp.csr_matrix(rna_mat[:49]),
            rna_var_names=rna_names,
            atac_matrix=sp.csr_matrix(atac_mat[:49]),
            atac_var_names=atac_names,
            gene_records=gene_records,
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="as_is",
            min_cells_detected=10,
            standard_chroms_only=True,
            background_null=False,
        )
        assert res49["status"] == "success"
        assert res49["background_match"] == "none"
        assert res49["background_method"] == "none"
        assert {lk["gene"] for lk in res49["links"]} == {"G_TARGET", "G_EDGE_10"}
        assert res49["n_tested"] == 6  # 2 genes x 3 cis peaks

        # Boundary: background_null=True with <5 valid peaks (4 peaks) vs ==5 valid peaks + zero-count peak in as_is
        res_lt5 = compute_peak_to_gene_links(
            rna_matrix=np.column_stack([g_valid, np.full(n_c, 3.0)]),
            rna_var_names=["G_TARGET", "G_CONST"],
            atac_matrix=np.column_stack([pk_strong, pk_weak, np.zeros(n_c), bg_cols[0]]),
            atac_var_names=["chr1:15000-15400", "chr1:12000-12401", "chr2:100-500", "chr2:2000-2400"],
            gene_records={
                "G_TARGET": {"gene": "G_TARGET", "chrom": "chr1", "canon_chrom": "chr1", "tss": 10000, "strand": "+"},
                "G_CONST": {"gene": "G_CONST", "chrom": "chr1", "canon_chrom": "chr1", "tss": 10000, "strand": "+"},
            },
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="as_is",
            background_null=True,
            seed=0,
        )
        assert res_lt5["status"] == "success"
        assert res_lt5["background_method"] == "none"
        assert res_lt5["background_match"] == "none"

        # Exactly 5 valid peaks (including zero-count peak in as_is) enables background_null!
        # With 3 trans peaks (< 5), background_method is initialized even though per-link pool falls short.
        res_eq5 = compute_peak_to_gene_links(
            rna_matrix=np.column_stack([g_valid, np.full(n_c, 3.0)]),
            rna_var_names=["G_TARGET", "G_CONST"],
            atac_matrix=np.column_stack([pk_strong, pk_weak, np.zeros(n_c), bg_cols[0], bg_cols[1]]),
            atac_var_names=[
                "chr1:15000-15400", "chr1:12000-12401",
                "chr2:100-500", "chr2:2000-2400", "chr2:4000-4400",
            ],
            gene_records={
                "G_TARGET": {"gene": "G_TARGET", "chrom": "chr1", "canon_chrom": "chr1", "tss": 10000, "strand": "+"},
                "G_CONST": {"gene": "G_CONST", "chrom": "chr1", "canon_chrom": "chr1", "tss": 10000, "strand": "+"},
                "G_CHR2": {"gene": "G_CHR2", "chrom": "chr2", "canon_chrom": "chr2", "tss": 1000, "strand": "+"},
            },
            distance_bp=50_000,
            min_corr=0.0,
            preprocess="as_is",
            background_null=True,
            min_bg_peaks=5,
            seed=0,
        )
        assert res_eq5["status"] == "success"
        assert res_eq5["background_method"] == "binned_matched_background"
        assert res_eq5["background_match"] == "accessibility_and_length"
