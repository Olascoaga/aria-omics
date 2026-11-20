"""Wave 3A regression tests: chromatin / epigenomics / TOBIAS scientific engines.

One test (or a small group) per audited engine bug. Every test is pure / light
(numpy, pandas, scipy, anndata, scanpy) -- heavy tools (pysam, pyBigWig, pydeseq2,
snapatac2, bedtools) are replaced by minimal fakes so the regression guards run in
the base ``aria-env``.
"""
from __future__ import annotations

import ast
import gzip
import json
import math
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


# ── 1. chromatin_motifs ──────────────────────────────────────────────────────

def test_motif_as_float_rejects_nonfinite_and_bool():
    from aria.scripts.chromatin_motifs import _as_float

    for bad in (float("nan"), float("inf"), float("-inf"), "inf", "-inf", "nan",
                True, False, None, "x"):
        assert _as_float(bad) is None, bad
    assert _as_float("1e-5") == pytest.approx(1e-5)
    assert _as_float(0) == 0.0


def _snap_records():
    # SnapATAC2 column names; "adjusted p-value" deliberately precedes "p-value".
    def rec(mid, name, fam, lfc, p, padj):
        return {"id": mid, "name": name, "family": fam,
                "adjusted p-value": padj, "log2(fold change)": lfc, "p-value": p}
    return [
        rec("MA0001.1", "DEPL_ZERO", None, float("-inf"), 1e-60, 1e-50),  # bound_fg == 0
        rec("MA0002.1", "DEPL", "bZIP", -2.0, 1e-12, 1e-10),
        rec("MA0003.1", "ENR_STRONG", "ETS", 1.5, 1e-7, 1e-5),
        rec("MA0004.1", "ENR_WEAK", "ETS", 0.5, 0.05, 0.2),
        rec("MA0005.1", "ENR_NANPADJ", "ETS", 0.9, float("nan"), float("nan")),
        rec("MA0006.1", "ENR_INF", "ETS", float("inf"), 1e-9, 1e-8),  # bound_bg == 0
    ]


def test_motif_depleted_never_ranked_or_counted_as_enriched():
    from aria.scripts.chromatin_motifs import _normalize_motif_records

    res = _normalize_motif_records(_snap_records(), padj_max=0.05, top_n=10)
    names = [m["name"] for m in res["top"]]
    assert "DEPL_ZERO" not in names and "DEPL" not in names
    assert names[0] == "ENR_STRONG"
    assert all(m["log2_enrichment"] is not None and m["log2_enrichment"] > 0
               for m in res["top"])
    assert res["n_enriched"] == 1          # only ENR_STRONG is enriched AND padj<0.05
    assert res["n_depleted"] == 2          # -inf and -2.0, both significant
    assert res["n_nonfinite"] >= 2         # -inf / +inf lfc, NaN padj
    # strict JSON: no NaN / Infinity can leak out of the normalized records
    json.dumps(res["top"], allow_nan=False)


def test_motif_pvalue_column_exact_match_and_no_none_string():
    from aria.scripts.chromatin_motifs import (
        _detect_pvalue_col, _normalize_motif_records)

    assert _detect_pvalue_col(["adjusted p-value", "p-value"]) == "p-value"
    assert _detect_pvalue_col(["adjusted p-value"]) is None
    res = _normalize_motif_records(_snap_records(), 0.05, 10)
    strong = next(m for m in res["top"] if m["name"] == "ENR_STRONG")
    assert strong["pvalue"] == pytest.approx(1e-7)   # raw p, not the adjusted one
    assert strong["padj"] == pytest.approx(1e-5)
    rec = [{"id": "MA0099.1 ZNF410", "name": None, "family": None,
            "log2(fold change)": 1.0, "p-value": 0.01, "adjusted p-value": 0.02}]
    out = _normalize_motif_records(rec, 0.05, 5)["top"][0]
    assert out["family"] is None and out["name"] == "ZNF410"
    assert out["motif_id"] == "MA0099.1"


# ── 2. chromatin_footprint_tobias ────────────────────────────────────────────

def _write_bindetect(path, rows, a="A", b="B"):
    lines = [f"name\t{a}_{b}_change\t{a}_{b}_pvalue\ttotal_tfbs"]
    lines += [f"{n}\t{c}\t{p}\t{s}" for n, c, p, s in rows]
    path.write_text("\n".join(lines) + "\n")


def test_tobias_summarize_bindetect_direction_strict_and_finite(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import summarize_bindetect

    txt = tmp_path / "bindetect_results.txt"
    _write_bindetect(txt, [
        ("TF1", 0.40, 1e-5, 100), ("TF2", 0.20, 1e-4, 100), ("TF3", 0.10, 1e-3, 100),
        ("TFNAN", "nan", 1e-6, 100), ("TFINF", "inf", 1e-6, 100),
        ("TFPINF", 0.30, "inf", 100),
    ])
    res = summarize_bindetect(str(txt), "A", "B", top_n=15)
    assert [r["tf"] for r in res["top_toward_A"]] == ["TF1", "TF2", "TF3"]
    # unidirectional perturbation: nothing may be listed toward B
    assert res["top_toward_B"] == []
    assert res["n_ranked_candidates"] == 3
    json.dumps(res, allow_nan=False)


def test_tobias_replicate_inference_direction_strict():
    from aria.scripts.chromatin_footprint_tobias import infer_replicate_footprints

    reps = {f"a{i}": "A" for i in range(3)} | {f"b{i}": "B" for i in range(3)}
    rng = np.random.default_rng(1)
    table = {}
    for k in range(30):
        shift = 0.5 if k < 6 else 0.0          # only "toward A" effects exist
        scores = {f"a{i}": 1.0 + shift + rng.normal(0, 0.01) for i in range(3)}
        scores |= {f"b{i}": 1.0 + rng.normal(0, 0.01) for i in range(3)}
        table[f"M{k}"] = {"tf": f"TF{k}", "motif_id": f"M{k}", "scores": scores,
                          "n_sites": 50}
    res = infer_replicate_footprints(table, reps, "A", "B", top_n=15)
    assert all(r["change"] > 0 for r in res["top_toward_A"])
    assert all(r["change"] < 0 for r in res["top_toward_B"])
    assert {r["tf"] for r in res["top_toward_A"]} <= {f"TF{k}" for k in range(6)}


def test_tobias_descriptive_fallback_direction_strict():
    from aria.scripts.chromatin_footprint_tobias import infer_replicate_footprints

    reps = {"a0": "A", "a1": "A", "b0": "B", "b1": "B"}  # < 3 reps -> descriptive
    table = {f"M{k}": {"tf": f"TF{k}", "scores": {"a0": 2.0 + k, "a1": 2.0 + k,
                                                  "b0": 1.0, "b1": 1.0}}
             for k in range(4)}
    res = infer_replicate_footprints(table, reps, "A", "B", top_n=15)
    assert res["inference"]["status"] == "descriptive_only"
    assert len(res["top_toward_A"]) == 4
    assert res["top_toward_B"] == []


# ── 3. chromatin_lsi_clustering ──────────────────────────────────────────────

def _lsi_input(tmp_path, n=300, p=400):
    import anndata as ad
    from scipy import sparse

    rng = np.random.default_rng(0)
    X = rng.poisson(0.3, size=(n, p)).astype(float)
    X[:10] = rng.poisson(3.0, size=(10, p))          # high depth + features
    a = ad.AnnData(sparse.csr_matrix(X))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.var_names = [f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(p)]
    path = tmp_path / "in.h5ad"
    a.write_h5ad(path)
    return path


@pytest.mark.parametrize("max_cells,sketched", [(250, True), (295, False)])
def test_lsi_sketch_uses_post_doublet_cell_count(tmp_path, max_cells, sketched):
    pytest.importorskip("scanpy")
    from aria.scripts.chromatin_lsi_clustering import chromatin_lsi_clustering

    path = _lsi_input(tmp_path)
    res = chromatin_lsi_clustering({
        "data_path": str(path), "output_dir": str(tmp_path / f"out{max_cells}"),
        "max_cells": max_cells, "n_components": 10,
    })
    assert res["status"] == "success", res
    assert res["doublets"]["removed"] == 10
    sketch_msgs = [w for w in res["warnings"] if w.startswith("Sketched")]
    if sketched:
        assert sketch_msgs and f"of 290 cells" in sketch_msgs[0]
    else:
        # 290 cells remain after doublet removal (<= 295): no sketch at all.
        assert not sketch_msgs


# ── 4. chromatin_qc ──────────────────────────────────────────────────────────

def test_mito_fraction_streams_full_sorted_file(tmp_path):
    from aria.scripts.chromatin_qc import _compute_mito_fraction

    frag = tmp_path / "fragments.tsv.gz"
    with gzip.open(frag, "wt") as fh:
        fh.write("# comment\n")
        for i in range(250_000):
            fh.write(f"chr1\t{i}\t{i + 100}\tBC{i % 50}\t1\n")
        for i in range(50_000):
            fh.write(f"chrM\t{i}\t{i + 100}\tBC{i % 50}\t1\n")
    assert _compute_mito_fraction(str(frag), "chrM") == pytest.approx(50_000 / 300_000)

    empty = tmp_path / "empty.tsv.gz"
    with gzip.open(empty, "wt") as fh:
        fh.write("# only a header\n")
    assert _compute_mito_fraction(str(empty), "chrM") is None


def test_mito_fraction_ensembl_mt_contig(tmp_path):
    from aria.scripts.chromatin_qc import _compute_mito_fraction

    frag = tmp_path / "fragments.tsv"
    frag.write_text("1\t0\t100\tA\t1\n1\t5\t100\tA\t1\nMT\t0\t50\tA\t1\n")
    assert _compute_mito_fraction(str(frag), "chrMT") == pytest.approx(1 / 3)


def test_frip_contig_normalization_ucsc_vs_ensembl(tmp_path):
    from aria.scripts.chromatin_qc import (
        _frip_from_intervals, _frip_per_barcode_from_intervals,
        _peaks_by_chrom_from_bed)

    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t100\t200\nchrX\t500\t600\n")
    peaks = _peaks_by_chrom_from_bed(str(bed))
    frip, n = _frip_from_intervals([("1", 150, 160), ("X", 550, 560), ("2", 1, 5)],
                                   peaks)
    assert n == 3 and frip == pytest.approx(2 / 3)
    dist = _frip_per_barcode_from_intervals(
        [("1", 150, 160, "bc"), ("1", 900, 950, "bc")], peaks, min_frags=1)
    assert dist == [pytest.approx(0.5)]


def test_bam_tss_enrichment_contig_normalized(monkeypatch, tmp_path):
    from contextlib import contextmanager

    from aria.scripts import chromatin_peak_counts, chromatin_qc

    gtf = tmp_path / "genes.gtf"
    gtf.write_text('1\tt\tgene\t3001\t4000\t.\t+\t.\tgene_id "g";\n')  # Ensembl "1"
    tss = chromatin_qc._read_local_tss(str(gtf))
    assert "1" in tss                       # caller-visible keys stay raw

    class _Read:
        is_unmapped = is_secondary = is_supplementary = is_duplicate = False
        is_reverse = is_paired = is_proper_pair = is_read1 = False
        template_length = 0
        query_length = 30

        def __init__(self, start):
            self.reference_name = "chr1"    # UCSC-named BAM
            self.reference_start = start
            self.reference_end = start + 30

    class _Header:
        def to_dict(self):
            return {}

    class _Bam:
        header = _Header()

        def fetch(self, until_eof=True):
            # cut = start + 4: two central cuts, two flank cuts (TSS 0-based 3000)
            return iter([_Read(2996), _Read(2990), _Read(1046), _Read(4946)])

    @contextmanager
    def fake_open(*a, **k):
        yield _Bam()

    monkeypatch.setattr(chromatin_peak_counts, "_open_alignment", fake_open)
    metrics, _ = chromatin_qc._scan_bam_metrics("x.bam", tss=tss)
    assert metrics["tss_enrichment"] == pytest.approx((2 / 101) / (2 / 200))
    assert chromatin_qc._norm_contig("chrM") == chromatin_qc._norm_contig("MT") == "MT"


def test_scatac_pass_qc_measured_failure_is_false():
    from aria.scripts.chromatin_qc import _scatac_pass_qc

    assert _scatac_pass_qc(["Mito 30% > 10%"], qc_complete=False) is False
    assert _scatac_pass_qc(["TSS < 4"], qc_complete=True) is False
    assert _scatac_pass_qc([], qc_complete=False) is None
    assert _scatac_pass_qc([], qc_complete=True) is True


def test_strand_cc_scores_stable_baseline():
    from aria.scripts.chromatin_qc import _strand_cc_scores

    n_bins = 1_000_000            # null SD of the CC ~ 1e-3
    cc = np.full(301, 0.02)
    cc[40] = 0.08                 # fragment peak (200 bp)
    cc[15] = 0.04                 # phantom (read length)
    nsc, rsc = _strand_cc_scores(cc, 40, 15, 70, 5, n_bins)
    assert nsc == pytest.approx(4.0) and rsc == pytest.approx(3.0)

    # Baseline indistinguishable from noise (2e-6) -> undefined, not NSC ~ 40,000.
    cc2 = cc.copy()
    cc2[200:] = 2e-6
    nsc2, rsc2 = _strand_cc_scores(cc2, 40, 15, 70, 5, n_bins)
    assert nsc2 is None
    assert rsc2 is None or math.isfinite(rsc2)

    # Phantom below the baseline -> no inverted (negative) RSC.
    cc3 = cc.copy()
    cc3[15] = 0.01
    _, rsc3 = _strand_cc_scores(cc3, 40, 15, 70, 5, n_bins)
    assert rsc3 is None


def _install_fake_pysam(monkeypatch):
    monkeypatch.setitem(sys.modules, "pysam", types.ModuleType("pysam"))


def test_bulk_qc_enforces_min_nrf_mito_none_and_unique_stems(monkeypatch, tmp_path):
    from aria.scripts import chromatin_qc

    _install_fake_pysam(monkeypatch)

    def fake_scan(path, *a, **k):
        return ({"file": path, "total_reads": 0, "mapped_reads": 0,
                 "mapped_nonduplicate_reads": 0, "mito_frac": None,
                 "dup_rate": None, "nrf": 0.2, "pbc1": None, "pbc2": None,
                 "nsc": None, "rsc": None, "estimated_fragment_length_bp": None,
                 "tss_enrichment": None, "tss_reads_sampled": 0}, [])

    monkeypatch.setattr(chromatin_qc, "_scan_bam_metrics", fake_scan)
    files = [str(tmp_path / "r1" / "rep.bam"), str(tmp_path / "r2" / "rep.bam")]
    res = chromatin_qc._bulk_chromatin_qc(files, "ChIP", "hg38", "tf", [], params={})
    assert res["status"] == "success"
    assert res["mito_fraction"] is None                 # never a fabricated 0.0
    assert res["pass_qc"] is False                      # NRF 0.2 < min_nrf 0.5
    assert any("NRF" in f for f in res["qc_failures"])
    assert len(res["nrf_per_sample"]) == 2              # no stem overwrite


def test_scan_bam_metrics_mito_none_without_mapped_reads(monkeypatch):
    from contextlib import contextmanager

    from aria.scripts import chromatin_peak_counts, chromatin_qc

    class _Header:
        def to_dict(self):
            return {}

    class _Bam:
        header = _Header()

        def fetch(self, until_eof=True):
            return iter(())

    @contextmanager
    def fake_open(*a, **k):
        yield _Bam()

    monkeypatch.setattr(chromatin_peak_counts, "_open_alignment", fake_open)
    metrics, _ = chromatin_qc._scan_bam_metrics("x.bam")
    assert metrics["mapped_reads"] == 0 and metrics["mito_frac"] is None


# ── 5. chromatin_bulk_diffacc ────────────────────────────────────────────────

def _bulk_da_inputs(tmp_path, levels=("trt/x", "ctrl")):
    samples = [f"s{i}" for i in range(6)]
    rng = np.random.default_rng(3)
    counts = pd.DataFrame(rng.poisson(50, size=(30, 6)),
                          index=[f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(30)],
                          columns=samples)
    counts.index.name = "peak"
    cpath = tmp_path / "counts.tsv"
    counts.to_csv(cpath, sep="\t")
    meta = pd.DataFrame({"sample_id": samples,
                         "condition": [levels[0]] * 3 + [levels[1]] * 3,
                         "replicate": [f"r{i}" for i in range(6)]})
    mpath = tmp_path / "meta.tsv"
    meta.to_csv(mpath, sep="\t", index=False)
    return cpath, mpath, counts


def test_bulk_da_spikein_requested_without_data_fails_closed(tmp_path):
    from aria.scripts.chromatin_bulk_diffacc import chromatin_bulk_diffacc

    cpath, mpath, _ = _bulk_da_inputs(tmp_path)
    res = chromatin_bulk_diffacc({
        "data_type": "CUT_AND_TAG", "counts_matrix_path": str(cpath),
        "sample_metadata_path": str(mpath),
        "comparisons": [{"test": "trt/x", "reference": "ctrl"}],
        "normalization_mode": "heterologous_spikein",
        "output_dir": str(tmp_path / "da"),
    })
    assert res["ran"] is False
    assert res["data_type"] == "CUT_AND_TAG"
    assert res["normalization_mode_requested"] == "heterologous_spikein"
    assert "median-of-ratios" in res["reason"]


def test_bulk_da_emits_normalization_mode_and_sanitized_paths(monkeypatch, tmp_path):
    from aria.scripts import chromatin_bulk_diffacc as mod
    from aria.scripts import rna_bulk_de

    cpath, mpath, counts = _bulk_da_inputs(tmp_path)

    def fake_deseq2(c, m, *a, **k):
        df = pd.DataFrame({"log2FoldChange": np.linspace(-2, 2, len(c)),
                           "padj": np.full(len(c), 0.01),
                           "pvalue": np.full(len(c), 0.001),
                           "baseMean": np.full(len(c), 50.0)}, index=c.index)
        return {"status": "success", "results": df, "n_replicates": 6}, []

    monkeypatch.setattr(rna_bulk_de, "_run_deseq2", fake_deseq2)
    monkeypatch.setattr(mod, "compute_counts_lobo_jackknife",
                        lambda *a, **k: {"status": "skipped"})
    out = tmp_path / "da"
    res = mod.chromatin_bulk_diffacc({
        "data_type": "ChIP", "counts_matrix_path": str(cpath),
        "sample_metadata_path": str(mpath),
        "comparisons": [{"test": "trt/x", "reference": "ctrl"}],
        "output_dir": str(out),
    })
    assert res["status"] == "success" and res["ran"] is True, res
    assert res["normalization_mode"] == "deseq2_median_of_ratios"
    csv_path = Path(res["comparisons"][0]["full_results_csv"])
    assert csv_path.parent == out and csv_path.is_file()
    assert "/" not in csv_path.name.replace("bulk_atac_da_", "")


def test_bulk_da_skip_carries_data_type(tmp_path):
    from aria.scripts.chromatin_bulk_diffacc import chromatin_bulk_diffacc

    cpath, mpath, counts = _bulk_da_inputs(tmp_path)
    (counts * 0).to_csv(cpath, sep="\t")
    res = chromatin_bulk_diffacc({
        "data_type": "CUT_AND_RUN", "counts_matrix_path": str(cpath),
        "sample_metadata_path": str(mpath),
        "comparisons": [{"test": "trt/x", "reference": "ctrl"}],
        "output_dir": str(tmp_path / "da"),
    })
    assert res["ran"] is False and res["data_type"] == "CUT_AND_RUN"


# ── 6. chromatin_diffacc ─────────────────────────────────────────────────────

def test_per_cluster_da_rejects_infinite_stats(monkeypatch):
    import anndata as ad
    import scanpy as sc

    from aria.scripts.chromatin_diffacc import _per_cluster_accessibility

    rng = np.random.default_rng(0)
    a = ad.AnnData(rng.poisson(1.0, size=(20, 4)).astype(float))
    a.var_names = ["p0", "p1", "p2", "p3"]
    a.obs["leiden"] = ["0"] * 10 + ["1"] * 10

    def fake_rgg(adata, groupby, method, key_added):
        dt2 = [("0", "O"), ("1", "O")]
        dtf = [("0", "f8"), ("1", "f8")]
        names = np.array([("p0", "p0"), ("p1", "p1"), ("p2", "p2")], dtype=dt2)
        lfc = np.array([(np.inf, 2.0), (3.0, -np.inf), (2.0, 2.0)], dtype=dtf)
        padj = np.array([(1e-5, np.inf), (1e-5, 1e-5), (np.nan, 1e-6)], dtype=dtf)
        adata.uns[key_added] = {"names": names, "logfoldchanges": lfc,
                                "pvals_adj": padj}

    monkeypatch.setattr(sc.tl, "rank_genes_groups", fake_rgg)
    res = _per_cluster_accessibility(a, "leiden", 0.05, 0.5, 10, [])
    for row in res["_rows"]:
        assert math.isfinite(row["log2fc"]) and math.isfinite(row["padj"])
    assert {(r["cluster"], r["peak"]) for r in res["_rows"]} == {
        ("0", "p1"), ("1", "p2")}


def test_pseudobulk_aggregation_with_duplicate_barcodes(monkeypatch):
    import anndata as ad
    from scipy import sparse

    from aria.scripts import rna_bulk_de
    from aria.scripts.chromatin_diffacc import _pseudobulk_da

    X = np.ones((8, 3))
    X[:, 0] = np.arange(1, 9)           # distinct per-cell counts in peak 0
    a = ad.AnnData(sparse.csr_matrix(X))
    a.obs_names = ["AAAC-1", "AAAC-1", "GGGT-1", "GGGT-1"] * 2   # duplicates
    a.var_names = ["p0", "p1", "p2"]
    a.obs["condition"] = ["t"] * 4 + ["c"] * 4
    a.obs["donor"] = ["d1", "d1", "d2", "d2", "d3", "d3", "d4", "d4"]
    captured = {}

    def fake_deseq2(counts, meta, *a_, **k):
        captured["counts"] = counts.copy()
        return {"status": "error", "error_type": "Fake", "details": "x"}, []

    monkeypatch.setattr(rna_bulk_de, "_run_deseq2", fake_deseq2)
    _pseudobulk_da(a, "condition", "donor", [{"test": "t", "reference": "c"}],
                   0.05, 0.5, 10, 1, 2, False, [])
    p0 = captured["counts"].loc["p0"].to_dict()
    assert p0 == {"d1": 3.0, "d2": 7.0, "d3": 11.0, "d4": 15.0}


# ── 7. chromatin_peaks + seacr_idr ───────────────────────────────────────────

@pytest.mark.parametrize("genome,expected", [
    ("hg38", "hs"), ("HG19", "hs"), ("GRCh38", "hs"), ("grch37", "hs"),
    ("mm10", "mm"), ("MM10", "mm"), ("mm9", "mm"), ("mm39", "mm"),
    ("GRCm39", "mm"), ("grcm38", "mm"), ("dm6", "dm"),
])
def test_genome_size_case_insensitive_aliases(genome, expected):
    from aria.scripts.chromatin_peaks import _get_genome_size

    assert _get_genome_size(genome) == expected


def test_effective_genome_size_bp():
    from aria.scripts.chromatin_peaks import _effective_genome_size_bp

    assert _effective_genome_size_bp("GRCm39") == 1_870_000_000
    assert _effective_genome_size_bp("hg38") == 2_700_000_000
    assert _effective_genome_size_bp("hg38", 2_913_022_398) == 2_913_022_398
    assert _effective_genome_size_bp("unknownAsm") is None


def test_peak_count_warning_not_inverted():
    src = (ROOT / "aria/scripts/chromatin_peaks.py").read_text()
    assert "Consider increasing q-value threshold" not in src
    assert "more stringent (smaller) q-value" in src


def test_peaks_cache_key_includes_all_caller_params(monkeypatch, tmp_path):
    from aria.scripts import chromatin_peaks

    bam = tmp_path / "ip.bam"
    bam.write_bytes(b"x")

    class _Captured(Exception):
        pass

    def capture(manifest, expected_inputs, expected_parameters):
        raise _Captured(expected_parameters)

    monkeypatch.setattr(chromatin_peaks, "load_valid_stage_result", capture)
    monkeypatch.setattr(chromatin_peaks, "file_record", lambda f: {"path": str(f)})
    with pytest.raises(_Captured) as exc:
        chromatin_peaks.chromatin_peaks({
            "data_type": "ChIP", "files": [str(bam)], "peak_caller": "broad_domain",
            "output_dir": str(tmp_path / "pk"), "control_scale_factors": [0.5],
            "domain_bin_size": 5000, "domain_min_log2fc": 1.0,
            "domain_max_gap_bp": 10_000, "domain_min_bp": 15_000,
        })
    params = exc.value.args[0]
    assert params["control_scale_factors"] == [0.5]
    assert params["domain_bin_size"] == 5000
    assert params["domain_min_log2fc"] == 1.0
    assert params["domain_max_gap_bp"] == 10_000
    assert params["domain_min_bp"] == 15_000


def test_overlap_merge_half_open_book_ended(tmp_path):
    from aria.scripts.chromatin_peaks import _merge_intervals_with_support

    r1 = tmp_path / "r1.bed"
    r2 = tmp_path / "r2.bed"
    r1.write_text("chr1\t100\t200\n")
    r2.write_text("chr1\t200\t300\n")
    merged = _merge_intervals_with_support([str(r1), str(r2)])
    assert len(merged) == 2 and all(len(m[3]) == 1 for m in merged)


def test_blacklist_refresh_never_overwrites_overlap_with_idr(tmp_path):
    from aria.scripts.chromatin_peaks import _refresh_consensus_counts

    overlap_bed = tmp_path / "reproducible_consensus.narrowPeak"
    idr_bed = tmp_path / "idr_reproducible_consensus.narrowPeak"
    overlap_bed.write_text("".join(f"chr1\t{i*10}\t{i*10+5}\n" for i in range(10)))
    idr_bed.write_text("chr1\t0\t5\nchr1\t10\t15\nchr1\t20\t25\n")
    repro = {
        "overlap": {"ran": True, "n_reproducible_regions": 10,
                    "n_candidate_regions": 20, "fraction_reproducible": 0.5,
                    "output_path": str(overlap_bed)},
        "idr": {"ran": True, "n_reproducible_peaks": 4, "output_path": str(idr_bed)},
    }
    _refresh_consensus_counts(repro, str(idr_bed))
    assert repro["idr"]["n_reproducible_peaks"] == 3
    assert repro["overlap"]["n_reproducible_regions"] == 10
    assert repro["overlap"]["fraction_reproducible"] == 0.5

    _refresh_consensus_counts(repro, str(overlap_bed))
    assert repro["overlap"]["n_reproducible_regions"] == 10


def test_bedgraph_value_format_keeps_tiny_scores():
    from aria.utils.seacr_idr import _fmt_bedgraph_value

    assert _fmt_bedgraph_value(3e-8) != "0"
    assert float(_fmt_bedgraph_value(3e-8)) == pytest.approx(3e-8)
    assert _fmt_bedgraph_value(3.0) == "3"
    assert float(_fmt_bedgraph_value(1.23456789012)) == pytest.approx(1.23456789, 1e-9)


def test_fragment_bedgraph_tiny_scale_not_zeroed(monkeypatch, tmp_path):
    from contextlib import contextmanager

    from aria.utils import seacr_idr

    class _Read:
        is_unmapped = is_secondary = is_supplementary = is_duplicate = False
        is_paired = False
        is_reverse = False
        mapping_quality = 60
        reference_name = "chr1"
        reference_start = 100
        reference_end = 150

    class _Header:
        references = ("chr1",)
        lengths = (10_000,)

    class _Bam:
        header = _Header()

        def fetch(self, until_eof=True):
            return iter([_Read(), _Read()])

    @contextmanager
    def fake_open(*a, **k):
        yield _Bam()

    monkeypatch.setattr(seacr_idr, "_open_alignment", fake_open)
    out = tmp_path / "x.bedgraph"
    seacr_idr.bam_to_fragment_bedgraph("x.bam", out, scale_factor=1e-8)
    vals = [float(l.split("\t")[3]) for l in out.read_text().splitlines()]
    assert vals and all(v == pytest.approx(2e-8) for v in vals)
    blocks = seacr_idr._integrate_bedgraph_signal_blocks(out)
    assert blocks, "tiny spike-in-scaled signal must survive as a SEACR block"


def test_scored_peaks_drop_nonfinite(tmp_path):
    from aria.utils.seacr_idr import _read_scored_peaks

    bed = tmp_path / "p.narrowPeak"
    bed.write_text(
        "chr1\t0\t100\tp1\t0\t.\t5.0\t-1\t-1\t50\n"
        "chr1\t200\t300\tp2\t0\t.\tnan\t-1\t-1\t50\n"
        "chr1\t400\t500\tp3\t0\t.\tinf\t-1\t-1\t50\n")
    peaks = _read_scored_peaks(bed)
    assert [p[1] for p in peaks] == [0]


def test_idr_merge_half_open_and_full_precision_local_idr(monkeypatch, tmp_path):
    from aria.utils import seacr_idr

    rows = "".join(
        f"chr1\t{i * 1000}\t{i * 1000 + 100}\tp\t0\t.\t{10 + i}\t-1\t-1\t50\n"
        f"chr1\t{i * 1000 + 100}\t{i * 1000 + 200}\tq\t0\t.\t{20 + i}\t-1\t-1\t50\n"
        for i in range(12))
    r1 = tmp_path / "r1.narrowPeak"
    r2 = tmp_path / "r2.narrowPeak"
    r1.write_text(rows)
    r2.write_text(rows)

    def fake_fit(s1, s2):
        n = len(s1)
        return {"global_idr": np.full(n, 1e-4), "local_idr": np.full(n, 1e-9),
                "converged": True, "n_iter": 1, "params": {}}

    monkeypatch.setattr(seacr_idr, "fit_copula_idr_pair", fake_fit)
    res = seacr_idr.compute_idr_reproducible_peaks(
        [str(r1), str(r2)], tmp_path / "idr.narrowPeak")
    assert res["ran"] is True
    # 24 book-ended peaks stay 24 (no merge of [a,b) with [b,c)).
    assert res["n_reproducible_peaks"] == 24
    assert res["pair_reports"][0]["median_local_idr_passing"] == pytest.approx(1e-9)


def test_broad_domains_reject_invalid_scale_factors():
    from aria.utils.seacr_idr import call_broad_domains

    for bad in ([float("nan")], [0.0], [-1.0], [float("inf")]):
        res = call_broad_domains(["ip.bam"], "/tmp/never.bed", ip_scale_factors=bad)
        assert res["status"] == "error" and res["error_type"] == "InvalidScaleFactors"
    res = call_broad_domains(["ip.bam"], "/tmp/never.bed",
                             control_bams=["c.bam"], control_scale_factors=[float("nan")])
    assert res["error_type"] == "InvalidScaleFactors"


def _fake_bigwig(monkeypatch, store):
    class _BW:
        def addHeader(self, header):
            store["header"] = header

        def addEntries(self, chroms, starts, ends=None, values=None):
            store["values"] = list(values)

        def close(self):
            pass

    mod = types.ModuleType("pyBigWig")
    mod.open = lambda path, mode="r": _BW()
    monkeypatch.setitem(sys.modules, "pyBigWig", mod)


def test_rpgc_bigwig_uses_effective_genome_size_and_real_fragment_bp(monkeypatch, tmp_path):
    from aria.utils import seacr_idr

    store = {}
    _fake_bigwig(monkeypatch, store)

    def fake_bg(bam, out, **k):
        # 2 fragments of 50 bp each -> 100 covered bp in total
        Path(out).write_text("chr1\t100\t150\t2\n")
        return {"n_fragments": 2, "chrom_lengths": {"chr1": 10_000}}

    monkeypatch.setattr(seacr_idr, "bam_to_fragment_bedgraph", fake_bg)
    info = seacr_idr.write_normalized_bigwig(
        "x.bam", tmp_path / "x.bw", normalization="RPGC", effective_genome_size=1000)
    assert info["normalization"] == "RPGC"
    assert info["scale_factor"] == pytest.approx(1000 / 100)
    assert info["mean_fragment_bp"] == pytest.approx(50.0)
    assert store["values"] == [pytest.approx(20.0)]

    with pytest.raises(ValueError, match="effective_genome_size"):
        seacr_idr.write_normalized_bigwig("x.bam", tmp_path / "y.bw",
                                          normalization="RPGC")


def test_chromatin_peaks_passes_effective_genome_size_to_bigwig():
    tree = ast.parse((ROOT / "aria/scripts/chromatin_peaks.py").read_text())
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
             and getattr(n.func, "id", None) == "write_normalized_bigwig"]
    assert calls and all(any(k.arg == "effective_genome_size" for k in c.keywords)
                         for c in calls)


# ── 8. chromatin_peak_annotation ─────────────────────────────────────────────

def _models(genes):
    by = {}
    for g in genes:
        by.setdefault("1", []).append(g)
    for c in by:
        by[c].sort(key=lambda g: g["tss"])
    return by, {c: [g["tss"] for g in gl] for c, gl in by.items()}


def test_annotation_exonic_beats_intronic_across_overlapping_genes():
    from aria.scripts.chromatin_peak_annotation import annotate_one_peak

    a = {"gene": "LONG", "strand": "+", "start": 1000, "end": 50_000, "tss": 1000,
         "exons": [(1000, 1100), (49_000, 50_000)]}
    b = {"gene": "INNER", "strand": "+", "start": 15_000, "end": 30_000,
         "tss": 15_000, "exons": [(15_000, 15_100), (19_900, 20_100)]}
    models, idx = _models([a, b])
    ann = annotate_one_peak("chr1", 19_950, 20_050, models, idx, 2000, 2000)
    assert ann["feature"] == "Exonic" and ann["nearest_gene"] == "INNER"


def test_annotation_promoter_picks_closest_tss():
    from aria.scripts.chromatin_peak_annotation import annotate_one_peak

    left = {"gene": "LEFT", "strand": "-", "start": 500, "end": 1000, "tss": 1000,
            "exons": [(500, 1000)]}
    right = {"gene": "RIGHT", "strand": "+", "start": 1300, "end": 5000,
             "tss": 1300, "exons": [(1300, 1400)]}
    models, idx = _models([left, right])
    ann = annotate_one_peak("1", 1200, 1300, models, idx, 2000, 2000)
    assert ann["feature"] == "Promoter" and ann["nearest_gene"] == "RIGHT"
    assert ann["distance_to_tss"] == -50


def _big_gtf(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#" + "x" * 120_000 + "\n")


def test_resolve_gtf_assembly_aliases_and_no_cross_species(monkeypatch, tmp_path):
    from aria.scripts.chromatin_peak_annotation import _resolve_gtf

    store = tmp_path / "genomes"
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("ARIA_GENOME_DIR", str(store))
    monkeypatch.delenv("ARIA_GTF", raising=False)
    monkeypatch.delenv("ARIA_HOME", raising=False)
    _big_gtf(store / "hg38" / "annotation.gtf")
    # Mouse explicitly requested, only a human GTF exists -> None (never hg38).
    assert _resolve_gtf({"genome": "mm10"}, tmp_path) is None
    assert _resolve_gtf({"genome": "GRCm39"}, tmp_path) is None
    _big_gtf(store / "mm39" / "annotation.gtf")
    got = _resolve_gtf({"genome": "GRCm39"}, tmp_path)
    assert got is not None and got.parent.name == "mm39"
    got = _resolve_gtf({"genome": "GRCh38"}, tmp_path)
    assert got is not None and got.parent.name == "hg38"


def test_select_comparison_rows_never_falls_back_to_other_contrasts():
    from aria.scripts.chromatin_peak_annotation import _select_comparison_rows

    multi = pd.DataFrame({"peak": ["a", "b"], "comparison": ["A_vs_B", "C_vs_D"]})
    assert _select_comparison_rows(multi, {"test": "X", "reference": "Y"}).empty
    sel = _select_comparison_rows(multi, {"test": "A", "reference": "B"})
    assert sel["peak"].tolist() == ["a"]
    single = pd.DataFrame({"peak": ["a"], "comparison": ["old_name"]})
    assert len(_select_comparison_rows(single, {"test": "X", "reference": "Y"})) == 1


# ── 9. chromatin_peak_counts ─────────────────────────────────────────────────

def test_peak_contigs_missing_from_genome_order(tmp_path):
    from aria.scripts.chromatin_peak_counts import _peak_contigs_missing

    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t0\t10\nchrUn_KI270\t0\t10\n# c\nchr1\t20\t30\n")
    assert _peak_contigs_missing(bed, ["chr1", "chr2"]) == ["chrUn_KI270"]
    assert _peak_contigs_missing(bed, ["chr1", "chrUn_KI270"]) == []


def test_peak_counts_guards_sorted_mode_on_missing_contigs():
    src = (ROOT / "aria/scripts/chromatin_peak_counts.py").read_text()
    assert "_peak_contigs_missing(peaks_path" in src


# ── 10. chromatin_regulatory ─────────────────────────────────────────────────

def test_regulatory_rna_obs_set_hoisted_out_of_comprehension():
    tree = ast.parse((ROOT / "aria/scripts/chromatin_regulatory.py").read_text())
    for comp in (n for n in ast.walk(tree) if isinstance(n, ast.ListComp)):
        for call in (n for n in ast.walk(comp) if isinstance(n, ast.Call)):
            assert getattr(call.func, "id", None) != "set", ast.unparse(comp)


# ── 11. epigenomic_spikein_de ────────────────────────────────────────────────

def test_spikein_all_nan_count_column_falls_through_to_scale_column():
    from aria.utils.epigenomic_spikein_de import resolve_replicate_spikein_size_factors

    meta = pd.DataFrame({"condition": ["a", "a", "b", "b"],
                         "replicate": ["1", "2", "1", "2"],
                         "spikein_reads": [np.nan] * 4,
                         "spikein_scale_factor": [1.0, 2.0, 1.0, 0.5]},
                        index=["s1", "s2", "s3", "s4"])
    rep_meta = pd.DataFrame({"condition": ["a", "a", "b", "b"]},
                            index=["a::1", "a::2", "b::1", "b::2"])
    sf, info = resolve_replicate_spikein_size_factors(meta, rep_meta, "condition",
                                                      "replicate")
    assert info["used"] is True
    assert info["source"] == "metadata_scale_column:spikein_scale_factor"
    assert np.isfinite(sf.to_numpy()).all()


@pytest.mark.parametrize("bad", [np.inf, -1.0, 0.0, np.nan])
def test_external_size_factors_must_be_finite_positive(bad):
    from aria.utils.epigenomic_spikein_de import run_deseq2_with_external_size_factors

    idx = [f"s{i}" for i in range(6)]
    counts = pd.DataFrame(np.random.default_rng(0).poisson(20, size=(10, 6)),
                          columns=idx, index=[f"g{i}" for i in range(10)])
    meta = pd.DataFrame({"condition": ["t"] * 3 + ["c"] * 3}, index=idx)
    sf = pd.Series([1.0, 1.0, bad, 1.0, 1.0, 1.0], index=idx)
    res, _ = run_deseq2_with_external_size_factors(
        counts, meta, sf, "condition", "t", "c", 0.05, 0.5)
    assert res["status"] == "error"
    assert res["error_type"] == "InvalidSpikeInSizeFactors"


# ── 12. chromatin_robustness ─────────────────────────────────────────────────

def test_row_features_ignore_explicit_sparse_zeros():
    from scipy import sparse

    from aria.utils.chromatin_robustness import _row_sums_and_features

    data = np.array([1.0, 0.0, 0.0, 2.0])
    x = sparse.csr_matrix((data, ([0, 0, 0, 1], [0, 1, 2, 0])), shape=(2, 3))
    assert x.nnz == 4                               # explicit zeros are stored
    totals, features = _row_sums_and_features(x)
    assert features == [1.0, 1.0]
    assert totals == [1.0, 2.0]


def test_doublets_nan_totals_and_zero_max_rate():
    from aria.utils.chromatin_robustness import detect_atac_doublets

    rng = np.random.default_rng(0)
    X = rng.poisson(0.3, size=(100, 200)).astype(float)
    X[:5] = rng.poisson(4.0, size=(5, 200))
    X[50, 0] = np.nan                               # one poisoned cell
    res = detect_atac_doublets(X, min_cells=50)
    assert res["ran"] and res["n_doublets"] >= 1
    assert all(math.isfinite(s) for s in res["_doublet_score"])

    res0 = detect_atac_doublets(X, min_cells=50, max_rate=0.0)
    assert res0["n_doublets"] == 0


def test_consensus_provenance_zero_reproducibility_is_not_missing():
    from aria.utils.chromatin_robustness import assess_consensus_peak_provenance

    res = assess_consensus_peak_provenance(
        ["chr1:1-100", "chr1:200-300"],
        metadata={"method": "consensus", "n_samples": 3, "overlap_fraction": 0.0,
                  "rare_peak_policy": "keep"})
    assert res["metadata"]["reproducibility"] == 0.0
    assert res["status"] == "verified"
