"""Wave 6C regression tests: Domain 2 (Chromatin + Epigenomics + TOBIAS) polish.

Covers all verified HIGH & MEDIUM remediation items in Wave 6C:
  1. D2-H6, D2-M2b, D2-M6, D2-M8c (chromatin_agent.py)
  2. D2-M1 (chromatin_qc.py)
  3. D2-M2a (chromatin_peaks.py)
  4. D2-M3 & CDX-07a (chromatin_peak_ora.py, chromatin_peak_annotation.py)
  5. D2-M4 & CDX-07b (chromatin_bulk_diffacc.py, epigenomic_spikein_de.py)
  6. D2-M5 (chromatin_robustness.py, chromatin_helpers.py)
  7. D2-M6b, CLD-05a/b/c, D2-M8a/b (chromatin_lsi_clustering.py,
     chromatin_fragments_to_matrix.py, atac_align.py, chromatin_regulatory.py)
  8. D2-M7 (narrators/chromatin.py, _narrative_chromatin.py)
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


def _bare_agent(env=None):
    from aria.agents.chromatin_agent import ChromatinAgent

    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = env
    agent.llm = None
    agent.published = []
    agent.publish_status = lambda *a, **k: None
    agent.publish_finding = lambda *a, **k: agent.published.append((a, k))
    agent.publish_escalation = lambda *a, **k: None
    agent._publish_qc_finding = lambda *a, **k: None
    agent._publish_peaks_finding = lambda *a, **k: None
    return agent


# ── 1. D2-H6, D2-M2b, D2-M6, D2-M8c (chromatin_agent.py) ─────────────────────

def test_run_motif_enrichment_returns_early_when_skipped(monkeypatch):
    import aria.agents.debate_council as dc_mod

    debate_called = []

    class _BoomCouncil:
        def __init__(self, *a, **k):
            debate_called.append(True)
            raise AssertionError("DebateCouncil must not be called when motifs skipped")

    monkeypatch.setattr(dc_mod, "DebateCouncil", _BoomCouncil)

    class _Env:
        def run_in_stack(self, **kw):
            return {"status": "success", "ran": False, "reason": "no_genome_fasta"}

    agent = _bare_agent(_Env())
    res = agent._run_motif_enrichment(
        "exp1", {"genome": "hg38"}, {}, {"peaks_path": "/tmp/p.bed"}, "scATAC"
    )
    assert res["ran"] is False
    assert debate_called == []
    assert agent.published == []


def test_select_peak_universe_skips_empty_consensus_in_bulk_and_epigenomic(tmp_path):
    from aria.agents.chromatin_agent import _select_peak_universe

    empty_cons = tmp_path / "empty_consensus.bed"
    empty_cons.write_text("", encoding="utf-8")
    nonempty_cons = tmp_path / "nonempty_consensus.bed"
    nonempty_cons.write_text("chr1\t10\t50\n", encoding="utf-8")
    pooled = tmp_path / "pooled.narrowPeak"
    pooled.write_text("chr1\t10\t50\n", encoding="utf-8")

    # Explicit n_consensus_peaks=0 falls back to peaks_path even if file exists
    assert _select_peak_universe({
        "consensus_peaks_path": str(nonempty_cons),
        "n_consensus_peaks": 0,
        "peaks_path": str(pooled),
    }) == str(pooled)

    # Absent n_consensus_peaks checks file size > 0
    assert _select_peak_universe({
        "consensus_peaks_path": str(empty_cons),
        "peaks_path": str(pooled),
    }) == str(pooled)
    assert _select_peak_universe({
        "consensus_peaks_path": str(nonempty_cons),
        "peaks_path": str(pooled),
    }) == str(nonempty_cons)

    # Explicit n_consensus_peaks > 0 selects consensus even on mock path
    assert _select_peak_universe({
        "consensus_peaks_path": "/virtual/consensus.bed",
        "n_consensus_peaks": 42,
        "peaks_path": str(pooled),
    }) == "/virtual/consensus.bed"

    # Wired into _run_bulk_atac and _run_epigenomic_bulk
    bam = tmp_path / "s1.bam"
    bam.write_bytes(b"BAM")
    seen_peak_paths = []

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            name = Path(script_path).name
            if name == "chromatin_qc.py":
                return {"status": "success", "data_type": params["data_type"], "n_samples": 1}
            if name == "chromatin_peaks.py":
                return {
                    "status": "success",
                    "data_type": params["data_type"],
                    "n_peaks": 10,
                    "peaks_path": str(pooled),
                    "consensus_peaks_path": str(empty_cons),
                    "n_consensus_peaks": 0,
                }
            if name == "chromatin_peak_counts.py":
                seen_peak_paths.append((params["data_type"], params["peaks_path"]))
                return {"status": "error", "error_type": "Stop"}
            return {"status": "skipped", "ran": False}

    agent = _bare_agent(_Env())
    ctx = {"genome": "hg38", "comparisons": [["A", "B"]]}
    agent._run_bulk_atac("exp_atac", ctx, {}, [str(bam)])
    agent._run_epigenomic_bulk("exp_chip", ctx, {}, [str(bam)], assay_type="ChIP")
    assert seen_peak_paths == [
        ("bulk_ATAC", str(pooled)),
        ("ChIP", str(pooled)),
    ]


def test_resolve_exp_gtf_consistent_across_annotation_ora_and_regulatory(tmp_path):
    from aria.agents.chromatin_agent import _resolve_exp_gtf

    assert _resolve_exp_gtf({"gtf_path": "/a.gtf", "gtf": "/b.gtf"}) == "/a.gtf"
    assert _resolve_exp_gtf({"gene_annotation_gtf": "/g.gtf"}) == "/g.gtf"
    assert _resolve_exp_gtf({"gtf": "/c.gtf"}) == "/c.gtf"
    assert _resolve_exp_gtf({"gtf_file": "/d.gtf"}) == "/d.gtf"
    assert _resolve_exp_gtf({"genome_config": {"gtf": "/cfg.gtf"}}) == "/cfg.gtf"
    assert _resolve_exp_gtf({}) is None

    seen_gtf = {}

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            name = Path(script_path).name
            if name == "chromatin_peak_annotation.py":
                seen_gtf["annotation"] = params.get("gtf")
                return {"status": "success", "ran": True}
            if name == "chromatin_peak_ora.py":
                seen_gtf["ora"] = params.get("gtf")
                return {"status": "success", "ran": True}
            if name == "chromatin_regulatory.py":
                seen_gtf["regulatory"] = params.get("gtf_path")
                return {"status": "success", "ran": True}
            if name == "chromatin_qc.py":
                return {"status": "success", "data_type": "scATAC"}
            if name == "chromatin_lsi_clustering.py":
                return {"status": "success", "output_path": str(tmp_path / "lsi.h5ad")}
            if name == "chromatin_diffacc.py":
                return {"status": "success", "per_cluster": {}, "pseudobulk": {}}
            if name == "chromatin_motifs.py":
                return {"status": "success", "ran": False}
            return {"status": "skipped", "ran": False}

    agent = _bare_agent(_Env())
    ctx = {
        "genome": "hg38",
        "genome_fasta": str(tmp_path / "ref.fa"),
        "genome_config": {"gtf": "/shared/genome_cfg.gtf"},
        "enable_footprinting": False,
    }
    da_res = {"comparisons": [{"test": "A", "reference": "B", "full_results_csv": "/da.csv"}]}
    agent._run_bulk_atac_annotation("e", ctx, da_res, [str(tmp_path / "a.bam")])
    agent._run_bulk_atac_ora("e", ctx, da_res, [str(tmp_path / "a.bam")])
    agent._run_scatac_matrix("e", ctx, {}, str(tmp_path / "in.h5ad"))

    assert seen_gtf == {
        "annotation": "/shared/genome_cfg.gtf",
        "ora": "/shared/genome_cfg.gtf",
        "regulatory": "/shared/genome_cfg.gtf",
    }


def test_run_scatac_cohort_matrices_forwards_structured_peak_provenance(tmp_path):
    captured_lsi = {}

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            name = Path(script_path).name
            if name == "chromatin_merge_matrices.py":
                return {
                    "status": "success",
                    "output_path": str(tmp_path / "merged.h5ad"),
                    "peak_provenance": "cohort_interval_union_consensus",
                    "n_consensus_peaks": 500,
                }
            if name == "chromatin_qc.py":
                return {"status": "success", "data_type": "scATAC"}
            if name == "chromatin_lsi_clustering.py":
                captured_lsi.update(params)
                return {"status": "error", "error_type": "Stop"}
            return {"status": "skipped"}

    agent = _bare_agent(_Env())
    f1, f2 = str(tmp_path / "s1.h5ad"), str(tmp_path / "s2.h5ad")
    agent._run_scatac_cohort_matrices("e", {"genome": "hg38"}, {}, [f1, f2])
    prov = captured_lsi["peak_provenance"]
    assert isinstance(prov, dict)
    assert prov["method"] == "overlap_unified"
    assert prov["source"] == "cohort_interval_union_consensus"
    assert prov["n_samples"] == 2
    assert prov["n_consensus_peaks"] == 500


# ── 2. D2-M1 (chromatin_qc.py) ───────────────────────────────────────────────

def test_chromatin_qc_merges_intervals_and_handles_multi_mb_broad_peaks():
    from aria.scripts.chromatin_qc import _frip_from_intervals, _norm_peak_keys

    # A 3 Mb broad domain containing a nested shorter peak at 1.5 Mb
    peaks = {"chr1": [(100, 3_100_000), (1_500_000, 1_500_500), (3_100_000, 3_100_500)]}
    norm = _norm_peak_keys(peaks)
    assert norm == {"1": [(100, 3_100_500)]}

    # Fragment at 2.8 Mb sits >1 Mb after both starts — previously missed by 1 Mb scan cap!
    frags = [
        ("chr1", 2_800_000, 2_800_150),  # inside 3 Mb broad peak
        ("1", 3_100_200, 3_100_300),     # inside merged adjacent peak
        ("chr1", 4_000_000, 4_000_100),  # outside
    ]
    frip, n = _frip_from_intervals(frags, peaks)
    assert n == 3
    assert frip == pytest.approx(2 / 3)


def test_chromatin_qc_discloses_frip_scan_truncation(tmp_path, monkeypatch):
    from aria.scripts import chromatin_qc

    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t100\t200\n", encoding="utf-8")
    frag = tmp_path / "fragments.tsv"
    lines = ["# comment header\n"]
    # First 5 fragments overlap peak; next 5 do not
    for i in range(5):
        lines.append(f"chr1\t120\t160\tBC1\t1\n")
    for i in range(5):
        lines.append(f"chr1\t500\t600\tBC1\t1\n")
    frag.write_text("".join(lines), encoding="utf-8")

    trunc = chromatin_qc._estimate_frip_stats(str(frag), str(bed), scan_limit=5)
    assert trunc["frip"] == pytest.approx(1.0)
    assert trunc["frip_scan_truncated"] is True
    assert trunc["frip_fragments_scanned"] == 5

    full = chromatin_qc._estimate_frip_stats(str(frag), str(bed), scan_limit=10)
    assert full["frip"] == pytest.approx(0.5)
    assert full["frip_scan_truncated"] is False
    assert full["frip_fragments_scanned"] == 10

    # Verify _scatac_qc propagates fields + warning
    fake_mu = types.ModuleType("muon")
    fake_mu.atac = SimpleNamespace(tl=SimpleNamespace(locate_fragments=lambda *a, **k: None))
    monkeypatch.setitem(sys.modules, "muon", fake_mu)
    monkeypatch.setattr(
        chromatin_qc,
        "_tss_qc",
        lambda *a, **k: {"median": 8.5, "reason": None, "tsse": [8.5], "log10_depth": [3.0]},
    )
    qc_res = chromatin_qc._scatac_qc(
        [str(frag)], "hg38", "Homo sapiens", [], peaks_bed=str(bed), frip_scan_limit=5
    )
    assert qc_res["frip_scan_truncated"] is True
    assert qc_res["frip_fragments_scanned"] == 5
    assert any("FRiP fragment scan capped" in w for w in qc_res["warnings"])


# ── 3. D2-M2a (chromatin_peaks.py) ───────────────────────────────────────────

def test_chromatin_peaks_dedups_bam_stems_and_clears_empty_post_blacklist_consensus(
    tmp_path, monkeypatch
):
    from aria.scripts import chromatin_peaks as pk
    from aria.utils import seacr_idr

    d1, d2 = tmp_path / "rep1", tmp_path / "rep2"
    d1.mkdir()
    d2.mkdir()
    bam1 = d1 / "aligned.bam"
    bam2 = d2 / "aligned.bam"
    bam1.write_bytes(b"b1")
    bam2.write_bytes(b"b2")

    # Blacklist wipes out the ONLY overlapping region (chr1:100-200) while leaving pooled chr1:500-600
    bl = tmp_path / "blacklist.bed"
    bl.write_text("chr1\t90\t250\n", encoding="utf-8")

    def fake_macs3(cmd, timeout=7200):
        name = cmd[cmd.index("-n") + 1]
        outdir = Path(cmd[cmd.index("--outdir") + 1])
        outdir.mkdir(parents=True, exist_ok=True)
        if name.startswith("rep1_"):
            text = "chr1\t100\t200\nchr1\t500\t600\n"
        elif name.startswith("rep2_"):
            text = "chr1\t100\t200\nchr1\t900\t1000\n"
        else:
            text = "chr1\t100\t200\nchr1\t500\t600\n"
        (outdir / f"{name}_peaks.narrowPeak").write_text(text, encoding="utf-8")
        return {"status": "success"}

    frip_vals = {str(bam1): 0.25, str(bam2): 0.45}
    monkeypatch.setattr(pk, "_run_macs3", fake_macs3)
    monkeypatch.setattr(pk, "_compute_frip", lambda p, *a, **k: frip_vals[str(p)])
    # Also simulate post-blacklist consensus wipe by letting rep calls keep chr1:100-200 before consensus filter
    orig_filter = pk._filter_blacklist
    calls = {"n": 0}

    def selective_filter(peaks_file, bl_bed):
        calls["n"] += 1
        # Only apply blacklist to pooled and consensus, so overlap consensus is created then emptied
        if "replicate_peaks" in str(peaks_file):
            return 0
        return orig_filter(peaks_file, bl_bed)

    monkeypatch.setattr(pk, "_filter_blacklist", selective_filter)
    monkeypatch.setattr(
        seacr_idr,
        "write_normalized_bigwig",
        lambda bam, out_bw, **kw: (
            Path(out_bw).parent.mkdir(parents=True, exist_ok=True),
            Path(out_bw).write_bytes(b"BW"),
            {"bigwig_path": str(out_bw)},
        )[2],
    )

    res = pk.chromatin_peaks({
        "data_type": "bulk_ATAC",
        "files": [str(bam1), str(bam2)],
        "genome": "hg38",
        "blacklist_bed": str(bl),
        "write_bigwig": True,
        "output_dir": str(tmp_path / "peaks"),
    })
    assert res["status"] == "success", res
    assert res["frip_per_sample"] == {"aligned": 0.25, "aligned_2": 0.45}
    assert res["frip"] == pytest.approx(0.35)
    assert set(res["bigwig_tracks"]) == {"aligned", "aligned_2"}
    assert res["bigwig_tracks"]["aligned"] != res["bigwig_tracks"]["aligned_2"]
    assert res["consensus_peaks_path"] is None
    assert res["n_consensus_peaks"] is None
    assert res["peak_reproducibility"]["consensus_peaks_path"] is None
    assert any("0 regions remaining" in w for w in res["warnings"])


# ── 4. D2-M3 & CDX-07a (peak_ora & peak_annotation) ──────────────────────────

def test_peak_ora_sanitizes_comp_key_in_output_paths(tmp_path, monkeypatch):
    import aria.scripts.rna_bulk.ora as ora_mod
    import aria.scripts.rna_pathway_viz as viz_mod
    from aria.scripts.chromatin_peak_ora import chromatin_peak_ora

    gtf = tmp_path / "ann.gtf"
    gtf.write_text(
        'chr1\tt\tgene\t1000\t5000\t.\t+\t.\tgene_id "G1"; gene_name "GENE1";\n'
        'chr1\tt\texon\t1000\t1200\t.\t+\t.\tgene_id "G1"; gene_name "GENE1";\n',
        encoding="utf-8",
    )
    da_csv = tmp_path / "da.csv"
    pd.DataFrame({
        "peak": ["chr1:1050-1150"],
        "log2FoldChange": [2.0],
        "padj": [0.001],
        "significant": [True],
    }).to_csv(da_csv, index=False)

    monkeypatch.setattr(
        ora_mod,
        "_run_pathway_enrichment",
        lambda *a, **k: (
            {"GO/BP": [{"term": "t1", "padj": 1e-4}]},
            [],
            {"method": "local_hypergeometric"},
        ),
    )

    def fake_dotplot(terms, db, title, out_path):
        p = Path(out_path)
        assert p.parent == tmp_path / "ora_out"
        p.write_bytes(b"PNG")
        p.with_suffix(".svg").write_text("<svg/>", encoding="utf-8")
        return str(p)

    monkeypatch.setattr(viz_mod, "make_ora_dotplot", fake_dotplot)
    res = chromatin_peak_ora({
        "data_type": "ChIP",
        "gtf": str(gtf),
        "genome": "hg38",
        "comparisons": [{
            "test": "T-cell/high",
            "reference": "ctrl/base",
            "full_results_csv": str(da_csv),
        }],
        "output_dir": str(tmp_path / "ora_out"),
    })
    assert res["status"] == "success"
    assert res["validation_level"] == "production"
    expected_key = "chip_peak_ora_GO_BP_T-cell_high_vs_ctrl_base"
    assert expected_key in res["figures"]
    assert expected_key in res["figures_svg"]
    assert Path(res["figures"][expected_key]).parent == tmp_path / "ora_out"


@pytest.mark.parametrize("dt,expected_level,expected_prefix", [
    ("bulk_ATAC", "production", "bulk_atac"),
    ("scATAC", "production", "scatac"),
    ("ChIP", "production", "chip"),
    ("CUT_AND_RUN", "production", "cut_and_run"),
    ("CUT_AND_TAG", "production", "cut_and_tag"),
])
def test_peak_annotation_modality_stamp_production_epigenomics(
    dt, expected_level, expected_prefix
):
    from aria.scripts.chromatin_peak_annotation import _modality_stamp

    got_dt, got_level, got_prefix = _modality_stamp({"data_type": dt})
    assert (got_dt, got_level, got_prefix) == (dt, expected_level, expected_prefix)


def test_annotate_peaks_prefers_overlapping_long_gene_over_closer_nonoverlapping_tss(tmp_path):
    from aria.scripts.chromatin_peak_annotation import (
        annotate_peaks,
        parse_gtf_gene_models,
    )

    lines = [
        # Long gene spanning 1,000 .. 2,000,000 with an exon at 1,500,000..1,500,500
        'chr1\tt\tgene\t1000\t2000000\t.\t+\t.\tgene_id "LONG"; gene_name "LONG_GENE";',
        'chr1\tt\texon\t1000\t1200\t.\t+\t.\tgene_id "LONG"; gene_name "LONG_GENE";',
        'chr1\tt\texon\t1500000\t1500500\t.\t+\t.\tgene_id "LONG"; gene_name "LONG_GENE";',
    ]
    # Insert 35 small non-overlapping genes between 1,000,000 and 1,490,000 so
    # LONG_GENE's TSS (1000) is far outside the +/-25 TSS neighbor window!
    for i in range(35):
        s = 1_000_000 + i * 14_000
        lines.append(
            f'chr1\tt\tgene\t{s}\t{s + 1000}\t.\t+\t.\tgene_id "S{i}"; gene_name "SHORT_{i}";'
        )
        lines.append(
            f'chr1\tt\texon\t{s}\t{s + 1000}\t.\t+\t.\tgene_id "S{i}"; gene_name "SHORT_{i}";'
        )
    # Closest non-overlapping gene just 5 kb downstream (outside 2 kb promoter window)
    lines.append(
        'chr1\tt\tgene\t1505000\t1507000\t.\t+\t.\tgene_id "NEAR"; gene_name "NEAR_SHORT";'
    )
    gtf = tmp_path / "long_gene.gtf"
    gtf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    models = parse_gtf_gene_models(str(gtf))

    res = annotate_peaks(
        [
            "chr1:1500100-1500400",  # midpoint inside LONG_GENE exon
            "chr1:1500450-1500700",  # interval overlaps LONG_GENE exon (midpoint 1500575 outside exon, in body)
            "chr1:1501000-1501200",  # intron of LONG_GENE (not overlapping NEAR_SHORT)
        ],
        models,
    )
    by_pk = {r["peak"]: r for r in res["rows"]}
    assert by_pk["chr1:1500100-1500400"]["nearest_gene"] == "LONG_GENE"
    assert by_pk["chr1:1500100-1500400"]["feature"] == "Exonic"
    assert by_pk["chr1:1500450-1500700"]["nearest_gene"] == "LONG_GENE"
    assert by_pk["chr1:1500450-1500700"]["feature"] == "Exonic"
    assert by_pk["chr1:1501000-1501200"]["nearest_gene"] == "LONG_GENE"
    assert by_pk["chr1:1501000-1501200"]["feature"] == "Intronic"


# ── 5. D2-M4 & CDX-07b (chromatin_bulk_diffacc & epigenomic_spikein_de) ──────

@pytest.mark.parametrize("opt_out_mode", [
    "deseq2_median_of_ratios",
    "median_of_ratios",
    "conventional",
])
def test_bulk_diffacc_opts_out_of_metadata_spikein_when_conventional_requested(
    tmp_path, monkeypatch, opt_out_mode
):
    from aria.scripts import chromatin_bulk_diffacc as mod
    from aria.scripts import rna_bulk_de

    samples = [f"s{i}" for i in range(6)]
    counts = pd.DataFrame(
        [[20, 22, 24, 90, 92, 94], [50, 52, 48, 51, 49, 50]],
        index=["chr1:100-200", "chr1:400-500"],
        columns=samples,
    )
    counts.index.name = "peak_id"
    cpath = tmp_path / "counts.tsv"
    counts.to_csv(cpath, sep="\t")
    # Metadata includes a spikein_reads column (even with a bad value!) — opting out must skip it.
    meta = pd.DataFrame({
        "sample_id": samples,
        "condition": ["ctrl"] * 3 + ["trt"] * 3,
        "replicate": ["r1", "r2", "r3", "r1", "r2", "r3"],
        "spikein_reads": [10_000, -5, 10_000, 10_000, None, 10_000],
    })
    mpath = tmp_path / "meta.tsv"
    meta.to_csv(mpath, sep="\t", index=False)

    used_conventional = []

    def fake_deseq2(c, m, *a, **k):
        used_conventional.append(True)
        df = pd.DataFrame(
            {"log2FoldChange": [2.0, 0.0], "padj": [0.001, 0.9], "baseMean": [50.0, 50.0]},
            index=c.index,
        )
        return {"status": "success", "results": df, "n_replicates": 6}, []

    monkeypatch.setattr(rna_bulk_de, "_run_deseq2", fake_deseq2)
    monkeypatch.setattr(mod, "compute_counts_lobo_jackknife", lambda *a, **k: {"status": "skipped"})

    res = mod.chromatin_bulk_diffacc({
        "data_type": "ChIP",
        "counts_matrix_path": str(cpath),
        "sample_metadata_path": str(mpath),
        "comparisons": [{"test": "trt", "reference": "ctrl"}],
        "normalization_mode": opt_out_mode,
        "output_dir": str(tmp_path / f"out_{opt_out_mode}"),
    })
    assert res["status"] == "success" and res["ran"] is True, res
    assert res["normalization_mode"] == "deseq2_median_of_ratios"
    assert res["spikein_normalization"]["used"] is False
    assert used_conventional == [True]


@pytest.mark.parametrize("kw_name", ["spikein_counts", "spikein_scale_factors"])
def test_spikein_all_none_dict_fails_closed(tmp_path, kw_name):
    from aria.scripts.chromatin_bulk_diffacc import chromatin_bulk_diffacc
    from aria.utils.epigenomic_spikein_de import resolve_replicate_spikein_size_factors

    samples = [f"s{i}" for i in range(6)]
    raw_meta = pd.DataFrame(
        {"condition": ["ctrl"] * 3 + ["trt"] * 3, "replicate": ["r1", "r2", "r3"] * 2},
        index=samples,
    )
    rep_meta = pd.DataFrame(
        {"condition": ["ctrl"] * 3 + ["trt"] * 3},
        index=[f"ctrl::r{i}" for i in (1, 2, 3)] + [f"trt::r{i}" for i in (1, 2, 3)],
    )
    all_none = {s: None for s in samples}
    sf, info = resolve_replicate_spikein_size_factors(
        raw_meta, rep_meta, "condition", "replicate", **{kw_name: all_none}
    )
    assert sf is None
    assert info["used"] is False
    assert info["explicit_requested"] is True
    assert info["invalid_reason"]

    # End-to-end chromatin_bulk_diffacc fails closed (even without normalization_mode set)
    cpath = tmp_path / "counts.tsv"
    pd.DataFrame(
        [[10, 12, 11, 50, 52, 48]], index=["chr1:1-100"], columns=samples
    ).rename_axis("peak_id").to_csv(cpath, sep="\t")
    mpath = tmp_path / "meta.tsv"
    raw_meta.rename_axis("sample_id").reset_index().to_csv(mpath, sep="\t", index=False)

    res = chromatin_bulk_diffacc({
        "data_type": "CUT_AND_RUN",
        "counts_matrix_path": str(cpath),
        "sample_metadata_path": str(mpath),
        "comparisons": [{"test": "trt", "reference": "ctrl"}],
        kw_name: all_none,
        "output_dir": str(tmp_path / f"da_{kw_name}"),
    })
    assert res["ran"] is False
    assert res["normalization_mode"] is None
    assert res["spikein_normalization"]["explicit_requested"] is True
    assert res["spikein_normalization"]["invalid_reason"]


# ── 6. D2-M5 (chromatin_robustness.py & chromatin_helpers.py) ────────────────

def test_chromatin_robustness_median_filters_nonfinite():
    from aria.utils.chromatin_robustness import _median, detect_atac_doublets

    assert _median([1.0, float("nan"), 3.0, float("inf"), None, 5.0]) == 3.0
    assert _median([float("nan"), float("inf"), None]) is None

    rng = np.random.default_rng(7)
    X = rng.poisson(2.0, size=(60, 20)).astype(float)
    X[0, 0] = float("nan")
    X[1, 0] = float("inf")
    res = detect_atac_doublets(X, min_cells=50)
    assert res["ran"] is True
    assert math.isfinite(res["metrics"]["median_fragments"])
    json.dumps({k: v for k, v in res.items() if not k.startswith("_")}, allow_nan=False)


def test_bulk_da_motif_regions_filters_nonfinite_lfc_and_padj(tmp_path):
    from aria.agents.chromatin_helpers import bulk_da_motif_regions

    csv_path = tmp_path / "da.csv"
    csv_path.write_text(
        "peak,log2FoldChange,padj,significant\n"
        "chr1:1-10,inf,1e-9,True\n"
        "chr1:20-30,-inf,1e-9,True\n"
        "chr1:40-50,nan,1e-9,True\n"
        "chr1:60-70,2.5,nan,True\n"
        "chr1:80-90,1.5,1e-4,True\n"
        "chr1:100-110,-2.0,inf,True\n"
        "chr1:120-130,-1.2,1e-3,True\n",
        encoding="utf-8",
    )
    regions, bg, warnings = bulk_da_motif_regions([{
        "test": "T",
        "reference": "C",
        "status": "success",
        "full_results_csv": str(csv_path),
    }])
    assert warnings == []
    # Non-finite LFC peaks excluded from foreground; finite padj ranks ahead of nan/inf padj
    assert regions["T_vs_C::up_in_T"] == ["chr1:80-90", "chr1:60-70"]
    assert regions["T_vs_C::up_in_C"] == ["chr1:120-130", "chr1:100-110"]
    assert len(bg) == 7


# ── 7. D2-M6b, CLD-05, D2-M8a/b ─────────────────────────────────────────────

def test_lsi_normalizes_string_provenance_and_runs_per_library_doublets(tmp_path):
    pytest.importorskip("scanpy")
    import anndata as ad
    from scipy import sparse
    from aria.scripts.chromatin_lsi_clustering import chromatin_lsi_clustering

    rng = np.random.default_rng(123)
    n_peaks = 120
    # Library A (100 cells, shallow: ~15 frags) + Library B (100 cells, deep: ~180 frags)
    # Plus 3 true doublets inside Library A and 3 inside Library B; plus tiny Library C (15 cells)
    X_a = rng.poisson(0.12, size=(100, n_peaks)).astype(float)
    X_a[:3] = rng.poisson(1.5, size=(3, n_peaks)).astype(float)  # A doublets

    X_b = rng.poisson(1.5, size=(100, n_peaks)).astype(float)
    X_b[:3] = rng.poisson(12.0, size=(3, n_peaks)).astype(float)  # B doublets

    X_c = rng.poisson(0.4, size=(15, n_peaks)).astype(float)  # tiny library (<50)

    X = np.vstack([X_a, X_b, X_c])
    adata = ad.AnnData(sparse.csr_matrix(X))
    adata.obs_names = [f"cell_{i}" for i in range(X.shape[0])]
    adata.var_names = [f"chr1:{i * 1000}-{i * 1000 + 400}" for i in range(n_peaks)]
    adata.obs["library_id"] = ["libA"] * 100 + ["libB"] * 100 + ["libC"] * 15
    adata.uns["peak_provenance"] = "cohort_interval_union_consensus"

    in_path = tmp_path / "cohort.h5ad"
    adata.write_h5ad(in_path)

    res = chromatin_lsi_clustering({
        "data_path": str(in_path),
        "output_dir": str(tmp_path / "lsi_out"),
        "n_components": 8,
        "doublet_min_cells": 50,
    })
    assert res["status"] == "success", res
    # String provenance normalized into structured dict -> partial (not unverified!)
    cp = res["consensus_peaks"]
    assert cp["status"] == "partial"
    assert cp["metadata"]["method"] == "overlap_unified"
    assert cp["metadata"]["n_samples"] == 3

    # Per-library doublet screen: deep libB must NOT be over-trimmed
    dbl = res["doublets"]
    assert dbl["scope"] == "per_library"
    assert dbl["groupby"] == "library_id"
    assert dbl["per_library"]["libA"]["scope"] == "per_library"
    assert dbl["per_library"]["libB"]["scope"] == "per_library"
    assert dbl["per_library"]["libC"]["scope"] == "global_fallback"
    assert dbl["per_library"]["libB"]["n_doublets"] <= 6  # true 3, never over-trimmed!
    assert dbl["per_library"]["libA"]["n_doublets"] >= 2


def test_fragments_to_matrix_synthesizes_manifest_when_explicit_none_or_empty(
    tmp_path, monkeypatch
):
    from aria.scripts import chromatin_fragments_to_matrix as br

    f1 = tmp_path / "s1.fragments.tsv"
    f2 = tmp_path / "s2.fragments.tsv"
    f1.write_text("chrA\t10\t20\tC1\t1\n", encoding="utf-8")
    f2.write_text("chrA\t30\t40\tC2\t1\n", encoding="utf-8")

    captured = {}

    class _Stop(Exception):
        pass

    def fake_stage_is_current(manifest_path, inputs, params, tool_version):
        captured["library_manifest"] = params["library_manifest"]
        raise _Stop()

    monkeypatch.setitem(sys.modules, "snapatac2", types.SimpleNamespace(__version__="2.0"))
    import aria.utils.stage_manifest as sm
    monkeypatch.setattr(sm, "stage_is_current", fake_stage_is_current)
    monkeypatch.setattr(br, "_matrix_valid", lambda *a, **k: True)

    for empty_manifest in (None, [], {"libraries": []}):
        captured.clear()
        with pytest.raises(_Stop):
            br.chromatin_fragments_to_matrix({
                "fragments_files": [str(f1), str(f2)],
                "sample_ids": ["s1", "s2"],
                "library_manifest": empty_manifest,
                "chrom_sizes": {"chrA": 1000},
                "output_dir": str(tmp_path / "out"),
            })
        assert len(captured["library_manifest"]) == 2
        assert [m["library_id"] for m in captured["library_manifest"]] == ["s1", "s2"]


def test_atac_align_chrmt_case_insensitive_full_precision_and_filtered_spikein(
    monkeypatch, tmp_path
):
    from aria.scripts import atac_align as aa

    assert "chrMT" in aa.DEFAULT_MITO_CONTIGS
    lines = [
        "@SQ\tSN:chr1\tLN:1000\n",
        "r1\t0\tchrMT\t1\t60\n",
        "r2\t0\tChRm\t1\t60\n",
        "r3\t0\tmt\t1\t60\n",
        "r4\t0\tchr1\t1\t60\n",
    ]
    assert list(aa._strip_mito_sam(lines)) == [lines[0], lines[4]]

    # _bam_metrics uses samtools view -c -q {mapq} -F {excl_flags} on spike-in contigs
    # and preserves full float precision (no underflow to 0.0)
    seen_cmds = []

    def fake_run(cmd, **kw):
        seen_cmds.append(cmd)
        if cmd[:2] == ["samtools", "flagstat"]:
            return subprocess.CompletedProcess(
                cmd, 0, stdout="100 + 0 in total\n95 + 0 mapped (95.0%)\n5 + 0 duplicates\n", stderr=""
            )
        if cmd[:2] == ["samtools", "idxstats"]:
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout="chr1\t1000\t900\t0\nchrMT\t100\t100\t0\nspike_chr1\t500\t20000000\t0\n*\t0\t0\t0\n",
                stderr="",
            )
        if cmd[:3] == ["samtools", "view", "-c"]:
            # Return filtered spike-in count (10M < raw 20M)
            return subprocess.CompletedProcess(cmd, 0, stdout="10000000\n", stderr="")
        return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    m = aa._bam_metrics(
        tmp_path / "markdup.bam",
        aa.DEFAULT_MITO_CONTIGS,
        [],
        spikein_contig_prefix="spike_",
        spikein_scale_constant=1.0,
        mapq=20,
        excl_flags=aa.SAMTOOLS_EXCLUDE_FLAGS_KEEP_DUP,
        paired=True,
    )
    assert m["pct_mito"] == pytest.approx(10.0)  # chrMT counted as mito
    assert m["spikein_reads"] == 10_000_000
    assert m["spikein_reads_filtered"] is True
    # Full float precision: 1.0 / 10_000_000 = 1e-7 (old round(..., 6) gave 0.0!)
    assert m["spikein_scale_factor"] == pytest.approx(1e-7)
    assert m["spikein_scale_factor"] > 0.0
    view_cmd = next(c for c in seen_cmds if c[:3] == ["samtools", "view", "-c"])
    assert "-q" in view_cmd and view_cmd[view_cmd.index("-q") + 1] == "20"
    assert "-F" in view_cmd and view_cmd[view_cmd.index("-F") + 1] == str(aa.SAMTOOLS_EXCLUDE_FLAGS_KEEP_DUP)
    assert "-f" in view_cmd and view_cmd[view_cmd.index("-f") + 1] == "2"
    assert "spike_chr1" in view_cmd

    # Clean fallback to idxstats when samtools view fails
    def fake_run_fallback(cmd, **kw):
        if cmd[:3] == ["samtools", "view", "-c"]:
            return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="no index")
        return fake_run(cmd, **kw)

    monkeypatch.setattr(subprocess, "run", fake_run_fallback)
    m_fb = aa._bam_metrics(
        tmp_path / "markdup.bam",
        aa.DEFAULT_MITO_CONTIGS,
        [],
        spikein_contig_prefix="spike_",
        spikein_scale_constant=1.0,
    )
    assert m_fb["spikein_reads"] == 20_000_000
    assert m_fb["spikein_reads_filtered"] is False
    assert m_fb["spikein_scale_factor"] == pytest.approx(5e-8)


def test_chromatin_regulatory_motif_activity_computes_depth_normalized_deviation(tmp_path):
    ad = pytest.importorskip("anndata")
    from aria.scripts.chromatin_regulatory import _motif_activity

    # 4 cells with depths [10, 10, 100, 100]; grand_total = 220
    # Peak 0 (MOTIF_HIGH_DEV): [5, 5, 10, 10] -> obs_frac=[0.5,0.5,0.1,0.1], mean=0.30, bg=30/220=0.136364 -> dev=+0.163636
    # Peak 1 (MOTIF_LOW_DEV):  [1, 1, 50, 50] -> obs_frac=[0.1,0.1,0.5,0.5], mean=0.30, bg=102/220=0.463636 -> dev=-0.163636
    # Peak 2 (bg filler):      [4, 4, 40, 40] -> obs_frac=[0.4,0.4,0.4,0.4], dev=0.0
    X = np.array([
        [5, 1, 4],
        [5, 1, 4],
        [10, 50, 40],
        [10, 50, 40],
    ], dtype=float)
    adata = ad.AnnData(X)
    adata.var_names = ["chr1:1-100", "chr1:200-300", "chr1:400-500"]
    adata.obs["n_fragments"] = X.sum(axis=1)

    mmap = tmp_path / "motifs.csv"
    mmap.write_text(
        "motif_id,peak\n"
        "MOTIF_HIGH_DEV,chr1:1-100\n"
        "MOTIF_LOW_DEV,chr1:200-300\n"
        "MOTIF_FLAT,chr1:400-500\n",
        encoding="utf-8",
    )
    res = _motif_activity(adata, {"motif_peak_map": str(mmap)}, tmp_path)
    assert res["ran"] is True
    top = res["top_motifs"]
    # Ranked by mean_accessibility_deviation (MOTIF_HIGH_DEV > MOTIF_FLAT > MOTIF_LOW_DEV)
    # even though MOTIF_LOW_DEV has a higher raw mean_accessibility_score (25.5 vs 7.5)!
    assert [r["motif_id"] for r in top] == ["MOTIF_HIGH_DEV", "MOTIF_FLAT", "MOTIF_LOW_DEV"]
    assert top[0]["mean_accessibility_deviation"] == pytest.approx(0.30 - 30 / 220, abs=1e-5)
    assert top[0]["mean_accessibility_score"] == pytest.approx(7.5)
    assert top[1]["mean_accessibility_deviation"] == pytest.approx(0.0, abs=1e-6)
    assert top[2]["mean_accessibility_deviation"] == pytest.approx(0.30 - 102 / 220, abs=1e-5)


# ── 8. D2-M7 (narrators/chromatin.py & _narrative_chromatin.py) ──────────────

def test_chromatin_narrator_tables_methods_svg_and_binding_titles(tmp_path, monkeypatch):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from aria.agents import _narrative_chromatin as nc
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.validators import find_causal_language

    captured_figs = []
    real_save = nc._save_dual

    def _cap(fig, png_path):
        captured_figs.append((Path(png_path).name, fig))
        return real_save(fig, png_path)

    monkeypatch.setattr(nc, "_save_dual", _cap)

    comp_csv = tmp_path / "chip_da.csv"
    pd.DataFrame({
        "peak": ["chr1:1-100", "chr1:200-300"],
        "log2FoldChange": [2.0, -1.8],
        "padj": [1e-4, 2e-3],
        "baseMean": [100.0, 80.0],
        "comparison": ["KO_vs_WT", "KO_vs_WT"],
        "significant": [True, True],
    }).to_csv(comp_csv, index=False)

    ora_png = tmp_path / "chip_ora.png"
    ora_svg = tmp_path / "chip_ora.svg"
    ora_png.write_bytes(b"PNG")
    ora_svg.write_text("<svg/>", encoding="utf-8")

    findings = {
        "differential_accessibility": {
            "status": "success",
            "ran": True,
            "data_type": "ChIP",
            "padj_max": 0.05,
            "comparisons": [{
                "test": "KO",
                "reference": "WT",
                "status": "success",
                "full_results_csv": str(comp_csv),
            }],
        },
        "peak_annotation": {
            "status": "success",
            "ran": True,
            "data_type": "ChIP",
            "gtf": "mm39.gtf",
            "promoter_upstream": 2000,
            "promoter_downstream": 2000,
            "feature_distribution_overall": {"Promoter": 15, "Exonic": 5},
            "comparisons": [{
                "test": "KO",
                "reference": "WT",
                "annotation_csv": str(tmp_path / "chip_peak_annotation_KO_vs_WT.csv"),
            }],
        },
        "peak_ora": {
            "status": "success",
            "ran": True,
            "data_type": "ChIP",
            "organism": "Mus musculus",
            "comparisons": [{
                "test": "KO",
                "reference": "WT",
                "ora_method": "local_hypergeometric",
                "pathways_summary": {"GO_BP": 4},
            }],
            "figures": {"chip_peak_ora_GO_BP_KO_vs_WT": str(ora_png)},
            "figures_svg": {"chip_peak_ora_GO_BP_KO_vs_WT": str(ora_svg)},
        },
    }

    nc.generate_figures(findings, None, tmp_path / "figs", env_manager=None)
    assert "da_volcano_KO_vs_WT" in findings["figures"]
    assert "da_volcano_KO_vs_WT" in findings["figures_svg"]
    assert "da_ma_KO_vs_WT" in findings["figures_svg"]
    assert "bulk_atac_peak_annotation" in findings["figures_svg"]
    assert findings["figures_svg"]["chip_peak_ora_GO_BP_KO_vs_WT"] == str(ora_svg)

    # Binding terminology in ChIP volcano, MA, and peak-annotation plots
    titles_by_name = {
        name: (fig.axes[0].get_title(), fig.axes[0].get_xlabel())
        for name, fig in captured_figs
    }
    assert "Differential binding" in titles_by_name["da_volcano_KO_vs_WT.png"][0]
    assert "(binding)" in titles_by_name["da_volcano_KO_vs_WT.png"][1]
    assert "Differential binding" in titles_by_name["da_ma_KO_vs_WT.png"][0]
    assert "mean binding" in titles_by_name["da_ma_KO_vs_WT.png"][1]
    assert "differential-binding peaks" in titles_by_name["bulk_atac_peak_annotation.png"][0]

    # ChromatinNarrator surfaces svg_path, peak_annotation CSVs, and peak_annotation/peak_ora methods
    narr = ChromatinNarrator()
    agent_res = {"status": "done", "findings": {"ChIP": {"status": "done", "findings": findings}}}

    figs = {f["id"]: f for f in narr.figures("chromatin_agent", agent_res)}
    assert figs["da_volcano_KO_vs_WT"]["svg_path"] == findings["figures_svg"]["da_volcano_KO_vs_WT"]
    assert figs["chip_peak_ora_GO_BP_KO_vs_WT"]["svg_path"] == str(ora_svg)

    tables = {t["id"]: t for t in narr.tables("chromatin_agent", agent_res)}
    assert "chromatin.peak_annotation.KO_vs_WT" in tables
    assert tables["chromatin.peak_annotation.KO_vs_WT"]["path"].endswith(
        "chip_peak_annotation_KO_vs_WT.csv"
    )

    methods_blob = " ".join(narr.methods("chromatin_agent", agent_res))
    assert "Genomic peak annotation" in methods_blob
    assert "differential binding peaks" in methods_blob
    assert "Peak-proximal functional enrichment" in methods_blob
    assert "GO_BP" in methods_blob
    assert find_causal_language(methods_blob) is None

    # scATAC pseudobulk full_results_csv + peak_annotation/peak_ora methods
    sc_res = {
        "status": "done",
        "findings": {
            "scATAC": {
                "status": "done",
                "findings": {
                    "differential_accessibility": {
                        "status": "success",
                        "pseudobulk": {
                            "ran": True,
                            "output_csv": "/tmp/pb_sig.csv",
                            "full_results_csv": "/tmp/pb_full.csv",
                        },
                    },
                    "peak_annotation": {
                        "status": "success",
                        "ran": True,
                        "data_type": "scATAC",
                        "gtf": "hg38.gtf",
                        "comparisons": [{
                            "test": "trt",
                            "reference": "ctrl",
                            "annotation_csv": "/tmp/scatac_ann.csv",
                        }],
                    },
                    "peak_ora": {
                        "status": "success",
                        "ran": True,
                        "data_type": "scATAC",
                        "organism": "Homo sapiens",
                        "comparisons": [{
                            "test": "trt",
                            "reference": "ctrl",
                            "ora_method": "local_hypergeometric",
                            "pathways_summary": {"Reactome": 3},
                        }],
                    },
                },
            }
        },
    }
    sc_tables = {t["id"]: t["path"] for t in narr.tables("chromatin_agent", sc_res)}
    assert sc_tables["chromatin.da_pseudobulk"] == "/tmp/pb_sig.csv"
    assert sc_tables["chromatin.da_pseudobulk_full"] == "/tmp/pb_full.csv"
    assert sc_tables["chromatin.peak_annotation.trt_vs_ctrl"] == "/tmp/scatac_ann.csv"

    sc_methods = " ".join(narr.methods("chromatin_agent", sc_res))
    assert "Genomic peak annotation" in sc_methods
    assert "Peak-proximal functional enrichment" in sc_methods
    assert find_causal_language(sc_methods) is None
