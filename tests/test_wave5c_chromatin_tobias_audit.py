"""Wave 5C regression tests: Domain 2 (Chromatin + Epigenomics + TOBIAS).

Pins every verified consensus/critical/regression bug fixed in Wave 5C:
  1. TRI-03 & BI-08: spike-in scale factor preservation (target + control),
     seacr_top_quantile / spikein_scale_constant parameter wiring & aliases,
     assay-specific MAPQ preservation, asymmetric spike-in guards in SEACR &
     call_broad_domains.
  2. TRI-04: _merge_peak_frip_into_qc delegates to chromatin_qc thresholds &
     failures, preserves NRF / prior pass_qc=False, fails bulk ATAC on low FRiP
     without TSS, keeps qc_complete=False on sample failure, populates
     frip_distribution (+ render_frip_distribution fallback to frip_per_sample).
  3. BI-06: strict _parse_bool_flag for TSV/CSV is_control strings + NA float
     guards on spike-in cells.
  4. D2-C3: samtools markdup omits -r so flagstat reports true pct_dup while
     samtools view -F 1804 removes duplicates; _count_spikein_index requires
     bwa-mem2 returncode == 0.
  5. D2-C4: chromap 0-based inclusive bc:start:end[:strand] slicing.
  6. BI-07 & D2-H2: _regions_from_csv directional splitting & positive-only
     per-cluster peaks; dash peak IDs + UCSC/Ensembl contig aliasing; explicit
     motif_file populates motif_source['n_motifs'].
  7. CDX-04 & D2-H7: paired _null_label_controls large-N sampling without OOM;
     paired batch-condition confounding check; full-precision change bucketing;
     TOBIAS filafy regex in _aggregate_plots.
  8. CDX-05: resolved data_path + size/mtime_ns in LSI _cache_params.
  9. CDX-06: drop cells with missing/empty condition or replicate before
     pseudobulk aggregation.
 10. D2-H1, D2-H3, D2-H4: gzipped/corrupt blacklist & contig normalization;
     copula IDR degenerate assignment; SEACR no-excess -> inf (0 peaks, no
     warning), cov > 0.0, _fmt_bedgraph_value, full-precision threshold.
 11. BI-09, BI-10, D2-H5, D2-H8, D2-H9: TOBIAS descriptive prose branching;
     _bulk_atac_methods *_alignment match; skipped bulk DA limitation block;
     spikein_normalization['used'] detection; multi-modality chromatin
     narration; render_motif_dotplot padj_max filter & small legend formatting.
"""

from __future__ import annotations

import gzip
import io
import json
import math
import subprocess
import sys
import types
import warnings
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _bare_agent(env=None):
    from aria.agents.chromatin_agent import ChromatinAgent

    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = env
    agent.publish_status = lambda *a, **k: None
    agent.publish_finding = lambda *a, **k: None
    agent.publish_escalation = lambda *a, **k: None
    agent._publish_qc_finding = lambda *a, **k: None
    agent._publish_peaks_finding = lambda *a, **k: None
    return agent


# ── 1. TRI-03 & BI-08 ────────────────────────────────────────────────────────

def test_agent_preprocessing_preserves_assay_mapq_and_forwards_spikein_constant(tmp_path):
    fq = tmp_path / "s_R1.fastq.gz"
    fq.write_bytes(b"")
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\nACGT\n")
    seen = {}

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            name = Path(script_path).name
            seen[name] = params
            if name == "rna_fastq_qc.py":
                return {
                    "status": "success",
                    "n_samples": 1,
                    "samples": [{"name": "s", "r1_trimmed": str(fq), "paired": False}],
                }
            return {
                "status": "success",
                "aligner": "bwa-mem2",
                "bam_files": [{"name": "s", "status": "success", "bam": str(tmp_path / "s.bam")}],
            }

    agent = _bare_agent(_Env())
    bams, pre = agent._run_bulk_atac_preprocessing(
        "exp",
        [str(fq)],
        {
            "genome": "hg38",
            "genome_fasta": str(fa),
            "spikein_scale_numerator": 50000.0,
        },
        data_type="ChIP",
    )
    assert bams == [str(tmp_path / "s.bam")]
    assert pre["analysis"] == "chip_alignment"
    align_p = seen["atac_align.py"]
    # Explicit mapq omitted -> atac_align uses assay default (mapq=20 for ChIP)
    assert "mapq" not in align_p
    assert align_p["assay_type"] == "ChIP"
    assert align_p["spikein_scale_constant"] == 50000.0

    # Explicit mapq IS forwarded
    seen.clear()
    agent._run_bulk_atac_preprocessing(
        "exp",
        [str(fq)],
        {"genome": "hg38", "genome_fasta": str(fa), "mapq": 25, "spikein_scale_constant": 12345.0},
        data_type="CUT_AND_RUN",
    )
    assert seen["atac_align.py"]["mapq"] == 25
    assert seen["atac_align.py"]["spikein_scale_constant"] == 12345.0


def test_agent_preserves_control_spikein_scale_factors_and_seacr_params(tmp_path):
    t1 = tmp_path / "t1.bam"
    t2 = tmp_path / "t2.bam"
    c1 = tmp_path / "igg.bam"
    for p in (t1, t2, c1):
        p.write_bytes(b"")

    meta_tsv = tmp_path / "meta.tsv"
    meta_tsv.write_text(
        "sample_id\tcondition\treplicate\tis_control\tspikein_scale_factor\tspikein_reads\n"
        "t1\ttreat\tr1\tFalse\t0.5\t20000\n"
        "t2\ttreat\tr2\t0\t0.25\t40000\n"
        "igg\tIgG\tr1\tTrue\t0.4\t25000\n",
        encoding="utf-8",
    )

    seen_peaks = {}

    class _Env:
        def run_in_stack(self, *, stack, script_path, params):
            name = Path(script_path).name
            if name == "chromatin_qc.py":
                return {"status": "success", "data_type": "CUT_AND_RUN", "n_samples": 2}
            if name == "chromatin_peaks.py":
                seen_peaks.update(params)
                return {
                    "status": "success",
                    "data_type": "CUT_AND_RUN",
                    "n_peaks": 10,
                    "peaks_path": str(tmp_path / "p.bed"),
                    "frip": 0.2,
                }
            return {"status": "success"}

    agent = _bare_agent(_Env())
    exp_ctx = {
        "genome": "hg38",
        "sample_metadata_path": str(meta_tsv),
        "seacr_top_fraction": 0.02,
    }
    resolved = agent._resolve_epigenomic_samples(exp_ctx, [str(t1), str(t2), str(c1)])
    assert resolved["target_files"] == [str(t1), str(t2)]
    assert resolved["control_files"] == [str(c1)]
    assert resolved["control_sample_ids"] == ["igg"]
    assert resolved["spikein_scale_factors"] == {"t1": 0.5, "t2": 0.25, "igg": 0.4}
    assert resolved["spikein_counts"] == {"t1": 20000.0, "t2": 40000.0, "igg": 25000.0}

    agent._run_epigenomic_bulk("exp", exp_ctx, {}, [str(t1), str(t2), str(c1)], "CUT_AND_RUN")
    assert seen_peaks["target_scale_factors"] == [0.5, 0.25]
    assert seen_peaks["control_scale_factors"] == [0.4]
    assert seen_peaks["seacr_top_quantile"] == pytest.approx(0.02)


def test_chromatin_peaks_legacy_aliases_enter_stage_cache(monkeypatch, tmp_path):
    from aria.scripts import chromatin_peaks

    bam = tmp_path / "ip.bam"
    bam.write_bytes(b"x")

    class _Stop(Exception):
        pass

    monkeypatch.setattr(
        chromatin_peaks,
        "load_valid_stage_result",
        lambda manifest, expected_inputs, expected_parameters: (_ for _ in ()).throw(
            _Stop(expected_parameters)
        ),
    )
    monkeypatch.setattr(chromatin_peaks, "file_record", lambda f: {"path": str(f)})
    with pytest.raises(_Stop) as exc:
        chromatin_peaks.chromatin_peaks({
            "data_type": "CUT_AND_RUN",
            "files": [str(bam)],
            "peak_caller": "seacr",
            "output_dir": str(tmp_path / "pk"),
            "spikein_scale_factors": [0.75],
            "seacr_top_fraction": 0.025,
        })
    params = exc.value.args[0]
    assert params["target_scale_factors"] == [0.75]
    assert params["seacr_top_quantile"] == pytest.approx(0.025)


def test_atac_align_spikein_alias_and_stage_manifest_invalidation(monkeypatch, tmp_path):
    from aria.scripts import atac_align as aa

    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT\n")
    sp_fa = tmp_path / "spike.fa"
    sp_fa.write_text(">sp1\nACGT\n")
    r1 = tmp_path / "R1.fq"
    r1.write_text("@r\nACGT\n+\nIIII\n")

    seen_constants = []

    def fake_align_sample(**kwargs):
        seen_constants.append(kwargs["spikein_scale_constant"])
        out_bam = Path(kwargs["output_dir"]) / kwargs["sample"]["name"] / f"{kwargs['sample']['name']}.filtered.bam"
        out_bam.parent.mkdir(parents=True, exist_ok=True)
        out_bam.write_bytes(b"bam")
        Path(str(out_bam) + ".bai").write_bytes(b"bai")
        return {"name": kwargs["sample"]["name"], "status": "success", "bam": str(out_bam), "n_final": 5}

    monkeypatch.setattr(aa, "_bwa_index_exists", lambda prefix: True)
    monkeypatch.setattr(aa, "_align_sample", fake_align_sample)

    base = {
        "samples": [{"name": "s1", "r1": str(r1), "paired": False}],
        "genome_fasta": str(fasta),
        "spikein_fasta": str(sp_fa),
        "spikein_scale_numerator": 25000.0,
        "output_dir": str(tmp_path / "out"),
    }
    res1 = aa.atac_align(base)
    assert res1["status"] == "success"
    assert seen_constants == [25000.0]
    assert res1["bam_files"][0]["reused_from_stage_manifest"] is False

    # Same params -> cache hit
    res2 = aa.atac_align(base)
    assert res2["bam_files"][0]["reused_from_stage_manifest"] is True

    # Changing spikein_scale_constant -> cache miss!
    res3 = aa.atac_align({**base, "spikein_scale_constant": 50000.0})
    assert res3["bam_files"][0]["reused_from_stage_manifest"] is False
    assert seen_constants[-1] == 50000.0

    # Changing spikein_fasta -> cache miss!
    sp_fa2 = tmp_path / "spike2.fa"
    sp_fa2.write_text(">sp2\nTGCA\n")
    res4 = aa.atac_align({**base, "spikein_scale_constant": 50000.0, "spikein_fasta": str(sp_fa2)})
    assert res4["bam_files"][0]["reused_from_stage_manifest"] is False


def test_seacr_asymmetric_spikein_keeps_norm_with_warning(tmp_path):
    from aria.utils.seacr_idr import call_seacr_peaks

    rng = np.random.default_rng(7)
    t_lines, c_lines = [], []
    pos = 1000
    for _ in range(100):
        t_lines.append(f"chr1\t{pos}\t{pos + 200}\t{float(rng.uniform(1.0, 2.0)):.3f}")
        c_lines.append(f"chr1\t{pos}\t{pos + 200}\t{float(rng.uniform(1.0, 2.0)):.3f}")
        pos += 1000
    for _ in range(20):
        t_lines.append(f"chr1\t{pos}\t{pos + 400}\t40.0")
        pos += 2000
    t_bg = tmp_path / "t.bedgraph"
    c_bg = tmp_path / "c.bedgraph"
    t_bg.write_text("\n".join(t_lines) + "\n")
    c_bg.write_text("\n".join(c_lines) + "\n")

    # Target scaled down 1000x without control scale factors: if norm="non" were
    # forced, all target blocks would sit below unscaled control and yield 0 peaks!
    res_asym = call_seacr_peaks(
        [str(t_bg)],
        tmp_path / "asym.bed",
        control_inputs=[str(c_bg)],
        target_scale_factors=[0.001],
        work_dir=tmp_path / "w_asym",
    )
    assert res_asym["norm"] == "norm"
    assert res_asym["already_spikein_scaled"] is False
    assert any("only one lane" in w for w in res_asym["warnings"])
    assert res_asym["n_peaks"] >= 18

    # Both lanes scaled -> norm="non" preserved
    res_both = call_seacr_peaks(
        [str(t_bg)],
        tmp_path / "both.bed",
        control_inputs=[str(c_bg)],
        target_scale_factors=[0.001],
        control_scale_factors=[0.001],
        work_dir=tmp_path / "w_both",
    )
    assert res_both["norm"] == "non"
    assert res_both["already_spikein_scaled"] is True
    assert res_both["n_peaks"] >= 18


def test_broad_domains_asymmetric_spikein_falls_back_to_cpm(monkeypatch, tmp_path):
    from aria.utils import seacr_idr

    class _Read:
        is_unmapped = is_secondary = is_supplementary = is_duplicate = False
        reference_name = "chr1"

        def __init__(self, pos):
            self.reference_start = pos

    class _Header:
        references = ("chr1",)
        lengths = (100_000,)

    class _Bam:
        header = _Header()

        def __init__(self, positions):
            self._pos = positions

        def fetch(self, until_eof=True):
            return iter(_Read(p) for p in self._pos)

    # IP enriched in 30kb..70kb; Control uniform
    ip_pos = list(range(0, 100_000, 1000)) + [p for p in range(30_000, 70_000, 200)]
    ctrl_pos = list(range(0, 100_000, 500))

    @contextmanager
    def fake_open(path, ref=None):
        yield _Bam(ip_pos if "ip" in str(path) else ctrl_pos)

    monkeypatch.setattr(seacr_idr, "_open_alignment", fake_open)
    res = seacr_idr.call_broad_domains(
        ["ip.bam"],
        tmp_path / "dom.bed",
        control_bams=["ctrl.bam"],
        bin_size=10_000,
        min_domain_bp=20_000,
        ip_scale_factors=[0.01],  # asymmetric: no control_scale_factors!
    )
    assert res["status"] == "success"
    assert res["normalization"] == "CPM"
    assert any("only one lane" in w for w in res["warnings"])
    assert res["n_peaks"] >= 1


# ── 2. TRI-04 ────────────────────────────────────────────────────────────────

def test_merge_peak_frip_preserves_nrf_failure_and_populates_distribution():
    from aria.agents.chromatin_agent import ChromatinAgent

    # ChIP failed NRF < 0.5 -> merging passing FRiP must NOT upgrade pass_qc to True
    findings = {
        "qc": {
            "status": "success",
            "data_type": "ChIP",
            "n_samples": 2,
            "nrf": 0.20,
            "metrics_not_computed": ["frip"],
            "qc_complete": False,
            "pass_qc": False,
            "qc_failures": ["NRF 0.200 < 0.5 (low library complexity)"],
            "warnings": ["NRF 0.200 < 0.5 (low library complexity)"],
        }
    }
    ChromatinAgent._merge_peak_frip_into_qc(
        findings,
        {"data_type": "ChIP", "frip": 0.12, "frip_per_sample": {"s1": 0.10, "s2": 0.14}},
    )
    qc = findings["qc"]
    assert qc["pass_qc"] is False
    assert any("NRF" in f for f in qc["qc_failures"])
    assert qc["frip_distribution"] == [0.10, 0.14]
    assert qc["qc_complete"] is True

    # Bulk ATAC low FRiP fails pass_qc even when TSS is missing
    atac_findings = {
        "qc": {
            "status": "success",
            "data_type": "bulk_ATAC",
            "n_samples": 1,
            "tss_enrichment": None,
            "metrics_not_computed": ["frip", "tss_enrichment"],
            "qc_complete": False,
            "pass_qc": None,
            "warnings": [],
        }
    }
    ChromatinAgent._merge_peak_frip_into_qc(
        atac_findings,
        {"data_type": "bulk_ATAC", "frip": 0.05, "frip_per_sample": {"s1": 0.05}},
    )
    assert atac_findings["qc"]["pass_qc"] is False
    assert atac_findings["qc"]["qc_complete"] is False

    # CUT&RUN uses 0.05 threshold (not stale 0.03) and keeps qc_complete=False on sample QC failure
    cr_findings = {
        "qc": {
            "status": "success",
            "data_type": "CUT_AND_RUN",
            "n_samples": 1,
            "metrics_not_computed": ["frip"],
            "qc_complete": False,
            "pass_qc": None,
            "warnings": ["Could not QC bad.bam: corrupt"],
        }
    }
    ChromatinAgent._merge_peak_frip_into_qc(
        cr_findings,
        {"data_type": "CUT_AND_RUN", "frip": 0.04, "frip_per_sample": {"s1": 0.04}},
    )
    assert cr_findings["qc"]["pass_qc"] is False
    assert cr_findings["qc"]["qc_complete"] is False
    assert any("FRiP" in f for f in cr_findings["qc"]["qc_failures"])


def test_render_frip_distribution_fallback_to_frip_per_sample(tmp_path):
    from aria.agents._narrative_chromatin import render_frip_distribution

    out = render_frip_distribution(
        {"data_type": "ChIP", "frip_per_sample": {"s1": 0.05, "s2": 0.12}},
        tmp_path / "frip.png",
    )
    assert out is not None and Path(out).is_file()
    assert Path(out).with_suffix(".svg").is_file()


# ── 3. BI-06 ─────────────────────────────────────────────────────────────────

def test_parse_bool_flag_and_na_spikein_cells(tmp_path):
    from aria.agents.chromatin_agent import ChromatinAgent
    from aria.scripts.chromatin_peak_counts import _is_control_meta, _parse_bool_flag

    for falsy in (False, "False", "false", "0", "no", "", "NA", None, 0, float("nan")):
        assert _parse_bool_flag(falsy) is False
        assert ChromatinAgent._parse_bool_flag(falsy) is False
        assert _is_control_meta({"is_control": falsy}) is False
    for truthy in (True, "True", "true", "1", "yes", "Y", 1):
        assert _parse_bool_flag(truthy) is True
        assert ChromatinAgent._parse_bool_flag(truthy) is True
        assert _is_control_meta({"is_control": truthy}) is True

    # Empty/NA spike-in cells in TSV metadata must not crash float()
    t1 = tmp_path / "s1.bam"
    t2 = tmp_path / "s2.bam"
    t1.write_bytes(b"")
    t2.write_bytes(b"")
    meta = tmp_path / "meta.tsv"
    meta.write_text(
        "sample_id\tis_control\tspikein_reads\tspikein_scale_factor\n"
        "s1\tFalse\t\tNA\n"
        "s2\t0\t15000\t0.66\n",
        encoding="utf-8",
    )
    res = _bare_agent()._resolve_epigenomic_samples(
        {"sample_metadata_path": str(meta)}, [str(t1), str(t2)]
    )
    assert res["target_files"] == [str(t1), str(t2)]
    assert res["control_files"] == []
    assert res["spikein_counts"] == {"s2": 15000.0}
    assert res["spikein_scale_factors"] == {"s2": 0.66}


# ── 4. D2-C3 ─────────────────────────────────────────────────────────────────

def test_atac_align_markdup_keeps_flags_and_view_strips_duplicates(monkeypatch, tmp_path):
    from aria.scripts import atac_align as aa

    commands = []

    def pipe(cmd1, cmd2, *a):
        commands.extend([cmd1, cmd2])
        return True

    def run(cmd, *a):
        commands.append(cmd)
        if cmd[:2] == ["samtools", "index"]:
            Path(cmd[-1] + ".bai").write_bytes(b"idx")
        return True

    def filt(cmd, out_bam, *a, **k):
        commands.append(cmd)
        Path(out_bam).write_bytes(b"bam")
        return True

    monkeypatch.setattr(aa, "_pipe", pipe)
    monkeypatch.setattr(aa, "_run", run)
    monkeypatch.setattr(aa, "_filter_bam", filt)
    monkeypatch.setattr(aa, "_bam_valid", lambda p: True)
    monkeypatch.setattr(aa, "_bam_metrics", lambda *a, **k: {"n_final": 10, "pct_dup": 12.5})
    monkeypatch.setattr(aa, "_bwa_index_exists", lambda p: True)

    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\nACGT\n")
    r1 = tmp_path / "R1.fq"
    r1.write_text("@r\nACGT\n+\nIIII\n")

    res = aa.atac_align({
        "samples": [{"name": "s", "r1": str(r1), "paired": False}],
        "genome_fasta": str(fa),
        "output_dir": str(tmp_path / "out"),
        "remove_duplicates": True,
    })
    assert res["status"] == "success"
    markdup_cmd = next(c for c in commands if c[:2] == ["samtools", "markdup"])
    assert "-r" not in markdup_cmd
    view_cmd = next(c for c in commands if c[:2] == ["samtools", "view"])
    assert view_cmd[view_cmd.index("-F") + 1] == str(aa.SAMTOOLS_EXCLUDE_FLAGS)


def test_count_spikein_index_fails_closed_when_bwa_errors(monkeypatch, tmp_path):
    from aria.scripts import atac_align as aa

    monkeypatch.setattr(aa, "_bwa_index_exists", lambda p: True)

    class _FailBwa:
        def __init__(self, cmd, stdout=None, stderr=None, **kw):
            self.stdout = io.BytesIO(b"")
            self.returncode = 1
            if stderr is not None and hasattr(stderr, "write"):
                stderr.write("bwa-mem2: corrupt index\n")

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def wait(self):
            return self.returncode

    monkeypatch.setattr(subprocess, "Popen", _FailBwa)
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **k: subprocess.CompletedProcess(a[0], 0, stdout="0\n", stderr=""),
    )
    warns = []
    got = aa._count_spikein_index("r1.fq", None, "sp_idx", "", 2, warns, "s1")
    assert got is None
    assert any("spike-in alignment failed" in w and "bwa-mem2 rc=1" in w for w in warns)


# ── 5. D2-C4 ─────────────────────────────────────────────────────────────────

def test_scatac_barcode_source_inclusive_coords_and_strand(tmp_path):
    from aria.scripts.chromatin_scatac_align import (
        _barcode_source,
        detect_barcode_whitelist_orientation,
    )

    assert _barcode_source("r1:0:-1,bc:0:15", ["r1.fq"], [], ["bc.fq"]) == (["bc.fq"], (0, 16))
    assert _barcode_source("bc:4:19", ["r1.fq"], [], []) == (["r1.fq"], (4, 20))
    assert _barcode_source("r1:16:-1,bc:0:15:-", ["r1.fq"], [], []) == (["r1.fq"], (0, 16))
    assert _barcode_source("bc:2:-1:+", ["r1.fq"], [], ["bc.fq"]) == (["bc.fq"], (2, None))

    wl = tmp_path / "wl.txt"
    bc_seq = "ACGTACGTACGTACGT"  # 16 bp
    wl.write_text(bc_seq + "\n")
    fq = tmp_path / "r1.fq"
    fq.write_text("".join(f"@r{i}\n{bc_seq}TTTTTTTT\n+\n{'I' * 24}\n" for i in range(25)))
    srcs, sl = _barcode_source("r1:16:-1,bc:0:15:-", [str(fq)], [], [])
    ev = detect_barcode_whitelist_orientation(srcs, str(wl), bc_slice=sl)
    assert ev["direct_rate"] == 1.0 and ev["orientation"] == "direct"


# ── 6. BI-07 & D2-H2 ─────────────────────────────────────────────────────────

def test_motifs_regions_from_csv_per_cluster_and_pseudobulk(tmp_path):
    from aria.scripts.chromatin_motifs import _regions_from_csv

    # Per-cluster CSV: negative LFC peaks excluded from cluster foreground
    pc_csv = tmp_path / "per_cluster.csv"
    pc_csv.write_text(
        "cluster,peak,log2fc,padj,significant\n"
        "0,chr1:10-20,2.0,1e-5,True\n"
        "0,chr1:30-40,-2.5,1e-6,True\n"
        "0,chr1:50-60,0.0,1e-4,True\n"
        "1,chr1:70-80,1.5,1e-3,True\n",
        encoding="utf-8",
    )
    pc_groups = _regions_from_csv(str(pc_csv), 5000, [])
    assert pc_groups == {"0": ["chr1:10-20"], "1": ["chr1:70-80"]}

    # Pseudobulk CSV: split by comparison::up and comparison::down
    pb_csv = tmp_path / "pseudobulk.csv"
    pb_csv.write_text(
        "comparison,peak,log2fc,padj,significant\n"
        "trt_vs_ctrl,chr1:100-200,1.8,1e-4,True\n"
        "trt_vs_ctrl,chr1:300-400,-1.4,1e-3,True\n"
        "trt_vs_ctrl,chr1:500-600,2.2,0.5,False\n",
        encoding="utf-8",
    )
    pb_groups = _regions_from_csv(str(pb_csv), 5000, [])
    assert pb_groups == {
        "trt_vs_ctrl::up": ["chr1:100-200"],
        "trt_vs_ctrl::down": ["chr1:300-400"],
    }


def test_motifs_dash_peaks_contig_aliasing_and_explicit_n_motifs(monkeypatch, tmp_path):
    from aria.scripts.chromatin_motifs import (
        _canonicalize_peak_for_genome,
        chromatin_motifs,
    )

    assert _canonicalize_peak_for_genome("chr1-100-500", {"chr1"}) == "chr1:100-500"
    assert _canonicalize_peak_for_genome("1-100-500", {"chr1"}) == "chr1:100-500"
    assert _canonicalize_peak_for_genome("chr1:100-500", {"1"}) == "1:100-500"
    assert _canonicalize_peak_for_genome("chrUn:100-500", {"chr1"}) is None

    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "ACGT" * 200 + "\n")
    (tmp_path / "g.fa.fai").write_text("chr1\t800\t6\t800\t801\n")
    meme = tmp_path / "custom.meme"
    meme.write_text("MEME version 4\nMOTIF MA0001.1 TF1\nMOTIF MA0002.1 TF2\n")

    captured = {}

    class _FakeDF:
        def to_dicts(self):
            return [{
                "id": "MA0001.1",
                "name": "TF1",
                "log2(fold change)": 1.5,
                "p-value": 1e-4,
                "adjusted p-value": 1e-3,
            }]

    fake_snap = types.SimpleNamespace(
        datasets=types.SimpleNamespace(read_motifs=lambda path: ["m1", "m2"]),
        tl=types.SimpleNamespace(
            motif_enrichment=lambda **kw: captured.update(kw) or {"g1": _FakeDF()}
        ),
    )
    monkeypatch.setitem(sys.modules, "snapatac2", fake_snap)

    res = chromatin_motifs({
        "motif_file": str(meme),
        "genome_fasta": str(fa),
        "regions": {"g1": ["1-100-200"]},
        "background": ["1-100-200", "chr1-300-400", "chrUn-10-50"],
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "success" and res["ran"] is True
    assert res["motif_source"]["n_motifs"] == 2
    assert captured["regions"] == {"g1": ["chr1:100-200"]}
    assert captured["background"] == ["chr1:100-200", "chr1:300-400"]


# ── 7. CDX-04 & D2-H7 ────────────────────────────────────────────────────────

def test_tobias_paired_null_controls_large_n_no_oom():
    from aria.scripts.chromatin_footprint_tobias import _null_label_controls

    rng = np.random.default_rng(0)
    n_pairs = 40  # 2**40 would OOM if list(range(...)) were materialized
    mat = rng.normal(0.0, 1.0, size=(5, 2 * n_pairs))
    res = _null_label_controls(
        mat, n_a=n_pairs, alpha=0.05, max_permutations=12, random_seed=42, paired=True
    )
    assert res["n_permutations"] == 12
    assert res["method"] == "within_donor_condition_label_swap_with_bh"


def test_tobias_paired_batch_confounding_and_tiny_significant_effect():
    from aria.scripts.chromatin_footprint_tobias import infer_replicate_footprints

    reps = {f"a{i}": "A" for i in range(4)} | {f"b{i}": "B" for i in range(4)}
    rep_ids = {f"a{i}": f"d{i}" for i in range(4)} | {f"b{i}": f"d{i}" for i in range(4)}
    table = {
        "M_POS": {
            "tf": "TF_TINY_POS",
            "motif_id": "M_POS",
            "n_sites": 100,
            "scores": {
                "a0": 1.00003, "a1": 2.00003, "a2": 3.00003, "a3": 4.00003,
                "b0": 1.00000, "b1": 2.00000, "b2": 3.00000, "b3": 4.00000,
            },
        },
        "M_NEG": {
            "tf": "TF_TINY_NEG",
            "motif_id": "M_NEG",
            "n_sites": 100,
            "scores": {
                "a0": 1.00000, "a1": 2.00000, "a2": 3.00000, "a3": 4.00000,
                "b0": 1.00003, "b1": 2.00003, "b2": 3.00003, "b3": 4.00003,
            },
        },
    }

    # Confounded batch even in paired design -> degrades to descriptive_only
    confounded_batches = {f"a{i}": "batch1" for i in range(4)} | {
        f"b{i}": "batch2" for i in range(4)
    }
    res_conf = infer_replicate_footprints(
        table, reps, "A", "B",
        replicate_ids=rep_ids,
        replicate_batches=confounded_batches,
    )
    assert res_conf["inference"]["status"] == "descriptive_only"
    assert res_conf["inference"]["reason"] == "batch_condition_confounded"
    assert res_conf["ranking_basis"]["pseudobulk"] is False
    assert [r["tf"] for r in res_conf["top_toward_A"]] == ["TF_TINY_POS"]
    assert [r["tf"] for r in res_conf["top_toward_B"]] == ["TF_TINY_NEG"]

    # Non-confounded batch (each donor assayed in one batch across both conditions) -> paired t-test
    matched_batches = {
        "a0": "b1", "b0": "b1", "a1": "b1", "b1": "b1",
        "a2": "b2", "b2": "b2", "a3": "b2", "b3": "b2",
    }
    res_ok = infer_replicate_footprints(
        table, reps, "A", "B",
        replicate_ids=rep_ids,
        replicate_batches=matched_batches,
    )
    assert res_ok["inference"]["status"] == "success"
    assert res_ok["n_significant"] == 2
    # Tiny significant effects (~3e-5) must NOT round to 0.0 and vanish from top lists
    assert [r["tf"] for r in res_ok["top_toward_A"]] == ["TF_TINY_POS"]
    assert [r["tf"] for r in res_ok["top_toward_B"]] == ["TF_TINY_NEG"]
    assert res_ok["top_toward_A"][0]["change"] == pytest.approx(3e-5)
    assert res_ok["top_toward_B"][0]["change"] == pytest.approx(-3e-5)


def test_tobias_aggregate_plots_filafy_matches_var2(monkeypatch, tmp_path):
    from aria.scripts import chromatin_footprint_tobias as mod

    bd = tmp_path / "bindetect"
    bed = bd / "FOSL2JUNDvar.2_MA1145.1" / "beds" / "FOSL2JUNDvar.2_MA1145.1_all.bed"
    bed.parent.mkdir(parents=True)
    bed.write_text("chr1\t10\t20\n")

    cmds = []
    monkeypatch.setattr(mod, "_run", lambda cmd, log: cmds.append(cmd))
    out = mod._aggregate_plots(bd, {"A": "a.bw", "B": "b.bw"}, ["FOSL2::JUND(var.2)"], tmp_path)
    assert "FOSL2::JUND(var.2)" in out
    assert out["FOSL2::JUND(var.2)"]["png"].endswith("FOSL2JUNDvar.2_aggregate.png")
    assert len(cmds) == 2


# ── 8. CDX-05 ────────────────────────────────────────────────────────────────

def test_lsi_cache_params_distinguishes_cohorts_and_mutations(tmp_path):
    from aria.scripts.chromatin_lsi_clustering import _cache_params

    a = tmp_path / "cohort_a.h5ad"
    b = tmp_path / "cohort_b.h5ad"
    a.write_bytes(b"cohort A")
    b.write_bytes(b"cohort B")

    pa = _cache_params({"data_path": str(a)})
    pb = _cache_params({"data_path": str(b)})
    assert pa["data_path"] == str(a.resolve())
    assert pb["data_path"] == str(b.resolve())
    assert pa != pb

    a.write_bytes(b"cohort A modified content")
    pa2 = _cache_params({"data_path": str(a)})
    assert pa2 != pa


# ── 9. CDX-06 ────────────────────────────────────────────────────────────────

def test_pseudobulk_da_drops_missing_donor_or_condition_cells(monkeypatch):
    import anndata as ad
    from scipy import sparse

    from aria.scripts import rna_bulk_de
    from aria.scripts.chromatin_diffacc import _pseudobulk_da

    X = np.array([
        [10, 1],  # valid d1 trt
        [20, 1],  # valid d1 trt
        [999, 9],  # missing donor (None) in trt -> MUST BE DROPPED
        [999, 9],  # missing donor (None) in ctrl -> MUST NOT form fake paired "None" donor!
        [30, 1],  # valid d2 ctrl
        [40, 1],  # valid d2 ctrl
        [888, 9],  # empty condition "" -> MUST BE DROPPED
        [777, 9],  # NaN donor -> MUST BE DROPPED
    ], dtype=float)
    a = ad.AnnData(sparse.csr_matrix(X))
    a.obs_names = [f"c{i}" for i in range(len(X))]
    a.var_names = ["p0", "p1"]
    a.obs["condition"] = ["trt", "trt", "trt", "ctrl", "ctrl", "ctrl", "", "ctrl"]
    a.obs["donor"] = ["d1", "d1", None, None, "d2", "d2", "d3", np.nan]

    captured = {}

    def fake_deseq2(counts, meta, *args, **kwargs):
        captured["counts"] = counts.copy()
        captured["meta"] = meta.copy()
        return {"status": "error", "error_type": "Fake", "details": "stop"}, []

    monkeypatch.setattr(rna_bulk_de, "_run_deseq2", fake_deseq2)
    warns = []
    _pseudobulk_da(
        a, "condition", "donor", [{"test": "trt", "reference": "ctrl"}],
        0.05, 0.5, 10, 2, 1, False, warns,
    )
    assert any("Dropped 4 cell(s)" in w for w in warns)
    assert set(captured["counts"].columns) == {"d1", "d2"}
    assert "None" not in captured["counts"].columns
    # Exact row alignment preserved after dropping rows!
    assert captured["counts"].loc["p0"].to_dict() == {"d1": 30.0, "d2": 70.0}


# ── 10. D2-H1, D2-H3, D2-H4 ──────────────────────────────────────────────────

def test_blacklist_gzip_corrupt_and_contig_aliasing(tmp_path):
    from aria.scripts.chromatin_peak_counts import _interval_index, _overlaps
    from aria.scripts.chromatin_peaks import _filter_blacklist, _read_peak_intervals

    # Contig normalization in _interval_index / _overlaps
    idx = _interval_index([("chr1", 100, 200), ("chrM", 10, 50)])
    assert _overlaps(idx, "1", 150, 160) is True
    assert _overlaps(idx, "chr1", 150, 160) is True
    assert _overlaps(idx, "MT", 20, 30) is True
    assert _overlaps(idx, "1", 200, 250) is False

    # Gzipped UCSC blacklist filters Ensembl-named peaks
    bl_gz = tmp_path / "hg38-blacklist.v2.bed.gz"
    with gzip.open(bl_gz, "wt", encoding="utf-8") as fh:
        fh.write("chr1\t100\t300\tHigh_Mappability\n")
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("1\t150\t250\tp1\n1\t500\t600\tp2\n", encoding="utf-8")
    removed = _filter_blacklist(str(peaks), str(bl_gz))
    assert removed == 1
    assert peaks.read_text(encoding="utf-8").strip() == "1\t500\t600\tp2"

    # Corrupt binary file does not raise UnicodeDecodeError
    bad = tmp_path / "corrupt.bed"
    bad.write_bytes(b"\xff\xfe\x00\x80\x81")
    assert _read_peak_intervals(str(bad)) == []


def test_copula_idr_degenerate_w_sum_updates_params(monkeypatch):
    from aria.utils import seacr_idr

    # Force gamma sum <= 1e-6 via strongly anti-correlated latent coordinates
    # (z1=+50, z2=-50 -> log_f1 << log_f0 since rho > 0, so gamma -> 1e-9).
    calls = {"n": 0}

    def fake_invert(u, pi1, mu, sigma):
        calls["n"] += 1
        return np.full_like(u, 50.0 if calls["n"] % 2 == 1 else -50.0)

    monkeypatch.setattr(seacr_idr, "_invert_mixture_cdf", fake_invert)
    x = np.arange(1, 26, dtype=float)
    fit = seacr_idr.fit_copula_idr_pair(x, x)
    assert fit["converged"] is True
    assert fit["n_iter"] == 1
    assert fit["params"]["pi1"] == pytest.approx(0.01)
    assert fit["params"]["rho"] == pytest.approx(0.0)


def test_seacr_no_excess_returns_zero_peaks_and_tiny_merged_signal_preserved(tmp_path):
    from aria.utils.seacr_idr import _select_seacr_empirical_threshold, call_seacr_peaks

    # Identical target & control: max excess <= 0 -> inf (0 peaks), zero RuntimeWarnings
    arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        thr, meta = _select_seacr_empirical_threshold(arr, arr, mode="stringent")
    runtime_warns = [w for w in caught if issubclass(w.category, RuntimeWarning)]
    assert runtime_warns == []
    assert math.isinf(thr) and thr > 0
    assert meta.get("no_signal_excess") is True

    bg_a = tmp_path / "a.bedgraph"
    bg_b = tmp_path / "b.bedgraph"
    bg_a.write_text("chr1\t0\t100\t10.0\nchr1\t500\t600\t2.0\n")
    bg_b.write_text("chr1\t0\t100\t10.0\nchr1\t500\t600\t2.0\n")
    # Identical target vs control -> 0 peaks and JSON-safe threshold=None
    res_null = call_seacr_peaks(
        [str(bg_a)], tmp_path / "null.bed", control_inputs=[str(bg_b)], work_dir=tmp_path / "wn"
    )
    assert res_null["n_peaks"] == 0
    assert res_null["threshold"] is None
    json.dumps(res_null, allow_nan=False)

    # Tiny merged scale factor 1e-10 survives cov > 0.0, _fmt_bedgraph_value, and full-precision threshold
    res_tiny = call_seacr_peaks(
        [str(bg_a), str(bg_b)],
        tmp_path / "tiny.bed",
        top_quantile=0.5,
        target_scale_factors=[1e-10, 1e-10],
        work_dir=tmp_path / "wt",
    )
    assert res_tiny["n_peaks"] == 1
    assert res_tiny["threshold"] == pytest.approx(1e-7)
    cols = (tmp_path / "tiny.bed").read_text().strip().split("\t")
    assert float(cols[3]) == pytest.approx(1e-7)
    assert float(cols[4]) == pytest.approx(1e-9)


# ── 11. BI-09, BI-10, D2-H5, D2-H8, D2-H9 ───────────────────────────────────

def test_narrator_tobias_descriptive_batch_confounded_prose():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    fp = {
        "ran": True,
        "group_a": "Treat",
        "group_b": "Ctrl",
        "group_label": "Conditions",
        "group_kind": "conditions",
        "differential_summary": {
            "parsed": True,
            "n_motifs_tested": 20,
            "n_ranked_candidates": 20,
            "ranking_basis": {
                "fdr_controlled": False,
                "replicate_inference": False,
                "pseudobulk": False,
            },
            "inference": {
                "status": "descriptive_only",
                "reason": "batch_condition_confounded",
            },
        },
    }
    narr = ChromatinNarrator()
    block = narr._footprint_block(fp)
    assert "pseudobulk" not in block.claim.lower()
    assert "uncorrected p" not in block.claim.lower()
    assert "confounded" in block.claim.lower()
    assert all("pseudobulk" not in c.text.lower() for c in block.caveats)

    methods = narr._footprint_methods_paragraph(fp)
    assert "pseudobulk" not in methods.lower()
    assert "replicate-mean" in methods.lower()
    assert "batch_condition_confounded" in methods


@pytest.mark.parametrize(
    "data_type,analysis",
    [
        ("ChIP", "chip_alignment"),
        ("CUT_AND_RUN", "cut_and_run_alignment"),
        ("CUT_AND_TAG", "cut_and_tag_alignment"),
    ],
)
def test_narrator_bulk_methods_matches_epigenomic_alignment_names(data_type, analysis):
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    findings = {
        "preprocessing": {
            "analysis": analysis,
            "data_type": data_type,
            "ran": True,
            "aligner": "bwa-mem2",
            "assembly": "mm39",
            "n_samples_aligned": 4,
            "filter": {"mapq": 20, "samtools_f": 2, "samtools_F": 780, "remove_duplicates": False},
        },
        "peaks": {"status": "success", "data_type": data_type, "n_peaks": 500, "genome": "mm39"},
    }
    txt = " ".join(ChromatinNarrator()._bulk_atac_methods(findings))
    assert "bwa-mem2" in txt
    assert "does not realign" not in txt


def test_narrator_bulk_da_skipped_emits_limitation_and_spikein_used_detected():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    narr = ChromatinNarrator()
    for da_payload in (
        {"status": "skipped", "data_type": "bulk_ATAC", "reason": "sample_metadata_missing"},
        {"status": "success", "ran": False, "data_type": "ChIP", "reason": "under_replicated"},
    ):
        blocks = narr.collect("chromatin_agent", {"findings": {"differential_accessibility": da_payload}})
        assert len(blocks) == 1
        assert blocks[0].block_type == "limitation"
        assert blocks[0].status == "limitation"
        assert da_payload["reason"] in blocks[0].caveats[0].text

    da_used = {
        "status": "success",
        "ran": True,
        "data_type": "CUT_AND_RUN",
        "n_comparisons_success": 1,
        "n_sig_total": 12,
        "spikein_normalization": {"used": True},
    }
    blk = narr._bulk_atac_diffacc_block(da_used)
    ev = {e.label: e.value for e in blk.evidence}
    assert ev["Normalization mode"] == "heterologous_spikein"
    assert any("spike-in sizeFactors" in c.text for c in blk.caveats)


def test_narrator_multimodality_collects_all_modalities_without_overwriting():
    from aria.agents.narrative.narrators.chromatin import (
        ChromatinNarrator,
        unwrap_chromatin_findings,
    )

    agent_res = {
        "status": "done",
        "findings": {
            "bulk_ATAC": {
                "status": "done",
                "findings": {
                    "qc": {"status": "success", "data_type": "bulk_ATAC", "n_samples": 4, "mito_fraction": 0.02},
                    "peaks": {"status": "success", "data_type": "bulk_ATAC", "n_peaks": 12000, "genome": "hg38", "peaks_path": "/a.bed"},
                    "motifs": {"status": "success", "ran": True, "per_group": {"g": {"n_enriched": 2}}, "motif_source": {"collection": "JASPAR", "n_motifs": 10}},
                    "figures": {"qc_frip": "/atac_frip.png"},
                },
            },
            "ChIP": {
                "status": "done",
                "findings": {
                    "qc": {"status": "success", "data_type": "ChIP", "n_samples": 6, "mito_fraction": 0.01},
                    "peaks": {"status": "success", "data_type": "ChIP", "n_peaks": 4500, "genome": "hg38", "peaks_path": "/c.bed"},
                    "motifs": {"status": "success", "ran": True, "per_group": {"g": {"n_enriched": 1}}, "motif_source": {"collection": "JASPAR", "n_motifs": 10}},
                    "figures": {"qc_frip": "/chip_frip.png"},
                },
            },
        },
    }
    narr = ChromatinNarrator()
    blocks = narr.collect("chromatin_agent", agent_res)
    ids = [b.id for b in blocks]
    assert len(ids) == len(set(ids)), f"duplicate block IDs: {ids}"
    assert "chromatin.qc.bulk_ATAC" in ids and "chromatin.qc.ChIP" in ids
    assert "chromatin.peak_calling.bulk_ATAC" in ids and "chromatin.peak_calling.ChIP" in ids
    assert sum(b.analysis == "motif_enrichment" for b in blocks) == 2

    methods = " ".join(narr.methods("chromatin_agent", agent_res))
    assert "Bulk ATAC" in methods and "ChIP-seq" in methods

    figs = narr.figures("chromatin_agent", agent_res)
    assert {f["path"] for f in figs} == {"/atac_frip.png", "/chip_frip.png"}
    assert len({f["id"] for f in figs}) == 2

    tables = narr.tables("chromatin_agent", agent_res)
    assert {t["path"] for t in tables} == {"/a.bed", "/c.bed"}
    assert len({t["id"] for t in tables}) == 2

    # Single-modality unwrap remains backwards-compatible
    single = {"findings": {"ChIP": {"findings": {"qc": {"n_samples": 3}}}}}
    assert unwrap_chromatin_findings(single) == {"qc": {"n_samples": 3}}


def test_render_motif_dotplot_padj_filter_and_small_legend(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt

    from aria.agents._narrative_chromatin import render_motif_dotplot

    # All motifs >= padj_max -> returns None
    nonsig = {
        "padj_max": 0.05,
        "per_group": {
            "0": {
                "top_motifs": [
                    {"name": "SPI1", "log2_enrichment": 1.8, "padj": 0.05},
                    {"name": "CTCF", "log2_enrichment": 2.5, "padj": 0.42},
                ]
            }
        },
    }
    assert render_motif_dotplot(nonsig, tmp_path / "none.png") is None

    # Small positive -log10(padj) (when padj_max is relaxed) formats cleanly without collapsing to "0"
    captured_figs = []
    orig_close = plt.close
    monkeypatch.setattr(plt, "close", lambda fig=None: captured_figs.append(fig))
    try:
        small = {
            "padj_max": 0.9,
            "per_group": {
                "0": {
                    "top_motifs": [
                        {"name": "SPI1", "log2_enrichment": 1.5, "padj": 0.8},
                    ]
                }
            },
        }
        out = render_motif_dotplot(small, tmp_path / "small.png")
        assert out is not None and Path(out).is_file()
        leg = captured_figs[-1].axes[0].get_legend()
        labels = [t.get_text() for t in leg.get_texts()]
        assert len(labels) == 3
        assert all(lbl != "0" and float(lbl) > 0 for lbl in labels)
    finally:
        for f in captured_figs:
            orig_close(f)
