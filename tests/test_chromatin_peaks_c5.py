from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace


def test_overlap_consensus_requires_replicate_support(tmp_path):
    from aria.scripts.chromatin_peaks import _write_overlap_consensus

    rep1 = tmp_path / "rep1.narrowPeak"
    rep2 = tmp_path / "rep2.narrowPeak"
    rep1.write_text("chr1\t100\t180\nchr1\t500\t600\n", encoding="utf-8")
    rep2.write_text("chr1\t120\t200\nchr1\t900\t1000\n", encoding="utf-8")

    out = tmp_path / "consensus.narrowPeak"
    res = _write_overlap_consensus([str(rep1), str(rep2)], str(out),
                                   min_support=2)

    assert res["ran"] is True
    assert res["method"] == "overlap_support_consensus"
    assert res["n_candidate_regions"] == 3
    assert res["n_reproducible_regions"] == 1
    assert res["fraction_reproducible"] == 0.3333
    assert out.read_text(encoding="utf-8").splitlines() == [
        "chr1\t100\t200\tsupport_2\t2"
    ]


def test_chromatin_peaks_runs_replicate_overlap_policy(tmp_path, monkeypatch):
    import subprocess
    import aria.scripts.chromatin_peaks as peaks_mod

    bam1 = tmp_path / "a.bam"
    bam2 = tmp_path / "b.bam"
    bam1.write_text("bam", encoding="utf-8")
    bam2.write_text("bam", encoding="utf-8")

    def fake_run(cmd, capture_output=True, text=True, timeout=None, **kwargs):
        assert cmd[:2] == ["macs3", "callpeak"]
        name = cmd[cmd.index("-n") + 1]
        outdir = Path(cmd[cmd.index("--outdir") + 1])
        outdir.mkdir(parents=True, exist_ok=True)
        if name.startswith("rep1_"):
            rows = ["chr1\t100\t180\n", "chr1\t500\t600\n"]
        elif name.startswith("rep2_"):
            rows = ["chr1\t120\t200\n", "chr1\t900\t1000\n"]
        else:
            rows = ["chr1\t100\t200\n", "chr1\t500\t600\n"]
        (outdir / f"{name}_peaks.narrowPeak").write_text(
            "".join(rows), encoding="utf-8")
        return SimpleNamespace(returncode=0, stderr="", stdout="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(peaks_mod, "_compute_frip", lambda *args, **kw: 0.42)

    res = peaks_mod.chromatin_peaks({
        "data_type": "bulk_ATAC",
        "files": [str(bam1), str(bam2)],
        "genome": "hg38",
        "macs3_params": {
            "format": "BAMPE",
            "nomodel": True,
            "extsize": 200,
            "keep_dup": "all",
        },
        "output_dir": str(tmp_path / "peaks"),
    })

    assert res["status"] == "success"
    assert res["peak_calling_strategy"] == (
        "pooled_macs3_with_replicate_overlap_qc")
    assert res["consensus_peaks_path"].endswith(
        "reproducible_consensus.narrowPeak")
    repro = res["peak_reproducibility"]
    assert repro["status"] == "verified"
    assert repro["replicate_peak_calling"]["ran"] is True
    assert repro["replicate_peak_calling"]["n_peak_files"] == 2
    assert repro["overlap"]["n_reproducible_regions"] == 1
    assert repro["idr"] == {
        "ran": False,
        "reason": "idr_not_run_overlap_reproducibility_policy",
    }
    assert "IDR was not run" in " ".join(res["warnings"])


def test_chromatin_narrator_surfaces_c5_reproducibility_policy():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.render_blocks import render_blocks

    repro = {
        "status": "verified",
        "strategy": "pooled_macs3_with_replicate_overlap_qc",
        "overlap": {
            "ran": True,
            "method": "overlap_support_consensus",
            "support_threshold": 2,
            "n_candidate_regions": 3,
            "n_reproducible_regions": 1,
        },
        "idr": {
            "ran": False,
            "reason": "idr_not_run_overlap_reproducibility_policy",
        },
    }
    findings = {
        "bulk_ATAC": {
            "status": "done",
            "findings": {
                "peaks": {
                    "status": "success",
                    "data_type": "bulk_ATAC",
                    "n_peaks": 2,
                    "genome": "hg38",
                    "peaks_path": "/tmp/pooled.narrowPeak",
                    "consensus_peaks_path": "/tmp/consensus.narrowPeak",
                    "frip": 0.42,
                    "peak_reproducibility": repro,
                    "macs3_cmd": "macs3 callpeak -t a.bam b.bam",
                    "warnings": [],
                }
            },
        }
    }

    blocks = ChromatinNarrator().collect(
        "chromatin_agent", {"status": "done", "findings": findings}, {})
    peak_block = next(b for b in blocks
                      if b.id == "chromatin.peak_calling.bulk_ATAC")
    assert peak_block.metrics["peak_reproducibility"]["status"] == "verified"
    html = render_blocks(blocks, strict=True)
    assert "Overlap reproducible peaks" in html
    assert "IDR was not run" in html

    methods = ChromatinNarrator().methods(
        "chromatin_agent", {"status": "done", "findings": findings}, {})
    blob = " ".join(methods)
    assert "support >= 2 replicate peak sets" in blob
    assert "IDR was not run" in blob
