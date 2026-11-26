"""Regression tests for D4 silent data misalignment fixes in chromatin / epigenomics."""
from __future__ import annotations

from pathlib import Path

import pytest


def _write_bam(path: Path, chrom: str, chrom_len: int, positions):
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": chrom, "LN": chrom_len}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (start, length) in enumerate(positions):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence, a.flag = f"r{i}", "A" * length, 0
            a.reference_id, a.reference_start, a.mapping_quality = 0, start, 60
            a.cigar = [(0, length)]
            a.query_qualities = pysam.qualitystring_to_array("I" * length)
            out.write(a)
    pysam.index(str(path))


def test_d4_007_seacr_and_broad_domains_fail_closed_on_scale_factor_length_mismatch(tmp_path: Path):
    """D4-007: Non-empty scale_factors shorter/longer than input BAMs must fail closed
    in call_seacr_peaks, call_broad_domains, and chromatin_peaks (never mix units)."""
    pytest.importorskip("pysam")
    from aria.scripts.chromatin_peaks import chromatin_peaks
    from aria.utils.seacr_idr import call_broad_domains, call_seacr_peaks

    pos = [(10_000 + i * 100, 50) for i in range(200)]
    bam1, bam2 = tmp_path / "rep1.bam", tmp_path / "rep2.bam"
    _write_bam(bam1, "chr1", 200_000, pos)
    _write_bam(bam2, "chr1", 200_000, pos)

    res_seacr = call_seacr_peaks(
        [str(bam1), str(bam2)],
        tmp_path / "seacr.bed",
        target_scale_factors=[0.001],
        work_dir=tmp_path / "work",
    )
    assert res_seacr.get("status") == "error" and res_seacr.get("error_type") == "InvalidScaleFactors"

    res_bd = call_broad_domains(
        [str(bam1), str(bam2)],
        tmp_path / "broad.bed",
        bin_size=10_000,
        ip_scale_factors=[0.01],
    )
    assert res_bd.get("status") == "error" and res_bd.get("error_type") == "InvalidScaleFactors"

    res_cp = chromatin_peaks({
        "data_type": "CUT_AND_RUN",
        "files": [str(bam1), str(bam2)],
        "peak_caller": "broad_domain",
        "target_scale_factors": [0.01],
        "output_dir": str(tmp_path / "cp_out"),
        "write_bigwig": True,
    })
    assert res_cp.get("status") == "error" and res_cp.get("error_type") == "InvalidScaleFactors"


def test_d4_008_broad_domains_normalizes_contigs_and_fails_on_zero_control_overlap(tmp_path: Path):
    """D4-008: call_broad_domains must match UCSC chr1 vs Ensembl 1 controls via
    canonical contig normalization, and fail closed if control_bams share 0 contigs."""
    pytest.importorskip("pysam")
    from aria.utils.seacr_idr import call_broad_domains

    pos = [
        (b * 10_000 + k * 50, 50)
        for b in range(20)
        for k in range(100 if 5 <= b <= 7 else 2)
    ]
    ip_bam = tmp_path / "ip.bam"
    ctrl_ucsc = tmp_path / "ctrl_ucsc.bam"
    ctrl_ens = tmp_path / "ctrl_ens.bam"
    ctrl_other = tmp_path / "ctrl_chr2.bam"
    _write_bam(ip_bam, "chr1", 200_000, pos)
    _write_bam(ctrl_ucsc, "chr1", 200_000, pos)
    _write_bam(ctrl_ens, "1", 200_000, pos)
    _write_bam(ctrl_other, "chr2", 200_000, pos)

    res_ok = call_broad_domains([str(ip_bam)], tmp_path / "ok.bed", control_bams=[str(ctrl_ucsc)])
    assert res_ok["status"] == "success" and res_ok["n_peaks"] == 0

    res_ens = call_broad_domains([str(ip_bam)], tmp_path / "ens.bed", control_bams=[str(ctrl_ens)])
    assert res_ens["status"] == "success" and res_ens["n_peaks"] == 0, (
        f"Expected chr1 vs 1 control match (0 peaks), got {res_ens}"
    )

    res_disjoint = call_broad_domains(
        [str(ip_bam)], tmp_path / "disjoint.bed", control_bams=[str(ctrl_other)]
    )
    assert res_disjoint["status"] == "error" and res_disjoint["error_type"] == "ControlContigMismatch"
