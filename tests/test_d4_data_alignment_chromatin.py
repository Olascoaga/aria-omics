"""Regression tests for D4 silent data misalignment fixes in chromatin / epigenomics."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse


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


def test_d4_010_label_transfer_rejects_duplicate_barcodes(tmp_path: Path):
    """D4-010: _label_transfer must reject non-unique barcodes in scATAC or scRNA
    instead of Cartesian-expanding cell counts."""
    ad = pytest.importorskip("anndata")
    from aria.scripts.chromatin_regulatory import _label_transfer

    bcs = [f"AAACAGCCAAG{i:05d}-1" for i in range(60)]
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.ones((60, 2), dtype=np.float32)),
        obs=pd.DataFrame({"leiden": ["0"] * 60}, index=bcs),
        var=pd.DataFrame(index=["chr1:1000-1500", "chr1:2000-2500"]),
    )
    rna_dup = ad.AnnData(
        X=sparse.csr_matrix(np.ones((120, 2), dtype=np.float32)),
        obs=pd.DataFrame({"cell_type": ["T_cell"] * 60 + ["B_cell"] * 60}, index=bcs + bcs),
        var=pd.DataFrame(index=["G0", "G1"]),
    )
    rna_dup_path = tmp_path / "rna_dup.h5ad"
    rna_dup.write_h5ad(rna_dup_path)

    res = _label_transfer(
        adata,
        {"rna_data_path": str(rna_dup_path), "rna_label_col": "cell_type"},
        tmp_path,
    )
    assert res["ran"] is False, f"Duplicate barcodes were not rejected: {res}"
    assert "duplicate" in res.get("reason", "").lower() or "unique" in res.get("reason", "").lower()

    rna_ok = ad.AnnData(
        X=sparse.csr_matrix(np.ones((60, 2), dtype=np.float32)),
        obs=pd.DataFrame({"cell_type": ["T_cell"] * 40 + ["B_cell"] * 20}, index=bcs),
        var=pd.DataFrame(index=["G0", "G1"]),
    )
    rna_ok_path = tmp_path / "rna_ok.h5ad"
    rna_ok.write_h5ad(rna_ok_path)
    res_ok = _label_transfer(
        adata,
        {"rna_data_path": str(rna_ok_path), "rna_label_col": "cell_type"},
        tmp_path,
    )
    assert res_ok["ran"] is True
    df_ok = pd.read_csv(res_ok["output_csv"])
    assert int(df_ok["n_cells"].sum()) == 60
