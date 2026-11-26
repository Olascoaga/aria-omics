"""Permanent regression tests for D1 silent-failure fixes in Hi-C compute kernels."""
from __future__ import annotations

import sys
import types

import pandas as pd
import pytest


def test_d1_001_hic_preprocess_must_not_drop_chrom_contacts(tmp_path, monkeypatch):
    """D1-001: Corrupt chromosome block during .hic -> .cool conversion must fail closed."""
    pytest.importorskip("cooler")
    pytest.importorskip("bioframe")
    from aria.scripts.hic_preprocess import _convert_hic_to_cooler

    class _C:
        def __init__(self, name, length):
            self.name, self.length = name, length

    class _R:
        def __init__(self, bx, by, c):
            self.binX, self.binY, self.counts = bx, by, c

    class _MZD:
        def __init__(self, c1, c2):
            self.c1, self.c2 = c1, c2

        def getRecords(self, *a):
            if self.c1 == "chr1" and self.c2 == "chr1":
                raise RuntimeError("corrupt compressed block on chr1")
            return [_R(0, 0, 50), _R(0, 10000, 20)] if self.c1 == self.c2 == "chr2" else []

    class _HiC:
        def __init__(self, p):
            pass

        def getResolutions(self):
            return [10000]

        def getChromosomes(self):
            return [_C("chr1", 30000), _C("chr2", 30000)]

        def getMatrixZoomData(self, c1, c2, *a):
            return _MZD(c1, c2)

    mod = types.ModuleType("hicstraw")
    mod.HiCFile = _HiC
    monkeypatch.setitem(sys.modules, "hicstraw", mod)
    with pytest.raises(RuntimeError, match="corrupt compressed block"):
        _convert_hic_to_cooler(tmp_path / "in.hic", tmp_path / "out.cool", 10000, "hg38")


def _make_compartment_cooler_and_gtf(tmp_path):
    """Build a small balanced cooler with strong A/B block signal and a matching GTF."""
    import cooler
    import numpy as np

    n_bins = 24
    res = 100_000
    bins = pd.DataFrame({
        "chrom": ["chr1"] * n_bins,
        "start": [i * res for i in range(n_bins)],
        "end": [(i + 1) * res for i in range(n_bins)],
    })
    rows = []
    # Two strong compartment blocks (first half A, second half B) + distance decay
    for i in range(n_bins):
        for j in range(i, n_bins):
            same_comp = (i < n_bins // 2) == (j < n_bins // 2)
            base = 200.0 / (1.0 + (j - i) * 0.35)
            cnt = int(round(base * (3.0 if same_comp else 0.5)))
            if cnt > 0:
                rows.append((i, j, cnt))
    cool_p = tmp_path / "ab.cool"
    cooler.create_cooler(
        str(cool_p),
        bins,
        pd.DataFrame(rows, columns=["bin1_id", "bin2_id", "count"]),
        assembly="hg38",
    )
    clr = cooler.Cooler(str(cool_p))
    cooler.balance_cooler(clr, store=True, mad_max=0, min_nnz=1, ignore_diags=1)

    # GTF with higher gene density in first half (A compartment) than second half
    gtf_lines = []
    for i in range(n_bins // 2):
        for k in range(4):
            pos = i * res + 1000 + k * 5000
            gtf_lines.append(
                f'chr1\tHAVANA\tgene\t{pos}\t{pos + 1000}\t.\t+\t.\tgene_id "G{i}_{k}";\n'
            )
    gtf_lines.append(
        f'chr1\tHAVANA\tgene\t{(n_bins - 2) * res + 1000}\t{(n_bins - 2) * res + 2000}\t.\t+\t.\tgene_id "GB1";\n'
    )
    gtf_p = tmp_path / "genes.gtf"
    gtf_p.write_text("".join(gtf_lines))
    return cool_p, gtf_p


def test_d1_020_corrupt_fasta_must_not_be_silently_swallowed(tmp_path):
    """D1-020: Corrupt FASTA sequence body detected during fetch after open must fail closed."""
    cooler = pytest.importorskip("cooler")
    pysam = pytest.importorskip("pysam")
    from aria.scripts.hic_topology import _build_phasing_track

    bins = pd.DataFrame({
        "chrom": ["chr1"] * 10,
        "start": [i * 100_000 for i in range(10)],
        "end": [(i + 1) * 100_000 for i in range(10)],
        "weight": [1.0] * 10,
    })
    cool_p = tmp_path / "c.cool"
    cooler.create_cooler(
        str(cool_p),
        bins,
        pd.DataFrame([(0, 1, 10)], columns=["bin1_id", "bin2_id", "count"]),
        assembly="hg38",
    )
    # Valid header + .fai, then truncated sequence body so open succeeds and fetch fails.
    bad_fa = tmp_path / "corrupt_trunc.fa"
    bad_fa.write_text(">chr1\n" + ("ACGT" * 20 + "\n") * 200)
    with pysam.FastaFile(str(bad_fa)):
        pass  # builds corrupt_trunc.fa.fai
    bad_fa.write_text(">chr1\nACGT\n")
    with pytest.raises((OSError, ValueError, RuntimeError)):
        _build_phasing_track(cooler.Cooler(str(cool_p)), ["chr1"], str(bad_fa), None)


def test_r3_d1_020_plain_gzip_fasta_falls_back_to_gtf_with_reason(tmp_path):
    """R3/D1-020 (a): Valid plain-gzip FASTA must fall back to GTF phasing with a recorded reason."""
    import gzip

    pytest.importorskip("cooler")
    pytest.importorskip("cooltools")
    pytest.importorskip("pysam")
    from aria.scripts.hic_topology import _compute_compartments

    cool_p, gtf_p = _make_compartment_cooler_and_gtf(tmp_path)
    gz_fa = tmp_path / "hg38.fa.gz"
    gz_fa.write_bytes(gzip.compress((">chr1\n" + "ACGT" * 100 + "\n").encode("utf-8")))

    res = _compute_compartments(
        str(cool_p), "hg38", "Homo sapiens", 100_000, ["chr1"], str(tmp_path / "out_gz"),
        genome_fasta=str(gz_fa), gtf_path=str(gtf_p),
    )
    assert res["status"] == "success", res
    assert res["phasing_method"] == "gene_density"
    assert res["pc1_validated"] is True
    reason = res.get("phasing_fallback_reason")
    assert isinstance(reason, str) and "hg38.fa.gz" in reason
    assert any("hg38.fa.gz" in w for w in (res.get("warnings") or []))


def test_r3_d1_020_readonly_dir_fasta_without_fai_falls_back_to_gtf_with_reason(tmp_path):
    """R3/D1-020 (b): Valid FASTA in a read-only dir without .fai must fall back to GTF with a recorded reason."""
    import os

    pytest.importorskip("cooler")
    pytest.importorskip("cooltools")
    pytest.importorskip("pysam")
    from aria.scripts.hic_topology import _compute_compartments

    if hasattr(os, "geteuid") and os.geteuid() == 0:
        pytest.skip("read-only dir permission check is bypassed by root")

    cool_p, gtf_p = _make_compartment_cooler_and_gtf(tmp_path)
    ro_dir = tmp_path / "ro_ref"
    ro_dir.mkdir()
    ro_fa = ro_dir / "ref.fa"
    ro_fa.write_text(">chr1\n" + "ACGT" * 100 + "\n")
    os.chmod(ro_dir, 0o555)
    try:
        res = _compute_compartments(
            str(cool_p), "hg38", "Homo sapiens", 100_000, ["chr1"], str(tmp_path / "out_ro"),
            genome_fasta=str(ro_fa), gtf_path=str(gtf_p),
        )
    finally:
        os.chmod(ro_dir, 0o755)
    assert res["status"] == "success", res
    assert res["phasing_method"] == "gene_density"
    assert res["pc1_validated"] is True
    reason = res.get("phasing_fallback_reason")
    assert isinstance(reason, str) and "ref.fa" in reason
    assert any("ref.fa" in w for w in (res.get("warnings") or []))


def test_d1_007_hic_qc_must_not_drop_corrupt_sample_silently(tmp_path):
    """D1-007: Corrupt Hi-C sample in a multi-sample cohort must fail closed."""
    cooler = pytest.importorskip("cooler")
    from aria.scripts.hic_qc_and_balance import hic_qc_and_balance

    bins = pd.DataFrame({
        "chrom": ["chr1"] * 12,
        "start": [i * 100_000 for i in range(12)],
        "end": [(i + 1) * 100_000 for i in range(12)],
        "weight": [1.0] * 12,
    })
    rows = [(i, j, max(1, 50 - 3 * (j - i))) for i in range(12) for j in range(i, 12)]
    good = tmp_path / "good.cool"
    cooler.create_cooler(
        str(good),
        bins,
        pd.DataFrame(rows, columns=["bin1_id", "bin2_id", "count"]),
        assembly="hg38",
    )
    bad = tmp_path / "corrupt.cool"
    bad.write_bytes(b"NOT_AN_HDF5_COOLER")
    res = hic_qc_and_balance({
        "files": [str(good), str(bad)],
        "sample_ids": ["good", "corrupt"],
        "resolution": 100_000,
        "min_valid_pairs": 1,
        "min_cis_trans_ratio": 0.0,
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "error"


def test_d1_004_phasing_valid_chroms_must_fail_closed_on_ice_error():
    """D1-004: ICE weight-read error must fail closed (return set()) rather than validating unmasked."""
    from aria.scripts.hic_topology import _phasing_valid_chroms

    phasing_df = pd.DataFrame({
        "chrom": ["chr1"] * 10,
        "start": [i * 100_000 for i in range(10)],
        "end": [(i + 1) * 100_000 for i in range(10)],
        "GC": [0.35 + 0.03 * i for i in range(10)],
    })

    class _BrokenBins:
        columns = ["chrom", "start", "end", "weight"]

        def __getitem__(self, key):
            raise OSError("corrupt HDF5 weight dataset")

    class _BrokenCooler:
        def bins(self):
            return _BrokenBins()

    assert _phasing_valid_chroms(phasing_df, ["chr1"], clr=_BrokenCooler()) == set()


def test_d1_019_hic_inspect_pairs_corrupt_must_report_error(tmp_path):
    """D1-019: Corrupt .pairs.gz must include an 'error' key in hic_inspect per_file entry."""
    from aria.scripts.hic_inspect import hic_inspect

    bad = tmp_path / "corrupt.pairs.gz"
    bad.write_bytes(b"NOT_GZIP_BYTES")
    res = hic_inspect({"files": [str(bad)]})
    assert "error" in res["per_file"][0] and res["per_file"][0]["chromosomes"] == []
