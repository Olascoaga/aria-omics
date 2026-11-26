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


def test_d1_020_corrupt_fasta_must_not_be_silently_swallowed(tmp_path):
    """D1-020: Corrupt user-supplied genome FASTA must fail closed instead of falling back to unphased."""
    cooler = pytest.importorskip("cooler")
    pytest.importorskip("pysam")
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
    bad_fa = tmp_path / "corrupt.fa"
    bad_fa.write_bytes(b"\x00\xff\xfeNOT_A_FASTA")
    with pytest.raises((OSError, ValueError, RuntimeError)):
        _build_phasing_track(cooler.Cooler(str(cool_p)), ["chr1"], str(bad_fa), None)


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
