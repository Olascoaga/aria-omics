"""scATAC P0 W0.5: real ENCODE QC from fragments — TSS enrichment (snapatac2 +
governed annotation) and FRiP (fraction of reads in peaks). The heavy snapatac2
compute is env/asset/egress-gated; this asserts the HONEST fallbacks (None + a
concrete reason, never a fabricated number) and the pure FRiP overlap computation.
"""

from __future__ import annotations

import sys
import types

import pytest


# ── TSS enrichment — honest fallbacks (no fabrication) ───────────────────────

def test_tss_enrichment_unknown_assembly_returns_reason():
    from aria.scripts.chromatin_qc import _compute_tss_enrichment

    value, reason = _compute_tss_enrichment("frags.tsv.gz", "ce11")
    assert value is None
    assert reason and "assembly" in reason.lower()


def test_tss_enrichment_air_gapped_skips_without_egress(monkeypatch):
    from aria.scripts import chromatin_qc
    from aria.utils import privacy

    monkeypatch.setattr(privacy, "egress_allowed", lambda: False)
    value, reason = chromatin_qc._compute_tss_enrichment("frags.tsv.gz", "hg38")
    assert value is None
    assert reason and "air-gapped" in reason.lower()


def test_tss_enrichment_missing_snapatac2_returns_reason(monkeypatch):
    from aria.scripts import chromatin_qc
    from aria.utils import privacy

    monkeypatch.setattr(privacy, "egress_allowed", lambda: True)
    # Ensure importing snapatac2 fails deterministically.
    monkeypatch.setitem(sys.modules, "snapatac2", None)
    value, reason = chromatin_qc._compute_tss_enrichment("frags.tsv.gz", "hg38")
    assert value is None
    assert reason and "snapatac2" in reason.lower()


def test_tss_enrichment_computes_median_with_fake_snapatac2(monkeypatch):
    """With snapatac2 present + egress allowed, the sample-level value is the
    MEDIAN of per-cell TSSe (not a fabricated constant)."""
    from aria.scripts import chromatin_qc
    from aria.utils import privacy

    monkeypatch.setattr(privacy, "egress_allowed", lambda: True)

    class _Obs(dict):
        pass

    fake_adata = types.SimpleNamespace(obs={"tsse": [3.0, 9.0, 6.0]})
    fake_snap = types.SimpleNamespace(
        genome=types.SimpleNamespace(hg38=object()),
        pp=types.SimpleNamespace(import_data=lambda *a, **k: fake_adata),
        metrics=types.SimpleNamespace(tsse=lambda adata, gobj: None),
    )
    monkeypatch.setitem(sys.modules, "snapatac2", fake_snap)
    value, reason = chromatin_qc._compute_tss_enrichment("frags.tsv.gz", "hg38")
    assert reason is None
    assert value == pytest.approx(6.0)   # median of [3, 9, 6]


# ── FRiP — real overlap computation ──────────────────────────────────────────

def test_frip_from_intervals_counts_overlaps():
    from aria.scripts.chromatin_qc import _frip_from_intervals

    peaks = {"chr1": [(100, 200), (500, 600)]}
    fragments = [
        ("chr1", 150, 160),   # inside peak 1 -> in
        ("chr1", 50, 60),     # before peak 1 -> out
        ("chr1", 550, 560),   # inside peak 2 -> in
        ("chr2", 150, 160),   # no peaks on chr2 -> out
        ("chr1", 195, 260),   # straddles peak 1 end -> in
    ]
    frip, n = _frip_from_intervals(fragments, peaks)
    assert n == 5
    assert frip == pytest.approx(3 / 5)


def test_frip_from_intervals_empty_is_none():
    from aria.scripts.chromatin_qc import _frip_from_intervals

    frip, n = _frip_from_intervals([], {"chr1": [(1, 2)]})
    assert frip is None and n == 0


def test_estimate_frip_none_without_peaks():
    from aria.scripts.chromatin_qc import _estimate_frip

    # FRiP is undefined before peak calling -> honest None, never fabricated.
    assert _estimate_frip("frags.tsv.gz") is None
    assert _estimate_frip("frags.tsv.gz", peaks_bed=None) is None


def test_estimate_frip_real_with_peaks_bed(tmp_path):
    from aria.scripts.chromatin_qc import _estimate_frip

    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t100\t200\nchr1\t500\t600\n")
    frags = tmp_path / "fragments.tsv"
    frags.write_text(
        "chr1\t150\t160\tBC1\t1\n"   # in
        "chr1\t50\t60\tBC1\t1\n"     # out
        "chr1\t550\t560\tBC2\t1\n"   # in
    )
    frip = _estimate_frip(str(frags), peaks_bed=str(bed))
    assert frip == pytest.approx(2 / 3)


# ── P4.1 (2): per-cell TSSe×depth arrays + per-barcode FRiP distribution ──────

def test_compute_tss_enrichment_still_returns_2tuple(monkeypatch):
    """P4.1 [A] is additive: _compute_tss_enrichment keeps the fixed (value,
    reason) contract on the success path (now a thin wrapper over _tss_qc)."""
    from aria.scripts import chromatin_qc
    from aria.utils import privacy

    monkeypatch.setattr(privacy, "egress_allowed", lambda: True)
    fake_adata = types.SimpleNamespace(obs={"tsse": [3.0, 9.0, 6.0]})
    fake_snap = types.SimpleNamespace(
        genome=types.SimpleNamespace(hg38=object()),
        pp=types.SimpleNamespace(import_data=lambda *a, **k: fake_adata),
        metrics=types.SimpleNamespace(tsse=lambda adata, gobj: None),
    )
    monkeypatch.setitem(sys.modules, "snapatac2", fake_snap)
    out = chromatin_qc._compute_tss_enrichment("frags.tsv.gz", "hg38")
    assert isinstance(out, tuple) and len(out) == 2
    value, reason = out
    assert reason is None and value == pytest.approx(6.0)


def test_tss_qc_emits_per_cell_arrays_with_depth(monkeypatch):
    """_tss_qc surfaces per-cell (TSSe, log10 depth) arrays from the SAME pass;
    the median is unchanged. Depth comes from the importer's per-cell count."""
    from aria.scripts import chromatin_qc
    from aria.utils import privacy

    monkeypatch.setattr(privacy, "egress_allowed", lambda: True)
    fake_adata = types.SimpleNamespace(
        obs={"tsse": [3.0, 9.0, 6.0], "n_fragment": [1000.0, 10000.0, 5000.0]})
    fake_snap = types.SimpleNamespace(
        genome=types.SimpleNamespace(hg38=object()),
        pp=types.SimpleNamespace(import_data=lambda *a, **k: fake_adata),
        metrics=types.SimpleNamespace(tsse=lambda adata, gobj: None),
    )
    monkeypatch.setitem(sys.modules, "snapatac2", fake_snap)
    r = chromatin_qc._tss_qc("frags.tsv.gz", "hg38")
    assert r["median"] == pytest.approx(6.0)
    assert r["tsse"] == [3.0, 9.0, 6.0]
    assert r["log10_depth"] == pytest.approx([3.0, 4.0, pytest.approx(3.69897, abs=1e-4)])


def test_tss_qc_arrays_none_on_skip():
    """Unknown assembly -> honest None arrays + reason, never fabricated."""
    from aria.scripts.chromatin_qc import _tss_qc

    r = _tss_qc("frags.tsv.gz", "ce11")
    assert r["median"] is None and r["tsse"] is None and r["log10_depth"] is None
    assert r["reason"] and "assembly" in r["reason"].lower()


def test_frip_per_barcode_distribution_counts_per_cell():
    from aria.scripts.chromatin_qc import _frip_per_barcode_from_intervals

    peaks = {"chr1": [(100, 200), (500, 600)]}
    frags = [
        ("chr1", 150, 160, "BC1"),   # in
        ("chr1", 50, 60, "BC1"),     # out  -> BC1 FRiP 1/2
        ("chr1", 550, 560, "BC2"),   # in
        ("chr1", 520, 530, "BC2"),   # in   -> BC2 FRiP 2/2
    ]
    dist = _frip_per_barcode_from_intervals(frags, peaks, min_frags=2)
    assert sorted(dist) == pytest.approx([0.5, 1.0])
    # min_frags filters low-count barcodes.
    assert _frip_per_barcode_from_intervals(frags, peaks, min_frags=3) == []


def test_estimate_frip_distribution_none_without_peaks(tmp_path):
    from aria.scripts.chromatin_qc import _estimate_frip_distribution

    frags = tmp_path / "fragments.tsv"
    frags.write_text("chr1\t150\t160\tBC1\t1\n")
    assert _estimate_frip_distribution(str(frags), peaks_bed=None) is None


def test_estimate_frip_distribution_real_with_peaks(tmp_path):
    from aria.scripts.chromatin_qc import _estimate_frip_distribution

    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t100\t200\nchr1\t500\t600\n")
    frags = tmp_path / "fragments.tsv"
    frags.write_text(
        "chr1\t150\t160\tBC1\t1\n"   # in   -> BC1 1/2
        "chr1\t50\t60\tBC1\t1\n"     # out
        "chr1\t550\t560\tBC2\t1\n"   # in   -> BC2 2/2
        "chr1\t520\t530\tBC2\t1\n"   # in
    )
    dist = _estimate_frip_distribution(str(frags), peaks_bed=str(bed), min_frags=2)
    assert sorted(dist) == pytest.approx([0.5, 1.0])
