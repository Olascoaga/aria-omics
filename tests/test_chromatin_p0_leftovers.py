"""scATAC P0 minor leftovers: FRiP wired through chromatin_qc (from a called-peak
BED) and the per-cell accessibility depth surfaced for the UMAP overlay.
"""

from __future__ import annotations

import inspect

import pytest


# ── FRiP wiring: chromatin_qc forwards a peaks BED to the FRiP computation ────

def test_chromatin_qc_forwards_peaks_bed_to_scatac_qc(tmp_path, monkeypatch):
    from aria.scripts import chromatin_qc

    frag = tmp_path / "fragments.tsv"
    frag.write_text("chr1\t150\t160\tBC1\t1\n")
    bed = tmp_path / "peaks.narrowPeak"
    bed.write_text("chr1\t100\t200\tp1\t0\t.\n")

    captured = {}

    def fake_scatac(files, genome, organism, warnings, allow_mocks=False,
                    peaks_bed=None):
        captured["peaks_bed"] = peaks_bed
        return {"status": "success", "frip": None}

    monkeypatch.setattr(chromatin_qc, "_scatac_qc", fake_scatac)
    chromatin_qc.chromatin_qc({
        "data_type": "scATAC",
        "files": [str(frag)],
        "peaks_bed": str(bed),
    })
    assert captured["peaks_bed"] == str(bed)


def test_chromatin_qc_accepts_peaks_path_alias(tmp_path, monkeypatch):
    from aria.scripts import chromatin_qc

    frag = tmp_path / "fragments.tsv"
    frag.write_text("chr1\t150\t160\tBC1\t1\n")
    captured = {}

    def fake_scatac(files, genome, organism, warnings, allow_mocks=False,
                    peaks_bed=None):
        captured["peaks_bed"] = peaks_bed
        return {"status": "success"}

    monkeypatch.setattr(chromatin_qc, "_scatac_qc", fake_scatac)
    # chromatin_peaks emits `peaks_path`; chromatin_qc accepts it as an alias.
    chromatin_qc.chromatin_qc({
        "data_type": "scATAC",
        "files": [str(frag)],
        "peaks_path": "/called/peaks.narrowPeak",
    })
    assert captured["peaks_bed"] == "/called/peaks.narrowPeak"


def test_estimate_frip_reads_narrowpeak(tmp_path):
    """The FRiP BED reader handles narrowPeak (extra columns) — cols 1-3 are
    chrom/start/end, the rest are ignored."""
    from aria.scripts.chromatin_qc import _estimate_frip

    bed = tmp_path / "peaks.narrowPeak"
    bed.write_text("chr1\t100\t200\tpeak1\t500\t.\t7.1\t-1\t-1\t50\n")
    frags = tmp_path / "fragments.tsv"
    frags.write_text("chr1\t150\t160\tBC1\t1\nchr1\t10\t20\tBC1\t1\n")
    frip = _estimate_frip(str(frags), peaks_bed=str(bed))
    assert frip == pytest.approx(0.5)


# ── UMAP depth overlay: lsi writes per-cell depth; color keys include it ──────

def test_umap_color_keys_include_depth():
    from aria.agents._narrative_chromatin import _umap_color_keys

    keys = _umap_color_keys({"differential_accessibility": {"groupby": "leiden"}})
    assert "leiden" in keys
    assert "log10_n_fragments" in keys   # per-cell depth overlay


def test_lsi_writes_per_cell_depth_to_obs():
    """lsi must persist the per-cell accessibility depth to obs so the UMAP can be
    coloured by it (W0.2 leftover)."""
    from aria.scripts import chromatin_lsi_clustering as mod

    src = inspect.getsource(mod)
    assert 'adata.obs["n_fragments"]' in src
    assert 'adata.obs["log10_n_fragments"]' in src
