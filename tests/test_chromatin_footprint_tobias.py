"""scATAC P4.3 — TOBIAS footprinting driver guards.

Pure helpers + the honest not-run gate are unit-tested in aria-env (no pysam, no
TOBIAS). The full ATACorrect->ScoreBigwig->BINDetect run is exercised separately in
the dedicated aria-tobias-env on real fragments. No fabrication: when TOBIAS or any
required asset is absent the driver returns ran=False with a concrete reason (ADR-002
/ W2.2: never an uncorrected footprint or an invented score).
"""

from __future__ import annotations

import pytest


def test_chrom_sizes_from_fai(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import _chrom_sizes_from_fai

    fai = tmp_path / "g.fa.fai"
    fai.write_text("chr1\t248956422\t112\t70\t71\nchr2\t242193529\t1\t70\t71\n",
                   encoding="utf-8")
    sizes = _chrom_sizes_from_fai(str(fai))
    assert sizes == {"chr1": 248956422, "chr2": 242193529}


def test_load_barcode_groups_skips_header(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import _load_barcode_groups

    tsv = tmp_path / "groups.tsv"
    tsv.write_text("barcode\tgroup\nBC1-1\tMonocyte\nBC2-1\tT_cell\nBC3-1\tMonocyte\n",
                   encoding="utf-8")
    groups = _load_barcode_groups(str(tsv))
    assert groups == {"BC1-1": "Monocyte", "BC2-1": "T_cell", "BC3-1": "Monocyte"}


def test_fragment_cut_reads_two_reads_at_tn5_sites():
    from aria.scripts.chromatin_footprint_tobias import _fragment_cut_reads

    reads = list(_fragment_cut_reads("chr1", 100, 200, "BC:0", read_len=50))
    assert len(reads) == 2
    fwd, rev = reads
    assert fwd["is_reverse"] is False and fwd["pos"] == 100
    assert rev["is_reverse"] is True and rev["pos"] == 150  # end - read_len
    # degenerate fragment (end <= start) yields nothing
    assert list(_fragment_cut_reads("chr1", 200, 200, "BC:1")) == []


def test_fragment_cut_reads_clamps_read_len_to_fragment():
    from aria.scripts.chromatin_footprint_tobias import _fragment_cut_reads

    reads = list(_fragment_cut_reads("chr1", 100, 130, "BC:0", read_len=50))
    # read_len clamped to the 30bp fragment so reads never exceed the fragment
    assert all(r["length"] == 30 for r in reads)
    assert reads[1]["pos"] == 100  # max(start, end - 30)


def test_summarize_bindetect_top_differential_per_group(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import summarize_bindetect

    res = tmp_path / "bindetect_results.txt"
    res.write_text(
        "name\ttotal_tfbs\tMonocyte_T_cell_change\tMonocyte_T_cell_pvalue\n"
        "CEBPB\t97\t1.10\t1e-100\n"        # strongly toward Monocyte
        "TCF7\t500\t-0.40\t1e-50\n"          # toward T_cell
        "NOISE\t5\t0.90\t1e-80\n"            # below min_sites -> excluded
        "FLAT\t300\t0.01\t0.4\n",            # not significant -> excluded
        encoding="utf-8")
    s = summarize_bindetect(str(res), "Monocyte", "T_cell")
    assert s["parsed"] is True
    assert s["n_motifs_tested"] == 4 and s["n_significant"] == 2
    assert s["top_toward_Monocyte"][0]["tf"] == "CEBPB"
    assert s["top_toward_T_cell"][0]["tf"] == "TCF7"


def test_summarize_bindetect_honest_when_columns_absent(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import summarize_bindetect

    res = tmp_path / "bad.txt"
    res.write_text("name\tother\nX\t1\n", encoding="utf-8")
    s = summarize_bindetect(str(res), "Monocyte", "T_cell")
    assert s["parsed"] is False and "missing column" in s["reason"]


def test_top_tfs_both_directions_deduped():
    from aria.scripts.chromatin_footprint_tobias import _top_tfs

    summary = {
        "top_toward_Monocyte": [{"tf": "CEBPB"}, {"tf": "CEBPA"}, {"tf": "SPIC"}],
        "top_toward_T_cell": [{"tf": "EGR1"}, {"tf": "CEBPB"}],  # CEBPB dup
    }
    tfs = _top_tfs(summary, "Monocyte", "T_cell", n=2)
    # top-2 each direction, deduped, order preserved (no biology dropped one-sided)
    assert tfs == ["CEBPB", "CEBPA", "EGR1"]


def test_aggregate_plots_skips_tfs_without_beds(tmp_path):
    from aria.scripts.chromatin_footprint_tobias import _aggregate_plots

    # No BINDetect bed dirs present -> honest empty result, no TOBIAS call needed.
    bindetect = tmp_path / "bindetect"
    bindetect.mkdir()
    out = _aggregate_plots(bindetect, {"A": "a.bw", "B": "b.bw"},
                           ["CEBPB", "EGR1"], tmp_path)
    assert out == {}


def test_driver_honest_not_run_without_tobias(tmp_path, monkeypatch):
    from aria.scripts import chromatin_footprint_tobias as mod

    # TOBIAS absent (the aria-env reality) -> honest ran=False, never a fake footprint.
    monkeypatch.setattr(mod.shutil, "which", lambda name: None)
    res = mod.chromatin_footprint_tobias({
        "fragments_file": str(tmp_path / "f.tsv.gz"),
        "genome_fasta": str(tmp_path / "g.fa"),
        "peaks_bed": str(tmp_path / "p.bed"),
        "motif_meme": str(tmp_path / "m.meme"),
        "barcode_groups": str(tmp_path / "g.tsv"),
        "group_a": "Monocyte", "group_b": "T_cell",
        "output_dir": str(tmp_path / "out"),
    })
    assert res["ran"] is False
    assert "TOBIAS not installed" in res["reason"]


def test_driver_honest_not_run_missing_asset(tmp_path, monkeypatch):
    from aria.scripts import chromatin_footprint_tobias as mod

    # TOBIAS present but a required asset is missing -> concrete reason, still no run.
    monkeypatch.setattr(mod.shutil, "which", lambda name: "/usr/bin/TOBIAS")
    res = mod.chromatin_footprint_tobias({
        "fragments_file": str(tmp_path / "absent_fragments.tsv.gz"),
        "genome_fasta": str(tmp_path / "g.fa"),
        "peaks_bed": str(tmp_path / "p.bed"),
        "motif_meme": str(tmp_path / "m.meme"),
        "barcode_groups": str(tmp_path / "g.tsv"),
        "group_a": "Monocyte", "group_b": "T_cell",
        "output_dir": str(tmp_path / "out"),
    })
    assert res["ran"] is False
    assert "fragments_file missing" in res["reason"]
