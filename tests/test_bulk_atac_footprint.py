"""B4 (bulk ATAC pre-print): condition-level TOBIAS footprinting.

The bulk entry builds the two pseudobulk BAMs by MERGING each condition's per-sample
BAMs (vs the scATAC fragments-by-barcode split) and reuses the same TOBIAS core. TOBIAS
lives in a dedicated env, so a missing tool/asset is an honest skip — footprinting never
reports uncorrected output, and a missing condition BAM is never fabricated."""

from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest


class _FakeAln:
    def __init__(self, mapped):
        self.mapped = mapped


def _install_fake_pysam(monkeypatch, mapped_by_bam):
    """Install a fake pysam whose merge writes a stub file and whose AlignmentFile
    reports a deterministic mapped count per merged path."""
    merged_counts = {}

    def merge(flag, out, *inputs):
        Path(out).write_bytes(b"BAM")
        merged_counts[str(out)] = sum(mapped_by_bam.get(Path(i).name, 0)
                                      for i in inputs)

    def index(path):
        Path(str(path) + ".bai").write_bytes(b"BAI")

    def alignment_file(path):
        return _FakeAln(merged_counts.get(str(path), 0))

    fake = types.SimpleNamespace(merge=merge, index=index,
                                 AlignmentFile=alignment_file)
    monkeypatch.setitem(sys.modules, "pysam", fake)


def test_merge_condition_bams_merges_and_counts(tmp_path, monkeypatch):
    from aria.scripts.chromatin_footprint_tobias import _merge_condition_bams

    for name in ("a1.bam", "a2.bam", "b1.bam"):
        (tmp_path / name).write_bytes(b"BAM")
    _install_fake_pysam(monkeypatch, {"a1.bam": 100, "a2.bam": 50, "b1.bam": 70})

    out = _merge_condition_bams(
        {"CondA": [str(tmp_path / "a1.bam"), str(tmp_path / "a2.bam")],
         "CondB": [str(tmp_path / "b1.bam")]},
        tmp_path)
    assert out["CondA"]["n_fragments"] == 150 and out["CondA"]["n_bams"] == 2
    assert out["CondB"]["n_fragments"] == 70 and out["CondB"]["n_bams"] == 1
    assert Path(out["CondA"]["bam"]).exists()


def test_merge_condition_bams_honest_when_no_valid_bam(tmp_path, monkeypatch):
    from aria.scripts.chromatin_footprint_tobias import _merge_condition_bams

    _install_fake_pysam(monkeypatch, {})
    out = _merge_condition_bams({"CondA": ["/nope/missing.bam"]}, tmp_path)
    assert out["CondA"]["bam"] is None and out["CondA"]["n_fragments"] == 0


def test_bulk_footprint_skips_without_tobias(tmp_path, monkeypatch):
    import aria.scripts.chromatin_footprint_tobias as mod

    monkeypatch.setattr(mod.shutil, "which", lambda _x: None)
    res = mod.chromatin_footprint_tobias_bulk({
        "condition_bams": {"A": ["a.bam"], "B": ["b.bam"]},
        "genome_fasta": str(tmp_path / "g.fa"), "group_a": "A", "group_b": "B",
        "output_dir": str(tmp_path)})
    assert res["ran"] is False
    assert "TOBIAS not installed" in res["reason"]


def test_bulk_footprint_skips_on_missing_assets(tmp_path, monkeypatch):
    import aria.scripts.chromatin_footprint_tobias as mod

    monkeypatch.setattr(mod.shutil, "which", lambda _x: "/usr/bin/TOBIAS")
    # genome FASTA does not exist -> honest skip (no uncorrected output).
    res = mod.chromatin_footprint_tobias_bulk({
        "condition_bams": {"A": ["a.bam"], "B": ["b.bam"]},
        "genome_fasta": str(tmp_path / "missing.fa"),
        "peaks_bed": str(tmp_path / "p.bed"), "motif_meme": str(tmp_path / "m.meme"),
        "group_a": "A", "group_b": "B", "output_dir": str(tmp_path)})
    assert res["ran"] is False
    assert "genome_fasta" in res["reason"]


def test_bulk_footprint_skips_when_group_missing(tmp_path, monkeypatch):
    import aria.scripts.chromatin_footprint_tobias as mod

    for n in ("g.fa", "g.fa.fai", "p.bed", "m.meme"):
        (tmp_path / n).write_bytes(b"x")
    monkeypatch.setattr(mod.shutil, "which", lambda _x: "/usr/bin/TOBIAS")
    res = mod.chromatin_footprint_tobias_bulk({
        "condition_bams": {"A": ["a.bam"]},  # B absent
        "genome_fasta": str(tmp_path / "g.fa"),
        "peaks_bed": str(tmp_path / "p.bed"), "motif_meme": str(tmp_path / "m.meme"),
        "group_a": "A", "group_b": "B", "output_dir": str(tmp_path)})
    assert res["ran"] is False
    assert "condition_bams must contain" in res["reason"]


def test_bulk_footprint_narrator_block_uses_condition_labels():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    footprinting = {
        "ran": True, "data_type": "bulk_ATAC",
        "group_a": "K562", "group_b": "GM12878",
        "group_label": "Conditions", "group_kind": "conditions",
        "differential_summary": {"parsed": True, "n_ranked_candidates": 12,
                                 "n_motifs_tested": 800,
                                 "ranking_basis": {"fdr_controlled": False},
                                 "top_group_a": [], "top_group_b": []},
    }
    blocks = ChromatinNarrator().collect(
        "chromatin_agent", {"findings": {"footprinting": footprinting}})
    block = next(b for b in blocks
                 if b.id == "chromatin.differential_tf_footprinting")
    # Condition framing, not cell-type, and counts-only claim (W-CLAIM safe).
    assert any(e.label == "Conditions" for e in block.evidence)
    assert "K562 and GM12878" in block.claim
    assert any("between conditions" in c.text for c in block.caveats)
    assert "regulates" not in block.claim.lower()
    # B7 interim: the block must NOT present pseudobulk footprints as significance.
    assert "not FDR-controlled significance" in block.claim
    assert not any(e.label.startswith("Significant") for e in block.evidence)
