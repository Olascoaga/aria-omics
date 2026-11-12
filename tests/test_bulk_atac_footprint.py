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

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def has_index(self):
        return True


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


@pytest.mark.parametrize("suffix,reference,expected_reference", [
    (".cram", "genome.fa", True), (".CRAM", "genome.fa", True),
    (".bam", "genome.fa", False), (".cram", None, False),
])
def test_merge_cram_forwards_reference_only_when_needed(
    tmp_path, monkeypatch, suffix, reference, expected_reference,
):
    from aria.scripts import chromatin_footprint_tobias as mod

    source = tmp_path / f"sample{suffix}"
    source.write_bytes(b"alignment")
    calls = []

    def merge(*args):
        calls.append(args)

    monkeypatch.setitem(sys.modules, "pysam", types.SimpleNamespace(
        merge=merge, index=lambda _: None))
    monkeypatch.setattr(mod, "_count_mapped", lambda *args: 123)
    result = mod._merge_condition_bams({"A": [str(source)]}, tmp_path,
                                      reuse=False, reference_fasta=reference)
    expected = ["-f"] + (["--reference", reference] if expected_reference else [])
    assert calls == [(*expected, str(tmp_path / "A.bam"), str(source))]
    assert result["A"]["n_fragments"] == 123


@pytest.mark.parametrize("indexed", [True, False])
def test_count_mapped_cram_opens_reference_and_counts_primary_only(tmp_path, monkeypatch, indexed):
    from aria.scripts.chromatin_footprint_tobias import _count_mapped

    calls = []
    reads = [types.SimpleNamespace(is_unmapped=u, is_secondary=s, is_supplementary=p)
             for u, s, p in [(False, False, False), (False, False, False),
                             (True, False, False), (False, True, False),
                             (False, False, True)]]

    class Alignment(_FakeAln):
        def has_index(self):
            return indexed

        def fetch(self, *, until_eof):
            assert until_eof is True
            return iter(reads)

    def open_alignment(path, **kwargs):
        calls.append((path, kwargs))
        return Alignment(7)

    monkeypatch.setitem(sys.modules, "pysam", types.SimpleNamespace(AlignmentFile=open_alignment))
    path = tmp_path / "input.CRAM"
    assert _count_mapped(path, "genome.fa") == (7 if indexed else 2)
    assert calls == [(str(path), {"reference_filename": "genome.fa"})]


def test_prepare_bulk_replicates_forwards_reference_and_reuse(tmp_path, monkeypatch):
    from aria.scripts import chromatin_footprint_tobias as mod

    calls = []

    def merge(inputs, out_dir, **kwargs):
        calls.append((inputs, out_dir, kwargs))
        return {}

    monkeypatch.setattr(mod, "_merge_condition_bams", merge)
    _, groups, ids = mod._prepare_bulk_replicate_bams(
        {"A": {"r1": ["a.cram", "a_lane2.cram"]}, "B": {"r2": "b.cram"}},
        "A", "B", tmp_path, reference_fasta="genome.fa", reuse=False)
    assert set(groups.values()) == {"A", "B"}
    assert set(ids.values()) == {"r1", "r2"}
    assert calls[0][1] == tmp_path / "replicate_bams"
    assert calls[0][2] == {"reference_fasta": "genome.fa", "reuse": False}
    assert list(calls[0][0].values()) == [["a.cram", "a_lane2.cram"], ["b.cram"]]
