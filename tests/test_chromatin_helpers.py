"""S15 guards for pure chromatin helper extraction."""

from __future__ import annotations


def test_chromatin_agent_reexports_pure_helpers_for_compatibility():
    from aria.agents import chromatin_agent
    from aria.agents import chromatin_helpers

    assert chromatin_agent._bulk_da_motif_regions is chromatin_helpers.bulk_da_motif_regions
    assert chromatin_agent._rank_bulk_da_motif_peaks is chromatin_helpers.rank_bulk_da_motif_peaks
    assert chromatin_agent._is_fastq is chromatin_helpers.is_fastq
    assert chromatin_agent._pick_read is chromatin_helpers.pick_read
    assert chromatin_agent._positive_int is chromatin_helpers.positive_int


def test_fastq_helpers_remain_case_insensitive():
    from aria.agents.chromatin_helpers import is_fastq, pick_read, positive_int

    files = ["sample_R1.FASTQ.GZ", "sample_R2.fastq.gz"]

    assert is_fastq(files) is True
    assert pick_read(files, ("_r1",)) == "sample_R1.FASTQ.GZ"
    assert positive_int("4", 1) == 4
    assert positive_int("0", 7) == 7
