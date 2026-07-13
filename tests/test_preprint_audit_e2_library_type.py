"""Preprint-readiness audit E2: library type is declared, filenames are a hint.

`DataAuditAgent._classify_files` fell back to filename regex over `SIGNATURES` in
dict order, and the generic `bulk_RNA_raw` paired-end rule (`.*_R[12]_.*fastq`)
came first — so an ATAC/scATAC/ChIP/HiC FASTQ (which also has R1/R2/R3) was
silently classified as bulk RNA, even when its name carried a modality keyword.

After E2: an explicit declared library type (manifest / context) is authoritative;
a modality keyword in the filename beats the generic paired-end rule (a non-binding
hint); and a truly generic R1/R2 FASTQ with no signal is classified as
`bulk_RNA_raw` only as a NON-BINDING hint, flagged `ambiguous_library_type` for
CHECKPOINT-1 confirmation — never silently bound.
"""
from __future__ import annotations

from pathlib import Path

import pytest


def _agent():
    pytest.importorskip("litellm")
    from aria.agents.data_audit_agent import DataAuditAgent
    agent = object.__new__(DataAuditAgent)
    return agent


def _classify(agent, names, declared=None):
    return agent._classify_files([Path(n) for n in names], declared)


# ── Filename keyword beats the generic paired-end rule (order bug) ─────────────

@pytest.mark.parametrize("name,expected", [
    ("sample_atac_R1.fastq.gz", "bulk_ATAC"),
    ("lib_scatac_R1.fastq.gz", "scATAC"),
    ("H3K27ac_rep1_R1.fastq.gz", "ChIP"),
    ("hic_sampleA_R1.fastq.gz", "HiC"),
])
def test_modality_keyword_beats_generic_bulk_rna(name, expected):
    agent = _agent()
    classified = _classify(agent, [name])
    assert name in classified.get(expected, []), classified
    assert name not in classified.get("bulk_RNA_raw", [])


# ── Generic R1/R2 FASTQ → non-binding bulk RNA hint, flagged ambiguous ─────────

def test_generic_paired_end_is_flagged_ambiguous_not_silently_bound():
    agent = _agent()
    names = ["GSM123_R1.fastq.gz", "GSM123_R2.fastq.gz"]
    classified = _classify(agent, names)
    # The hint is still bulk RNA (least-surprise default)...
    assert set(names) <= set(classified.get("bulk_RNA_raw", []))
    # ...but it is recorded ambiguous, so CP1 confirms instead of silent dispatch.
    ambiguous = getattr(agent, "_ambiguous_library_types", [])
    assert set(names) <= set(ambiguous)
    assert agent._ambiguous_library_type_warnings()  # a CP1 warning is emitted


# ── An explicit declared library type is authoritative ────────────────────────

def test_global_declared_library_type_overrides_filename():
    agent = _agent()
    # A user asserts the whole run is scATAC; a generic R1/R2 name must not win.
    classified = _classify(
        agent, ["GSM123_R1.fastq.gz"], {"global": "scATAC"})
    assert classified.get("scATAC") == ["GSM123_R1.fastq.gz"]
    assert "bulk_RNA_raw" not in classified
    # Authoritatively declared → not ambiguous.
    assert not getattr(agent, "_ambiguous_library_types", [])


def test_per_file_declared_library_type_overrides_filename():
    agent = _agent()
    classified = _classify(
        agent,
        ["tumor_R1.fastq.gz", "tumor_R2.fastq.gz"],
        {"per_file": {"tumor_R1.fastq.gz": "bulk_ATAC",
                      "tumor_R2.fastq.gz": "bulk_ATAC"}},
    )
    assert set(classified.get("bulk_ATAC", [])) == {
        "tumor_R1.fastq.gz", "tumor_R2.fastq.gz"}
    assert "bulk_RNA_raw" not in classified


# ── Control: non-FASTQ classification is unchanged ────────────────────────────

def test_processed_counts_matrix_still_bulk_rna():
    agent = _agent()
    classified = _classify(agent, ["experiment_counts.tsv"])
    assert classified.get("bulk_RNA") == ["experiment_counts.tsv"]
    assert not getattr(agent, "_ambiguous_library_types", [])


def test_declared_atac_fastq_is_not_flagged_ambiguous():
    agent = _agent()
    classified = _classify(agent, ["sample_atac_R1.fastq.gz"])
    # A keyword-typed FASTQ is a confident hint, not ambiguous.
    assert not getattr(agent, "_ambiguous_library_types", [])
