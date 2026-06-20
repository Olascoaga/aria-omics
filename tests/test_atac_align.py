"""B2a (bulk ATAC raw entry): FASTQ -> filtered BAM alignment.

atac_align is the ATAC analogue of rna_align (STAR): bwa-mem2 + ENCODE-style
filtering. It is binary-gated — a missing aligner is an honest MissingDependency
skip, never a fabricated BAM (ADR-002) unless mocks are explicitly enabled. The
ChromatinAgent prepends it to the bulk ATAC lane only for raw FASTQ inputs;
scATAC FASTQ is intentionally not handled.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from aria.scripts import atac_align as aa
from aria.agents.chromatin_agent import _is_fastq


# ── pure-input validation ─────────────────────────────────────────────────

def test_atac_align_requires_samples():
    out = aa.atac_align({"genome_fasta": "/tmp/g.fa"})
    assert out["status"] == "error"
    assert out["error_type"] == "NoSamples"


def test_atac_align_requires_genome():
    out = aa.atac_align({"samples": [{"name": "s", "r1": "/tmp/s_R1.fq.gz"}]})
    assert out["status"] == "error"
    assert out["error_type"] == "MissingGenomeFasta"


def test_atac_align_missing_index_and_fasta():
    out = aa.atac_align({
        "samples": [{"name": "s", "r1": "/tmp/s_R1.fq.gz"}],
        "bwa_index": "/tmp/nonexistent_index",
    })
    assert out["status"] == "error"
    assert out["error_type"] == "MissingGenomeIndex"


# ── honest-skip when the aligner is absent ────────────────────────────────

def test_atac_align_honest_skip_when_bwa_missing(monkeypatch, tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")

    def _boom(*a, **k):
        raise FileNotFoundError("bwa-mem2")

    monkeypatch.setattr(subprocess, "run", _boom)
    out = aa.atac_align({
        "samples": [{"name": "s", "r1": str(tmp_path / "s_R1.fq.gz"),
                     "paired": False}],
        "genome_fasta": str(fasta),
        "output_dir": str(tmp_path / "out"),
    })
    assert out["status"] == "error"
    assert out["error_type"] == "MissingDependency"


def test_atac_align_mock_mode_when_explicit(monkeypatch, tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")

    def _boom(*a, **k):
        raise FileNotFoundError("bwa-mem2")

    monkeypatch.setattr(subprocess, "run", _boom)
    out = aa.atac_align({
        "samples": [{"name": "s", "r1": str(tmp_path / "s_R1.fq.gz")}],
        "genome_fasta": str(fasta),
        "output_dir": str(tmp_path / "out"),
        "allow_mock": True,
    })
    assert out["status"] == "success"
    assert out["bam_files"][0]["note"].startswith("mock")


# ── ATAC filter parameters are reported (provenance) ──────────────────────

def test_filter_constants():
    assert aa.SAMTOOLS_EXCLUDE_FLAGS == 1804
    assert aa.SAMTOOLS_REQUIRE_FLAGS == 2
    assert aa.DEFAULT_MAPQ == 30


# ── FASTQ detection used by the agent prepend ─────────────────────────────

def test_is_fastq_helper():
    assert _is_fastq(["/data/s_R1.fastq.gz"]) is True
    assert _is_fastq(["/data/s_R1.fq"]) is True
    assert _is_fastq(["/data/sample.bam"]) is False
    assert _is_fastq([]) is False


# ── agent prepend: FASTQ triggers alignment, then the validated lane ──────

class _FakeFastqEnv:
    """Returns trim -> align -> qc -> peaks success for a FASTQ bulk ATAC run."""

    def __init__(self):
        self.calls = []

    def run_in_stack(self, *, stack, script_path, params):
        self.calls.append((stack, Path(script_path).name))
        if script_path.endswith("rna_fastq_qc.py"):
            return {"status": "success", "n_samples": 1,
                    "samples": [{"name": "s", "r1_trimmed": "/tmp/s_R1.fq.gz",
                                 "paired": False}],
                    "multiqc_report": None}
        if script_path.endswith("atac_align.py"):
            return {"status": "success", "aligner": "bwa-mem2",
                    "filter": {"mapq": 30, "samtools_f": 2, "samtools_F": 1804,
                               "remove_mito": True},
                    "bam_files": [{"name": "s", "status": "success",
                                   "bam": "/tmp/aligned/s/s.filtered.bam"}]}
        if script_path.endswith("chromatin_qc.py"):
            return {"status": "success", "data_type": "bulk_ATAC",
                    "n_samples": 1, "warnings": []}
        if script_path.endswith("chromatin_peaks.py"):
            return {"status": "success", "data_type": "bulk_ATAC",
                    "n_peaks": 100, "genome": params.get("genome"),
                    "peaks_path": "/tmp/p.narrowPeak", "warnings": []}
        return {"status": "success"}


def _agent(env):
    from aria.agents.chromatin_agent import ChromatinAgent
    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = env
    agent.publish_status = lambda *a, **k: None
    agent.publish_finding = lambda *a, **k: None
    return agent


def test_bulk_atac_fastq_triggers_alignment(tmp_path):
    fq = tmp_path / "s_R1.fastq.gz"
    fq.write_bytes(b"")
    env = _FakeFastqEnv()
    agent = _agent(env)

    result = agent._run_bulk_atac(
        "exp",
        {"genome": "hg38", "genome_fasta": str(tmp_path / "genome.fa")},
        {},
        [str(fq)],
    )

    assert result["status"] == "done"
    pre = result["findings"]["preprocessing"]
    assert pre["ran"] is True and pre["aligner"] == "bwa-mem2"
    # The alignment ran before QC, and QC/peaks consumed the produced BAM(s).
    stacks = [s for s, _ in env.calls]
    scripts = [n for _, n in env.calls]
    assert scripts[:2] == ["rna_fastq_qc.py", "atac_align.py"]
    assert stacks[:2] == ["rnaseq", "atacseq"]
    assert "chromatin_qc.py" in scripts and "chromatin_peaks.py" in scripts


def test_bulk_atac_fastq_alignment_failure_skips_honestly(tmp_path):
    fq = tmp_path / "s_R1.fastq.gz"
    fq.write_bytes(b"")

    class _FailEnv(_FakeFastqEnv):
        def run_in_stack(self, *, stack, script_path, params):
            self.calls.append((stack, Path(script_path).name))
            if script_path.endswith("rna_fastq_qc.py"):
                return {"status": "success", "n_samples": 1,
                        "samples": [{"name": "s", "r1_trimmed": "x",
                                     "paired": False}]}
            if script_path.endswith("atac_align.py"):
                return {"status": "error", "error_type": "MissingDependency",
                        "details": "bwa-mem2 absent"}
            return {"status": "success"}

    env = _FailEnv()
    agent = _agent(env)
    result = agent._run_bulk_atac(
        "exp", {"genome": "hg38", "genome_fasta": str(tmp_path / "g.fa")},
        {}, [str(fq)])

    assert result["status"] == "failed"
    assert result["reason"] == "bulk_atac_preprocessing_failed"
    # The validated lane never ran on un-aligned FASTQ.
    assert "chromatin_peaks.py" not in [n for _, n in env.calls]


def test_bulk_atac_fastq_no_genome_skips(tmp_path, monkeypatch):
    fq = tmp_path / "s_R1.fastq.gz"
    fq.write_bytes(b"")
    # No genome_fasta in exp_ctx and none staged locally.
    from aria.utils import genomes
    monkeypatch.setattr(genomes, "resolve_local_genome_fasta",
                        lambda a: (None, None))
    env = _FakeFastqEnv()
    agent = _agent(env)
    result = agent._run_bulk_atac("exp", {"genome": "hg38"}, {}, [str(fq)])
    assert result["status"] == "failed"
    assert result["findings"]["preprocessing"]["reason"] == \
        "no_reference_genome_fasta"


# ── narrator methods discloses the alignment when it ran ──────────────────

def test_methods_reports_alignment_when_preprocessed():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    narr = ChromatinNarrator.__new__(ChromatinNarrator)
    findings = {
        "preprocessing": {"analysis": "bulk_atac_alignment", "ran": True,
                          "aligner": "bwa-mem2", "assembly": "hg38",
                          "n_samples_aligned": 2,
                          "filter": {"mapq": 30, "samtools_f": 2,
                                     "samtools_F": 1804, "remove_mito": True}},
        "qc": {"data_type": "bulk_ATAC", "status": "success"},
        "peaks": {"data_type": "bulk_ATAC", "status": "success",
                  "n_peaks": 100, "genome": "hg38"},
    }
    lines = narr._bulk_atac_methods(findings)
    joined = " ".join(lines)
    assert "bwa-mem2" in joined
    assert "markdup" in joined
    assert "MAPQ >= 30" in joined
    # The "does not realign" claim must NOT appear once ARIA aligned the reads.
    assert "does not realign" not in joined


def test_methods_keeps_no_realign_claim_without_preprocessing():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    narr = ChromatinNarrator.__new__(ChromatinNarrator)
    findings = {
        "qc": {"data_type": "bulk_ATAC", "status": "success"},
        "peaks": {"data_type": "bulk_ATAC", "status": "success",
                  "n_peaks": 100, "genome": "hg38"},
    }
    joined = " ".join(narr._bulk_atac_methods(findings))
    assert "does not realign" in joined
    assert "bwa-mem2" not in joined
