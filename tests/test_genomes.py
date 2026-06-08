"""v4.6 — reference-genome resolver (aria/utils/genomes.py).

ARIA must resolve a reference genome automatically from the inferred assembly,
not by asking the user for an env var. These guards cover the LIGHT resolution
policy (env + managed store + assembly→snapatac2 attribute); the governed
auto-fetch path is exercised in test_chromatin_motifs (snapatac2-gated).
"""

import pytest

from aria.utils import genomes


def test_genome_dir_respects_env(tmp_path, monkeypatch):
    monkeypatch.setenv("ARIA_GENOME_DIR", str(tmp_path / "g"))
    assert genomes.genome_dir() == tmp_path / "g"


def test_env_fallback_resolves_existing_file(tmp_path, monkeypatch):
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\nACGT\n")
    monkeypatch.setenv("ARIA_GENOME_FASTA", str(fa))
    assert genomes.genome_fasta_from_env() == str(fa)
    path, src = genomes.resolve_local_genome_fasta("hg38")
    assert path == str(fa) and src == "env:ARIA_GENOME_FASTA"


def test_managed_store_resolves_by_assembly_and_alias(tmp_path, monkeypatch):
    monkeypatch.delenv("ARIA_GENOME_FASTA", raising=False)
    monkeypatch.setenv("ARIA_GENOME_DIR", str(tmp_path))
    # stage under the GENCODE name; a run asking for hg38 must still find it
    sub = tmp_path / "grch38"
    sub.mkdir()
    fa = sub / "GRCh38.primary_assembly.genome.fa"
    fa.write_text(">chr1\nACGT\n")
    path, src = genomes.resolve_local_genome_fasta("hg38")
    assert path == str(fa) and src == "managed:ARIA_GENOME_DIR"


def test_unknown_assembly_resolves_to_nothing(tmp_path, monkeypatch):
    monkeypatch.delenv("ARIA_GENOME_FASTA", raising=False)
    monkeypatch.setenv("ARIA_GENOME_DIR", str(tmp_path))
    path, src = genomes.resolve_local_genome_fasta("not_an_assembly")
    assert path is None and src is None


def test_snapatac2_attr_maps_known_assemblies():
    assert genomes.snapatac2_attr("hg38") == "hg38"
    assert genomes.snapatac2_attr("GRCh38") == "GRCh38"   # case-insensitive
    assert genomes.snapatac2_attr("mm10") == "mm10"
    assert genomes.snapatac2_attr("hg19") == "hg19"
    assert genomes.snapatac2_attr("unknown") is None


def test_motifs_reexports_genome_env_helper_for_back_compat():
    # aria.utils.motifs.genome_fasta_from_env must still import (re-exported).
    from aria.utils.motifs import genome_fasta_from_env, GENOME_FASTA_ENV
    assert genome_fasta_from_env is genomes.genome_fasta_from_env
    assert GENOME_FASTA_ENV == "ARIA_GENOME_FASTA"
