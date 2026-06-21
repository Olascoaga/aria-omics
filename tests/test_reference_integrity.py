"""S13 reference checksum integrity guards."""

from __future__ import annotations

import hashlib
import json

from aria.agents.modality_audit import build_capability_matrix
from aria.utils import genomes, motifs, ora
from aria.utils.reference_integrity import verify_reference_file


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def test_verify_reference_file_matches_declared_sha256(tmp_path):
    ref = tmp_path / "ref.gmt"
    text = "TERM\t\tGENE1\n"
    ref.write_text(text, encoding="utf-8")
    (tmp_path / "manifest.json").write_text(
        json.dumps({"sha256": _sha(text), "release": "test"}),
        encoding="utf-8",
    )

    result = verify_reference_file(ref)

    assert result["status"] == "ok"
    assert result["expected_sha256"] == _sha(text)
    assert result["observed_sha256"] == _sha(text)


def test_verify_reference_file_detects_declared_sha256_mismatch(tmp_path):
    ref = tmp_path / "ref.gmt"
    ref.write_text("TERM\t\tGENE1\n", encoding="utf-8")
    (tmp_path / "manifest.json").write_text(
        json.dumps({"sha256": "0" * 64}),
        encoding="utf-8",
    )

    result = verify_reference_file(ref)

    assert result["status"] == "checksum_mismatch"
    assert result["expected_sha256"] == "0" * 64
    assert result["observed_sha256"] != result["expected_sha256"]


def test_gmt_loader_rejects_checksum_mismatch(tmp_path, monkeypatch):
    monkeypatch.setenv(ora.GMT_DIR_ENV, str(tmp_path))
    base = tmp_path / "LIB"
    base.mkdir()
    (base / "LIB.gmt").write_text("TERM\t\tGENE1\n", encoding="utf-8")
    (base / "manifest.json").write_text(
        json.dumps({"library": "LIB", "sha256": "0" * 64}),
        encoding="utf-8",
    )

    assert ora.load_local_library("LIB") is None


def test_motif_loader_rejects_checksum_mismatch(tmp_path, monkeypatch):
    monkeypatch.setenv(motifs.MOTIF_DIR_ENV, str(tmp_path))
    base = tmp_path / "JASPAR_TEST"
    base.mkdir()
    (base / "JASPAR_TEST.meme").write_text(
        "MEME version 4\n\nMOTIF M1\nletter-probability matrix: alength= 4 w= 1\n"
        "0.25 0.25 0.25 0.25\n",
        encoding="utf-8",
    )
    (base / "manifest.json").write_text(
        json.dumps({"collection": "JASPAR_TEST", "sha256": "0" * 64}),
        encoding="utf-8",
    )

    assert motifs.load_local_motif_collection("JASPAR_TEST") is None


def test_managed_genome_rejects_checksum_mismatch(tmp_path, monkeypatch):
    monkeypatch.delenv(genomes.GENOME_FASTA_ENV, raising=False)
    monkeypatch.setenv(genomes.GENOME_DIR_ENV, str(tmp_path))
    base = tmp_path / "hg38"
    base.mkdir()
    fasta = base / "hg38.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    (base / "manifest.json").write_text(
        json.dumps({"sha256": "0" * 64}),
        encoding="utf-8",
    )

    path, source = genomes.resolve_local_genome_fasta("hg38")
    path_i, _source_i, integrity = genomes.resolve_local_genome_fasta_with_integrity("hg38")

    assert path is None and source is None
    assert path_i is None
    assert integrity["status"] == "checksum_mismatch"


def test_capability_matrix_surfaces_reference_checksum_mismatch(tmp_path, monkeypatch):
    monkeypatch.setenv(ora.GMT_DIR_ENV, str(tmp_path))
    base = tmp_path / "LIB"
    base.mkdir()
    (base / "LIB.gmt").write_text("TERM\t\tGENE1\n", encoding="utf-8")
    (base / "manifest.json").write_text(
        json.dumps({"library": "LIB", "sha256": "0" * 64}),
        encoding="utf-8",
    )

    matrix = build_capability_matrix({"modalities": {}})

    ref = matrix["preflight"]["reference_integrity"]
    assert ref["status"] == "red"
    assert ref["checks"][0]["status"] == "checksum_mismatch"
    assert any(
        f["check"] == "reference_checksum_mismatch"
        and f["severity"] == "blocking"
        for f in matrix["findings"]
    )
