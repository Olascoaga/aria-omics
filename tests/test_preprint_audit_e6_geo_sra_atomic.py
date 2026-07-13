"""E6 guards for atomic, content-validated GEO/SRA retrieval."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import tarfile
import urllib.error
import zipfile
from pathlib import Path

import pytest


def _fastq(path: Path, name: str = "READ_1") -> None:
    path.write_text(f"@{name}\nACGT\n+\nIIII\n", encoding="utf-8")


def test_e6_atomic_download_verifies_checksum_and_preserves_prior_file(tmp_path):
    from aria.utils.atomic_retrieval import RetrievalError, download_atomic

    source = tmp_path / "source.tsv"
    source.write_text("gene\tsample\nGENE_1\t4\n", encoding="utf-8")
    destination = tmp_path / "published.tsv"
    destination.write_text("prior valid bytes\n", encoding="utf-8")

    expected = hashlib.md5(source.read_bytes()).hexdigest()
    record = download_atomic(source.as_uri(), destination, expected_md5=expected)
    assert destination.read_bytes() == source.read_bytes()
    assert record["md5"] == expected
    assert record["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()

    destination.write_text("prior valid bytes\n", encoding="utf-8")
    with pytest.raises(RetrievalError, match="MD5"):
        download_atomic(source.as_uri(), destination, expected_md5="0" * 32)
    assert destination.read_text(encoding="utf-8") == "prior valid bytes\n"
    assert not list(tmp_path.glob(".*.part-*"))


def test_e6_safe_archive_extraction_supports_tar_and_zip(tmp_path):
    from aria.utils.atomic_retrieval import safe_extract_archive

    payload = b"gene\tsample\nGENE_1\t9\n"
    tar_path = tmp_path / "bundle.tar.gz"
    with tarfile.open(tar_path, "w:gz") as archive:
        info = tarfile.TarInfo("nested/counts.tsv")
        info.size = len(payload)
        archive.addfile(info, io.BytesIO(payload))
    zip_path = tmp_path / "bundle.zip"
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("matrix/matrix.mtx", "%%MatrixMarket matrix coordinate integer general\n")

    tar_files = safe_extract_archive(tar_path, tmp_path / "tar_out")
    zip_files = safe_extract_archive(zip_path, tmp_path / "zip_out")

    assert [p.relative_to(tmp_path / "tar_out").as_posix() for p in tar_files] == [
        "nested/counts.tsv"
    ]
    assert [p.relative_to(tmp_path / "zip_out").as_posix() for p in zip_files] == [
        "matrix/matrix.mtx"
    ]


def test_e6_geo_supplement_archive_is_extracted_and_classified(tmp_path, monkeypatch):
    from aria.connectors.geo_connector import GEOConnector

    payload = b"gene\tsample\nGENE_1\t9\n"
    archive_bytes = io.BytesIO()
    with tarfile.open(fileobj=archive_bytes, mode="w:gz") as archive:
        info = tarfile.TarInfo("nested/counts.tsv")
        info.size = len(payload)
        archive.addfile(info, io.BytesIO(payload))
    html = b'<a href="GSE000001_processed.tar.gz">payload</a>'

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    def fake_urlopen(url, *args, **kwargs):
        return _Response(html if str(url).endswith("/suppl/") else archive_bytes.getvalue())

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    files = GEOConnector(cache_dir=str(tmp_path))._download_supplementary(
        "GSE000001", tmp_path / "stage", None
    )

    assert len(files["counts"]) == 1
    assert Path(files["counts"][0]).read_bytes() == payload


def test_e6_compressed_h5ad_is_published_as_analyzable_h5ad(
    tmp_path, monkeypatch
):
    from aria.connectors.geo_connector import GEOConnector

    hdf5 = b"\x89HDF\r\n\x1a\n" + b"public-h5ad-payload"
    compressed = io.BytesIO()
    with gzip.GzipFile(fileobj=compressed, mode="wb") as handle:
        handle.write(hdf5)
    html = b'<a href="GSE000001_processed.h5ad.gz">payload</a>'

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    def fake_urlopen(url, *args, **kwargs):
        return _Response(html if str(url).endswith("/suppl/") else compressed.getvalue())

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    files = GEOConnector(cache_dir=str(tmp_path))._download_supplementary(
        "GSE000001", tmp_path / "stage", None
    )

    assert len(files["h5ad"]) == 1
    published = Path(files["h5ad"][0])
    assert published.name == "GSE000001_processed.h5ad"
    assert published.read_bytes() == hdf5


@pytest.mark.parametrize("kind", ["tar", "zip"])
def test_e6_archive_extraction_rejects_path_traversal(tmp_path, kind):
    from aria.utils.atomic_retrieval import RetrievalError, safe_extract_archive

    archive_path = tmp_path / ("escape.tar" if kind == "tar" else "escape.zip")
    if kind == "tar":
        with tarfile.open(archive_path, "w") as archive:
            info = tarfile.TarInfo("../escape.tsv")
            info.size = 1
            archive.addfile(info, io.BytesIO(b"x"))
    else:
        with zipfile.ZipFile(archive_path, "w") as archive:
            archive.writestr("../escape.tsv", "x")

    with pytest.raises(RetrievalError, match="unsafe archive member"):
        safe_extract_archive(archive_path, tmp_path / "out")
    assert not (tmp_path / "escape.tsv").exists()


def test_e6_content_validation_rejects_corrupt_gzip(tmp_path):
    from aria.utils.atomic_retrieval import RetrievalError, validate_payload

    corrupt = tmp_path / "counts.tsv.gz"
    corrupt.write_bytes(b"not gzip")
    with pytest.raises(RetrievalError, match="gzip"):
        validate_payload(corrupt, "counts")


def test_e6_retrieval_manifest_detects_post_publish_tampering(tmp_path):
    from aria.utils.atomic_retrieval import (
        validate_retrieval_manifest,
        write_retrieval_manifest,
    )

    root = tmp_path / "published"
    root.mkdir()
    payload = root / "counts.tsv"
    payload.write_text("gene\tsample\nGENE_1\t1\n", encoding="utf-8")
    write_retrieval_manifest(root, accession="GSE000001", payloads=[payload])

    assert validate_retrieval_manifest(root)["status"] == "valid"
    payload.write_bytes(b"x" * payload.stat().st_size)
    result = validate_retrieval_manifest(root)
    assert result["status"] == "invalid"
    assert any("sha256" in error for error in result["errors"])


def test_e6_retrieval_manifest_rejects_unmanifested_payload(tmp_path):
    from aria.utils.atomic_retrieval import (
        validate_retrieval_manifest,
        write_retrieval_manifest,
    )

    root = tmp_path / "published"
    root.mkdir()
    payload = root / "counts.tsv"
    payload.write_text("gene\tsample\nGENE_1\t1\n", encoding="utf-8")
    write_retrieval_manifest(root, accession="GSE000001", payloads=[payload])
    (root / "unexpected.fastq").write_text("@R1\nACGT\n+\nIIII\n", encoding="utf-8")

    result = validate_retrieval_manifest(root)
    assert result["status"] == "invalid"
    assert any("unmanifested payload" in error for error in result["errors"])


def test_e6_retrieval_manifest_binds_the_requested_accession(tmp_path):
    from aria.utils.atomic_retrieval import (
        validate_retrieval_manifest,
        write_retrieval_manifest,
    )

    root = tmp_path / "published"
    root.mkdir()
    payload = root / "counts.tsv"
    payload.write_text("gene\tsample\nGENE_1\t1\n", encoding="utf-8")
    write_retrieval_manifest(root, accession="GSE000001", payloads=[payload])

    result = validate_retrieval_manifest(root, expected_accession="GSE000002")
    assert result["status"] == "invalid"
    assert any("accession mismatch" in error for error in result["errors"])


def test_e6_network_primitive_respects_air_gap(monkeypatch):
    from aria.utils.atomic_retrieval import open_url_with_retry
    from aria.utils.privacy import EgressBlocked

    called = []
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    monkeypatch.setattr(
        "urllib.request.urlopen", lambda *args, **kwargs: called.append(args)
    )

    with pytest.raises(EgressBlocked):
        open_url_with_retry("https://example.invalid/payload")
    assert called == []


def test_e6_content_validation_accepts_real_bam_and_cram_signatures(tmp_path):
    from aria.utils.atomic_retrieval import validate_payload

    bam = tmp_path / "reads.bam"
    with gzip.open(bam, "wb") as handle:
        handle.write(b"BAM\x01minimal-header")
    cram = tmp_path / "reads.cram"
    cram.write_bytes(b"CRAMminimal-header")

    assert validate_payload(bam, "bam")["kind"] == "bam"
    assert validate_payload(cram, "bam")["kind"] == "bam"


def test_e6_soft_parser_preserves_sra_relations():
    from aria.connectors.geo_connector import _parse_soft_text

    parsed = _parse_soft_text(
        "!Series_title = Public study\n"
        "!Series_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRP123456\n"
        "^SAMPLE = GSM000001\n"
        "!Sample_relation = BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA999999\n"
    )
    assert parsed["sra_accessions"] == ["SRP123456", "PRJNA999999"]


def test_e6_sra_runinfo_is_parsed_without_pysradb(tmp_path, monkeypatch):
    from aria.connectors.geo_connector import GEOConnector

    runinfo = (
        "Run,ScientificName,LibraryStrategy,LibraryLayout,SampleName,BioSample\n"
        "SRR000001,Homo sapiens,RNA-Seq,PAIRED,sample_a,SAMN1\n"
    ).encode()
    search = json.dumps({
        "esearchresult": {
            "count": "1", "querykey": "1", "webenv": "test-history"
        }
    }).encode()

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    responses = iter((_Response(search), _Response(runinfo)))
    monkeypatch.setattr("urllib.request.urlopen", lambda *args, **kwargs: next(responses))
    rows = GEOConnector(cache_dir=str(tmp_path))._fetch_sra_runinfo("SRP000001")

    assert rows[0]["Run"] == "SRR000001"
    assert rows[0]["LibraryStrategy"] == "RNA-Seq"


def test_e6_sra_run_accession_does_not_expand_to_sibling_runs(
    tmp_path, monkeypatch
):
    from aria.connectors.geo_connector import GEOConnector

    runinfo = (
        "Run,ScientificName,LibraryStrategy\n"
        "SRR000001,Homo sapiens,RNA-Seq\n"
        "SRR000002,Homo sapiens,RNA-Seq\n"
    ).encode()
    search = json.dumps({
        "esearchresult": {"count": "2", "querykey": "1", "webenv": "history"}
    }).encode()

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    responses = iter((_Response(search), _Response(runinfo)))
    monkeypatch.setattr("urllib.request.urlopen", lambda *args, **kwargs: next(responses))
    rows = GEOConnector(cache_dir=str(tmp_path))._fetch_sra_runinfo("SRR000001")

    assert [row["Run"] for row in rows] == ["SRR000001"]


def test_e6_sra_runinfo_retries_ncbi_rate_limit(tmp_path, monkeypatch):
    from aria.connectors.geo_connector import GEOConnector

    search = json.dumps({
        "esearchresult": {
            "count": "1", "querykey": "1", "webenv": "test-history"
        }
    }).encode()
    runinfo = (
        "Run,ScientificName,LibraryStrategy\n"
        "SRR000001,Homo sapiens,RNA-Seq\n"
    ).encode()

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    limited = urllib.error.HTTPError(
        "https://example.test", 429, "Too Many Requests", {}, None
    )
    responses = iter((_Response(search), limited, _Response(runinfo)))

    def fake_urlopen(*args, **kwargs):
        response = next(responses)
        if isinstance(response, Exception):
            raise response
        return response

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    monkeypatch.setattr("time.sleep", lambda *_: None)
    rows = GEOConnector(cache_dir=str(tmp_path))._fetch_sra_runinfo("SRP000001")

    assert rows[0]["Run"] == "SRR000001"


def test_e6_sra_fetch_downloads_every_run_and_preserves_mixed_modalities(
    tmp_path, monkeypatch
):
    from aria.connectors.geo_connector import GEOConnector

    rows = [
        {
            "Run": "SRR000001", "ScientificName": "Homo sapiens",
            "LibraryStrategy": "RNA-Seq", "LibraryLayout": "PAIRED",
            "SampleName": "rna_sample", "BioSample": "SAMN1",
        },
        {
            "Run": "SRR000002", "ScientificName": "Homo sapiens",
            "LibraryStrategy": "ATAC-seq", "LibraryLayout": "PAIRED",
            "SampleName": "atac_sample", "BioSample": "SAMN2",
        },
    ]
    connector = GEOConnector(cache_dir=str(tmp_path))
    monkeypatch.setattr(connector, "_fetch_sra_runinfo", lambda accession: rows)

    calls = []

    def fake_toolkit(row, destination, status_cb=None):
        calls.append(row["Run"])
        destination.mkdir(parents=True, exist_ok=True)
        outputs = []
        for mate in (1, 2):
            path = destination / f"{row['Run']}_{mate}.fastq"
            _fastq(path, f"{row['Run']}/{mate}")
            outputs.append(path)
        return outputs

    monkeypatch.setattr(connector, "_retrieve_sra_run", fake_toolkit)
    result = connector.fetch("SRP000001")

    assert calls == ["SRR000001", "SRR000002"]
    assert result["retrieval_status"] == "complete"
    assert result["data_type"] == "mixed"
    assert len(result["files"]["fastq"]) == 4
    assert {result["file_modalities"][path] for path in result["files"]["fastq"]} == {
        "bulk_RNA_raw", "bulk_ATAC"
    }
    assert "fastq_pending" not in result["files"]
    assert Path(result["retrieval_manifest"]).is_file()

    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    cached = connector.fetch("SRP000001")
    assert cached == result
    assert calls == ["SRR000001", "SRR000002"]


def test_e6_sra_partial_failure_never_publishes_accession(tmp_path, monkeypatch):
    from aria.connectors.geo_connector import GEOConnector
    from aria.utils.atomic_retrieval import RetrievalError

    connector = GEOConnector(cache_dir=str(tmp_path))
    rows = [
        {"Run": "SRR000001", "LibraryStrategy": "RNA-Seq"},
        {"Run": "SRR000002", "LibraryStrategy": "RNA-Seq"},
    ]
    monkeypatch.setattr(connector, "_fetch_sra_runinfo", lambda accession: rows)

    def partial(row, destination, status_cb=None):
        destination.mkdir(parents=True, exist_ok=True)
        if row["Run"] == "SRR000002":
            raise RetrievalError("toolkit failed")
        path = destination / "SRR000001.fastq"
        _fastq(path)
        return [path]

    monkeypatch.setattr(connector, "_retrieve_sra_run", partial)
    with pytest.raises(RetrievalError, match="toolkit failed"):
        connector.fetch("SRP000001")

    assert not (tmp_path / "SRP000001").exists()
    assert not list(tmp_path.glob(".SRP000001.*.staging"))


def test_e6_gse_falls_back_to_related_sra_when_no_processed_payload(
    tmp_path, monkeypatch
):
    from aria.connectors.geo_connector import GEOConnector

    connector = GEOConnector(cache_dir=str(tmp_path))
    metadata = {
        "title": "Study", "organism": "Homo sapiens", "samples": [],
        "suppl_files": [], "library_strategy": "RNA-Seq",
        "sra_accessions": ["SRP000001"],
    }
    monkeypatch.setattr(connector, "_parse_soft", lambda *args: metadata)
    monkeypatch.setattr(connector, "_download_supplementary", lambda *args: {
        "counts": [], "h5ad": [], "h5": [], "mtx": [], "fragments": [],
        "peaks": [], "bam": [], "fastq": [],
    })
    monkeypatch.setattr(connector, "_fetch_sra_runinfo", lambda accession: [{
        "Run": "SRR000001", "ScientificName": "Homo sapiens",
        "LibraryStrategy": "RNA-Seq", "LibraryLayout": "SINGLE",
        "SampleName": "sample_a", "BioSample": "SAMN1",
    }])

    def fake_toolkit(row, destination, status_cb=None):
        destination.mkdir(parents=True, exist_ok=True)
        path = destination / "SRR000001.fastq"
        _fastq(path)
        return [path]

    monkeypatch.setattr(connector, "_retrieve_sra_run", fake_toolkit)
    result = connector.fetch("GSE000001")

    assert result["source_accessions"] == ["GSE000001", "SRP000001"]
    assert result["data_type"] == "bulk_RNA"
    assert len(result["files"]["fastq"]) == 1
    assert result["file_modalities"][result["files"]["fastq"][0]] == "bulk_RNA_raw"


def test_e6_gse_fallback_unions_sample_level_sra_relations(tmp_path, monkeypatch):
    from aria.connectors.geo_connector import GEOConnector

    connector = GEOConnector(cache_dir=str(tmp_path))
    metadata = {
        "title": "Study", "organism": "Homo sapiens", "samples": [],
        "suppl_files": [], "library_strategy": "mixed",
        "sra_accessions": ["SRX000001", "SRX000002"],
    }
    monkeypatch.setattr(connector, "_parse_soft", lambda *args: metadata)
    monkeypatch.setattr(connector, "_download_supplementary", lambda *args: {
        "counts": [], "h5ad": [], "h5": [], "mtx": [], "fragments": [],
        "peaks": [], "bam": [], "fastq": [],
    })

    def fake_runinfo(accession):
        suffix = "1" if accession == "SRX000001" else "2"
        return [{
            "Run": f"SRR00000{suffix}", "ScientificName": "Homo sapiens",
            "LibraryStrategy": "RNA-Seq" if suffix == "1" else "ATAC-seq",
            "LibraryLayout": "SINGLE", "SampleName": f"sample_{suffix}",
            "BioSample": f"SAMN{suffix}",
        }]

    monkeypatch.setattr(connector, "_fetch_sra_runinfo", fake_runinfo)

    def fake_toolkit(row, destination, status_cb=None):
        destination.mkdir(parents=True, exist_ok=True)
        path = destination / f"{row['Run']}.fastq"
        _fastq(path, row["Run"])
        return [path]

    monkeypatch.setattr(connector, "_retrieve_sra_run", fake_toolkit)
    result = connector.fetch("GSE000001")

    assert result["source_accessions"] == [
        "GSE000001", "SRX000001", "SRX000002"
    ]
    assert {Path(path).name for path in result["files"]["fastq"]} == {
        "SRR000001.fastq", "SRR000002.fastq"
    }
    assert result["data_type"] == "mixed"


def test_e6_data_audit_uses_authoritative_per_fastq_geo_modalities(tmp_path):
    from aria.agents.data_audit_agent import DataAuditAgent
    from aria.utils.atomic_retrieval import write_retrieval_manifest

    rna = tmp_path / "SRR000001_1.fastq"
    atac = tmp_path / "SRR000002_1.fastq"
    _fastq(rna)
    _fastq(atac)
    manifest = write_retrieval_manifest(
        tmp_path, accession="SRP000001", payloads=[rna, atac]
    )
    geo = {
        "accession": "SRP000001",
        "retrieval_status": "complete",
        "retrieval_manifest": str(manifest),
        "data_type": "mixed",
        "files": {"fastq": [str(rna), str(atac)]},
        "file_modalities": {
            str(rna): "bulk_RNA_raw",
            str(atac): "bulk_ATAC",
        },
        "inferred_design": {"organism": "Homo sapiens", "genome": "hg38"},
    }

    class _Memory:
        create_wing = lambda *args, **kwargs: None
        create_hall = lambda *args, **kwargs: None

    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent.memory = _Memory()
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_escalation = lambda *args, **kwargs: None
    result = agent.run("exp", {
        "data_dir": str(tmp_path), "user_question": "compare", "geo_metadata": geo,
    })

    modalities = result["exp_context"]["modalities"]
    assert modalities["bulk_RNA_raw"] == [str(rna)]
    assert modalities["bulk_ATAC"] == [str(atac)]


def test_e6_data_audit_excludes_manifested_transport_containers(tmp_path):
    from aria.agents.data_audit_agent import DataAuditAgent
    from aria.utils.atomic_retrieval import write_retrieval_manifest

    compressed = tmp_path / "processed.h5ad.gz"
    with gzip.open(compressed, "wb") as handle:
        handle.write(b"\x89HDF\r\n\x1a\ntransport")
    analyzable = tmp_path / "processed.h5ad"
    analyzable.write_bytes(b"\x89HDF\r\n\x1a\ndownstream")
    manifest = write_retrieval_manifest(
        tmp_path,
        accession="GSE000001",
        payloads=[compressed, analyzable],
    )
    geo = {
        "accession": "GSE000001", "retrieval_status": "complete",
        "retrieval_manifest": str(manifest), "data_type": "scRNA",
        "files": {"h5ad": [str(analyzable)]},
        "inferred_design": {"organism": "Homo sapiens", "genome": "hg38"},
    }

    class _Memory:
        create_wing = lambda *args, **kwargs: None
        create_hall = lambda *args, **kwargs: None

    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent.memory = _Memory()
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_escalation = lambda *args, **kwargs: None
    result = agent.run("exp", {
        "data_dir": str(tmp_path), "user_question": "compare", "geo_metadata": geo,
    })

    assert result["exp_context"]["modalities"]["scRNA"] == [str(analyzable)]


def test_e6_incomplete_explicit_geo_result_is_not_audited(tmp_path):
    from aria.agents.data_audit_agent import DataAuditAgent

    path = tmp_path / "partial.tsv"
    path.write_text("partial", encoding="utf-8")
    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent.publish_status = lambda *args, **kwargs: None
    result = agent.run("exp", {
        "data_dir": str(tmp_path), "user_question": "compare",
        "geo_metadata": {
            "retrieval_status": "incomplete",
            "files": {"counts": [str(path)]},
        },
    })

    assert result["status"] == "failed"
    assert "incomplete" in result["error"].lower()


def test_e6_complete_geo_result_is_revalidated_before_data_audit(tmp_path):
    from aria.agents.data_audit_agent import DataAuditAgent
    from aria.utils.atomic_retrieval import write_retrieval_manifest

    root = tmp_path / "published"
    root.mkdir()
    path = root / "counts.tsv"
    path.write_text("gene\tsample\nGENE_1\t1\n", encoding="utf-8")
    manifest = write_retrieval_manifest(
        root, accession="GSE000001", payloads=[path]
    )
    path.write_bytes(b"x" * path.stat().st_size)

    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent.publish_status = lambda *args, **kwargs: None
    result = agent.run("exp", {
        "data_dir": str(root), "user_question": "compare",
        "geo_metadata": {
            "accession": "GSE000001",
            "retrieval_status": "complete",
            "retrieval_manifest": str(manifest),
            "files": {"counts": [str(path)]},
        },
    })

    assert result["status"] == "failed"
    assert "manifest" in result["error"].lower()


def test_e6_geo_result_cannot_reference_a_payload_outside_its_manifest(tmp_path):
    from aria.agents.data_audit_agent import DataAuditAgent
    from aria.utils.atomic_retrieval import write_retrieval_manifest

    root = tmp_path / "published"
    root.mkdir()
    inside = root / "counts.tsv"
    inside.write_text("gene\tsample\nGENE_1\t1\n", encoding="utf-8")
    manifest = write_retrieval_manifest(
        root, accession="GSE000001", payloads=[inside]
    )
    outside = tmp_path / "outside.tsv"
    outside.write_text("gene\tsample\nGENE_2\t2\n", encoding="utf-8")

    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent.publish_status = lambda *args, **kwargs: None
    result = agent.run("exp", {
        "data_dir": str(root), "user_question": "compare",
        "geo_metadata": {
            "accession": "GSE000001",
            "retrieval_status": "complete",
            "retrieval_manifest": str(manifest),
            "files": {"counts": [str(outside)]},
        },
    })

    assert result["status"] == "failed"
    assert "outside" in result["error"].lower()


def test_e6_tui_handles_a_connector_result_without_payloads(monkeypatch, tmp_path):
    from aria import tui

    result = {
        "accession": "GSE000001", "title": "Empty public study",
        "organism": "Unknown", "genome": "unknown", "data_type": "unknown",
        "n_samples": 0, "local_dir": str(tmp_path), "files": {},
        "inferred_design": {"groups": {}},
    }

    class _Connector:
        def fetch(self, accession, status_callback=None):
            return result

    monkeypatch.setattr("aria.connectors.geo_connector.GEOConnector", _Connector)
    assert tui._resolve_geo_accession("GSE000001") == result
