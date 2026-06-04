"""Real-run bug (2026-06-04): bulk RNA-seq FASTQs (B1_1.fq.gz, ...) were
classified as bulk_ATAC. Root cause: AssayDetector._detect_alignment treated ANY
gzipped file as alignment content and scanned its decompressed bytes for assay
keywords — but a gzipped FASTQ decompresses to read sequences, and ACGT reads
contain the 4-mer "atac" by chance, so it confidently returned bulk_ATAC.
"""

import gzip

from aria.utils.assay_detector import AssayDetector


def _write_fastq_gz(path, n=50):
    # Sequences deliberately contain the "atac" 4-mer (as real reads do).
    seq = "GATACACGGATACAGTTACGT"
    rec = "".join(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i in range(n))
    with gzip.open(path, "wt") as fh:
        fh.write(rec)


def test_gzipped_fastq_is_not_content_detected_as_alignment(tmp_path):
    p = tmp_path / "B1_1.fq.gz"
    _write_fastq_gz(p)
    det = AssayDetector().detect_file(p)
    # A FASTQ is not a recognizable content type for this detector; it must defer
    # (return None) so the filename/path signal decides — NOT claim bulk_ATAC.
    assert det is None, f"FASTQ.gz misdetected as {det and det.modality}"


def test_data_audit_classifies_paired_rna_fastq_as_bulk_rna_raw(tmp_path):
    from aria.agents.data_audit_agent import DataAuditAgent

    files = []
    for cond in ("B1", "R1", "WT1"):
        for mate in ("1", "2"):
            p = tmp_path / f"{cond}_{mate}.fq.gz"
            _write_fastq_gz(p)
            files.append(p)

    agent = DataAuditAgent.__new__(DataAuditAgent)
    classified = agent._classify_files(files)
    assert set(classified) == {"bulk_RNA_raw"}, (
        f"paired RNA FASTQs must classify as bulk_RNA_raw, got {dict((k, len(v)) for k, v in classified.items())}")
    assert len(classified["bulk_RNA_raw"]) == 6


def test_real_bam_like_gzip_still_detected_as_alignment(tmp_path):
    # Guard the fix doesn't break real BAM detection: a gzip whose decompressed
    # content begins with the BAM magic must still be treated as alignment.
    p = tmp_path / "aln.bam"
    payload = b"BAM\x01" + (4).to_bytes(4, "little") + b"@HD\tVN:1.6\n" + b"\x00" * 32
    with gzip.open(p, "wb") as fh:
        fh.write(payload)
    det = AssayDetector().detect_file(p)
    assert det is not None
    assert det.evidence.get("format") == "bam_or_sam"
