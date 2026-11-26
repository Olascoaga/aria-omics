"""Permanent regression tests for D1 silent-failure fixes in bulk RNA compute kernels."""
from __future__ import annotations

import gzip

import pytest


def test_d1_003_gtf_truncated_must_not_return_partial_biotype(tmp_path):
    """D1-003: Truncated GTF must fail closed instead of returning a partial biotype_map."""
    from aria.scripts.rna_bulk.gtf_io import _parse_gtf_biotype_and_length

    gtf_gz = tmp_path / "trunc.gtf.gz"
    lines = [
        f'chr1\tHAVANA\tgene\t{1000 + i * 100}\t{1050 + i * 100}\t.\t+\t.\t'
        f'gene_id "ENSG{i:011d}"; gene_biotype "protein_coding";\n'
        for i in range(600)
    ]
    raw_gz = gzip.compress("".join(lines).encode("utf-8"))
    gtf_gz.write_bytes(raw_gz[: len(raw_gz) // 2])
    with pytest.raises((EOFError, OSError, ValueError)):
        _parse_gtf_biotype_and_length(str(gtf_gz))
