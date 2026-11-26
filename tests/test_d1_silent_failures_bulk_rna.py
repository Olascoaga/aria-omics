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


def test_d1_006_bulk_counts_multifile_must_fail_closed_on_corrupt_sample(tmp_path):
    """D1-006: Multi-file bulk counts must fail closed when a declared sample file is corrupt/empty."""
    pytest.importorskip("pydeseq2")
    import numpy as np
    from aria.scripts.rna_bulk_de import bulk_rna_de

    rng = np.random.default_rng(7)
    genes = [f"GENE{i}" for i in range(60)]
    specs = [
        ("ctrl_1", "ctrl"), ("ctrl_2", "ctrl"), ("ctrl_3", "ctrl"), ("ctrl_4", "ctrl"),
        ("trt_1", "treat"), ("trt_2", "treat"), ("trt_3", "treat"), ("trt_4", "treat"),
    ]
    files = []
    for sid, cond in specs:
        p = tmp_path / f"{sid}.tsv"
        if sid == "trt_4":
            p.write_text("")
        else:
            vals = rng.poisson(50 if cond == "ctrl" else 200, size=len(genes))
            p.write_text("\n".join(f"{g}\t{v}" for g, v in zip(genes, vals)) + "\n")
        files.append(str(p))
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tcondition\n" + "\n".join(f"{s}\t{c}" for s, c in specs) + "\n")
    res = bulk_rna_de({
        "files": files,
        "metadata_file": str(meta),
        "design_factor": "condition",
        "contrasts": [{"numerator": "treat", "denominator": "ctrl", "name": "t_vs_c"}],
        "run_pathways": False,
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "error"
