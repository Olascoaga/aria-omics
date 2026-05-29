"""Bulk DE raw-count guard (audit 2026-05-29, B10).

`_load_counts` must hard-refuse a non-raw matrix (TPM/CPM/log-normalized) rather
than silently round it into pseudo-counts for DESeq2, and must coerce only when
`allow_nonraw` is explicit, recording provenance either way. Needs pandas only.
"""

import numpy as np
import pandas as pd
import pytest

from aria.scripts.rna_bulk_de import _load_counts


def _write_matrix(tmp_path, mat, name="counts.tsv"):
    genes = [f"G{i}" for i in range(mat.shape[0])]
    samples = [f"S{j}" for j in range(mat.shape[1])]
    df = pd.DataFrame(mat, index=genes, columns=samples)
    path = tmp_path / name
    df.to_csv(path, sep="\t")
    return str(path)


def test_raw_counts_load_and_are_tagged_raw(tmp_path):
    rng = np.random.default_rng(0)
    mat = rng.poisson(40, size=(300, 6))
    path = _write_matrix(tmp_path, mat)
    counts, warnings, meta = _load_counts([path])
    assert counts is not None
    assert meta["count_source"] == "raw_counts"
    assert str(counts.dtypes.iloc[0]).startswith("int")


def test_lognorm_matrix_is_hard_refused_by_default(tmp_path):
    rng = np.random.default_rng(1)
    mat = rng.uniform(0, 8, size=(300, 6))  # non-integer, bounded
    path = _write_matrix(tmp_path, mat)
    counts, warnings, meta = _load_counts([path])
    assert counts is None
    assert meta["refused"] is True
    assert meta["error_type"] == "NonRawCounts"


def test_nonraw_matrix_coerced_only_when_allowed(tmp_path):
    rng = np.random.default_rng(2)
    mat = rng.uniform(0, 5000, size=(300, 6))  # TPM-like continuous
    path = _write_matrix(tmp_path, mat)
    counts, warnings, meta = _load_counts([path], allow_nonraw=True)
    assert counts is not None
    assert meta["count_source"] == "coerced_nonraw"
    assert any("LOW CONFIDENCE" in w for w in warnings)
