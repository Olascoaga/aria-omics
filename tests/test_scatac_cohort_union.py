"""Cohort unions preserve raw counts, interval overlap, and donor provenance."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from aria.utils.scatac_cohort_union import merge_scatac_matrices, resolve_sample_ids
from aria.scripts.chromatin_merge_matrices import chromatin_merge_matrices


def _write_matrix(path, peaks, counts, *, obs=None, raw_counts=None, coordinates=None):
    matrix = ad.AnnData(sparse.csr_matrix(counts),
                        obs=pd.DataFrame(obs or {}, index=["CELL"]),
                        var=pd.DataFrame(index=peaks))
    if raw_counts is not None:
        matrix.layers["counts"] = sparse.csr_matrix(raw_counts)
    if coordinates is not None:
        matrix.var[["chrom", "start", "end"]] = pd.DataFrame(coordinates, index=matrix.var_names)
    if path.suffix == ".h5mu":
        mudata = pytest.importorskip("mudata")
        mudata.MuData({"atac": matrix}).write_h5mu(path)
    else:
        matrix.write_h5ad(path)
    return str(path)


@pytest.mark.parametrize("suffixes", [(".h5ad", ".h5ad"), (".h5mu", ".h5ad"), (".h5mu", ".h5mu")])
def test_shifted_interval_union_and_count_projection(tmp_path, suffixes):
    paths = [
        _write_matrix(tmp_path / ("first" + suffixes[0]),
                      ["chr1:100-200", "chr1-180-240", "chr1_300_400", "chr2:0-10"],
                      [[1, 2, 3, 4]], obs={"donor_id": ["existing_donor"], "batch": ["run1"]}),
        _write_matrix(tmp_path / ("second" + suffixes[1]),
                      ["chr1:220-280", "chr1:400-450", "chr2:5-15"], [[5, 6, 7]]),
    ]
    output = tmp_path / "cohort.h5ad"
    result = merge_scatac_matrices(paths, str(output), sample_ids=["first", "second"],
                                   sample_metadata={"first": {"condition": "group_a"},
                                                    "second": {"donor_id": "donor_b", "condition": "group_b"}})
    assert result["status"] == "success", result
    assert result["n_cells"] == 2 and result["n_consensus_peaks"] == 4
    assert result["peak_provenance"] == "cohort_interval_union_consensus"
    merged = ad.read_h5ad(output)
    assert list(merged.var_names) == ["chr1:100-280", "chr1:300-400", "chr1:400-450", "chr2:0-15"]
    assert sparse.issparse(merged.X)
    np.testing.assert_array_equal(merged.X.toarray(), [[3, 3, 0, 4], [5, 0, 6, 7]])
    np.testing.assert_array_equal(merged.layers["counts"].toarray(), merged.X.toarray())
    assert list(merged.obs_names) == ["first#CELL", "second#CELL"]
    assert list(merged.obs["donor_id"]) == ["existing_donor", "donor_b"]
    assert list(merged.obs["condition"]) == ["group_a", "group_b"]
    assert list(merged.obs["sample_id"]) == ["first", "second"]
    assert merged.obs.loc["first#CELL", "batch"] == "run1"
    assert list(merged.var["start"]) == [100, 300, 400, 0]


def test_counts_layer_and_coordinate_columns(tmp_path):
    first = _write_matrix(tmp_path / "a.h5ad", ["peak_a", "peak_b"], [[0.2, 0.3]],
                          raw_counts=[[10, 20]], coordinates=[("chrA", 1, 10), ("chrA", 5, 20)])
    second = _write_matrix(tmp_path / "b.h5ad", ["chrA_15_25"], [[4]])
    result = chromatin_merge_matrices({"files": [first, second], "output_dir": str(tmp_path / "out")})
    assert result["status"] == "success", result
    merged = ad.read_h5ad(result["output_path"])
    np.testing.assert_array_equal(merged.X.toarray(), [[30], [4]])
    assert list(merged.var_names) == ["chrA:1-25"]
    assert list(merged.obs["donor_id"]) == ["a", "b"]


def test_non_genomic_features_rejected(tmp_path):
    paths = [_write_matrix(tmp_path / "a.h5ad", ["feature_a"], [[1]]),
             _write_matrix(tmp_path / "b.h5ad", ["chrA:0-10"], [[1]])]
    result = merge_scatac_matrices(paths, str(tmp_path / "out.h5ad"))
    assert result["status"] == "error"
    assert result["error_type"] == "NonGenomicPeakFeatures"
    assert not (tmp_path / "out.h5ad").exists()


def test_normalized_only_counts_rejected(tmp_path):
    paths = [_write_matrix(tmp_path / "a.h5ad", ["chrA:0-10"], [[0.25]]),
             _write_matrix(tmp_path / "b.h5ad", ["chrA:0-10"], [[1]])]
    result = merge_scatac_matrices(paths, str(tmp_path / "out.h5ad"))
    assert result["error_type"] == "NonIntegerCounts"


def test_noninteger_layer_falls_back_to_raw_x(tmp_path):
    paths = [_write_matrix(tmp_path / "a.h5ad", ["chrA:0-10"], [[3]], raw_counts=[[0.2]]),
             _write_matrix(tmp_path / "b.h5ad", ["chrA:5-15"], [[4]])]
    result = merge_scatac_matrices(paths, str(tmp_path / "out.h5ad"))
    assert result["status"] == "success", result
    assert ad.read_h5ad(result["output_path"]).X.sum() == 7


def test_prefix_idempotent_and_duplicate_stems(tmp_path):
    assert resolve_sample_ids(["/a/sample.fragments.tsv.gz", "/b/sample.tsv.gz", "/c/sample_2.bed.gz"]) == ["sample", "sample_2", "sample_2_2"]
    first = _write_matrix(tmp_path / "a.h5ad", ["chrA:0-10"], [[1]])
    matrix = ad.read_h5ad(first)
    matrix.obs_names = ["a#CELL"]
    matrix.write_h5ad(first)
    second = _write_matrix(tmp_path / "b.h5ad", ["chrA:5-15"], [[2]])
    result = merge_scatac_matrices([first, second], str(tmp_path / "out.h5ad"))
    assert result["status"] == "success", result
    assert list(ad.read_h5ad(result["output_path"]).obs_names) == ["a#CELL", "b#CELL"]


def test_input_output_validation(tmp_path):
    result = merge_scatac_matrices([], str(tmp_path / "out.h5ad"))
    assert result["error_type"] == "InvalidInputs"
    result = merge_scatac_matrices(["a.h5ad", "b.h5ad"], str(tmp_path / "out.h5mu"))
    assert result["error_type"] == "InvalidOutputPath"
    result = merge_scatac_matrices(["a.h5ad", "b.h5ad"], str(tmp_path / "out.h5ad"))
    assert result["error_type"] == "MissingInput"
