"""scATAC P2 regulatory-layer guards.

Synthetic coordinates only. These tests verify that P2 runs from explicit
technical inputs and otherwise skips honestly; no biology is inferred from names.
"""

from pathlib import Path

import pytest


def test_regulatory_helpers_parse_coordinates_and_gtf(tmp_path):
    from aria.scripts.chromatin_regulatory import (
        _genes_to_peak_indices,
        _parse_gtf_tss,
        _peak_coord,
    )

    assert _peak_coord("chr1:10-20") == ("chr1", 10, 20)
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\ttest\tgene\t100\t150\t.\t+\t.\tgene_id "g1"; gene_name "GENE1";\n'
        'chr2\ttest\tgene\t500\t550\t.\t-\t.\tgene_id "g2"; gene_name "GENE2";\n',
        encoding="utf-8",
    )
    tss = _parse_gtf_tss(str(gtf))
    assert tss["GENE1"] == ("chr1", 100)
    assert tss["GENE2"] == ("chr2", 550)
    mapping = _genes_to_peak_indices(
        ["chr1:90-120", "chr2:10-20"], tss, window_bp=50)
    assert mapping == {"GENE1": [0]}


def test_vectorized_corr_matches_per_pair_corr():
    """P4.2 refactor guard: ``_vectorized_corr`` must reproduce the per-pair
    ``_corr`` it replaced in ``_peak_to_gene`` — same value (rounded to 6) on
    defined columns, and the same undefined cases (zero variance / n < 3), where
    NaN from the vectorized path corresponds to ``None`` from the scalar path."""
    import math

    import numpy as np

    from aria.scripts.chromatin_regulatory import _corr, _vectorized_corr

    x = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    cols = np.column_stack([
        2.0 * x + 0.5,                       # perfect positive
        -x + 7.0,                            # perfect negative
        np.array([1.0, 0.0, 3.0, 1.0, 4.0, 2.0]),  # partial
        np.array([2.0, 2.0, 2.0, 2.0, 2.0, 2.0]),  # zero variance -> undefined
    ])

    rvec = _vectorized_corr(x, cols)
    assert rvec.shape == (cols.shape[1],)
    for j in range(cols.shape[1]):
        scalar = _corr(x, cols[:, j])
        vec = rvec[j]
        if scalar is None:
            # undefined for the scalar path -> NaN in the vectorized path
            assert math.isnan(float(vec))
        else:
            assert round(float(vec), 6) == pytest.approx(scalar)


def test_vectorized_corr_undefined_when_too_few_cells():
    import math

    import numpy as np

    from aria.scripts.chromatin_regulatory import _corr, _vectorized_corr

    x = np.array([0.0, 1.0])
    cols = np.array([[0.0], [1.0]])
    assert _corr(x, cols[:, 0]) is None
    rvec = _vectorized_corr(x, cols)
    assert rvec.shape == (1,)
    assert math.isnan(float(rvec[0]))


def _write_atac(tmp_path):
    ad = pytest.importorskip("anndata")
    import numpy as np

    cells = [f"c{i}" for i in range(6)]
    peaks = ["chr1:90-120", "chr1:180-220", "chr2:10-40"]
    x = np.array([
        [0, 1, 0],
        [1, 1, 0],
        [2, 0, 0],
        [3, 0, 1],
        [4, 0, 1],
        [5, 0, 1],
    ], dtype="float32")
    adata = ad.AnnData(x)
    adata.obs_names = cells
    adata.var_names = peaks
    adata.obs["n_fragments"] = x.sum(axis=1)
    path = tmp_path / "atac.h5ad"
    adata.write_h5ad(path)
    return path


def _write_rna(tmp_path):
    ad = pytest.importorskip("anndata")
    import numpy as np

    cells = [f"c{i}" for i in range(6)]
    x = np.array([[0], [1], [2], [3], [4], [5]], dtype="float32")
    adata = ad.AnnData(x)
    adata.obs_names = cells
    adata.var_names = ["GENE1"]
    adata.obs["cell_type"] = ["A", "A", "B", "B", "B", "B"]
    path = tmp_path / "rna.h5ad"
    adata.write_h5ad(path)
    return path


def _write_gtf(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\ttest\tgene\t100\t150\t.\t+\t.\tgene_id "g1"; gene_name "GENE1";\n',
        encoding="utf-8",
    )
    return gtf


def test_regulatory_script_skips_missing_optional_inputs(tmp_path):
    pytest.importorskip("anndata")
    from aria.scripts.chromatin_regulatory import chromatin_regulatory

    atac = _write_atac(tmp_path)
    res = chromatin_regulatory({
        "data_path": str(atac),
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "success"
    assert res["motif_activity"]["ran"] is False
    assert res["gene_scores"]["ran"] is False
    assert res["footprinting"]["ran"] is False
    assert res["peak_to_gene"]["ran"] is False
    assert res["label_transfer"]["trusted_for_inference"] is False


def test_regulatory_script_runs_with_explicit_maps_gtf_and_paired_rna(tmp_path):
    pytest.importorskip("anndata")
    from aria.scripts.chromatin_regulatory import chromatin_regulatory

    atac = _write_atac(tmp_path)
    rna = _write_rna(tmp_path)
    gtf = _write_gtf(tmp_path)
    motif_map = tmp_path / "motif_peak_map.csv"
    motif_map.write_text(
        "motif_id,peak\nMOTIF_A,chr1:90-120\nMOTIF_A,chr1:180-220\n",
        encoding="utf-8",
    )

    res = chromatin_regulatory({
        "data_path": str(atac),
        "output_dir": str(tmp_path / "out"),
        "motif_peak_map": str(motif_map),
        "gtf_path": str(gtf),
        "rna_data_path": str(rna),
        "rna_label_col": "cell_type",
        "min_shared_cells": 3,
        "min_peak_gene_corr": 0.8,
        "peak_gene_distance_bp": 1000,
    })

    assert res["status"] == "success"
    assert res["motif_activity"]["ran"] is True
    assert res["motif_activity"]["n_motifs"] == 1
    assert Path(res["motif_activity"]["output_csv"]).exists()
    assert res["gene_scores"]["ran"] is True
    assert res["gene_scores"]["n_genes_scored"] == 1
    assert res["peak_to_gene"]["ran"] is True
    assert res["peak_to_gene"]["n_links"] >= 1
    assert res["label_transfer"]["ran"] is True
    assert res["label_transfer"]["trusted_for_inference"] is False
    assert res["footprinting"]["ran"] is False

    # P4.2/ADR-050 scoped promotion: ONLY peak-to-gene is beta-grade; the script
    # umbrella stays alpha and the other regulatory layers are NOT over-promoted.
    assert res["peak_to_gene"]["validation_level"] == "beta"
    assert res["validation_level"] == "alpha"
    for other in ("motif_activity", "gene_scores", "footprinting", "label_transfer"):
        assert res[other].get("validation_level") != "beta"
