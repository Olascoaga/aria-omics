"""D8 convergence regression tests for count readers, multiome loaders, HDF5 obs decoding, and 10x barcodes."""

from __future__ import annotations

from pathlib import Path


def test_d8_010_bulk_count_readers_converge_on_featurecounts_and_multifile(tmp_path: Path):
    """D8-010: all bulk count readers must strip featureCounts metadata, default to tab, and keep multi-file samples."""
    import pytest

    pytest.importorskip("litellm")
    from aria.agents.audit_agent import _load_count_matrix
    from aria.scripts.chromatin_figure_bulk import _load_counts as chrom_fig_load
    from aria.scripts.rna_bulk.counts import _load_counts

    fc = tmp_path / "counts.featureCounts"
    rows = ["Geneid\tChr\tStart\tEnd\tStrand\tLength\tctrl_1\tctrl_2\ttreat_1\ttreat_2"]
    for i in range(30):
        rows.append(
            f"ENSG{i:011d}\tchr1\t{1000 + i * 2000}\t{1500 + i * 2000}\t+\t500\t"
            f"{100 + i}\t{110 + i}\t{200 + i}\t{210 + i}"
        )
    fc.write_text("\n".join(rows) + "\n", encoding="utf-8")

    c_bulk, _, _ = _load_counts([str(fc)])
    c_audit = _load_count_matrix([str(fc)])
    c_fig, _ = chrom_fig_load(str(fc))

    expected_cols = ["ctrl_1", "ctrl_2", "treat_1", "treat_2"]
    assert c_bulk is not None and list(c_bulk.columns) == expected_cols
    assert c_audit is not None and list(c_audit.columns) == expected_cols
    assert c_fig is not None and list(c_fig.columns) == expected_cols

    # Multi-file headered TSVs sharing the column name 'count' must keep every sample uniquely.
    s1 = tmp_path / "sample_A.tsv"
    s2 = tmp_path / "sample_B.tsv"
    s1.write_text("gene_id\tcount\nG1\t10\nG2\t20\n", encoding="utf-8")
    s2.write_text("gene_id\tcount\nG1\t30\nG2\t40\n", encoding="utf-8")
    m_bulk, _, _ = _load_counts([str(s1), str(s2)])
    m_audit = _load_count_matrix([str(s1), str(s2)])
    assert m_bulk is not None and list(m_bulk.columns) == ["sample_A", "sample_B"]
    assert m_audit is not None and list(m_audit.columns) == ["sample_A", "sample_B"]


def test_d8_013_multiome_10x_atac_loaders_converge(tmp_path: Path):
    """D8-013: WNN, MOFA+, Peak2Gene, and eRegulons must all load 10x Multiome .h5 Peaks and RNA."""
    import pytest

    h5py = pytest.importorskip("h5py")
    pytest.importorskip("anndata")
    import numpy as np
    from scipy import sparse

    from aria.scripts.integration_eregulons import _load_atac_adata as ereg_atac
    from aria.scripts.integration_mofa import _load_modality
    from aria.scripts.integration_peak2gene import _load_atac_adata as p2g_atac
    from aria.scripts.integration_wnn import _load_atac as wnn_atac, _load_rna as wnn_rna
    from aria.utils.mudata_io import _read_10x_matrix

    h5_path = tmp_path / "filtered_feature_bc_matrix.h5"
    mat = sparse.csc_matrix(np.array([[1, 2], [3, 0], [5, 1]], dtype=np.int32))
    with h5py.File(h5_path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("data", data=mat.data)
        g.create_dataset("indices", data=mat.indices)
        g.create_dataset("indptr", data=mat.indptr)
        g.create_dataset("shape", data=np.array(mat.shape, dtype=np.int64))
        g.create_dataset("barcodes", data=np.array([b"AAACCCAAGAAACACT-1", b"AAACCCAAGAAACCAT-1"]))
        fg = g.create_group("features")
        fg.create_dataset("id", data=np.array([b"ENSG0001", b"chr1:100-500", b"chr1:1000-1500"]))
        fg.create_dataset("name", data=np.array([b"TP53", b"chr1:100-500", b"chr1:1000-1500"]))
        fg.create_dataset("feature_type", data=np.array([b"Gene Expression", b"Peaks", b"Peaks"]))
        fg.create_dataset("genome", data=np.array([b"GRCh38", b"GRCh38", b"GRCh38"]))

    _read_10x_matrix(str(h5_path), require_paired=True)
    p2g_ad, _ = p2g_atac([str(h5_path)])
    ereg_ad, _ = ereg_atac([str(h5_path)])
    mofa_ad = _load_modality("atac", str(h5_path))
    wnn_ad = wnn_atac([str(h5_path)], "hg38", "Homo sapiens")
    wnn_rna_ad = wnn_rna(str(h5_path))

    assert p2g_ad is not None and p2g_ad.n_vars == 2
    assert ereg_ad is not None and ereg_ad.n_vars == 2
    assert mofa_ad is not None and mofa_ad.n_vars == 2
    assert wnn_ad is not None and wnn_ad.n_vars == 2
    assert wnn_rna_ad is not None and wnn_rna_ad.n_vars == 1

    # Raw fragments must still raise NotImplementedError in WNN (P0-8).
    with pytest.raises(NotImplementedError):
        wnn_atac([str(tmp_path / "atac_fragments.tsv.gz")], "hg38", "Homo sapiens")
