"""D8 convergence regression tests for count readers, multiome loaders, HDF5 obs decoding, and 10x barcodes."""

from __future__ import annotations

from pathlib import Path


def test_d8_010_bulk_count_readers_converge_on_featurecounts_and_multifile(tmp_path: Path):
    """D8-010: all bulk count readers must strip featureCounts metadata, default to tab, and keep multi-file samples."""
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
