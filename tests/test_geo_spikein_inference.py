"""Deferred WIP closeout: GEO multi-organism (spike-in) handling.

Pure helper for spike-in GEO series, where the SOFT metadata lists several
organisms:

- geo_connector._organism_from_gene_symbols: infer the experimental organism from
  gene-ID style (a technical species-detection, like the ADR-011 `human_markers`
  exception — it picks a reference genome, it makes no biological claim).

The former ``BulkRNAAgent._infer_col_groups`` column-name group recovery was
removed by preprint-readiness audit B4 (no post-confirmation column-name
inference; an unmappable confirmed design fails closed). See
``tests/test_preprint_audit_b4_no_name_fallback.py``.
"""

import pytest


def _write_counts(path, genes):
    path.write_text("gene\ts1\ts2\ts3\ts4\n"
                    + "".join(f"{g}\t1\t2\t3\t4\n" for g in genes))
    return str(path)


def test_organism_from_human_hgnc_symbols(tmp_path):
    from aria.connectors.geo_connector import _organism_from_gene_symbols
    p = _write_counts(tmp_path / "h.tsv",
                      ["CXCL1", "IL6", "TP53", "GAPDH", "ACTB", "STAT1"])
    assert _organism_from_gene_symbols(p) == "homo sapiens"


def test_organism_from_mouse_mgi_symbols(tmp_path):
    from aria.connectors.geo_connector import _organism_from_gene_symbols
    p = _write_counts(tmp_path / "m.tsv",
                      ["Cxcl1", "Il6", "Trp53", "Gapdh", "Actb", "Stat1"])
    assert _organism_from_gene_symbols(p) == "mus musculus"


def test_organism_from_ensembl_prefixes(tmp_path):
    from aria.connectors.geo_connector import _organism_from_gene_symbols
    human = _write_counts(tmp_path / "eh.tsv",
                          [f"ENSG00000{i:06d}" for i in range(6)])
    mouse = _write_counts(tmp_path / "em.tsv",
                          [f"ENSMUSG00000{i:06d}" for i in range(6)])
    assert _organism_from_gene_symbols(human) == "homo sapiens"
    assert _organism_from_gene_symbols(mouse) == "mus musculus"


def test_organism_inference_is_uncertain_for_mixed_ids(tmp_path):
    from aria.connectors.geo_connector import _organism_from_gene_symbols
    p = _write_counts(tmp_path / "x.tsv", ["1", "2", "foo", "bar"])
    assert _organism_from_gene_symbols(p) == ""


def test_infer_col_groups_helper_is_removed():
    """B4: the column-name group-inference back door no longer exists."""
    pytest.importorskip("litellm")  # importing the agent pulls aria.llm.provider
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    assert not hasattr(BulkRNAAgent, "_infer_col_groups")
