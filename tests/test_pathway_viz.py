from __future__ import annotations


def test_ranked_signature_keeps_symbol_ids_without_symbol_map():
    pd = __import__("pandas")
    from aria.scripts.rna_pathway_viz import _ranked_signature_frame

    de = pd.DataFrame(
        {"log2FoldChange": [1.2, -2.5, 0.4]},
        index=["GENE_A", "GENE_B", "GENE_A"],
    )

    signature = _ranked_signature_frame(de, symbol_map={})

    assert signature.to_dict("records") == [
        {"gene": "GENE_B", "score": -2.5},
        {"gene": "GENE_A", "score": 1.2},
    ]


def test_ranked_signature_drops_unmapped_ensembl_ids():
    pd = __import__("pandas")
    from aria.scripts.rna_pathway_viz import _ranked_signature_frame

    de = pd.DataFrame(
        {"log2FoldChange": [1.2, -2.5]},
        index=["ENSG000001.1", "GENE_B"],
    )

    signature = _ranked_signature_frame(de, symbol_map={})

    assert signature.to_dict("records") == [{"gene": "GENE_B", "score": -2.5}]
