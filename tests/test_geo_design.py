from __future__ import annotations


def test_geo_design_prefers_experimental_keys_and_uses_geo_ids():
    from aria.connectors.geo_connector import _infer_design

    metadata = {
        "samples": [
            {
                "id": "GSM0001",
                "title": "Control replicate 1",
                "characteristics": {
                    "sex": "female",
                    "genotype": "WT",
                    "sample id": "donor_1",
                },
            },
            {
                "id": "GSM0002",
                "title": "Control replicate 2",
                "characteristics": {
                    "sex": "male",
                    "genotype": "WT",
                    "sample id": "donor_2",
                },
            },
            {
                "id": "GSM0003",
                "title": "Knockout replicate 1",
                "characteristics": {
                    "sex": "female",
                    "genotype": "KO",
                    "sample id": "donor_3",
                },
            },
            {
                "id": "GSM0004",
                "title": "Knockout replicate 2",
                "characteristics": {
                    "sex": "male",
                    "genotype": "KO",
                    "sample id": "donor_4",
                },
            },
        ]
    }

    design = _infer_design(metadata)

    assert design["condition_col"] == "genotype"
    assert design["main_factor"] == "genotype"
    assert design["groups"] == {
        "KO": ["GSM0003", "GSM0004"],
        "WT": ["GSM0001", "GSM0002"],
    }
    assert design["sample_aliases"]["GSM0001"] == [
        "GSM0001",
        "Control replicate 1",
    ]


def test_bulk_design_maps_geo_ids_and_title_aliases_to_count_columns(tmp_path):
    from aria.agents.bulk_rna_agent import BulkRNAAgent

    counts = tmp_path / "counts.tsv"
    counts.write_text(
        "gene_id\tGSM0001\tGSM0002\tKnockout replicate 1\tKnockout replicate 2\n"
        "GENE_1\t10\t11\t80\t75\n",
        encoding="utf-8",
    )
    design = {
        "groups": {
            "WT": ["GSM0001", "GSM0002"],
            "KO": ["GSM0003", "GSM0004"],
        },
        "main_factor": "genotype",
        "plan_contrasts": [{"numerator": "KO", "denominator": "WT"}],
        "sample_aliases": {
            "GSM0001": ["GSM0001", "Control replicate 1"],
            "GSM0002": ["GSM0002", "Control replicate 2"],
            "GSM0003": ["GSM0003", "Knockout replicate 1"],
            "GSM0004": ["GSM0004", "Knockout replicate 2"],
        },
    }

    agent = BulkRNAAgent.__new__(BulkRNAAgent)
    sample_names, group_labels, factor, contrasts = agent._apply_design(
        design,
        [str(counts)],
        "exp_geo",
    )

    assert sample_names == [
        "GSM0001",
        "GSM0002",
        "Knockout replicate 1",
        "Knockout replicate 2",
    ]
    assert group_labels == {
        "GSM0001": "WT",
        "GSM0002": "WT",
        "Knockout replicate 1": "KO",
        "Knockout replicate 2": "KO",
    }
    assert factor == "genotype"
    assert contrasts == [{
        "numerator": "KO",
        "denominator": "WT",
        "name": "KO vs WT",
    }]
