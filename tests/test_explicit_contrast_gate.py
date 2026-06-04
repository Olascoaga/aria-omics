"""P0-5 regression: DE must not choose reference levels implicitly.

Bulk and pseudobulk DE are publication-facing tests. If ARIA runs them without
an explicit numerator/test and denominator/reference, the first sorted group can
become the reference by accident. These tests lock the production behavior:
missing contrasts produce suggestions and a checkpoint, not DE execution.
"""

from __future__ import annotations


def test_bulk_script_requires_explicit_contrasts(tmp_path):
    import numpy as np
    import pandas as pd
    from aria.scripts.rna_bulk_de import bulk_rna_de

    rng = np.random.default_rng(5)
    counts = pd.DataFrame(
        rng.poisson(80, size=(80, 6)),
        index=[f"GENE_{i:03d}" for i in range(80)],
        columns=["A_1", "A_2", "A_3", "B_1", "B_2", "B_3"],
    )
    counts_path = tmp_path / "counts.tsv"
    counts.to_csv(counts_path, sep="\t")

    result = bulk_rna_de({
        "files": [str(counts_path)],
        "design_factor": "condition",
        "run_pathways": False,
        "output_dir": str(tmp_path / "out"),
    })

    assert result["status"] == "error"
    assert result["error_type"] == "ExplicitContrastRequired"
    assert result["available_groups"] == ["A", "B"]
    assert result["suggested_contrasts"] == [{
        "numerator": "B",
        "denominator": "A",
        "name": "B vs A",
    }]


def test_bulk_agent_uses_only_confirmed_plan_contrasts(tmp_path):
    from aria.agents.bulk_rna_agent import BulkRNAAgent

    counts = tmp_path / "counts.tsv"
    counts.write_text(
        "gene_id\tS1\tS2\tS3\tS4\n"
        "GENE_1\t10\t11\t80\t75\n",
        encoding="utf-8",
    )
    base_design = {
        "groups": {"A": ["S1", "S2"], "B": ["S3", "S4"]},
        "main_factor": "condition",
    }
    agent = BulkRNAAgent.__new__(BulkRNAAgent)

    _, _, _, contrasts = agent._apply_design(base_design, [str(counts)], "exp")
    assert contrasts == []

    confirmed = {
        **base_design,
        "plan_contrasts": [{"numerator": "B", "denominator": "A"}],
    }
    _, _, _, contrasts = agent._apply_design(confirmed, [str(counts)], "exp")
    assert contrasts == [{
        "numerator": "B",
        "denominator": "A",
        "name": "B vs A",
    }]


def test_bulk_agent_writes_confirmed_design_metadata_for_raw_fastq_counts(tmp_path):
    from aria.agents.bulk_rna_agent import BulkRNAAgent

    counts = tmp_path / "counts.tsv"
    counts.write_text(
        "gene_id\tB1\tB2\tB3\tR1\tR2\tR3\tWT1\tWT2\tWT3\n"
        "GENE_1\t10\t11\t12\t80\t75\t77\t20\t22\t19\n",
        encoding="utf-8",
    )
    design = {
        "groups": {
            "BMAL1_KO": ["B1", "B2", "B3"],
            "REV-ERBa_KO": ["R1", "R2", "R3"],
            "WT": ["WT1", "WT2", "WT3"],
        },
        "main_factor": "genotype",
        "plan_contrasts": [
            {"numerator": "BMAL1_KO", "denominator": "WT"},
            # LLM plan text may drop punctuation; this must resolve to the
            # confirmed design level rather than being discarded.
            {"numerator": "REVERBa_KO", "denominator": "WT"},
            {"numerator": "BMAL1_KO", "denominator": "REVERBa_KO"},
        ],
    }

    agent = BulkRNAAgent.__new__(BulkRNAAgent)
    sample_names, group_labels, factor, contrasts = agent._apply_design(
        design, [str(counts)], "exp-h9"
    )
    metadata_path = agent._write_design_metadata(
        group_labels, factor, tmp_path / "de_out"
    )

    assert sample_names == ["B1", "B2", "B3", "R1", "R2", "R3", "WT1", "WT2", "WT3"]
    assert group_labels == {
        "B1": "BMAL1_KO",
        "B2": "BMAL1_KO",
        "B3": "BMAL1_KO",
        "R1": "REV-ERBa_KO",
        "R2": "REV-ERBa_KO",
        "R3": "REV-ERBa_KO",
        "WT1": "WT",
        "WT2": "WT",
        "WT3": "WT",
    }
    assert contrasts == [
        {"numerator": "BMAL1_KO", "denominator": "WT", "name": "BMAL1_KO vs WT"},
        {"numerator": "REV-ERBa_KO", "denominator": "WT", "name": "REV-ERBa_KO vs WT"},
        {
            "numerator": "BMAL1_KO",
            "denominator": "REV-ERBa_KO",
            "name": "BMAL1_KO vs REV-ERBa_KO",
        },
    ]
    assert metadata_path.read_text(encoding="utf-8").splitlines() == [
        "sample\tgenotype",
        "B1\tBMAL1_KO",
        "B2\tBMAL1_KO",
        "B3\tBMAL1_KO",
        "R1\tREV-ERBa_KO",
        "R2\tREV-ERBa_KO",
        "R3\tREV-ERBa_KO",
        "WT1\tWT",
        "WT2\tWT",
        "WT3\tWT",
    ]


def test_scrna_pseudobulk_missing_comparison_requests_confirmation(tmp_path):
    from aria.agents.scrna_agent import scRNAAgent

    h5ad_path = tmp_path / "input.h5ad"
    h5ad_path.write_text("placeholder", encoding="utf-8")
    findings = []
    escalations = []

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.publish_finding = lambda *args, **kwargs: findings.append(args)
    agent.publish_escalation = lambda **kwargs: escalations.append(kwargs)

    result = agent._run_pseudobulk(
        "exp",
        str(h5ad_path),
        {
            "design": {
                "groups": {"A": ["r1", "r2", "r3"], "B": ["r4", "r5", "r6"]},
                "main_factor": "condition",
                "pseudobulk": {
                    "from_obs": True,
                    "condition_col": "condition",
                    "replicate_col": "donor",
                    "groupby_col": "cell_type",
                },
            },
        },
        {},
        {},
    )

    assert result["status"] == "skipped"
    assert result["reason"] == "explicit_comparison_required"
    assert result["suggested_comparisons"] == [["B", "A"]]
    assert findings
    assert escalations
    assert escalations[0]["checkpoint"] == "scrna.pseudobulk.contrast"
