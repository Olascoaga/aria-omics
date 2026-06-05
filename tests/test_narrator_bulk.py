def _bulk_findings():
    return {
        "sample_qc": {
            "n_samples": 6,
            "size_ratio": 1.4,
            "outliers": [],
        },
        "design_used": "~condition",
        "padj_threshold": 0.05,
        "lfc_threshold": 0.5,
        "contrasts": [{
            "name": "treat_vs_ctrl",
            "status": "success",
            "n_significant": 223,
            "n_upregulated": 120,
            "n_downregulated": 103,
            "power_estimate_at_lfc_min": 0.71,
            "top_genes": [
                {"gene": "GENE_UP_1", "log2fc": 2.2},
                {"gene": "GENE_UP_2", "log2fc": 1.8},
            ],
            "pathways": {
                "GO_BP": [{
                    "term": "pathway_alpha_response",
                    "adjusted_p": 1e-4,
                }]
            },
        }, {
            "name": "low_power_vs_ctrl",
            "status": "success",
            "n_significant": 8,
            "n_upregulated": 5,
            "n_downregulated": 3,
            "power_estimate_at_lfc_min": 0.31,
            "low_power_warning": True,
            "low_power_reason": "n<=2 replicates on one side.",
        }],
    }


def test_bulk_narrator_generates_qc_contrast_pathway_and_power_blocks():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
    from aria.agents.narrative.validators import validate_blocks

    agent_result = {"status": "done", "findings": _bulk_findings()}
    blocks = validate_blocks(BulkRnaNarrator().collect("bulk_rna_agent", agent_result))
    ids = {block.id for block in blocks}

    assert "bulk.qc" in ids
    assert "bulk.contrast.treat_vs_ctrl" in ids
    assert "bulk.pathway.treat_vs_ctrl" in ids
    assert "bulk.contrast.low_power_vs_ctrl" in ids
    assert "bulk.power" in ids

    contrast = next(b for b in blocks if b.id == "bulk.contrast.treat_vs_ctrl")
    assert contrast.claim == "Bulk contrast treat_vs_ctrl had 223 DE genes."
    assert any("GENE_UP_1" in ev.label for ev in contrast.evidence)

    low_power = next(b for b in blocks
                     if b.id == "bulk.contrast.low_power_vs_ctrl")
    assert any("n<=2" in caveat.text for caveat in low_power.caveats)

    power = next(b for b in blocks if b.id == "bulk.power")
    assert power.metrics["min_power"] == 0.31
    assert power.metrics["max_power"] == 0.71
    assert any(ev.value == "31%" for ev in power.evidence)
    assert any(ev.value == "71%" for ev in power.evidence)


def test_bulk_power_block_renders_under_strict_evidence_gate(tmp_path):
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
    from aria.agents.narrative.render_blocks import render_blocks
    from aria.agents.narrative.validators import validate_blocks

    agent_result = {"status": "done", "findings": _bulk_findings()}
    blocks = validate_blocks(BulkRnaNarrator().collect("bulk_rna_agent", agent_result))
    power = [block for block in blocks if block.id == "bulk.power"]

    html = render_blocks(power, report_dir=tmp_path)

    assert "31%" in html
    assert "71%" in html


def test_bulk_narrator_methods_are_generic_and_auditable():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    agent_result = {"status": "done", "findings": _bulk_findings()}
    methods = BulkRnaNarrator().methods("bulk_rna_agent", agent_result)
    assert len(methods) == 1
    assert "design ~condition" in methods[0]
    assert "adjusted p-value < 0.05" in methods[0]


def test_bulk_narrator_methods_disclose_outlier_sensitivity():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    findings = _bulk_findings()
    findings["outlier_sensitivity"] = {
        "status": "success",
        "removed_samples": ["ctrl_1"],
        "conclusion_robust": False,
    }
    agent_result = {"status": "done", "findings": findings}
    methods = BulkRnaNarrator().methods("bulk_rna_agent", agent_result)
    assert any("retained in the primary DE analysis" in line for line in methods)
    assert any("not robust" in line for line in methods)


def test_bulk_narrator_surfaces_gsea_as_narrative_block(tmp_path):
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
    from aria.agents.narrative.validators import validate_blocks

    gsea = tmp_path / "gsea_results.csv"
    gsea.write_text(
        ",nes,fdr\n"
        "pathway_alpha_response,1.9,0.04\n"
        "pathway_beta_response,-1.4,0.2\n"
        "pathway_gamma_response,1.1,0.4\n",
        encoding="utf-8",
    )
    findings = _bulk_findings()
    findings["contrasts"][0]["plots"] = {
        "gsea_table": str(gsea),
        "gsea_running_sums": [str(tmp_path / "running.png")],
    }

    agent_result = {"status": "done", "findings": findings}
    blocks = validate_blocks(BulkRnaNarrator().collect("bulk_rna_agent", agent_result))
    gsea_block = next(block for block in blocks if block.id == "bulk.gsea.treat_vs_ctrl")

    assert gsea_block.analysis == "gsea_preranked"
    assert gsea_block.metrics["n_pathways"] == 2
    assert gsea_block.metrics["top_pathways"][0]["term"] == "pathway_alpha_response"
    assert any(ev.label == "FDR<0.25 pathways" for ev in gsea_block.evidence)


def test_bulk_gsea_nonfinite_values_are_not_narrated_as_top_pathways(tmp_path):
    from aria.agents.narrative.compose_prose import compose_block_prose
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
    from aria.agents.narrative.validators import validate_blocks

    gsea = tmp_path / "gsea_results.csv"
    gsea.write_text(
        ",nes,fdr\n"
        "unstable_pathway,inf,0.0\n"
        "stable_pathway,1.9,0.04\n",
        encoding="utf-8",
    )
    findings = _bulk_findings()
    findings["contrasts"][0]["plots"] = {
        "gsea_table": str(gsea),
        "gsea_running_sums": [str(tmp_path / "running.png")],
    }

    blocks = validate_blocks(
        BulkRnaNarrator().collect(
            "bulk_rna_agent", {"status": "done", "findings": findings}
        )
    )
    gsea_block = next(block for block in blocks if block.id == "bulk.gsea.treat_vs_ctrl")
    prose = compose_block_prose(gsea_block)

    assert gsea_block.metrics["n_numeric_unstable"] == 1
    assert gsea_block.metrics["top_pathways"][0]["term"] == "stable_pathway"
    assert "unstable_pathway" not in prose
    assert "NES=inf" not in prose
    assert "FDR=0.0)" not in prose
    assert "non-finite" in prose
