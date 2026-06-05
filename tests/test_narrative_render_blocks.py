import json


def test_render_blocks_strict_false_withholds_bad_block_without_aborting():
    # Real-run bug: a single block failing strict verification aborted the WHOLE
    # report ("termina pero no genera el reporte html"). With strict=False the bad
    # block is withheld and every other block still renders.
    from aria.agents.narrative.render_blocks import render_blocks
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock

    good = NarrativeBlock(
        id="bulk.qc", modality="bulk RNA-seq", analysis="sample_qc",
        block_type="qc", title="QC", status="success", confidence="high",
        claim="Sample QC evaluated the inputs.",
        evidence=[EvidenceItem(label="samples", value=6)])
    bad = NarrativeBlock(
        id="bulk.gsea.bad", modality="bulk RNA-seq", analysis="gsea_preranked",
        block_type="exploratory", title="Bad GSEA", status="success",
        confidence="low", claim="GSEA had 999 unsupported signals.",
        evidence=[EvidenceItem(label="x", value=1)])
    html = render_blocks([good, bad], strict=False)
    assert "QC" in html                       # the good block rendered
    assert "withheld" in html                 # the bad block was withheld
    assert "999 unsupported" not in html      # its unverified prose is NOT shown


def test_gsea_prose_nes_matches_evidence_after_rounding_fix():
    # The GSEA composer rounded NES to :.3g (5.43) while the evidence held the raw
    # value (5.4321), so strict verification rejected the prose and killed the
    # report. The composer now uses the verifier's :.6g normalization.
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
    from aria.agents.narrative.compose_prose import compose_block_prose
    from aria.agents.narrative.evidence_verifier import verify_block_claim_support

    block = NarrativeBlock(
        id="bulk.gsea.X", modality="bulk RNA-seq", analysis="gsea_preranked",
        block_type="exploratory", title="Preranked GSEA X", status="success",
        confidence="low",
        claim="Preranked GSEA generated ranked pathway support for X.",
        evidence=[
            EvidenceItem(label="FDR<0.25 pathways", value=1, source="gsea_preranked"),
            EvidenceItem(label="GSEA term RNA Polymerase II",
                         value=5.43209876, source="gsea_preranked"),
        ],
        metrics={"n_pathways": 1, "top_pathways": [
            {"term": "RNA Polymerase II", "nes": 5.43209876, "fdr": 0.000312}]})
    prose = compose_block_prose(block)
    assert "NES=5.43" in prose                 # the value is shown
    verify_block_claim_support(block, prose, strict=True)   # and it verifies


def test_render_blocks_shows_claim_evidence_caveats_and_validates_files(tmp_path):
    import pytest
    from aria.agents.narrative.render_blocks import render_blocks
    from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock
    from aria.agents.narrative.validators import NarrativeValidationError

    fig = tmp_path / "figures" / "plot.png"
    table = tmp_path / "tables" / "de.tsv"
    fig.parent.mkdir()
    table.parent.mkdir()
    fig.write_bytes(b"png")
    table.write_text("gene\nGENE_UP_1\n", encoding="utf-8")

    block = NarrativeBlock(
        id="scrna.pseudobulk.GroupA.condition_a_vs_condition_b",
        modality="scRNA-seq",
        analysis="pseudobulk_de",
        block_type="result",
        title="GroupA condition_a_vs_condition_b",
        status="success",
        confidence="medium",
        claim="GroupA had 140 global-FDR DE genes.",
        evidence=[EvidenceItem("global-FDR DE genes", 140, "pseudobulk_de")],
        caveats=[Caveat("Composition covariate was included.", "info")],
        figures=[{"path": "figures/plot.png", "caption": "DE plot"}],
        tables=[{"path": "tables/de.tsv", "label": "DE genes"}],
    )

    html = render_blocks([block], report_dir=tmp_path)
    assert "data-block-id=\"scrna.pseudobulk.GroupA.condition_a_vs_condition_b\"" in html
    assert "GroupA condition_a_vs_condition_b contributed 140 global-FDR DE genes" in html
    assert "Structured evidence" in html
    assert "global-FDR DE genes" in html
    assert "Composition covariate" in html
    assert "DE genes" in html
    assert "data:image/png;base64" in html

    block.figures = [{"path": "figures/missing.png"}]
    with pytest.raises(NarrativeValidationError, match="referenced figure"):
        render_blocks([block], report_dir=tmp_path)


def test_render_blocks_fails_on_unsupported_claim_sentence():
    import pytest
    from aria.agents.narrative.claim_compiler import annotate_claim_tiers
    from aria.agents.narrative.render_blocks import render_blocks
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
    from aria.agents.narrative.validators import NarrativeValidationError

    block = NarrativeBlock(
        id="scrna.pseudobulk.Monocytes.STIM_vs_CTRL",
        modality="scRNA-seq",
        analysis="pseudobulk_de",
        block_type="result",
        title="Monocytes STIM_vs_CTRL",
        status="success",
        confidence="medium",
        claim="Monocytes had 140 global-FDR DE genes and IFNG drives response.",
        evidence=[EvidenceItem("global-FDR DE genes", 12, "pseudobulk_de")],
        metrics={"n_significant_global": 12},
    )
    annotate_claim_tiers([block], exp_ctx={})

    with pytest.raises(NarrativeValidationError, match="unsupported claim sentence"):
        render_blocks([block])


def test_render_blocks_stores_claim_verification_metadata():
    from aria.agents.narrative.render_blocks import render_blocks
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock

    block = NarrativeBlock(
        id="scrna.pseudobulk.Monocytes.STIM_vs_CTRL",
        modality="scRNA-seq",
        analysis="pseudobulk_de",
        block_type="result",
        title="Monocytes STIM_vs_CTRL",
        status="success",
        confidence="medium",
        claim="Monocytes had 12 global-FDR DE genes.",
        evidence=[EvidenceItem("global-FDR DE genes", 12, "pseudobulk_de")],
        metrics={"n_significant_global": 12},
    )

    render_blocks([block])
    verification = block.metadata["claim_verification"]
    assert verification["status"] == "supported"
    assert verification["evidence_card"]["evidence_card_id"] == (
        "scrna.pseudobulk.Monocytes.STIM_vs_CTRL#evidence"
    )


def test_render_blocks_does_not_duplicate_associative_badge_label():
    from aria.agents.narrative.render_blocks import render_blocks
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock

    block = NarrativeBlock(
        id="bulk.pathway.treat_vs_ctrl",
        modality="bulk RNA-seq",
        analysis="pathway_enrichment",
        block_type="result",
        title="Pathway enrichment treat_vs_ctrl",
        status="success",
        confidence="medium",
        claim="Pathway enrichment found 3 term(s) for treat_vs_ctrl.",
        evidence=[EvidenceItem("enriched terms", 3, "bulk_pathways")],
        metrics={"n_terms": 3},
        metadata={
            "claim": {
                "tier": "associative",
                "licensed_language": "associative",
                "rationale": "observational omics evidence",
            }
        },
    )

    html = render_blocks([block])

    assert "Evidence scope: association only" in html
    assert "associative · associative" not in html


def test_narrative_agent_composes_scrna_from_blocks_and_persists_json(tmp_path):
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    agent.memory = type("M", (), {"db_path": ":memory:"})()
    agent_results = {
        "scrna_agent": {
            "status": "done",
            "findings": {
                "pseudobulk_de": {
                    "groupby": "cluster",
                    "condition_col": "condition",
                    "n_groups": 1,
                    "per_group": {
                        "GroupA": {
                            "per_comparison": {
                                "condition_a_vs_condition_b": {
                                    "status": "success",
                                    "n_significant_global": 140,
                                    "n_significant_local": 180,
                                    "n_up_global": 80,
                                    "n_down_global": 60,
                                    "corrected_for_composition": True,
                                    "top_genes": [
                                        {"gene": "GENE_UP_1", "log2fc": 2.4},
                                        {"gene": "GENE_DOWN_1", "log2fc": -1.1},
                                    ],
                                }
                            }
                        }
                    },
                }
            },
        }
    }

    report = agent._render_html_report(
        experiment_id="exp_blocks",
        exp_ctx={"organism": "Homo sapiens", "genome": "GRCh38"},
        intent={"summary": "Compare condition A versus condition B."},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={
            "high": [],
            "medium": [],
            "low": [],
            "insufficient": [],
        },
        methods="methods",
        decisions=[],
        agent_results=agent_results,
        report_dir=tmp_path / "report",
    )
    html = report.read_text(encoding="utf-8")
    assert "data-block-id=\"scrna.pseudobulk.GroupA.condition_a_vs_condition_b\"" in html
    assert "GroupA condition_a_vs_condition_b contributed 140 global-FDR DE genes" in html
    assert "GENE_UP_1" in html
    assert "GENE_DOWN_1" in html

    methodology = json.loads((report.parent / "methodology.json").read_text())
    ids = [block["id"] for block in methodology["narrative_blocks"]]
    assert "scrna.pseudobulk.GroupA.condition_a_vs_condition_b" in ids


def test_rna_narrative_adapter_persists_input_hashes(tmp_path):
    from aria.scripts.rna_narrative_adapter import adapt

    h5ad = tmp_path / "synthetic.h5ad"
    h5ad.write_bytes(b"not really h5ad but hashable")
    report = {
        "data": [str(h5ad)],
        "organism": "Homo sapiens",
        "tissue_hint": "synthetic",
        "workspace": str(tmp_path),
        "mode": "pseudobulk",
        "stages": {},
        "pseudobulk_inputs": {
            "condition": "condition",
            "replicate": "sample_id",
            "groupby": "cluster",
            "comparisons": [["condition_a", "condition_b"]],
        },
    }

    bundle = adapt(report, workspace=tmp_path)
    inputs = bundle["exp_context"]["input_files"]

    assert inputs[0]["path"] == str(h5ad)
    assert inputs[0]["size_bytes"] == h5ad.stat().st_size
    assert len(inputs[0]["sha256"]) == 64
