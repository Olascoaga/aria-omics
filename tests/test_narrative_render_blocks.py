import json


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
    table.write_text("gene\nISG15\n", encoding="utf-8")

    block = NarrativeBlock(
        id="scrna.pseudobulk.Monocytes.STIM_vs_CTRL",
        modality="scRNA-seq",
        analysis="pseudobulk_de",
        block_type="result",
        title="Monocytes STIM_vs_CTRL",
        status="success",
        confidence="medium",
        claim="Monocytes had 140 global-FDR DE genes.",
        evidence=[EvidenceItem("global-FDR DE genes", 140, "pseudobulk_de")],
        caveats=[Caveat("Composition covariate was included.", "info")],
        figures=[{"path": "figures/plot.png", "caption": "DE plot"}],
        tables=[{"path": "tables/de.tsv", "label": "DE genes"}],
    )

    html = render_blocks([block], report_dir=tmp_path)
    assert "data-block-id=\"scrna.pseudobulk.Monocytes.STIM_vs_CTRL\"" in html
    assert "Monocytes had 140 global-FDR DE genes" in html
    assert "global-FDR DE genes" in html
    assert "Composition covariate" in html
    assert "DE genes" in html
    assert "data:image/png;base64" in html

    block.figures = [{"path": "figures/missing.png"}]
    with pytest.raises(NarrativeValidationError, match="referenced figure"):
        render_blocks([block], report_dir=tmp_path)


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
                    "condition_col": "stim",
                    "n_groups": 1,
                    "per_group": {
                        "Monocytes": {
                            "per_comparison": {
                                "STIM_vs_CTRL": {
                                    "status": "success",
                                    "n_significant_global": 140,
                                    "n_significant_local": 180,
                                    "n_up_global": 80,
                                    "n_down_global": 60,
                                    "corrected_for_composition": True,
                                    "top_genes": [
                                        {"gene": "ISG15", "log2fc": 2.4},
                                        {"gene": "CCR2", "log2fc": -1.1},
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
        intent={"summary": "Compare STIM versus CTRL."},
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
    assert "data-block-id=\"scrna.pseudobulk.Monocytes.STIM_vs_CTRL\"" in html
    assert "Monocytes STIM_vs_CTRL had 140 global-FDR DE genes" in html
    assert "ISG15" in html
    assert "CCR2" in html

    methodology = json.loads((report.parent / "methodology.json").read_text())
    ids = [block["id"] for block in methodology["narrative_blocks"]]
    assert "scrna.pseudobulk.Monocytes.STIM_vs_CTRL" in ids


def test_rna_narrative_adapter_persists_input_hashes(tmp_path):
    from aria.scripts.rna_narrative_adapter import adapt

    h5ad = tmp_path / "pbmc.h5ad"
    h5ad.write_bytes(b"not really h5ad but hashable")
    report = {
        "data": [str(h5ad)],
        "organism": "Homo sapiens",
        "tissue_hint": "pbmc",
        "workspace": str(tmp_path),
        "mode": "pseudobulk",
        "stages": {},
        "pseudobulk_inputs": {
            "condition": "stim",
            "replicate": "Donor",
            "groupby": "cluster",
            "comparisons": [["STIM", "CTRL"]],
        },
    }

    bundle = adapt(report, workspace=tmp_path)
    inputs = bundle["exp_context"]["input_files"]

    assert inputs[0]["path"] == str(h5ad)
    assert inputs[0]["size_bytes"] == h5ad.stat().st_size
    assert len(inputs[0]["sha256"]) == 64
