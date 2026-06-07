import pytest

pytest.importorskip("litellm")

from aria.agents.narrative_agent import NarrativeAgent


def _agent():
    return NarrativeAgent.__new__(NarrativeAgent)


def _bulk_agent_results():
    return {
        "bulk_rna_agent": {
            "status": "done",
            "findings": {
                "sample_qc": {
                    "n_samples": 9,
                    "outliers": [],
                },
                "design_used": "~ genotype",
                "padj_threshold": 0.05,
                "lfc_threshold": 0.58,
                "overlap": {
                    "BMAL1_KO vs WT intersect REV-ERBa_KO vs WT": {
                        "n_shared": 268,
                        "jaccard": 0.312,
                    }
                },
                "contrasts": [
                    {
                        "name": "BMAL1_KO vs WT",
                        "status": "success",
                        "n_significant": 481,
                        "n_upregulated": 225,
                        "n_downregulated": 256,
                        "power_estimate_at_lfc_min": 0.76,
                        "pathways": {
                            "GO_BP": [
                                {
                                    "term": "extracellular matrix organization",
                                    "adjusted_p": 0.001,
                                }
                            ]
                        },
                        "pathway_background": {
                            "background_size": 18791,
                            "background_source": "dataset_expressed_genes",
                        },
                        "pathway_ora": {
                            "method": "local_hypergeometric",
                            "gene_set_versions": {
                                "GO_BP": {
                                    "library": "GO_Biological_Process_2023",
                                    "release": "2023",
                                }
                            },
                        },
                    },
                    {
                        "name": "REV-ERBa_KO vs WT",
                        "status": "success",
                        "n_significant": 646,
                        "n_upregulated": 220,
                        "n_downregulated": 426,
                        "power_estimate_at_lfc_min": 0.7825,
                        "pathways": {},
                    },
                    {
                        "name": "BMAL1_KO vs REV-ERBa_KO",
                        "status": "success",
                        "n_significant": 198,
                        "n_upregulated": 183,
                        "n_downregulated": 15,
                        "power_estimate_at_lfc_min": 0.8669,
                        "pathways": {},
                    },
                ],
            },
        }
    }


def _exp_ctx():
    return {
        "organism": "Homo sapiens",
        "genome": "hg38",
        "user_question": (
            "H9 cells BMAL1 knockout, REV-ERBa knockout, and wildtype"
        ),
        "design": {
            "replicates": {
                "BMAL1_KO": 3,
                "REV-ERBa_KO": 3,
                "WT": 3,
            }
        },
    }


def test_bulk_executive_summary_is_deterministic_and_grounded():
    summary = _agent()._write_executive_summary(
        exp_ctx=_exp_ctx(),
        intent={"summary": "H9 knockout differential expression"},
        grouped={"high": [], "medium": [], "low": [], "insufficient": []},
        agent_results=_bulk_agent_results(),
    )

    assert "481 DE genes" in summary
    assert "646 DE genes" in summary
    assert "198 DE genes" in summary
    assert "BMAL1_KO n=3" in summary
    assert "absence of reported replicate counts" not in summary
    assert not summary.rstrip().endswith("ChIP-")
    assert "suppresses" not in summary


def test_bulk_methods_describe_local_ora_not_enrichr_endpoint():
    methods = _agent()._write_methods_section(
        exp_ctx=_exp_ctx(),
        agent_results=_bulk_agent_results(),
        decisions=[],
    )

    assert "local" in methods.lower()
    assert "hypergeometric" in methods.lower()
    assert "gene lists were not sent to enrichr" in methods.lower()
    assert "gseapy (Enrichr endpoint)" not in methods
    assert "instead of Enrichr's default universe" not in methods


def test_bulk_report_layout_and_single_modality_wording(tmp_path):
    agent = _agent()
    agent.reports_dir = tmp_path
    agent.memory = type("M", (), {"db_path": ":memory:"})()
    exp_ctx = _exp_ctx()
    agent_results = _bulk_agent_results()
    agent_results["raw_ingestion_agent"] = {
        "records": [{
            "mode": "fastq_kb_plan",
            "status": "blocked",
            "source_directory": "/data/raw_fastq",
            "reason": "raw_ingestion_kb_incomplete",
            "blockers": ["FASTQ ingestion requires explicit chemistry."],
            "missing_fields": ["chemistry", "index_path"],
        }]
    }

    report = agent._render_html_report(
        experiment_id="bulk_layout",
        exp_ctx=exp_ctx,
        intent={"summary": "H9 knockout differential expression"},
        executive_summary="grounded summary",
        findings_sections={
            "conflicts": agent._summarize_conflicts(agent_results, {
                "high": [], "medium": [], "low": [], "insufficient": [],
            })
        },
        grouped_findings={"high": [], "medium": [], "low": [], "insufficient": []},
        methods="methods",
        decisions=[],
        agent_results=agent_results,
        report_dir=tmp_path / "report",
    )
    html = report.read_text(encoding="utf-8")

    assert html.index("<h2>Executive Summary</h2>") < html.index("<h2>Provenance</h2>")
    assert "FASTQ-to-h5ad/kb ingestion routes" in html
    assert "STAR/featureCounts execution is reported separately" in html
    assert "raw_ingestion_kb_incomplete" in html
    assert "Missing fields: chemistry, index_path" in html
    assert "Cross-modal conflict analysis: not applicable; single-modality report." in html
    assert "No cross-modal conflicts identified." not in html
