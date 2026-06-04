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
