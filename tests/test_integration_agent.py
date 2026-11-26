"""
ARIA IntegrationAgent Tests (native pytest — converted from legacy script harness).

Validates: strategy selection, WNN weight evaluation, peak-to-gene conflict
detection, MOFA+ no-hardcoded-biology checks, cross-modal conflict resolution,
script error contracts, and system prompt invariants.
"""

from __future__ import annotations

import pytest


class _MockLLM:
    def complete(self, *a, **kw):
        return (
            "ALTERNATIVE_HYPOTHESIS: Weight imbalance could reflect "
            "low ATAC data quality rather than biology\n"
            "CHALLENGES: Check FRiP score for ATAC\n"
            "EVIDENCE_REQUESTED: ATAC FRiP > 0.2 required\n"
            "VERDICT: ACCEPT_REVISED\n"
            "REVISED_CLAIM: WNN dominated by RNA (weight=0.91); "
            "ATAC contribution limited — verify ATAC QC before reporting"
        )

    def complete_heavy(self, *a, **kw):
        return self.complete()

    def complete_medium(self, *a, **kw):
        return self.complete()

    def complete_light(self, *a, **kw):
        return "ok"


class _MockMemory:
    def store_decision(self, *a, **kw):
        pass

    def create_tunnel(self, *a, **kw):
        pass

    def get_decisions(self, *a, **kw):
        return []

    def list_wings(self):
        return []


@pytest.fixture
def integration_agent():
    pytest.importorskip("litellm")
    from aria.agents.integration_agent import IntegrationAgent

    agent = IntegrationAgent.__new__(IntegrationAgent)
    agent.llm = _MockLLM()
    agent.memory = _MockMemory()
    agent._findings = []
    agent._escalations = []
    agent.publish_finding = lambda exp_id, payload, conf: agent._findings.append(
        {"payload": payload, "conf": conf}
    )
    agent.publish_escalation = lambda **kw: agent._escalations.append(kw)
    agent.publish_status = lambda *a, **kw: None
    return agent


def test_integration_agent_exposes_core_methods(integration_agent):
    for name in (
        "_select_strategy",
        "_run_wnn",
        "_run_peak_to_gene",
        "_run_mofa",
        "_resolve_conflicts",
    ):
        assert hasattr(integration_agent, name)


@pytest.mark.parametrize(
    "kwargs,expected",
    [
        (
            dict(
                has_scrna=True,
                has_scatac=True,
                has_hic=False,
                n_modalities=2,
                intent={"user_question": "What cell types are in my PBMC multiome?"},
            ),
            "wnn",
        ),
        (
            dict(
                has_scrna=True,
                has_scatac=True,
                has_hic=True,
                n_modalities=3,
                intent={"user_question": "What factors drive variation across modalities?"},
            ),
            "wnn_then_mofa",
        ),
        (
            dict(
                has_scrna=True,
                has_scatac=False,
                has_hic=False,
                n_modalities=2,
                intent={"user_question": "Run MOFA decomposition on RNA and ChIP"},
            ),
            "mofa",
        ),
        (
            dict(
                has_scrna=False,
                has_scatac=True,
                has_hic=False,
                n_modalities=2,
                intent={"user_question": "Link chromatin peaks to genes"},
            ),
            "peak_to_gene_only",
        ),
    ],
)
def test_integration_strategy_selection(integration_agent, kwargs, expected):
    assert integration_agent._select_strategy(**kwargs) == expected


def test_wnn_weight_evaluation_and_debate_downgrade(integration_agent):
    from aria.bus.message_bus import Confidence

    wnn_balanced = {
        "status": "success",
        "mean_rna_weight": 0.62,
        "mean_atac_weight": 0.38,
        "n_cells": 5000,
        "warnings": [],
    }
    integration_agent._evaluate_wnn_weights(
        "exp_001", wnn_balanced, {"user_question": "test"}
    )
    assert integration_agent._findings[-1]["conf"] == Confidence.HIGH

    integration_agent._findings = []
    wnn_imbalanced = {
        "status": "success",
        "mean_rna_weight": 0.91,
        "mean_atac_weight": 0.09,
        "n_cells": 5000,
        "atac_frip": 0.12,
        "warnings": ["RNA dominates WNN"],
    }
    integration_agent._evaluate_wnn_weights(
        "exp_002", wnn_imbalanced, {"user_question": "PBMC multiome analysis"}
    )
    assert integration_agent._findings[-1]["conf"] in (Confidence.LOW, Confidence.MEDIUM)


def test_resolve_conflicts_peak2gene_wnn_and_eregulon_flags(integration_agent):
    p2g_result = {
        "status": "success",
        "n_links": 847,
        "n_positive_correlations": 721,
        "n_negative_correlations": 126,
        "top_links": [
            {
                "gene": "PDCD1",
                "peak": "chr2:241800000-241802000",
                "correlation": -0.44,
                "distance_kb": 22.1,
                "direction": "negative",
            },
        ],
        "ctcf_validated": False,
        "hic_corroborated": False,
        "warnings": ["126 negative correlations detected"],
    }
    c_p2g = integration_agent._resolve_conflicts(
        "exp_003",
        {"genome": "hg38", "organism": "Homo sapiens"},
        {"user_question": "regulatory landscape of T cells"},
        {"peak_to_gene": p2g_result},
    )
    assert any(
        c.get("type") == "accessibility_expression_discordance"
        for c in c_p2g.get("conflicts", [])
    )

    c_wnn = integration_agent._resolve_conflicts(
        "exp_disc",
        {"genome": "hg38"},
        {"user_question": "test"},
        {
            "wnn": {
                "status": "success",
                "mean_rna_weight": 0.60,
                "mean_atac_weight": 0.40,
                "n_discordant_clusters": 3,
                "n_cells": 5000,
                "warnings": [],
            }
        },
    )
    assert any("discordance" in c.get("type", "") for c in c_wnn.get("conflicts", []))

    c_ereg = integration_agent._resolve_conflicts(
        "exp_ereg",
        {"genome": "hg38"},
        {"user_question": "eregulon analysis"},
        {
            "eregulons": {
                "status": "success",
                "n_eregulons": 12,
                "n_neg_tf_pos_p2g": 3,
                "n_paralog_ambiguous": 2,
                "n_depth_confounded": 1,
            }
        },
    )
    ctypes = {c.get("type") for c in c_ereg.get("conflicts", [])}
    assert {
        "eregulon_negative_tf_gene_discordance",
        "eregulon_paralog_motif_ambiguity",
        "eregulon_depth_confounded",
    } <= ctypes


def test_mofa_no_hardcoded_biology_and_factor_advice(integration_agent):
    import aria.scripts.integration_mofa as mofa
    from aria.scripts.integration_mofa import _mock_mofa

    assert not hasattr(mofa, "CELL_CYCLE_GENES")
    result = _mock_mofa(n_factors=10, reason="test")
    assert result["technical_factor_flag"] is False
    assert result["technical_factor_id"] is None
    assert all(str(x).startswith("feature_") for x in result["factor1_top_features"])
    assert not any("cell cycle" in str(w).lower() for w in result.get("warnings", []))

    mods_2 = {"scRNA": ["/data/rna.h5ad"], "scATAC": ["/data/atac.h5ad"]}
    mods_3 = {**mods_2, "HiC": ["/data/hic.mcool"]}
    n2 = integration_agent._advise_mofa_factors(
        "exp_2mod", {"user_question": "find shared variation in PBMC"}, mods_2
    )
    n3 = integration_agent._advise_mofa_factors(
        "exp_3mod", {"user_question": "atlas of complex tissue with many cell types"}, mods_3
    )
    assert n3 >= n2


def test_integration_scripts_fail_closed_on_missing_inputs():
    from aria.scripts.integration_eregulons import integration_eregulons
    from aria.scripts.integration_mofa import integration_mofa
    from aria.scripts.integration_peak2gene import integration_peak2gene
    from aria.scripts.integration_wnn import integration_wnn

    r_wnn = integration_wnn({"rna_files": ["/nonexistent/rna.h5ad"], "atac_files": ["/nonexistent/atac.h5ad"]})
    assert r_wnn["status"] == "error" and "Missing" in r_wnn["error_type"]

    r_p2g = integration_peak2gene({"rna_files": ["/nonexistent/rna.h5ad"], "atac_files": []})
    assert r_p2g["status"] == "error" and r_p2g["error_type"] == "MissingRNA"

    r_mofa = integration_mofa({"modalities": {"scRNA": ["/only/one"]}, "n_factors": 10})
    assert r_mofa["status"] == "error" and "Insufficient" in r_mofa["error_type"]

    r_ereg = integration_eregulons({"rna_files": ["/nonexistent/rna.h5ad"], "atac_files": []})
    assert r_ereg["status"] == "error" and r_ereg["error_type"] == "MissingRNA"


def test_integration_system_prompt_guardrails():
    pytest.importorskip("litellm")
    from aria.agents.integration_agent import INTEGRATION_SYSTEM

    low = INTEGRATION_SYSTEM.lower()
    assert "WNN" in INTEGRATION_SYSTEM and "MOFA" in INTEGRATION_SYSTEM and "peak" in low
    assert "negative correlation" in low and "poised" not in low
    assert "0.85" in INTEGRATION_SYSTEM or "imbalance" in low
    assert "feature names alone" in low and "cell cycle" not in low
    assert "independent" in low
