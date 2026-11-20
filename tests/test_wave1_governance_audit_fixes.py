"""Wave 1 regression tests: C1-C4 governance, falsification, reporting, synthesis, adapters."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from aria.agents.narrative.claim_compiler import (
    annotate_claim_tiers,
    classify_claim,
    compile_public_claims,
)
from aria.agents.narrative.compose_prose import compose_block_prose
from aria.agents.narrative.evidence_verifier import (
    _as_nonnegative_number,
    _numbers_in,
    verify_block_claim_support,
)
from aria.agents.narrative.falsification_certificates import (
    build_falsification_certificates,
)
from aria.agents.narrative.hypothesis.adapters.hic import hic_evidence
from aria.agents.narrative.render_blocks import render_blocks
from aria.agents.narrative.reporting._base import _executive_summary_numbers
from aria.agents.narrative.reporting.render import _sanitize_json_value
from aria.agents.narrative.reporting.sections import ReportSectionsMixin
from aria.agents.narrative.run_ledger import _status_from_finding
from aria.agents.narrative.synthesis.cross_modal_synthesis import (
    compose_cross_modal_discussion_blocks,
    detect_cross_modal_patterns,
)
from aria.agents.narrative.synthesis.cross_modal_triangulation import (
    build_cross_modal_evidence_graph,
)
from aria.agents.narrative.types import EvidenceItem, NarrativeBlock


def test_scientific_notation_verification_and_gene_symbol_safety():
    assert "1.2e-09" in _numbers_in("FDR=1.2e-09")
    assert "1e-05" in _numbers_in("p = 1e-05")
    assert any(tok in _numbers_in("N = 1.2e5") for tok in ("1.2e+05", "1.2e5", "120000"))
    # Gene symbol E2F1 must not be tokenized as a number.
    assert _numbers_in("Top gene E2F1 was enriched.") == set()
    # Scientific notation in executive summary numbers.
    exec_nums = _executive_summary_numbers("Pathway reached p = 1.2e-09 across 4 samples.")
    assert "4" in exec_nums
    assert any("1.2e-09" in n or "1.2e-9" in n for n in exec_nums)

    block = NarrativeBlock(
        id="bulk.gsea.ko_vs_wt",
        modality="bulk_RNA",
        analysis="gsea_preranked",
        block_type="result",
        title="GSEA KO vs WT",
        status="success",
        confidence="medium",
        claim="Found 1 pathway: TP53 Regulates Transcription (NES=2.1, FDR=9.999999e-05).",
        evidence=[
            EvidenceItem(label="FDR<0.25 pathways", value=1, source="gsea"),
            EvidenceItem(
                label="top GSEA pathway",
                value="TP53 Regulates Transcription (NES=2.1, FDR=9.999999e-05)",
                source="gsea",
            ),
        ],
        metrics={
            "n_pathways": 1,
            "top_pathways": [
                {"term": "TP53 Regulates Transcription", "nes": 2.1, "fdr": 9.999999e-5}
            ],
        },
    )
    prose = compose_block_prose(block)
    res = verify_block_claim_support(block, prose, strict=True)
    assert res["status"] == "supported"


def test_as_nonnegative_number_rejects_bool_and_nonfinite():
    assert _as_nonnegative_number(True) is None
    assert _as_nonnegative_number(False) is None
    assert _as_nonnegative_number(float("inf")) is None
    assert _as_nonnegative_number(float("nan")) is None
    assert _as_nonnegative_number(-1.0) is None
    assert _as_nonnegative_number(0.0) == 0.0


def test_claim_compiler_da_causal_strong_mechanistic_and_nan_gate():
    exp_ctx = {
        "design": {
            "causal_estimands": [
                {
                    "id": "est_da",
                    "contrast": "KO_vs_WT",
                    "interventional": True,
                    "randomized": True,
                    "confounding_verified": True,
                }
            ]
        }
    }
    da_block = NarrativeBlock(
        id="chromatin.da.ko_vs_wt",
        modality="bulk_ATAC",
        analysis="differential_accessibility",
        block_type="result",
        title="Bulk ATAC DA KO_vs_WT",
        status="success",
        confidence="high",
        claim="KO_vs_WT showed 120 differentially accessible peaks.",
        evidence=[EvidenceItem(label="DA peaks", value=120, source="diffacc")],
        metrics={"n_significant": 120, "contrast": "KO_vs_WT"},
        metadata={"estimand": {"id": "est_da", "contrast": "KO_vs_WT"}},
    )
    annotated_da = annotate_claim_tiers([da_block], exp_ctx)
    assert annotated_da[0].metadata["claim"]["tier"] == "causal_experimental"

    b_de = NarrativeBlock(
        id="bulk.de.c1",
        modality="bulk_RNA",
        analysis="differential_expression",
        block_type="result",
        title="DE",
        status="success",
        confidence="high",
        claim="Found 50 DE genes.",
        evidence=[EvidenceItem(label="DE genes", value=50, source="deseq2")],
        metrics={"n_significant": 50},
        metadata={"subject": "shared_axis"},
    )
    b_m1 = NarrativeBlock(
        id="chrom.motif.c1",
        modality="bulk_ATAC",
        analysis="motif_enrichment",
        block_type="result",
        title="Motifs",
        status="success",
        confidence="medium",
        claim="Found 10 enriched motifs.",
        evidence=[EvidenceItem(label="enriched motifs", value=10, source="motifs")],
        metrics={"n_significant": 10},
        metadata={"subject": "shared_axis"},
    )
    b_m2 = NarrativeBlock(
        id="integ.p2g.c1",
        modality="integration",
        analysis="peak_to_gene",
        block_type="result",
        title="P2G",
        status="success",
        confidence="medium",
        claim="Found 25 cis links.",
        evidence=[EvidenceItem(label="cis links", value=25, source="p2g")],
        metrics={"n_links": 25},
        metadata={"subject": "shared_axis"},
    )
    annotated = annotate_claim_tiers([b_de, b_m1, b_m2], {})
    assert any(b.metadata["claim"]["tier"] == "strong_mechanistic" for b in annotated)

    nan_block = NarrativeBlock(
        id="bulk.de.nan",
        modality="bulk_RNA",
        analysis="differential_expression",
        block_type="result",
        title="NaN DE",
        status="success",
        confidence="high",
        claim="Found 0 genes.",
        evidence=[EvidenceItem(label="DE genes", value=0, source="deseq2")],
        metrics={"n_significant": float("nan"), "power_estimate": float("nan")},
    )
    nan_annotated = annotate_claim_tiers([nan_block], {})
    assert nan_annotated[0].metadata["claim"]["tier"] == "descriptive"


def test_claim_compiler_report_dir_withholds_missing_figure(tmp_path: Path):
    block = NarrativeBlock(
        id="bulk.de.missing_fig",
        modality="bulk_RNA",
        analysis="differential_expression",
        block_type="result",
        title="DE Missing Fig",
        status="success",
        confidence="high",
        claim="Found 10 DE genes.",
        evidence=[EvidenceItem(label="DE genes", value=10, source="deseq2")],
        metrics={"n_significant": 10},
        figures=[{"id": "volcano", "path": "figures/does_not_exist.png"}],
    )
    compiled = compile_public_claims([block], {}, report_dir=tmp_path)
    assert len(compiled.blocks) == 0
    assert len(compiled.withheld) == 1


def test_run_ledger_ran_false_precedes_success():
    status, reason = _status_from_finding(
        {"status": "success", "ran": False, "reason": "no_peaks"}
    )
    assert status == "skipped"
    assert reason == "no_peaks"


def test_falsification_certificates_cross_modal_and_replicates():
    claims = [
        {
            "claim_id": "integration.cross_modal_de_da.c1",
            "modality": "cross_modal",
            "analysis": "cross_modal_concordance",
            "tier": "associative",
        },
        {
            "claim_id": "integration.cross_modal_sample_links",
            "modality": "cross_modal",
            "analysis": "cross_modal_concordance",
            "tier": "associative",
        },
    ]
    agent_results_no_perm = {
        "cross_modal_concordance": {
            "status": "done",
            "findings": {"sample_matched_links": {"ran": False}},
        }
    }
    fc1 = build_falsification_certificates(
        claims, agent_results_no_perm, {"design": {"groups": {"A": ["s1"], "B": ["s2"]}}}
    )
    for cert in fc1["certificates"]:
        assert cert["block_preserving_permutation_null"]["status"] == "parametric_bh_only"
        assert cert["design_aliasing_and_unit_audit"]["low_power_n2_warning"] is True

    agent_results_with_perm = {
        "cross_modal_concordance": {
            "status": "done",
            "findings": {"sample_matched_links": {"ran": True}},
        }
    }
    fc2 = build_falsification_certificates(claims, agent_results_with_perm, {})
    by_id = {c["claim_id"]: c for c in fc2["certificates"]}
    assert (
        by_id["integration.cross_modal_de_da.c1"]["block_preserving_permutation_null"]["status"]
        == "parametric_bh_only"
    )
    assert (
        by_id["integration.cross_modal_sample_links"]["block_preserving_permutation_null"]["status"]
        == "calibrated"
    )


def test_sanitize_json_value_forbids_nan_and_inf():
    dirty = {"a": float("nan"), "b": [1.0, float("inf"), -float("inf")], "c": True}
    clean = _sanitize_json_value(dirty)
    encoded = json.dumps(clean, allow_nan=False)
    assert json.loads(encoded) == {"a": None, "b": [1.0, None, None], "c": True}


def test_render_blocks_relativizes_staged_table_links(tmp_path: Path):
    src_dir = tmp_path / "work"
    src_dir.mkdir()
    src_csv = src_dir / "cross_modal_de_da_concordance.csv"
    src_csv.write_text("gene,r\nG1,0.8\n", encoding="utf-8")

    report_dir = tmp_path / "report"
    (report_dir / "tables").mkdir(parents=True)
    (report_dir / "tables" / src_csv.name).write_text("gene,r\nG1,0.8\n", encoding="utf-8")

    block = NarrativeBlock(
        id="integration.cross_modal_de_da",
        modality="cross_modal",
        analysis="cross_modal_concordance",
        block_type="integration",
        title="Cross-Modal DE x DA",
        status="success",
        confidence="medium",
        claim="Evaluated 1 concordant gene G1.",
        evidence=[
            EvidenceItem(label="concordant genes", value=1, source="cm"),
            EvidenceItem(label="top concordant genes", value="G1", source="cm"),
        ],
        tables=[{"id": "de_da", "path": str(src_csv), "label": "DE x DA"}],
    )
    html_out = render_blocks([block], report_dir=report_dir, strict=True)
    assert "href='tables/cross_modal_de_da_concordance.csv'" in html_out
    assert str(src_dir) not in html_out


def test_sections_qc_guards_zero_peaks_nan_and_png_staging(tmp_path: Path):
    mixin = ReportSectionsMixin()
    grouped = {"high": [], "medium": [], "low": [], "insufficient": []}
    agent_results = {
        "bulk_rna_agent": {
            "status": "done",
            "findings": {
                "sample_qc": {
                    "n_samples": 6,
                    "size_ratio": float("nan"),
                    "pca_plot": str(tmp_path / "pca.png"),
                },
                "contrasts": [],
            },
        },
        "chromatin_agent": {
            "status": "done",
            "findings": {
                "bulk_ATAC": {
                    "findings": {
                        "peaks": {
                            "n_peaks": 0,
                            "frip": float("nan"),
                            "tss_enrichment": float("inf"),
                        }
                    }
                }
            },
        },
    }
    (tmp_path / "pca.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    html_qc = mixin._build_qc_section(grouped, agent_results=agent_results, exp_ctx={})
    assert "0 peaks" in html_qc
    assert "nan" not in html_qc.lower()
    assert "inf" not in html_qc.lower()

    report_dir = tmp_path / "rep"
    (report_dir / "figures").mkdir(parents=True)
    (report_dir / "tables").mkdir(parents=True)
    mixin._stage_artifacts(agent_results, report_dir)
    assert (report_dir / "figures" / "pca_all_samples.png").is_file()
    assert not (report_dir / "figures" / "pca_all_samples.svg").exists()


def test_compose_prose_never_interpolates_none_or_nan():
    b_de = NarrativeBlock(
        id="bulk.de.x",
        modality="bulk_RNA",
        analysis="differential_expression",
        block_type="result",
        title="DE",
        status="success",
        confidence="medium",
        claim="Adversarial hype claim must not leak.",
        evidence=[EvidenceItem(label="contrast", value="KO_vs_WT", source="deseq2")],
    )
    de_prose = compose_block_prose(b_de)
    assert "None" not in de_prose
    assert "Adversarial hype" not in de_prose

    b_pow = NarrativeBlock(
        id="bulk.power.x",
        modality="bulk_RNA",
        analysis="power",
        block_type="qc",
        title="Power",
        status="success",
        confidence="medium",
        claim="Power unavailable.",
        evidence=[EvidenceItem(label="status", value="unavailable", source="power")],
        metrics={"min_power": float("nan"), "max_power": True},
    )
    assert compose_block_prose(b_pow) == "Power unavailable."

    b_mofa = NarrativeBlock(
        id="integ.mofa.x",
        modality="integration",
        analysis="mofa",
        block_type="result",
        title="MOFA",
        status="success",
        confidence="medium",
        claim="MOFA fallback.",
        evidence=[EvidenceItem(label="MOFA cells evaluated", value=100, source="mofa")],
        metrics={
            "n_cells": 100,
            "n_factors": 5,
            "n_active_factors": 3,
            "n_shared_factors": 2,
            "n_modality_specific_factors": 1,
            "n_excluded_factors": 0,
        },
    )
    prose = compose_block_prose(b_mofa)
    assert "None" not in prose


def test_cross_modal_synthesis_perfect_concordance_and_held_out_fdr():
    agent_results = {
        "cross_modal_concordance": {
            "status": "done",
            "findings": {
                "ran": True,
                "rna_modality": "bulk_RNA",
                "atac_modality": "bulk_ATAC",
                "matched_contrasts": [
                    {
                        "contrast_id": "KO_vs_WT",
                        "de_da_summary": {
                            "n_pairs_tested": 20,
                            "n_both_significant": 10,
                            "n_concordant_total": 10,
                            "n_concordant_up": 6,
                            "n_concordant_down": 4,
                            "n_discordant_total": 0,
                            "fisher_odds_ratio": None,
                            "fisher_pvalue": 1.5e-6,
                        },
                    }
                ],
                "sample_matched_links": {
                    "ran": True,
                    "n_shared_samples": 8,
                    "n_pairs_tested": 50,
                    "n_significant_links": 4,
                    "n_positive_links": 4,
                    "n_negative_links": 0,
                    "permutation_null": {
                        "n_permutations": 100,
                        "held_out_discovery_fraction": 0.02,
                    },
                },
            },
        }
    }
    patterns = detect_cross_modal_patterns(agent_results)
    blocks = compose_cross_modal_discussion_blocks(patterns, agent_results)
    de_da = next(b for b in blocks if "cross_modal_de_da" in b.id)
    assert any(ev.label == "directional Fisher p-value" for ev in de_da.evidence)
    links = next(b for b in blocks if "cross_modal_sample_links" in b.id)
    assert any(ev.label == "empirical permutation FDR" for ev in links.evidence)


def test_cross_modal_triangulation_preserves_padj_mofa_views_and_wnn_axes():
    agent_results = {
        "integration_agent": {
            "status": "done",
            "findings": {
                "wnn": {
                    "status": "success",
                    "guardrail_status": "warn",
                    "n_discordant_clusters": 0,
                    "anti_degradation": {
                        "issues": [{"code": "cluster_modality_collapse", "severity": "warn"}],
                        "modality_preservation": {"evaluable": True},
                        "rare_population_retention": {"evaluable": True},
                        "spurious_bridge_screen": {"anchor_evaluable": True},
                        "per_cluster_weight_audit": {"evaluable": True},
                        "held_out_check": {"evaluable": False},
                    },
                },
                "mofa": {
                    "status": "success",
                    "views": {
                        "scRNA": {"kind": "rna"},
                        "scATAC": {"kind": "atac"},
                    },
                    "top_factors": [
                        {
                            "factor_id": 1,
                            "active": True,
                            "category": "shared",
                            "r2_by_view": {"scRNA": 0.14, "scATAC": 0.09},
                        }
                    ],
                },
                "peak_to_gene": {
                    "status": "success",
                    "top_links": [
                        {
                            "gene": "G1",
                            "peak": "chr1:100-500",
                            "corr": 0.72,
                            "padj": 1.2e-12,
                        }
                    ],
                },
            },
        },
        "chromatin_agent": {
            "status": "done",
            "findings": {
                "scATAC": {
                    "findings": {
                        "diffacc": {
                            "per_cluster": {
                                "ran": True,
                                "top_peaks_by_cluster": {
                                    "0": [{"peak": "chr1:100-500", "log2fc": None}]
                                },
                            }
                        }
                    }
                }
            },
        },
    }
    graph = build_cross_modal_evidence_graph(agent_results)
    nodes = {n["id"]: n for n in graph["nodes"]}
    guard = nodes["guardrail:wnn_anti_degradation"]
    assert guard["axis_statuses"]["modality_preservation"] == "pass"
    assert guard["axis_statuses"]["per_cluster_weight_audit"] == "warn"
    assert guard["axis_statuses"]["held_out_check"] == "not_evaluable"

    mofa_node = nodes["mofa:Factor_1"]
    assert mofa_node["rna_r2"] == 0.14
    assert mofa_node["atac_r2"] == 0.09

    p2g_edge = next(e for e in graph["edges"] if e["relation"] == "peak_to_gene_association")
    assert p2g_edge["padj"] == pytest.approx(1.2e-12)
    assert p2g_edge["padj"] > 0.0


def test_hic_adapter_preserves_small_qval_precision():
    agent_results = {
        "genome_arch_agent": {
            "status": "done",
            "findings": {
                "qc": {"n_samples": 3},
                "loops": {
                    "status": "success",
                    "ran": True,
                    "fdr_threshold": 0.05,
                    "top_loops": [
                        {
                            "chrom1": "chr1",
                            "start1": 1000,
                            "end1": 6000,
                            "chrom2": "chr1",
                            "start2": 50000,
                            "end2": 55000,
                            "qval": 2.5e-11,
                            "anchor1_genes": "MYC",
                        }
                    ],
                },
            },
        }
    }
    ledger = {"entries": [{"node_id": "ledger://hic/loops", "status": "ran"}]}
    signals = hic_evidence(agent_results, run_ledger=ledger)
    assert len(signals) == 1
    assert signals[0].value == pytest.approx(2.5e-11)
    assert signals[0].value > 0.0
