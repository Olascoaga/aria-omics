"""Wave 5A regression suite: Cross-Modal, Core Stats, Governance C1-C4, Capsules, and Top-Level Agents."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import fisher_exact, pearsonr, spearmanr

from aria.agents.data_audit_agent import DataAuditAgent
from aria.agents.narrative.claim_compiler import (
    annotate_claim_tiers,
    compile_public_claims,
)
from aria.agents.narrative.compose_prose import compose_block_prose
from aria.agents.narrative.evidence_verifier import (
    EvidenceVerificationError,
    build_evidence_card,
    verify_block_claim_support,
)
from aria.agents.narrative.falsification_certificates import (
    build_falsification_certificates,
)
from aria.agents.narrative.hypothesis import (
    DiscriminatingExperiment,
    EvidenceSignal,
    Hypothesis,
    verify_hypothesis_grounding,
)
from aria.agents.narrative.ledger_export import (
    _portable_methodology,
    verify_reproducible_capsule,
    write_reproducible_capsule,
)
from aria.agents.narrative.render_blocks import render_blocks
from aria.agents.narrative.report_builder import ReportBuilderMixin
from aria.agents.narrative.run_ledger import (
    LedgerLinkageError,
    _entries_for,
    _status_from_finding,
    link_claims_to_ledger,
    verify_blocks_against_ledger,
)
from aria.agents.narrative.synthesis.cross_modal_synthesis import (
    compose_cross_modal_discussion_blocks,
    detect_cross_modal_patterns,
)
from aria.agents.narrative.synthesis.cross_modal_triangulation import (
    build_cross_modal_evidence_graph,
)
from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
from aria.agents.narrative.validators import (
    _unnegated_causal_hit,
    find_causal_language,
    validate_block,
)
from aria.headless import _find_report_on_disk, default_answer_policy
from aria.utils import cross_modal_concordance as cm
from aria.utils.decision_manifest import (
    build_decision_manifest,
    commit_decision_manifest,
    verify_decision_manifest,
)
from aria.utils.design_power import _assess_contrast, assess_design_power
from aria.utils.lobo_jackknife import _unpaired_fold_stats, compute_lobo_jackknife
from aria.utils.sensitivity import classify_sensitivity


# ── 1. TRI-01: Empirical FDR vs held-out per-pair FPR ────────────────────────


def test_tri01_sample_matched_empirical_fdr_and_honest_synthesis_label(tmp_path):
    rng = np.random.default_rng(7)
    n = 8
    cols = [f"D{i}" for i in range(n)]
    x = rng.normal(5.0, 1.5, size=(1, n))
    y = x + rng.normal(0.0, 0.005, size=(1, n))
    rna = pd.DataFrame(np.exp2(x), index=["G1"], columns=cols)
    atac = pd.DataFrame(np.exp2(y), index=["chr1:100-200"], columns=cols)
    models = {"chr1": [{"gene": "G1", "tss": 150, "strand": "+"}]}
    params = {
        "cis_distance_bp": 50000,
        "min_matched_samples": 6,
        "seed": 3,
        "n_permutations": 100,
        "padj_threshold": 0.05,
    }
    _, res = cm.sample_matched_layer(rna, atac, models, params, [])
    perm = res["permutation_null"]
    assert res["n_significant_links"] == 1
    expected_efdr = min(1.0, float(perm["held_out_n_significant"]) / 1.0)
    assert perm["empirical_fdr_at_padj"] == pytest.approx(expected_efdr)

    # Zero significant links -> empirical_fdr_at_padj must be None (not 0/0 or FPR).
    _, res_zero = cm.sample_matched_layer(
        rna, atac, models, {**params, "padj_threshold": 0.0}, []
    )
    assert res_zero["n_significant_links"] == 0
    assert res_zero["permutation_null"]["empirical_fdr_at_padj"] is None

    # Synthesis must NEVER call held_out_discovery_fraction "empirical FDR".
    base_ar = {
        "cross_modal_concordance": {
            "status": "done",
            "findings": {
                "ran": True,
                "rna_modality": "bulk_RNA",
                "atac_modality": "bulk_ATAC",
                "matched_contrasts": [],
                "sample_matched_links": {
                    "ran": True,
                    "n_shared_samples": 8,
                    "n_pairs_tested": 500,
                    "n_significant_links": 10,
                    "n_positive_links": 8,
                    "n_negative_links": 2,
                    "permutation_null": {
                        "n_permutations": 100,
                        "empirical_fdr_at_padj": None,
                        "held_out_discovery_fraction": 0.004,
                    },
                },
            },
        }
    }
    blocks = compose_cross_modal_discussion_blocks(detect_cross_modal_patterns(base_ar), base_ar)
    link_block = next(b for b in blocks if b.id == "integration.cross_modal_sample_links")
    assert "empirical FDR" not in link_block.claim
    assert "held-out null discovery rate (per-pair FPR)" in link_block.claim
    assert any(
        ev.label == "held-out null discovery rate (per-pair FPR)"
        for ev in link_block.evidence
    )
    assert not any(ev.label == "empirical permutation FDR" for ev in link_block.evidence)
    assert verify_block_claim_support(link_block, compose_block_prose(link_block), strict=True)["status"] == "supported"
    assert find_causal_language(link_block.claim) is None

    # When empirical_fdr_at_padj IS finite, report it as empirical FDR.
    base_ar["cross_modal_concordance"]["findings"]["sample_matched_links"]["permutation_null"]["empirical_fdr_at_padj"] = 0.05
    blocks2 = compose_cross_modal_discussion_blocks(detect_cross_modal_patterns(base_ar), base_ar)
    link_block2 = next(b for b in blocks2 if b.id == "integration.cross_modal_sample_links")
    assert "empirical FDR of 0.05" in link_block2.claim
    assert any(ev.label == "empirical permutation FDR" for ev in link_block2.evidence)


# ── 2. CLD-01 & D5-H9: Gene-level DE×DA stats, depleted motifs, Panel C sort ─


def test_cld01_d5h9_gene_level_de_da_depleted_motifs_and_panel_c_sort(tmp_path):
    # 1 gene (G_MULTI) has 20 concordant_up peaks; 4 genes have 1 peak each (3 discordant, 1 concordant_down).
    rows = []
    for k in range(20):
        rows.append({
            "contrast_id": "T_vs_C",
            "gene": "G_MULTI",
            "gene_key": "G_MULTI",
            "peak": f"chr1:{1000 + k * 50}-{1020 + k * 50}",
            "distance_to_tss": k * 10,
            "rna_log2fc": 3.0,
            "rna_padj": 0.001,
            "atac_log2fc": 2.5 - k * 0.01,
            "atac_padj": 0.001 + k * 0.0001,
            "classification": "concordant_up",
            "rna_significant": True,
            "atac_significant": True,
            "best_cell_group": None,
        })
    extras = [
        ("G_CD", -2.0, -2.0, "concordant_down"),
        ("G_D1", 2.0, -1.5, "discordant_rna_up_atac_down"),
        ("G_D2", 1.5, -2.5, "discordant_rna_up_atac_down"),
        ("G_D3", -1.5, 2.0, "discordant_rna_down_atac_up"),
    ]
    for i, (g, rlfc, alfc, cls) in enumerate(extras):
        rows.append({
            "contrast_id": "T_vs_C",
            "gene": g,
            "gene_key": g,
            "peak": f"chr2:{100 + i * 1000}-{200 + i * 1000}",
            "distance_to_tss": 50,
            "rna_log2fc": rlfc,
            "rna_padj": 0.002,
            "atac_log2fc": alfc,
            "atac_padj": 0.002,
            "classification": cls,
            "rna_significant": True,
            "atac_significant": True,
            "best_cell_group": None,
        })
    frame = pd.DataFrame(rows)
    summary = cm.summarize_de_da(frame, promoter_bp=2000)

    # Pair counts stay backwards-compatible; statistical inference is strictly gene-level (5 genes, NOT 24 pairs).
    assert summary["n_pairs_tested"] == 24
    assert summary["n_gene_peak_pairs_tested"] == 24
    assert summary["n_genes_tested"] == 5
    # 2 concordant genes out of 5 both-significant genes -> 0.4, NOT 21/24 = 0.875!
    assert summary["concordance_fraction"] == pytest.approx(2 / 5)
    exp_or, exp_p = fisher_exact([[1, 2], [1, 1]], alternative="two-sided")
    assert summary["fisher_odds_ratio"] == pytest.approx(exp_or)
    assert summary["fisher_pvalue"] == pytest.approx(exp_p)

    # Gene-level primary Spearman/Pearson correlations (5 points, not 24 pseudoreplicated points)
    primary = frame.drop_duplicates("gene_key")
    assert summary["spearman_r"] == pytest.approx(
        spearmanr(primary["rna_log2fc"], primary["atac_log2fc"]).statistic
    )
    assert summary["pearson_r"] == pytest.approx(
        pearsonr(primary["rna_log2fc"], primary["atac_log2fc"]).statistic
    )
    assert summary["promoter_spearman_r"] == pytest.approx(summary["spearman_r"])

    # D5-H9: depleted motifs (log2_enrichment <= 0) must be skipped, not wrapped with abs().
    rna_df, _ = cm.load_rna_table(
        {"all_sig": [{"gene": "DEPLETED_TF", "log2fc": 2.0, "padj": 0.001},
                     {"gene": "ENRICHED_TF", "log2fc": 2.0, "padj": 0.001}]},
        {},
        [],
    )
    atac_payload = {
        "motifs": {
            "per_group": {
                "up": {
                    "top_motifs": [
                        {"name": "DEPLETED_TF", "padj": 0.001, "log2_enrichment": -2.5},
                        {"name": "ZERO_TF", "padj": 0.001, "log2_enrichment": 0.0},
                        {"name": "ENRICHED_TF", "padj": 0.001, "log2_enrichment": 1.8},
                    ]
                }
            }
        }
    }
    tf_frame, _ = cm.tf_layer(atac_payload, {"name": "A_vs_B"}, rna_df, False, 0.05, [])
    assert tf_frame["tf"].tolist() == ["ENRICHED_TF"]
    assert tf_frame["chromatin_score"].iloc[0] == pytest.approx(1.8)

    # Panel C uses _pathway_sort_key so zero-padj (0, 0.049) ranks conservatively after tiny (1e-200, 1e-200).
    p_rows = pd.DataFrame([
        {"safe_label": "ZeroPrecision", "term_key": "zeroprecision", "rna_padj": 0.0, "atac_padj": 0.049, "geometric_mean_padj": 0.0},
        {"safe_label": "TinyBoth", "term_key": "tinyboth", "rna_padj": 1e-200, "atac_padj": 1e-200, "geometric_mean_padj": 1e-200},
    ])
    _, svg_path = cm.render_starburst(pd.DataFrame(), pd.DataFrame(), p_rows, {}, tmp_path)
    svg_text = Path(svg_path).read_text(encoding="utf-8")
    # Both rendered; top bar (last in reversed iloc[::-1] / higher y) is TinyBoth.
    assert "TinyBoth" in svg_text and "ZeroPrecision" in svg_text


# ── 3. D5-C3 & D5-H7: Triangulation motif names, DA comparisons, MOFA quarantine, distance_bp ──


def test_d5c3_d5h7_triangulation_da_motifs_mofa_quarantine_and_distance_guard():
    agent_results = {
        "integration_agent": {
            "status": "done",
            "findings": {
                "mofa": {
                    "status": "success",
                    "batch_factor_ids": [2],
                    "cell_cycle_factor_ids": ["Factor_3"],
                    "technical_factor_id": 4,
                    "top_factors": [
                        {"factor_id": 1, "category": "shared", "rna_r2": 0.12, "atac_r2": 0.08},
                        {"factor_id": 2, "category": "shared", "rna_r2": 0.20, "atac_r2": 0.15},
                        {"factor_id": 3, "category": "shared", "rna_r2": 0.10, "atac_r2": 0.09},
                        {"factor_id": 4, "category": "shared", "rna_r2": 0.11, "atac_r2": 0.07},
                        {"factor_id": 5, "category": "shared", "batch_flag": True, "rna_r2": 0.09, "atac_r2": 0.06},
                        {"factor_id": 6, "category": "shared", "technical_depth_flag": True, "rna_r2": 0.08, "atac_r2": 0.05},
                        {"factor_id": 7, "category": "shared", "cell_cycle_flag": True, "rna_r2": 0.07, "atac_r2": 0.04},
                        {"factor_id": 8, "category": "shared", "confounder_flag": True, "rna_r2": 0.06, "atac_r2": 0.03},
                    ],
                    "excluded_factors": [
                        {"factor_id": 2, "exclusion_reasons": ["batch_r2>=0.70"]},
                    ],
                },
                "peak_to_gene": {
                    "status": "success",
                    "top_links": [
                        {"gene": "SOX10", "peak": "chr1:100-200", "corr": 0.5, "padj": 0.001, "distance_bp": float("nan")},
                        {"gene": "MBP", "peak": "chr1:300-400", "corr": 0.4, "padj": 0.002, "distance_bp": True},
                        {"gene": "PLP1", "peak": "chr1:500-600", "corr": 0.45, "padj": 0.002, "distance_bp": float("inf")},
                    ],
                },
                "eregulons": {
                    "status": "success",
                    "top_eregulons": [
                        {"tf": "SOX10", "sign_class": "pos_tf_pos_p2g", "top_target_genes": ["MBP"]},
                    ],
                },
            },
        },
        "chromatin_agent": {
            "status": "done",
            "findings": {
                "differential_accessibility": {
                    "pseudobulk": {
                        "ran": True,
                        "comparisons": [
                            {"status": "success", "top_peaks": [{"peak": "chr1:100-200", "log2fc": 1.8}]},
                        ],
                    },
                    "comparisons": [
                        {"status": "success", "top_peaks": [{"peak": "chr1:300-400", "log2FoldChange": -1.4}]},
                    ],
                },
                "motifs": {
                    "ran": True,
                    "per_group": {
                        "up": {"top_motifs": [{"name": "SOX10(HMG)", "padj": 0.001}]}
                    },
                },
            },
        },
    }
    graph = build_cross_modal_evidence_graph(agent_results)
    nodes = {n["id"]: n for n in graph["nodes"]}

    # Only Factor 1 is active/shared; Factors 2..8 are ALL quarantined as excluded_confounder.
    assert graph["mofa_summary"]["n_shared_factors"] == 1
    assert graph["mofa_summary"]["n_excluded_confounder_factors"] == 7
    assert nodes["mofa:Factor_1"]["status"] == "active"
    for fid in range(2, 9):
        assert nodes[f"mofa:Factor_{fid}"]["status"] == "excluded_confounder"
        assert nodes[f"mofa:Factor_{fid}"]["classification"] == "excluded_confounder"

    # Non-finite / bool distance_bp guarded to None cleanly (strict JSON serializable).
    assert nodes["peak:chr1:100-200"]["distance_bp"] is None
    assert nodes["peak:chr1:300-400"]["distance_bp"] is None
    assert nodes["peak:chr1:500-600"]["distance_bp"] is None
    json.dumps(graph, allow_nan=False)

    # Both scATAC pseudobulk.comparisons and top-level bulk ATAC comparisons produced DA edges.
    da_targets = {e["target"] for e in graph["edges"] if e["relation"] == "atac_da_association"}
    assert {"peak:chr1:100-200", "peak:chr1:300-400"} <= da_targets

    # Motif row with `name` (instead of `tf`/`motif_name`) produced atac_motif_enrichment edge.
    motif_targets = {e["target"] for e in graph["edges"] if e["relation"] == "atac_motif_enrichment"}
    assert "tf:SOX10" in motif_targets


# ── 4. TRI-02 / D5-H4: Falsification permutation axis verifies real calibration ─


def test_tri02_d5h4_falsification_permutation_axis_no_unconditional_calibrated():
    claims = [
        {"claim_id": "c.fp", "modality": "chromatin", "analysis": "differential_tf_footprinting", "tier": "associative"},
        {"claim_id": "c.p2g", "modality": "integration", "analysis": "peak_to_gene", "tier": "associative"},
        {"claim_id": "c.ereg", "modality": "integration", "analysis": "eregulons", "tier": "associative"},
    ]

    # Empty agent_results -> NONE may be marked "calibrated"!
    fc_empty = build_falsification_certificates(claims, {}, {})
    by_id = {c["claim_id"]: c["block_preserving_permutation_null"] for c in fc_empty["certificates"]}
    assert by_id["c.fp"]["status"] == "not_calibrated_descriptive_only"
    assert by_id["c.p2g"]["status"] == "not_calibrated"
    assert by_id["c.ereg"]["status"] == "not_calibrated"

    # Uncalibrated / descriptive / zero-permutation runs -> still NOT calibrated.
    uncal_results = {
        "chromatin_agent": {
            "findings": {
                "footprinting": {
                    "ran": True,
                    "status": "success",
                    "differential_summary": {
                        "inference": {"status": "success", "null_label_controls": {"n_permutations": 0}},
                    },
                }
            }
        },
        "integration_agent": {
            "findings": {
                "peak_to_gene": {"ran": True, "status": "success", "background_method": "none"},
                "eregulons": {"ran": True, "status": "success", "n_null_permutations": 0},
            }
        },
    }
    fc_uncal = build_falsification_certificates(claims, uncal_results, {})
    by_id_uncal = {c["claim_id"]: c["block_preserving_permutation_null"] for c in fc_uncal["certificates"]}
    assert by_id_uncal["c.fp"]["status"] == "parametric_bh_only"
    assert by_id_uncal["c.p2g"]["status"] == "not_calibrated"
    assert by_id_uncal["c.ereg"]["status"] == "not_calibrated"

    # Genuinely calibrated runs -> "calibrated".
    cal_results = {
        "chromatin_agent": {
            "findings": {
                "footprinting": {
                    "ran": True,
                    "status": "success",
                    "differential_summary": {
                        "inference": {"status": "success", "null_label_controls": {"n_permutations": 50}},
                    },
                }
            }
        },
        "integration_agent": {
            "findings": {
                "peak_to_gene": {
                    "ran": True,
                    "status": "success",
                    "background_method": "trans_chromosome_gc_accessibility",
                },
                "eregulons": {
                    "ran": True,
                    "status": "success",
                    "n_null_permutations": 500,
                    "p2g_summary": {"background_match": "trans_gc_acc"},
                },
            }
        },
    }
    fc_cal = build_falsification_certificates(claims, cal_results, {})
    by_id_cal = {c["claim_id"]: c["block_preserving_permutation_null"] for c in fc_cal["certificates"]}
    assert by_id_cal["c.fp"]["status"] == "calibrated"
    assert by_id_cal["c.p2g"]["status"] == "calibrated"
    assert by_id_cal["c.ereg"]["status"] == "calibrated"


# ── 5. BI-11 / D5-H3 / CDX-03 / D5-C2: C4 bypass removal, negation, pseudobulk ─


def test_bi11_d5h3_no_causal_evidence_bypass_multiclause_negation_and_idempotency():
    # Forging metadata["causal_evidence"] = True without a C4 estimand contract must NOT bypass!
    forged = NarrativeBlock(
        id="bulk.forged",
        modality="bulk_rna",
        analysis="differential_expression",
        block_type="result",
        title="Forged causal",
        status="success",
        confidence="high",
        claim="GeneX drives 10 downstream changes.",
        evidence=[
            EvidenceItem("significant features", 10, "deseq2"),
            EvidenceItem("lead gene", "GeneX", "deseq2"),
        ],
        metrics={"n_significant": 10},
        metadata={"causal_evidence": True},
    )
    ledger = {
        "entries": [{
            "node_id": "ledger://bulk/differential_expression",
            "node_type": "analysis_run",
            "provenance": {"kind": "structured_agent_result"},
            "status": "ran",
        }]
    }
    compiled = compile_public_claims([forged], exp_ctx={}, run_ledger=ledger)
    assert compiled.blocks == []
    assert len(compiled.withheld) == 1
    assert forged.metadata["claim"]["language_violation"] == "drives"

    # Multi-clause sentence: negated first clause + unnegated second clause MUST be caught!
    assert _unnegated_causal_hit("Factor A does not cause X; Factor B drives Z.") == "drives"
    # Properly negated verbs must NOT be false positives!
    assert _unnegated_causal_hit("Factor A never drives X and cannot cause Y.") is None

    # evidence_verifier catches multi-clause causal smuggling and accepts "never drives"/"cannot cause"
    b_multi = NarrativeBlock(
        id="bulk.multi",
        modality="bulk_rna",
        analysis="custom",
        block_type="result",
        title="Multi-clause",
        status="success",
        confidence="high",
        claim="This signal does not cause change; it drives progression.",
        evidence=[EvidenceItem("signal", "present", "test")],
        metadata={"claim": {"tier": "associative", "licensed_language": "associative"}},
    )
    with pytest.raises(EvidenceVerificationError, match="causal language"):
        verify_block_claim_support(b_multi, strict=True)

    b_neg = NarrativeBlock(
        id="bulk.neg",
        modality="bulk_rna",
        analysis="custom",
        block_type="result",
        title="Negated causal",
        status="success",
        confidence="high",
        claim="This signal never drives change and cannot cause progression.",
        evidence=[EvidenceItem("signal", "present", "test")],
        metadata={"claim": {"tier": "associative", "licensed_language": "associative"}},
    )
    assert verify_block_claim_support(b_neg, strict=True)["status"] == "supported"

    # validators._apply_causal_guard is idempotent (no double-downgrade) & respects licensed_language=="causal"
    b_idem = NarrativeBlock(
        id="bulk.idem",
        modality="bulk_rna",
        analysis="freeform",
        block_type="result",
        title="Idempotent guard",
        status="success",
        confidence="high",
        claim="Module drives proliferation.",
        evidence=[EvidenceItem("n", 5, "test")],
    )
    validate_block(b_idem)
    assert b_idem.confidence == "medium"
    validate_block(b_idem)
    assert b_idem.confidence == "medium"  # NOT downgraded to "low" on second pass!

    b_licensed = NarrativeBlock(
        id="bulk.licensed",
        modality="bulk_rna",
        analysis="differential_expression",
        block_type="result",
        title="Licensed causal",
        status="success",
        confidence="high",
        claim="Drug treatment causes 10 DE genes.",
        evidence=[EvidenceItem("DE genes", 10, "deseq2")],
        metrics={"n_significant": 10},
        metadata={"estimand": {"id": "drug", "contrast": "drug_vs_ctrl"}},
    )
    ctx_c4 = {
        "design": {
            "causal_estimands": [{
                "id": "drug",
                "contrast": "drug_vs_ctrl",
                "interventional": True,
                "randomized": True,
                "confounding_verified": True,
            }]
        }
    }
    annotate_claim_tiers([b_licensed], ctx_c4)
    assert b_licensed.metadata["claim"]["licensed_language"] == "causal"
    validate_block(b_licensed)
    assert b_licensed.confidence == "high"
    html_out = render_blocks([b_licensed], strict=True)
    assert "correlation, not causation" not in html_out


def test_cdx03_d5c2_pseudobulk_prose_and_semantic_facts_align_with_fdr_strategy():
    # Case 1: fdr_strategy == "per_cluster", global == 0, local == 42 (25 up, 17 down)
    b_local = NarrativeBlock(
        id="scrna.pseudobulk.T.drug_vs_ctrl",
        modality="scRNA-seq",
        analysis="pseudobulk_de",
        block_type="result",
        title="T drug_vs_ctrl",
        status="success",
        confidence="medium",
        claim="T drug_vs_ctrl had 42 per-cluster FDR DE genes.",
        evidence=[
            EvidenceItem("global-FDR DE genes", 0, "pseudobulk_de"),
            EvidenceItem("local-FDR DE genes", 42, "pseudobulk_de"),
            EvidenceItem("up genes", 25, "pseudobulk_de"),
            EvidenceItem("down genes", 17, "pseudobulk_de"),
            EvidenceItem("top gene GeneAlpha", 2.1, "pseudobulk_de"),
        ],
        metrics={
            "fdr_strategy": "per_cluster",
            "n_significant": 42,
            "n_significant_global": 0,
            "n_significant_local": 42,
        },
    )
    prose_local = compose_block_prose(b_local)
    assert "0 global-FDR DE genes (25 higher, 17 lower)" not in prose_local
    assert "42 per-cluster local-FDR DE genes (25 higher, 17 lower)" in prose_local
    card_local = build_evidence_card(b_local)
    star_polarities = {f.polarity for f in card_local.semantic_facts if f.subject == "*"}
    assert star_polarities == {"affirmed"}
    assert verify_block_claim_support(b_local, prose_local, strict=True)["status"] == "supported"

    # Case 2: fdr_strategy == "global", global == 0, local == 42 (25 up, 17 down)
    b_global_zero = NarrativeBlock(
        id="scrna.pseudobulk.B.drug_vs_ctrl",
        modality="scRNA-seq",
        analysis="pseudobulk_de",
        block_type="result",
        title="B drug_vs_ctrl",
        status="success",
        confidence="medium",
        claim="No significant genes were detected under global FDR.",
        evidence=[
            EvidenceItem("global-FDR DE genes", 0, "pseudobulk_de"),
            EvidenceItem("local-FDR DE genes", 42, "pseudobulk_de"),
            EvidenceItem("up genes", 25, "pseudobulk_de"),
            EvidenceItem("down genes", 17, "pseudobulk_de"),
        ],
        metrics={
            "fdr_strategy": "global",
            "n_significant": 0,
            "n_significant_global": 0,
            "n_significant_local": 42,
        },
    )
    prose_global = compose_block_prose(b_global_zero)
    assert "0 global-FDR DE genes (25 higher, 17 lower)" not in prose_global
    assert "0 global-FDR DE genes; 42 passed the per-block local FDR threshold (25 higher, 17 lower)" in prose_global
    card_global = build_evidence_card(b_global_zero)
    star_polarities_g = {f.polarity for f in card_global.semantic_facts if f.subject == "*"}
    assert star_polarities_g == {"negated"}
    assert verify_block_claim_support(b_global_zero, prose_global, strict=True)["status"] == "supported"


# ── 6. BI-15 / D5-H5: RunLedger status allowlist & multi-key reconciliation ──


def test_bi15_d5h5_run_ledger_status_allowlist_and_multi_key_preference():
    assert _status_from_finding({"status": "partial"})[0] == "not_run"
    assert _status_from_finding({"status": "timeout"})[0] == "error"
    assert _status_from_finding({"status": "cancelled"})[0] == "skipped"
    assert _status_from_finding({"status": "missing"})[0] == "not_run"
    assert _status_from_finding({"status": "custom_unknown"})[0] == "not_run"
    assert _status_from_finding({"status": "success"})[0] == "ran"

    # _entries_for prefers a ran sibling key over an earlier skipped/error key
    specs = [{"key": "pathway_enrichment", "label": "Pathways", "plan_kw": ["ora"],
              "finding_keys": ["pseudobulk_pathways", "pathways"]}]
    entries = _entries_for(
        "scRNA",
        specs,
        {"pseudobulk_pathways": {"status": "skipped", "reason": "no pb"}, "pathways": {"status": "success"}},
        [],
    )
    assert entries[0]["status"] == "ran"

    # Allowlist rejects "partial" / "timeout" / None across link_claims, verify_blocks, and hypothesis grounding
    for bad_status in ("partial", "timeout", "unknown_state", None):
        bad_ledger = {
            "entries": [{
                "node_id": "ledger://bulk/differential_expression",
                "node_type": "analysis_run",
                "provenance": {"kind": "structured_agent_result"},
                "status": bad_status,
            }]
        }
        claims = [{
            "claim_id": "bulk.de.1",
            "modality": "bulk",
            "analysis": "differential_expression",
            "tier": "associative",
            "evidence_card_id": "bulk.de.1#evidence",
        }]
        linkage = link_claims_to_ledger(claims, bad_ledger)
        assert linkage["n_violations"] == 1

        blk = NarrativeBlock(
            id="bulk.de.1",
            modality="bulk",
            analysis="differential_expression",
            block_type="result",
            title="DE",
            status="success",
            confidence="medium",
            claim="Found 5 DE genes.",
            evidence=[EvidenceItem("DE genes", 5, "deseq2")],
            metadata={
                "claim": {"tier": "associative"},
                "claim_verification": {
                    "status": "supported",
                    "evidence_card": {
                        "evidence_card_id": "bulk.de.1#evidence",
                        "node_type": "evidence_card",
                        "provenance": {"block_id": "bulk.de.1", "sources": ["deseq2"]},
                    },
                },
            },
        )
        with pytest.raises(LedgerLinkageError):
            verify_blocks_against_ledger([blk], bad_ledger, strict=True)

        sig = EvidenceSignal(
            entity="GATA1",
            entity_kind="gene",
            modality="bulk",
            measure="log2fc",
            audited_node_ref="ledger://bulk/differential_expression",
            value=2.0,
            direction="up",
        )
        hyp = Hypothesis(
            id="h_bad",
            mechanism="GATA1 may act",
            entities=["GATA1"],
            observation_refs=["ledger://bulk/differential_expression"],
            observed_claims=[{"signal_id": sig.signal_id, "stated_direction": "up"}],
            experiment=DiscriminatingExperiment("KD", "qPCR", "down", "unchanged"),
        )
        gr = verify_hypothesis_grounding(hyp, [sig], bad_ledger)
        assert gr.grounded is False
        assert gr.not_run_refs


# ── 7. CLD-02: Withheld claims close speculative verification rail #1 ─────────


def test_cld02_withheld_claims_close_speculative_verification_gate():
    rb = ReportBuilderMixin.__new__(ReportBuilderMixin)

    class _Block:
        metadata = {"claim_verification": {"status": "supported"}}

    clean_ledger = {"claim_ledger_verification": {"n_violations": 0}}
    r_ok = rb._speculative_verification_state(clean_ledger, [_Block()], n_withheld=0)
    assert r_ok.w_claim_passed is True and r_ok.gate_open is True

    r_withheld = rb._speculative_verification_state(clean_ledger, [_Block()], n_withheld=2)
    assert r_withheld.complete is True
    assert r_withheld.w_claim_passed is False
    assert r_withheld.gate_open is False
    assert "2 withheld" in r_withheld.w_claim_evidence


# ── 8. BI-12 / D5-H14 / D5-C4 / CLD-03 / D5-M1: DecisionManifest, replay, capsules, headless ──


def test_bi12_d5h14_decision_manifest_contrasts_none_thresholds_posthoc_and_portability(tmp_path):
    base_ctx = {
        "organism": "Mus musculus",
        "genome": "mm39",
        "reproducible_mode": True,
        "global_padj": None,
        "global_lfc": None,
        "data_dir": "/private/lab/study42",
        "input_files": [
            {"modality": "bulk_ATAC", "path": "/private/lab/study42/KO1.bam", "size_bytes": 100, "sha256": "a" * 64},
        ],
        "design": {
            "organism": "Mus musculus",
            "genome": "mm39",
            "main_factor": "condition",
            "contrasts": [
                {"test": "KO", "reference": "WT"},
                {"case": "Disease", "baseline": "Healthy"},
            ],
            "replicate_handling": {"source_path": "/private/lab/study42/meta.tsv"},
        },
    }
    intent = {
        "analysis_type": "differential",
        "raw_user_question": "Analyze /private/lab/study42/KO1.bam KO vs WT",
    }
    m1 = build_decision_manifest(
        experiment_id="exp1", exp_context=base_ctx, plan={}, biological_intent=intent
    )
    assert m1["thresholds"]["global_padj"] == 0.05
    assert m1["thresholds"]["global_lfc"] == 1.0
    assert m1["design"]["contrasts"] == [
        {"numerator": "KO", "denominator": "WT", "factor": "condition"},
        {"numerator": "Disease", "denominator": "Healthy", "factor": "condition"},
    ]
    assert "/private/lab/study42" not in json.dumps(m1)
    assert verify_decision_manifest(m1) == (True, "verified")

    # _portable_methodology MUST preserve cryptographic verification of decision_manifest!
    portable = _portable_methodology({"decision_manifest": m1})
    assert verify_decision_manifest(portable["decision_manifest"]) == (True, "verified")

    # Re-running a post_hoc_fork must preserve commitment_status == "post_hoc_fork"
    ledger_path = tmp_path / "fork.jsonl"
    _, fl1 = commit_decision_manifest(m1, ledger_path=ledger_path)
    assert fl1["fork_index"] == 0 and fl1["commitment_status"] == "locked_before_execution"

    fork_ctx = {**base_ctx, "global_lfc": 0.58}
    m2 = build_decision_manifest(experiment_id="exp2", exp_context=fork_ctx, plan={}, biological_intent=intent)
    _, fl2 = commit_decision_manifest(m2, ledger_path=ledger_path)
    assert fl2["fork_index"] == 1 and fl2["commitment_status"] == "post_hoc_fork"

    # Re-run exact forked spec with default commitment_status -> stays post_hoc_fork!
    m2_rerun = build_decision_manifest(experiment_id="exp3", exp_context=fork_ctx, plan={}, biological_intent=intent)
    _, fl2_rerun = commit_decision_manifest(m2_rerun, ledger_path=ledger_path)
    assert fl2_rerun["fork_index"] == 1
    assert fl2_rerun["is_exact_spec_replay"] is True
    assert fl2_rerun["commitment_status"] == "post_hoc_fork"


def test_d5c4_orchestrator_replay_preserves_specification_sha_and_fork_index(monkeypatch, tmp_path):
    pytest.importorskip("litellm")
    from aria.agents.orchestrator_agent import OrchestratorAgent
    from aria.memory.memory import ARIAMemory

    mem = ARIAMemory(db_path=str(tmp_path / "replay_mem.db"))
    orch = OrchestratorAgent(mem, llm=None)
    ledger_path = tmp_path / "replay_mem_fork_ledger.jsonl"

    base_ctx = {
        "organism": "Mus musculus",
        "genome": "mm39",
        "reproducible_mode": True,
        "modalities": {},
        "input_files": [{"modality": "bulk_RNA", "path": "/data/a.tsv", "size_bytes": 10, "sha256": "c" * 64}],
        "design": {"organism": "Mus musculus", "genome": "mm39", "main_factor": "cond", "groups": {"A": ["a1"], "B": ["b1"]}},
    }
    plan = {"steps": []}
    intent = {"analysis_type": "differential", "raw_user_question": "A vs B", "user_question": "A vs B"}
    m1 = build_decision_manifest(experiment_id="orig", exp_context=base_ctx, plan=plan, biological_intent=intent)
    m1_committed, fl1 = commit_decision_manifest(m1, ledger_path=ledger_path)
    assert fl1["fork_index"] == 0
    assert m1_committed["locked_parameter_decisions"] == {}

    from aria.utils.decision_manifest import seal_realized_parameter_decisions
    # Simulate post-run realized_parameter_decisions cryptographically sealed in methodology.json
    replayed_dm = seal_realized_parameter_decisions(
        m1_committed, {"leiden_clustering.resolution": 0.6}
    )

    captured: dict = {}

    def _stop_after_commit(exp_id, p, ctx):
        captured["ctx"] = ctx
        captured["intent"] = orch._experiment_plans[exp_id]["intent"]

    monkeypatch.setattr(orch, "_blocked_modalities", lambda m: {})
    monkeypatch.setattr(orch, "_experimental_modalities", lambda m: {})
    # Stub Setup/Narrative/downstream by making _dispatch_agents stop via empty modalities + mocked NarrativeAgent
    from aria.agents import narrative_agent as na_mod

    class _StubNA:
        name = "narrative_agent"

        def __init__(self, *a, **k):
            pass

        def run(self, exp_id, ctx):
            captured["exp_context"] = ctx["exp_context"]
            captured["biological_intent"] = ctx["biological_intent"]
            return {"status": "done", "report_path": str(tmp_path / "r.html")}

    monkeypatch.setattr(na_mod, "NarrativeAgent", _StubNA)
    orch._experiment_plans["replay_run"] = {"intent": dict(intent)}
    orch._dispatch_agents(
        "replay_run",
        plan,
        {**base_ctx, "replayed_decision_manifest": replayed_dm, "zero_llm_replay": True},
    )
    replay_dm_out = captured["exp_context"]["decision_manifest"]
    replay_fl = captured["exp_context"]["analytical_fork_ledger"]
    assert replay_dm_out["specification_sha256"] == m1_committed["specification_sha256"]
    assert replay_fl["fork_index"] == 0
    assert replay_fl["total_analytical_forks"] == 1
    assert replay_fl["is_exact_spec_replay"] is True
    # Runtime ParameterAdvisor still receives realized_parameter_decisions!
    assert captured["biological_intent"]["locked_parameter_decisions"] == {"leiden_clustering.resolution": 0.6}
    mem.close()


def test_cld03_d5m1_no_zero_byte_memory_snapshot_and_safe_headless(tmp_path):
    # 1. _write_memory_snapshot skips creating 0-byte file for :memory: and on export failure
    class _Host(ReportBuilderMixin):
        def __init__(self, memory):
            self.memory = memory

    rep_dir = tmp_path / "aria_reproducible_exp999"
    rep_dir.mkdir()
    _Host(type("M", (), {"db_path": ":memory:"})())._write_memory_snapshot(rep_dir, "exp999")
    assert not (rep_dir / "memory_snapshot.sqlite").exists()

    def _boom(*_a, **_k):
        raise RuntimeError("disk full")

    _Host(type("M", (), {"db_path": "/tmp/x.db", "export_experiment_snapshot": _boom})())._write_memory_snapshot(rep_dir, "exp999")
    assert not (rep_dir / "memory_snapshot.sqlite").exists()

    # 2. Capsule write & verify succeed without crashing when memory_snapshot.sqlite is absent
    (rep_dir / "report.html").write_text("<html></html>", encoding="utf-8")
    (rep_dir / "methodology.json").write_text(
        json.dumps({"provenance": {"version": "4.9.1", "git_commit": "abc"}, "inputs": [], "claims": [], "run_ledger": {"entries": []}}),
        encoding="utf-8",
    )
    (rep_dir / "ro-crate-metadata.json").write_text("{}", encoding="utf-8")
    cap = write_reproducible_capsule(rep_dir, repo_root=tmp_path / "empty_repo")
    ver = verify_reproducible_capsule(cap, repo_root=tmp_path / "empty_repo")
    assert ver["status"] == "pass"
    assert ver["memory_snapshot"]["present"] is False and ver["memory_snapshot"]["valid"] is True

    # 3. _find_report_on_disk does NOT hang when memory_snapshot.sqlite is absent in :memory: reproducible mode,
    #    yet still waits if a 0-byte in-progress file exists.
    assert _find_report_on_disk("exp999", reproducible_mode=True, reports_root=tmp_path) == str(rep_dir / "report.html")
    (rep_dir / "memory_snapshot.sqlite").write_bytes(b"")
    assert _find_report_on_disk("exp999", reproducible_mode=True, reports_root=tmp_path) is None

    # 4. default_answer_policy never raises IndexError on empty options
    for cp in (1, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2, 3.5, 4):
        assert isinstance(default_answer_policy(cp, "Question?", []), str)


# ── 9. BI-13 / D5-H1: Design power alpha, negative log2fc, two-sided, batch_map ──


def test_bi13_d5h1_design_power_alpha_two_sided_and_batch_prefix():
    ctx_005 = {
        "global_padj": 0.05,
        "design": {"groups": {"ctrl": [f"c{i}" for i in range(6)], "trt": [f"t{i}" for i in range(6)]}},
    }
    ctx_001 = {**ctx_005, "global_padj": 0.01}
    ctx_0001 = {**ctx_005, "global_padj": 0.001}
    p_005 = assess_design_power(ctx_005)["contrasts"][0]["power_estimate_at_target_log2fc"]
    p_001 = assess_design_power(ctx_001)["contrasts"][0]["power_estimate_at_target_log2fc"]
    p_0001 = assess_design_power(ctx_0001)["contrasts"][0]["power_estimate_at_target_log2fc"]
    assert p_005 > p_001 > p_0001

    # Negative target_log2fc is symmetric with positive target_log2fc
    p_pos = assess_design_power({**ctx_005, "target_log2fc": 1.2})["contrasts"][0]["power_estimate_at_target_log2fc"]
    p_neg = assess_design_power({**ctx_005, "target_log2fc": -1.2})["contrasts"][0]["power_estimate_at_target_log2fc"]
    assert p_pos == p_neg

    # Two-sided power at effect=0 equals alpha (~0.05), not alpha/2 (0.025)
    zero_res = _assess_contrast(
        "trt",
        "ctrl",
        {"ctrl": ["c1", "c2", "c3", "c4"], "trt": ["t1", "t2", "t3", "t4"]},
        {"alpha": 0.05, "target_log2fc": 0.0, "assumed_log2_sigma": 1.0},
    )
    assert zero_res["power_estimate_at_target_log2fc"] == pytest.approx(0.05, abs=0.002)

    # batch_map exact/longest-key matching: "S1" must NOT steal "S10"
    batch_res = assess_design_power({
        "design": {
            "groups": {"ctrl": ["S1", "S2", "S3"], "trt": ["S10", "S11", "S12"]},
            "batch_map": {"S1": "B1", "S2": "B1", "S3": "B1", "S10": "B2", "S11": "B2", "S12": "B2"},
        }
    })
    assert batch_res["checks"]["batch_condition_confounding"]["condition_batches"] == {
        "ctrl": ["B1"],
        "trt": ["B2"],
    }
    assert batch_res["checks"]["batch_condition_confounding"]["confounded"] is True


# ── 10. BI-14 / D5-H10: LOBO NaN-preserving BH and partially paired designs ──


def test_bi14_d5h10_lobo_nan_preserved_in_bh_and_partially_paired_fallback():
    # Adding zero-variance constant features (NaN raw p-values in Welch t-test) must NOT dilute BH padj!
    rng = np.random.default_rng(11)
    real_mat = np.array([
        [6.0, 6.2, 5.9, 1.0, 1.1, 0.9],
        [4.0, 4.3, 3.8, 2.0, 2.2, 1.9],
    ])
    const_zeros = np.zeros((10, 6))
    _, _, padj_real_only = _unpaired_fold_stats(real_mat, [0, 1, 2], [3, 4, 5], moderate=False)
    _, _, padj_with_const = _unpaired_fold_stats(
        np.vstack([real_mat, const_zeros]), [0, 1, 2], [3, 4, 5], moderate=False
    )
    np.testing.assert_allclose(padj_with_const[:2], padj_real_only)
    assert np.all(padj_with_const[2:] == 1.0)

    # Partially paired with n_pairs=1 (<3) but 4v4 arms -> falls back to unpaired LOBO!
    units = ["T1", "T2", "T3", "T4", "C1", "C2", "C3", "C4"]
    groups = {u: ("T" if u.startswith("T") else "C") for u in units}
    mat = rng.normal(2.0, 0.2, size=(5, 8))
    mat[0, :4] += 4.0
    res_fallback = compute_lobo_jackknife(
        [f"G{i}" for i in range(5)],
        mat,
        units,
        groups,
        numerator="T",
        denominator="C",
        unit_to_donor={"T1": "D1", "C1": "D1", "T2": "D2", "T3": "D3", "T4": "D4", "C2": "D5", "C3": "D6", "C4": "D7"},
    )
    assert res_fallback["status"] == "measured"
    assert res_fallback["design_mode"] == "unpaired_biological_replicate"
    assert res_fallback["n_jackknife_folds"] == 8

    # Partially paired with n_pairs=3 + 2 extra unpaired units -> paired_donor with n_unpaired_units_excluded == 2
    res_partial_paired = compute_lobo_jackknife(
        [f"G{i}" for i in range(5)],
        mat,
        units,
        groups,
        numerator="T",
        denominator="C",
        unit_to_donor={"T1": "D1", "C1": "D1", "T2": "D2", "C2": "D2", "T3": "D3", "C3": "D3", "T4": "D4", "C4": "D5"},
    )
    assert res_partial_paired["status"] == "measured"
    assert res_partial_paired["design_mode"] == "paired_donor"
    assert res_partial_paired["n_unpaired_units_excluded"] == 2
    assert res_partial_paired["n_jackknife_folds"] == 3


# ── 11. D5-H2: Sensitivity tokenizer letter-digit splits & clinical plurals ──


@pytest.mark.parametrize(
    "field",
    ["icd10", "icd9", "patient01", "patients", "subjects", "biopsies", "specimens", "patientid", "subjectid", "donorid"],
)
def test_d5h2_sensitivity_detects_clinical_plurals_and_alphanumeric_compounds(field):
    res = classify_sensitivity("human", [field])
    assert res["level"] == "high", (field, res)
    assert res["phi_signals"], (field, res)


@pytest.mark.parametrize("field", ["mRNA", "pseudobulk", "neutrophil"])
def test_d5h2_sensitivity_preserves_non_phi_scientific_tokens(field):
    res = classify_sensitivity("human", [field])
    assert res["level"] == "elevated"
    assert res["phi_signals"] == []


# ── 12. D5-C5, D5-C6, CDX-02: DataAudit json/boundaries, SetupAgent boundaries, DesignAgent CP2.2 ──


def test_d5c5_d5c6_cdx02_species_boundaries_manifest_json_and_mouse_cp22(tmp_path):
    # 1. Word-boundary matching in DataAuditAgent and SetupAgent prevents "rat" false positives
    da = DataAuditAgent.__new__(DataAuditAgent)
    g, org = da._infer_genome_organism(
        [Path("regulatory_ratio_integrative.h5ad")],
        tmp_path / "integrative_regulatory_ratio",
        "Analyze regulatory and integrative expression ratio across clusters",
    )
    assert (g, org) == ("unknown", "unknown")

    g_rat, org_rat = da._infer_genome_organism(
        [Path("rat_cortex_rn7.bam")], tmp_path, "Compare rat cortex"
    )
    assert (g_rat, org_rat) == ("rn7", "Rattus norvegicus")

    pytest.importorskip("litellm")
    from aria.agents.design_agent import DesignAgent
    from aria.agents.setup_agent import SetupAgent
    from aria.bus.message_bus import bus
    from aria.memory.memory import ARIAMemory

    sa = SetupAgent.__new__(SetupAgent)
    assert sa._detect_genome("", {"summary": "Integrative regulatory ratio analysis"}, audited_genome="unknown") is None
    assert sa._detect_genome("", {"summary": "Compare rat liver rn7"}, audited_genome="unknown") == "rn7"

    # 2. DataAuditAgent parses multiome_fastq_manifest.json without NameError (json imported)
    mm_dir = tmp_path / "mm_data"
    mm_dir.mkdir()
    (mm_dir / "multiome_fastq_manifest.json").write_text(
        json.dumps({"samples": [{"rna_r1": "s1_rna_R1.fastq.gz", "atac_r1": "s1_atac_R1.fastq.gz"}]}),
        encoding="utf-8",
    )
    (mm_dir / "s1_rna_R1.fastq.gz").write_bytes(b"")
    (mm_dir / "s1_atac_R1.fastq.gz").write_bytes(b"")
    mem = ARIAMemory(db_path=str(tmp_path / "audit_mem.db"))
    audit = DataAuditAgent(mem)
    res = audit.run("exp_mm_json", {"data_dir": str(mm_dir), "user_question": "Analyze integrative regulatory ratio"})
    assert res["status"] == "awaiting_checkpoint"
    assert isinstance(res["exp_context"]["multiome_fastq_manifest"], dict)
    assert res["exp_context"]["organism"] == "unknown"  # NOT falsely "Rattus norvegicus"!

    # 3. DesignAgent._publish_organism_checkpoint puts detected Mus musculus (mm39) as option #1
    design = DesignAgent(mem, llm=None)
    design._experiment_id = "exp_mouse_cp22"
    design._organism = "Mus musculus"
    design._genome = "mm39"
    design._publish_organism_checkpoint()
    pending = [
        m for m in bus.get_pending_checkpoints(experiment_id="exp_mouse_cp22")
        if m.payload.get("checkpoint") == 2.2
    ]
    assert pending
    opts = pending[-1].payload["options"]
    assert opts[0] == "Mus musculus (mm39)"
    assert default_answer_policy(2.2, pending[-1].payload["question"], opts) == "Mus musculus (mm39)"
    mem.close()
