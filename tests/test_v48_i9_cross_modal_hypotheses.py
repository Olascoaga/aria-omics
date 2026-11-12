"""V48 Slice I9 (S10, ADR-072): Cross-modal HypothesisAgent adapter & governance guards.

Covers:
  1. ``integration_evidence`` contract across ``peak_to_gene``, ``eregulons``,
     ``mofa`` (strict confounder exclusion + inline/disk parity + latent sign
     indeterminacy), and ``wnn`` anti-degradation guardrail propagation.
  2. Identifier & numerical hygiene (NaN/Inf/range rejection, peak-coordinate
     filtering, composite ``::`` motif rejection, sign-class consistency).
  3. Non-circular competitive ranking (ADR-070 + ADR-072): citing dependent
     paired cis representations (``scATAC`` ``peak2gene`` + ``integration``
     ``peak_to_gene`` + ``integration`` ``eregulons``) preserves every audited
     signal while counting as ONE independent lineage line without inflating
     ``mean_effect_norm``; genuinely distinct contrasts preserve independence.
  4. Full 4-gate + causal-gate + quarantine + opt-in-off enforcement on
     cross-modal hypotheses.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from aria.agents.hypothesis_agent import HypothesisAgent
from aria.agents.narrative.hypothesis import (
    CAVEAT_CODES,
    DiscriminatingExperiment,
    EvidenceSignal,
    Hypothesis,
    LLMProposer,
    SpeculativePromotionError,
    VerificationReceipt,
    assert_no_speculative_promotion,
    build_proposer_prompt,
    build_speculative_manifest,
    build_speculative_section,
    caveat_gloss,
    gather_evidence,
    integration_evidence,
    persist_speculative_manifest,
    rank_hypotheses,
    render_speculative_section_html,
)
from aria.agents.narrative.hypothesis.adapters.integration import (
    PAIRED_CIS_REGULATORY_KEY,
    PAIRED_MOFA_LATENT_KEY,
)
from aria.agents.narrative.hypothesis.caveats import (
    DEPTH_CONFOUNDED,
    MOFA_LATENT_ASSOCIATIVE,
    MOTIF_NOT_BINDING,
    PARALOG_MOTIF_AMBIGUITY,
    PEAK2GENE_ASSOCIATIVE,
    WNN_GUARDRAIL_DEGRADED,
)


def _exp(
    perturbation: str = "CRISPRi knockdown of ETV6",
    readout: str = "SOX6 transcript abundance by RT-qPCR",
    predicted_direction: str = "decrease",
    refuting_outcome: str = "SOX6 transcript abundance increases or stays unchanged",
) -> DiscriminatingExperiment:
    return DiscriminatingExperiment(
        perturbation=perturbation,
        readout=readout,
        predicted_direction=predicted_direction,
        refuting_outcome=refuting_outcome,
    )


def _integration_ledger(
    *,
    pqc: str = "ran",
    wnn: str = "ran",
    mofa: str = "ran",
    p2g: str = "ran",
    ereg: str = "ran",
    chrom_motif: str = "ran",
    chrom_reg: str = "ran",
    scrna_de: str = "ran",
) -> dict:
    return {
        "entries": [
            {"node_id": "ledger://integration/paired_qc", "status": pqc},
            {"node_id": "ledger://integration/wnn", "status": wnn},
            {"node_id": "ledger://integration/mofa", "status": mofa},
            {"node_id": "ledger://integration/peak_to_gene", "status": p2g},
            {"node_id": "ledger://integration/eregulons", "status": ereg},
            {"node_id": "ledger://chromatin/motif_enrichment", "status": chrom_motif},
            {"node_id": "ledger://chromatin/regulatory_layers", "status": chrom_reg},
            {"node_id": "ledger://scRNA/pseudobulk_de", "status": scrna_de},
        ]
    }


def _sample_integration_results(*, wnn_guardrail: str = "pass") -> dict:
    return {
        "integration_agent": {
            "status": "done",
            "findings": {
                "paired_qc": {
                    "status": "success",
                    "n_cells_before": 3143,
                    "n_cells_after": 2517,
                },
                "wnn": {
                    "status": "success",
                    "guardrail_status": wnn_guardrail,
                    "n_clusters": 7,
                    "clusters": ["0", "1", "2"],
                },
                "peak_to_gene": {
                    "status": "success",
                    "max_padj": 0.05,
                    "top_links": [
                        {
                            "gene": "SOX6",
                            "peak": "chr11:16001000-16001500",
                            "correlation": 0.61,
                            "padj": 1.2e-6,
                            "region_class": "distal_cis",
                        },
                        {
                            "gene": "GFAP",
                            "peak": "chr17:44910000-44910500",
                            "correlation": -0.44,
                            "padj": 3.4e-4,
                            "region_class": "promoter_proximal",
                        },
                        # Peak coordinate mistakenly placed in gene slot -> rejected
                        {
                            "gene": "chr1:10000-10500",
                            "peak": "chr1:10000-10500",
                            "correlation": 0.75,
                            "padj": 1e-5,
                        },
                        # Non-significant link -> rejected
                        {
                            "gene": "NOISEG",
                            "peak": "chr2:20000-20500",
                            "correlation": 0.25,
                            "padj": 0.22,
                        },
                    ],
                },
                "eregulons": {
                    "status": "success",
                    "max_padj": 0.05,
                    "top_eregulons": [
                        {
                            "regulon_id": "ETV6_pos_tf_pos_p2g",
                            "tf": "ETV6",
                            "sign_class": "pos_tf_pos_p2g",
                            "mean_r_tf_gene": 0.48,
                            "regulon_padj": 2.0e-5,
                            "target_genes": ["SOX6", "MBP", "ETV6", "PLP1", "MOG"],
                            "enhancer_peaks": ["chr11:16001000-16001500"],
                            "motif_family_ambiguity": False,
                            "depth_confounded": False,
                        },
                        {
                            "regulon_id": "NFIA_neg_tf_pos_p2g",
                            "tf": "NFIA",
                            "sign_class": "neg_tf_pos_p2g",
                            "mean_r_tf_gene": -0.39,
                            "regulon_padj": 4.5e-4,
                            "target_genes": ["GFAP", "AQP4"],
                            "motif_family_ambiguity": False,
                            "depth_confounded": False,
                        },
                        {
                            "regulon_id": "SPI1_pos_tf_pos_p2g",
                            "tf": "SPI1",
                            "sign_class": "pos_tf_pos_p2g",
                            "mean_r_tf_gene": 0.36,
                            "regulon_padj": 1.1e-3,
                            "target_genes": ["CSF1R", "CX3CR1"],
                            "motif_family_ambiguity": True,
                            "ambiguous_paralog_tfs": [{"tf": "IKZF1"}],
                            "depth_confounded": True,
                        },
                    ],
                },
                "mofa": {
                    "status": "success",
                    "views": {
                        "rna": {"kind": "rna"},
                        "atac": {"kind": "atac"},
                    },
                    "top_factors": [
                        {
                            "factor_id": 1,
                            "active": True,
                            "category": "shared",
                            "r2_by_view": {"rna": 0.14, "atac": 0.09},
                            "excluded_from_biological_interpretation": False,
                        },
                        # Technical depth factor — even if active/shared, MUST be excluded!
                        {
                            "factor_id": 5,
                            "active": True,
                            "category": "shared",
                            "r2_by_view": {"rna": 0.22, "atac": 0.31},
                            "excluded_from_biological_interpretation": False,
                        },
                        # Modality-specific factor -> omitted from cross-modal shared signals
                        {
                            "factor_id": 2,
                            "active": True,
                            "category": "modality_specific",
                            "r2_by_view": {"rna": 0.11, "atac": 0.002},
                            "excluded_from_biological_interpretation": False,
                        },
                    ],
                    "excluded_factors": [
                        {
                            "factor_id": 5,
                            "exclusion_reasons": ["technical_depth(|r|=0.941)"],
                        }
                    ],
                    "technical_factor_flag": True,
                    "technical_factor_id": 5,
                    "factor1_top_features": ["UNVERIFIED_BARE_FEATURE"],
                    "top_loadings": {
                        "factors": {
                            "Factor1": {
                                "rna": [
                                    {"feature": "MBP", "weight": 0.84},
                                    {"feature": "GAD1", "weight": -0.73},
                                    {"feature": "chr11:16001000-16001500", "weight": 0.99},
                                ],
                                "atac": [
                                    {"feature": "chr11:16001000-16001500", "weight": 0.95},
                                ],
                            },
                            "Factor5": {
                                "rna": [
                                    {"feature": "DEPTHGENE", "weight": 0.99},
                                ]
                            },
                            "Factor2": {
                                "rna": [
                                    {"feature": "RNAONLYGENE", "weight": 0.88},
                                ]
                            },
                        }
                    },
                },
            },
        }
    }


# ── 1. Adapter basics, ledger gates & fail-closed envelopes ─────────────────


def test_integration_adapter_empty_when_absent_or_blocked():
    assert integration_evidence({}) == []
    blocked = {
        "integration_agent": {
            "status": "blocked",
            "reason": "multiome_same_cell_gate_blocked",
            "findings": _sample_integration_results()["integration_agent"]["findings"],
        }
    }
    assert integration_evidence(blocked, _integration_ledger()) == []


def test_integration_adapter_emits_all_three_layers_with_exact_contracts():
    signals = integration_evidence(
        _sample_integration_results(), _integration_ledger()
    )
    assert signals

    # No peak coordinates, WNN numeric clusters, excluded Factor5, or modality-specific Factor2
    entities = {s.entity for s in signals}
    assert "chr1:10000-10500" not in entities
    assert "chr11:16001000-16001500" not in entities
    assert "0" not in entities
    assert "NOISEG" not in entities
    assert "DEPTHGENE" not in entities
    assert "RNAONLYGENE" not in entities
    assert "UNVERIFIED_BARE_FEATURE" not in entities

    # Layer 1: peak_to_gene
    p2g = {
        s.entity: s
        for s in signals
        if s.audited_node_ref == "ledger://integration/peak_to_gene"
    }
    assert set(p2g) == {"SOX6", "GFAP"}
    assert p2g["SOX6"].direction == "up"
    assert p2g["SOX6"].value == 0.61
    assert p2g["GFAP"].direction == "down"
    assert p2g["GFAP"].value == -0.44
    assert p2g["SOX6"].caveats_inherited == [PEAK2GENE_ASSOCIATIVE]
    assert p2g["SOX6"].independence_key == PAIRED_CIS_REGULATORY_KEY

    # Layer 2: eregulons (TFs + max 3 non-self target memberships with value=None)
    ereg_tfs = {
        s.entity: s
        for s in signals
        if s.audited_node_ref == "ledger://integration/eregulons"
        and s.entity_kind == "tf_motif"
    }
    assert set(ereg_tfs) == {"ETV6", "NFIA", "SPI1"}
    assert ereg_tfs["ETV6"].direction == "up"
    assert ereg_tfs["ETV6"].value == 0.48
    assert ereg_tfs["ETV6"].caveats_inherited == [
        MOTIF_NOT_BINDING,
        PEAK2GENE_ASSOCIATIVE,
    ]
    # Negative TF-gene correlation preserved faithfully as "down"
    assert ereg_tfs["NFIA"].direction == "down"
    assert ereg_tfs["NFIA"].value == -0.39
    # Conditional caveats attached ONLY to SPI1 where flagged
    assert set(ereg_tfs["SPI1"].caveats_inherited) == {
        MOTIF_NOT_BINDING,
        PEAK2GENE_ASSOCIATIVE,
        PARALOG_MOTIF_AMBIGUITY,
        DEPTH_CONFOUNDED,
    }

    etv6_targets = [
        s
        for s in signals
        if s.audited_node_ref == "ledger://integration/eregulons"
        and s.entity_kind == "gene"
        and "ETV6_pos_tf_pos_p2g" in s.context
    ]
    # Self-target "ETV6" skipped; capped at 3 targets; value=None & direction="na"
    assert [s.entity for s in etv6_targets] == ["SOX6", "MBP", "PLP1"]
    assert all(
        s.value is None
        and s.direction == "na"
        and s.measure == "eregulon_target_membership"
        and s.independence_key == PAIRED_CIS_REGULATORY_KEY
        for s in etv6_targets
    )

    # Layer 3: MOFA+ shared active non-excluded RNA loadings
    mofa_sigs = {
        s.entity: s
        for s in signals
        if s.audited_node_ref == "ledger://integration/mofa"
    }
    assert set(mofa_sigs) == {"MBP", "GAD1"}
    assert mofa_sigs["MBP"].direction == "na"
    assert mofa_sigs["MBP"].value == 0.84
    assert "positive_loading" in mofa_sigs["MBP"].context
    assert mofa_sigs["GAD1"].direction == "na"
    assert mofa_sigs["GAD1"].value == -0.73
    assert "negative_loading" in mofa_sigs["GAD1"].context
    assert mofa_sigs["MBP"].caveats_inherited == [MOFA_LATENT_ASSOCIATIVE]
    assert mofa_sigs["MBP"].independence_key == PAIRED_MOFA_LATENT_KEY


def test_per_stage_ledger_gating_is_independent():
    res = _sample_integration_results()
    sigs_no_ereg = integration_evidence(res, _integration_ledger(ereg="skipped"))
    assert not any(
        s.audited_node_ref == "ledger://integration/eregulons"
        for s in sigs_no_ereg
    )
    assert any(
        s.audited_node_ref == "ledger://integration/peak_to_gene"
        for s in sigs_no_ereg
    )
    assert any(
        s.audited_node_ref == "ledger://integration/mofa" for s in sigs_no_ereg
    )

    sigs_no_mofa = integration_evidence(res, _integration_ledger(mofa="error"))
    assert not any(
        s.audited_node_ref == "ledger://integration/mofa" for s in sigs_no_mofa
    )


def test_mofa_disk_json_equivalence_and_missing_file_safety(tmp_path: Path):
    res_inline = _sample_integration_results()
    inline_sigs = [
        s.to_dict()
        for s in integration_evidence(res_inline, _integration_ledger())
        if s.audited_node_ref == "ledger://integration/mofa"
    ]

    # Write identical loadings JSON to disk and remove inline dict
    res_disk = json.loads(json.dumps(res_inline))
    mofa_block = res_disk["integration_agent"]["findings"]["mofa"]
    loadings_payload = mofa_block.pop("top_loadings")
    loadings_file = tmp_path / "mofa_loadings.json"
    loadings_file.write_text(json.dumps(loadings_payload), encoding="utf-8")
    mofa_block["factor_loadings_path"] = str(loadings_file)

    disk_sigs = [
        s.to_dict()
        for s in integration_evidence(res_disk, _integration_ledger())
        if s.audited_node_ref == "ledger://integration/mofa"
    ]
    assert disk_sigs == inline_sigs

    # Missing file on disk -> omits MOFA+ cleanly, NEVER falls back to bare factor1_top_features
    mofa_block["factor_loadings_path"] = str(tmp_path / "does_not_exist.json")
    missing_sigs = [
        s
        for s in integration_evidence(res_disk, _integration_ledger())
        if s.audited_node_ref == "ledger://integration/mofa"
    ]
    assert missing_sigs == []


def test_numerical_hygiene_and_composite_tf_rejection():
    dirty = {
        "integration_agent": {
            "status": "done",
            "findings": {
                "peak_to_gene": {
                    "status": "success",
                    "top_links": [
                        {"gene": "NAN1", "peak": "chr1:1-2", "correlation": float("nan")},
                        {"gene": "INF1", "peak": "chr1:1-2", "correlation": float("inf")},
                        {"gene": "OOB1", "peak": "chr1:1-2", "correlation": 1.35},
                        {"gene": "ZERO1", "peak": "chr1:1-2", "correlation": 0.0},
                        {"gene": "chr3-100-500", "peak": "chr3:100-500", "correlation": 0.5},
                    ],
                },
                "eregulons": {
                    "status": "success",
                    "top_eregulons": [
                        # Composite :: heterodimer motif string -> rejected, not split!
                        {
                            "regulon_id": "JUN::FOS_pos",
                            "tf": "JUN::FOS",
                            "sign_class": "pos_tf_pos_p2g",
                            "mean_r_tf_gene": 0.52,
                            "target_genes": ["FOSL1"],
                        },
                        # Sign mismatch: pos_tf_pos_p2g with negative r -> rejected!
                        {
                            "regulon_id": "BADSIGN_pos",
                            "tf": "BADSIGN",
                            "sign_class": "pos_tf_pos_p2g",
                            "mean_r_tf_gene": -0.41,
                            "target_genes": ["SOX6"],
                        },
                    ],
                },
            },
        }
    }
    assert integration_evidence(dirty, _integration_ledger()) == []


def test_wnn_degraded_guardrail_propagates_caveat():
    signals = integration_evidence(
        _sample_integration_results(wnn_guardrail="degraded"),
        _integration_ledger(),
    )
    assert signals
    assert all(
        WNN_GUARDRAIL_DEGRADED in s.caveats_inherited for s in signals
    )


# ── 2. Decisive Non-Circular Ranking Proof (ADR-070 + ADR-072) ──────────────


def test_non_circular_ranking_collapses_shared_cis_lineage_but_keeps_contrasts():
    combined_results = _sample_integration_results()
    # Add co-run ChromatinAgent (scATAC peak2gene on same paired container +
    # cluster motif enrichment) AND independent scRNA pseudobulk DE contrast.
    combined_results["chromatin_agent"] = {
        "status": "done",
        "findings": {
            "motifs": {
                "ran": True,
                "padj_max": 0.05,
                "per_group": {
                    "Oligo": {
                        "top_motifs": [
                            {
                                "name": "ETV6",
                                "log2_enrichment": 2.1,
                                "padj": 1e-4,
                            }
                        ]
                    }
                },
            },
            "regulatory": {
                "ran": True,
                "peak_to_gene": {
                    "ran": True,
                    "top_links": [
                        {
                            "gene": "SOX6",
                            "correlation": 0.61,
                            "direction": "positive",
                        }
                    ],
                },
            },
        },
    }
    combined_results["scrna_agent"] = {
        "status": "done",
        "findings": {
            "pseudobulk_de": {
                "per_group": {
                    "Oligo": {
                        "per_comparison": {
                            "aged_vs_young": {
                                "status": "success",
                                "top_genes": [
                                    {"gene": "SOX6", "log2fc": 1.25},
                                ],
                            }
                        }
                    }
                }
            }
        },
    }

    ledger = _integration_ledger()
    signals = gather_evidence(combined_results, ledger, {})

    # 1. No destructive dedup: all 3 cis representations exist with distinct signal_ids
    scatac_p2g = next(
        s
        for s in signals
        if s.modality == "scATAC" and s.measure == "peak2gene" and s.entity == "SOX6"
    )
    integ_p2g = next(
        s
        for s in signals
        if s.modality == "integration"
        and s.measure == "peak2gene_correlation"
        and s.entity == "SOX6"
    )
    integ_ereg_tf = next(
        s
        for s in signals
        if s.modality == "integration"
        and s.measure == "eregulon_mean_tf_gene_correlation"
        and s.entity == "ETV6"
    )
    scrna_de = next(
        s
        for s in signals
        if s.modality == "scRNA" and s.measure == "log2fc" and s.entity == "SOX6"
    )
    atac_motif = next(
        s
        for s in signals
        if s.modality == "scATAC"
        and s.audited_node_ref == "ledger://chromatin/motif_enrichment"
        and s.entity == "ETV6"
    )

    assert (
        len({scatac_p2g.signal_id, integ_p2g.signal_id, integ_ereg_tf.signal_id})
        == 3
    )
    # Single-modality DE signal keeps empty independence_key & clean to_dict()
    assert scrna_de.independence_key == ""
    assert "independence_key" not in scrna_de.to_dict()

    # 2. Hypothesis citing ONLY the 3 dependent cis signals -> 1 line, max effect 0.61!
    h_dependent_three = Hypothesis(
        id="dep3",
        mechanism="ETV6 may co-vary with SOX6 cis accessibility",
        entities=["ETV6", "SOX6"],
        observation_refs=[
            "ledger://chromatin/regulatory_layers",
            "ledger://integration/peak_to_gene",
            "ledger://integration/eregulons",
        ],
        observed_claims=[
            {"signal_id": scatac_p2g.signal_id, "stated_direction": "up"},
            {"signal_id": integ_p2g.signal_id, "stated_direction": "up"},
            {"signal_id": integ_ereg_tf.signal_id, "stated_direction": "up"},
        ],
        experiment=_exp(),
    )
    # 3. Hypothesis citing 1 cis signal + independent DE contrast -> check effect invariance!
    h_one_cis_plus_de = Hypothesis(
        id="cis1_de",
        mechanism="ETV6 and SOX6 may vary across aging",
        entities=["ETV6", "SOX6"],
        observation_refs=[
            "ledger://integration/peak_to_gene",
            "ledger://scRNA/pseudobulk_de",
        ],
        observed_claims=[
            {"signal_id": integ_p2g.signal_id, "stated_direction": "up"},
            {"signal_id": scrna_de.signal_id, "stated_direction": "up"},
        ],
        experiment=_exp(),
    )
    h_padded_cis_plus_de = Hypothesis(
        id="cis3_de",
        mechanism="ETV6 and SOX6 may vary across aging via enhancer",
        entities=["ETV6", "SOX6"],
        observation_refs=[
            "ledger://chromatin/regulatory_layers",
            "ledger://integration/peak_to_gene",
            "ledger://integration/eregulons",
            "ledger://scRNA/pseudobulk_de",
        ],
        observed_claims=[
            {"signal_id": scatac_p2g.signal_id, "stated_direction": "up"},
            {"signal_id": integ_p2g.signal_id, "stated_direction": "up"},
            {"signal_id": integ_ereg_tf.signal_id, "stated_direction": "up"},
            {"signal_id": scrna_de.signal_id, "stated_direction": "up"},
        ],
        experiment=_exp(),
    )
    # 4. Hypothesis with 3 genuinely distinct lines (paired cis + scRNA DE + cluster motif)
    h_three_genuine = Hypothesis(
        id="true3",
        mechanism="ETV6 motif enrichment and aging DE may converge on SOX6",
        entities=["ETV6", "SOX6"],
        observation_refs=[
            "ledger://integration/eregulons",
            "ledger://scRNA/pseudobulk_de",
            "ledger://chromatin/motif_enrichment",
        ],
        observed_claims=[
            {"signal_id": integ_ereg_tf.signal_id, "stated_direction": "up"},
            {"signal_id": scrna_de.signal_id, "stated_direction": "up"},
            {"signal_id": atac_motif.signal_id, "stated_direction": "up"},
        ],
        experiment=_exp(),
    )

    ranked = rank_hypotheses(
        [h_dependent_three, h_one_cis_plus_de, h_padded_cis_plus_de, h_three_genuine],
        signals,
    )
    assert [h.id for h in ranked][0] == "true3"
    assert h_dependent_three.rank_evidence["n_independent_lines"] == 1
    assert h_dependent_three.rank_evidence["mean_effect_norm"] == 0.61
    # Padding 2 extra dependent cis representations changes NEITHER line count NOR mean_effect_norm!
    assert h_one_cis_plus_de.rank_evidence["n_independent_lines"] == 2
    assert h_padded_cis_plus_de.rank_evidence["n_independent_lines"] == 2
    expected_two_line_eff = round((0.61 + math.tanh(1.25)) / 2.0, 4)
    assert h_one_cis_plus_de.rank_evidence["mean_effect_norm"] == expected_two_line_eff
    assert (
        h_padded_cis_plus_de.rank_evidence["mean_effect_norm"]
        == expected_two_line_eff
    )
    assert h_three_genuine.rank_evidence["n_independent_lines"] == 3


# ── 3. Causal gate, 4 publication gates & quarantine wall ───────────────────


def test_new_caveat_codes_are_registered_with_glosses_and_prompt_tags():
    for code in (
        PARALOG_MOTIF_AMBIGUITY,
        DEPTH_CONFOUNDED,
        MOFA_LATENT_ASSOCIATIVE,
        WNN_GUARDRAIL_DEGRADED,
    ):
        assert code in CAVEAT_CODES
        assert caveat_gloss(code) != code

    sigs = integration_evidence(
        _sample_integration_results(wnn_guardrail="warn"),
        _integration_ledger(),
    )
    prompt = build_proposer_prompt(
        sigs, {"biological_question": "Candidate oligodendrocyte regulators?"}, 3
    )
    assert "[paralog_motif_ambiguity]" in prompt
    assert "[depth_confounded]" in prompt
    assert "[mofa_latent_associative]" in prompt
    assert "[wnn_guardrail_degraded]" in prompt


def test_closed_verification_receipt_never_calls_proposer():
    sigs = integration_evidence(
        _sample_integration_results(), _integration_ledger()
    )
    called = {"n": 0}

    def _spy_proposer(s, c):
        called["n"] += 1
        return []

    closed_receipt = VerificationReceipt(
        w_claim_passed=True,
        w_ledger_passed=False,
        complete=True,
        w_claim_evidence="ok",
        w_ledger_evidence="n_violations=1",
    )
    out = HypothesisAgent(proposer=_spy_proposer).generate(
        sigs, _integration_ledger(), verification=closed_receipt
    )
    assert out["ran"] is False
    assert called["n"] == 0


def test_four_gates_and_quarantine_on_cross_modal_hypotheses(tmp_path: Path):
    res = _sample_integration_results(wnn_guardrail="warn")
    ledger = _integration_ledger()
    sigs = gather_evidence(res, ledger, {})

    by_key = {(s.entity, s.measure): s for s in sigs}
    etv6_tf = by_key[("ETV6", "eregulon_mean_tf_gene_correlation")]
    sox6_p2g = by_key[("SOX6", "peak2gene_correlation")]
    nfia_tf = by_key[("NFIA", "eregulon_mean_tf_gene_correlation")]
    spi1_tf = by_key[("SPI1", "eregulon_mean_tf_gene_correlation")]
    mbp_mofa = by_key[("MBP", "mofa_shared_rna_loading")]

    candidates_payload = [
        # 1. VALID positive cross-modal hypothesis (ETV6 + SOX6) owning all visible caveats
        {
            "id": "valid_etv6",
            "mechanism": (
                "ETV6 expression and cis enhancer accessibility may co-vary with "
                "SOX6 transcription in paired cells"
            ),
            "entities": ["ETV6", "SOX6"],
            "observation_refs": [
                "ledger://integration/eregulons",
                "ledger://integration/peak_to_gene",
            ],
            "observed_claims": [
                {"signal_id": etv6_tf.signal_id, "stated_direction": "up"},
                {"signal_id": sox6_p2g.signal_id, "stated_direction": "up"},
            ],
            "experiment": _exp().to_dict(),
            "devils_advocate": {
                "simpler_explanation": (
                    "Co-expression within a shared lineage state could explain "
                    "the correlation without ETV6 acting on SOX6"
                ),
                "confounds": [
                    "motif_not_binding",
                    "peak2gene_associative",
                    "wnn_guardrail_degraded",
                ],
            },
        },
        # 2. VALID negative eRegulon + MOFA+ hypothesis owning MOFA + WNN caveats
        {
            "id": "valid_nfia_mbp",
            "mechanism": (
                "NFIA inverse target correlation alongside MBP shared latent "
                "loading might reflect opposing glial state programs"
            ),
            "entities": ["NFIA", "MBP"],
            "observation_refs": [
                "ledger://integration/eregulons",
                "ledger://integration/mofa",
            ],
            "observed_claims": [
                {"signal_id": nfia_tf.signal_id, "stated_direction": "down"},
                {"signal_id": mbp_mofa.signal_id, "stated_direction": "na"},
            ],
            "experiment": {
                "perturbation": "NFIA overexpression in primary glia",
                "readout": "MBP mRNA by RT-qPCR",
                "predicted_direction": "decrease",
                "refuting_outcome": "MBP mRNA increases or remains unchanged",
            },
            "devils_advocate": {
                "simpler_explanation": "Independent lineage maturation gradient",
                "confounds": [
                    "motif_not_binding",
                    "peak2gene_associative",
                    "mofa_latent_associative",
                    "wnn guardrail warning",
                ],
            },
        },
        # 3. REJECTED by grounding (H15): NFIA is neg_tf_pos_p2g ("down") misread as "up"
        {
            "id": "bad_direction_nfia",
            "mechanism": "NFIA may elevate GFAP transcription",
            "entities": ["NFIA", "GFAP"],
            "observation_refs": ["ledger://integration/eregulons"],
            "observed_claims": [
                {"signal_id": nfia_tf.signal_id, "stated_direction": "up"}
            ],
            "experiment": _exp().to_dict(),
            "devils_advocate": {
                "simpler_explanation": "Shared astrocyte state",
                "confounds": [
                    "motif_not_binding",
                    "peak2gene_associative",
                    "wnn_guardrail_degraded",
                ],
            },
        },
        # 4. REJECTED by devils_advocate: cites SPI1 but hides paralog_motif_ambiguity & depth_confounded
        {
            "id": "bad_hides_paralog",
            "mechanism": "SPI1 may modulate CSF1R in microglia",
            "entities": ["SPI1", "CSF1R"],
            "observation_refs": ["ledger://integration/eregulons"],
            "observed_claims": [
                {"signal_id": spi1_tf.signal_id, "stated_direction": "up"}
            ],
            "experiment": {
                "perturbation": "SPI1 CRISPRi knockdown",
                "readout": "CSF1R transcript level by RT-qPCR",
                "predicted_direction": "decrease",
                "refuting_outcome": "CSF1R unchanged",
            },
            "devils_advocate": {
                "simpler_explanation": "Microglial abundance shift",
                "confounds": [
                    "motif_not_binding",
                    "peak2gene_associative",
                    "wnn_guardrail_degraded",
                ],
            },
        },
        # 5. REJECTED by language gate: causal verb "drives"
        {
            "id": "bad_causal_verb",
            "mechanism": "ETV6 drives SOX6 expression via enhancer binding",
            "entities": ["ETV6", "SOX6"],
            "observation_refs": ["ledger://integration/eregulons"],
            "observed_claims": [
                {"signal_id": etv6_tf.signal_id, "stated_direction": "up"}
            ],
            "experiment": _exp().to_dict(),
            "devils_advocate": {
                "simpler_explanation": "Co-expression",
                "confounds": [
                    "motif_not_binding",
                    "peak2gene_associative",
                    "wnn_guardrail_degraded",
                ],
            },
        },
        # 6. REJECTED by grounding: unmeasured TF FOXP3 in prose/entities
        {
            "id": "bad_ungrounded_tf",
            "mechanism": "FOXP3 may cooperate with ETV6 at SOX6 enhancers",
            "entities": ["ETV6", "SOX6"],
            "observation_refs": ["ledger://integration/eregulons"],
            "observed_claims": [
                {"signal_id": etv6_tf.signal_id, "stated_direction": "up"}
            ],
            "experiment": _exp().to_dict(),
            "devils_advocate": {
                "simpler_explanation": "Co-expression",
                "confounds": [
                    "motif_not_binding",
                    "peak2gene_associative",
                    "wnn_guardrail_degraded",
                ],
            },
        },
    ]

    proposer = LLMProposer(
        lambda prompt, system: json.dumps(candidates_payload),
        model_label="deterministic-test-proposer",
    )
    section = build_speculative_section(
        res,
        ledger,
        {"reproducible_mode": True},
        proposer=proposer,
        w_claim_passed=True,
        w_ledger_passed=True,
    )
    assert section is not None
    assert section["ran"] is True
    assert section["honest_null"] is False

    accepted_ids = {h["id"] for h in section["hypotheses"]}
    assert accepted_ids == {"valid_etv6", "valid_nfia_mbp"}

    rejected_by_id = {
        r["hypothesis_id"]: {f["gate"]: f for f in r["failures"]}
        for r in section["rejected"]
    }
    assert set(rejected_by_id) == {
        "bad_direction_nfia",
        "bad_hides_paralog",
        "bad_causal_verb",
        "bad_ungrounded_tf",
    }
    assert "grounding" in rejected_by_id["bad_direction_nfia"]
    assert rejected_by_id["bad_direction_nfia"]["grounding"]["contradicting_claims"]
    assert "devils_advocate" in rejected_by_id["bad_hides_paralog"]
    assert "paralog_motif_ambiguity" in (
        rejected_by_id["bad_hides_paralog"]["devils_advocate"]["reason"] or ""
    )
    assert "depth_confounded" in (
        rejected_by_id["bad_hides_paralog"]["devils_advocate"]["reason"] or ""
    )
    assert "language" in rejected_by_id["bad_causal_verb"]
    assert "grounding" in rejected_by_id["bad_ungrounded_tf"]
    assert "FOXP3" in (
        rejected_by_id["bad_ungrounded_tf"]["grounding"]["ungrounded_prose_entities"]
    )

    # Quarantine & non-promotion wall
    manifest_path = persist_speculative_manifest(
        section, tmp_path, reproducible=True
    )
    assert manifest_path is not None and manifest_path.is_file()
    manifest = build_speculative_manifest(
        section, generated_utc="<timestamp redacted for byte-identity>"
    )
    assert manifest["promotable"] is False
    assert manifest["tier"] == "SPECULATIVE"
    for h in section["hypotheses"]:
        assert h["tier"] == "SPECULATIVE"
        assert str(h["ledger_node"]).startswith("hypothesis://")
        with pytest.raises(SpeculativePromotionError):
            assert_no_speculative_promotion(
                [{"claim_id": "leak", "ledger_node": h["ledger_node"]}]
            )

    html = render_speculative_section_html(section)
    assert "SPECULATIVE" in html
    assert "not part of the audited claim manifest" in html
    assert "MOFA+ factor loading reflects undirected latent co-variation" in html
    assert "anti-degradation guardrail warnings" in html


def test_paralog_and_depth_caveats_pass_when_explicitly_acknowledged():
    res = _sample_integration_results()
    ledger = _integration_ledger()
    sigs = gather_evidence(res, ledger, {})
    spi1_tf = next(
        s
        for s in sigs
        if s.entity == "SPI1"
        and s.measure == "eregulon_mean_tf_gene_correlation"
    )
    csf1r_tg = next(
        s
        for s in sigs
        if s.entity == "CSF1R" and s.measure == "eregulon_target_membership"
    )

    payload = [
        {
            "id": "spi1_honest",
            "mechanism": (
                "SPI1 motif-bearing cis enhancers may co-vary with CSF1R "
                "transcript levels in microglial cells"
            ),
            "entities": ["SPI1", "CSF1R"],
            "observation_refs": ["ledger://integration/eregulons"],
            "observed_claims": [
                {"signal_id": spi1_tf.signal_id, "stated_direction": "up"},
                {"signal_id": csf1r_tg.signal_id, "stated_direction": "na"},
            ],
            "experiment": {
                "perturbation": "SPI1 CRISPRi knockdown in human microglia",
                "readout": "CSF1R mRNA by RT-qPCR",
                "predicted_direction": "decrease",
                "refuting_outcome": "CSF1R mRNA unchanged or increased",
            },
            "devils_advocate": {
                "simpler_explanation": (
                    "Co-expressed paralog transcription factors or library-size "
                    "variation could explain the shared accessibility-expression profile"
                ),
                "confounds": [
                    "motif_not_binding",
                    "peak2gene_associative",
                    "paralog motif family sharing",
                    "sequencing depth covariation",
                ],
            },
        },
        # Falsifiability failure: vague direction "changes"
        {
            "id": "spi1_vague_exp",
            "mechanism": "SPI1 may modulate CSF1R expression",
            "entities": ["SPI1", "CSF1R"],
            "observation_refs": ["ledger://integration/eregulons"],
            "observed_claims": [
                {"signal_id": spi1_tf.signal_id, "stated_direction": "up"}
            ],
            "experiment": {
                "perturbation": "SPI1 knockdown",
                "readout": "CSF1R mRNA",
                "predicted_direction": "changes significantly",
                "refuting_outcome": "no effect",
            },
            "devils_advocate": {
                "simpler_explanation": "Paralog TF activity",
                "confounds": [
                    "motif_not_binding",
                    "peak2gene_associative",
                    "paralog_motif_ambiguity",
                    "depth_confounded",
                ],
            },
        },
    ]
    proposer = LLMProposer(lambda p, s: json.dumps(payload))
    section = build_speculative_section(
        res,
        ledger,
        {},
        proposer=proposer,
        w_claim_passed=True,
        w_ledger_passed=True,
    )
    assert [h["id"] for h in section["hypotheses"]] == ["spi1_honest"]
    assert set(section["hypotheses"][0]["inherited_caveats"]) == {
        "motif_not_binding",
        "peak2gene_associative",
        "paralog_motif_ambiguity",
        "depth_confounded",
    }
    rej_gates = {f["gate"] for f in section["rejected"][0]["failures"]}
    assert "falsifiability" in rej_gates


def test_multiple_shared_mofa_factors_count_as_one_latent_line_and_cell_cycle_excluded():
    res = _sample_integration_results()
    mofa = res["integration_agent"]["findings"]["mofa"]
    mofa["top_factors"].extend(
        [
            {
                "factor_id": 3,
                "active": True,
                "category": "shared",
                "r2_by_view": {"rna": 0.08, "atac": 0.06},
                "excluded_from_biological_interpretation": False,
            },
            {
                "factor_id": 4,
                "active": True,
                "category": "shared",
                "r2_by_view": {"rna": 0.10, "atac": 0.09},
                "excluded_from_biological_interpretation": False,
            },
        ]
    )
    mofa["cell_cycle_factor_ids"] = [4]
    mofa["top_loadings"]["factors"]["Factor3"] = {
        "rna": [{"feature": "OLIG2", "weight": 0.67}]
    }
    mofa["top_loadings"]["factors"]["Factor4"] = {
        "rna": [{"feature": "MKI67", "weight": 0.91}]
    }

    sigs = integration_evidence(res, _integration_ledger())
    mofa_by_ent = {
        s.entity: s
        for s in sigs
        if s.audited_node_ref == "ledger://integration/mofa"
    }
    assert "MKI67" not in mofa_by_ent  # cell-cycle Factor4 excluded!
    assert {"MBP", "GAD1", "OLIG2"} <= set(mofa_by_ent)

    h_two_factors = Hypothesis(
        id="mofa2",
        mechanism="MBP and OLIG2 may co-vary along shared glial latent axes",
        entities=["MBP", "OLIG2"],
        observation_refs=["ledger://integration/mofa"],
        observed_claims=[
            {"signal_id": mofa_by_ent["MBP"].signal_id, "stated_direction": "na"},
            {"signal_id": mofa_by_ent["OLIG2"].signal_id, "stated_direction": "na"},
        ],
        experiment=_exp(),
    )
    rank_hypotheses([h_two_factors], sigs)
    assert h_two_factors.rank_evidence["n_independent_lines"] == 1
    assert h_two_factors.rank_evidence["mean_effect_norm"] == 0.0


def test_narrative_render_opt_in_default_off_vs_on_for_integration(tmp_path: Path):
    import sys
    import types

    litellm_stub = types.ModuleType("litellm")
    litellm_stub.completion = lambda *args, **kwargs: None
    sys.modules.setdefault("litellm", litellm_stub)

    from aria.agents.narrative.reporting.render import ReportRenderMixin

    class _DummyRenderer(ReportRenderMixin):
        llm = None

    renderer = _DummyRenderer()
    res = _sample_integration_results()
    ledger = _integration_ledger()
    ledger["claim_ledger_verification"] = {"n_violations": 0}

    # Default OFF: returns "" and writes NO speculative_hypotheses.json
    off_dir = tmp_path / "off"
    off_dir.mkdir()
    html_off = renderer._build_speculative_section_html(
        res, ledger, {}, narrative_blocks=[], report_dir=off_dir
    )
    assert html_off == ""
    assert not (off_dir / "speculative_hypotheses.json").exists()

    # Opt-in ON (with default honest-null proposer when llm=None): writes manifest
    on_dir = tmp_path / "on"
    on_dir.mkdir()
    html_on = renderer._build_speculative_section_html(
        res,
        ledger,
        {"enable_hypotheses": True, "reproducible_mode": True},
        narrative_blocks=[],
        report_dir=on_dir,
    )
    assert "SPECULATIVE" in html_on
    assert (on_dir / "speculative_hypotheses.json").is_file()

