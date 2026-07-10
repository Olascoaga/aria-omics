"""scATAC P4.3 governance tie-in: TOBIAS differential TF binding becomes a GOVERNED
associative narrative block with footprint<->RNA cross-evidence. The claim is capped
associative (no causal language, counts-only so W-CLAIM is satisfied), and the
cross-modal concordance is honest (unavailable RNA is not assumed concordant).
"""

from __future__ import annotations

import pytest


def _summary():
    return {
        "parsed": True, "n_motifs_tested": 879, "n_ranked_candidates": 873,
        "ranking_basis": {"fdr_controlled": False, "replicate_inference": False,
                          "pseudobulk": True},
        "top_toward_Monocyte": [
            {"tf": "CEBPB", "change": 1.01, "pvalue": 1e-195, "n_sites": 8987},
            {"tf": "FOSL2::JUND", "change": 0.48, "pvalue": 1e-188, "n_sites": 25363},
        ],
        "top_toward_T_cell": [
            {"tf": "EGR1", "change": -0.27, "pvalue": 1e-180, "n_sites": 56877},
            {"tf": "GHOST_TF", "change": -0.10, "pvalue": 1e-9, "n_sites": 30},
        ],
    }


def test_footprint_rna_concordance_scores_cross_modal():
    from aria.agents.narrative.synthesis.footprint_rna import footprint_rna_concordance

    rna = {
        "CEBPB": {"Monocyte": 3.0, "T_cell": 0.5},   # concordant (footprint+RNA Mono)
        "JUND":  {"Monocyte": 2.0, "T_cell": 1.0},   # dimer component, concordant Mono
        "EGR1":  {"Monocyte": 0.2, "T_cell": 1.5},   # concordant (footprint+RNA T_cell)
        # GHOST_TF has no RNA -> not checkable
    }
    res = footprint_rna_concordance(_summary(), rna, "Monocyte", "T_cell")
    assert res["n_checked"] == 3          # CEBPB, FOSL2::JUND(JUND), EGR1
    assert res["n_concordant"] == 3
    assert res["concordance_rate"] == pytest.approx(1.0)
    ghost = res["toward_T_cell"][1]
    assert ghost["rna_checked"] is False and ghost["rna_concordant"] is None


def test_footprint_rna_concordance_marks_discordant():
    from aria.agents.narrative.synthesis.footprint_rna import footprint_rna_concordance

    rna = {"CEBPB": {"Monocyte": 0.1, "T_cell": 2.0}}  # footprint Mono but RNA T_cell
    summary = {"parsed": True, "n_significant": 1, "n_motifs_tested": 1,
               "top_toward_Monocyte": [{"tf": "CEBPB", "change": 1.0,
                                        "pvalue": 1e-9, "n_sites": 50}],
               "top_toward_T_cell": []}
    res = footprint_rna_concordance(summary, rna, "Monocyte", "T_cell")
    assert res["n_checked"] == 1 and res["n_concordant"] == 0
    assert res["toward_Monocyte"][0]["rna_concordant"] is False


def test_footprint_block_is_governed_associative():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.validators import find_causal_language
    from aria.agents.narrative.synthesis.footprint_rna import footprint_rna_concordance

    rna = {"CEBPB": {"Monocyte": 3.0, "T_cell": 0.5},
           "JUND": {"Monocyte": 2.0, "T_cell": 1.0},
           "EGR1": {"Monocyte": 0.2, "T_cell": 1.5}}
    cross = footprint_rna_concordance(_summary(), rna, "Monocyte", "T_cell")
    footprinting = {"ran": True, "group_a": "Monocyte", "group_b": "T_cell",
                    "differential_summary": _summary(), "rna_cross_evidence": cross}

    block = ChromatinNarrator()._footprint_block(footprinting)
    assert block.status == "success" and block.confidence == "medium"
    assert "873" in block.claim and "879" in block.claim
    # cross-modal counts surfaced honestly
    assert "3" in block.claim  # n_concordant
    # associative only: no causal language in the claim; mandatory caveat present
    assert not find_causal_language(block.claim)
    assert any("associative" in c.text.lower() for c in block.caveats)
    # no specific TF name in the claim (W-CLAIM named-entity safety)
    assert "CEBPB" not in block.claim and "EGR1" not in block.claim


def test_footprint_block_surfaces_aggregate_figures():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    footprinting = {"ran": True, "group_a": "Monocyte", "group_b": "T_cell",
                    "differential_summary": _summary(),
                    "aggregate_plots": {
                        "CEBPB": {"png": "/tmp/x/CEBPB_aggregate.png",
                                  "svg": "/tmp/x/CEBPB_aggregate.svg"},
                        "EGR1": {"png": "/tmp/x/EGR1_aggregate.png"}}}
    block = ChromatinNarrator()._footprint_block(footprinting)
    paths = {f["path"] for f in block.figures}
    assert "/tmp/x/CEBPB_aggregate.png" in paths
    assert len(block.figures) == 2
    assert all("footprint" in f["caption"].lower() for f in block.figures)


def test_collect_emits_footprint_block_when_present():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    footprinting = {"ran": True, "group_a": "Monocyte", "group_b": "T_cell",
                    "differential_summary": _summary()}
    blocks = ChromatinNarrator().collect(
        "ChromatinAgent", {"findings": {"footprinting": footprinting}})
    ids = [b.id for b in blocks]
    assert "chromatin.differential_tf_footprinting" in ids
