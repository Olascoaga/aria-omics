"""BiologicalSynthesisAgent (Slice 1): evidence-governed integrated discussion.

The agent must organize evidence, never invent it. These guards mirror the
design spec: no fabricated pathways, no causal language, cross-modal only when
the modality is present, mandatory limitations, conflicts surfaced honestly, and
every composed claim mapped to its evidence card (strict verification).
"""

from aria.agents.narrative.synthesis.pattern_detector import (
    detect_bulk_patterns,
    detect_scrna_patterns,
)
from aria.agents.narrative.synthesis.discussion_composer import (
    compose_discussion_blocks,
    compose_scrna_discussion_blocks,
)
from aria.agents.narrative.claim_compiler import annotate_claim_tiers
from aria.agents.narrative.evidence_verifier import verify_block_claim_support
from aria.agents.narrative.validators import find_causal_language
from aria.agents.biological_synthesis_agent import BiologicalSynthesisAgent


def _contrast(name, num, den, up, down, terms, symbols=None):
    ids = list(up) + list(down)
    syms = symbols if symbols is not None else ids
    return {
        "name": name, "numerator": num, "denominator": den, "status": "success",
        "n_significant": len(ids), "all_sig_gene_ids": ids, "all_sig_genes": syms,
        "up_gene_ids": list(up), "down_gene_ids": list(down),
        "pathways": {"GO_BP": [{"term": t} for t in terms]},
        "power_estimate_at_lfc_min": 0.8,
    }


def _two_contrasts():
    # A: up g1-6, down g7-10 ; B: up g5-9, down g10-14 ; shared = g5..g10 (6)
    A = _contrast("KOa vs WT", "KOa", "WT",
                  [f"g{i}" for i in range(1, 7)], [f"g{i}" for i in range(7, 11)],
                  ["GO:1", "GO:2"])
    B = _contrast("KOb vs WT", "KOb", "WT",
                  [f"g{i}" for i in range(5, 10)], [f"g{i}" for i in range(10, 15)],
                  ["GO:2", "GO:3"])
    return [A, B]


def _scrna_agent_result():
    return {
        "status": "done",
        "findings": {
            "scRNA": {
                "findings": {
                    "qc": {"n_cells_after": 1200},
                    "cell_types": {
                        "cell_types": {
                            "GroupA": {"cell_type": "GroupA"},
                            "GroupB": {"cell_type": "GroupB"},
                        }
                    },
                    "differential_abundance": {
                        "per_comparison": {
                            "condition_a_vs_condition_b": {
                                "status": "success",
                                "n_significant": 1,
                                "per_cell_type": [{
                                    "name": "GroupA",
                                    "significant": True,
                                    "direction": "up",
                                }],
                            }
                        }
                    },
                    "pseudobulk_de": {
                        "per_group": {
                            "GroupA": {
                                "per_comparison": {
                                    "condition_a_vs_condition_b": {
                                        "status": "success",
                                        "n_significant": 140,
                                        "n_up": 90,
                                        "n_down": 50,
                                        "corrected_for_composition": True,
                                        "power_estimate_at_lfc_min": 0.82,
                                        "top_genes": [
                                            {"gene": "GENE_A", "log2fc": 2.4},
                                            {"gene": "GENE_B", "log2fc": -1.1},
                                        ],
                                    }
                                }
                            }
                        }
                    },
                    "pseudobulk_pathways": {
                        "per_cluster": {
                            "GroupA::condition_a_vs_condition_b": {
                                "n_significant": 2,
                                "results": {
                                    "GO_BP": [{
                                        "term": "pathway_alpha",
                                        "adjusted_p": 1e-5,
                                    }]
                                },
                            }
                        }
                    },
                    "cell_communication": {
                        "status": "success",
                        "n_interactions": 3,
                        "n_cell_types": 2,
                        "n_autocrine_dropped": 1,
                    },
                    "trajectory": {
                        "status": "success",
                        "paga": {"n_connections": 6, "n_strong": 1},
                        "pseudotime": {"computed": True},
                    },
                }
            }
        },
    }


# ── detector math ────────────────────────────────────────────────────────────

def test_detector_counts_shared_specific_direction_and_terms():
    p = detect_bulk_patterns(_two_contrasts())
    x = p["cross_contrast"][0]
    assert x["shared_reference"] == "WT"
    assert x["n_shared_genes"] == 6                 # g5..g10
    assert x["n_specific_a"] == 4 and x["n_specific_b"] == 4
    # concordant: up∩up {g5,g6} + down∩down {g10} = 3 ; discordant {g7,g8,g9} = 3
    assert x["n_direction_concordant"] == 3
    assert x["n_direction_discordant"] == 3
    assert x["n_shared_terms"] == 1                 # GO:2


def test_no_pairs_when_single_contrast():
    p = detect_bulk_patterns([_two_contrasts()[0]])
    assert p["cross_contrast"] == []


# ── governance guards ────────────────────────────────────────────────────────

def test_every_claim_passes_strict_evidence_verification():
    blocks = compose_discussion_blocks(detect_bulk_patterns(_two_contrasts()))
    annotate_claim_tiers(blocks, {})
    assert blocks
    for b in blocks:                                 # Test 6: evidence mapping
        verify_block_claim_support(b, strict=True)


def test_no_causal_language_in_any_claim():
    blocks = compose_discussion_blocks(detect_bulk_patterns(_two_contrasts()))
    for b in blocks:                                 # Test 2: no causality
        assert find_causal_language(b.claim) is None, b.claim


def test_claims_are_capped_at_associative():
    blocks = compose_discussion_blocks(detect_bulk_patterns(_two_contrasts()))
    annotate_claim_tiers(blocks, {})
    for b in blocks:
        tier = b.metadata.get("claim", {}).get("tier")
        assert tier in {"descriptive", "associative"}, (b.title, tier)


def test_convergent_claim_names_the_shared_processes_not_just_counts():
    # Integrative, not a summary: the claim must say what the convergence points
    # to (the shared enriched processes), and those names must be evidence-backed.
    A = _contrast("KOa vs WT", "KOa", "WT", ["g1", "g2", "g3"], ["g4"],
                  ["response to zinc ion", "extracellular matrix organization"])
    B = _contrast("KOb vs WT", "KOb", "WT", ["g2", "g3"], ["g4", "g5"],
                  ["response to zinc ion", "cholesterol biosynthetic process"])
    blocks = compose_discussion_blocks(detect_bulk_patterns([A, B]))
    conv = next(b for b in blocks if b.analysis == "convergent_evidence")
    assert "points to" in conv.claim
    assert "response to zinc ion" in conv.claim          # the shared process named
    verify_block_claim_support(conv, strict=True)        # still evidence-backed
    div = next(b for b in blocks if b.analysis == "divergent_evidence")
    assert "uniquely engages" in div.claim
    # each contrast's specific process is named
    assert "extracellular matrix organization" in div.claim
    assert "cholesterol biosynthetic process" in div.claim
    verify_block_claim_support(div, strict=True)


def test_main_pattern_block_leads_and_names_genes_and_processes():
    # The synthesis must open with a headline "main pattern" and name the top
    # shared genes by symbol (reviewers expect the gene names).
    A = _contrast("KOa vs WT", "KOa", "WT", ["i1", "i2", "i3"], ["i4"],
                  ["response to zinc ion"],
                  symbols=["IGFBP5", "COL3A1", "ANKRD1", "RFX3"])
    B = _contrast("KOb vs WT", "KOb", "WT", ["i1", "i2"], ["i4", "i5"],
                  ["response to zinc ion"],
                  symbols=["IGFBP5", "COL3A1", "KITLG", "GABRP"])
    blocks = compose_discussion_blocks(detect_bulk_patterns([A, B]))
    annotate_claim_tiers(blocks, {})
    assert blocks[0].id == "integration.main_pattern"        # headline first
    assert "dominant integrated pattern" in blocks[0].claim.lower()
    verify_block_claim_support(blocks[0], strict=True)
    conv = next(b for b in blocks if b.analysis == "convergent_evidence")
    # shared genes IGFBP5/COL3A1 (i1,i2) named by symbol, ranked
    assert "IGFBP5" in conv.claim and "COL3A1" in conv.claim
    assert "KITLG" not in conv.claim                          # not shared -> not named
    verify_block_claim_support(conv, strict=True)


def test_no_pathways_means_no_enrichment_claim():
    # Test 1: pathways empty -> no integrated-signal block, no "enriched term" prose
    A = _contrast("KOa vs WT", "KOa", "WT", ["g1", "g2"], ["g3"], [])
    B = _contrast("KOb vs WT", "KOb", "WT", ["g2", "g3"], ["g4"], [])
    blocks = compose_discussion_blocks(detect_bulk_patterns([A, B]))
    ids = {b.id for b in blocks}
    assert "integration.signal" not in ids
    assert not any("enriched term" in b.claim.lower() for b in blocks)


def test_discordant_evidence_is_surfaced_not_smoothed():
    # Test 5: shared genes moving opposite directions are called discordant.
    blocks = compose_discussion_blocks(detect_bulk_patterns(_two_contrasts()))
    div = next(b for b in blocks if b.analysis == "divergent_evidence")
    assert "discordant evidence" in div.claim.lower()


def test_limitations_block_is_mandatory_when_claims_made():
    # Test 4: a limitations block is always present alongside integrated claims.
    blocks = compose_discussion_blocks(detect_bulk_patterns(_two_contrasts()))
    assert any(b.id == "integration.limitations" for b in blocks)


def test_cross_modal_not_claimed_without_the_modality():
    # Test 3: Slice 1 reads only bulk; nothing about chromatin/accessibility.
    blocks = compose_discussion_blocks(detect_bulk_patterns(_two_contrasts()))
    blob = " ".join(b.claim.lower() for b in blocks)
    for forbidden in ("accessibility", "chromatin", "atac", "motif", "peak"):
        assert forbidden not in blob


def test_agent_returns_empty_without_bulk_results():
    assert BiologicalSynthesisAgent().synthesize({}, {}) == []
    assert BiologicalSynthesisAgent().synthesize(
        {"bulk_rna_agent": {"findings": {"contrasts": []}}}, {}) == []


def test_scrna_synthesis_detects_measured_layers_only():
    p = detect_scrna_patterns(_scrna_agent_result())
    assert p["modalities_present"] == ["scRNA-seq"]
    assert p["n_cells"] == 1200
    assert p["strongest_pseudobulk"]["group"] == "GroupA"
    assert p["strongest_pseudobulk"]["n_de"] == 140
    assert p["strongest_pseudobulk"]["n_pathway_terms"] == 2
    assert p["n_abundance_shifts"] == 1
    assert p["cellcomm"]["n_interactions"] == 3
    assert p["trajectory"]["n_connections"] == 6


def test_scrna_synthesis_blocks_pass_strict_verification():
    blocks = compose_scrna_discussion_blocks(
        detect_scrna_patterns(_scrna_agent_result())
    )
    annotate_claim_tiers(blocks, {})
    ids = {b.id for b in blocks}
    assert "integration.scrna.main_pattern" in ids
    assert "integration.scrna.context_layers" in ids
    assert "integration.scrna.limitations" in ids
    for block in blocks:
        assert find_causal_language(block.claim) is None, block.claim
        verify_block_claim_support(block, strict=True)


def test_agent_synthesizes_scRNA_without_bulk():
    blocks = BiologicalSynthesisAgent().synthesize(
        {"scrna_agent": _scrna_agent_result()}, {}
    )
    assert any(b.id == "integration.scrna.main_pattern" for b in blocks)
    blob = " ".join(b.claim.lower() for b in blocks)
    for forbidden in ("accessibility", "chromatin", "atac", "motif", "peak"):
        assert forbidden not in blob


# ── real-data golden (runs on Samael's machine; skips in CI) ──────────────────

def test_real_report_overlap_golden():
    import csv
    import os
    import pytest

    T = ("/home/medusa/.aria/reports/"
         "aria_20260604_170909_bmal1_reverba_h9cells_-39a/tables")
    if not os.path.isdir(T):
        pytest.skip("real report tables not present")

    def contrast(name, num, den, de, pw):
        ids, up, down = [], [], []
        for r in csv.DictReader(open(f"{T}/{de}"), delimiter="\t"):
            g = r["gene_id"]
            ids.append(g)
            (up if float(r["log2FoldChange"]) > 0 else down).append(g)
        paths = {"GO_BP": [{"term": r.get("term") or list(r.values())[0]}
                           for r in csv.DictReader(open(f"{T}/{pw}"), delimiter="\t")]}
        return {"name": name, "numerator": num, "denominator": den,
                "status": "success", "n_significant": len(ids),
                "all_sig_gene_ids": ids, "up_gene_ids": up, "down_gene_ids": down,
                "pathways": paths, "power_estimate_at_lfc_min": 0.76}

    contrasts = [
        contrast("BMAL1_KO vs WT", "BMAL1_KO", "WT",
                 "bmal1_ko_vs_wt_de_genes.tsv", "bmal1_ko_vs_wt_pathways.tsv"),
        contrast("REV-ERBa_KO vs WT", "REV-ERBa_KO", "WT",
                 "rev_erba_ko_vs_wt_de_genes.tsv", "rev_erba_ko_vs_wt_pathways.tsv"),
    ]
    x = detect_bulk_patterns(contrasts)["cross_contrast"][0]
    assert x["n_shared_genes"] == 268
    assert x["n_direction_concordant"] == 259
    assert x["n_direction_discordant"] == 9
    assert x["n_shared_terms"] == 16
