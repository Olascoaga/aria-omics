"""ADR-080 (ADR-028 Slice 2): non-same-cell cross-modal synthesis, RunLedger C3
linkage, narrative reporting and the non-circular HypothesisAgent adapter.

Covers both the concordance-engine-backed path (``cross_modal_concordance``
ran) and the deterministic in-memory fallback (bulk RNA + bulk ATAC agents),
under the same W-CLAIM strict verification, C3 ledger linkage and C4 causal
guard the renderer enforces.
"""

from __future__ import annotations

import copy
import csv
from pathlib import Path

import pytest

from aria.agents.biological_synthesis_agent import BiologicalSynthesisAgent
from aria.agents.narrative.claim_compiler import compile_public_claims
from aria.agents.narrative.compose_prose import compose_block_prose
from aria.agents.narrative.evidence_verifier import verify_block_claim_support
from aria.agents.narrative.hypothesis import (
    CROSS_MODAL_UNPAIRED_COHORT,
    Hypothesis,
    build_signals_by_entity,
    check_devils_advocate,
    cross_modal_evidence,
    gather_evidence,
    rank_hypotheses,
)
from aria.agents.narrative.hypothesis.caveats import (
    CAVEAT_CODES,
    FOOTPRINT_ASSOCIATIVE,
    MOTIF_NOT_BINDING,
    PARALOG_MOTIF_AMBIGUITY,
    PEAK2GENE_ASSOCIATIVE,
    caveat_gloss,
)
from aria.agents.narrative.hypothesis.ranking import _normalized_effect
from aria.agents.narrative.reporting.sections import ReportSectionsMixin
from aria.agents.narrative.run_ledger import (
    build_run_ledger,
    ensure_report_ledger_nodes,
    link_claims_to_ledger,
    node_id_for,
    verify_blocks_against_ledger,
)
from aria.agents.narrative.synthesis.cross_modal_synthesis import (
    compose_cross_modal_discussion_blocks,
    detect_cross_modal_patterns,
)
from aria.agents.narrative.validators import find_causal_language
from aria.agents.narrative_agent import NarrativeAgent

CM_NODE = "ledger://cross_modal/concordance"
SYNTH_NODE = "ledger://report/biological_synthesis"


# ── fixtures ──────────────────────────────────────────────────────────────────

def _touch(path: Path, text: str = "x\n") -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return str(path)


def _engine_result(tmp_path: Path, *, flipped: bool = False) -> dict:
    out = tmp_path / "cross_modal"
    tables = {
        "de_da_concordance_csv": _touch(out / "cross_modal_de_da_concordance.csv"),
        "tf_triangulation_csv": _touch(out / "cross_modal_tf_triangulation.csv"),
        "shared_pathways_csv": _touch(out / "cross_modal_shared_pathways.csv"),
        "sample_matched_links_csv": _touch(
            out / "cross_modal_matched_sample_links.csv"),
    }
    return {
        "status": "success",
        "ran": True,
        "reason": None,
        "analysis": "cross_modal_concordance",
        "validation_level": "production",
        "evidence_tier": "associative",
        "associative_only": True,
        "rna_modality": "bulk_RNA",
        "atac_modality": "bulk_ATAC",
        "n_matched_contrasts": 1,
        "matched_contrasts": [{
            "contrast_id": "treated_vs_control",
            "rna_contrast_name": "treated_vs_control",
            "atac_contrast_name": (
                "control_vs_treated" if flipped else "treated_vs_control"),
            "atac_sign_flipped": flipped,
            "matching_granularity": "bulk_contrast",
            "de_da_summary": {
                "n_genes_tested_with_cis_peaks": 420,
                "n_gene_peak_pairs_tested": 650,
                "n_promoter_pairs_tested": 310,
                "n_both_significant": 28,
                "n_concordant_up": 15,
                "n_concordant_down": 9,
                "n_concordant_total": 24,
                "n_discordant_rna_up_atac_down": 2,
                "n_discordant_rna_down_atac_up": 2,
                "n_discordant_total": 4,
                "n_atac_da_rna_unchanged": 35,
                "n_rna_de_atac_unchanged": 41,
                "concordance_fraction": 0.8571,
                "fisher_odds_ratio": 33.75,
                "fisher_pvalue": 0.00012,
                "spearman_r": 0.54,
                "spearman_pvalue": 1.2e-12,
                "pearson_r": 0.51,
                "pearson_pvalue": 4.5e-11,
                "top_concordant_genes": [
                    {"gene": "GATA1", "peak": "chrX:48786000-48786800",
                     "feature": "Promoter", "distance_to_tss": -120,
                     "rna_log2fc": 3.2, "rna_padj": 1e-6, "atac_log2fc": 2.4,
                     "atac_padj": 2e-5, "quadrant": "concordant_up"},
                    {"gene": "ENSG00000102145", "peak": "chrX:1-2"},
                    {"gene": "TAL1", "peak": "chr1:47216000-47216500"},
                    {"gene": "KLF1", "peak": "chr19:12884000-12884400"},
                ],
                "top_discordant_genes": [{"gene": "MYB"}],
            },
            "tf_triangulation": {
                "ran": True,
                "n_tfs_evaluated": 12,
                "n_concordant_tfs": 4,
                "n_discordant_tfs": 1,
                "n_chromatin_shift_rna_unchanged": 5,
                "n_paralog_ambiguous": 2,
                "top_concordant_tfs": [
                    {"tf": "GATA1", "motif_id": "MA0035.4",
                     "source_layer": "motif_and_footprint",
                     "chromatin_effect": 2.1, "chromatin_padj": 0.001,
                     "chromatin_direction": "up", "rna_log2fc": 3.2,
                     "rna_padj": 1e-6, "rna_direction": "up",
                     "status": "concordant", "paralog_motif_ambiguity": True,
                     "paralog_family": "GATA"},
                    {"tf": "KLF1", "source_layer": "motif",
                     "chromatin_effect": -1.3, "status": "concordant",
                     "paralog_motif_ambiguity": False},
                    {"tf": "TAL1::TCF3", "chromatin_effect": 1.0,
                     "status": "concordant"},
                ],
            },
            "pathway_concordance": {
                "ran": True,
                "n_rna_terms": 18,
                "n_atac_terms": 11,
                "n_shared_terms": 5,
                "jaccard_index": 0.2083,
                "shared_terms": [
                    {"term": "erythrocyte differentiation (GO:0030218)",
                     "clean_term": "erythrocyte differentiation",
                     "database": "GO_BP", "rna_padj": 0.0004,
                     "atac_padj": 0.0031},
                    {"term": "heme biosynthetic process (GO:0006783)",
                     "clean_term": "heme biosynthetic process"},
                    # Causal vocabulary in a DB term name must never be named.
                    {"term": "process that activates erythropoiesis (GO:1)",
                     "clean_term": "process that activates erythropoiesis"},
                ],
            },
        }],
        "de_da_concordance": {},
        "tf_triangulation": {},
        "pathway_concordance": {},
        "sample_matched_links": {
            "ran": True,
            "reason": None,
            "n_shared_samples": 8,
            "shared_sample_ids": [f"D{i}" for i in range(1, 9)],
            "cis_distance_bp": 100000,
            "n_pairs_tested": 1200,
            "n_significant_links": 19,
            "n_positive_links": 17,
            "n_negative_links": 2,
            "permutation_null": {"n_permutations": 200,
                                 "empirical_fdr_at_padj": 0.02},
            "top_links": [
                {"gene": "GATA1", "peak": "chrX:48786000-48786800",
                 "distance_to_tss": -120, "spearman_r": 0.88,
                 "pearson_r": 0.86, "pvalue": 0.001, "padj": 0.012,
                 "direction": "positive"},
                {"gene": "chrX:1-2", "spearman_r": 0.9, "padj": 0.01},
                {"gene": "SPI1", "spearman_r": -0.71, "padj": 0.03},
                {"gene": "MYC", "spearman_r": 0.95, "padj": 0.4},
            ],
        },
        "tables": tables,
        "figures": {"cross_modal_starburst": _touch(
            out / "cross_modal_starburst.png")},
        "figures_svg": {"cross_modal_starburst": _touch(
            out / "cross_modal_starburst.svg", "<svg/>")},
        "stage_manifest_path": _touch(
            out / "cross_modal_concordance.stage_manifest.json", "{}"),
        "warnings": [],
    }


def _write_csv(path: Path, header: list[str], rows: list[list]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    return str(path)


def _bulk_rna_agent() -> dict:
    return {
        "status": "done",
        "findings": {
            "contrasts": [
                {
                    "name": "treated_vs_control", "numerator": "treated",
                    "denominator": "control", "status": "success",
                    "n_significant": 4,
                    "top_genes": [
                        {"symbol": "GATA1", "log2fc": 3.2, "padj": 1e-6},
                        {"symbol": "TAL1", "log2fc": 2.0, "padj": 1e-4},
                        {"symbol": "SPI1", "log2fc": -1.5, "padj": 1e-3},
                        {"symbol": "MYB", "log2fc": 1.1, "padj": 1e-2},
                    ],
                    "pathways": {"GO_BP": [
                        {"term": "erythrocyte differentiation (GO:0030218)"},
                        {"term": "heme biosynthetic process (GO:0006783)"},
                    ]},
                },
                {
                    # Incompatible arm (vehicle is not a chromatin level).
                    "name": "treated_vs_vehicle", "numerator": "treated",
                    "denominator": "vehicle", "status": "success",
                    "n_significant": 1,
                    "top_genes": [{"symbol": "GATA1", "log2fc": 1.0}],
                },
            ],
        },
    }


def _chromatin_agent(tmp_path: Path) -> dict:
    """Bulk ATAC run in the REVERSE orientation (control_vs_treated)."""
    da_csv = _write_csv(
        tmp_path / "atac" / "bulk_atac_da_control_vs_treated_full.csv",
        ["peak", "log2FoldChange", "padj", "significant"],
        [
            ["chr1:100-200", -2.0, 1e-5, "True"],   # GATA1, up in treated
            ["chr1:300-400", -1.0, 1e-3, "True"],   # TAL1, up in treated
            ["chr2:100-200", 1.5, 1e-4, "True"],    # SPI1, up in control
            ["chr3:100-200", 1.0, 1e-3, "True"],    # MYB: discordant
            ["chr4:100-200", -0.8, 1e-2, "True"],   # HBB: not RNA-DE
            ["chr5:100-200", 0.1, 0.9, "False"],    # not significant
        ],
    )
    ann_csv = _write_csv(
        tmp_path / "atac" / "bulk_atac_peak_annotation_control_vs_treated.csv",
        ["peak", "nearest_gene", "distance_to_tss", "feature"],
        [
            ["chr1:100-200", "GATA1", -120, "Promoter"],
            ["chr1:300-400", "TAL1", 5000, "Intronic"],
            ["chr2:100-200", "SPI1", -300, "Promoter"],
            ["chr3:100-200", "MYB", 20000, "Distal Intergenic"],
            ["chr4:100-200", "HBB", 100, "Promoter"],
            ["chr5:100-200", "ACTB", 100, "Promoter"],
        ],
    )
    return {
        "status": "done",
        "findings": {"bulk_ATAC": {"status": "done", "findings": {
            "differential_accessibility": {
                "data_type": "bulk_ATAC", "status": "success",
                "comparisons": [{
                    "test": "control", "reference": "treated",
                    "status": "success", "n_sig": 5,
                    "full_results_csv": da_csv,
                }],
            },
            "peak_annotation": {
                "ran": True, "status": "success", "data_type": "bulk_ATAC",
                "comparisons": [{"test": "control", "reference": "treated",
                                 "annotation_csv": ann_csv}],
            },
            "peak_ora": {
                "ran": True, "status": "success", "data_type": "bulk_ATAC",
                "comparisons": [{
                    "test": "control", "reference": "treated",
                    "pathways": {"GO_BP": [
                        {"term": "erythrocyte differentiation (GO:0030218)"},
                        {"term": "myeloid cell activation (GO:0002274)"},
                    ]},
                }],
            },
            "motifs": {
                "ran": True, "padj_max": 0.05,
                "per_group": {
                    "control_vs_treated::up_in_treated": {"top_motifs": [
                        {"name": "GATA1", "padj": 0.001,
                         "log2_enrichment": 1.4},
                        {"name": "TAL1::TCF3", "padj": 0.01,
                         "log2_enrichment": 0.9},
                    ]},
                    "control_vs_treated::up_in_control": {"top_motifs": [
                        {"name": "SPI1", "padj": 0.002, "log2_enrichment": 1.1},
                        {"name": "CEBPA", "padj": 0.01, "log2_enrichment": 0.7},
                        {"name": "NOTSIG", "padj": 0.2, "log2_enrichment": 0.1},
                    ]},
                },
            },
        }}},
    }


def _memory_results(tmp_path: Path) -> dict:
    return {"bulk_rna_agent": _bulk_rna_agent(),
            "chromatin_agent": _chromatin_agent(tmp_path)}


def _cross_modal_blocks(blocks) -> list:
    return [b for b in blocks if str(b.id).startswith("integration.cross_modal")]


def _assert_governed(blocks) -> None:
    assert blocks
    for b in blocks:
        prose = compose_block_prose(b)
        manifest = verify_block_claim_support(b, prose, strict=True)
        assert manifest["status"] == "supported", (b.id, manifest["unsupported"])
        texts = [b.claim, prose]
        texts += [f"{e.label} {e.value}" for e in b.evidence]
        texts += [c.text for c in b.caveats]
        for text in texts:
            assert find_causal_language(text) is None, (b.id, text)
    kept = NarrativeAgent._drop_unverifiable_synthesis_blocks(list(blocks))
    assert [b.id for b in kept] == [b.id for b in blocks]


# ── 1. RunLedger C3 ───────────────────────────────────────────────────────────

def test_ledger_reconciles_ran_cross_modal_node_with_provenance(tmp_path):
    cm = _engine_result(tmp_path)
    for payload in (cm, {"status": "done", "findings": cm}):
        ledger = build_run_ledger(
            {"design_intelligence": {"recommended": [
                "RNA ATAC concordance (DE x DA)"]}},
            {"cross_modal_concordance": payload},
        )
        assert "cross_modal" in ledger["modalities"]
        node = next(e for e in ledger["entries"] if e["node_id"] == CM_NODE)
        assert node["modality"] == "cross_modal"
        assert node["analysis"] == "concordance"
        assert node["node_type"] == "analysis_run"
        assert node["status"] == "ran"
        assert node["planned"] is True and node["divergence"] is False
        prov = node["provenance"]
        assert prov["kind"] == "ipc_script"
        assert prov["script"] == "aria/scripts/cross_modal_concordance.py"
        assert prov["rna_modality"] == "bulk_RNA"
        assert prov["atac_modality"] == "bulk_ATAC"
        assert prov["stage_manifest_path"] == cm["stage_manifest_path"]
        assert set(prov["tables"]) == set(cm["tables"])
        assert ledger["n_divergences"] == 0


@pytest.mark.parametrize("patch,expected", [
    ({"ran": False, "status": "skipped",
      "reason": "no_matched_contrasts"}, "skipped"),
    ({"ran": False, "status": "success"}, "skipped"),
    ({"ran": None, "status": "error", "reason": "boom"}, "error"),
    ({"ran": None, "status": "success"}, "not_run"),
])
def test_ledger_non_run_statuses_are_honest(tmp_path, patch, expected):
    cm = _engine_result(tmp_path)
    cm.update(patch)
    ledger = build_run_ledger(
        {"design_intelligence": {"recommended": ["cross-modal concordance"]}},
        {"cross_modal_concordance": cm},
    )
    node = next(e for e in ledger["entries"] if e["node_id"] == CM_NODE)
    assert node["status"] == expected
    assert node["divergence"] is True
    # A non-run engine never produces public cross-modal claims.
    blocks = BiologicalSynthesisAgent().synthesize(
        {"cross_modal_concordance": cm}, {})
    assert _cross_modal_blocks(blocks) == []


def test_node_id_routing_for_cross_modal_and_in_memory_aliases():
    assert node_id_for("cross_modal", "concordance") == CM_NODE
    assert node_id_for("cross-modal", "cross_modal_concordance") == CM_NODE
    assert node_id_for("integrated", "cross_modal_synthesis") == SYNTH_NODE
    assert node_id_for("integrated", "cross_modal_limitations") == SYNTH_NODE
    # Generic RNA/ATAC fallbacks are unaffected.
    assert node_id_for("bulk", "differential_expression") == (
        "ledger://bulk/differential_expression")
    assert node_id_for("chromatin", "motif_enrichment") == (
        "ledger://chromatin/motif_enrichment")


# ── 2-3. Engine-backed Slice 2 ────────────────────────────────────────────────

def test_engine_backed_blocks_are_governed_and_link_to_ipc_node(tmp_path):
    cm = _engine_result(tmp_path, flipped=True)
    ar = {"cross_modal_concordance": {"status": "done", "findings": cm}}
    blocks = BiologicalSynthesisAgent().synthesize(ar, {})
    cm_blocks = _cross_modal_blocks(blocks)
    ids = [b.id for b in cm_blocks]
    assert ids == [
        "integration.cross_modal_de_da",
        "integration.cross_modal_tf",
        "integration.cross_modal_pathways",
        "integration.cross_modal_sample_links",
        "integration.cross_modal_limitations",
    ]
    assert all(b.modality == "cross_modal" and b.analysis == "concordance"
               for b in cm_blocks)
    _assert_governed(cm_blocks)

    de_da = cm_blocks[0]
    for token in ("650", "24", "4", "35", "0.54", "33.8", "GATA1", "TAL1",
                  "KLF1"):
        assert token in de_da.claim
    assert "ENSG" not in de_da.claim and "chrX" not in de_da.claim
    fig = de_da.figures[0]
    assert fig["path"].endswith("cross_modal_starburst.png")
    assert fig["svg_path"].endswith("cross_modal_starburst.svg")
    assert {Path(t["path"]).name for t in de_da.tables} == {
        Path(p).name for p in cm["tables"].values()}

    tf = cm_blocks[1]
    assert "12" in tf.claim and "paralog" in tf.claim
    assert "TAL1::TCF3" not in tf.claim  # composite motif never named
    pw = cm_blocks[2]
    assert "erythrocyte differentiation" in pw.claim
    assert "GO:" not in pw.claim and "activates" not in pw.claim
    links = cm_blocks[3]
    assert "8" in links.claim and "1200" in links.claim and "19" in links.claim
    lim = cm_blocks[4]
    assert lim.block_type == "limitation"
    assert "separate libraries" in lim.claim
    assert "control_vs_treated" in lim.claim  # reversed-contrast disclosure

    ledger = build_run_ledger({}, ar)
    ensure_report_ledger_nodes(ledger, cm_blocks)
    compiled = compile_public_claims(cm_blocks, {}, run_ledger=ledger)
    assert compiled.withheld == []
    assert len(compiled.claims) == len(cm_blocks)
    linkage = link_claims_to_ledger(copy.deepcopy(compiled.claims), ledger)
    assert linkage["n_violations"] == 0
    assert {c["ledger_node_id"] for c in compiled.claims} == {CM_NODE}
    assert verify_blocks_against_ledger(compiled.blocks, ledger)["status"] == "ok"


def test_engine_claims_withheld_when_ledger_node_skipped(tmp_path):
    cm = _engine_result(tmp_path)
    blocks = _cross_modal_blocks(BiologicalSynthesisAgent().synthesize(
        {"cross_modal_concordance": cm}, {}))
    skipped = dict(cm, ran=False, status="skipped")
    ledger = build_run_ledger({}, {"cross_modal_concordance": skipped})
    compiled = compile_public_claims(blocks, {}, run_ledger=ledger)
    assert compiled.claims == []
    assert len(compiled.withheld) == len(blocks)


def test_engine_honest_skip_reason_and_missing_layers(tmp_path):
    cm = _engine_result(tmp_path)
    cm["sample_matched_links"] = {
        "ran": False,
        "reason": "fewer than 6 matched samples across disjoint namespaces",
    }
    mc = cm["matched_contrasts"][0]
    mc["tf_triangulation"] = {"ran": False, "reason": "no motifs"}
    mc["pathway_concordance"] = {"ran": True, "n_rna_terms": 0,
                                 "n_atac_terms": 3, "n_shared_terms": 0}
    for key in ("figures", "figures_svg"):
        cm[key] = {"cross_modal_starburst": str(tmp_path / "missing.png")}
    blocks = _cross_modal_blocks(BiologicalSynthesisAgent().synthesize(
        {"cross_modal_concordance": cm}, {}))
    assert [b.id for b in blocks] == ["integration.cross_modal_de_da",
                                      "integration.cross_modal_limitations"]
    assert blocks[0].figures == []  # never attaches a missing file
    assert "fewer than 6 matched samples" in blocks[1].claim
    _assert_governed(blocks)


def test_multiple_matched_contrasts_get_distinct_block_ids(tmp_path):
    cm = _engine_result(tmp_path)
    second = copy.deepcopy(cm["matched_contrasts"][0])
    second["contrast_id"] = "day7_vs_day0"
    cm["matched_contrasts"].append(second)
    blocks = _cross_modal_blocks(BiologicalSynthesisAgent().synthesize(
        {"cross_modal_concordance": cm}, {}))
    ids = [b.id for b in blocks]
    assert len(ids) == len(set(ids))
    assert "integration.cross_modal_de_da.treated_vs_control" in ids
    assert "integration.cross_modal_de_da.day7_vs_day0" in ids
    assert sum(1 for b in blocks if b.figures) == 1
    _assert_governed(blocks)


# ── 2-3. In-memory fallback ───────────────────────────────────────────────────

def test_in_memory_fallback_sign_alignment_and_incompatible_rejection(tmp_path):
    ar = _memory_results(tmp_path)
    patterns = detect_cross_modal_patterns(ar, {})
    assert patterns["active"] is True
    assert patterns["backed_by_concordance_engine"] is False
    assert patterns["rna_modality"] == "bulk_RNA"
    assert patterns["atac_modality"] == "bulk_ATAC"
    assert [r["rna_contrast_name"] for r in patterns["rejected_contrasts"]] == [
        "treated_vs_vehicle"]
    (mc,) = patterns["contrasts"]
    assert mc["atac_sign_flipped"] is True
    assert mc["atac_contrast_name"] == "control_vs_treated"
    d = mc["de_da"]
    assert (d["n_da_peaks_annotated"], d["n_both_significant"]) == (5, 4)
    assert (d["n_concordant_up"], d["n_concordant_down"]) == (2, 1)
    assert d["n_discordant_total"] == 1
    assert [g["gene"] for g in d["top_concordant_genes"]] == [
        "GATA1", "SPI1", "TAL1"]
    t = mc["tf"]
    assert t["n_tfs_evaluated"] == 5
    assert t["n_concordant_tfs"] == 3
    assert t["n_discordant_tfs"] == 0
    assert t["n_chromatin_shift_rna_unchanged"] == 2
    assert t["n_paralog_ambiguous"] == 2  # TAL1::TCF3 dimer subunits
    p = mc["pathways"]
    assert (p["n_rna_terms"], p["n_atac_terms"], p["n_shared_terms"]) == (2, 2, 1)

    blocks = _cross_modal_blocks(BiologicalSynthesisAgent().synthesize(ar, {}))
    assert [b.id for b in blocks] == [
        "integration.cross_modal_de_da",
        "integration.cross_modal_tf",
        "integration.cross_modal_pathways",
        "integration.cross_modal_limitations",
    ]
    assert all(b.modality == "integrated" for b in blocks)
    assert {b.analysis for b in blocks[:-1]} == {"cross_modal_synthesis"}
    assert blocks[-1].analysis == "cross_modal_limitations"
    assert all(not b.figures for b in blocks)
    assert "control_vs_treated" in blocks[-1].claim
    assert "1 RNA contrast(s) had no orientation-compatible" in blocks[-1].claim
    _assert_governed(blocks)

    ledger = build_run_ledger({}, ar)
    assert CM_NODE not in {e["node_id"] for e in ledger["entries"]}
    ensure_report_ledger_nodes(ledger, blocks)
    node = next(e for e in ledger["entries"] if e["node_id"] == SYNTH_NODE)
    assert node["status"] == "ran"
    assert node["node_type"] == "derived_report_analysis"
    compiled = compile_public_claims(blocks, {}, run_ledger=ledger)
    assert compiled.withheld == []
    linkage = link_claims_to_ledger(copy.deepcopy(compiled.claims), ledger)
    assert linkage["n_violations"] == 0
    assert {c["ledger_node_id"] for c in compiled.claims} == {SYNTH_NODE}


def test_in_memory_same_orientation_is_not_flipped(tmp_path):
    ar = _memory_results(tmp_path)
    chrom = ar["chromatin_agent"]["findings"]["bulk_ATAC"]["findings"]
    chrom["differential_accessibility"]["comparisons"][0].update(
        {"test": "treated", "reference": "control"})
    for rec in (chrom["peak_annotation"]["comparisons"]
                + chrom["peak_ora"]["comparisons"]):
        rec.update({"test": "treated", "reference": "control"})
    chrom["motifs"]["per_group"] = {
        k.replace("control_vs_treated", "treated_vs_control"): v
        for k, v in chrom["motifs"]["per_group"].items()
    }
    (mc,) = detect_cross_modal_patterns(ar)["contrasts"]
    assert mc["atac_sign_flipped"] is False
    # Same raw CSV signs read in the forward orientation invert the agreement.
    assert mc["de_da"]["n_concordant_total"] == 1
    assert mc["de_da"]["n_discordant_total"] == 3
    blocks = _cross_modal_blocks(BiologicalSynthesisAgent().synthesize(ar, {}))
    assert "sign-aligned" not in blocks[-1].claim
    _assert_governed(blocks)


def test_paired_same_cell_integration_does_not_duplicate_slice2(tmp_path):
    ar = _memory_results(tmp_path)
    ar["integration_agent"] = {"status": "done", "findings": {
        "paired_qc": {"status": "success"}}}
    patterns = detect_cross_modal_patterns(ar)
    assert patterns["active"] is False
    assert patterns["skip_reason"] == (
        "paired_same_cell_integration_owns_triangulation")
    assert _cross_modal_blocks(BiologicalSynthesisAgent().synthesize(ar, {})) == []
    # An explicit engine run still synthesizes (no duplicate of multiome.*).
    ar["cross_modal_concordance"] = _engine_result(tmp_path)
    assert _cross_modal_blocks(BiologicalSynthesisAgent().synthesize(ar, {}))


def test_single_modality_runs_never_fire(tmp_path):
    for ar in ({"bulk_rna_agent": _bulk_rna_agent()},
               {"chromatin_agent": _chromatin_agent(tmp_path)}, {}):
        assert detect_cross_modal_patterns(ar)["active"] is False
        assert _cross_modal_blocks(
            BiologicalSynthesisAgent().synthesize(ar, {})) == []
    assert compose_cross_modal_discussion_blocks({"active": False}) == []


def test_slice2_failure_is_isolated(monkeypatch, tmp_path):
    import aria.agents.biological_synthesis_agent as mod

    def _boom(*_a, **_k):
        raise RuntimeError("boom")

    monkeypatch.setattr(mod, "detect_cross_modal_patterns", _boom)
    blocks = BiologicalSynthesisAgent().synthesize(_memory_results(tmp_path), {})
    assert _cross_modal_blocks(blocks) == []
    assert any(b.id.startswith("integration.") for b in blocks)  # Slice 1a kept


# ── 5. Reporting: staging, QC row, conflict notice, Methods ───────────────────

class _Sections(ReportSectionsMixin):
    pass


def test_stage_artifacts_copies_cross_modal_tables_and_starburst(tmp_path):
    cm = _engine_result(tmp_path)
    report_dir = tmp_path / "report"
    _Sections()._stage_artifacts(
        {"cross_modal_concordance": {"status": "done", "findings": cm}},
        report_dir,
    )
    tables = sorted(p.name for p in (report_dir / "tables").iterdir())
    assert tables == sorted([
        "cross_modal_de_da_concordance.csv",
        "cross_modal_tf_triangulation.csv",
        "cross_modal_shared_pathways.csv",
        "cross_modal_matched_sample_links.csv",
    ])
    figures = sorted(p.name for p in (report_dir / "figures").iterdir())
    assert figures == ["cross_modal_starburst.png", "cross_modal_starburst.svg"]


def test_qc_row_and_conflict_notice(tmp_path):
    cm = _engine_result(tmp_path)
    grouped = {"high": [], "medium": [], "low": [], "insufficient": []}
    html = _Sections()._build_qc_section(grouped, {"cross_modal_concordance": cm})
    assert "Cross-Modal Concordance (RNA × ATAC)" in html
    assert "1 matched contrast(s)" in html
    assert "24 concordant DE×DA pair(s)" in html
    assert "8 matched sample(s)" in html
    notice = ReportSectionsMixin._public_conflict_notice(
        {"cross_modal_concordance": cm})
    assert "not applicable" not in notice
    skipped = dict(cm, ran=False, status="skipped")
    assert "not applicable" in ReportSectionsMixin._public_conflict_notice(
        {"cross_modal_concordance": skipped})
    html_skipped = _Sections()._build_qc_section(
        grouped, {"cross_modal_concordance": skipped})
    assert "Cross-Modal Concordance" not in html_skipped


def test_methods_section_is_html_safe_and_non_causal(tmp_path):
    cm = _engine_result(tmp_path, flipped=True)
    cm["rna_modality"] = "<b>bulk_RNA</b>"
    agent = NarrativeAgent.__new__(NarrativeAgent)
    text = agent._write_methods_section(
        {"organism": "human", "genome": "GRCh38"},
        {"cross_modal_concordance": cm}, [])
    assert "**Cross-modal concordance (RNA × ATAC)**" in text
    paragraph = text.split("**Cross-modal concordance (RNA × ATAC)**", 1)[1]
    assert "<" not in paragraph and ">" not in paragraph
    assert find_causal_language(paragraph) is None
    for phrase in ("100000 bp", "Fisher", "Spearman", "Pearson", "paralog",
                   "hypergeometric", "log2(CPM + 1)", "8 biological units",
                   "200 sample-label permutations", "Benjamini-Hochberg",
                   "1 sign-aligned contrast(s)", "multiplied by -1"):
        assert phrase in paragraph, phrase

    no_links = dict(cm, sample_matched_links={"ran": False, "reason": "x"})
    text2 = agent._write_methods_section({}, {"cross_modal_concordance": no_links}, [])
    assert "log2(CPM + 1)" not in text2
    skipped = dict(cm, ran=False)
    assert "Cross-modal concordance" not in agent._write_methods_section(
        {}, {"cross_modal_concordance": skipped}, [])


# ── 6. HypothesisAgent non-circular adapter ───────────────────────────────────

def test_caveat_code_registered_with_gloss():
    assert CROSS_MODAL_UNPAIRED_COHORT == "cross_modal_unpaired_cohort"
    assert CROSS_MODAL_UNPAIRED_COHORT in CAVEAT_CODES
    gloss = caveat_gloss(CROSS_MODAL_UNPAIRED_COHORT)
    assert "non-same-cell" in gloss
    assert find_causal_language(gloss) is None


def test_adapter_gated_on_ledger_node(tmp_path):
    cm = _engine_result(tmp_path)
    ar = {"cross_modal_concordance": cm}
    assert cross_modal_evidence({}, None) == []
    skipped_ledger = build_run_ledger(
        {}, {"cross_modal_concordance": dict(cm, ran=False, status="skipped")})
    assert cross_modal_evidence(ar, skipped_ledger) == []
    assert cross_modal_evidence(ar, {"entries": []}) == []  # node not_run
    not_run = dict(cm, ran=None)
    assert cross_modal_evidence({"cross_modal_concordance": not_run}) == []

    signals = cross_modal_evidence(ar, build_run_ledger({}, ar))
    links = [s for s in signals if s.measure == "sample_matched_peak2gene_corr"]
    tfs = [s for s in signals if s.measure == "cross_modal_tf_concordance"]
    # Peak coordinate rejected; padj above alpha dropped.
    assert [s.entity for s in links] == ["GATA1", "SPI1"]
    assert links[0].direction == "up" and links[1].direction == "down"
    assert all(s.audited_node_ref == CM_NODE and s.modality == "cross_modal"
               for s in signals)
    assert links[0].context == "matched_samples_n=8"
    assert links[0].caveats_inherited == [PEAK2GENE_ASSOCIATIVE,
                                          CROSS_MODAL_UNPAIRED_COHORT]
    assert not links[0].independence_key
    assert _normalized_effect(links[0]) == pytest.approx(0.88)
    # Composite motif labels are never invented into component genes.
    assert [s.entity for s in tfs] == ["GATA1", "KLF1"]
    gata1 = tfs[0]
    assert gata1.entity_kind == "tf_motif"
    assert gata1.value == 2.1 and gata1.direction == "up"
    assert gata1.context == "treated_vs_control"
    assert gata1.independence_key == "chromatin_tf_lineage::gata1"
    assert gata1.caveats_inherited == [MOTIF_NOT_BINDING,
                                       CROSS_MODAL_UNPAIRED_COHORT,
                                       PARALOG_MOTIF_AMBIGUITY,
                                       FOOTPRINT_ASSOCIATIVE]
    assert tfs[1].direction == "down"
    assert tfs[1].caveats_inherited == [MOTIF_NOT_BINDING,
                                        CROSS_MODAL_UNPAIRED_COHORT]


def _full_results(tmp_path: Path) -> dict:
    ar = _memory_results(tmp_path)
    ar["cross_modal_concordance"] = _engine_result(tmp_path)
    return ar


def _sig(signals, *, entity, measure_prefix):
    return next(s for s in signals
                if s.entity == entity and s.measure.startswith(measure_prefix))


def test_gather_evidence_never_double_counts_derived_tf_concordance(tmp_path):
    ar = _full_results(tmp_path)
    ledger = build_run_ledger({}, ar)
    signals = gather_evidence(ar, ledger, {})
    rna = _sig(signals, entity="GATA1", measure_prefix="log2fc")
    atac = _sig(signals, entity="GATA1", measure_prefix="motif_enrich")
    cm_tf = _sig(signals, entity="GATA1",
                 measure_prefix="cross_modal_tf_concordance")
    assert rna.modality == "bulk_RNA" and atac.modality == "bulk_ATAC"
    assert atac.independence_key == cm_tf.independence_key == (
        "chromatin_tf_lineage::gata1")
    assert not rna.independence_key
    # A chromatin TF without a cross-modal concordance keeps its own line.
    cebpa = _sig(signals, entity="CEBPA", measure_prefix="motif_enrich")
    assert not cebpa.independence_key

    def _hyp(hid, cited):
        return Hypothesis(
            id=hid, mechanism="GATA1 co-varies across assays.",
            entities=["GATA1"],
            observed_claims=[{"signal_id": s.signal_id,
                              "stated_direction": s.direction}
                             for s in cited],
        )

    h_all = _hyp("h_all", [rna, atac, cm_tf])
    h_cm = _hyp("h_cm", [cm_tf])
    h_rna_cm = _hyp("h_rna_cm", [rna, cm_tf])
    h_atac_cm = _hyp("h_atac_cm", [atac, cm_tf])
    rank_hypotheses([h_all, h_cm, h_rna_cm, h_atac_cm], signals)
    assert h_all.rank_evidence["n_independent_lines"] == 2  # NOT 3
    assert h_cm.rank_evidence["n_independent_lines"] == 1
    assert h_rna_cm.rank_evidence["n_independent_lines"] == 2
    assert h_atac_cm.rank_evidence["n_independent_lines"] == 1

    # Without the cross-modal engine the single-assay keys are untouched.
    plain = gather_evidence(_memory_results(tmp_path),
                            build_run_ledger({}, _memory_results(tmp_path)), {})
    assert not _sig(plain, entity="GATA1",
                    measure_prefix="motif_enrich").independence_key


def test_devils_advocate_requires_unpaired_cohort_acknowledgement(tmp_path):
    ar = _full_results(tmp_path)
    ledger = build_run_ledger({}, ar)
    signals = gather_evidence(ar, ledger, {})
    index = build_signals_by_entity(signals)
    cm_tf = _sig(signals, entity="GATA1",
                 measure_prefix="cross_modal_tf_concordance")
    base_confounds = [
        "motif_not_binding", "paralog_motif_ambiguity",
        "footprint_associative", "peak2gene_associative",
    ]
    hyp = Hypothesis(
        id="h1", mechanism="GATA1 co-varies across assays.",
        entities=["GATA1"],
        observed_claims=[{"signal_id": cm_tf.signal_id,
                          "stated_direction": "up"}],
        devils_advocate={
            "simpler_explanation": "Shared condition effects in both assays.",
            "confounds": list(base_confounds),
        },
    )
    result = check_devils_advocate(hyp, index)
    assert result.passed is False
    assert CROSS_MODAL_UNPAIRED_COHORT in result.reason

    for acknowledgement in (
        CROSS_MODAL_UNPAIRED_COHORT,
        "the assays are non-same-cell (separate libraries)",
    ):
        hyp.devils_advocate["confounds"] = base_confounds + [acknowledgement]
        assert check_devils_advocate(hyp, index).passed is True
