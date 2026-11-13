"""V48 Slice I7 tests: Cross-Modal Evidence Graph, IntegrationNarrator, C1-C4/Claim 6, and Dual PNG+SVG Figures."""

from __future__ import annotations

import copy
import csv
import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from aria.agents import _narrative_integration
from aria.agents.narrative.claim_compiler import compile_public_claims
from aria.agents.narrative.compose_prose import compose_block_prose
from aria.agents.narrative.devils_advocate import build_devils_advocate
from aria.agents.narrative.evidence_verifier import (
    EvidenceVerificationError,
    verify_block_claim_support,
)
from aria.agents.narrative.narrators.integration import IntegrationNarrator
from aria.agents.narrative.render_blocks import render_blocks
from aria.agents.narrative.run_ledger import (
    build_run_ledger,
    ensure_report_ledger_nodes,
    node_id_for,
    verify_blocks_against_ledger,
)
from aria.agents.narrative.synthesis.cross_modal_triangulation import (
    build_cross_modal_evidence_graph,
    build_triangulation_block,
)
from aria.agents.narrative.validators import find_causal_language
from aria.agents.narrative_agent import NarrativeAgent
from aria.memory.memory import ARIAMemory


def _write_synthetic_wnn_h5mu(path: Path, n_cells: int = 60) -> Path:
    """Create a realistic HDF5 .h5mu with shuffled modality barcodes and categorical -1 NA codes."""
    rng = np.random.default_rng(42)
    barcodes = [f"cell_{i:03d}" for i in range(n_cells)]
    wnn_umap = rng.normal(0.0, 1.5, size=(n_cells, 2)).astype(np.float32)

    # Permute modality barcode order to verify explicit barcode alignment.
    rna_perm = rng.permutation(n_cells)
    atac_perm = rng.permutation(n_cells)
    rna_barcodes = [barcodes[i] for i in rna_perm]
    atac_barcodes = [barcodes[i] for i in atac_perm]
    rna_umap = (wnn_umap[rna_perm] * 0.9 + 0.2).astype(np.float32)
    atac_umap = (wnn_umap[atac_perm] * 1.1 - 0.3).astype(np.float32)

    categories = np.array(["Oligo", "Astro", "Microglia"], dtype="S16")
    codes = np.array([i % 3 for i in range(n_cells)], dtype=np.int8)
    # Inject AnnData NA categorical code -1; must NEVER decode to categories[-1] ("Microglia").
    codes[0] = -1
    codes[5] = -1

    leiden_wnn = np.array([str(i % 4) for i in range(n_cells)], dtype="S4")
    rna_w = np.clip(rng.normal(0.54, 0.12, size=n_cells), 0.05, 0.95)
    atac_w = 1.0 - rna_w

    str_dt = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as f:
        obs = f.create_group("obs")
        obs.attrs["_index"] = "_index"
        obs.create_dataset("_index", data=np.asarray(barcodes, dtype=object), dtype=str_dt)
        sub_grp = obs.create_group("subclass")
        sub_grp.create_dataset("categories", data=categories)
        sub_grp.create_dataset("codes", data=codes)
        obs.create_dataset("leiden_wnn", data=leiden_wnn)
        obs.create_dataset("rna_weight", data=rna_w)
        obs.create_dataset("atac_weight", data=atac_w)

        obsm = f.create_group("obsm")
        obsm.create_dataset("X_wnn_umap", data=wnn_umap)
        obsm.create_dataset("X_umap", data=wnn_umap)

        mod = f.create_group("mod")
        for mod_name, mod_bcs, mod_coords in (
            ("rna", rna_barcodes, rna_umap),
            ("atac", atac_barcodes, atac_umap),
        ):
            mg = mod.create_group(mod_name)
            m_obs = mg.create_group("obs")
            m_obs.attrs["_index"] = "_index"
            m_obs.create_dataset("_index", data=np.asarray(mod_bcs, dtype=object), dtype=str_dt)
            m_obsm = mg.create_group("obsm")
            m_obsm.create_dataset("X_umap", data=mod_coords)

    return path


def _synthetic_integration_result(tmp_path: Path, guardrail_status: str = "pass") -> dict:
    h5mu_path = _write_synthetic_wnn_h5mu(tmp_path / "wnn_integrated.h5mu")
    p2g_csv = tmp_path / "peak_gene_links.csv"
    with open(p2g_csv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["gene", "peak", "corr", "padj", "distance_bp", "region_class"],
        )
        writer.writeheader()
        writer.writerow({
            "gene": "SOX10",
            "peak": "chr22:38380000-38380500",
            "corr": "0.4821",
            "padj": "0.00012",
            "distance_bp": "1400",
            "region_class": "promoter_proximal",
        })
        writer.writerow({
            "gene": "MBP",
            "peak": "chr18:74690000-74690500",
            "corr": "0.3915",
            "padj": "0.0014",
            "distance_bp": "42500",
            "region_class": "distal_cis",
        })
        writer.writerow({
            "gene": "GFAP",
            "peak": "chr17:44910000-44910500",
            "corr": "-0.2610",
            "padj": "0.0120",
            "distance_bp": "18200",
            "region_class": "distal_cis",
        })

    ereg_csv = tmp_path / "eregulons.csv"
    ereg_csv.write_text("tf,sign_class,n_target_genes\nSOX10,pos_tf_pos_p2g,12\n", encoding="utf-8")
    trip_csv = tmp_path / "eregulon_triplets.csv"
    trip_csv.write_text("tf,peak,gene\nSOX10,chr22:38380000-38380500,MBP\n", encoding="utf-8")
    mofa_csv = tmp_path / "mofa_factor_loadings.csv"
    mofa_csv.write_text("factor,modality,feature,loading\n1,rna,SOX10,0.82\n", encoding="utf-8")

    return {
        "status": "done",
        "strategy": "multiome_wnn_mofa",
        "findings": {
            "multiome_gate": {"passed": True, "shared_barcodes": 2517},
            "paired_qc": {
                "status": "success",
                "n_cells_before": 3143,
                "n_cells_after": 2517,
                "pct_removed": 19.92,
                "n_genes": 38606,
                "n_peaks": 46200,
                "rna_qc_summary": {"min_genes": 200, "max_pct_mt": 20.0},
                "atac_qc_summary": {"min_peaks": 200},
                "label_join": {"n_matched": 2490},
                "metrics_not_computed": [
                    "fragment_length_periodicity",
                    "nucleosome_signal",
                ],
                "output_path": str(tmp_path / "paired_qc_filtered.h5mu"),
            },
            "wnn": {
                "status": "success",
                "n_cells": 2517,
                "k": 20,
                "resolution": 0.5,
                "n_joint_clusters": 11,
                "n_rna_only_clusters": 12,
                "n_atac_only_clusters": 9,
                "mean_rna_weight": 0.5312,
                "mean_atac_weight": 0.4688,
                "wnn_vs_rna_ari": 0.6142,
                "wnn_vs_atac_ari": 0.4025,
                "wnn_vs_rna_nmi": 0.7215,
                "wnn_vs_atac_nmi": 0.5380,
                "n_discordant_clusters": 1 if guardrail_status != "pass" else 0,
                "guardrail_status": guardrail_status,
                "per_cluster_modality_weights": {
                    "0": {"n_cells": 25, "mean_rna_weight": 0.58, "mean_atac_weight": 0.42, "dominant_modality": "rna"},
                    "1": {"n_cells": 15, "mean_rna_weight": 0.51, "mean_atac_weight": 0.49, "dominant_modality": "balanced"},
                    "2": {"n_cells": 10, "mean_rna_weight": 0.44, "mean_atac_weight": 0.56, "dominant_modality": "atac"},
                    "3": {"n_cells": 10, "mean_rna_weight": 0.60, "mean_atac_weight": 0.40, "dominant_modality": "rna"},
                },
                "anti_degradation": {
                    "guardrail_status": guardrail_status,
                    "axes": {
                        "modality_balance": {"status": "pass"},
                        "neighborhood_preservation": {"status": guardrail_status},
                    },
                    "issues": (
                        [{"code": "knn_jaccard_drop", "severity": guardrail_status, "detail": "synthetic check"}]
                        if guardrail_status != "pass"
                        else []
                    ),
                },
                "output_path": str(h5mu_path),
            },
            "mofa": {
                "status": "success",
                "n_cells": 2517,
                "n_cells_trained": 2517,
                "subsampled": False,
                "n_factors": 10,
                "n_active_factors": 9,
                "top_factors": [
                    {"factor": 1, "rna_r2": 0.0842, "atac_r2": 0.0415, "classification": "shared"},
                    {"factor": 2, "rna_r2": 0.0510, "atac_r2": 0.0040, "classification": "rna_specific"},
                    {"factor": 3, "rna_r2": 0.0035, "atac_r2": 0.0390, "classification": "atac_specific"},
                ],
                "shared_factors": [{"factor": 1}],
                "modality_specific_factors": [{"factor": 2}, {"factor": 3}],
                "excluded_factors": [
                    {
                        "factor": 5,
                        "rna_r2": 0.0215,
                        "atac_r2": 0.0680,
                        "exclusion_reasons": ["technical_depth_correlation_|r|>=0.70"],
                    }
                ],
                "variance_explained": {"rna": 0.2145, "atac": 0.1480},
                "variance_explained_additive": {"rna": 0.2310, "atac": 0.1620},
                "factor_loadings_path": str(mofa_csv),
            },
            "peak_to_gene": {
                "status": "success",
                "n_cells": 2517,
                "distance_kb": 150,
                "min_corr": 0.10,
                "max_padj": 0.05,
                "background_method": "trans_chromosome_gc_accessibility",
                "n_genes_tested": 1200,
                "n_peaks_tested": 8500,
                "n_tested": 14200,
                "n_links": 184,
                "n_positive_correlations": 162,
                "n_negative_correlations": 22,
                "n_promoter_proximal": 48,
                "n_distal_cis": 136,
                "top_links": [
                    {
                        "gene": "SOX10",
                        "peak": "chr22:38380000-38380500",
                        "corr": 0.4821,
                        "padj": 0.00012,
                        "distance_bp": 1400,
                        "region_class": "promoter_proximal",
                    },
                    {
                        "gene": "MBP",
                        "peak": "chr18:74690000-74690500",
                        "corr": 0.3915,
                        "padj": 0.0014,
                        "distance_bp": 42500,
                        "region_class": "distal_cis",
                    },
                    {
                        "gene": "GFAP",
                        "peak": "chr17:44910000-44910500",
                        "corr": -0.2610,
                        "padj": 0.0120,
                        "distance_bp": 18200,
                        "region_class": "distal_cis",
                    },
                ],
                "output_path": str(p2g_csv),
            },
            "eregulons": {
                "status": "success",
                "n_shared_cells": 2517,
                "scan_window_bp": 500,
                "motif_pvalue": 1e-4,
                "scanner_backend": "lightweight_pwm",
                "motif_source": "JASPAR2024_CORE_vertebrates",
                "n_candidate_modules_tested": 92,
                "n_eregulons": 8,
                "n_triplets": 64,
                "n_unique_tfs": 7,
                "n_unique_enhancers": 51,
                "n_unique_target_genes": 43,
                "n_pos_tf_pos_p2g": 6,
                "n_neg_tf_pos_p2g": 2,
                "n_paralog_ambiguous": 2,
                "n_depth_confounded": 0,
                "top_eregulons": [
                    {
                        "tf": "SOX10",
                        "motif_id": "MA0442.2",
                        "sign_class": "pos_tf_pos_p2g",
                        "n_target_genes": 12,
                        "n_enhancer_peaks": 15,
                        "tf_gene_corr": 0.4412,
                        "padj": 0.00021,
                        "rna_atac_aucell_corr": 0.3620,
                        "motif_family_ambiguity": False,
                        "depth_confounded": False,
                        "top_target_genes": ["MBP", "PLP1"],
                    },
                    {
                        "tf": "ETV6",
                        "motif_id": "MA0645.1",
                        "sign_class": "neg_tf_pos_p2g",
                        "n_target_genes": 7,
                        "n_enhancer_peaks": 8,
                        "tf_gene_corr": -0.3105,
                        "padj": 0.0038,
                        "rna_atac_aucell_corr": -0.2150,
                        "motif_family_ambiguity": True,
                        "depth_confounded": False,
                        "top_target_genes": ["OLIG1"],
                    },
                ],
                "output_path": str(ereg_csv),
                "triplets_path": str(trip_csv),
            },
        },
    }


def test_ledger_routing_prevents_peak2gene_alias_hijack_and_failed_status_promotion(tmp_path: Path):
    assert node_id_for("integration", "peak_to_gene") == "ledger://integration/peak_to_gene"
    assert node_id_for("integration", "peak2gene") == "ledger://integration/peak_to_gene"
    assert node_id_for("chromatin", "peak2gene") == "ledger://chromatin/regulatory_layers"
    assert (
        node_id_for("integrated synthesis", "integrated_signal")
        == "ledger://report/biological_synthesis"
    )

    agent_results = {
        "integration_agent": {
            "status": "partial",
            "findings": {
                "paired_qc": {"status": "success", "n_cells_before": 100, "n_cells_after": 90, "pct_removed": 10.0, "n_genes": 500, "n_peaks": 800},
                "wnn": {"status": "not_run", "reason": "skipped by test"},
                "mofa": {"status": "failed", "details": "solver diverged"},
                "peak_to_gene": {"status": "blocked", "reason": "missing gtf"},
            },
        }
    }
    exp_ctx = {"design_intelligence": {"recommended": ["Run MOFA factor analysis and peak2gene"]}}
    ledger = build_run_ledger(exp_ctx, agent_results)
    by_id = {e["node_id"]: e for e in ledger["entries"]}
    assert by_id["ledger://integration/paired_qc"]["status"] == "ran"
    assert by_id["ledger://integration/wnn"]["status"] == "not_run"
    assert by_id["ledger://integration/mofa"]["status"] == "error"
    assert by_id["ledger://integration/mofa"]["divergence"] is True
    assert by_id["ledger://integration/peak_to_gene"]["status"] == "skipped"
    assert by_id["ledger://integration/peak_to_gene"]["divergence"] is True


def test_cross_modal_evidence_graph_non_circular_and_portable(tmp_path: Path):
    integ = _synthetic_integration_result(tmp_path)
    agent_results = {
        "integration_agent": integ,
        "bulk_rna_agent": {
            "status": "done",
            "findings": {
                "contrasts": [{
                    "name": "diff_vs_prolif",
                    "status": "success",
                    "n_significant": 2,
                    "top_genes": [
                        {"symbol": "SOX10", "log2fc": 2.1},
                        {"symbol": "MBP", "log2fc": 3.4},
                    ],
                }]
            },
        },
    }
    graph = build_cross_modal_evidence_graph(
        agent_results, {"organism": "Homo sapiens", "genome": "GRCh38"}
    )
    assert graph["ran"] is True
    assert graph["non_causal_contract"] is True
    assert graph["circular_convergence_prevented"] is True

    # Verify zero absolute paths leaked into the serialized graph.
    serialized = json.dumps(graph)
    assert str(tmp_path) not in serialized
    assert "/home/" not in serialized

    # Verify eRegulon edges carry explicit shared_upstream=["peak_to_gene"].
    ereg_edges = [
        e for e in graph["edges"] if e["relation"] == "eregulon_tf_to_gene"
    ]
    assert ereg_edges
    assert all(e["shared_upstream"] == ["peak_to_gene"] for e in ereg_edges)
    assert all(e["non_causal"] is True for e in graph["edges"])

    # SOX10 corroborated by external bulk RNA DE -> convergent_cross_contrast;
    # ETV6 neg_tf_pos_p2g -> discordant_signal;
    # GFAP negative peak2gene corr -> discordant_signal;
    # Factor 5 depth confounder -> caveated_ambiguous;
    # ATAC DA was not run -> tracked in absent_not_tested, NOT discordant.
    axes_by_id = {a["axis_id"]: a for a in graph["regulatory_axes"]}
    assert axes_by_id["eregulon:SOX10:pos_tf_pos_p2g"]["classification"] == "convergent_cross_contrast"
    assert axes_by_id["eregulon:ETV6:neg_tf_pos_p2g"]["classification"] == "discordant_signal"
    assert axes_by_id["peak_to_gene:chr17:44910000-44910500->GFAP"]["classification"] == "discordant_signal"
    assert axes_by_id["mofa_confounder:Factor_5"]["classification"] == "caveated_ambiguous"
    assert "modality.atac_differential_accessibility" in graph["absent_not_tested"]


def test_integration_narrator_c1_c4_ledger_and_claim6_monotone_invariance(tmp_path: Path):
    integ = _synthetic_integration_result(tmp_path, guardrail_status="pass")
    fig_dir = tmp_path / "figures"
    _narrative_integration.generate_figures(integ["findings"], fig_dir)

    agent_results = {"integration_agent": integ}
    narrator = IntegrationNarrator()
    assert narrator.accepts("integration_agent", integ)

    blocks = narrator.collect(
        "integration_agent",
        integ,
        {"exp_context": {"organism": "Homo sapiens"}, "agent_results": agent_results},
    )
    block_ids = [b.id for b in blocks]
    assert block_ids == [
        "multiome.paired_qc",
        "multiome.wnn",
        "multiome.mofa",
        "multiome.peak_to_gene",
        "multiome.eregulons",
        "multiome.triangulation",
    ]

    # Zero causal patterns across claims, caveats, and rendered prose.
    for block in blocks:
        prose = compose_block_prose(block)
        assert find_causal_language(block.claim) is None
        assert find_causal_language(prose) is None
        for cav in block.caveats:
            assert find_causal_language(cav.text) is None

    # Full C1-C4 public compilation + C3 strict ledger verification.
    ledger = build_run_ledger({}, agent_results)
    ensure_report_ledger_nodes(ledger, blocks)
    compiled = compile_public_claims(blocks, exp_ctx={}, run_ledger=ledger)
    assert compiled.withheld == []
    assert len(compiled.claims) == 6
    # Non-circular tier cap: observational paired multiome claims never exceed associative.
    for claim in compiled.claims:
        assert claim["tier"] in {"descriptive", "associative"}
        assert claim["licensed_language"] in {"descriptive", "associative"}
        assert claim["language_violation"] is None
        assert claim["ledger_node_id"].startswith("ledger://integration/")

    ledger_check = verify_blocks_against_ledger(compiled.blocks, ledger, strict=True)
    assert ledger_check["status"] == "ok"
    assert ledger_check["n_violations"] == 0

    # Devil's advocate challenges attach cleanly without causal-lexicon collisions.
    da_manifest = build_devils_advocate(compiled.blocks, agent_results, {})
    assert len(da_manifest) >= 4
    for block in compiled.blocks:
        for cav in block.caveats:
            assert find_causal_language(cav.text) is None

    # Claim 6 monotone invariance: benign claim text variation does not alter rendered_prose;
    # overreaching causal or ungrounded entity variation is withheld.
    ereg_orig = next(b for b in blocks if b.id == "multiome.eregulons")
    prose_orig = compose_block_prose(ereg_orig)
    ereg_paraphrase = copy.deepcopy(ereg_orig)
    ereg_paraphrase.claim = (
        f"Across 2517 shared cells, 92 modules were tested and 8 BH-FDR "
        f"candidate eRegulons (6 pos_tf_pos_p2g, 2 neg_tf_pos_p2g) across "
        f"7 TFs and 64 triplets were retained (2 paralog-ambiguous, 0 depth-confounded)."
    )
    assert compose_block_prose(ereg_paraphrase) == prose_orig

    ereg_causal = copy.deepcopy(ereg_orig)
    ereg_causal.claim = "SOX10 directly regulates MBP across 2517 shared cells."
    causal_compiled = compile_public_claims([ereg_causal], exp_ctx={}, run_ledger=ledger)
    assert len(causal_compiled.claims) == 0
    assert len(causal_compiled.withheld) == 1

    # Single-family guard stays active on individual stage blocks.
    p2g_bad_family = copy.deepcopy(next(b for b in blocks if b.id == "multiome.peak_to_gene"))
    p2g_bad_family.claim += " Pathway enrichment GSEA found 184 terms."
    with pytest.raises(EvidenceVerificationError):
        verify_block_claim_support(p2g_bad_family, strict=True)


def test_wnn_degraded_guardrail_and_skipped_stages_remain_honest_and_visible(tmp_path: Path):
    integ = _synthetic_integration_result(tmp_path, guardrail_status="degraded")
    integ["findings"]["eregulons"] = {
        "status": "skipped",
        "reason": "genome_fasta not supplied",
    }
    agent_results = {"integration_agent": integ}
    blocks = IntegrationNarrator().collect(
        "integration_agent", integ, {"agent_results": agent_results}
    )
    wnn_block = next(b for b in blocks if b.id == "multiome.wnn")
    ereg_block = next(b for b in blocks if b.id == "multiome.eregulons")
    assert wnn_block.confidence == "low"
    assert any("LOW/INSUFFICIENT" in w for w in NarrativeAgent._narrative_registry(None).collect_blocks(agent_results)[1].warnings)
    assert ereg_block.status == "skipped"
    assert ereg_block.claim == ""

    ledger = build_run_ledger({}, agent_results)
    compiled = compile_public_claims(blocks, exp_ctx={}, run_ledger=ledger)
    assert compiled.withheld == []
    assert verify_blocks_against_ledger(compiled.blocks, ledger, strict=True)["n_violations"] == 0
    rendered_html = render_blocks(compiled.blocks, report_dir=tmp_path, strict=True)
    assert "genome_fasta not supplied" in rendered_html


def test_dual_png_svg_figures_and_h5py_categorical_minus_one_safety(tmp_path: Path):
    integ = _synthetic_integration_result(tmp_path)
    bundle = _narrative_integration._read_wnn_h5mu(integ["findings"]["wnn"]["output_path"])
    assert bundle is not None
    # Verify code == -1 decoded to "NA" and NEVER categories[-1] ("Microglia").
    subclass = bundle["obs"]["subclass"]
    assert subclass[0] == "NA"
    assert subclass[5] == "NA"
    assert subclass[1] == "Astro"
    # Verify shuffled RNA/ATAC UMAPs were barcode-aligned back to /obs/_index.
    assert np.allclose(bundle["rna_umap"], bundle["wnn_umap"] * 0.9 + 0.2, atol=1e-5)
    assert np.allclose(bundle["atac_umap"], bundle["wnn_umap"] * 1.1 - 0.3, atol=1e-5)

    fig_dir = tmp_path / "figures"
    _narrative_integration.generate_figures(integ["findings"], fig_dir)
    expected = (
        "wnn_joint_vs_modality_umap",
        "wnn_modality_weights_by_cluster",
        "mofa_variance_explained",
        "regulatory_triangulation_summary",
    )
    for key in expected:
        png_path = Path(integ["findings"]["figures"][key])
        svg_path = Path(integ["findings"]["figures_svg"][key])
        assert png_path.exists() and png_path.stat().st_size > 1000
        assert svg_path.exists() and svg_path.stat().st_size > 1000
        svg_text = svg_path.read_text(encoding="utf-8")
        assert "<svg" in svg_text and "<text" in svg_text


def test_methods_html_safety_and_full_html_report_e2e(tmp_path: Path):
    integ = _synthetic_integration_result(tmp_path)
    # Inject adversarial HTML tag into motif_source to verify tag stripping + single escaping.
    integ["findings"]["eregulons"]["motif_source"] = "JASPAR2024<script>alert(1)</script>"
    # Mark MOFA skipped to verify methods() omits unexecuted stages.
    integ["findings"]["mofa"] = {"status": "skipped", "reason": "disabled"}

    paragraphs = IntegrationNarrator().methods("integration_agent", integ)
    joined = "\n".join(paragraphs)
    assert "<script>" not in joined
    assert "MOFA+ latent factor decomposition" not in joined
    assert "Paired multiome quality control" in joined
    assert "Candidate eRegulon inference" in joined

    mem = ARIAMemory(db_path=":memory:")
    agent = NarrativeAgent(memory=mem, llm=None)
    agent.reports_dir = tmp_path / "reports"
    agent.reports_dir.mkdir(parents=True, exist_ok=True)

    res = agent.run(
        experiment_id="exp_v48_i7_e2e",
        context={
            "exp_context": {
                "organism": "Homo sapiens",
                "genome": "GRCh38",
                "user_question": "Characterize paired RNA+ATAC regulatory programs",
                "reproducible_mode": True,
            },
            "biological_intent": {"summary": "Paired multiome regulatory analysis"},
            "agent_results": {"integration_agent": integ},
            "findings": [],
        },
    )
    assert res["status"] == "done"
    report_html_path = Path(res["report_path"])
    assert report_html_path.exists()
    html_text = report_html_path.read_text(encoding="utf-8")
    assert "Multi-Omics Integration (Paired RNA + ATAC)" in html_text
    assert "<script>alert(1)</script>" not in html_text
    assert "&amp;lt;" not in html_text

    meth_json = json.loads((report_html_path.parent / "methodology.json").read_text(encoding="utf-8"))
    assert meth_json["cross_modal_evidence_graph"]["ran"] is True
    assert meth_json["run_ledger"]["claim_ledger_verification"]["n_violations"] == 0
    assert (report_html_path.parent / "ro-crate-metadata.json").exists()


def test_integration_validation_promoted_to_beta_in_slice_i8():
    from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION

    assert INTEGRATION_VALIDATION["level"] == "production"
    assert INTEGRATION_VALIDATION["dispatch_enabled"] is True
    assert INTEGRATION_VALIDATION["requires_ack"] is False
    assert INTEGRATION_VALIDATION["production_promoted"] is True
    assert INTEGRATION_VALIDATION["scope"] == "paired_same_cell_rna_atac"


def test_real_ipc_schema_mofa_dict_p2g_correlation_and_json_loadings(tmp_path: Path):
    """Verify exact real I4/I5/I6 IPC payload schemas stage CSVs, figures, and claims cleanly."""
    integ = _synthetic_integration_result(tmp_path)
    mofa_json_path = tmp_path / "mofa_loadings.json"
    mofa_json_path.write_text(
        json.dumps(
            {
                "seed": 0,
                "n_factors": 3,
                "views": ["rna", "atac"],
                "factors": {
                    "Factor1": {
                        "rna": [
                            {
                                "feature": "MBP",
                                "weight": 0.812,
                                "abs_weight": 0.812,
                                "direction": "positive",
                            }
                        ],
                        "atac": [
                            {
                                "feature": "chr18:76980000-76980500",
                                "weight": 0.541,
                                "abs_weight": 0.541,
                                "direction": "positive",
                            }
                        ],
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    f1 = {
        "factor_id": 1,
        "active": True,
        "category": "shared",
        "r2_by_view": {"rna": 0.1515, "atac": 0.0110},
        "excluded_from_biological_interpretation": False,
        "exclusion_reasons": [],
    }
    f2 = {
        "factor_id": 2,
        "active": True,
        "category": "modality_specific",
        "r2_by_view": {"rna": 0.1105, "atac": 0.0006},
        "excluded_from_biological_interpretation": False,
        "exclusion_reasons": [],
    }
    f3_ex = {
        "factor_id": 3,
        "active": True,
        "category": "modality_specific",
        "r2_by_view": {"rna": 0.0019, "atac": 0.0164},
        "excluded_from_biological_interpretation": True,
        "exclusion_reasons": ["sequencing_depth_correlated:atac(|r|=0.941)"],
    }
    integ["findings"]["mofa"] = {
        "status": "success",
        "n_cells": 2517,
        "n_cells_trained": 2517,
        "n_factors": 3,
        "n_active_factors": 2,
        "subsampled": False,
        "variance_explained": {"rna": 0.2620, "atac": 0.0280},
        "variance_explained_additive": {"rna": 0.2639, "atac": 0.0280},
        "top_factors": [f1, f2, f3_ex],
        "shared_factors": [f1],
        "modality_specific_factors": {"rna": [f2], "atac": [f3_ex]},
        "excluded_factors": [f3_ex],
        "factor_loadings_path": str(mofa_json_path),
    }
    integ["findings"]["peak_to_gene"]["top_links"] = [
        {
            "gene": "PDE4B",
            "peak": "chr1:66258000-66258500",
            "correlation": 0.3412,
            "padj": 0.0004,
            "distance_bp": 14200,
            "region_class": "distal_cis",
        }
    ]
    integ["findings"]["eregulons"]["motif_source"] = {
        "collection": "JASPAR2024_CORE_vertebrates",
        "release": "2024",
        "n_motifs": 879,
    }
    integ["findings"]["eregulons"]["top_eregulons"] = [
        {
            "tf": "TCF12",
            "best_motif_id": "MA1648.2",
            "sign_class": "pos_tf_pos_p2g",
            "n_target_genes": 18,
            "n_enhancer_peaks": 21,
            "mean_r_tf_gene": 0.3821,
            "regulon_padj": 0.00012,
            "rna_atac_activity_corr": 0.2940,
            "motif_family_ambiguity": True,
            "depth_confounded": False,
            "target_genes": ["PDE4B", "MBP"],
            "enhancer_peaks": ["chr1:66258000-66258500"],
        }
    ]

    fig_dir = tmp_path / "figs_real_schema"
    tables_dir = tmp_path / "tables_real_schema"
    _narrative_integration.generate_figures(integ["findings"], fig_dir)
    staged = _narrative_integration.stage_integration_tables(integ["findings"], tables_dir)
    assert staged["mofa_factor_loadings"] == "tables/mofa_factor_loadings.csv"
    csv_lines = (tables_dir / "mofa_factor_loadings.csv").read_text(encoding="utf-8").splitlines()
    assert csv_lines[0] == "factor,modality,feature,weight,abs_weight,direction"
    assert any("MBP" in ln for ln in csv_lines[1:])

    blocks = IntegrationNarrator().collect(
        "integration_agent",
        integ,
        {"exp_context": {"organism": "Homo sapiens"}, "agent_results": {"integration_agent": integ}},
    )
    by_id = {b.id: b for b in blocks}
    assert by_id["multiome.mofa"].metrics["n_shared_factors"] == 1
    assert by_id["multiome.mofa"].metrics["n_modality_specific_factors"] == 1
    assert by_id["multiome.mofa"].metrics["n_excluded_factors"] == 1
    assert "PDE4B" in by_id["multiome.peak_to_gene"].metrics["top_genes"]
    assert "TCF12" in by_id["multiome.eregulons"].metrics["top_tfs"]
    methods_txt = "\n".join(IntegrationNarrator().methods("integration_agent", integ))
    assert "JASPAR2024_CORE_vertebrates" in methods_txt
    assert "{'collection'" not in methods_txt
