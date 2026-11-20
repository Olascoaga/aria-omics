"""Wave 4A regression tests: narrator contracts & integration hypothesis adapter.

Each test pins one audit finding (bulk RNA, scRNA, chromatin, integration,
Hi-C narrators + the integration hypothesis adapter) so the contract/key
mismatches and NaN/inf crashes cannot silently return.
"""
from __future__ import annotations

import math

import pytest


# ═══════════════════════════════════════════════════════════════════════════
# Bulk RNA narrator
# ═══════════════════════════════════════════════════════════════════════════

def _bulk_result(contrast: dict) -> dict:
    return {"status": "done", "findings": {"contrasts": [contrast]}}


def test_bulk_ora_evidence_reads_padj():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    contrast = {
        "name": "treated_vs_ctrl", "status": "success",
        "pathways": {"GO_BP": [{"term": "immune response", "padj": 1.2e-6,
                                "pvalue": 1e-8}]},
    }
    blocks = BulkRnaNarrator()._pathway_blocks_for_contrast(contrast)
    labels = {ev.label: ev.value for ev in blocks[0].evidence}
    assert labels["immune response FDR"] == pytest.approx(1.2e-6)


def test_bulk_ora_evidence_legacy_adjusted_p_alias():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    contrast = {"name": "c", "pathways": {"KEGG": [
        {"Term": "glycolysis", "adjusted_p": 0.003}]}}
    blocks = BulkRnaNarrator()._pathway_blocks_for_contrast(contrast)
    labels = {ev.label: ev.value for ev in blocks[0].evidence}
    assert labels["glycolysis FDR"] == pytest.approx(0.003)


def test_bulk_gsea_fdr_keeps_full_precision(tmp_path):
    from aria.agents.narrative.narrators.bulk_rna import (
        _parse_gsea_row, _read_gsea_summary,
    )

    row = _parse_gsea_row({"Term": "T1", "NES": "2.123456", "FDR q-val": "1e-05"})
    assert row["fdr"] == pytest.approx(1e-5)
    assert row["fdr"] != 0.0
    assert row["nes"] == pytest.approx(2.1235)

    csv = tmp_path / "gsea.csv"
    csv.write_text(
        "Term,NES,FDR q-val\n"
        "weaker,3.0,4e-5\n"
        "stronger,1.5,1e-15\n",
        encoding="utf-8",
    )
    summary = _read_gsea_summary(str(csv))
    # 1e-15 must sort ahead of 4e-5 (rounding both to 0.0 broke this).
    assert summary["top_pathways"][0]["term"] == "stronger"
    assert summary["top_pathways"][0]["fdr"] == pytest.approx(1e-15)


def test_bulk_figures_only_images_and_running_sums():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    contrast = {
        "name": "c1", "status": "success",
        "plots": {
            "volcano": "/r/c1/volcano.png",
            "heatmap": "/r/c1/heatmap_padj.png",
            "heatmap_padj": "/r/c1/heatmap_padj.png",
            "gsea_top_table": "/r/c1/gsea_top.png",
            "gsea_table": "/r/c1/gsea.csv",
            "error": "Pathway visualization failed: boom",
            "tables": {"de_genes": "/r/c1/de_genes.tsv"},
            "ora_dotplots": {"GO_BP": "/r/c1/ora_go.png"},
            "gsea_running_sums": ["/r/c1/rs1.png", "/r/c1/rs2.svg"],
        },
    }
    figs = BulkRnaNarrator().figures("rna_agent", _bulk_result(contrast))
    paths = [f["path"] for f in figs]
    assert "/r/c1/rs1.png" in paths and "/r/c1/rs2.svg" in paths
    assert "/r/c1/ora_go.png" in paths
    assert "/r/c1/gsea_top.png" in paths
    assert not any(p.endswith((".csv", ".tsv")) for p in paths)
    assert not any("failed" in p for p in paths)
    assert paths.count("/r/c1/heatmap_padj.png") == 1


def test_bulk_tables_include_gsea_and_full_results():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    contrast = {
        "name": "c1", "status": "success",
        "full_results_csv": "/r/c1/full_de.tsv",
        "plots": {"tables": {"de_genes": "/r/c1/de.tsv"},
                  "gsea_table": "/r/c1/gsea.csv"},
    }
    tables = BulkRnaNarrator().tables("rna_agent", _bulk_result(contrast))
    paths = {t["path"] for t in tables}
    assert paths == {"/r/c1/de.tsv", "/r/c1/gsea.csv", "/r/c1/full_de.tsv"}


def test_bulk_failed_contrast_reports_error_key():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator

    findings = {"contrasts": [{"name": "c1", "status": "error",
                               "error": "DESeq2 failed: singular design"}]}
    blocks = BulkRnaNarrator()._contrast_blocks(findings)
    assert blocks[0].error == "DESeq2 failed: singular design"

    findings = {"contrasts": [{"name": "c2", "status": "skipped",
                               "reason": "n<2 replicates"}]}
    blocks = BulkRnaNarrator()._contrast_blocks(findings)
    assert blocks[0].error == "n<2 replicates"


# ═══════════════════════════════════════════════════════════════════════════
# scRNA tables / text / _common
# ═══════════════════════════════════════════════════════════════════════════

def _per_cluster_comp() -> dict:
    """Local BH found 120 genes; global BH found 0 (rna_pseudobulk_de default)."""
    return {
        "status": "success",
        "fdr_strategy": "per_cluster",
        "n_significant": 120, "n_significant_local": 120,
        "n_significant_global": 0,
        "n_up": 70, "n_up_local": 70, "n_up_global": 0,
        "n_down": 50, "n_down_local": 50, "n_down_global": 0,
        "top_genes": [
            {"gene": "GeneA", "log2fc": 2.0, "padj": 1e-4,
             "padj_local": 1e-4, "padj_global": 0.42},
            {"gene": "GeneZ", "log2fc": -1.5, "padj": 2e-3,
             "padj_local": 2e-3, "padj_global": 0.6},
        ],
    }


def _pb_findings(comp: dict | None = None) -> dict:
    return {"pseudobulk_de": {
        "groupby": "cell_type", "n_groups": 1,
        "multiple_testing": {"fdr_strategy": "per_cluster"},
        "per_group": {"Tcell": {"n_pseudosamples": 8, "per_comparison": {
            "trt_vs_ctrl": comp or _per_cluster_comp()}}},
    }}


def test_scrna_pseudobulk_table_uses_primary_fdr_family():
    from aria.agents.narrative.scrna.tables import extract_pseudobulk_de_table

    html_rows = extract_pseudobulk_de_table(_pb_findings())
    assert "<strong>120</strong>" in html_rows
    assert "70 ↑" in html_rows and "50 ↓" in html_rows
    assert ">0 ↑" not in html_rows and ">0 ↓" not in html_rows
    assert "per-cluster" in html_rows


def test_scrna_prose_uses_primary_fdr_family():
    from aria.agents.narrative.scrna.text import _describe_top_de_blocks

    lines = _describe_top_de_blocks(_pb_findings())
    txt = " ".join(lines)
    assert "120 per-cluster FDR DE genes" in txt
    assert "70 up" in txt and "50 down" in txt
    assert "0 global-FDR DE genes" not in txt
    # top gene FDR is the primary (per-cluster) padj, not padj_global=0.42
    assert "GeneA (log2FC=+2.00, FDR=1.00e-04)" in txt


def test_scrna_legacy_global_fields_still_supported():
    from aria.agents.narrative.scrna._common import _primary_de_counts

    legacy = {"n_significant_global": 12, "n_significant_local": 20,
              "n_up_global": 8, "n_down_global": 4}
    c = _primary_de_counts(legacy)
    assert (c["n_significant"], c["n_up"], c["n_down"]) == (12, 8, 4)
    assert c["label"] == "global-FDR"


def test_scrna_dpt_nonfinite_does_not_crash():
    from aria.agents.narrative.scrna.tables import (
        export_supplementary_tables, extract_trajectory_tables,
    )

    findings = {"trajectory": {"pseudotime": {"pseudotime_by_group": {
        "A": 0.1, "B": 0.5, "C": float("inf"), "D": float("nan")}}}}
    out = extract_trajectory_tables(findings)
    rows = out["pseudotime_rows"]
    assert "disconnected" in rows
    assert rows.index("<strong>A</strong>") < rows.index("<strong>B</strong>")


def test_scrna_paga_nan_connectivity_and_zero_threshold(tmp_path):
    from aria.agents.narrative.scrna.tables import (
        export_supplementary_tables, extract_trajectory_tables,
    )

    findings = {"trajectory": {"paga": {
        "top_connections": {"0->1": 0.01, "1->2": float("nan")},
        "strong_threshold": 0.0,
    }}}
    rows = extract_trajectory_tables(findings)["paga_rows"]
    # explicit 0.0 threshold honoured: 0.01 >= 0.0 is strong
    assert rows.split("</tr>")[0].count("strong") == 1
    tables = export_supplementary_tables(findings, tmp_path)
    text = (tmp_path / "scrna_paga_connections.tsv").read_text()
    assert "True" in text
    assert tables["paga_connections"]


def test_scrna_cellphone_pval_zero_is_rendered():
    from aria.agents.narrative.scrna.tables import extract_cellcomm_table

    findings = {"cell_communication": {"top_interactions": [
        {"source": "A", "target": "B", "ligand": "L", "receptor": "R",
         "score": 0.1, "rank": 1, "cellphone_pval": 0.0},
        {"source": "B", "target": "A", "ligand": "L2", "receptor": "R2",
         "score": 0.2, "rank": 2, "cellphone_pval": float("nan")},
    ]}}
    rows = extract_cellcomm_table(findings).split("</tr>")
    assert "<td><code>0</code></td>" in rows[0]
    assert rows[1].rstrip().endswith("<td>—</td>")


def test_scrna_pathway_export_reads_padj_and_pvalue(tmp_path):
    from aria.agents.narrative.scrna.tables import export_supplementary_tables

    findings = {"pseudobulk_pathways": {"per_cluster": {"T::c": {"results": {
        "GO_BP": [{"term": "x", "padj": 0.01, "pvalue": 0.0001}]}}}}}
    export_supplementary_tables(findings, tmp_path)
    lines = (tmp_path / "scrna_pathway_enrichment.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["adjusted_p"] == "0.01"
    assert row["p_value"] == "0.0001"


def test_scrna_ora_method_wording_local_vs_enrichr():
    from aria.agents.narrative.scrna.text import (
        build_scrna_methods, summarize_scrna_text,
    )

    pwp = {"per_cluster": {"T::c": {"n_significant": 1, "results": {}}},
           "ora_method": "local_hypergeometric", "background_size": 1000}
    f = {"pseudobulk_pathways": pwp}
    assert "Enrichr" not in build_scrna_methods(f)
    assert "Enrichr" not in summarize_scrna_text(f)
    assert "local hypergeometric" in build_scrna_methods(f)

    pwp_legacy = dict(pwp, ora_method="enrichr")
    assert "Enrichr" in build_scrna_methods({"pseudobulk_pathways": pwp_legacy})


def test_scrna_composition_clause_respects_any_significant():
    from aria.agents.narrative.scrna.text import summarize_scrna_text

    comp = _per_cluster_comp()
    comp["corrected_for_composition"] = False
    comp["composition_skipped_reason"] = "collinear with condition"
    f = _pb_findings(comp)
    f["differential_abundance"] = {"any_significant": True, "per_comparison": {}}
    txt = summarize_scrna_text(f)
    assert "found no significant shifts" not in txt

    f["differential_abundance"] = {"any_significant": False, "per_comparison": {}}
    assert "found no significant shifts" in summarize_scrna_text(f)


def test_scrna_methods_mt_threshold_used_and_plain_text_operators():
    from aria.agents.narrative.scrna.text import build_scrna_methods

    f = {"qc": {"n_cells_before": 10, "n_cells_after": 9,
                "mt_threshold_used": 12.5}}
    f.update(_pb_findings())
    f["pseudobulk_de"]["thresholds"] = {"padj_max": 0.05, "lfc_min": 0.5}
    f["pseudobulk_pathways"] = {"per_cluster": {"x": {}},
                                "ora_method": "local_hypergeometric"}
    txt = build_scrna_methods(f)
    assert "mt-fraction cap at 12.5%" in txt
    assert "&lt;" not in txt and "&gt;" not in txt
    assert "< 0.05" in txt


def test_scrna_find_pathway_block_no_substring_collision():
    from aria.agents.narrative.scrna._common import _find_pathway_block

    pwp = {"per_cluster": {"CD8_T::trt_vs_ctrl": {"n_significant": 3},
                           "10::trt_vs_ctrl": {"n_significant": 2}}}
    assert _find_pathway_block(pwp, "T", "trt_vs_ctrl") == ("", {})
    assert _find_pathway_block(pwp, "1", "trt_vs_ctrl") == ("", {})
    key, _ = _find_pathway_block(pwp, "CD8_T", "trt_vs_ctrl")
    assert key == "CD8_T::trt_vs_ctrl"


# ═══════════════════════════════════════════════════════════════════════════
# scRNA narrator
# ═══════════════════════════════════════════════════════════════════════════

def test_scrna_narrator_dotplot_exact_block_id():
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock

    blocks = [
        NarrativeBlock(id=f"scrna.pathway.{k}", modality="scRNA-seq",
                       analysis="pseudobulk_pathways", block_type="result",
                       title=k, status="success", confidence="medium",
                       claim=f"ORA found 1 enriched term(s) for {k}.",
                       evidence=[EvidenceItem("enriched terms", 1,
                                              "pseudobulk_pathways")])
        for k in ("1", "10", "11")
    ]
    findings = {"figures": {"pathway_dotplots": {
        "1": ["/f/c1.png"], "10": ["/f/c10.png"], "11": ["/f/c11.png"]}}}
    ScrnaNarrator()._attach_artifacts(blocks, findings)
    by_id = {b.id: [f["path"] for f in b.figures] for b in blocks}
    assert by_id["scrna.pathway.1"] == ["/f/c1.png"]
    assert by_id["scrna.pathway.10"] == ["/f/c10.png"]
    assert by_id["scrna.pathway.11"] == ["/f/c11.png"]


def test_scrna_narrator_keeps_null_pseudobulk_blocks():
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator

    comp = {"status": "success", "fdr_strategy": "per_cluster",
            "n_significant": 0, "n_significant_global": 0,
            "n_significant_local": 0, "power_estimate_at_lfc_min": 0.4,
            "top_genes": []}
    blocks = ScrnaNarrator()._pseudobulk_blocks(_pb_findings(comp))
    assert len(blocks) == 1
    b = blocks[0]
    assert b.status == "success"
    assert b.metrics["n_significant"] == 0
    assert b.metrics["null_result"] is True
    assert b.metrics["power_estimate_at_lfc_min"] == 0.4
    assert "0 per-cluster FDR DE genes" in b.claim
    assert any("absence of significance" in c.text for c in b.caveats)


def test_scrna_narrator_pathway_padj_evidence():
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator

    f = {"pseudobulk_pathways": {"per_cluster": {"T::c": {
        "n_significant": 1,
        "results": {"GO_BP": [{"term": "t cell activation", "padj": 3e-7}]},
    }}}}
    blocks = ScrnaNarrator()._pathway_blocks(f)
    labels = {ev.label: ev.value for ev in blocks[0].evidence}
    assert labels["t cell activation FDR"] == pytest.approx(3e-7)


def test_scrna_annotation_confidence_ignores_nan():
    from aria.agents.narrative.narrators.scrna import (
        _annotation_confidence_by_label,
    )

    out = _annotation_confidence_by_label({
        "0": {"label": "B", "mean_confidence": 0.3, "n_cells": 100},
        "1": {"label": "B", "mean_confidence": float("nan"), "n_cells": 100},
    })
    assert out["B"]["mean_confidence"] == pytest.approx(0.3)
    assert out["B"]["uncertain"] is True


# ═══════════════════════════════════════════════════════════════════════════
# Chromatin narrator
# ═══════════════════════════════════════════════════════════════════════════

def test_chromatin_qc_reads_jsd_and_nested_fragment_fractions():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    qc = {"status": "success", "data_type": "CUT_AND_RUN", "n_samples": 2,
          "mito_fraction": None, "jsd_vs_control": 0.31,
          "fragment_sizes": {"subnucleosomal_fraction": 0.42,
                             "mononucleosomal_fraction": 0.2,
                             "cutrun_mononucleosomal_fraction": 0.33}}
    block = ChromatinNarrator()._qc_block(qc)
    labels = {ev.label: ev.value for ev in block.evidence}
    assert labels["IP vs control JSD"] == 0.31
    assert labels["Subnucleosomal fraction (<120 bp)"] == 0.42
    assert labels["Mononucleosomal fraction (180-247 bp)"] == 0.2
    assert labels["CUT&RUN mononucleosomal fraction (140-250 bp)"] == 0.33
    assert "None" not in block.claim


def test_chromatin_scatac_qc_claim_uses_fragments_not_unknown_peaks():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    qc = {"status": "success", "data_type": "scATAC", "n_cells": 900,
          "n_fragments": 123456, "mito_fraction": 0.02}
    block = ChromatinNarrator()._qc_block(qc)
    assert "? peaks" not in block.claim
    assert "123456 fragments" in block.claim


def test_chromatin_idr_reads_n_reproducible_peaks():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    peaks = {"status": "success", "data_type": "ChIP", "n_peaks": 5000,
             "genome": "hg38", "peak_reproducibility": {
                 "strategy": "idr",
                 "idr": {"ran": True, "n_reproducible_peaks": 4321}}}
    block = ChromatinNarrator()._peaks_block(peaks)
    labels = {ev.label: ev.value for ev in block.evidence}
    assert labels["IDR reproducible peaks"] == 4321

    lines = ChromatinNarrator().methods(
        "chromatin_agent", {"findings": {"peaks": peaks}})
    txt = " ".join(lines)
    assert "retaining 4321 reproducible peaks" in txt
    assert "? reproducible" not in txt


def test_chromatin_scatac_methods_qc_sentence_honours_computed_metrics():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    computed = {"findings": {"qc": {"status": "success", "data_type": "scATAC",
                                    "tss_enrichment": 7.1, "frip": 0.4,
                                    "qc_complete": True}}}
    txt = " ".join(ChromatinNarrator().methods("chromatin_agent", computed))
    assert "remain uncomputed" not in txt
    assert "TSS enrichment" in txt and "FRiP" in txt

    missing = {"findings": {"qc": {"status": "success", "data_type": "scATAC",
                                   "tss_enrichment": None, "frip": None}}}
    txt = " ".join(ChromatinNarrator().methods("chromatin_agent", missing))
    assert "TSS enrichment and FRiP remain uncomputed" in txt


def test_chromatin_footprint_heterodimer_figure_id_sanitised():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    fp = {"ran": True, "group_a": "A", "group_b": "B",
          "differential_summary": {"parsed": True},
          "aggregate_plots": {"FOS::JUN": {"png": "/f/fosjun.png"}}}
    block = ChromatinNarrator()._footprint_block(fp)
    ids = [f["id"] for f in block.figures]
    assert "footprint_FOS__JUN" in ids
    assert not any("::" in i for i in ids)


def test_chromatin_bulk_diffacc_normalization_mode_reported():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    da = {"status": "success", "ran": True, "data_type": "bulk_ATAC",
          "normalization_mode": "deseq2_median_of_ratios",
          "n_comparisons_success": 1, "n_sig_total": 10}
    block = ChromatinNarrator()._bulk_atac_diffacc_block(da)
    labels = {ev.label: ev.value for ev in block.evidence}
    assert labels["Normalization mode"] == "deseq2_median_of_ratios"
    assert any("median-of-ratios" in c.text for c in block.caveats)

    da["normalization_mode"] = "heterologous_spikein"
    block = ChromatinNarrator()._bulk_atac_diffacc_block(da)
    assert any("spike-in sizeFactors" in c.text for c in block.caveats)


# ═══════════════════════════════════════════════════════════════════════════
# Integration narrator
# ═══════════════════════════════════════════════════════════════════════════

def _integ(findings: dict) -> dict:
    return {"status": "done", "findings": findings}


def test_integration_gate_blocked_via_allowed_and_status():
    from aria.agents.narrative.narrators.integration import IntegrationNarrator

    for gate in ({"allowed": False, "status": "blocked", "reason": "unpaired"},
                 {"status": "blocked"},
                 {"passed": False}):
        blocks = IntegrationNarrator().collect(
            "integration_agent", _integ({"multiome_gate": gate}))
        assert any(b.id == "multiome.gate" for b in blocks), gate

    blocks = IntegrationNarrator().collect(
        "integration_agent",
        _integ({"multiome_gate": {"allowed": True, "status": "allowed"}}))
    assert not any(b.id == "multiome.gate" for b in blocks)


def test_integration_batch_qc_reads_check_key():
    from aria.agents.narrative.narrators.integration import _batch_integration

    wnn = {"batch_integration": {"applied": True, "quality_assessment": {
        "issues": [{"severity": "warning", "check": "possible_overcorrection",
                    "message": "m"},
                   {"severity": "warning", "code": "residual_batch_effect"}]}}}
    out = _batch_integration(wnn)
    assert "possible_overcorrection" in out["qc_warnings"]
    assert "residual_batch_effect" in out["qc_warnings"]


def _mofa_payload() -> dict:
    return {
        "status": "success", "n_cells": 1000, "n_factors": 5,
        "n_active_factors": 4,
        "views": {"scRNA": {"kind": "rna"}, "scATAC": {"kind": "atac"}},
        "variance_explained": {"scRNA": 0.42, "scATAC": 0.23},
        "variance_explained_additive": {"scRNA": 0.5, "scATAC": 0.3},
        "top_factors": [
            {"factor_id": 1, "active": True, "category": "shared"},
            {"factor_id": 2, "active": True, "category": "shared",
             "batch_flag": True},
            {"factor_id": 3, "active": True, "category": "shared",
             "technical_depth_flag": True,
             "excluded_from_biological_interpretation": True},
            {"factor_id": 4, "active": True, "category": "modality_specific"},
            {"factor_id": 5, "active": False, "category": "shared"},
        ],
        "shared_factors": [
            {"factor_id": 1}, {"factor_id": 2},
            {"factor_id": 3, "excluded_from_biological_interpretation": True},
        ],
        "modality_specific_factors": {"scRNA": [{"factor_id": 4}]},
        "excluded_factors": [{"factor_id": 3}],
        "batch_factor_ids": [2],
    }


def test_integration_mofa_view_keys_case_insensitive_and_bio_count():
    from aria.agents.narrative.narrators.integration import IntegrationNarrator

    block = IntegrationNarrator()._mofa_block(_mofa_payload(), {})
    assert "42.0% in RNA and 23.0% in ATAC" in block.claim
    assert "0.0%" not in block.claim
    # factors 1 and 4 are the only active non-quarantined factors
    assert block.metrics["n_active_biological_factors"] == 2
    assert "retained 2 active biological factors" in block.claim
    assert block.metrics["rna_r2_pct"] == 42.0

    mofa = _mofa_payload()
    mofa.pop("views")
    mofa["variance_explained"] = {"RNA": 0.1, "ATAC": 0.2}
    block = IntegrationNarrator()._mofa_block(mofa, {})
    assert "10.0% in RNA and 20.0% in ATAC" in block.claim


def test_integration_methods_mofa_bio_count_and_qc_thresholds():
    from aria.agents.narrative.narrators.integration import IntegrationNarrator

    findings = {
        "mofa": _mofa_payload(),
        "paired_qc": {"status": "success", "n_cells_before": 100,
                      "n_cells_after": 80, "n_genes": 5, "n_peaks": 6,
                      "rna_qc_summary": {}, "atac_qc_summary": {},
                      "thresholds_used": {"min_genes": 500, "max_pct_mt": None,
                                          "min_peaks": 1000}},
    }
    txt = " ".join(IntegrationNarrator().methods(
        "integration_agent", _integ(findings)))
    assert "min_genes >= 500" in txt and "min_peaks >= 1000" in txt
    assert "mitochondrial fraction <=" not in txt
    assert "min_genes >= 200" not in txt
    assert "yielding 2 active biological factors" in txt

    findings["paired_qc"]["thresholds_used"]["max_pct_mt"] = 15.0
    txt = " ".join(IntegrationNarrator().methods(
        "integration_agent", _integ(findings)))
    assert "mitochondrial fraction <= 15.0%" in txt


def test_integration_methods_peak2gene_dynamic_wording():
    from aria.agents.narrative.narrators.integration import IntegrationNarrator

    p2g = {"status": "success", "atac_normalization": "raw_counts",
           "background_match": "accessibility_and_length",
           "background_method": "trans_chromosome"}
    txt = " ".join(IntegrationNarrator().methods(
        "integration_agent", _integ({"peak_to_gene": p2g})))
    assert "raw peak counts" in txt
    assert "two-sided normal tail p-values" in txt
    assert "one-sided" not in txt
    assert "reference-FASTA" not in txt
    assert "no reference FASTA" in txt

    p2g.update(atac_normalization="signac_tfidf",
               background_match="gc_and_accessibility")
    txt = " ".join(IntegrationNarrator().methods(
        "integration_agent", _integ({"peak_to_gene": p2g})))
    assert "log-TF-IDF" in txt and "reference-FASTA GC-fraction" in txt


# ═══════════════════════════════════════════════════════════════════════════
# Integration hypothesis adapter
# ═══════════════════════════════════════════════════════════════════════════

def test_adapter_mofa_excludes_batch_factors():
    from aria.agents.narrative.hypothesis.adapters.integration import (
        _mofa_excluded_ids, _mofa_signals,
    )

    loadings = {
        f"Factor{i}": {"scRNA": [{"feature": f"GENE{i}", "weight": 0.9}]}
        for i in (1, 2, 6)
    }
    mofa = {
        "views": {"scRNA": {"kind": "rna"}},
        "top_loadings": loadings,
        "batch_factor_ids": [2],
        "top_factors": [
            {"factor_id": 1, "active": True, "category": "shared",
             "r2_by_view": {"scRNA": 0.1}},
            {"factor_id": 2, "active": True, "category": "shared",
             "r2_by_view": {"scRNA": 0.5}},
            {"factor_id": 6, "active": True, "category": "shared",
             "batch_flag": True, "r2_by_view": {"scRNA": 0.4}},
        ],
    }
    assert "2" in _mofa_excluded_ids(mofa)
    signals = _mofa_signals(mofa, wnn_caveats=[], max_factors=5,
                            max_genes_per_factor=5)
    genes = {s.entity for s in signals}
    assert genes == {"GENE1"}


# ═══════════════════════════════════════════════════════════════════════════
# Hi-C narrator
# ═══════════════════════════════════════════════════════════════════════════

def _hic_findings() -> dict:
    return {
        "preprocessing": {"status": "success", "n_samples": 2,
                          "cooler_files": ["/c/a.mcool", "/c/b.mcool"],
                          "samples": [{"resolutions": [5000, 10000]},
                                      {"resolutions": [5000]}]},
        "balancing": {"status": "success", "n_samples": 0,
                      "mean_replicate_scc": float("nan"),
                      "n_valid_pairs": float("inf")},
        "compartments": {"status": "success", "pct_A": 55.0, "pct_B": 45.0,
                         "resolution": 100000, "pc1_validated": True,
                         "phasing_method": "gc_content"},
        "tads": {"status": "success", "n_tads": 100, "median_size_kb": 400.0,
                 "window_size": 50000},
        "loops": {"status": "success", "n_loops": 20, "resolution": 10000,
                  "loops_path": "/t/loops.bedpe",
                  "annotated_loops_path": "/t/loops_annot.tsv",
                  "loop_class_counts": {"promoter_enhancer": 7,
                                        "promoter_promoter": 3,
                                        "distal_distal": 10}},
        "differential_3d": {
            "status": "success", "ran": True,
            "comparisons": [{
                "contrast": "KO_vs_WT",
                "differential_compartments_csv": "/t/dc.csv",
                "differential_insulation_csv": "/t/di.csv",
                "n_differential_compartment_bins": float("nan"),
                "n_ab_switches": float("inf"),
                "n_differential_insulation_bins": 4,
            }],
        },
        "figures": {"hic_contact_map": "/f/map.png"},
        "figure_bundles": {"hic_contact_map": {"png": "/f/map.png",
                                               "svg": "/f/map.svg"}},
    }


def test_hic_as_int_nonfinite():
    from aria.agents.narrative.narrators.hic import _as_int

    assert _as_int(float("inf")) is None
    assert _as_int(float("nan")) is None
    assert _as_int(3.0) == 3


def test_hic_svg_from_figure_bundles():
    from aria.agents.narrative.narrators.hic import HiCNarrator

    figs = HiCNarrator().figures("genome_arch_agent",
                                 {"findings": _hic_findings()})
    rec = next(f for f in figs if f["id"] == "hic.hic_contact_map")
    assert rec["svg_path"] == "/f/map.svg"


def test_hic_blocks_canonical_keys_and_numeric_edges():
    from aria.agents.narrative.narrators.hic import HiCNarrator

    blocks = {b.id: b for b in HiCNarrator().collect(
        "genome_arch_agent", {"findings": _hic_findings()})}

    pre = blocks["hic.preprocessing"]
    assert pre.metrics["base_resolution"] == 5000

    qc = blocks["hic.qc"]
    assert qc.metrics["n_samples"] == 0  # explicit 0 is not coerced to 1

    tads = blocks["hic.tads"]
    assert tads.metrics["window_bp"] == 50000

    loops = blocks["hic.loops"]
    assert loops.metrics["n_promoter_enhancer_candidates"] == 7
    paths = {ev.path for ev in loops.evidence if ev.path}
    assert {"/t/loops.bedpe", "/t/loops_annot.tsv"} <= paths

    diff = blocks["hic.differential_3d"]  # int(nan)/int(inf) no longer crash
    assert diff.metrics["n_significant_compartment_bins"] == 0
    assert diff.metrics["n_significant_boundaries"] == 4

    comp = blocks["hic.compartments"]
    assert "55.0% active" in comp.claim


def test_hic_unphased_compartments_not_called_active():
    from aria.agents.narrative.narrators.hic import HiCNarrator

    f = _hic_findings()
    f["compartments"]["pc1_validated"] = False
    blocks = {b.id: b for b in HiCNarrator().collect(
        "genome_arch_agent", {"findings": f})}
    comp = blocks["hic.compartments"]
    assert "active" not in comp.claim.replace("not active/inactive", "")
    assert "positive-E1" in comp.claim and "unphased" in comp.claim
    labels = [ev.label for ev in comp.evidence]
    assert "Active compartment fraction" not in labels
    assert "Positive-E1 bin fraction (unphased sign)" in labels


def test_hic_tables_canonical_and_legacy_keys():
    from aria.agents.narrative.narrators.hic import HiCNarrator

    tables = HiCNarrator().tables("genome_arch_agent",
                                  {"findings": _hic_findings()})
    paths = {t["path"] for t in tables}
    assert {"/t/loops.bedpe", "/t/loops_annot.tsv",
            "/t/dc.csv", "/t/di.csv"} <= paths

    legacy = _hic_findings()
    legacy["loops"] = {"status": "success", "loops_bedpe": "/l/old.bedpe",
                       "annotated_loops_tsv": "/l/old.tsv"}
    legacy["differential_3d"]["comparisons"] = [
        {"contrast": "x", "compartments_csv": "/l/c.csv",
         "insulation_csv": "/l/i.csv"}]
    paths = {t["path"] for t in HiCNarrator().tables(
        "genome_arch_agent", {"findings": legacy})}
    assert {"/l/old.bedpe", "/l/old.tsv", "/l/c.csv", "/l/i.csv"} <= paths


def test_hic_methods_no_placeholders():
    from aria.agents.narrative.narrators.hic import HiCNarrator

    f = _hic_findings()
    f["balancing"]["mean_replicate_scc"] = float("nan")
    txt = " ".join(HiCNarrator().methods("genome_arch_agent", {"findings": f}))
    assert "base bin 5000 bp" in txt
    assert "window 50000 bp" in txt
    assert "phased against gc_content" in txt
    assert "?" not in txt
    assert "nan" not in txt.lower().replace("financ", "")
    assert "mean SCC" not in txt

    f["balancing"]["mean_replicate_scc"] = 0.912
    txt = " ".join(HiCNarrator().methods("genome_arch_agent", {"findings": f}))
    assert "mean SCC 0.912" in txt

    legacy = _hic_findings()
    legacy["compartments"].pop("phasing_method")
    legacy["compartments"]["phasing_track"] = "tss_density"
    legacy["tads"] = {"status": "success", "window_bp": 25000}
    legacy["preprocessing"]["base_resolution"] = 1000
    txt = " ".join(HiCNarrator().methods(
        "genome_arch_agent", {"findings": legacy}))
    assert "phased against tss_density" in txt
    assert "window 25000 bp" in txt
    assert "base bin 1000 bp" in txt


def test_hic_qc_per_file_fallback():
    from aria.agents.narrative.narrators.hic import HiCNarrator

    bal = {"status": "success", "per_file": [
        {"status": "success"}, {"status": "success"}, {"status": "error"}]}
    block = HiCNarrator()._qc_block(bal, {})
    assert block.metrics["n_samples"] == 2
    assert not math.isnan(block.metrics["n_samples"])
