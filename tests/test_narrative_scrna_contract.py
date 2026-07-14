"""A7 characterization: pin the public contract of ``aria.agents._narrative_scrna``.

FASE 7 slice A7 refactors giant modules (1.8-2.3k lines) that mix orchestration,
statistics, policy and rendering. The tracker rule is: pin the typed contracts
with characterization tests BEFORE extracting anything, so a later split into a
subpackage is provably behavior-preserving.

``_narrative_scrna`` is a module of module-level functions consumed by
``narrative_agent``, ``devils_advocate`` and ``_narrative_chromatin``. This file
locks the input->output contract of its public surface (and the "private"
helpers other modules already import) against a representative findings fixture,
plus the honest-null behavior on empty findings. These assertions describe what
the module does TODAY; they are the safety net the extraction must keep green.
"""
from __future__ import annotations

import pytest

from aria.agents import _narrative_scrna as ns


# ── Representative findings fixture (one scrna_agent.run() direct emit) ────────

def _findings() -> dict:
    return {
        "qc": {
            "n_cells_before": 3000,
            "n_cells_after": 2600,
            "pct_removed": 13,
            "n_samples": 2,
            "mt_threshold": 15,
        },
        "clustering": {"n_clusters": 7, "resolution": 0.5},
        "cell_types": {
            "label_col": "cell_type",
            "cell_types": {
                "0": {"cell_type": "T cell", "annotation_source": "celltypist"},
                "1": {"cell_type": "B cell", "annotation_source": "celltypist"},
                "2": {"cell_type": "annotation_failed",
                      "annotation_source": "celltypist"},
            },
        },
        "pseudobulk_de": {
            "groupby": "cell_type",
            "per_group": {
                "T cell": {
                    "per_comparison": {
                        "COND_A_vs_COND_B": {
                            "status": "success",
                            "n_significant": 12,
                            "top_genes": [
                                {"gene": "GeneUp", "symbol": "GeneUp",
                                 "log2fc": 2.1, "padj_global": 1e-4},
                                {"gene": "GeneDn", "symbol": "GeneDn",
                                 "log2fc": -1.7, "padj_global": 1e-3},
                            ],
                        }
                    }
                },
                "B cell": {
                    "per_comparison": {
                        "COND_A_vs_COND_B": {
                            "status": "success",
                            "n_significant": 3,
                            "top_genes": [
                                {"gene": "GeneX", "symbol": "GeneX",
                                 "log2fc": 1.2, "padj_global": 2e-3},
                            ],
                        }
                    }
                },
            },
        },
        "cell_communication": {
            "top_interactions": [
                {"source": "T cell", "target": "B cell",
                 "ligand": "LIG1", "receptor": "REC1"},
                {"source": "B cell", "target": "T cell",
                 "ligand": "LIG2", "receptor": "REC2"},
            ]
        },
    }


# ── unwrap_scrna_findings: robust to both envelope shapes ──────────────────────

def test_unwrap_accepts_multimodal_wrapped_shape():
    inner = {"qc": {"n_cells_after": 10}}
    wrapped = {"findings": {"scRNA": {"findings": inner}}}
    assert ns.unwrap_scrna_findings(wrapped) == inner


def test_unwrap_accepts_direct_run_shape():
    inner = {"qc": {"n_cells_after": 10}, "clustering": {"n_clusters": 3}}
    assert ns.unwrap_scrna_findings({"findings": inner}) == inner


def test_unwrap_empty_envelope_returns_empty_dict():
    assert ns.unwrap_scrna_findings({}) == {}
    assert ns.unwrap_scrna_findings({"findings": {}}) == {}


# ── _group_label: crisp obs-column label mapping ──────────────────────────────

@pytest.mark.parametrize(
    "groupby,n,expected",
    [
        ("leiden", 5, "Leiden clusters"),
        ("leiden", 1, "Leiden cluster"),
        ("cell_type", 3, "cell types"),
        ("celltype", 3, "cell types"),
        ("cell_type_celltypist", 3, "cell types"),
        ("condition", 4, "condition groups"),
        ("condition", 1, "condition group"),
        (None, 2, "groups"),
        ("", 2, "groups"),
    ],
)
def test_group_label_mapping(groupby, n, expected):
    assert ns._group_label(groupby, n) == expected


# ── _top_de_blocks: filter to success, sort by n_significant desc, limit ──────

def test_top_de_blocks_filters_sorts_and_limits():
    pb = _findings()["pseudobulk_de"]
    # add a non-success comparison that must be dropped
    pb["per_group"]["T cell"]["per_comparison"]["OTHER"] = {
        "status": "error", "n_significant": 999,
    }
    blocks = ns._top_de_blocks(pb, limit=5)
    assert [(g, c) for g, c, _ in blocks] == [
        ("T cell", "COND_A_vs_COND_B"),
        ("B cell", "COND_A_vs_COND_B"),
    ]
    assert ns._top_de_blocks(pb, limit=1)[0][0] == "T cell"


# ── _annotation_state: valid-label filtering + marker-fallback flag ───────────

def test_annotation_state_filters_invalid_and_reports_source():
    state = ns._annotation_state(_findings())
    assert state["has_valid"] is True
    assert set(state["labels"]) == {"T cell", "B cell"}  # annotation_failed dropped
    assert state["n_unique"] == 2
    assert state["source"] == "celltypist"
    assert state["is_marker_fallback"] is False
    assert state["label_col"] == "cell_type"


def test_annotation_state_flags_marker_fallback():
    findings = {
        "cell_types": {
            "cell_types": {
                "0": {"cell_type": "T cell",
                      "annotation_source": "marker_fallback"},
            }
        }
    }
    state = ns._annotation_state(findings)
    assert state["is_marker_fallback"] is True


# ── summarize_scrna_text: multi-line summary + honest null ────────────────────

def test_summarize_reports_qc_clustering_and_celltypes():
    text = ns.summarize_scrna_text(_findings())
    assert isinstance(text, str)
    assert "2,600 of 3,000 cells were retained" in text
    assert "7 clusters" in text
    assert "resolution 0.5" in text
    # only the two valid labels are surfaced
    assert "T cell" in text and "B cell" in text
    assert "annotation_failed" not in text


def test_summarize_empty_findings_returns_generic_completion_no_fabrication():
    # Empty findings yield a neutral completion notice, never fabricated numbers.
    assert ns.summarize_scrna_text({}) == (
        "scRNA analysis completed. See findings table for details."
    )


# ── build_scrna_methods: deterministic methods prose + honest null ────────────

def test_methods_describe_scanpy_pipeline():
    methods = ns.build_scrna_methods(_findings())
    assert isinstance(methods, str)
    assert "scanpy" in methods
    assert "Scrublet" in methods


def test_methods_empty_findings_is_empty():
    assert ns.build_scrna_methods({}) == ""


# ── build_scrna_integrated_interpretation: deterministic synthesis + null ─────

def test_integrated_interpretation_is_string_and_null_safe():
    # Empty findings yield a generic framing sentence, not fabricated results.
    assert ns.build_scrna_integrated_interpretation({}) == (
        "Integrated interpretation: ARIA addressed the submitted single-cell "
        "RNA-seq question using the structured scRNA outputs available in this "
        "run."
    )
    synthesis = ns.build_scrna_integrated_interpretation(_findings())
    assert isinstance(synthesis, str)
    assert synthesis  # a populated fixture yields a non-empty synthesis


# ── HTML table extractors: rows on data, empty string on no data ──────────────

def test_pseudobulk_de_table_rows_and_null():
    assert ns.extract_pseudobulk_de_table({}) == ""
    rows = ns.extract_pseudobulk_de_table(_findings())
    assert isinstance(rows, str) and rows
    assert "T cell" in rows


def test_cellcomm_table_rows_and_null():
    assert ns.extract_cellcomm_table({}) == ""
    rows = ns.extract_cellcomm_table(_findings())
    assert isinstance(rows, str) and rows
    assert "LIG1" in rows and "REC1" in rows
