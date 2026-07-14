"""A7 characterization: pin the scRNAAgent contract before/through extraction.

``scrna_agent.py`` (2.3k lines) is the second A7 giant. Unlike ``_narrative_scrna``
it is one stateful ``scRNAAgent`` class, so the extraction splits it into concern
mixins (qc/annotation/de/pseudobulk/advanced) that the composed class inherits,
behind the ``aria/agents/scrna_agent.py`` facade. This file locks two things a
mixin split must preserve:

1. Structural contract — the composed class still exposes EVERY method (nothing is
   dropped or renamed when methods move into mixins).
2. Behavior of the pure routing/policy methods that decide which analyses run and
   how comparisons are normalised (the parts a split is most likely to break).

The agent never imports scanpy (subprocess-only), so these run in the standard
env. Heavy pipeline behavior stays covered by the scanpy/e2e suites.
"""
from __future__ import annotations

from aria.agents.scrna_agent import scRNAAgent


# ── 1. Structural contract: every method survives the split ───────────────────

_EXPECTED_METHODS = [
    # entry / core
    "run", "receive", "_workspace", "_log_decision",
    "_design_intelligence_blocks", "_design_intelligence_optional_selected",
    # focus routing
    "_prepare_focused_h5ads", "_scrna_focus_workspace", "_design_groupby_col",
    "_infer_cell_focus_values", "_cell_focus_text", "_available_groupby_values",
    "_sample_id_from_path",
    # qc / integration / clustering
    "_run_qc", "_qc_single", "_resolve_batch_column", "_run_integration",
    "_run_clustering",
    # annotation
    "_infer_tissue_hint", "_allow_default_immune_model", "_predefined_celltype_col",
    "_annotation_from_obs", "_annotate_cell_types", "_marker_based_annotation",
    "_parse_annotation_json", "_trusted_annotation_groupby",
    "_annotation_is_report_only",
    # de / pathway
    "_differential_expression", "_run_pathway_per_cluster",
    "_expressed_gene_background",
    # pseudobulk
    "_needs_pseudobulk", "_inject_condition_obs", "_integration_qc_has_blocking",
    "_integration_qc_blocking_reason", "_run_pseudobulk",
    "_normalise_pseudobulk_comparisons", "_suggest_pseudobulk_comparisons",
    # advanced
    "_needs_trajectory", "_run_trajectory", "_needs_cell_communication",
    "_run_cell_communication",
]


def test_composed_agent_exposes_every_method():
    missing = [m for m in _EXPECTED_METHODS if not hasattr(scRNAAgent, m)]
    assert missing == [], f"methods lost in the split: {missing}"


def test_class_constants_and_identity_preserved():
    assert scRNAAgent.name == "scrna_agent"
    assert isinstance(scRNAAgent.PSEUDOBULK_KEYWORDS, tuple)
    assert "clustered" in scRNAAgent._ARIA_INTERMEDIATE_STEMS


def _bare():
    """A scRNAAgent instance without __init__ (pure methods need no state)."""
    return scRNAAgent.__new__(scRNAAgent)


# ── 2. Pure routing/policy behavior ───────────────────────────────────────────

def test_needs_trajectory_keyword_routing():
    agent = _bare()
    assert agent._needs_trajectory({"summary": "study the differentiation lineage"})
    assert agent._needs_trajectory(
        {"summary": "x", "biological_entities": ["pseudotime progression"]}
    )
    assert not agent._needs_trajectory({"summary": "compare disease vs healthy"})


def test_needs_cell_communication_keyword_routing():
    agent = _bare()
    assert agent._needs_cell_communication(
        {"summary": "ligand receptor crosstalk in the niche"}
    )
    assert not agent._needs_cell_communication(
        {"summary": "cluster the cells and annotate"}
    )


def test_normalise_pseudobulk_comparisons_shapes():
    norm = scRNAAgent._normalise_pseudobulk_comparisons
    assert norm(None) == []
    assert norm([{"test": "A", "ref": "B"}]) == [["A", "B"]]
    assert norm([{"case": "A", "control": "B"}]) == [["A", "B"]]
    assert norm([("A", "B", "extra")]) == [["A", "B"]]
    assert norm([{"test": "A"}]) == []  # incomplete pair dropped


def test_suggest_pseudobulk_comparisons_prefers_control_reference():
    sug = scRNAAgent._suggest_pseudobulk_comparisons
    assert sug({}) == []
    assert sug({"only": {}}) == []
    # a control-like level becomes the shared reference
    assert sug({"treated": {}, "control": {}}) == [["treated", "control"]]
    # otherwise all ordered pairs, each as [later, earlier] over sorted names
    assert sug({"A": {}, "B": {}, "C": {}}) == [["B", "A"], ["C", "A"], ["C", "B"]]
