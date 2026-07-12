"""Preprint-readiness audit B2: unresolved confounding must block, not dispatch.

Two silent failures are closed here:

  * A batch (or any covariate) perfectly aliased with the biological condition
    is non-identifiable.  The headless answer policy defaults to "ignore batch",
    but ignoring a confounded factor does not make the biology estimable — it
    hides the confound.  The assembled design must fail closed regardless of the
    ignore-batch choice.
  * GEO/SRA metadata inference must not collapse a multifactorial study down to a
    single "best" characteristic and silently discard donor/batch/covariates.
    Secondary design factors are preserved in an explicit sample sheet, and a
    secondary factor completely confounded with the chosen condition is flagged
    so the design phase blocks instead of dispatching on the single best feature.
"""
from __future__ import annotations


# ── Shared pure helper ──────────────────────────────────────────────────────

def test_confounding_helper_detects_bijective_alias():
    from aria.utils.design_matrix import factors_confounded_with_condition

    # batch maps 1:1 onto condition — completely confounded.
    sheet = [
        {"sample": "s1", "condition": "ctrl", "batch": "b1", "sex": "F"},
        {"sample": "s2", "condition": "ctrl", "batch": "b1", "sex": "M"},
        {"sample": "s3", "condition": "stim", "batch": "b2", "sex": "F"},
        {"sample": "s4", "condition": "stim", "batch": "b2", "sex": "M"},
    ]
    confounded = factors_confounded_with_condition(
        sheet, "condition", ["batch", "sex"]
    )
    assert confounded == ["batch"]


def test_confounding_helper_ignores_balanced_covariate():
    from aria.utils.design_matrix import factors_confounded_with_condition

    # donor crosses condition (balanced block) — identifiable, not confounded.
    sheet = [
        {"sample": "s1", "condition": "ctrl", "donor": "d1"},
        {"sample": "s2", "condition": "stim", "donor": "d1"},
        {"sample": "s3", "condition": "ctrl", "donor": "d2"},
        {"sample": "s4", "condition": "stim", "donor": "d2"},
    ]
    assert factors_confounded_with_condition(sheet, "condition", ["donor"]) == []


# ── GEO metadata preserves secondary factors and flags confounding ──────────

def _confounded_geo_metadata():
    # condition (genotype) is perfectly aliased with the processing batch, and a
    # donor id is present.  Neither may be discarded.
    return {
        "samples": [
            {"id": "GSM01", "title": "WT rep1",
             "characteristics": {"genotype": "WT", "batch": "day1",
                                 "donor": "p1", "sex": "F"}},
            {"id": "GSM02", "title": "WT rep2",
             "characteristics": {"genotype": "WT", "batch": "day1",
                                 "donor": "p2", "sex": "M"}},
            {"id": "GSM03", "title": "KO rep1",
             "characteristics": {"genotype": "KO", "batch": "day2",
                                 "donor": "p3", "sex": "F"}},
            {"id": "GSM04", "title": "KO rep2",
             "characteristics": {"genotype": "KO", "batch": "day2",
                                 "donor": "p4", "sex": "M"}},
        ]
    }


def test_geo_infer_design_preserves_secondary_factors():
    from aria.connectors.geo_connector import _infer_design

    design = _infer_design(_confounded_geo_metadata())

    # The primary condition is still chosen.
    assert design["condition_col"] == "genotype"

    # A per-sample sample sheet carries every characteristic — nothing dropped.
    sheet = {row["sample"]: row for row in design["sample_sheet"]}
    assert set(sheet) == {"GSM01", "GSM02", "GSM03", "GSM04"}
    assert sheet["GSM01"]["batch"] == "day1"
    assert sheet["GSM01"]["donor"] == "p1"
    assert sheet["GSM03"]["genotype"] == "KO"

    # Secondary multi-level factors are surfaced, not silently discarded.
    assert "batch" in design["covariates"]
    assert "donor" in design["covariates"]


def test_geo_infer_design_flags_confounded_covariate():
    from aria.connectors.geo_connector import _infer_design

    design = _infer_design(_confounded_geo_metadata())

    # batch is perfectly aliased with genotype -> flagged as unresolved.
    assert "batch" in design["confounded_covariates"]
    assert design["unresolved_confounding"] is True


# ── DesignAgent fails closed on confounded metadata ─────────────────────────

def _design_agent(*, groups, inferred_design, batch_covariate=None):
    from aria.agents.design_agent import DesignAgent

    agent = object.__new__(DesignAgent)
    agent._confirmed_groups = groups
    agent._main_factor = "genotype"
    agent._batch_covariate = batch_covariate
    agent._organism = "homo sapiens"
    agent._genome = "hg38"
    agent._pseudobulk_design = {}
    agent._inferred_design = inferred_design
    agent._replicate_handling = None
    agent._publish_confirm_checkpoint = lambda: None
    return agent


def test_build_design_blocks_confounded_metadata():
    from aria.connectors.geo_connector import _infer_design

    inferred = _infer_design(_confounded_geo_metadata())
    agent = _design_agent(
        groups={"WT": ["GSM01", "GSM02"], "KO": ["GSM03", "GSM04"]},
        inferred_design=inferred,
    )
    design = agent._build_design()

    assert design["design_status"] == "blocking"
    assert "condition_covariate_confounding" in design["blocking_reasons"]


def test_ignore_batch_choice_cannot_bypass_confounding():
    """Even with the headless 'ignore batch' default (batch_covariate=None), a
    metadata-level confound must still block the assembled design."""
    from aria.connectors.geo_connector import _infer_design

    inferred = _infer_design(_confounded_geo_metadata())
    agent = _design_agent(
        groups={"WT": ["GSM01", "GSM02"], "KO": ["GSM03", "GSM04"]},
        inferred_design=inferred,
        batch_covariate=None,   # headless default: "No — ignore batch"
    )
    design = agent._build_design()

    assert design["design_status"] == "blocking"
    assert "condition_covariate_confounding" in design["blocking_reasons"]


def test_confirm_response_reports_confounding_reason():
    """The confirm handler must fail closed with the confounding reason, not the
    B1 replicate label, so headless/TUI cancel with an accurate cause."""
    from aria.connectors.geo_connector import _infer_design

    inferred = _infer_design(_confounded_geo_metadata())
    agent = _design_agent(
        groups={"WT": ["GSM01", "GSM02"], "KO": ["GSM03", "GSM04"]},
        inferred_design=inferred,
    )
    result = agent._handle_confirm_response("Yes — proceed")
    assert result["status"] == "cancelled"
    assert result["reason"] == "condition_covariate_confounding"


def test_build_design_allows_balanced_block():
    """A donor that crosses condition is a valid block, not a confound."""
    balanced = {
        "sample_sheet": [
            {"sample": "GSM01", "genotype": "WT", "donor": "d1"},
            {"sample": "GSM02", "genotype": "KO", "donor": "d1"},
            {"sample": "GSM03", "genotype": "WT", "donor": "d2"},
            {"sample": "GSM04", "genotype": "KO", "donor": "d2"},
        ],
        "covariates": ["donor"],
        "confounded_covariates": [],
        "unresolved_confounding": False,
        "condition_col": "genotype",
    }
    agent = _design_agent(
        groups={"WT": ["GSM01", "GSM03"], "KO": ["GSM02", "GSM04"]},
        inferred_design=balanced,
    )
    design = agent._build_design()

    assert design["design_status"] == "ready"
    assert "condition_covariate_confounding" not in design.get(
        "blocking_reasons", []
    )
