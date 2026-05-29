"""scRNA pseudobulk dispatch gate honors the approved plan.

Regression guard for the 2026-05-28 PBMC blocker: `_needs_pseudobulk` used to
re-gate purely on PSEUDOBULK_KEYWORDS in the free-text question, silently
overriding a DesignIntelligence recommendation that the user had approved at
CP2. An IFN-β "response programs / signaling networks" question carries no
keyword, so DE was dropped without a logged decision or reported skip.
"""

from aria.agents.scrna_agent import scRNAAgent


def _agent():
    return scRNAAgent.__new__(scRNAAgent)


_OBS_DESIGN = {
    "groups": {"CTRL": ["d1", "d2", "d3"], "STIM": ["d1", "d2", "d3"]},
    "pseudobulk": {
        "condition_col": "stim",
        "replicate_col": "Donor",
        "groupby_col": "cluster",
        "comparisons": [["STIM", "CTRL"]],
    },
}

# A descriptive question with NO comparison keyword (the blocker's phrasing).
_DESCRIPTIVE_INTENT = {
    "summary": "Map cell-type-resolved interferon-beta transcriptional response "
               "programs and infer myeloid-lymphoid signaling networks.",
    "biological_entities": [],
}


def test_pseudobulk_runs_when_design_intelligence_recommended_it():
    agent = _agent()
    exp_ctx = {
        "design": _OBS_DESIGN,
        "run_optional_supported": True,
        "user_question": _DESCRIPTIVE_INTENT["summary"],
        "design_intelligence": {
            "recommended": [
                "QC using available h5ad obs metrics or count data.",
                "Donor/sample-level pseudobulk DE: condition=stim, "
                "replicate=Donor, groupby=cluster.",
            ],
        },
    }
    assert agent._needs_pseudobulk(_DESCRIPTIVE_INTENT, exp_ctx) is True


def test_pseudobulk_falls_back_to_keywords_without_di_recommendation():
    """No DI recommendation + keyword-free descriptive question -> do not run."""
    agent = _agent()
    exp_ctx = {
        "design": _OBS_DESIGN,
        "run_optional_supported": True,
        "user_question": _DESCRIPTIVE_INTENT["summary"],
        "design_intelligence": {"recommended": []},
    }
    assert agent._needs_pseudobulk(_DESCRIPTIVE_INTENT, exp_ctx) is False


def test_pseudobulk_keyword_fallback_still_runs_on_comparison_question():
    agent = _agent()
    intent = {"summary": "Compare STIM versus CTRL: differential expression.",
              "biological_entities": []}
    exp_ctx = {
        "design": _OBS_DESIGN,
        "run_optional_supported": True,
        "user_question": intent["summary"],
        "design_intelligence": {"recommended": []},
    }
    assert agent._needs_pseudobulk(intent, exp_ctx) is True
