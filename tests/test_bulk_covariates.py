"""P0-4 regression: bulk DESeq2 must honor confirmed covariates/batch.

`_run_deseq2` hardcoded `design="~ {factor}"` and called the design-matrix
validator with `covariates=[]`, so a batch covariate confirmed by the user at
DesignAgent CHECKPOINT 2.4 (recorded as `design_formula = "~ batch + condition"`)
was silently dropped — the fitted model adjusted for nothing and the report did
not disclose it. These tests lock covariate propagation: the fitted formula
includes the covariate, the validator receives it, an unavailable covariate is
disclosed (not silently ignored), and the agent forwards the confirmed batch.
"""

import pandas as pd
import pytest


# ── Pure helpers (no pydeseq2) ───────────────────────────────────────────────

def test_build_design_formula_places_factor_of_interest_last():
    from aria.scripts.rna_bulk_de import _build_design_formula
    assert _build_design_formula("condition", []) == "~ condition"
    assert _build_design_formula("condition", ["batch"]) == "~ batch + condition"
    # Dedup + never repeat the factor of interest as a covariate.
    assert _build_design_formula(
        "condition", ["batch", "batch", "condition"]
    ) == "~ batch + condition"


def test_resolve_covariates_keeps_usable_and_discloses_dropped():
    from aria.scripts.rna_bulk_de import _resolve_covariates
    meta = pd.DataFrame({
        "condition": ["A", "A", "B", "B"],
        "batch":     ["b1", "b2", "b1", "b2"],   # varies -> usable
        "lane":      ["L1", "L1", "L1", "L1"],   # constant -> dropped
    }, index=["s1", "s2", "s3", "s4"])

    usable, dropped = _resolve_covariates(meta, "condition",
                                          ["batch", "lane", "missing"])
    assert usable == ["batch"]
    dropped_names = {name for name, _ in dropped}
    assert dropped_names == {"lane", "missing"}   # constant + absent disclosed


def test_resolve_covariates_never_returns_the_design_factor():
    from aria.scripts.rna_bulk_de import _resolve_covariates
    meta = pd.DataFrame({"condition": ["A", "A", "B", "B"]},
                        index=["s1", "s2", "s3", "s4"])
    usable, _ = _resolve_covariates(meta, "condition", ["condition"])
    assert usable == []


# ── Agent forwards the confirmed batch covariate ─────────────────────────────

def test_agent_design_covariates_from_confirmed_design():
    pytest.importorskip("litellm")  # importing the agent pulls aria.llm.provider
    from aria.agents.bulk_rna_agent import BulkRNAAgent
    assert BulkRNAAgent._design_covariates(
        {"batch_covariate": "batch"}) == ["batch"]
    assert BulkRNAAgent._design_covariates({"batch_covariate": None}) == []
    assert BulkRNAAgent._design_covariates({}) == []
    # An explicit covariates list is merged with the batch covariate, deduped.
    assert BulkRNAAgent._design_covariates(
        {"batch_covariate": "batch", "covariates": ["batch", "sex"]}
    ) == ["batch", "sex"]


# ── End-to-end DESeq2 with a real covariate (pydeseq2-gated) ─────────────────

def test_deseq2_fits_and_reports_the_covariate_adjusted_formula():
    pytest.importorskip("pydeseq2")
    import numpy as np
    from aria.scripts.rna_bulk_de import _run_deseq2

    rng = np.random.default_rng(0)
    samples = [f"s{i}" for i in range(8)]
    metadata = pd.DataFrame({
        "condition": ["ctrl"] * 4 + ["treat"] * 4,
        # batch crosses condition (not confounded): adjustable, not rank-deficient
        "batch":     ["b1", "b2", "b1", "b2", "b1", "b2", "b1", "b2"],
    }, index=samples)
    counts = pd.DataFrame(
        rng.poisson(100, size=(200, 8)),
        index=[f"GENE_{i:04d}" for i in range(200)],
        columns=samples,
    )

    result, warnings = _run_deseq2(
        counts, metadata, "condition", "treat", "ctrl",
        padj_thr=0.05, lfc_thr=1.0, covariates=["batch"],
    )
    assert result["status"] == "success", result
    assert result["fitted_design_formula"] == "~ batch + condition"
    assert result["covariates_adjusted"] == ["batch"]


# ── Report Methods state the actually-fitted formula ─────────────────────────

def test_bulk_methods_report_the_fitted_covariate_formula():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
    agent_result = {"findings": {"contrasts": [{
        "name": "treat vs ctrl",
        "status": "success",
        "fitted_design_formula": "~ batch + condition",
        "covariates_adjusted": ["batch"],
        "covariates_dropped": [],
    }]}}
    methods = BulkRnaNarrator().methods("bulk_rna_agent", agent_result)
    assert any("~ batch + condition" in m for m in methods)


def test_bulk_methods_disclose_a_dropped_covariate():
    from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
    agent_result = {"findings": {"contrasts": [{
        "name": "treat vs ctrl",
        "status": "success",
        "fitted_design_formula": "~ condition",
        "covariates_adjusted": [],
        "covariates_dropped": [{"covariate": "batch",
                                "reason": "not present in the sample metadata"}],
    }]}}
    methods = BulkRnaNarrator().methods("bulk_rna_agent", agent_result)
    assert any("not adjusted" in m and "batch" in m for m in methods)
