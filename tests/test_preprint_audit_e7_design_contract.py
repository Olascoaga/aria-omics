"""Preprint-readiness audit E7: one shared RNA/ATAC design contract.

The confirmed design (DesignAgent._build_design) is the canonical contract:
main_factor / batch_covariate / pseudobulk.{condition_col,replicate_col} /
covariates. Before E7 the ATAC lanes re-derived those fields ad-hoc:

* the scATAC lane walked a scattered exp_ctx/design/pseudobulk fallback fan with
  key-name drift (batch_col / batch_factor vs the contract's batch_covariate) and
  forwarded only condition/replicate/comparisons to DA — dropping the batch
  covariate entirely;
* the bulk ATAC lane ignored exp_ctx["design"] and defaulted condition_col /
  replicate_col to the literals "condition" / "replicate", never adding the
  confirmed batch.

After E7 a single pure resolver (resolve_design_contract) normalizes the design
for RNA and both ATAC lanes, and the confirmed batch covariate reaches the scATAC
pseudobulk DESeq2 fit like it does for RNA.
"""
from __future__ import annotations

import pytest


# ── The shared resolver (pure, dependency-light) ──────────────────────────────

def _resolve():
    from aria.utils.design_matrix import resolve_design_contract
    return resolve_design_contract


def _design():
    # Shape emitted by DesignAgent._build_design for a donor-paired, batched design.
    return {
        "main_factor": "genotype",
        "batch_covariate": "batch",
        "pseudobulk": {"condition_col": "genotype", "replicate_col": "donor"},
        "groups": {"KO": ["d1", "d2"], "WT": ["d3", "d4"]},
    }


def test_resolver_reads_condition_replicate_from_contract():
    contract = _resolve()(_design(), {})
    assert contract["condition_col"] == "genotype"
    assert contract["replicate_col"] == "donor"


def test_resolver_carries_batch_as_covariate():
    contract = _resolve()(_design(), {})
    # The confirmed batch covariate must survive as a modelled covariate.
    assert "batch" in contract["covariates"]


def test_resolver_never_lists_condition_or_replicate_as_covariate():
    contract = _resolve()(_design(), {})
    assert "genotype" not in contract["covariates"]
    assert "donor" not in contract["covariates"]


def test_resolver_dedupes_and_preserves_order():
    design = _design()
    design["covariates"] = ["batch", "sex"]  # batch already the batch_covariate
    contract = _resolve()(design, {"covariates": ["sex", "prep"]})
    assert contract["covariates"] == ["batch", "sex", "prep"]


def test_resolver_exp_ctx_overrides_design():
    contract = _resolve()(
        _design(),
        {"condition_col": "treatment", "replicate_col": "subject"},
    )
    assert contract["condition_col"] == "treatment"
    assert contract["replicate_col"] == "subject"


def test_resolver_tolerates_missing_design():
    contract = _resolve()(None, {"condition_col": "cond"})
    assert contract["condition_col"] == "cond"
    assert contract["covariates"] == []


# ── Both ATAC lanes resolve the DA design through the shared contract ──────────

def _chromatin_agent():
    pytest.importorskip("litellm")  # importing the agent pulls aria.llm.provider
    from aria.agents.chromatin_agent import ChromatinAgent
    return object.__new__(ChromatinAgent)


def test_chromatin_da_design_uses_the_contract():
    agent = _chromatin_agent()
    exp_ctx = {"design": _design()}
    resolved = agent._resolve_da_design(exp_ctx)
    # No hardcoded "condition"/"replicate" defaults; batch not dropped.
    assert resolved["condition_col"] == "genotype"
    assert resolved["replicate_col"] == "donor"
    assert "batch" in resolved["covariates"]


def test_chromatin_da_design_no_literal_default_when_contract_present():
    agent = _chromatin_agent()
    # A confirmed design whose factor is not literally "condition"/"replicate"
    # must not silently fall back to those names.
    exp_ctx = {"design": _design()}
    resolved = agent._resolve_da_design(exp_ctx)
    assert resolved["condition_col"] != "condition"
    assert resolved["replicate_col"] != "replicate"


# ── The confirmed batch reaches the scATAC pseudobulk DESeq2 fit ──────────────

def test_scatac_pseudobulk_models_confirmed_batch():
    pytest.importorskip("pydeseq2")
    pytest.importorskip("anndata")
    import numpy as np
    import anndata as ad
    from aria.scripts.chromatin_diffacc import _pseudobulk_da

    rng = np.random.default_rng(0)
    # 4 donors × 2 conditions, batch confounded-free (2 batches crossing condition),
    # enough cells per donor for the min-cells filter.
    donors = ["d1", "d2", "d3", "d4"]
    cond = {"d1": "KO", "d2": "KO", "d3": "WT", "d4": "WT"}
    batch = {"d1": "b1", "d2": "b2", "d3": "b1", "d4": "b2"}
    obs_rows = []
    for d in donors:
        for _ in range(40):
            obs_rows.append((d, cond[d], batch[d]))
    import pandas as pd
    obs = pd.DataFrame(obs_rows, columns=["donor", "condition", "batch"])
    X = rng.poisson(5, size=(len(obs), 30)).astype(float)
    adata = ad.AnnData(X=X, obs=obs.reset_index(drop=True))

    result = _pseudobulk_da(
        adata, "condition", "donor",
        comparisons=[{"test": "KO", "reference": "WT"}],
        padj_max=0.1, lfc_min=0.0, top_n=10, min_cells=5, min_reps=2,
        allow_mock=False, warnings=[], covariates=["batch"],
    )
    assert result.get("ran") is True
    comp = result["comparisons"][0]
    # The batch term must appear in the fitted model, not be dropped.
    assert "batch" in (comp.get("fitted_design_formula") or "")
