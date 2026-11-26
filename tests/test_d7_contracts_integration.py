"""Regression tests for D7 layer-contract fixes in multiome integration."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


def test_d7_018_warning_router_surfaces_mofa_batch_condition_confounded(tmp_path: Path):
    """D7-018: integration_mofa detects 1:1 batch/donor-to-condition confounding
    and sets batch_condition_confounded=True. Through NarrativeRegistry.collect_blocks,
    warning_router must surface a structured warning caveat and cap confidence
    from high to medium even when free-text warnings are absent."""
    integ_py = Path("/home/medusa/anaconda3/envs/aria-integration-env/bin/python")
    if not integ_py.is_file():
        pytest.skip("aria-integration-env not available")
    pytest.importorskip("litellm")
    pytest.importorskip("pydantic")
    from aria.agents.narrative.narrators.integration import IntegrationNarrator
    from aria.agents.narrative.registry import registry_with
    from aria.utils.environment_manager import env_manager

    h5mu_path = tmp_path / "confounded.h5mu"
    build_code = f"""
import anndata as ad, mudata as md, numpy as np, pandas as pd
rng = np.random.default_rng(29)
n_cells = 120
donor = ['D1'] * 60 + ['D2'] * 60
cond = ['ctrl'] * 60 + ['treat'] * 60
z = rng.normal(size=(n_cells, 2))
z[60:, 0] += 2.0
rna_x = np.clip(np.round(np.exp(1.5 + 0.4 * (z @ rng.normal(size=(2, 30))))), 0, 500).astype(np.int32)
atac_x = np.clip(np.round(np.exp(1.0 + 0.4 * (z @ rng.normal(size=(2, 30))))), 0, 500).astype(np.int32)
barcodes = [f'cell_{{i}}' for i in range(n_cells)]
obs = pd.DataFrame({{'donor': donor, 'condition': cond}}, index=barcodes)
mdata = md.MuData({{
    'rna': ad.AnnData(X=rna_x, obs=obs.copy(), var=pd.DataFrame(index=[f'G{{i}}' for i in range(30)])),
    'atac': ad.AnnData(X=atac_x, obs=obs.copy(), var=pd.DataFrame(index=[f'chr1:{{i*1000}}-{{i*1000+500}}' for i in range(30)])),
}})
mdata.obs['donor'] = donor
mdata.obs['condition'] = cond
mdata.write_h5mu({str(h5mu_path)!r})
"""
    subprocess.run([str(integ_py), "-c", build_code], check=True)
    mofa_res = env_manager.run_in_stack(
        "integration",
        "aria/scripts/integration_mofa.py",
        {
            "h5mu_path": str(h5mu_path),
            "output_dir": str(tmp_path / "mofa_out"),
            "n_factors": 3,
            "max_iter": 30,
            "batch_key": "donor",
            "condition_col": "condition",
            "seed": 0,
        },
    )
    assert mofa_res["status"] == "success", mofa_res
    assert mofa_res["batch_condition_confounded"] is True

    # Strip free-text warnings and raise n_cells>=500 to verify the STRUCTURED
    # batch_condition_confounded signal both adds a warning caveat AND caps
    # confidence from 'high' to 'medium'.
    mofa_structured_only = {**mofa_res, "warnings": [], "n_cells": 600}
    reg = registry_with((IntegrationNarrator(),))
    blocks = reg.collect_blocks(
        {"integration_agent": {"status": "done", "findings": {"mofa": mofa_structured_only}}}
    )
    mofa_block = next(b for b in blocks if b.id == "multiome.mofa")
    assert mofa_block.confidence == "medium"
    warn_caveats = [c for c in mofa_block.caveats if c.severity == "warning"]
    assert any(
        "1:1" in c.text and "confounded" in c.text.lower()
        for c in warn_caveats
    ), f"Expected structured batch_condition_confounded caveat, got {warn_caveats}"
    sig_fields = {s["field"] for s in mofa_block.metadata.get("degradation_signals", [])}
    assert "batch_condition_confounded" in sig_fields


def test_d7_021_integration_agent_mofa_resolves_design_main_factor_and_batch_key(tmp_path: Path):
    """D7-021: DesignAgent._build_design stores condition column under main_factor
    and blocking covariate under batch_covariate (never top-level condition_col).
    IntegrationAgent._run_mofa must resolve condition_col (and fallback batch_key)
    via resolve_design_contract."""
    pytest.importorskip("litellm")
    from aria.agents.design_agent import DesignAgent
    from aria.agents.integration_agent import IntegrationAgent
    from aria.memory.memory import ARIAMemory

    da = DesignAgent(ARIAMemory(":memory:"), llm=None)
    da._organism = "Homo sapiens"
    da._genome = "hg38"
    da._confirmed_groups = {"armA": ["s1", "s2"], "armB": ["s3", "s4"]}
    da._main_factor = "perturbation_arm"
    da._batch_covariate = "donor"
    da._inferred_design = {"source": "sample_sheet"}
    design = da._build_design()
    assert design["main_factor"] == "perturbation_arm"
    assert "condition_col" not in design

    captured: dict = {}

    class _CapEnv:
        def run_in_stack(self, stack, script_path, params):
            captured["params"] = dict(params)
            return {"status": "skipped", "reason": "captured"}

    agent = IntegrationAgent(ARIAMemory(":memory:"), llm=None)
    agent.env = _CapEnv()
    agent._advise_mofa_factors = lambda exp_id, intent, mods: 5
    fake_h5mu = tmp_path / "paired.h5mu"
    fake_h5mu.write_bytes(b"")
    agent._run_mofa(
        "exp_mofa",
        {"design": design, "h5mu_path": str(fake_h5mu), "mofa_output_dir": str(tmp_path / "out")},
        intent={},
        mods={"scRNA": [str(fake_h5mu)], "scATAC": [str(fake_h5mu)]},
    )
    assert captured["params"].get("condition_col") == "perturbation_arm", captured["params"]
    assert captured["params"].get("batch_key") == "donor", captured["params"]

