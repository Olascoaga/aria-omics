"""Regression tests for D7 layer-contract fixes in scRNA."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest


def test_d7_017_warning_router_surfaces_scrna_partial_pairing_dropped_donors(tmp_path):
    """D7-017: When a paired scRNA cohort has partial pairing, rna_pseudobulk_de
    drops unpaired donors and emits pairing_status='partial' and
    excluded_unpaired_samples=[...]. Through NarrativeRegistry.collect_blocks,
    warning_router must attach a structured warning caveat naming the dropped donor."""
    pytest.importorskip("litellm")
    pytest.importorskip("pydantic")
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator
    from aria.agents.narrative.registry import registry_with
    from aria.utils.environment_manager import env_manager

    rng = np.random.default_rng(17)
    donors_conds = [
        ("D1", "ctrl"), ("D1", "treat"),
        ("D2", "ctrl"), ("D2", "treat"),
        ("D3", "ctrl"), ("D3", "treat"),
        ("D4", "ctrl"),
    ]
    n_per = 25
    n_genes = 40
    xs, obs_rows = [], []
    for donor, cond in donors_conds:
        counts = rng.poisson(20, size=(n_per, n_genes))
        if cond == "treat":
            counts[:, :5] = rng.poisson(120, size=(n_per, 5))
        xs.append(counts)
        for j in range(n_per):
            obs_rows.append({
                "cell_id": f"{donor}_{cond}_{j}",
                "donor": donor,
                "condition": cond,
                "cell_type": "T_cell",
            })
    obs = pd.DataFrame(obs_rows).set_index("cell_id")
    adata = ad.AnnData(
        X=np.vstack(xs).astype(np.float32),
        obs=obs,
        var=pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes)]),
    )
    h5ad_path = tmp_path / "pb.h5ad"
    adata.write_h5ad(str(h5ad_path))

    pb_res = env_manager.run_in_stack(
        "rna",
        "aria/scripts/rna_pseudobulk_de.py",
        {
            "data_path": str(h5ad_path),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["treat", "ctrl"]],
            "min_cells_per_sample": 5,
            "min_replicates_per_condition": 3,
            "output_dir": str(tmp_path / "pb_out"),
            "padj_max": 0.05,
            "lfc_min": 0.5,
        },
    )
    assert pb_res["status"] == "success", pb_res
    comp = pb_res["per_group"]["T_cell"]["per_comparison"]["treat_vs_ctrl"]
    assert comp["pairing_status"] == "partial"
    assert len(comp["excluded_unpaired_samples"]) == 1

    reg = registry_with((ScrnaNarrator(),))
    blocks = reg.collect_blocks(
        {"scrna_agent": {"status": "done", "findings": {"pseudobulk_de": pb_res}}}
    )
    pb_blocks = [b for b in blocks if b.analysis == "pseudobulk_de"]
    assert pb_blocks
    pb_block = pb_blocks[0]
    warn_caveats = [c for c in pb_block.caveats if c.severity == "warning"]
    assert any("unpaired" in c.text.lower() and "d4" in c.text.lower() for c in warn_caveats), (
        f"Dropped unpaired donor(s) {comp['excluded_unpaired_samples']} not surfaced: {warn_caveats}"
    )
    sig_fields = {s["field"] for s in pb_block.metadata.get("degradation_signals", [])}
    assert "pairing_status" in sig_fields
