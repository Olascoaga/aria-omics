"""Permanent regression tests for D1 silent-failure fixes in integration compute kernels."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


def test_d1_013_mofa_must_not_silently_drop_corrupt_third_modality(tmp_path):
    """D1-013: Corrupt 3rd modality in MOFA+ must fail closed (must run in aria-integration-env)."""
    pytest.importorskip("mofapy2")
    ad = pytest.importorskip("anndata")
    pytest.importorskip("scanpy")
    import scipy.sparse as sp
    from aria.scripts.integration_mofa import integration_mofa

    rng = np.random.default_rng(1)
    bcs = [f"cell_{i}" for i in range(110)]
    rna_p = tmp_path / "rna.h5ad"
    atac_p = tmp_path / "atac.h5ad"
    prot_p = tmp_path / "protein.h5ad"
    ad.AnnData(
        X=sp.csr_matrix(rng.poisson(5, size=(110, 25)).astype(np.int32)),
        obs=pd.DataFrame(index=bcs),
        var=pd.DataFrame(index=[f"G{i}" for i in range(25)]),
    ).write_h5ad(rna_p)
    ad.AnnData(
        X=sp.csr_matrix(rng.poisson(2, size=(110, 25)).astype(np.int32)),
        obs=pd.DataFrame(index=bcs),
        var=pd.DataFrame(index=[f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(25)]),
    ).write_h5ad(atac_p)
    prot_p.write_bytes(b"CORRUPT_H5AD_FILE")

    res = integration_mofa({
        "modalities": {
            "rna": [str(rna_p)],
            "atac": [str(atac_p)],
            "protein": [str(prot_p)],
        },
        "n_factors": 3,
        "max_iter": 10,
        "output_dir": str(tmp_path / "out"),
    })
    assert res["status"] == "error" and res.get("error_type") != "MissingDependency"
