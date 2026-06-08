"""v4.6 scATAC step 3 — chromatin differential accessibility guards.

Covers `aria/scripts/chromatin_diffacc.py`:

- IPC contract registration;
- structured error routing (FileNotFound);
- honest gating of the pseudobulk lane (no condition / no replicate / no
  explicit comparison => not-run with a concrete reason, never fabricated);
- per-cluster Wilcoxon marker-peak accessibility on a synthetic peak matrix;
- pseudobulk DA through the SHARED validated DESeq2 core (pydeseq2-gated).

Synthetic labels only (ADR-011): no biological gene/peak names.
"""

import pytest


# ── Contract + light error routing ───────────────────────────────────────────

def test_diffacc_has_ipc_contract():
    from aria.utils.script_contracts import SCRIPT_CONTRACTS
    key = "aria/scripts/chromatin_diffacc.py"
    assert key in SCRIPT_CONTRACTS
    assert "data_path" in {f.name for f in SCRIPT_CONTRACTS[key].inputs}


def test_diffacc_file_not_found(tmp_path):
    from aria.scripts.chromatin_diffacc import chromatin_diffacc
    res = chromatin_diffacc({"data_path": str(tmp_path / "missing.h5ad")})
    assert res["status"] == "error" and res["error_type"] == "FileNotFound"


# ── Pseudobulk lane gating (honest skips) ────────────────────────────────────

def _clustered_adata(tmp_path, with_condition=False, with_replicate=False,
                     seed=0):
    ad = pytest.importorskip("anndata")
    import numpy as np

    rng = np.random.default_rng(seed)
    n_per, n_peaks = 30, 120
    blocks, leiden = [], []
    for grp in range(2):
        base = np.zeros((n_per, n_peaks))
        lo, hi = (0, 60) if grp == 0 else (60, 120)
        base[:, lo:hi] = 1.0
        blocks.append(rng.poisson(base * 3.0 + 0.2))
        leiden += [str(grp)] * n_per
    X = np.vstack(blocks).astype("float32")
    adata = ad.AnnData(X)
    adata.obs["leiden"] = leiden
    if with_condition:
        # two conditions, balanced across both clusters
        adata.obs["condition"] = (["A", "B"] * (adata.n_obs // 2 + 1))[:adata.n_obs]
    if with_replicate:
        reps = []
        for i in range(adata.n_obs):
            cond = adata.obs["condition"].iloc[i] if with_condition else "X"
            reps.append(f"{cond}_rep{i % 3}")
        adata.obs["replicate"] = reps
    p = tmp_path / "lsi_clustered.h5ad"
    adata.write_h5ad(str(p))
    return p


def test_pseudobulk_skips_without_condition(tmp_path):
    pytest.importorskip("scanpy")
    from aria.scripts.chromatin_diffacc import chromatin_diffacc
    p = _clustered_adata(tmp_path)
    res = chromatin_diffacc({"data_path": str(p), "output_dir": str(tmp_path)})
    assert res["status"] == "success"
    assert res["pseudobulk"]["ran"] is False
    assert "condition" in res["pseudobulk"]["reason"]


def test_pseudobulk_skips_without_replicate(tmp_path):
    pytest.importorskip("scanpy")
    from aria.scripts.chromatin_diffacc import chromatin_diffacc
    p = _clustered_adata(tmp_path, with_condition=True)
    res = chromatin_diffacc({
        "data_path": str(p), "condition_col": "condition",
        "comparisons": [{"test": "A", "reference": "B"}],
        "output_dir": str(tmp_path),
    })
    assert res["pseudobulk"]["ran"] is False
    assert "replicate" in res["pseudobulk"]["reason"]


def test_pseudobulk_requires_explicit_comparison(tmp_path):
    pytest.importorskip("scanpy")
    from aria.scripts.chromatin_diffacc import chromatin_diffacc
    p = _clustered_adata(tmp_path, with_condition=True, with_replicate=True)
    res = chromatin_diffacc({
        "data_path": str(p), "condition_col": "condition",
        "replicate_col": "replicate", "output_dir": str(tmp_path),
    })
    # P0-5: no DE without an explicit reference contrast.
    assert res["pseudobulk"]["ran"] is False
    assert "explicit comparison" in res["pseudobulk"]["reason"]


# ── Per-cluster accessibility (scanpy-gated) ──────────────────────────────────

def test_per_cluster_accessibility_runs(tmp_path):
    pytest.importorskip("scanpy")
    from aria.scripts.chromatin_diffacc import chromatin_diffacc
    p = _clustered_adata(tmp_path)
    res = chromatin_diffacc({"data_path": str(p), "output_dir": str(tmp_path)})
    assert res["status"] == "success"
    assert res["n_clusters"] == 2
    pc = res["per_cluster"]
    assert pc["ran"] is True
    # planted block structure should yield differential peaks for the clusters
    assert pc["n_da_total"] > 0
    # records are descriptive (peak coordinate + stats + accessibility fraction)
    for recs in pc["da_peaks_by_cluster"].values():
        for r in recs:
            assert {"peak", "log2fc", "padj", "pct_in", "pct_out"} <= set(r)
            assert 0.0 <= r["pct_in"] <= 1.0


def test_per_cluster_skips_single_cluster(tmp_path):
    pytest.importorskip("scanpy")
    ad = pytest.importorskip("anndata")
    import numpy as np
    from aria.scripts.chromatin_diffacc import chromatin_diffacc

    adata = ad.AnnData(np.abs(np.ones((20, 30))).astype("float32"))
    adata.obs["leiden"] = ["0"] * 20
    p = tmp_path / "lsi_clustered.h5ad"
    adata.write_h5ad(str(p))
    res = chromatin_diffacc({"data_path": str(p), "output_dir": str(tmp_path)})
    assert res["status"] == "success"
    assert res["per_cluster"]["ran"] is False


# ── Pseudobulk DA through the shared DESeq2 core (pydeseq2-gated) ──────────────

def test_pseudobulk_da_runs_with_replicates(tmp_path):
    pytest.importorskip("scanpy")
    pytest.importorskip("pydeseq2")
    ad = pytest.importorskip("anndata")
    import numpy as np
    from aria.scripts.chromatin_diffacc import chromatin_diffacc

    rng = np.random.default_rng(0)
    n_peaks = 200
    rows, cond, rep = [], [], []
    # 2 conditions x 3 replicates x 25 cells; peaks 0..49 up in condition B.
    for c, cname in enumerate(["A", "B"]):
        for r in range(3):
            for _ in range(25):
                base = np.full(n_peaks, 2.0)
                if cname == "B":
                    base[:50] += 8.0
                rows.append(rng.poisson(base))
                cond.append(cname)
                rep.append(f"{cname}_rep{r}")
    X = np.vstack(rows).astype("float32")
    adata = ad.AnnData(X)
    adata.obs["leiden"] = ["0"] * (adata.n_obs // 2) + ["1"] * (adata.n_obs - adata.n_obs // 2)
    adata.obs["condition"] = cond
    adata.obs["replicate"] = rep
    p = tmp_path / "lsi_clustered.h5ad"
    adata.write_h5ad(str(p))

    res = chromatin_diffacc({
        "data_path": str(p), "condition_col": "condition",
        "replicate_col": "replicate",
        "comparisons": [{"test": "B", "reference": "A"}],
        "min_replicates_per_condition": 3,
        "output_dir": str(tmp_path),
    })
    assert res["status"] == "success"
    pb = res["pseudobulk"]
    assert pb["ran"] is True
    comp = pb["comparisons"][0]
    assert comp["status"] == "success"
    assert comp["n_sig"] > 0
    assert comp["n_up"] > 0   # the planted peaks are up in B vs A
