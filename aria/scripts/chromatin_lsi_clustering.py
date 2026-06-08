"""
ARIA Chromatin LSI + Clustering Script (v4.6 scATAC, step 2)
-----------------------------------------------------------
Data-only dimensionality reduction and clustering for single-cell ATAC peak
matrices. The pipeline is the standard scATAC topic-model spine:

  1. TF-IDF on the peak x cell count matrix (signac-style log-TF-log-IDF);
  2. SVD / LSI (truncated SVD on the sparse TF-IDF matrix);
  3. drop LSI component(s) correlated with sequencing depth (the first LSI
     component on ATAC typically captures total fragment count, not biology);
  4. neighbors / UMAP / Leiden on the depth-corrected LSI embedding.

Executed inside aria-chromatin-env by EnvironmentManager.

Honesty contract (ADR-002 / ADR-011): this is a STRUCTURAL analysis only. It
produces clusters and an embedding from the accessibility matrix; it assigns NO
biological identity to clusters and infers nothing from peak/feature names. If
the chromatin stack (scanpy + scipy + sklearn + mudata/anndata) is unavailable,
it returns a structured MissingDependency error instead of fabricating an
embedding. Nothing that is not genuinely computed is reported.

Resume-by-file-validity: a prior `lsi_clustered.h5ad` with a matching parameter
signature and a readable file is reused instead of recomputed.

Input params:
    data_path:          str   — ATAC peak matrix: a `.h5mu` (ATAC modality is
                                 selected automatically) or a peak `.h5ad`.
    n_components:       int   (optional) — LSI/SVD components (default: 50).
    n_neighbors:        int   (optional) — kNN neighbors (default: 15).
    resolution:         float (optional) — Leiden resolution (default: 1.0).
    depth_corr_cutoff:  float (optional) — |Pearson r| with log depth above which
                                 an LSI component is dropped (default: 0.9).
    max_cells:          int   (optional) — random sketch size for very large
                                 matrices (default: 200000).
    seed:               int   (optional) — random_state (default: 0).
    output_dir:         str   (optional) — where to write the output .h5ad.

Output:
    {
      "status":               "success",
      "input_kind":           "h5mu" | "h5ad",
      "atac_modality":        str | None,
      "n_cells_total":        int,
      "n_cells_used":         int,
      "n_peaks":              int,
      "n_components_computed": int,
      "depth_correlations":   [float],   — |r| of each LSI comp vs log depth
      "dropped_components":   [int],     — LSI comps removed as depth-driven
      "n_components_used":    int,
      "rep_used":             "X_lsi",
      "n_clusters":           int,
      "cluster_sizes":        {cluster_id: n_cells},
      "resolution":           float,
      "seed":                 int,
      "sketch_used":          bool,
      "output_path":          str,       — path to lsi_clustered .h5ad
      "warnings":             [str],
      "cache_version":        int,
      "cache_params":         dict,
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script, mocks_allowed

CACHE_VERSION = 1


def _cache_params(params: dict) -> dict:
    return {
        "n_components": int(params.get("n_components", 50)),
        "n_neighbors": int(params.get("n_neighbors", 15)),
        "resolution": float(params.get("resolution", 1.0)),
        "depth_corr_cutoff": float(params.get("depth_corr_cutoff", 0.9)),
        "max_cells": int(params.get("max_cells", 200_000)),
        "seed": int(params.get("seed", 0)),
    }


def _cache_matches(cached: dict, expected: dict) -> bool:
    return (cached.get("cache_version") == CACHE_VERSION
            and cached.get("cache_params") == expected)


def _load_atac_adata(data_path: str):
    """
    Load the ATAC peak matrix as an AnnData. Routes a `.h5mu` through the real
    MuData reader (ATAC modality auto-selected) and a `.h5ad` through the safe
    reader. Returns (adata, input_kind, atac_modality) on success or a
    structured error dict on failure.
    """
    from pathlib import Path

    p = Path(data_path)
    if not p.exists():
        return {"status": "error", "error_type": "FileNotFound",
                "details": f"data_path does not exist: {data_path}"}

    suffix = p.suffix.lower()
    if suffix == ".h5mu":
        from aria.utils.mudata_io import read_h5mu_atac
        res = read_h5mu_atac(str(p))
        if res.get("status") != "success":
            if res.get("error_type") == "MuDataToolingMissing":
                return {"status": "error", "error_type": "MissingDependency",
                        "details": (res["details"] + " Install the chromatin "
                                    "stack (aria-chromatin-env) to read .h5mu.")}
            return res
        return res["adata"], "h5mu", res.get("modality")

    if suffix in (".h5ad", ".h5"):
        from aria.utils.safe_h5ad import h5ad_is_readable, read_h5ad
        if not h5ad_is_readable(str(p)):
            return {"status": "error", "error_type": "UnreadableInput",
                    "details": f"AnnData file is not readable: {data_path}"}
        return read_h5ad(str(p)), "h5ad", None

    return {"status": "error", "error_type": "UnsupportedInput",
            "details": (f"chromatin LSI clustering expects a .h5mu or peak "
                        f".h5ad, got: {p.name}")}


def _tfidf(counts):
    """
    Signac-style log-TF-log-IDF on a sparse peak x cell count matrix laid out as
    (cells x peaks). Returns a sparse CSR TF-IDF matrix. Pure linear algebra;
    no fabrication. Empty cells/peaks are guarded against division by zero.
    """
    import numpy as np
    import scipy.sparse as sp

    X = sp.csr_matrix(counts).astype(np.float64)
    n_cells = X.shape[0]

    cell_totals = np.asarray(X.sum(axis=1)).ravel()       # fragments per cell
    peak_totals = np.asarray(X.sum(axis=0)).ravel()       # cells per peak

    # Term frequency: per-cell L1 normalization.
    inv_cell = np.zeros_like(cell_totals)
    nz = cell_totals > 0
    inv_cell[nz] = 1.0 / cell_totals[nz]
    tf = sp.diags(inv_cell) @ X

    # Inverse document frequency: log(1 + n_cells / peak_total).
    idf = np.log1p(n_cells / np.maximum(peak_totals, 1.0))
    tfidf = tf @ sp.diags(idf)

    # log1p scaling (scale factor 1e4, the signac default).
    tfidf = tfidf.tocsr()
    tfidf.data = np.log1p(tfidf.data * 1e4)
    return tfidf, cell_totals


def chromatin_lsi_clustering(params: dict) -> dict:
    import json
    from pathlib import Path

    data_path = params["data_path"]
    cache_params = _cache_params(params)
    n_components = cache_params["n_components"]
    n_neighbors = cache_params["n_neighbors"]
    resolution = cache_params["resolution"]
    depth_corr_cutoff = cache_params["depth_corr_cutoff"]
    max_cells = cache_params["max_cells"]
    seed = cache_params["seed"]
    output_dir = params.get("output_dir")
    allow_mocks = mocks_allowed(params)

    input_path = Path(data_path)
    out_dir = Path(output_dir) if output_dir else input_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / "lsi_clustered.h5ad"
    summary_path = output_path.with_suffix(".summary.json")

    # ── Resume by file validity ──────────────────────────────────────────────
    if (output_path.exists() and input_path.exists()
            and output_path.stat().st_mtime >= input_path.stat().st_mtime
            and summary_path.exists()):
        try:
            with open(summary_path) as f:
                cached = json.load(f)
            if _cache_matches(cached, cache_params):
                from aria.utils.safe_h5ad import h5ad_is_readable
                if h5ad_is_readable(str(output_path)):
                    cached["resumed"] = True
                    cached["warnings"] = cached.get("warnings", []) + [
                        f"[resume] LSI clustering skipped "
                        f"(valid {output_path.name} exists)"
                    ]
                    return cached
        except Exception:
            pass

    # ── Dependencies ─────────────────────────────────────────────────────────
    try:
        import numpy as np
        import scipy.sparse as sp  # noqa: F401  (used in _tfidf)
        import scanpy as sc
        from sklearn.decomposition import TruncatedSVD
    except ImportError as e:
        if not allow_mocks:
            return {
                "status": "error",
                "error_type": "MissingDependency",
                "details": (
                    "chromatin LSI clustering requires scanpy + scipy + "
                    "scikit-learn (aria-chromatin-env). Install the chromatin "
                    f"stack or pass allow_mock=true. ImportError: {e}"
                ),
            }
        return {
            "status": "success",
            "mocked": True,
            "warnings": [f"Mock LSI clustering (deps missing: {e})."],
            "n_clusters": 0,
            "cluster_sizes": {},
        }

    loaded = _load_atac_adata(data_path)
    if isinstance(loaded, dict):   # structured error
        return loaded
    adata, input_kind, atac_modality = loaded

    warnings: list[str] = []
    n_cells_total = int(adata.n_obs)
    n_peaks = int(adata.n_vars)

    if n_cells_total < 2 or n_peaks < 2:
        return {
            "status": "error",
            "error_type": "InsufficientData",
            "details": (f"peak matrix too small for LSI: {n_cells_total} cells "
                        f"x {n_peaks} peaks."),
        }

    # ── Optional random sketch for very large matrices ──────────────────────
    sketch_used = False
    if n_cells_total > max_cells:
        rng = np.random.default_rng(seed)
        keep = np.sort(rng.choice(n_cells_total, size=max_cells, replace=False))
        adata = adata[keep, :].copy()
        sketch_used = True
        warnings.append(
            f"Sketched {max_cells} of {n_cells_total} cells (seed={seed}) "
            f"before LSI; cluster assignments cover the sketch only."
        )

    # Cap components below the matrix rank.
    max_comp = min(n_components, adata.n_obs - 1, adata.n_vars - 1)
    if max_comp < n_components:
        warnings.append(
            f"Reduced n_components {n_components} -> {max_comp} "
            f"(bounded by matrix dimensions)."
        )
    n_components = int(max_comp)

    # ── TF-IDF → SVD/LSI ─────────────────────────────────────────────────────
    tfidf, cell_totals = _tfidf(adata.X)
    svd = TruncatedSVD(n_components=n_components, random_state=seed)
    lsi = svd.fit_transform(tfidf)   # (n_cells, n_components)

    # ── Drop depth-associated component(s) ──────────────────────────────────
    # The first LSI component on ATAC typically tracks sequencing depth rather
    # than biology. We measure |Pearson r| between each component and log10 total
    # fragments per cell and drop components above the cutoff (signac practice).
    log_depth = np.log10(np.maximum(cell_totals, 1.0))
    depth_corr = []
    for k in range(lsi.shape[1]):
        comp = lsi[:, k]
        if np.std(comp) == 0 or np.std(log_depth) == 0:
            depth_corr.append(0.0)
        else:
            r = np.corrcoef(comp, log_depth)[0, 1]
            depth_corr.append(float(abs(r)) if np.isfinite(r) else 0.0)

    dropped = [k for k, r in enumerate(depth_corr) if r >= depth_corr_cutoff]
    keep_cols = [k for k in range(lsi.shape[1]) if k not in dropped]
    if not keep_cols:
        # Never drop every component: keep all and warn instead of returning an
        # empty embedding.
        warnings.append(
            f"All {lsi.shape[1]} LSI components correlate with depth at "
            f"|r|>={depth_corr_cutoff}; kept all components instead of dropping."
        )
        keep_cols = list(range(lsi.shape[1]))
        dropped = []
    elif dropped:
        warnings.append(
            "Dropped depth-associated LSI component(s) "
            f"{dropped} (|r| vs log depth >= {depth_corr_cutoff})."
        )

    adata.obsm["X_lsi"] = lsi[:, keep_cols]

    # ── Neighbors / UMAP / Leiden on the depth-corrected LSI rep ────────────
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_lsi",
                    random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph",
                 n_iterations=2, directed=False, random_state=seed)

    n_clusters = int(adata.obs["leiden"].nunique())
    cluster_sizes = {
        str(c): int((adata.obs["leiden"] == c).sum())
        for c in adata.obs["leiden"].unique()
    }

    adata.uns["aria_n_cells_total"] = n_cells_total
    adata.write_h5ad(str(output_path))

    result = {
        "status": "success",
        "input_kind": input_kind,
        "atac_modality": atac_modality,
        "n_cells_total": n_cells_total,
        "n_cells_used": int(adata.n_obs),
        "n_peaks": n_peaks,
        "n_components_computed": int(lsi.shape[1]),
        "depth_correlations": [round(r, 4) for r in depth_corr],
        "dropped_components": dropped,
        "n_components_used": len(keep_cols),
        "rep_used": "X_lsi",
        "n_clusters": n_clusters,
        "cluster_sizes": cluster_sizes,
        "resolution": resolution,
        "seed": seed,
        "sketch_used": sketch_used,
        "output_path": str(output_path),
        "warnings": warnings,
        "cache_version": CACHE_VERSION,
        "cache_params": cache_params,
    }
    try:
        with open(summary_path, "w") as f:
            json.dump(result, f, indent=2)
    except Exception:
        pass
    return result


if __name__ == "__main__":
    run_script(chromatin_lsi_clustering)
