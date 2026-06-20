"""
ARIA Chromatin Differential Accessibility Script (v4.6 scATAC, step 3)
---------------------------------------------------------------------
Differential accessibility (DA) on a clustered scATAC peak matrix, in two
honest lanes:

  1. Per-cluster accessibility (descriptive): Wilcoxon rank-sum marker peaks,
     cluster vs rest, with the fraction of cells accessible in/out of the
     cluster (pct_in / pct_out). This is the analog of per-cluster marker
     genes — it ranks the peaks that distinguish each cluster. It assigns NO
     biological identity to clusters or peaks.

  2. Per-condition pseudobulk DA (inferential, replicate-gated): when the obs
     carries a condition column AND a biological-replicate column AND the caller
     supplies an explicit comparison (numerator vs reference), peak counts are
     summed per biological replicate into pseudobulk samples and the SHARED,
     validated DESeq2 core (`rna_bulk_de._run_deseq2`) is run on the peaks. This
     is the only statistically defensible cross-condition DA: cells from one
     replicate are not independent observations.

Executed inside aria-chromatin-env by EnvironmentManager.

Honesty contract (ADR-002 / ADR-005 / ADR-011):
  - No peak/feature identity is inferred from names; peaks are coordinates.
  - The pseudobulk lane runs ONLY with an explicit comparison (P0-5: no DE
    without an explicit reference) and ONLY with enough replicates per group
    (`min_replicates_per_condition`, default from AnalysisThresholds). Otherwise
    it is reported as not-run with a concrete reason — never fabricated.
  - Significance thresholds come from the confirmed experiment context
    (`AnalysisThresholds`, P0-7), not loose literals.
  - Missing scientific tooling returns a structured error, not a placeholder.

Input params:
    data_path:        str   — clustered scATAC .h5ad (peak matrix + cluster obs).
    groupby:          str   (optional) — cluster obs column (default: "leiden").
    condition_col:    str   (optional) — obs column with the condition/group.
    replicate_col:    str   (optional) — obs column identifying biological
                              replicates (donor/sample); REQUIRED for pseudobulk.
    comparisons:      list  (optional) — [{"test": lvl, "reference": lvl}, ...]
                              explicit contrasts for the pseudobulk lane.
    padj_max:         float (optional) — adjusted p-value cutoff.
    lfc_min:          float (optional) — minimum |log2FC|.
    top_n:            int   (optional) — top peaks per cluster/contrast (def 20).
    min_cells_per_pseudosample:   int (optional, default 10).
    min_replicates_per_condition: int (optional).
    exp_context:      dict  (optional) — confirmed thresholds (global_padj/lfc).
    output_dir:       str   (optional).

Output:
    {
      "status":              "success" | "error",
      "groupby":             str,
      "n_clusters":          int,
      "padj_max":            float,
      "lfc_min":             float,
      "per_cluster": {
        "ran":               bool,
        "reason":            str | None,
        "n_da_total":        int,
        "n_da_by_cluster":   {cluster_id: int},
        "da_peaks_by_cluster": {cluster_id: [{peak, log2fc, padj,
                                              pct_in, pct_out}, ...]},
        "output_csv":        str | None,
      },
      "pseudobulk": {
        "ran":               bool,
        "reason":            str | None,
        "comparisons":       [ {test, reference, status, n_sig, n_up, n_down,
                                top_peaks, low_power_warning, ...} ],
        "output_csv":        str | None,
      },
      "warnings":            [str],
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script, mocks_allowed


def _binary_pct(X, mask, col_idx):
    """Fraction of cells in `mask` with a nonzero count at peak `col_idx`."""
    import numpy as np
    col = X[mask, col_idx]
    if hasattr(col, "toarray"):
        col = col.toarray().ravel()
    return float((np.asarray(col).ravel() > 0).mean())


def _per_cluster_accessibility(adata, groupby, padj_max, lfc_min, top_n,
                               warnings):
    """Wilcoxon marker peaks per cluster (cluster vs rest). Descriptive only."""
    import scanpy as sc
    import numpy as np
    import pandas as pd

    n_clusters = int(adata.obs[groupby].astype(str).nunique())
    if n_clusters < 2:
        return {
            "ran": False,
            "reason": f"only {n_clusters} cluster(s) — per-cluster DA needs >=2",
            "n_da_total": 0,
            "n_da_by_cluster": {},
            "da_peaks_by_cluster": {},
            "output_csv": None,
        }

    # Wilcoxon is computed on a log-normalized COPY so a few high-count peaks do
    # not dominate the ranking; pct_in/pct_out are computed on the raw matrix
    # (a peak is "accessible" in a cell when its count is > 0). The normalized
    # copy is not persisted and makes no biological claim.
    work = adata.copy()
    raw_X = adata.X
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    work.obs[groupby] = work.obs[groupby].astype(str)

    sc.tl.rank_genes_groups(
        work, groupby=groupby, method="wilcoxon", key_added="aria_da",
    )
    rgg = work.uns["aria_da"]

    var_names = list(adata.var_names)
    var_idx = {p: i for i, p in enumerate(var_names)}
    membership = {
        str(c): np.asarray(adata.obs[groupby].astype(str) == str(c))
        for c in adata.obs[groupby].astype(str).unique()
    }

    da_by_cluster: dict = {}
    n_by_cluster: dict = {}
    rows = []
    n_total = 0

    for cl in rgg["names"].dtype.names:
        cl_s = str(cl)
        mask_in = membership.get(cl_s)
        mask_out = ~mask_in if mask_in is not None else None
        names = rgg["names"][cl]
        lfc = rgg["logfoldchanges"][cl]
        padj = rgg["pvals_adj"][cl]

        sig = []
        for peak, l, p in zip(names, lfc, padj):
            peak = str(peak)
            try:
                p_val = float(p)
                l_val = float(l)
            except (TypeError, ValueError):
                continue
            if (np.isnan(p_val) or np.isnan(l_val)
                    or peak in ("nan", "") or peak not in var_idx):
                continue
            if p_val >= padj_max or abs(l_val) < lfc_min:
                continue
            pi = var_idx[peak]
            rec = {
                "cluster": cl_s,
                "peak": peak,
                "log2fc": round(l_val, 3),
                "padj": round(p_val, 6),
                "pct_in": round(_binary_pct(raw_X, mask_in, pi), 3)
                if mask_in is not None else None,
                "pct_out": round(_binary_pct(raw_X, mask_out, pi), 3)
                if mask_out is not None else None,
            }
            sig.append(rec)
            rows.append(rec)

        n_by_cluster[cl_s] = len(sig)
        n_total += len(sig)
        da_by_cluster[cl_s] = sig[:top_n]

    return {
        "ran": True,
        "reason": None,
        "n_da_total": n_total,
        "n_da_by_cluster": n_by_cluster,
        "da_peaks_by_cluster": da_by_cluster,
        "_rows": rows,
        "n_clusters": n_clusters,
    }


def _gate_da_by_convergence(results_df, padj_max):
    """T15 / scATAC P0 W0.4: a peak is reported SIGNIFICANT only if its padj clears
    the threshold AND its apeGLM LFC fit converged.

    On non-converged fits pydeseq2 can report a near-zero shrunk log2FoldChange with
    a still-significant Wald p (the contradiction ADR-042/043 disclosed). When the
    convergence flag is available (`lfc_converged` column, exposed opt-in by
    `_run_deseq2`), those peaks are excluded from the significant set and counted.
    When it is NOT available (older pydeseq2 / mock), NO gating is applied — we do
    not fabricate a convergence verdict. Returns ``(sig_df, info)`` where ``info``
    carries ``convergence_available`` + ``n_nonconverged_excluded``.
    """
    padj_sig = results_df[results_df["padj"] < padj_max]
    if "lfc_converged" in results_df.columns:
        converged = padj_sig["lfc_converged"].fillna(False).astype(bool)
        sig = padj_sig[converged]
        n_excluded = int((~converged).sum())
        available = True
    else:
        sig = padj_sig
        n_excluded = 0
        available = False
    return sig.sort_values("padj"), {
        "convergence_available": available,
        "n_nonconverged_excluded": n_excluded,
    }


def _pseudobulk_da(adata, condition_col, replicate_col, comparisons,
                   padj_max, lfc_min, top_n, min_cells, min_reps,
                   allow_mock, warnings, n_cpus=None):
    """Per-condition pseudobulk DA on peaks via the shared validated DESeq2 core.

    Returns a structured result. The lane is gated: it requires the condition +
    replicate columns and explicit comparisons, and each contrast needs enough
    replicates per group. Anything missing is reported as not-run with a reason,
    never fabricated.
    """
    import numpy as np
    import pandas as pd
    from scipy import sparse

    if not condition_col or condition_col not in adata.obs.columns:
        return {"ran": False,
                "reason": f"no usable condition column "
                          f"({condition_col!r} not in obs)",
                "comparisons": [], "output_csv": None}
    if not replicate_col or replicate_col not in adata.obs.columns:
        return {"ran": False,
                "reason": (f"no biological-replicate column "
                           f"({replicate_col!r} not in obs); pseudobulk DA "
                           f"needs replicates, single cells are not replicates"),
                "comparisons": [], "output_csv": None}
    if not comparisons:
        return {"ran": False,
                "reason": ("no explicit comparison supplied; ARIA does not "
                           "invent a reference contrast (P0-5)"),
                "comparisons": [], "output_csv": None}

    # Import the SHARED DESeq2 core (no module-level heavy deps; safe to import
    # in the chromatin env). This keeps one validated DE engine across modalities.
    from aria.scripts.rna_bulk_de import _run_deseq2

    obs = adata.obs[[condition_col, replicate_col]].copy()
    obs[condition_col] = obs[condition_col].astype(str)
    obs[replicate_col] = obs[replicate_col].astype(str)

    # Pseudosample = biological replicate. One replicate can only map to one
    # condition; if a replicate spans conditions, key by replicate × condition.
    rep_cond = obs.groupby(replicate_col)[condition_col].nunique()
    if (rep_cond > 1).any():
        obs["_ps"] = obs[replicate_col] + "::" + obs[condition_col]
    else:
        obs["_ps"] = obs[replicate_col]

    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    X = X.tocsr()

    # Aggregate raw peak counts per pseudosample, dropping pseudosamples with
    # too few cells (unstable aggregates).
    ps_rows = []
    ps_meta = []
    for ps, idx in obs.groupby("_ps").groups.items():
        positions = obs.index.get_indexer(idx)
        if len(positions) < min_cells:
            warnings.append(
                f"Dropped pseudosample '{ps}' with {len(positions)} cells "
                f"(< min_cells_per_pseudosample={min_cells})."
            )
            continue
        summed = np.asarray(X[positions, :].sum(axis=0)).ravel()
        ps_rows.append(summed)
        cond = obs[condition_col].iloc[positions[0]]
        ps_meta.append({"pseudosample": str(ps), condition_col: str(cond)})

    if len(ps_rows) < 2:
        return {"ran": False,
                "reason": (f"only {len(ps_rows)} usable pseudosample(s) after "
                           f"the min-cells filter; need >=2"),
                "comparisons": [], "output_csv": None}

    counts_mat = np.vstack(ps_rows).astype(float)          # pseudosamples x peaks
    meta_df = pd.DataFrame(ps_meta).set_index("pseudosample")
    peak_names = list(adata.var_names)

    # _run_deseq2 expects counts as features x samples (columns indexed by the
    # metadata index) and transposes internally.
    counts_df = pd.DataFrame(counts_mat.T, index=peak_names, columns=meta_df.index)

    # Drop peaks that are all-zero across the usable pseudosamples (untestable).
    nonzero = counts_df.sum(axis=1) > 0
    dropped_peaks = int((~nonzero).sum())
    counts_df = counts_df.loc[nonzero]
    if dropped_peaks:
        warnings.append(
            f"Excluded {dropped_peaks} peaks with zero pseudobulk counts "
            f"across all replicates from DA testing."
        )

    comp_results = []
    rows_csv = []
    full_rows = []   # W0.3: every tested peak per comparison (for volcano/MA)
    for comp in comparisons:
        test = str(comp.get("test", "")).strip()
        ref = str(comp.get("reference", comp.get("ref", ""))).strip()
        if not test or not ref:
            comp_results.append({"test": test, "reference": ref,
                                 "status": "skipped",
                                 "reason": "incomplete comparison"})
            continue
        levels = set(meta_df[condition_col].unique())
        if test not in levels or ref not in levels:
            comp_results.append({
                "test": test, "reference": ref, "status": "skipped",
                "reason": (f"contrast levels not present in pseudobulk "
                           f"(have {sorted(levels)})")})
            continue

        de, de_warn = _run_deseq2(
            counts_df, meta_df, condition_col, test, ref,
            padj_max, lfc_min, allow_mock=allow_mock,
            min_replicates_per_condition=min_reps,
            n_cpus=n_cpus,
            expose_convergence=True,   # T15 / W0.4
        )
        warnings.extend(de_warn)
        if de.get("status") != "success":
            comp_results.append({
                "test": test, "reference": ref, "status": "error",
                "error_type": de.get("error_type"),
                "reason": de.get("details", "")[:300],
            })
            continue

        res_df = de["results"]
        comp_key = f"{test}_vs_{ref}"
        # T15 / W0.4: gate the significant set by apeGLM LFC convergence, then
        # recompute the DA counts from the gated set (the shared _run_deseq2 n_sig
        # stays convergence-agnostic for bulk RNA; the scATAC lane is stricter).
        sig, conv_info = _gate_da_by_convergence(res_df, padj_max)
        n_excluded = conv_info["n_nonconverged_excluded"]
        if n_excluded:
            warnings.append(
                f"[{comp_key}] excluded {n_excluded} non-converged peak(s) from "
                f"the DA significant set (apeGLM LFC≈0 with a significant Wald p; "
                f"T15 convergence gate)."
            )
        n_up = int((sig["log2FoldChange"] > 0).sum())
        n_down = int((sig["log2FoldChange"] < 0).sum())

        top_peaks = []
        for peak in list(sig.index)[:top_n]:
            r = sig.loc[peak]
            rec = {
                "comparison": comp_key,
                "peak": str(peak),
                "log2fc": round(float(r["log2FoldChange"]), 3),
                "padj": round(float(r["padj"]), 6),
            }
            top_peaks.append(rec)
            rows_csv.append(rec)

        # W0.3: persist the FULL per-peak table (all tested peaks) so the volcano /
        # MA figures have the complete cloud, with the convergence flag carried
        # through. Non-converged peaks are kept here (flagged) but never counted as
        # significant nor plotted as real effects.
        for peak in list(res_df.index):
            r = res_df.loc[peak]
            conv = r.get("lfc_converged") if "lfc_converged" in res_df.columns else None
            base_mean = r.get("baseMean") if "baseMean" in res_df.columns else None
            full_rows.append({
                "comparison": comp_key,
                "peak": str(peak),
                "log2fc": float(r["log2FoldChange"]),
                "padj": float(r["padj"]),
                "base_mean": (float(base_mean) if base_mean is not None
                              and base_mean == base_mean else None),
                "lfc_converged": (bool(conv) if conv is not None and conv == conv
                                  else None),
                "significant": bool(peak in sig.index),
            })

        comp_results.append({
            "test": test, "reference": ref, "status": "success",
            "n_sig": int(len(sig)),
            "n_up": n_up,
            "n_down": n_down,
            "n_replicates": de.get("n_replicates"),
            "low_power_warning": de.get("low_power_warning", False),
            "fitted_design_formula": de.get("fitted_design_formula"),
            "convergence_available": conv_info["convergence_available"],
            "n_nonconverged_excluded": n_excluded,
            "top_peaks": top_peaks,
        })

    return {
        "ran": any(c.get("status") == "success" for c in comp_results),
        "reason": None,
        "n_pseudosamples": int(len(ps_rows)),
        "n_peaks_tested": int(counts_df.shape[0]),
        "comparisons": comp_results,
        "_rows": rows_csv,
        "_full_rows": full_rows,
    }


def chromatin_diffacc(params: dict) -> dict:
    import json
    from pathlib import Path

    data_path = params["data_path"]
    groupby = params.get("groupby", "leiden")
    condition_col = params.get("condition_col")
    replicate_col = params.get("replicate_col")
    comparisons = params.get("comparisons") or []
    top_n = int(params.get("top_n", 20))
    output_dir = params.get("output_dir")
    allow_mock = mocks_allowed(params)

    from aria.utils.thresholds import AnalysisThresholds
    thr = AnalysisThresholds.from_exp_context(
        params.get("exp_context"), modality="scATAC")
    padj_max = float(params.get("padj_max", thr.padj))
    lfc_min = float(params.get("lfc_min", thr.log2fc))
    min_cells = int(params.get("min_cells_per_pseudosample", thr.min_cells))
    min_reps = int(params.get("min_replicates_per_condition", thr.min_replicates))
    n_cpus = params.get("n_cpus")
    n_cpus = int(n_cpus) if n_cpus else None
    # Optional: skip the DESCRIPTIVE per-cluster Wilcoxon lane and go straight to
    # the inferential pseudobulk DA. The per-cluster lane is single-threaded scipy
    # CSR slicing (pct_in/pct_out) that scales pathologically with cluster count ×
    # significant peaks on a large multi-sample pool; it adds no pseudobulk
    # evidence. Default False keeps the production path unchanged.
    skip_per_cluster = bool(params.get("skip_per_cluster", False))

    in_path = Path(data_path)
    if not in_path.exists():
        return {"status": "error", "error_type": "FileNotFound",
                "details": f"data_path does not exist: {data_path}"}
    out_dir = Path(output_dir) if output_dir else in_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        import scanpy  # noqa: F401
        import numpy  # noqa: F401
        import pandas as pd
        from scipy import sparse  # noqa: F401
        from aria.utils.safe_h5ad import h5ad_is_readable, read_h5ad
    except ImportError as e:
        if not allow_mock:
            return {"status": "error", "error_type": "MissingDependency",
                    "details": (
                        "chromatin DA requires scanpy + scipy + pandas "
                        f"(aria-chromatin-env). ImportError: {e}")}
        return {"status": "success", "mocked": True,
                "warnings": [f"Mock DA (deps missing: {e})."],
                "per_cluster": {"ran": False, "reason": "mock"},
                "pseudobulk": {"ran": False, "reason": "mock"}}

    if not h5ad_is_readable(str(in_path)):
        return {"status": "error", "error_type": "UnreadableInput",
                "details": f"clustered AnnData is not readable: {data_path}"}

    adata = read_h5ad(str(in_path))
    if groupby not in adata.obs.columns:
        return {"status": "error", "error_type": "NoGroupBy",
                "details": (f"cluster column '{groupby}' not in obs. "
                            f"Run chromatin_lsi_clustering first. "
                            f"Available: {list(adata.obs.columns)[:15]}")}

    warnings: list = []

    if skip_per_cluster:
        per_cluster = {
            "ran": False,
            "reason": "skipped (skip_per_cluster=True): descriptive lane "
                      "disabled to run the inferential pseudobulk DA directly",
            "n_da_total": 0, "n_da_by_cluster": {}, "da_peaks_by_cluster": {},
        }
        n_clusters = int(adata.obs[groupby].astype(str).nunique())
    else:
        per_cluster = _per_cluster_accessibility(
            adata, groupby, padj_max, lfc_min, top_n, warnings)
        n_clusters = per_cluster.pop("n_clusters", None) or int(
            adata.obs[groupby].astype(str).nunique())

    # Persist per-cluster table. F9: a failed write is disclosed, not swallowed.
    from aria.utils.csv_artifacts import persist_csv
    pc_rows = per_cluster.pop("_rows", [])
    pc_csv = None
    if pc_rows:
        pc_path = str(out_dir / "chromatin_da_per_cluster.csv")
        pc_csv = persist_csv(
            pc_path, "chromatin_da_per_cluster.csv",
            lambda: pd.DataFrame(pc_rows).to_csv(pc_path, index=False), warnings)
    per_cluster["output_csv"] = pc_csv

    pseudobulk = _pseudobulk_da(
        adata, condition_col, replicate_col, comparisons,
        padj_max, lfc_min, top_n, min_cells, min_reps, allow_mock, warnings,
        n_cpus=n_cpus)
    pb_rows = pseudobulk.pop("_rows", [])
    pb_csv = None
    if pb_rows:
        pb_path = str(out_dir / "chromatin_da_pseudobulk.csv")
        pb_csv = persist_csv(
            pb_path, "chromatin_da_pseudobulk.csv",
            lambda: pd.DataFrame(pb_rows).to_csv(pb_path, index=False), warnings)
    pseudobulk["output_csv"] = pb_csv

    # W0.3: persist the FULL per-peak pseudobulk table (all tested peaks, with the
    # convergence flag + significance) so the volcano/MA figures have the complete
    # cloud, not just the top-N significant rows.
    pb_full_rows = pseudobulk.pop("_full_rows", [])
    pb_full_csv = None
    if pb_full_rows:
        pb_full_path = str(out_dir / "chromatin_da_pseudobulk_full.csv")
        pb_full_csv = persist_csv(
            pb_full_path, "chromatin_da_pseudobulk_full.csv",
            lambda: pd.DataFrame(pb_full_rows).to_csv(pb_full_path, index=False),
            warnings)
    pseudobulk["full_results_csv"] = pb_full_csv

    return {
        "status": "success",
        "groupby": groupby,
        "n_clusters": n_clusters,
        "padj_max": padj_max,
        "lfc_min": lfc_min,
        "per_cluster": per_cluster,
        "pseudobulk": pseudobulk,
        "warnings": warnings,
    }


if __name__ == "__main__":
    run_script(chromatin_diffacc)
