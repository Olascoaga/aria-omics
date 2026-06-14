"""
ARIA RNA Trajectory Script
---------------------------
Pseudotime and trajectory inference for scRNA-seq data.
Executed inside aria-rna-env by EnvironmentManager (standalone entry point).

Methods (in order of application):
  1. PAGA         — graph abstraction of cluster connectivity (always)
  2. DPT          — diffusion pseudotime when a defensible root is available
  3. scVelo       — RNA velocity (only if spliced/unspliced layers present)

Input params:
    data_path:       str  — path to clustered .h5ad
    root_cell_type:  str  (optional) — cell type to use as trajectory root
    cell_type_col:   str  (optional) — obs column with cell type labels (default: "cell_type")
    output_dir:      str  (optional)

Output:
    {
      "status":  "success" | "error",
      "paga":    {"top_connections": {"A→B": float, ...}},
      "pseudotime": {
          "computed": bool,
          "pseudotime_by_group": {group: mean_pt, ...},
          "root_used": str | None,
          "root_selection_policy": str,
          "reason": str  (only when not computed),
      },
      "velocity": {"computed": bool, "method": str | "reason": str},
      "groupby":  str,
      "output_path": str,
    }
"""

from __future__ import annotations
import re
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


# T14 (tri-audit 2026-06-14): PAGA "strong" cluster adjacency threshold. `n_strong`
# counts off-diagonal connectivities at or above this value. (top_connections is a
# separate display list — the top 25 edges by connectivity, no floor applied.)
_PAGA_STRONG_CONNECTIVITY = 0.05

# T13 (tri-audit 2026-06-14): progenitor-like cell-type tokens for DPT root
# auto-detection. Matched as WHOLE WORDS (with optional plural) so "stem" does NOT
# fire on "brainstem"; a plain `tok in label` substring test rooted pseudotime on
# unrelated populations.
_PROGENITOR_TOKENS = ("stem", "progenitor", "precursor")
_PROGENITOR_RE = re.compile(
    r"\b(?:" + "|".join(_PROGENITOR_TOKENS) + r")s?\b", re.IGNORECASE
)


def _is_progenitor_label(label: str) -> bool:
    """True iff a cell-type label names a progenitor-like population by a generic
    token as a WHOLE WORD ("neural stem cell", "Radial glia progenitors") — never
    by substring ("brainstem", "stemness")."""
    return bool(_PROGENITOR_RE.search(str(label or "")))


def rna_trajectory(params: dict) -> dict:
    import scanpy as sc
    import numpy as np
    from pathlib import Path
    from aria.utils.safe_h5ad import read_h5ad

    data_path      = params["data_path"]
    root_cell_type = params.get("root_cell_type")
    cell_type_col  = params.get("cell_type_col", "cell_type")
    output_dir     = params.get("output_dir", str(Path(data_path).parent))

    adata = read_h5ad(data_path)

    # Pick grouping column
    if "leiden" in adata.obs.columns:
        groupby = "leiden"
    elif cell_type_col in adata.obs.columns:
        groupby = cell_type_col
    else:
        return {"status": "error", "error_type": "NoClusters",
                "details": "Need 'leiden' or cell_type column for PAGA."}

    # Ensure neighbor graph
    if "neighbors" not in adata.uns:
        rep = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"
        sc.pp.neighbors(adata, use_rep=rep, n_pcs=30)

    # ── 1. PAGA ──────────────────────────────────────────────────────────
    sc.tl.paga(adata, groups=groupby)

    # PAGA-initialised UMAP for downstream figure generation. Skipped if
    # the AnnData already carries a precomputed UMAP we should respect.
    if "X_umap" not in adata.obsm:
        try:
            sc.tl.umap(adata, init_pos="paga")
        except Exception as e:
            log_msg = f"sc.tl.umap(init_pos=paga) skipped: {e}"
            adata.uns.setdefault("paga_log", []).append(log_msg)

    paga_conn: dict = {}
    max_off_diag = 0.0
    try:
        conn_mat = adata.uns["paga"]["connectivities"]
        if hasattr(conn_mat, "toarray"):
            conn_mat = conn_mat.toarray()
        cats = (list(adata.obs[groupby].cat.categories)
                if hasattr(adata.obs[groupby], "cat")
                else sorted(adata.obs[groupby].unique()))
        # We record EVERY off-diagonal connectivity (top_connections below keeps
        # the strongest 25 for display, no floor applied) and separately count the
        # "strong" adjacencies at or above `_PAGA_STRONG_CONNECTIVITY` (0.05) as
        # `n_strong` — see the constant for the single source of truth.
        for i, g1 in enumerate(cats):
            for j, g2 in enumerate(cats):
                if i < j:
                    v = float(conn_mat[i, j])
                    max_off_diag = max(max_off_diag, v)
                    paga_conn[f"{g1}→{g2}"] = round(v, 5)
    except Exception:
        pass

    # Keep all edges, then sort & take top 25 (the report uses ≤15).
    top_connections = dict(
        sorted(paga_conn.items(), key=lambda x: -x[1])[:25]
    )
    n_strong = sum(
        1 for v in paga_conn.values() if v >= _PAGA_STRONG_CONNECTIVITY
    )

    # ── 2. Diffusion Pseudotime ───────────────────────────────────────────
    dpt_result: dict = {"computed": False}
    try:
        # Root selection
        #
        # Order of preference:
        #   1. Precomputed `adata.uns["iroot"]` from an upstream curated file.
        #   2. Explicit `root_cell_type` from the caller (user / intent).
        #   3. Auto-detect a progenitor-like cell type from generic wording
        #      ("stem", "progenitor", "precursor") in the annotation column,
        #      matched as WHOLE WORDS (so "stem" never fires on "brainstem").
        #
        # No low-complexity fallback: low n_genes can indicate low-quality cells,
        # so using it as DPT root gives a false sense of biological direction.
        root_used = None
        root_policy = "explicit_or_progenitor_label_required"
        if "iroot" in adata.uns:
            root_used = "precomputed iroot"
            root_policy = "precomputed_iroot"
        elif root_cell_type and cell_type_col in adata.obs.columns:
            mask = adata.obs[cell_type_col] == root_cell_type
            if mask.sum() > 0:
                adata.uns["iroot"] = int(np.where(mask)[0][0])
                root_used = root_cell_type
                root_policy = "explicit_root_cell_type"
        if root_used is None and cell_type_col in adata.obs.columns:
            labels = [
                str(x) for x in adata.obs[cell_type_col].unique() if x
            ]
            # Prefer the largest progenitor cluster so root pseudotime is
            # statistically meaningful — small contaminants are skipped.
            best_label = None
            best_size  = -1
            for lab in labels:
                if _is_progenitor_label(lab):
                    n = int((adata.obs[cell_type_col] == lab).sum())
                    if n > best_size:
                        best_size = n
                        best_label = lab
            if best_label is not None:
                mask = adata.obs[cell_type_col] == best_label
                adata.uns["iroot"] = int(np.where(mask)[0][0])
                root_used = f"auto (progenitor match: {best_label})"
                root_policy = "progenitor_label"

        if root_used is None:
            dpt_result = {
                "computed": False,
                "reason": "root_unresolved",
                "root_used": None,
                "root_selection_policy": root_policy,
            }
        else:
            sc.tl.diffmap(adata)

            sc.tl.dpt(adata)

            # Mean pseudotime per group
            col = cell_type_col if cell_type_col in adata.obs.columns else groupby
            pt_by_group = (
                adata.obs.groupby(col)["dpt_pseudotime"]
                .mean()
                .sort_values()
                .round(4)
                .to_dict()
            )
            dpt_result = {
                "computed":             True,
                "pseudotime_by_group":  pt_by_group,
                "root_used":            root_used,
                "root_selection_policy": root_policy,
            }
    except Exception as e:
        dpt_result = {"computed": False, "reason": str(e)}

    # ── 3. RNA Velocity (optional) ────────────────────────────────────────
    velocity_result: dict
    if "spliced" in adata.layers and "unspliced" in adata.layers:
        try:
            import scvelo as scv
            scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
            scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
            scv.tl.velocity(adata)
            scv.tl.velocity_graph(adata)
            velocity_result = {"computed": True, "method": "scvelo_stochastic"}
        except ImportError:
            velocity_result = {"computed": False, "reason": "scvelo not installed"}
        except Exception as e:
            velocity_result = {"computed": False, "reason": str(e)}
    else:
        velocity_result = {
            "computed": False,
            "reason":   "no spliced/unspliced layers — provide loom/raw data for velocity",
        }

    output_path = str(Path(output_dir) / "trajectory.h5ad")
    adata.write_h5ad(output_path)

    return {
        "status":      "success",
        "paga":        {"top_connections":   top_connections,
                        "n_connections":     len(paga_conn),
                        "n_strong":          n_strong,
                        "max_connectivity":  round(max_off_diag, 5),
                        "strong_threshold":  0.05},
        "pseudotime":  dpt_result,
        "velocity":    velocity_result,
        "groupby":     groupby,
        "output_path": output_path,
    }


if __name__ == "__main__":
    run_script(rna_trajectory)
