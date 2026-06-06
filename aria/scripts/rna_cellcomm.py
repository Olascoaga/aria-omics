"""
ARIA RNA Cell-Cell Communication Script
-----------------------------------------
Infers ligand-receptor interactions between cell types in scRNA-seq data.
Executed inside aria-rna-env by EnvironmentManager (standalone entry point).

Primary method: LIANA rank_aggregate (if liana-py installed)
Fallback:       mean-expression scoring with a curated L-R resource

Input params:
    data_path:      str  — path to annotated .h5ad
    groupby:        str  (optional) — obs column with cell types (default: "cell_type").
                         Legacy alias `cell_type_col` is still accepted (P0-1).
    organism:       str  (optional) — "Homo sapiens" | "Mus musculus" (default: "Homo sapiens")
    output_dir:     str  (optional)
    n_perms:        int  (optional) — LIANA permutations (default: 1000)

Output:
    {
      "status":           "success" | "skipped" | "error",
      "method":           str,
      "n_cell_types":     int,
      "n_interactions":   int,
      "top_interactions": [{"source", "target", "ligand", "receptor", "score"}, ...],
      "top_pairs":        ["TypeA→TypeB", ...],
      "output_path":      str | None,
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def _resolve_groupby(params: dict) -> str:
    """Resolve the grouping obs column from the IPC params.

    Canonical key is `groupby`; `cell_type_col` is the accepted legacy alias
    (P0-1). Falls back to "cell_type" when neither is provided. Empty/blank
    values are treated as absent so they do not override the alias/default.
    """
    for key in ("groupby", "cell_type_col"):
        value = params.get(key)
        if isinstance(value, str) and value.strip():
            return value
    return "cell_type"


def rna_cellcomm(params: dict) -> dict:
    import pandas as pd
    from pathlib import Path
    from aria.utils.safe_h5ad import read_h5ad

    data_path     = params["data_path"]
    cell_type_col = _resolve_groupby(params)
    output_dir    = params.get("output_dir", str(Path(data_path).parent))
    # C4: 100 permutations was too low for stable LIANA ranks in publication
    # reports. Keep it configurable, but default to 1000.
    n_perms       = int(params.get("n_perms", 1000))

    adata = read_h5ad(data_path)

    # Resolve cell type column
    if cell_type_col not in adata.obs.columns:
        if "leiden" in adata.obs.columns:
            cell_type_col = "leiden"
        else:
            return {"status": "error", "error_type": "NoCellTypes",
                    "details": f"Column '{cell_type_col}' not found and no leiden clusters."}

    n_types = int(adata.obs[cell_type_col].nunique())
    if n_types < 2:
        return {"status": "skipped",
                "reason":  "need ≥2 cell types for communication analysis"}

    # ── Primary: LIANA ────────────────────────────────────────────────────
    interactions: list[dict] = []
    method = "liana_rank_aggregate"
    autocrine_count = 0
    try:
        import liana as li
        li.mt.rank_aggregate(
            adata,
            groupby=cell_type_col,
            use_raw=False,
            verbose=False,
            n_perms=n_perms,
        )
        liana_df = adata.uns["liana_res"].copy()

        # Drop self-interactions (autocrine) — they dominate the output
        # because expression overlap with self is trivially perfect, but
        # they're rarely the biological question for cell-cell comm.
        autocrine_mask = liana_df["source"] == liana_df["target"]
        autocrine_count = int(autocrine_mask.sum())
        liana_df = liana_df[~autocrine_mask].copy()

        # Pick a usable rank column. Recent LIANA versions emit
        # magnitude_rank as NaN when the aggregate doesn't include any
        # magnitude-scoring method, so we fall back to specificity_rank
        # (RRA-aggregated p-value rank, lower = more specific).
        # Prefer specificity_rank for rank_aggregate. magnitude_rank can be
        # present but all-NaN, and in some LIANA versions it is less useful
        # for aggregate outputs. All rank metrics are lower-is-better.
        rank_cols = [
            c for c in ("specificity_rank", "magnitude_rank", "lrscore")
            if c in liana_df.columns
        ]
        chosen_rank = None
        for c in rank_cols:
            if liana_df[c].notna().any():
                chosen_rank = c
                break
        if chosen_rank is None:
            raise RuntimeError(
                f"LIANA returned only NaN ranks across {rank_cols}"
            )

        # specificity_rank / magnitude_rank: lower = better.
        # lrscore: higher = better.
        ascending = chosen_rank != "lrscore"
        ranked = (liana_df
                  .dropna(subset=[chosen_rank])
                  .sort_values(chosen_rank, ascending=ascending)
                  .head(50))
        for rank_idx, (_, row) in enumerate(ranked.iterrows(), start=1):
            raw_score = float(row.get(chosen_rank, 0))
            raw_pval = row.get("cellphone_pvals") if "cellphone_pvals" in row else None
            try:
                pval = float(raw_pval) if raw_pval is not None else None
                if pval <= 0:
                    pval = None
            except Exception:
                pval = None
            interactions.append({
                "source":   str(row.get("source", "")),
                "target":   str(row.get("target", "")),
                "ligand":   str(row.get("ligand_complex",
                                         row.get("ligand", ""))),
                "receptor": str(row.get("receptor_complex",
                                         row.get("receptor", ""))),
                "score":    raw_score,
                "rank":     rank_idx,
                "rank_metric": chosen_rank,
                "score_direction": (
                    "lower is better" if ascending else "higher is better"
                ),
                "cellphone_pval": pval,
            })
        method = f"liana_rank_aggregate ({chosen_rank})"

    except ImportError:
        return {
            "status": "skipped",
            "reason": "liana_not_installed",
            "details": (
                "Cell-cell communication requires LIANA. ARIA will not emit "
                "interactions from an embedded ligand-receptor list; install "
                "LIANA in the RNA environment or skip this analysis."
            ),
            "method": method,
            "n_cell_types": n_types,
            "n_interactions": 0,
            "n_autocrine_dropped": 0,
            "top_interactions": [],
            "top_pairs": [],
            "output_path": None,
        }

    # Summarize by sender→receiver pair
    pair_counts: dict = {}
    for ia in interactions:
        k = f"{ia['source']}→{ia['target']}"
        pair_counts[k] = pair_counts.get(k, 0) + 1
    top_pairs = [p for p, _ in sorted(pair_counts.items(), key=lambda x: -x[1])[:10]]

    # Save CSV
    result_path: str | None = None
    try:
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        result_path = str(out_dir / "cellcomm_interactions.csv")
        pd.DataFrame(interactions).to_csv(result_path, index=False)
    except Exception as e:
        result_path = None
        # Surface the reason in the result for downstream debugging
        method = f"{method} (csv_write_failed: {e})"

    return {
        "status":             "success",
        "method":             method,
        "n_cell_types":       n_types,
        "n_interactions":     len(interactions),
        "n_autocrine_dropped": autocrine_count,
        "top_interactions":   interactions[:20],
        "top_pairs":          top_pairs,
        "output_path":        result_path,
    }


if __name__ == "__main__":
    run_script(rna_cellcomm)
