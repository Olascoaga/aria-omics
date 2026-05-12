"""
ARIA RNA Cell-Cell Communication Script
-----------------------------------------
Infers ligand-receptor interactions between cell types in scRNA-seq data.
Executed inside aria-rna-env by EnvironmentManager (standalone entry point).

Primary method: LIANA rank_aggregate (if liana-py installed)
Fallback:       mean-expression scoring with a curated L-R resource

Input params:
    data_path:      str  — path to annotated .h5ad
    cell_type_col:  str  (optional) — obs column with cell types (default: "cell_type")
    organism:       str  (optional) — "Homo sapiens" | "Mus musculus" (default: "Homo sapiens")
    output_dir:     str  (optional)
    n_perms:        int  (optional) — LIANA permutations (default: 100)

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

# Curated ligand-receptor pairs (conserved, high-confidence)
_LR_PAIRS_HUMAN = [
    ("TGFB1",  "TGFBR1"), ("TGFB1",  "TGFBR2"),
    ("VEGFA",  "KDR"),    ("VEGFA",  "FLT1"),
    ("MIF",    "CD44"),   ("SPP1",   "CD44"),
    ("SPP1",   "ITGAV"),  ("IL6",    "IL6R"),
    ("IL6",    "IL6ST"),  ("CXCL12", "CXCR4"),
    ("CCL5",   "CCR5"),   ("CCL5",   "CCR1"),
    ("EGF",    "EGFR"),   ("HGF",    "MET"),
    ("PDGFB",  "PDGFRB"), ("FGF2",   "FGFR1"),
    ("WNT5A",  "ROR2"),   ("NOTCH1", "JAG1"),
    ("DLL4",   "NOTCH1"), ("SEMA3A", "NRP1"),
    ("ANGPT1", "TEK"),    ("BMP2",   "BMPR1A"),
    ("BMP4",   "BMPR2"),  ("KITLG",  "KIT"),
    ("THBS1",  "CD36"),   ("FN1",    "ITGA5"),
    ("LAMB1",  "ITGB1"),  ("CDH1",   "ITGAE"),
    ("EFNB2",  "EPHB4"),  ("EFNA1",  "EPHA2"),
    ("IGF1",   "IGF1R"),  ("PDGFA",  "PDGFRA"),
    ("TNF",    "TNFRSF1A"),("IFNG",  "IFNGR1"),
    ("IL10",   "IL10RA"), ("IL4",    "IL4R"),
]


def rna_cellcomm(params: dict) -> dict:
    import scanpy as sc
    import numpy as np
    import pandas as pd
    from pathlib import Path

    data_path     = params["data_path"]
    cell_type_col = params.get("cell_type_col", "cell_type")
    organism      = params.get("organism", "Homo sapiens").lower()
    output_dir    = params.get("output_dir", str(Path(data_path).parent))
    n_perms       = int(params.get("n_perms", 100))

    adata = sc.read_h5ad(data_path)

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

    interactions: list[dict] = []
    method = "mean_expression_fallback"

    # ── Primary: LIANA ────────────────────────────────────────────────────
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
        rank_cols = [
            c for c in ("magnitude_rank", "specificity_rank", "lrscore")
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
        for _, row in ranked.iterrows():
            interactions.append({
                "source":   str(row.get("source", "")),
                "target":   str(row.get("target", "")),
                "ligand":   str(row.get("ligand_complex",
                                         row.get("ligand", ""))),
                "receptor": str(row.get("receptor_complex",
                                         row.get("receptor", ""))),
                "score":    round(float(row.get(chosen_rank, 0)), 4),
                "rank_metric": chosen_rank,
                "cellphone_pval": float(row.get("cellphone_pvals"))
                                   if "cellphone_pvals" in row else None,
            })
        method = f"liana_rank_aggregate ({chosen_rank})"

    except (ImportError, Exception):
        # ── Fallback: curated mean-expression scoring ─────────────────────
        lr_pairs = _LR_PAIRS_HUMAN
        if "musculus" in organism:
            lr_pairs = [(l.capitalize(), r.capitalize()) for l, r in lr_pairs]

        var_names = set(adata.var_names)
        valid_pairs = [(l, r) for l, r in lr_pairs
                       if l in var_names and r in var_names]

        if not valid_pairs:
            return {"status": "skipped",
                    "reason": "no curated L-R pairs found in gene names"}

        cell_types = adata.obs[cell_type_col].unique()
        var_idx    = {g: i for i, g in enumerate(adata.var_names)}

        mean_expr: dict = {}
        for ct in cell_types:
            mask = adata.obs[cell_type_col] == ct
            X_sub = adata.X[mask]
            if hasattr(X_sub, "toarray"):
                X_sub = X_sub.toarray()
            mean_expr[ct] = np.asarray(X_sub, dtype=float).mean(axis=0)

        for l_gene, r_gene in valid_pairs:
            li_idx, ri_idx = var_idx[l_gene], var_idx[r_gene]
            for src in cell_types:
                for tgt in cell_types:
                    if src == tgt:
                        continue
                    score = float(mean_expr[src][li_idx]) * float(mean_expr[tgt][ri_idx])
                    if score > 0:
                        interactions.append({
                            "source":   str(src),  "target":   str(tgt),
                            "ligand":   l_gene,    "receptor": r_gene,
                            "score":    round(score, 4),
                        })

        interactions.sort(key=lambda x: -x["score"])
        interactions = interactions[:50]
        method = "mean_expression_fallback"

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
