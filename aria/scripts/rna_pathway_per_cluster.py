"""
ARIA RNA Per-Cluster Pathway Enrichment
----------------------------------------
ORA (over-representation analysis) per Leiden cluster via gseapy/Enrichr.
Same database choices as rna_bulk_de.py so the user gets consistent
pathway vocabulary across modalities. Executed inside aria-rna-env.

Input params:
    de_genes_by_cluster: {cluster_id: [{gene, log2fc, padj, ...}, ...]}
                         Typically the output of rna_de_per_cluster.py.
    organism:            str — "Homo sapiens" | "Mus musculus" | ...
    top_genes_per_cluster: int (optional) — max genes submitted per cluster
                                              (default: 200)
    background_genes:     list[str] (optional) — expressed/detectable genes
                          in the analyzed dataset. Used as the ORA universe
                          instead of Enrichr's default all-database universe.
    padj_db_max:         float (optional) — filter Enrichr hits (default: 0.05)
    output_dir:          str (optional) — CSV destination

Output:
    {
      "status":      "success" | "error",
      "organism":    str,
      "databases":   {db_label: db_enrichr_name, ...},
      "per_cluster": {
          cluster_id: {
              "n_input_genes": int,
              "results": {
                  db_label: [
                      {term, padj, overlap, odds_ratio, combined_score, genes},
                      ...
                  ],
                  ...
              },
              "n_significant": int   # across all databases
          },
          ...
      },
      "output_csv":  str | None
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def _gseapy_organism(organism: str) -> str:
    """gseapy 1.1.13+ requires lowercase organism strings."""
    o = organism.lower()
    if "sapiens" in o or o in ("human", "hsapiens", "h. sapiens", "hs"):
        return "human"
    if "musculus" in o or o in ("mouse", "mus musculus", "mm"):
        return "mouse"
    if "rerio" in o or "zebrafish" in o:
        return "zebrafish"
    if "norvegicus" in o or "rat" in o:
        return "rat"
    if "melanogaster" in o or "drosophila" in o:
        return "fly"
    return "human"  # safe default for Enrichr


def _get_gene_sets(organism: str) -> dict:
    o = organism.lower()
    if "sapiens" in o:
        return {
            "GO_BP":    "GO_Biological_Process_2021",
            "KEGG":     "KEGG_2021_Human",
            "Reactome": "Reactome_2022",
        }
    if "musculus" in o:
        return {
            "GO_BP":    "GO_Biological_Process_2021",
            "KEGG":     "KEGG_2019_Mouse",
            "Reactome": "Reactome_2022",
        }
    return {"GO_BP": "GO_Biological_Process_2021"}


def rna_pathway_per_cluster(params: dict) -> dict:
    import time
    import pandas as pd
    from pathlib import Path

    de_by_cluster        = params.get("de_genes_by_cluster") or {}
    organism             = params.get("organism", "Homo sapiens")
    top_n                = int(params.get("top_genes_per_cluster", 200))
    background_genes_in  = params.get("background_genes") or []
    padj_db_max          = float(params.get("padj_db_max", 0.05))
    output_dir           = params.get("output_dir") or "."

    if not de_by_cluster:
        return {"status": "skipped",
                "reason": "no de_genes_by_cluster provided"}

    try:
        import gseapy as gp
    except ImportError:
        return {
            "status":     "error",
            "error_type": "GseapyMissing",
            "details":    "gseapy not installed in aria-rna-env.",
        }

    gene_sets   = _get_gene_sets(organism)
    enrichr_org = _gseapy_organism(organism)
    background_genes = sorted({
        str(g) for g in background_genes_in
        if g and str(g).lower() != "nan"
    })
    background_set = set(background_genes)
    background_source = (
        "dataset_expressed_genes"
        if background_genes else
        "enrichr_default_universe"
    )
    per_cluster: dict = {}
    csv_rows:    list = []

    # Rate-limit: Enrichr blocks aggressive callers. Same pattern as
    # rna_bulk_de.py — 8s between calls, applied across the (cluster, db)
    # cartesian product.
    first_call = True
    for cluster_id, gene_records in de_by_cluster.items():
        if not gene_records:
            per_cluster[cluster_id] = {
                "n_input_genes": 0,
                "results":       {},
                "n_significant": 0,
                "background_size": len(background_genes),
                "background_source": background_source,
            }
            continue

        # Order genes by significance × magnitude so top_n captures the most
        # informative signal — Enrichr's overlap metric is sensitive to list
        # size and the noise floor rises if you submit too many marginal hits.
        sorted_records = sorted(
            gene_records,
            key=lambda r: (-abs(float(r.get("log2fc", 0))),
                           float(r.get("padj", 1.0))),
        )
        symbols = [
            str(r["gene"])
            for r in sorted_records[:top_n]
            if r.get("gene") and str(r["gene"]) not in ("nan", "")
        ]
        if background_set:
            symbols = [g for g in symbols if g in background_set]
        if not symbols:
            per_cluster[cluster_id] = {
                "n_input_genes": 0,
                "results":       {},
                "n_significant": 0,
                "background_size": len(background_genes),
                "background_source": background_source,
            }
            continue

        cluster_results: dict = {}
        cluster_n_sig = 0
        for db_label, db_name in gene_sets.items():
            if not first_call:
                time.sleep(8)
            first_call = False
            try:
                enrichr_kwargs = {
                    "gene_list": symbols,
                    "gene_sets": db_name,
                    "organism": enrichr_org,
                    "outdir": None,
                    "verbose": False,
                }
                if background_genes:
                    enrichr_kwargs["background"] = background_genes
                try:
                    enr = gp.enrichr(**enrichr_kwargs)
                except TypeError as exc:
                    if "background" not in str(exc):
                        raise
                    enrichr_kwargs.pop("background", None)
                    enr = gp.enrichr(**enrichr_kwargs)
                if enr.results is None or enr.results.empty:
                    cluster_results[db_label] = []
                    continue
                sig = enr.results[enr.results["Adjusted P-value"] < padj_db_max]
                sig = sig.sort_values("Adjusted P-value").head(20)
                entries = []
                for _, row in sig.iterrows():
                    entry = {
                        "term":           str(row["Term"]),
                        "padj":           round(float(row["Adjusted P-value"]), 5),
                        "overlap":        str(row.get("Overlap", "")),
                        "odds_ratio":     round(float(row.get("Odds Ratio", 1)), 2),
                        "combined_score": round(float(row.get("Combined Score", 0)), 1),
                        "genes":          str(row.get("Genes", "")).split(";")[:10],
                    }
                    entries.append(entry)
                    csv_rows.append({
                        "cluster":  cluster_id, "database": db_label,
                        **{k: (v if not isinstance(v, list) else ";".join(v))
                           for k, v in entry.items()}
                    })
                cluster_results[db_label] = entries
                cluster_n_sig += len(entries)
            except Exception as e:
                cluster_results[db_label] = []
                csv_rows.append({
                    "cluster":  cluster_id, "database": db_label,
                    "term":     f"ERROR: {type(e).__name__}: {str(e)[:120]}",
                    "padj":     None, "overlap": None,
                    "odds_ratio": None, "combined_score": None, "genes": None,
                })

        per_cluster[cluster_id] = {
            "n_input_genes": len(symbols),
            "results":       cluster_results,
            "n_significant": cluster_n_sig,
            "background_size": len(background_genes),
            "background_source": background_source,
        }

    csv_path = None
    try:
        csv_path = str(Path(output_dir) / "pathways_per_cluster.csv")
        pd.DataFrame(csv_rows).to_csv(csv_path, index=False)
    except Exception:
        csv_path = None

    return {
        "status":       "success",
        "organism":     organism,
        "databases":    gene_sets,
        "background_size": len(background_genes),
        "background_source": background_source,
        "per_cluster":  per_cluster,
        "output_csv":   csv_path,
    }


if __name__ == "__main__":
    run_script(rna_pathway_per_cluster)
