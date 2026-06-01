"""
ARIA RNA Per-Cluster Pathway Enrichment
----------------------------------------
ORA (over-representation analysis) per Leiden cluster. P1-7/W-PRIV: the default
engine is a LOCAL hypergeometric test against versioned GMT libraries and each
cell type's expressed-gene background, so gene lists never leave the machine.
Enrichr (network) is used only for databases lacking a local GMT and only when
opted in (ARIA_ALLOW_ENRICHR=1, air-gapped off). Same database choices as
rna_bulk_de.py so the user gets consistent pathway vocabulary across modalities.
Executed inside aria-rna-env.

Input params:
    de_genes_by_cluster: {cluster_id: [{gene, log2fc, padj, ...}, ...]}
                         Typically the output of rna_de_per_cluster.py.
    organism:            str — "Homo sapiens" | "Mus musculus" | ...
    top_genes_per_cluster: int (optional) — max genes submitted per cluster
                                              (default: 200)
    background_genes:     list[str] (optional) — expressed/detectable genes
                          in the analyzed dataset. Used as the ORA universe
                          instead of Enrichr's default all-database universe.
                          Acts as the fallback universe when a cluster has no
                          entry in background_genes_by_cluster.
    background_genes_by_cluster: {cluster_id: [str]} (optional) — per-cluster
                          ORA universe (C2, audit 2026-05-29). When present for
                          a cluster, the genes tested in THAT cell type's
                          pseudobulk are the universe, instead of the global
                          dataset background, which inflates per-cluster
                          enrichment.
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
    background_by_cluster_in = params.get("background_genes_by_cluster") or {}
    padj_db_max          = float(params.get("padj_db_max", 0.05))
    output_dir           = params.get("output_dir") or "."

    if not de_by_cluster:
        return {"status": "skipped",
                "reason": "no de_genes_by_cluster provided"}

    # P1-7 / W-PRIV: default to LOCAL hypergeometric ORA against versioned GMT
    # libraries (offline — gene lists never leave the machine). Enrichr is used
    # only for databases lacking a local GMT and only when opted in
    # (ARIA_ALLOW_ENRICHR=1) AND air-gapped mode is off.
    from aria.utils import ora as _ora
    from aria.utils.privacy import egress_allowed

    gene_sets   = _get_gene_sets(organism)   # {label: enrichr/GMT library name}
    enrichr_org = _gseapy_organism(organism)

    # Load local versioned libraries once (shared across all clusters).
    local_libs: dict = {}          # label -> (gene_sets_dict, version)
    gene_set_versions: dict = {}   # label -> version manifest
    for db_label, lib_name in gene_sets.items():
        loaded = _ora.load_local_library(lib_name)
        if loaded:
            local_libs[db_label] = loaded
            gene_set_versions[db_label] = loaded[1]
    missing = [lbl for lbl in gene_sets if lbl not in local_libs]

    use_enrichr = bool(missing) and _ora.enrichr_opt_in() and egress_allowed()
    gp = None
    if use_enrichr:
        try:
            import gseapy as gp  # noqa: F401
        except ImportError:
            gp = None
            use_enrichr = False

    # Nothing to enrich with: skip honestly (no fabrication) and explain how to
    # enable ORA, instead of silently leaking via Enrichr.
    if not local_libs and not use_enrichr:
        if _ora.enrichr_opt_in() and not egress_allowed():
            reason = "air_gapped_egress_blocked"
        elif _ora.enrichr_opt_in():
            reason = "gseapy_missing"
        else:
            reason = "no_local_gene_sets"
        return {
            "status":  "skipped",
            "reason":  reason,
            "details": (
                "Per-cluster ORA was skipped: no local versioned GMT library "
                "is available for the requested databases and Enrichr is "
                "opt-in. Provision GMTs (python scripts/fetch_genesets.py) or "
                "set ARIA_ALLOW_ENRICHR=1 (with air-gapped mode off) to use "
                "Enrichr."
            ),
            "databases_skipped": list(gene_sets.keys()),
        }

    def _clean_genes(genes) -> list:
        return sorted({
            str(g) for g in (genes or [])
            if g and str(g).lower() != "nan"
        })

    # Global/default universe — fallback when a cluster has no per-cluster
    # universe in background_genes_by_cluster.
    background_genes = _clean_genes(background_genes_in)
    default_background_source = (
        "dataset_expressed_genes"
        if background_genes else
        "all_geneset_genes"
    )
    per_cluster: dict = {}
    csv_rows:    list = []

    # Enrichr (opt-in fallback) is rate-limited (~8s between calls); the local
    # path has no such limit. Track Enrichr calls only.
    enrichr_first_call = True

    def _enrichr_cluster(symbols, db_name, cl_background):
        nonlocal enrichr_first_call
        if not enrichr_first_call:
            time.sleep(8)
        enrichr_first_call = False
        enrichr_kwargs = {
            "gene_list": symbols,
            "gene_sets": db_name,
            "organism": enrichr_org,
            "outdir": None,
            "verbose": False,
        }
        if cl_background:
            enrichr_kwargs["background"] = cl_background
        try:
            enr = gp.enrichr(**enrichr_kwargs)
        except TypeError as exc:
            if "background" not in str(exc):
                raise
            enrichr_kwargs.pop("background", None)
            enr = gp.enrichr(**enrichr_kwargs)
        if enr.results is None or enr.results.empty:
            return []
        sig = enr.results[enr.results["Adjusted P-value"] < padj_db_max]
        sig = sig.sort_values("Adjusted P-value").head(20)
        return [
            {
                "term":           str(row["Term"]),
                "padj":           round(float(row["Adjusted P-value"]), 5),
                "overlap":        str(row.get("Overlap", "")),
                "odds_ratio":     round(float(row.get("Odds Ratio", 1)), 2),
                "combined_score": round(float(row.get("Combined Score", 0)), 1),
                "genes":          str(row.get("Genes", "")).split(";")[:10],
            }
            for _, row in sig.iterrows()
        ]

    for cluster_id, gene_records in de_by_cluster.items():
        # C2 (audit 2026-05-29): resolve THIS cluster's ORA universe. A
        # per-cluster universe (genes tested in that cell type's pseudobulk)
        # avoids the enrichment inflation a global background causes.
        cl_bg_in = background_by_cluster_in.get(cluster_id)
        if cl_bg_in:
            cl_background = _clean_genes(cl_bg_in)
            cl_background_source = "cluster_expressed_genes"
        else:
            cl_background = background_genes
            cl_background_source = default_background_source
        cl_background_set = set(cl_background)

        if not gene_records:
            per_cluster[cluster_id] = {
                "n_input_genes": 0,
                "results":       {},
                "n_significant": 0,
                "background_size": len(cl_background),
                "background_source": cl_background_source,
            }
            continue

        # Order genes by significance × magnitude so top_n captures the most
        # informative signal (the ORA noise floor rises with marginal hits).
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
        if cl_background_set:
            symbols = [g for g in symbols if g in cl_background_set]
        if not symbols:
            per_cluster[cluster_id] = {
                "n_input_genes": 0,
                "results":       {},
                "n_significant": 0,
                "background_size": len(cl_background),
                "background_source": cl_background_source,
            }
            continue

        cluster_results: dict = {}
        cluster_n_sig = 0
        for db_label, db_name in gene_sets.items():
            try:
                if db_label in local_libs:
                    gene_sets_dict, _ = local_libs[db_label]
                    entries = _ora.run_ora(
                        symbols, gene_sets_dict, cl_background,
                        padj_max=padj_db_max, top=20,
                    )
                elif use_enrichr:
                    entries = _enrichr_cluster(symbols, db_name, cl_background)
                else:
                    entries = []
                for entry in entries:
                    csv_rows.append({
                        "cluster":  cluster_id, "database": db_label,
                        **{k: (v if not isinstance(v, list) else ";".join(map(str, v)))
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
            "background_size": len(cl_background),
            "background_source": cl_background_source,
        }

    csv_path = None
    try:
        csv_path = str(Path(output_dir) / "pathways_per_cluster.csv")
        pd.DataFrame(csv_rows).to_csv(csv_path, index=False)
    except Exception:
        csv_path = None

    if local_libs and use_enrichr:
        ora_method = "mixed_local_enrichr"
    elif local_libs:
        ora_method = "local_hypergeometric"
    else:
        ora_method = "enrichr"

    # Top-level disclosure: when per-cluster universes were supplied, the ORA
    # universe is cell-type-specific, so report that rather than a single global
    # size. The per_cluster entries carry each cluster's exact universe.
    used_per_cluster_bg = bool(background_by_cluster_in)
    return {
        "status":       "success",
        "organism":     organism,
        "databases":    gene_sets,
        # P1-7/W-PRIV: ORA engine + exact versioned gene-set release per database.
        "ora_method":   ora_method,
        "gene_set_versions": gene_set_versions,
        "databases_skipped": ([] if use_enrichr else list(missing)),
        "background_size": (None if used_per_cluster_bg
                            else len(background_genes)),
        "background_source": ("per_cluster_expressed_genes"
                              if used_per_cluster_bg
                              else default_background_source),
        "per_cluster":  per_cluster,
        "output_csv":   csv_path,
    }


if __name__ == "__main__":
    run_script(rna_pathway_per_cluster)
