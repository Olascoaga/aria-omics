"""Bulk RNA-seq over-representation analysis / pathway helpers (P2-8 split).

Behavior-preserving extraction; re-exported from `aria.scripts.rna_bulk_de`."""

from __future__ import annotations

import math

from aria.scripts.rna_bulk.gtf_io import _normalize_species, _to_symbols


def _finite_float(value, default=None):
    """Coerce external numeric fields without emitting JSON NaN/Infinity."""
    try:
        value = float(value)
        return value if math.isfinite(value) else default
    except (TypeError, ValueError):
        return default


def _enrichr_entries(results, padj_max=0.05, top=20):
    """Normalize Enrichr rows; retain full precision for p-values and FDR."""
    entries = []
    for _, row in results.iterrows():
        padj = _finite_float(row.get("Adjusted P-value"))
        if padj is None or not 0 <= padj < padj_max:
            continue
        genes = row.get("Genes")
        genes = genes.split(";") if isinstance(genes, str) else (
            genes if isinstance(genes, (list, tuple)) else [])
        odds = _finite_float(row.get("Odds Ratio"))
        score = _finite_float(row.get("Combined Score"))
        entries.append({
            "term": str(row["Term"]),
            "padj": padj,
            "adjusted_p": padj,
            "pvalue": _finite_float(row.get("P-value")),
            "overlap": row.get("Overlap") if isinstance(row.get("Overlap"), str) else "",
            "odds_ratio": round(odds, 2) if odds is not None else None,
            "combined_score": round(score, 1) if score is not None else None,
            "genes": [g.strip() for g in genes if isinstance(g, str)
                      and g.strip() and g.strip().lower() != "nan"][:10],
        })
    entries.sort(key=lambda r: (
        r["padj"], _finite_float(r["pvalue"], math.inf),
        -_finite_float(r["combined_score"], -math.inf), r["term"]))
    return entries[:top]


def _enrichr_enrichment(sig_symbols: list, db_labels: list, gene_sets: dict,
                        organism: str, background_symbols: list,
                        allow_mock: bool, original_sig_genes: list,
                        failed_databases: dict | None = None) -> tuple:
    """Enrichr ORA for the given database labels (opt-in, network egress).

    Only invoked for databases that have no local versioned GMT, and only when
    the caller already confirmed ARIA_ALLOW_ENRICHR + egress is allowed. The
    [:500] submission cap is an Enrichr server-side practicality (its URL/GET
    payload is bounded); the local hypergeometric path has no such cap."""
    pathways: dict = {}
    warnings:  list = []
    try:
        import gseapy as gp
        import time as _time

        for i, db in enumerate(db_labels):
            db_name = gene_sets[db]
            # Rate-limit: Enrichr blocks aggressive callers (~5-15s between calls).
            if i > 0:
                _time.sleep(8)
            try:
                enrichr_kwargs = {
                    "gene_list": sig_symbols[:500],
                    "gene_sets": db_name,
                    "organism": _gseapy_organism(organism),
                    "outdir": None,
                    "verbose": False,
                }
                if background_symbols:
                    enrichr_kwargs["background"] = background_symbols
                try:
                    enr = gp.enrichr(**enrichr_kwargs)
                except TypeError as exc:
                    if "background" not in str(exc):
                        raise
                    enrichr_kwargs.pop("background", None)
                    enr = gp.enrichr(**enrichr_kwargs)
                results = enr.results

                if results is not None and not results.empty:
                    pathways[db] = _enrichr_entries(results)
                    if not pathways[db]:
                        warnings.append(
                            f"No significant {db} pathways at padj < 0.05 "
                            f"({len(results)} terms tested via Enrichr)."
                        )
                else:
                    warnings.append(
                        f"{db} returned empty results from Enrichr "
                        f"(0 of {len(sig_symbols)} symbols matched)."
                    )
            except Exception as e:  # aria: allow-broad-except auxiliary Enrichr network ORA records per-database failure
                err_txt = str(e)[:200]
                if isinstance(failed_databases, dict):
                    failed_databases[db] = err_txt
                warnings.append(
                    f"{db} Enrichr ERROR: {err_txt}. "
                    f"(Sent {len(sig_symbols)} symbols, organism="
                    f"{_gseapy_organism(organism)}.)"
                )
    except ImportError:
        if allow_mock:
            warnings.append(
                "gseapy not installed — returning mock pathway enrichment "
                "(explicit mock mode)."
            )
            pathways = _mock_pathways(original_sig_genes)
        else:
            warnings.append(
                "gseapy not installed — Enrichr pathway enrichment skipped. "
                "No simulated pathway terms were generated."
            )
    return pathways, warnings

def _run_pathway_enrichment(sig_genes: list, up_genes: list,
                              down_genes: list, organism: str,
                              output_dir: str,
                              symbol_map: dict = None,
                              background_genes: list | None = None,
                              allow_mock: bool = False) -> tuple:
    """Pathway over-representation analysis (ORA).

    P1-7 / W-PRIV: the DEFAULT engine is a LOCAL hypergeometric test against
    versioned GMT libraries and the dataset's expressed-gene background, so the
    DE gene list never leaves the machine. Enrichr (which ships the list to an
    external server) is used ONLY for databases lacking a local GMT and ONLY
    when the user opts in (ARIA_ALLOW_ENRICHR=1) AND air-gapped mode is off;
    otherwise those databases are skipped with an explicit caveat (no
    fabrication). The local path imposes no [:500] submission cap.

    Returns (pathways, warnings, ora_meta). `ora_meta` records the method and
    the exact gene-set release per database for `methodology.json`.
    """
    from aria.utils import ora as _ora
    from aria.utils.privacy import egress_allowed

    warnings  = []
    pathways  = {}
    gene_sets = _get_gene_sets(organism)   # {label: enrichr/GMT library name}
    ora_meta  = {
        "method": "none",
        "gene_set_versions": {},
        "databases_local": [],
        "databases_enrichr": [],
        "databases_skipped": [],
        "failed_databases": {},
    }

    if not sig_genes:
        return {}, ["No significant genes for pathway enrichment."], ora_meta

    # Convert IDs → symbols if we have a mapping
    sig_symbols  = _to_symbols(sig_genes,  symbol_map or {})
    up_symbols   = _to_symbols(up_genes,   symbol_map or {})
    down_symbols = _to_symbols(down_genes, symbol_map or {})
    background_symbols = sorted({
        str(g) for g in (background_genes or [])
        if g and str(g).lower() != "nan"
    })
    background_set = {g.upper() for g in background_symbols}
    if background_set:
        sig_symbols = [g for g in sig_symbols if g.upper() in background_set]
        up_symbols = [g for g in up_symbols if g.upper() in background_set]
        down_symbols = [g for g in down_symbols if g.upper() in background_set]

    if symbol_map and len(sig_symbols) < len(sig_genes) * 0.3:
        warnings.append(
            f"Symbol mapping yielded only {len(sig_symbols)}/{len(sig_genes)} "
            f"genes (< 30%). Pathway enrichment may be underpowered."
        )

    if not sig_symbols:
        warnings.append(
            "No valid gene symbols to enrich — check that GTF gene_name "
            "annotations are present."
        )
        return {}, warnings, ora_meta

    # 1) Local hypergeometric ORA (default; offline; works even air-gapped).
    local_results, versions, missing = _ora.local_ora_for_databases(
        sig_symbols, gene_sets, background_symbols, padj_max=0.05, top=20,
    )
    for db, terms in local_results.items():
        pathways[db] = terms
        ora_meta["databases_local"].append(db)
        ora_meta["gene_set_versions"][db] = versions[db]
        if not terms:
            warnings.append(
                f"No significant {db} pathways at padj < 0.05 (local ORA)."
            )
    if local_results:
        ora_meta["method"] = "local_hypergeometric"

    # 2) Databases with no local GMT: Enrichr only when opted in + egress allowed.
    if missing:
        if _ora.enrichr_opt_in() and egress_allowed():
            failed_dbs: dict[str, str] = {}
            enr_pw, enr_warn = _enrichr_enrichment(
                sig_symbols, missing, gene_sets, organism,
                background_symbols, allow_mock, sig_genes,
                failed_databases=failed_dbs,
            )
            pathways.update(enr_pw)
            warnings.extend(enr_warn)
            ora_meta["databases_enrichr"] = list(enr_pw.keys())
            if failed_dbs:
                ora_meta["failed_databases"] = failed_dbs
            if enr_pw:
                ora_meta["method"] = (
                    "mixed_local_enrichr" if local_results else "enrichr"
                )
        else:
            ora_meta["databases_skipped"] = list(missing)
            if _ora.enrichr_opt_in() and not egress_allowed():
                reason = "ARIA_AIR_GAPPED is enabled (Enrichr egress refused)"
            else:
                reason = ("no local versioned GMT is available and Enrichr is "
                          "opt-in")
            warnings.append(
                f"Pathway ORA skipped for {', '.join(missing)}: {reason}. "
                f"Provision versioned GMTs (python scripts/fetch_genesets.py) "
                f"or set ARIA_ALLOW_ENRICHR=1 to use Enrichr."
            )

    # 3) Directional KEGG (up vs down) — local-only, complementary to combined.
    kegg_lib = gene_sets.get("KEGG")
    if kegg_lib and (up_symbols or down_symbols):
        loaded = _ora.load_local_library(kegg_lib)
        if loaded:
            kegg_sets, _ = loaded
            for direction, genes in [("up", up_symbols), ("down", down_symbols)]:
                if not genes:
                    continue
                terms = _ora.run_ora(
                    genes, kegg_sets, background_symbols, padj_max=0.1, top=10,
                )
                if terms:
                    pathways[f"KEGG_{direction}"] = terms

    return pathways, warnings, ora_meta

def _get_gene_sets(organism: str) -> dict:
    """Return gene set databases appropriate for organism."""
    species = _normalize_species(organism)
    if species == "Homo sapiens":
        return {
            "GO_BP":    "GO_Biological_Process_2021",
            "KEGG":     "KEGG_2021_Human",
            "Reactome": "Reactome_2022",
        }
    elif species == "Mus musculus":
        return {
            "GO_BP":    "GO_Biological_Process_2021",
            "KEGG":     "KEGG_2019_Mouse",
            "Reactome": "Reactome_2022",
        }
    else:
        return {"GO_BP": "GO_Biological_Process_2021"}

def _gseapy_organism(organism: str) -> str:
    """
    Return the organism string in the format gseapy/Enrichr expects.

    gseapy 1.1.13+ explicitly validates against a lowercase whitelist:
        human, hsapiens, h. sapiens, hs, mouse, mus musculus, ...
    Capitalized strings like "Human" or "Homo sapiens" are REJECTED with:
        "Invalid organism 'Human'. Valid options are: ..."
    """
    s = _normalize_species(organism).lower()
    if "sapiens" in s or s in ("human", "hs", "hsapiens", "h. sapiens"):
        return "human"
    if "musculus" in s or s in ("mouse", "mm"):
        return "mouse"
    if "rerio" in s or s in ("zebrafish", "danio rerio"):
        return "fish"
    if "melanogaster" in s or s in ("fly", "drosophila"):
        return "fly"
    if "norvegicus" in s or s == "rat":
        return "rat"
    if "elegans" in s:
        return "worm"
    if "cerevisiae" in s or s == "yeast":
        return "yeast"
    return "human"   # safest default for the most common use case

def _mock_pathways(sig_genes: list) -> dict:
    """Neutral mock pathways for explicit mock mode only."""
    return {
        "GO_BP": [
            {"term": "mock pathway A", "padj": 0.001,
             "genes": sig_genes[:5], "overlap": "12/234"},
            {"term": "mock pathway B", "padj": 0.003,
             "genes": sig_genes[:3], "overlap": "8/156"},
        ],
        "KEGG": [
            {"term": "mock pathway C", "padj": 0.005,
             "genes": sig_genes[:4], "overlap": "6/87"},
        ],
        "note": "mock — install gseapy for real enrichment",
    }
