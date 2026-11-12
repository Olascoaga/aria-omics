"""B3 (bulk ATAC pre-print): functional over-representation analysis (ORA) of the
genes near differential-accessibility peaks.

State-of-the-art ATAC interpretation asks *what biology* the differential peaks point
to. ARIA reuses two validated cores, adding no new statistics:
- the B2 peak→gene assignment (`chromatin_peak_annotation`): each peak's nearest-TSS
  gene;
- the bulk RNA ORA engine (`rna_bulk.ora._run_pathway_enrichment`): a LOCAL
  hypergeometric test against versioned GMT libraries (GO BP / KEGG / Reactome), with
  Enrichr as an opt-in fallback only.

Foreground = genes near the significant peaks of a comparison. Background = genes near
ALL tested peaks of that comparison (this controls for the assay's gene-accessibility
bias rather than testing against the whole genome). Enrichment is ASSOCIATIVE — that a
pathway is over-represented among genes near DA peaks is a hypothesis, not evidence that
the peaks regulate those genes. Honest-skip without a GTF, without significant genes, or
when no versioned GMT is provisioned (no fabrication).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from aria.scripts.chromatin_peak_annotation import (
    _modality_stamp, _resolve_gtf, _select_comparison_rows, annotate_peaks,
    parse_gtf_gene_models)


def _organism_from_genome(genome) -> str:
    g = str(genome or "").lower()
    if g.startswith(("hg", "grch")) or "sapiens" in g or "human" in g:
        return "Homo sapiens"
    if g.startswith(("mm", "grcm")) or "musculus" in g or "mouse" in g:
        return "Mus musculus"
    return "Homo sapiens"


def _genes_for_peaks(peaks, models) -> list:
    """Unique nearest-TSS gene symbols for an iterable of peak ids (drops peaks with
    no assigned gene; never fabricates)."""
    rows = annotate_peaks(peaks, models)["rows"]
    seen: dict = {}
    for r in rows:
        g = r.get("nearest_gene")
        if g:
            seen[g] = None
    return list(seen)


def chromatin_peak_ora(params: dict) -> dict:
    """Dispatchable entry: functional ORA of genes near each ATAC comparison's DA
    peaks (bulk ATAC DA, or scATAC pseudobulk DA when ``data_type`` is
    ``"scATAC"``), against a peak-proximal background. Honest-skip without a GTF /
    sig genes."""
    import pandas as pd

    data_type, validation_level, _prefix = _modality_stamp(params)
    # Bulk ATAC keeps its historical (unprefixed) figure keys/filenames.
    fig_prefix = "scatac_" if data_type == "scATAC" else ""

    out_dir = Path(params.get("output_dir", "."))
    out_dir.mkdir(parents=True, exist_ok=True)

    gtf_path = _resolve_gtf(params, out_dir)
    if not gtf_path:
        return {"status": "skipped", "ran": False, "analysis": "peak_ora",
                "data_type": data_type, "validation_level": validation_level,
                "reason": "no_gtf_for_peak_ora",
                "details": ("No GTF found. Place one at "
                            "~/.aria/genomes/<assembly>/annotation.gtf[.gz], pass "
                            "gtf=<path>, or set ARIA_GTF=<path>.")}
    models = parse_gtf_gene_models(str(gtf_path))
    if not models:
        return {"status": "skipped", "ran": False, "analysis": "peak_ora",
                "data_type": data_type, "validation_level": validation_level,
                "reason": "gtf_parsed_zero_gene_models",
                "gtf": Path(gtf_path).name}

    organism = (params.get("organism")
                or _organism_from_genome(params.get("genome")))
    from aria.scripts.rna_bulk.ora import _run_pathway_enrichment
    try:
        from aria.scripts.rna_pathway_viz import make_ora_dotplot
    except Exception:  # seaborn/matplotlib optional — enrichment still runs
        make_ora_dotplot = None

    warnings: list = []
    comp_results: list = []
    figures: dict = {}
    for comp in params.get("comparisons") or []:
        csv_path = (comp or {}).get("full_results_csv")
        if not csv_path or not Path(csv_path).exists():
            continue
        df = _select_comparison_rows(pd.read_csv(csv_path), comp)
        peak_col = "peak" if "peak" in df.columns else df.columns[0]
        sig = df
        if "significant" in df.columns:
            sig = df[df["significant"].fillna(False).astype(bool)]
        if sig.empty:
            continue

        sig_peaks = [str(p) for p in sig[peak_col].tolist()]
        fg_genes = _genes_for_peaks(sig_peaks, models)
        if not fg_genes:
            continue
        # Directional foreground (up/down by accessibility log2FC).
        up_genes: list = []
        down_genes: list = []
        # bulk ATAC DA tables carry `log2FoldChange`; the scATAC pseudobulk
        # full table carries the same estimate as `log2fc`.
        lfc_col = next((c for c in ("log2FoldChange", "log2fc")
                        if c in sig.columns), None)
        if lfc_col is not None:
            up_peaks = [str(p) for p in
                        sig.loc[sig[lfc_col] > 0, peak_col].tolist()]
            down_peaks = [str(p) for p in
                          sig.loc[sig[lfc_col] < 0, peak_col].tolist()]
            up_genes = _genes_for_peaks(up_peaks, models)
            down_genes = _genes_for_peaks(down_peaks, models)
        # Peak-proximal background: genes near ALL tested peaks.
        bg_genes = _genes_for_peaks(
            [str(p) for p in df[peak_col].tolist()], models)

        comp_key = f"{comp.get('test')}_vs_{comp.get('reference')}"
        pathways, pw_warn, ora_meta = _run_pathway_enrichment(
            fg_genes, up_genes, down_genes, organism, str(out_dir),
            symbol_map=None, background_genes=bg_genes, allow_mock=False)
        warnings.extend(f"[{comp_key}] {w}" for w in pw_warn)

        comp_figs: dict = {}
        if make_ora_dotplot is not None:
            for db, terms in pathways.items():
                if not terms:
                    continue
                key = f"{fig_prefix}peak_ora_{db}_{comp_key}"
                path = make_ora_dotplot(
                    terms, db, comp_key.replace("_", " "),
                    str(out_dir / f"{key}.png"))
                if path:
                    comp_figs[key] = path
                    figures[key] = path

        comp_results.append({
            "test": comp.get("test"), "reference": comp.get("reference"),
            "n_foreground_genes": len(fg_genes),
            "n_background_genes": len(bg_genes),
            "ora_method": ora_meta.get("method"),
            "gene_set_versions": ora_meta.get("gene_set_versions"),
            "databases_skipped": ora_meta.get("databases_skipped"),
            "pathways_summary": {db: len(terms)
                                 for db, terms in pathways.items()},
            "pathways": pathways,
            "figures": comp_figs,
        })

    if not comp_results:
        return {"status": "skipped", "ran": False, "analysis": "peak_ora",
                "data_type": data_type, "validation_level": validation_level,
                "reason": "no_enrichable_genes_near_da_peaks",
                "gtf": Path(gtf_path).name, "warnings": warnings}

    return {
        "status": "success", "ran": True, "analysis": "peak_ora",
        "data_type": data_type, "validation_level": validation_level,
        "method": ("functional ORA (local hypergeometric over versioned GMTs) of "
                   "genes near DA peaks, against a peak-proximal background"),
        "organism": organism,
        "gtf": Path(gtf_path).name,
        "comparisons": comp_results,
        "figures": figures,
        "warnings": warnings,
    }


if __name__ == "__main__":
    import json

    if len(sys.argv) == 2 and sys.argv[1].lstrip().startswith("{"):
        print(json.dumps(chromatin_peak_ora(json.loads(sys.argv[1]))))
    else:
        from aria.scripts._base import run_script

        run_script(chromatin_peak_ora)
