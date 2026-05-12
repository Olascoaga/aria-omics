"""
ARIA NarrativeAgent — scRNA / pseudobulk extension
---------------------------------------------------
Helpers that produce the scRNA-specific text, figures, and HTML blocks for
the NarrativeAgent. Kept in a separate module so the main agent stays
focused on cross-modal orchestration and the scRNA logic can grow without
bloating it.

The functions here are pure (no LLM, no bus). They:

  - read the same `agent_results["scrna_agent"]["findings"]["scRNA"]["findings"]`
    dict that NarrativeAgent already consumes
  - produce text summaries, methods blocks, and HTML cards
  - optionally call out to env_manager to render UMAPs in the rna stack
  - render pathway dotplots in-process via aria.scripts.rna_pathway_viz

Shape consumed (see aria/scripts/rna_narrative_adapter.py for construction):

    findings: {
      qc, integration, clustering, clustering_decision, cell_types,
      differential_expression, pathways,
      pseudobulk_de:        { groupby, condition_col, replicate_col,
                              covariates, thresholds, n_groups, per_group },
      pseudobulk_pathways:  { organism, databases, per_cluster },
      figures:              { umap_<key>: png_path, pathway_dotplots: {...},
                              per_celltype_de_bar: png_path },
    }
"""

from __future__ import annotations

import base64
import html
import logging
from pathlib import Path
from typing import Optional

log = logging.getLogger("aria.narrative.scrna")


# ── Text summaries ────────────────────────────────────────────────────────

def summarize_scrna_text(findings: dict) -> str:
    """Multi-line text summary for findings_sections['scrna']."""
    lines = []

    qc = findings.get("qc") or {}
    if qc:
        n_b = qc.get("n_cells_before")
        n_a = qc.get("n_cells_after")
        if n_b and n_a:
            lines.append(
                f"After QC, {n_a:,} of {n_b:,} cells were retained "
                f"({qc.get('pct_removed', 0)}% removed across "
                f"{qc.get('n_samples', '?')} samples)."
            )
        elif n_a:
            lines.append(f"After QC, {n_a:,} cells were retained.")

    integ = findings.get("integration") or {}
    if integ.get("status") in ("done", "success"):
        s_b = integ.get("silhouette_before")
        s_a = integ.get("silhouette_after")
        method = integ.get("method", "harmony")
        if s_b is not None and s_a is not None:
            lines.append(
                f"Batch correction ({method}) across "
                f"{integ.get('n_batches', '?')} batches: "
                f"silhouette {s_b:+.3f} → {s_a:+.3f} "
                f"(Δ={integ.get('batch_correction_delta', 0):+.3f}; "
                f"lower silhouette indicates better mixing)."
            )

    clu = findings.get("clustering") or {}
    if clu.get("n_clusters"):
        lines.append(
            f"Leiden clustering identified {clu['n_clusters']} clusters "
            f"at resolution {clu.get('resolution', '?')}."
        )

    ct = (findings.get("cell_types") or {}).get("cell_types", {}) or {}
    if ct:
        unique = sorted({str(v) for v in ct.values() if v})
        if unique:
            lines.append(
                f"Cell-type annotation labelled {len(unique)} unique types "
                f"(top: {', '.join(unique[:5])}{'…' if len(unique) > 5 else ''})."
            )

    pb = findings.get("pseudobulk_de") or {}
    if pb:
        n_groups = pb.get("n_groups", 0)
        per_group = pb.get("per_group", {}) or {}
        n_with_de = sum(
            1 for g in per_group.values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "success" and c.get("n_significant", 0) > 0
        )
        thr = pb.get("thresholds", {}) or {}
        lines.append(
            f"Pseudobulk DE (DESeq2 on pseudosamples) ran across "
            f"{n_groups} {pb.get('groupby', 'group')}s. "
            f"{n_with_de} (group × comparison) blocks yielded significant DE "
            f"at padj < {thr.get('padj_max', 0.05)} and "
            f"|log2FC| > {thr.get('lfc_min', 0.5)}."
        )

    pwp = findings.get("pseudobulk_pathways") or {}
    if pwp.get("per_cluster"):
        n_blocks = len(pwp["per_cluster"])
        n_sig_blocks = sum(
            1 for b in pwp["per_cluster"].values()
            if b.get("n_significant", 0) > 0
        )
        lines.append(
            f"Pathway over-representation (Enrichr) on top-200 DE genes per "
            f"(group × comparison): {n_sig_blocks}/{n_blocks} blocks "
            f"with significant enrichment."
        )

    return "\n".join(lines) if lines else (
        "scRNA analysis completed. See findings table for details."
    )


# ── Methods block ─────────────────────────────────────────────────────────

def build_scrna_methods(findings: dict) -> str:
    """Methods section text for the report (scRNA + optional pseudobulk)."""
    lines = []

    qc = findings.get("qc") or {}
    if qc:
        mt = qc.get("mt_threshold")
        mt_str = f" with mt-fraction cap at {mt}%" if mt else ""
        lines.append(
            f"Raw count matrices were processed using scanpy. "
            f"Cells were filtered by adaptive MAD-based thresholds on "
            f"total_counts, n_genes, and percent.mt{mt_str}. Doublets were "
            f"flagged with Scrublet. Counts were normalised to 10,000 "
            f"per cell and log1p-transformed."
        )

    integ = findings.get("integration") or {}
    if integ.get("status") in ("done", "success"):
        lines.append(
            f"Batch correction was performed with "
            f"{integ.get('method', 'Harmony')} on the {integ.get('rep_used', 'X_pca')} "
            f"representation across the '{integ.get('batch_col', 'batch')}' "
            f"covariate. Mixing quality was assessed by silhouette score on "
            f"the corrected embedding."
        )

    clu = findings.get("clustering") or {}
    cdec = findings.get("clustering_decision") or {}
    if clu.get("n_clusters"):
        lines.append(
            f"Dimensionality reduction used PCA (50 components) followed by "
            f"k-NN graph construction (k=15) and UMAP visualisation. "
            f"Leiden clustering at resolution={cdec.get('recommended', '?')} "
            f"(selected by silhouette across "
            f"{cdec.get('n_candidates', '?')} candidates) yielded "
            f"{clu['n_clusters']} clusters."
        )

    ct = findings.get("cell_types") or {}
    if ct.get("model_used"):
        lines.append(
            f"Cell-type annotation used CellTypist with model "
            f"'{ct['model_used']}', assigning a majority label per cluster."
        )

    pb = findings.get("pseudobulk_de") or {}
    if pb:
        thr = pb.get("thresholds", {}) or {}
        cov = ", ".join(pb.get("covariates", []) or []) or "none"
        lines.append(
            f"Between-condition differential expression was performed by "
            f"pseudobulk aggregation: raw counts were summed per "
            f"({pb.get('groupby', 'cell_type')} × "
            f"{pb.get('replicate_col', 'replicate')}) and fitted with "
            f"pyDESeq2 (design ~ {pb.get('condition_col', 'condition')} "
            f"+ {cov if cov != 'none' else 'no covariates'}). "
            f"Pseudosamples with < "
            f"{thr.get('min_cells_per_pseudosample', 10)} cells were dropped; "
            f"groups requiring ≥ "
            f"{thr.get('min_replicates_per_condition', 2)} replicates per "
            f"condition. Significance: padj &lt; {thr.get('padj_max', 0.05)}, "
            f"|log2FC| &gt; {thr.get('lfc_min', 0.5)}. "
            f"For Seurat-derived h5ads with log-normalised raw.X, counts were "
            f"recovered as expm1(x) × nCount_RNA / 10000 prior to aggregation."
        )

    pwp = findings.get("pseudobulk_pathways") or {}
    if pwp.get("per_cluster"):
        dbs = list((pwp.get("databases") or {}).keys()) or [
            "GO_BP", "KEGG", "Reactome"
        ]
        lines.append(
            f"Over-representation analysis (gseapy / Enrichr endpoint) was "
            f"run on the top-200 DE genes per (group × comparison) against "
            f"{', '.join(dbs)}. Significance: adjusted p &lt; 0.05."
        )

    return "\n\n".join(lines)


# ── Pseudobulk DE table (HTML rows) ───────────────────────────────────────

def extract_pseudobulk_de_table(findings: dict,
                                  top_genes_per_row: int = 4) -> str:
    """
    Build an HTML <tbody> with one row per (group × comparison) summarising
    pseudobulk DE. Returns the inner rows (caller wraps in <table>).
    """
    pb = findings.get("pseudobulk_de") or {}
    per_group = pb.get("per_group", {}) or {}
    if not per_group:
        return ""

    rows = []
    # Sort groups by max n_significant across their comparisons (desc)
    def _max_sig(group_info):
        return max(
            (c.get("n_significant", 0)
             for c in (group_info.get("per_comparison", {}) or {}).values()),
            default=-1,
        )
    ordered = sorted(per_group.items(),
                     key=lambda kv: _max_sig(kv[1]),
                     reverse=True)
    for group, info in ordered:
        n_ps = info.get("n_pseudosamples", "?")
        comps = info.get("per_comparison", {}) or {}
        if not comps:
            rows.append(
                f"<tr><td>{html.escape(str(group))}</td>"
                f"<td>{n_ps}</td><td colspan='4'>"
                f"<em>{html.escape(str(info.get('reason', 'no comparison')))}</em>"
                f"</td></tr>"
            )
            continue
        for comp_key, comp in comps.items():
            if comp.get("status") == "skipped":
                rows.append(
                    f"<tr><td>{html.escape(str(group))}</td>"
                    f"<td>{n_ps}</td><td>{html.escape(str(comp_key))}</td>"
                    f"<td colspan='3' style='color:var(--muted)'>"
                    f"<em>skipped: "
                    f"{html.escape(str(comp.get('reason', '')))}</em></td></tr>"
                )
                continue
            if comp.get("status") != "success":
                continue
            n_sig = comp.get("n_significant", 0)
            n_up   = comp.get("n_up", 0)
            n_down = comp.get("n_down", 0)
            top_genes = comp.get("top_genes", []) or []
            up_tops = [g["gene"] for g in top_genes
                       if g.get("log2fc", 0) > 0][:top_genes_per_row]
            dn_tops = [g["gene"] for g in top_genes
                       if g.get("log2fc", 0) < 0][:top_genes_per_row]
            rows.append(
                f"<tr>"
                f"<td><strong>{html.escape(str(group))}</strong></td>"
                f"<td>{n_ps}</td>"
                f"<td>{html.escape(str(comp_key))}</td>"
                f"<td><strong>{n_sig}</strong></td>"
                f"<td style='color:var(--red)'>{n_up} ↑ "
                f"<span style='color:var(--muted);font-size:0.85em'>"
                f"{html.escape(', '.join(up_tops))}</span></td>"
                f"<td style='color:var(--blue)'>{n_down} ↓ "
                f"<span style='color:var(--muted);font-size:0.85em'>"
                f"{html.escape(', '.join(dn_tops))}</span></td>"
                f"</tr>"
            )
    return "\n".join(rows)


# ── Figure rendering ──────────────────────────────────────────────────────

def render_pathway_dotplots(findings: dict,
                              output_dir: Path,
                              top_n_blocks: int = 8,
                              top_n_terms: int = 12) -> dict:
    """
    Render ORA dotplots for the top-K (group × comparison) blocks of the
    pseudobulk_pathways stage. Uses aria.scripts.rna_pathway_viz directly
    (in-process; needs matplotlib + seaborn + pandas).

    Returns dict mapping block_key → list[png_path] (one per database).
    """
    from aria.scripts import rna_pathway_viz as pviz

    # Prefer pseudobulk_pathways (between-condition); fall back to per-cluster
    # marker pathways from the standard scRNA pipeline. Both share schema.
    pwp = (findings.get("pseudobulk_pathways")
           or findings.get("pathways") or {})
    per_cluster = pwp.get("per_cluster", {}) or {}
    if not per_cluster:
        return {}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Rank blocks by total n_significant pathways across all databases
    ranked = sorted(
        per_cluster.items(),
        key=lambda kv: kv[1].get("n_significant", 0),
        reverse=True,
    )[:top_n_blocks]

    figures: dict = {}
    for block_key, block in ranked:
        results = block.get("results", {}) or {}
        per_db = []
        for db_name, terms in results.items():
            if not terms:
                continue
            safe_block = (
                str(block_key).replace("::", "__").replace(" ", "_")
                .replace("/", "_")
            )
            out_path = output_dir / f"pathway_{safe_block}__{db_name}.png"
            png = pviz.make_ora_dotplot(
                pathways_list=terms,
                db_name=db_name,
                contrast_name=block_key.replace("::", " — "),
                output_path=str(out_path),
                top_n=top_n_terms,
            )
            if png:
                per_db.append(png)
        if per_db:
            figures[block_key] = per_db
    return figures


def render_per_celltype_de_bar(findings: dict, output_path: Path) -> Optional[str]:
    """
    Stacked bar chart of n_up / n_down DE genes per (group × comparison)
    from the pseudobulk stage. Returns PNG path on success, None if no data.
    """
    pb = findings.get("pseudobulk_de") or {}
    per_group = pb.get("per_group", {}) or {}
    rows = []
    for group, info in per_group.items():
        for comp_key, comp in (info.get("per_comparison", {}) or {}).items():
            if comp.get("status") != "success":
                continue
            rows.append((
                f"{group}\n({comp_key})",
                int(comp.get("n_up", 0)),
                int(comp.get("n_down", 0)),
            ))
    if not rows:
        return None

    rows.sort(key=lambda r: r[1] + r[2], reverse=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    labels = [r[0] for r in rows]
    n_up   = np.array([r[1] for r in rows])
    n_down = np.array([r[2] for r in rows])

    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(rows) + 2), 4.2),
                           dpi=160)
    x = np.arange(len(rows))
    ax.bar(x, n_up,   color="#991b1b", label="up",   width=0.7)
    ax.bar(x, -n_down, color="#1d4ed8", label="down", width=0.7)
    ax.axhline(0, color="#1e293b", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=8)
    ax.set_ylabel("DE genes  (up ↑ / down ↓)", fontsize=9)
    ax.set_title("Pseudobulk DE — per cell type",
                 fontsize=10, fontweight="bold")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.tick_params(labelsize=7)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    return str(output_path)


# ── HTML embedding ────────────────────────────────────────────────────────

def _embed_png(path: str) -> str:
    """Inline a PNG as a base64 data URI (returns '' on failure)."""
    try:
        p = Path(path)
        if not p.exists():
            return ""
        data = p.read_bytes()
        b64 = base64.b64encode(data).decode("ascii")
        return f"data:image/png;base64,{b64}"
    except Exception as e:
        log.warning(f"_embed_png failed for {path}: {e}")
        return ""


def build_scrna_html_section(findings: dict,
                              max_pathway_blocks: int = 8) -> str:
    """
    Build the inner-HTML for the scRNA findings card: UMAPs, per-cell-type
    DE bar, DE table, and pathway dotplot mosaic. Returns concatenated
    HTML snippet (no outer <div class="card"> — caller wraps).
    """
    parts: list[str] = []
    figs = findings.get("figures") or {}

    # 1. UMAP figures ─────────────────────────────────────────────────────
    umaps = {k: v for k, v in figs.items() if k.startswith("umap_")}
    if umaps:
        parts.append('<h4 style="margin-top:1rem">Embedding</h4>')
        parts.append('<div style="display:flex;flex-wrap:wrap;gap:1rem">')
        for key, path in sorted(umaps.items()):
            uri = _embed_png(path)
            if not uri:
                continue
            caption = html.escape(key.replace("umap_", "UMAP — "))
            parts.append(
                f'<figure style="flex:1 1 320px;min-width:300px;max-width:480px">'
                f'<img src="{uri}" alt="{caption}">'
                f'<figcaption>{caption}</figcaption>'
                f'</figure>'
            )
        parts.append('</div>')

    # 2. Per-cell-type DE summary bar ─────────────────────────────────────
    de_bar = figs.get("per_celltype_de_bar")
    if de_bar:
        uri = _embed_png(de_bar)
        if uri:
            parts.append('<h4>Pseudobulk DE — counts per cell type</h4>')
            parts.append(
                f'<figure><img src="{uri}" '
                f'alt="Per cell-type DE bar"></figure>'
            )

    # 3. Pseudobulk DE table ──────────────────────────────────────────────
    table_rows = extract_pseudobulk_de_table(findings)
    if table_rows:
        parts.append('<h4>DE summary by (cell type × comparison)</h4>')
        parts.append(
            '<table style="width:100%;font-size:0.85em">'
            '<thead><tr>'
            '<th>Cell type</th>'
            '<th>n<sub>pseudo</sub></th>'
            '<th>Comparison</th>'
            '<th>Sig.</th>'
            '<th>Up (top genes)</th>'
            '<th>Down (top genes)</th>'
            '</tr></thead>'
            f'<tbody>{table_rows}</tbody>'
            '</table>'
        )

    # 4. Pathway dotplots ─────────────────────────────────────────────────
    pw_figs = figs.get("pathway_dotplots") or {}
    if pw_figs:
        parts.append(
            '<h4 style="margin-top:1.4rem">'
            'Pathway enrichment — top cell types</h4>'
        )
        # Render up to N blocks, two-column grid
        n = 0
        parts.append(
            '<div style="display:grid;'
            'grid-template-columns:repeat(auto-fit, minmax(340px, 1fr));'
            'gap:1rem">'
        )
        for block_key, png_list in pw_figs.items():
            if n >= max_pathway_blocks:
                break
            for png in png_list:
                uri = _embed_png(png)
                if not uri:
                    continue
                # Database label = filename suffix between __ and .png
                db_label = (
                    Path(png).stem.rsplit("__", 1)[-1]
                    if "__" in Path(png).stem else "ORA"
                )
                caption = html.escape(
                    f"{block_key.replace('::', ' — ')}  ·  {db_label}"
                )
                parts.append(
                    f'<figure><img src="{uri}" alt="{caption}">'
                    f'<figcaption>{caption}</figcaption></figure>'
                )
            n += 1
        parts.append('</div>')

    return "\n".join(parts)


# ── Top-level orchestration helper ────────────────────────────────────────

def generate_figures(findings: dict,
                      h5ad_path: Optional[str],
                      output_dir: Path,
                      env_manager=None,
                      umap_color_keys: Optional[list] = None) -> dict:
    """
    Generate all scRNA figures and write their paths back into
    findings['figures']. Mutates `findings` and returns it.

    Args:
        findings:        scRNA findings dict (mutated in place).
        h5ad_path:       path to AnnData with UMAP. Required for UMAP figs.
        output_dir:      where to write PNGs.
        env_manager:     ARIA EnvironmentManager (so UMAP runs in rna stack).
                         If None, UMAP rendering is skipped.
        umap_color_keys: obs columns to color the UMAP by. If None, picks
                         sensible defaults from the design.

    Returns:
        findings dict with `figures` populated:
            {umap_<key>: png_path, per_celltype_de_bar: png_path,
             pathway_dotplots: {block_key: [png_paths]}}
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figs = findings.setdefault("figures", {})

    # 1. UMAP figures via the rna stack ─────────────────────────────────
    if h5ad_path and env_manager is not None:
        if umap_color_keys is None:
            keys: list = []
            pb = findings.get("pseudobulk_de") or {}
            if pb.get("groupby"):
                keys.append(pb["groupby"])
            if pb.get("condition_col"):
                keys.append(pb["condition_col"])
            if not keys:
                # Standard mode: pick a sensible set in priority order. The
                # rna_figure_umap script silently skips missing columns, so
                # listing redundant fallbacks is safe.
                ct = findings.get("cell_types") or {}
                for candidate in (
                    "cell_type_celltypist",
                    ct.get("label_col"),
                    "leiden",
                    findings.get("integration", {}).get("batch_col"),
                    "batch",
                    "sample_id",
                ):
                    if candidate and candidate not in keys:
                        keys.append(candidate)
            umap_color_keys = keys
        if umap_color_keys:
            try:
                res = env_manager.run_in_stack(
                    stack="rna",
                    script_path="aria/scripts/rna_figure_umap.py",
                    params={
                        "h5ad_path":  str(h5ad_path),
                        "color_by":   umap_color_keys,
                        "output_dir": str(output_dir),
                    },
                )
                if res.get("status") == "success":
                    for key, path in (res.get("figures") or {}).items():
                        figs[f"umap_{key}"] = path
                else:
                    log.warning(
                        f"UMAP figure generation failed: "
                        f"{res.get('error_type')} — {res.get('details', '')[:200]}"
                    )
            except Exception as e:
                log.warning(f"UMAP figure subprocess crashed: {e}")

    # 2. Per-cell-type DE summary bar ───────────────────────────────────
    bar_path = render_per_celltype_de_bar(
        findings, output_dir / "pseudobulk_de_per_celltype_bar.png"
    )
    if bar_path:
        figs["per_celltype_de_bar"] = bar_path

    # 3. Pathway dotplots ───────────────────────────────────────────────
    pw_figs = render_pathway_dotplots(findings, output_dir / "pathways")
    if pw_figs:
        figs["pathway_dotplots"] = pw_figs

    return findings
