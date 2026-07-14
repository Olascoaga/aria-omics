"""scRNA narrative figures: matplotlib renders, HTML section, figure generation.

Extracted verbatim from aria/agents/_narrative_scrna.py (A7); behavior pinned by
tests/test_narrative_scrna_contract.py and tests/test_scrna_e2e.py.
"""
from __future__ import annotations

import base64
import html
import logging
from pathlib import Path
from typing import Optional

from aria.agents.narrative.scrna._common import *  # noqa: F401,F403
from aria.agents.narrative.scrna.tables import (
    extract_cellcomm_table,
    extract_pseudobulk_de_table,
    extract_trajectory_tables,
)

log = logging.getLogger("aria.narrative.scrna")


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


def render_cellcomm_heatmap(findings: dict,
                              output_path: Path) -> Optional[str]:
    """
    Heatmap of n_interactions per (source × target). Top_interactions only.
    """
    ccc = findings.get("cell_communication") or {}
    top = ccc.get("top_interactions") or []
    if not top:
        return None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    df = pd.DataFrame(top)
    if df.empty or "source" not in df.columns or "target" not in df.columns:
        return None

    counts = (df.groupby(["source", "target"]).size()
              .reset_index(name="n"))
    cell_types = sorted(set(counts["source"]).union(counts["target"]))
    mat = pd.DataFrame(0, index=cell_types, columns=cell_types, dtype=int)
    for _, r in counts.iterrows():
        mat.loc[r["source"], r["target"]] = int(r["n"])

    n = len(cell_types)
    fig, ax = plt.subplots(figsize=(max(5, 0.55 * n + 2),
                                     max(4, 0.55 * n + 1.5)),
                           dpi=160)
    im = ax.imshow(mat.values, cmap="magma_r", aspect="auto",
                   vmin=0, vmax=max(1, mat.values.max()))
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(cell_types, rotation=55, ha="right", fontsize=8)
    ax.set_yticklabels(cell_types, fontsize=8)
    ax.set_xlabel("Receiver", fontsize=9)
    ax.set_ylabel("Sender", fontsize=9)
    ax.set_title(
        f"Cell-cell communication — interactions among top {len(top)} "
        f"(autocrine excluded)",
        fontsize=10, fontweight="bold",
    )
    # Cell-level labels
    for i in range(n):
        for j in range(n):
            v = mat.values[i, j]
            if v > 0:
                ax.text(j, i, str(v),
                        ha="center", va="center", fontsize=7,
                        color="white" if v >= mat.values.max() * 0.5
                              else "#1e293b")
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02,
                 label="n interactions").ax.tick_params(labelsize=7)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    return str(output_path)


def render_cellcomm_top_pairs_bar(findings: dict,
                                    output_path: Path,
                                    top_n: int = 15) -> Optional[str]:
    """Horizontal barplot of top-N L-R interactions by rank/score."""
    ccc = findings.get("cell_communication") or {}
    top = ccc.get("top_interactions") or []
    if not top:
        return None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Ranks: lower = better for spec/mag. Prefer explicit rank order when
    # present so tied/underflowed LIANA scores do not render as all-zero bars.
    rows = []
    for ia in top[:top_n]:
        label = (f"{ia.get('source', '?')[:18]} → "
                 f"{ia.get('target', '?')[:18]}  "
                 f"({ia.get('ligand', '?')}-{ia.get('receptor', '?')})")
        rows.append((label, float(ia.get("rank", ia.get("score", 0)))))
    if not rows:
        return None
    labels = [r[0] for r in rows]
    scores = np.array([r[1] for r in rows], dtype=float)

    is_rank = (ccc.get("method") or "").startswith("liana")
    if is_rank:
        # Invert: higher bar = lower rank = better.
        max_score = scores.max() if scores.max() > 0 else 1.0
        plot_vals = max_score - scores + (max_score * 0.05)
        xlabel = f"strength (inverted rank; metric: {ccc.get('method')})"
    else:
        plot_vals = scores
        xlabel = "score"

    fig, ax = plt.subplots(figsize=(8, max(3, 0.32 * len(rows) + 1.2)),
                           dpi=160)
    y = np.arange(len(rows))
    ax.barh(y, plot_vals, color="#0d9488", edgecolor="#0f172a",
            linewidth=0.4, height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()  # first = top
    ax.set_xlabel(xlabel, fontsize=8.5)
    ax.set_title("Top ligand-receptor interactions",
                 fontsize=10, fontweight="bold")
    ax.tick_params(labelsize=7)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    return str(output_path)


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
    ann = _annotation_state(findings)
    label_unit = "cell type" if ann["has_valid"] else "cluster"
    if ann["is_marker_fallback"]:
        parts.append(
            '<div class="warning">'
            'Cell labels in this report are unresolved fallback labels because '
            'CellTypist did not complete. Treat UMAP, '
            'trajectory, and communication labels as curation targets, not '
            'final cell identities.'
            '</div>'
        )
    elif not ann["has_valid"]:
        parts.append(
            '<div class="warning">'
            'Cell-type annotation did not produce usable biological labels. '
            'Embedding, trajectory, and communication results are reported at '
            'Leiden-cluster resolution.'
            '</div>'
        )

    # 1. UMAP figures ─────────────────────────────────────────────────────
    # Exclude trajectory-specific UMAPs (e.g. dpt_pseudotime) — those are
    # rendered inside the Trajectory section to keep the narrative
    # contiguous.
    umaps = {
        k: v for k, v in figs.items()
        if k.startswith("umap_") and k != "umap_dpt_pseudotime"
    }
    if umaps:
        embedding_label = findings.get("embedding_label") or "UMAP"
        parts.append('<h4 style="margin-top:1rem">Embedding</h4>')
        parts.append('<div style="display:flex;flex-wrap:wrap;gap:1rem">')
        for key, path in sorted(umaps.items()):
            uri = _embed_png(path)
            if not uri:
                continue
            pretty = key.replace("umap_", "")
            pretty = {
                "cell_type_marker": "marker-based cell label",
                "cell_type_celltypist": "CellTypist cell label",
                "leiden": "Leiden cluster",
                "batch": "batch",
                "sample_id": "sample",
            }.get(pretty, pretty)
            caption = html.escape(f"{embedding_label} — {pretty}")
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
            parts.append(f'<h4>Pseudobulk DE — counts per {label_unit}</h4>')
            parts.append(
                f'<figure><img src="{uri}" '
                f'alt="Per cell-type DE bar"></figure>'
            )

    # 3. Pseudobulk DE table ──────────────────────────────────────────────
    table_rows = extract_pseudobulk_de_table(findings)
    if table_rows:
        parts.append(f'<h4>DE summary by ({label_unit} × comparison)</h4>')
        parts.append(
            '<table style="width:100%;font-size:0.85em">'
            '<thead><tr>'
                f'<th>{html.escape(label_unit.title())}</th>'
            '<th>n<sub>pseudo</sub></th>'
            '<th>Comparison</th>'
            '<th>Sig.</th>'
            '<th>Up (top genes)</th>'
            '<th>Down (top genes)</th>'
            '</tr></thead>'
            f'<tbody>{table_rows}</tbody>'
            '</table>'
        )

    # 4. Trajectory section (PAGA + DPT) ──────────────────────────────────
    traj = findings.get("trajectory") or {}
    if traj.get("status") in ("done", "success"):
        parts.append('<h4 style="margin-top:1.4rem">'
                     f'Trajectory — PAGA + DPT by {label_unit}</h4>')
        parts.append(
            '<p style="color:var(--muted);font-size:0.88em">'
            'This section reports graph connectivity and DPT ordering for '
            f'{html.escape(label_unit)} groups. It is an exploratory manifold '
            'summary, not causal evidence of differentiation by itself.'
            '</p>'
        )

        # PAGA + DPT-coloured UMAP figures, side by side
        traj_figs = []
        for fkey in ("paga_graph", "paga_log10_graph",
                     "umap_dpt_pseudotime"):
            p = figs.get(fkey)
            if p:
                uri = _embed_png(p)
                if uri:
                    cap = html.escape({
                        "paga_graph":         "PAGA — group connectivity",
                        "paga_log10_graph":   "PAGA — log-scaled edges",
                        "umap_dpt_pseudotime": "UMAP — DPT pseudotime",
                    }[fkey])
                    traj_figs.append(
                        f'<figure style="flex:1 1 280px;min-width:260px;'
                        f'max-width:430px">'
                        f'<img src="{uri}" alt="{cap}">'
                        f'<figcaption>{cap}</figcaption>'
                        f'</figure>'
                    )
        if traj_figs:
            parts.append('<div style="display:flex;flex-wrap:wrap;'
                         'gap:1rem">')
            parts.extend(traj_figs)
            parts.append('</div>')

        # Tables: PAGA top connections + DPT pseudotime by group
        tables = extract_trajectory_tables(findings)
        if tables["paga_rows"]:
            paga_meta = traj.get("paga", {}) or {}
            max_c = paga_meta.get("max_connectivity", 0)
            n_str = paga_meta.get("n_strong", 0)
            thr_  = paga_meta.get("strong_threshold", 0.05)
            note = (
                f'<p style="color:var(--muted);font-size:0.82em;'
                f'margin-top:0.4rem;font-style:italic">'
                f'Max connectivity = {max_c:.4f}. '
                f'{n_str} edge(s) above the {thr_} threshold. '
                f'In mature / non-developmental populations, absolute '
                f'connectivities are typically &lt; 0.01 — interpret '
                f'rankings rather than absolute magnitudes.</p>'
            )
            parts.append(
                '<h4 style="margin-top:1.2rem">PAGA — top connections</h4>'
                '<table style="width:100%;font-size:0.88em">'
                '<thead><tr><th>Edge</th><th>Connectivity</th>'
                '<th>Strength</th><th>Visual</th></tr></thead>'
                f'<tbody>{tables["paga_rows"]}</tbody>'
                '</table>'
                + note
            )
        if tables["pseudotime_rows"]:
            pt = traj.get("pseudotime", {}) or {}
            root_str = html.escape(str(pt.get("root_used", "auto")))
            parts.append(
                f'<h4 style="margin-top:1.2rem">DPT pseudotime by '
                f'{html.escape(label_unit)} group '
                f'(root: {root_str})</h4>'
                '<table style="width:100%;font-size:0.88em">'
                '<thead><tr><th>Rank</th><th>Group</th>'
                '<th>Mean DPT</th><th>Visual</th></tr></thead>'
                f'<tbody>{tables["pseudotime_rows"]}</tbody>'
                '</table>'
            )

    # 4b. Cell-cell communication section ─────────────────────────────────
    ccc = findings.get("cell_communication") or {}
    if ccc.get("status") in ("done", "success"):
        parts.append('<h4 style="margin-top:1.4rem">'
                     f'Cell-cell communication by {label_unit}</h4>')
        parts.append(
            '<p style="color:var(--muted);font-size:0.88em">'
            'Ligand-receptor scores are summarized between observed '
            f'{html.escape(label_unit)} groups. These results require '
            'manual review of sender and receiver labels before biological '
            'interpretation.'
            '</p>'
        )

        ccc_figs = []
        for fkey, caption in (
            ("cellcomm_heatmap",   "Sender → receiver interaction count"),
            ("cellcomm_top_pairs", "Top ligand-receptor interactions"),
        ):
            p = figs.get(fkey)
            if p:
                uri = _embed_png(p)
                if uri:
                    ccc_figs.append(
                        f'<figure style="flex:1 1 320px;min-width:300px;'
                        f'max-width:520px"><img src="{uri}" '
                        f'alt="{html.escape(caption)}">'
                        f'<figcaption>{html.escape(caption)}</figcaption>'
                        f'</figure>'
                    )
        if ccc_figs:
            parts.append('<div style="display:flex;flex-wrap:wrap;'
                         'gap:1rem">')
            parts.extend(ccc_figs)
            parts.append('</div>')

        cc_rows = extract_cellcomm_table(findings)
        if cc_rows:
            method = html.escape(str(ccc.get("method", "?")))
            n_ct = ccc.get("n_cell_types", "?")
            n_int = ccc.get("n_interactions", "?")
            n_auto = ccc.get("n_autocrine_dropped", 0)
            parts.append(
                f'<h4 style="margin-top:1rem">'
                f'Top L-R interactions  '
                f'<span style="color:var(--muted);font-weight:400;'
                f'font-size:0.85em">'
                f'({method} · {n_int} interactions across {n_ct} '
                f'{html.escape(label_unit)} groups · {n_auto} autocrine '
                f'pairs excluded)</span></h4>'
                '<table style="width:100%;font-size:0.85em">'
                '<thead><tr><th>Sender</th><th>Receiver</th>'
                '<th>Ligand</th><th>Receptor</th>'
                '<th>Rank</th><th>Metric value</th><th>CellPhone p</th>'
                '</tr></thead>'
                f'<tbody>{cc_rows}</tbody></table>'
            )

    # 5. Pathway dotplots ─────────────────────────────────────────────────
    pw_figs = figs.get("pathway_dotplots") or {}
    if pw_figs:
        parts.append(
            '<h4 style="margin-top:1.4rem">'
            f'Pathway enrichment — top {label_unit} groups</h4>'
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

    # 6. Supplementary table links ────────────────────────────────────────
    table_links = findings.get("tables") or {}
    if table_links:
        labels = {
            "qc_per_sample": "QC per sample",
            "cell_types": "Cell types",
            "cluster_markers": "Cluster markers",
            "pseudobulk_de_summary": "Pseudobulk DE summary",
            "pseudobulk_de_genes": "Pseudobulk DE genes",
            "pathway_enrichment": "Pathway enrichment",
            "cellcomm_interactions": "Cell-cell communication",
            "paga_connections": "PAGA connections",
            "pseudotime_by_group": "Pseudotime by group",
        }
        links = []
        for key, path in table_links.items():
            name = labels.get(key, key.replace("_", " ").title())
            rel = f"tables/{Path(path).name}"
            links.append(
                f'<a href="{html.escape(rel)}" style="color:var(--blue);'
                f'text-decoration:underline">{html.escape(name)}</a>'
            )
        if links:
            parts.append(
                '<h4 style="margin-top:1.4rem">Supplementary tables</h4>'
                '<p style="font-size:0.85rem;color:var(--muted)">'
                + " &middot; ".join(links)
                + '</p>'
            )

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
            ann = _annotation_state(findings)
            if ann.get("label_col"):
                keys.append(ann["label_col"])
            pb = findings.get("pseudobulk_de") or {}
            if pb.get("groupby") and pb.get("groupby") not in keys:
                keys.append(pb["groupby"])
            if pb.get("condition_col") and pb.get("condition_col") not in keys:
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
            else:
                for candidate in (
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
                    if res.get("embedding_label"):
                        findings["embedding_label"] = res.get("embedding_label")
                    if res.get("embedding_key"):
                        findings["embedding_key"] = res.get("embedding_key")
                    if res.get("embedding_was_computed"):
                        findings["embedding_was_computed"] = True
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

    # 3b. Cell-cell communication figures ───────────────────────────────
    heat_path = render_cellcomm_heatmap(
        findings, output_dir / "cellcomm_heatmap.png"
    )
    if heat_path:
        figs["cellcomm_heatmap"] = heat_path
    bar_path_ccc = render_cellcomm_top_pairs_bar(
        findings, output_dir / "cellcomm_top_pairs.png"
    )
    if bar_path_ccc:
        figs["cellcomm_top_pairs"] = bar_path_ccc

    # 4. Trajectory figures (PAGA graph + DPT UMAP) ─────────────────────
    traj = findings.get("trajectory") or {}
    if traj.get("status") in ("done", "success") and env_manager is not None:
        traj_h5ad = traj.get("output_path") or h5ad_path
        if traj_h5ad:
            try:
                paga_res = env_manager.run_in_stack(
                    stack="rna",
                    script_path="aria/scripts/rna_figure_paga.py",
                    params={
                        "h5ad_path":  str(traj_h5ad),
                        "output_dir": str(output_dir / "trajectory"),
                        "groupby":    traj.get("groupby"),
                    },
                )
                if paga_res.get("status") == "success":
                    figs.update(paga_res.get("figures") or {})
                else:
                    log.warning(
                        f"PAGA figure failed: "
                        f"{paga_res.get('error_type')} — "
                        f"{paga_res.get('details', '')[:200]}"
                    )
            except Exception as e:
                log.warning(f"PAGA figure subprocess crashed: {e}")

            # DPT-coloured UMAP — only if dpt_pseudotime obs col exists.
            pt = traj.get("pseudotime", {}) or {}
            if pt.get("computed"):
                try:
                    dpt_res = env_manager.run_in_stack(
                        stack="rna",
                        script_path="aria/scripts/rna_figure_umap.py",
                        params={
                            "h5ad_path":  str(traj_h5ad),
                            "color_by":   ["dpt_pseudotime"],
                            "output_dir": str(output_dir / "trajectory"),
                        },
                    )
                    if dpt_res.get("status") == "success":
                        path = (dpt_res.get("figures") or {}).get(
                            "dpt_pseudotime"
                        )
                        if path:
                            figs["umap_dpt_pseudotime"] = path
                except Exception as e:
                    log.warning(f"DPT UMAP subprocess crashed: {e}")

    return findings
