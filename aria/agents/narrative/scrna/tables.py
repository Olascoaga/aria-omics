"""scRNA narrative tables: HTML table bodies and supplementary TSV export.

Extracted verbatim from aria/agents/_narrative_scrna.py (A7); behavior pinned by
tests/test_narrative_scrna_contract.py.
"""
from __future__ import annotations

import csv
import html
import math
from pathlib import Path
from typing import Optional

from aria.agents.narrative.scrna._common import *  # noqa: F401,F403


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
            if comp.get("status") == "error":
                err_msg = (
                    comp.get("error")
                    or comp.get("reason")
                    or comp.get("details")
                    or "failed"
                )
                rows.append(
                    f"<tr><td>{html.escape(str(group))}</td>"
                    f"<td>{n_ps}</td><td>{html.escape(str(comp_key))}</td>"
                    f"<td colspan='3' style='color:var(--red)'>"
                    f"<em>error: "
                    f"{html.escape(str(err_msg))}</em></td></tr>"
                )
                continue
            if comp.get("status") != "success":
                continue
            # Wave 4A: counts follow the PRIMARY FDR family (fdr_strategy) —
            # the same family that populated top_genes — so the row can never
            # read "0 ↑ GeneA" when per-cluster BH found the up genes.
            counts = _primary_de_counts(comp, pb)
            n_sig = counts["n_significant"]
            n_up = counts["n_up"]
            n_down = counts["n_down"]
            primary_tag = ("per-cluster" if counts["strategy"] == "per_cluster"
                           else "global")
            secondary = counts["n_secondary"]
            secondary_html = (
                f" / {secondary} {counts['secondary_label']}"
                if secondary is not None else ""
            )
            top_genes = [g for g in (comp.get("top_genes", []) or [])
                         if isinstance(g, dict)]
            up_tops = [str(g.get("gene", "?")) for g in top_genes
                       if (_finite_number(g.get("log2fc")) or 0) > 0
                       ][:top_genes_per_row]
            dn_tops = [str(g.get("gene", "?")) for g in top_genes
                       if (_finite_number(g.get("log2fc")) or 0) < 0
                       ][:top_genes_per_row]
            rows.append(
                f"<tr>"
                f"<td><strong>{html.escape(str(group))}</strong></td>"
                f"<td>{n_ps}</td>"
                f"<td>{html.escape(str(comp_key))}</td>"
                f"<td><strong>{n_sig}</strong> "
                f"<span style='color:var(--muted);font-size:0.85em'>"
                f"{primary_tag}{html.escape(secondary_html)}</span></td>"
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

def extract_cellcomm_table(findings: dict, top_n: int = 20) -> str:
    """HTML <tbody> rows for top L-R interactions."""
    ccc = findings.get("cell_communication") or {}
    top = (ccc.get("top_interactions") or [])[:top_n]
    if not top:
        return ""
    rows = []
    for ia in top:
        src = html.escape(str(ia.get("source", "")))
        tgt = html.escape(str(ia.get("target", "")))
        lig = html.escape(str(ia.get("ligand", "")))
        rec = html.escape(str(ia.get("receptor", "")))
        score = ia.get("score", "?")
        rank = ia.get("rank")
        metric = ia.get("rank_metric") or (
            (ccc.get("method", "").split("(")[-1].rstrip(")").strip())
            if ccc.get("method") else ""
        )
        pval = ia.get("cellphone_pval")
        # Wave 4A: p == 0.0 is a valid finite permutation p (no permutation
        # exceeded the observed statistic) — render it; NaN/inf/negative -> "—".
        pval_f = _finite_number(pval)
        pval_str = (f"<code>{_fmt_stat(pval_f)}</code>"
                    if pval_f is not None and pval_f >= 0 else "—")
        rank_str = f"#{int(rank)}" if isinstance(rank, (int, float)) else "—"
        score_str = _fmt_stat(score)
        metric_str = html.escape(str(metric or "score"))
        rows.append(
            f"<tr><td>{src}</td><td>{tgt}</td>"
            f"<td><strong>{lig}</strong></td><td>{rec}</td>"
            f"<td>{rank_str}</td><td><code>{score_str}</code><br>"
            f"<span style='color:var(--muted);font-size:0.82em'>"
            f"{metric_str}</span></td><td>{pval_str}</td></tr>"
        )
    return "\n".join(rows)


def extract_trajectory_tables(findings: dict) -> dict:
    """
    Return {paga_rows, pseudotime_rows} HTML strings for the trajectory
    section. Each is a <tbody> inner snippet (caller wraps in <table>).
    """
    traj = findings.get("trajectory") or {}
    paga = traj.get("paga", {}) or {}
    pt = traj.get("pseudotime", {}) or {}

    paga_rows = ""
    top_conn = paga.get("top_connections", {}) or {}
    if top_conn:
        finite_conn = [
            _finite_number(v) for v in top_conn.values()
            if _finite_number(v) is not None
        ]
        max_c = _finite_number(paga.get("max_connectivity")) or max(
            finite_conn, default=0,
        )
        thr = paga.get("strong_threshold")
        thr = 0.05 if _finite_number(thr) is None else float(thr)
        rows = []
        for edge, val in top_conn.items():
            val_f = _finite_number(val)
            strong = val_f is not None and val_f >= thr
            badge = (
                '<span style="background:#dcfce7;color:var(--green);'
                'padding:2px 6px;border-radius:3px;font-size:0.75em;'
                'font-weight:600">strong</span>'
                if strong else
                '<span style="color:var(--muted);font-size:0.85em">weak</span>'
            )
            rel = (val_f / max_c) if (max_c and val_f is not None) else 0
            bar_w = max(2, int(rel * 100)) if math.isfinite(rel) else 2
            rows.append(
                f"<tr><td>{html.escape(str(edge))}</td>"
                f"<td><code>{val}</code></td>"
                f"<td>{badge}</td>"
                f"<td><div style='background:#e2e8f0;width:120px;"
                f"height:8px;border-radius:3px;overflow:hidden'>"
                f"<div style='background:var(--teal);width:{bar_w}px;"
                f"height:100%;'></div></div></td></tr>"
            )
        paga_rows = "\n".join(rows)

    pseudotime_rows = ""
    pt_by = pt.get("pseudotime_by_group", {}) or {}
    if pt_by:
        # Wave 4A: disconnected DPT components give inf/NaN group means; rank
        # and scale bars on finite values only and list non-finite groups last
        # (never int(nan * 140) -> ValueError).
        finite_items = [(g, _finite_number(v)) for g, v in pt_by.items()
                        if _finite_number(v) is not None]
        nonfinite_items = [(g, v) for g, v in pt_by.items()
                           if _finite_number(v) is None]
        ordered = sorted(finite_items, key=lambda kv: kv[1])
        max_pt = max((v for _, v in finite_items), default=0)
        rows = []
        for rank, (group, val) in enumerate(ordered, 1):
            rel = (val / max_pt) if max_pt else 0
            bar_w = max(2, int(rel * 140)) if math.isfinite(rel) else 2
            rows.append(
                f"<tr><td>{rank}</td>"
                f"<td><strong>{html.escape(str(group))}</strong></td>"
                f"<td><code>{val:.4f}</code></td>"
                f"<td><div style='background:#e2e8f0;width:160px;"
                f"height:8px;border-radius:3px;overflow:hidden'>"
                f"<div style='background:var(--blue);width:{bar_w}px;"
                f"height:100%'></div></div></td></tr>"
            )
        for group, val in nonfinite_items:
            rows.append(
                f"<tr><td>—</td>"
                f"<td><strong>{html.escape(str(group))}</strong></td>"
                f"<td><code>{html.escape(str(val))}</code></td>"
                f"<td><em style='color:var(--muted)'>not reached from root "
                f"(disconnected)</em></td></tr>"
            )
        pseudotime_rows = "\n".join(rows)

    return {"paga_rows": paga_rows,
            "pseudotime_rows": pseudotime_rows}


# ── Supplementary table export ────────────────────────────────────────────

def _write_tsv(path: Path, rows: list[dict]) -> Optional[str]:
    """Write rows as TSV. Returns path on success, None if no rows."""
    if not rows:
        return None
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return str(path)


def export_supplementary_tables(findings: dict, output_dir: Path) -> dict:
    """
    Materialize scRNA result objects into report/tables/*.tsv.

    The analytical scripts often return rich in-memory structures for the
    NarrativeAgent but the report staging layer only copied bulk RNA tables.
    This exporter keeps the report directory self-contained for scRNA runs.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables: dict[str, str] = {}

    qc = findings.get("qc") or {}
    qc_rows = []
    for sample in qc.get("per_sample", []) or []:
        if isinstance(sample, dict):
            qc_rows.append(sample)
    p = _write_tsv(output_dir / "scrna_qc_per_sample.tsv", qc_rows)
    if p:
        tables["qc_per_sample"] = p

    ct = (findings.get("cell_types") or {}).get("cell_types", {}) or {}
    ct_rows = []
    for cluster, value in ct.items():
        row = {"cluster": cluster, "label": _label_cell_type(value)}
        if isinstance(value, dict):
            for key, val in value.items():
                if isinstance(val, (str, int, float, bool)) or val is None:
                    row[key] = val
                elif isinstance(val, list):
                    row[key] = ", ".join(map(str, val))
        ct_rows.append(row)
    p = _write_tsv(output_dir / "scrna_cell_types.tsv", ct_rows)
    if p:
        tables["cell_types"] = p

    # Standard per-cluster marker DE.
    de = findings.get("differential_expression") or {}
    marker_rows = []
    for cluster, genes in (de.get("de_genes_by_cluster", {}) or {}).items():
        for gene in genes or []:
            if isinstance(gene, dict):
                row = {"cluster": cluster}
                row.update(gene)
            else:
                row = {"cluster": cluster, "gene": gene}
            marker_rows.append(row)
    p = _write_tsv(output_dir / "scrna_cluster_markers.tsv", marker_rows)
    if p:
        tables["cluster_markers"] = p

    pb = findings.get("pseudobulk_de") or {}
    pb_summary_rows = []
    pb_gene_rows = []
    for group, info in (pb.get("per_group", {}) or {}).items():
        n_ps = info.get("n_pseudosamples")
        for comp_key, comp in (info.get("per_comparison", {}) or {}).items():
            pb_summary_rows.append({
                "group": group,
                "comparison": comp_key,
                "status": comp.get("status"),
                "n_pseudosamples": n_ps,
                "n_significant": comp.get("n_significant", 0),
                "n_significant_local": comp.get(
                    "n_significant_local", comp.get("n_significant", 0)
                ),
                "n_significant_global": comp.get(
                    "n_significant_global", comp.get("n_significant", 0)
                ),
                "n_up": comp.get("n_up", 0),
                "n_up_local": comp.get("n_up_local", comp.get("n_up", 0)),
                "n_up_global": comp.get("n_up_global", comp.get("n_up", 0)),
                "n_down": comp.get("n_down", 0),
                "n_down_local": comp.get(
                    "n_down_local", comp.get("n_down", 0)
                ),
                "n_down_global": comp.get(
                    "n_down_global", comp.get("n_down", 0)
                ),
                "reason": comp.get("reason", ""),
            })
            records = comp.get("all_sig") or comp.get("top_genes") or []
            for rec in records:
                if not isinstance(rec, dict):
                    continue
                row = {
                    "group": group,
                    "comparison": comp_key,
                    "gene": rec.get("gene"),
                    "log2fc": rec.get("log2fc"),
                    "padj": rec.get("padj"),
                    "padj_local": rec.get("padj_local"),
                    "padj_global": rec.get("padj_global"),
                    "pvalue": rec.get("pvalue"),
                }
                for key, val in rec.items():
                    if key not in row:
                        row[key] = val
                pb_gene_rows.append(row)
    p = _write_tsv(output_dir / "scrna_pseudobulk_de_summary.tsv",
                   pb_summary_rows)
    if p:
        tables["pseudobulk_de_summary"] = p
    p = _write_tsv(output_dir / "scrna_pseudobulk_de_genes.tsv",
                   pb_gene_rows)
    if p:
        tables["pseudobulk_de_genes"] = p

    def _pathway_rows(container: dict, mode: str) -> list[dict]:
        rows = []
        for block_key, block in (container.get("per_cluster", {}) or {}).items():
            results = block.get("results", {}) or {}
            for db_name, terms in results.items():
                for term in terms or []:
                    if not isinstance(term, dict):
                        continue
                    raw_genes = _term_value(
                        term, "genes", "Genes", "lead_genes"
                    )
                    if isinstance(raw_genes, (list, tuple)):
                        genes_val = ";".join(
                            str(g) for g in raw_genes
                            if g is not None and str(g).strip()
                        )
                    else:
                        genes_val = raw_genes
                    rows.append({
                        "mode": mode,
                        "block": block_key,
                        "database": db_name,
                        "term": _term_value(term, "term", "Term"),
                        "adjusted_p": _term_value(
                            term, "padj", "adjusted_p", "Adjusted P-value",
                            "adj_p",
                        ),
                        "p_value": _term_value(
                            term, "pvalue", "p_value", "P-value",
                        ),
                        "overlap": _term_value(term, "overlap", "Overlap"),
                        "odds_ratio": _term_value(
                            term, "odds_ratio", "Odds Ratio"
                        ),
                        "combined_score": _term_value(
                            term, "combined_score", "Combined Score"
                        ),
                        "genes": genes_val,
                    })
        return rows

    pathway_rows = []
    pathway_rows.extend(_pathway_rows(findings.get("pathways") or {},
                                      "cluster_markers"))
    pathway_rows.extend(_pathway_rows(findings.get("pseudobulk_pathways") or {},
                                      "pseudobulk_de"))
    p = _write_tsv(output_dir / "scrna_pathway_enrichment.tsv",
                   pathway_rows)
    if p:
        tables["pathway_enrichment"] = p

    ccc = findings.get("cell_communication") or {}
    cc_rows = []
    for rec in ccc.get("top_interactions", []) or []:
        if isinstance(rec, dict):
            cc_rows.append(rec)
    p = _write_tsv(output_dir / "scrna_cellcomm_interactions.tsv", cc_rows)
    if p:
        tables["cellcomm_interactions"] = p

    traj = findings.get("trajectory") or {}
    paga = traj.get("paga", {}) or {}
    paga_rows = []
    for edge, val in (paga.get("top_connections", {}) or {}).items():
        if "->" in str(edge):
            source, target = str(edge).split("->", 1)
        elif "→" in str(edge):
            source, target = str(edge).split("→", 1)
        else:
            source, target = "", ""
        paga_rows.append({
            "edge": edge,
            "source": source.strip(),
            "target": target.strip(),
            "connectivity": val,
            "strong_threshold": paga.get("strong_threshold"),
            "is_strong": (
                _finite_number(val) is not None
                and val >= (
                    0.05 if _finite_number(paga.get("strong_threshold")) is None
                    else float(paga.get("strong_threshold"))
                )
            ),
        })
    p = _write_tsv(output_dir / "scrna_paga_connections.tsv", paga_rows)
    if p:
        tables["paga_connections"] = p

    pt = traj.get("pseudotime", {}) or {}
    pt_rows = [
        {"group": group, "mean_dpt": val}
        for group, val in (pt.get("pseudotime_by_group", {}) or {}).items()
    ]
    pt_rows.sort(key=lambda r: (
        _finite_number(r["mean_dpt"]) is None,
        _finite_number(r["mean_dpt"]) or 0.0,
    ))
    for i, row in enumerate(pt_rows, 1):
        row["rank"] = i
    p = _write_tsv(output_dir / "scrna_pseudotime_by_group.tsv", pt_rows)
    if p:
        tables["pseudotime_by_group"] = p

    if tables:
        findings["tables"] = tables
    return tables


