"""scRNA narrative leaf helpers: findings normalisation, formatting, selectors.

Shared by the text/tables/figures modules of this subpackage. Extracted verbatim
from the former monolithic aria/agents/_narrative_scrna.py (A7); behavior is
unchanged and pinned by tests/test_narrative_scrna_contract.py.
"""
from __future__ import annotations

__all__ = [
    "unwrap_scrna_findings", "_fmt_int", "_fdr_primary_clause",
    "_lfc_shrinkage_clause", "_fmt_stat", "_group_label", "_label_cell_type",
    "_annotation_state", "_top_de_blocks", "_top_pathway_blocks", "_gene_name",
    "_gene_brief", "_top_directional_genes", "_find_pathway_block",
    "_top_pathway_terms", "_term_value",
]


def _term_value(term: dict, *keys, default=""):
    """Return the first present key's value from a term dict. Leaf helper shared
    by the pathway selectors here and the supplementary-table export."""
    for key in keys:
        if key in term:
            return term.get(key)
    return default


def unwrap_scrna_findings(agent_result: dict) -> dict:
    """
    Return the scRNA findings dict from a scrna_agent envelope, robust to
    both shapes that exist in the codebase:

        - Adapter / multimodal-wrapped:  {findings: {scRNA: {findings: {...}}}}
        - scrna_agent.run() direct:      {findings: {qc, clustering, ...}}

    Without this helper the TUI / Orchestrator path silently returns empty
    findings (the inner scRNA wrapper does not exist on the direct emit).
    """
    f = agent_result.get("findings", {}) or {}
    wrapped = (f.get("scRNA", {}) or {}).get("findings", {}) or {}
    return wrapped or f


# ── Text summaries ────────────────────────────────────────────────────────

def _fmt_int(value) -> str:
    try:
        return f"{int(value):,}"
    except Exception:
        return str(value)


def _fdr_primary_clause(pb: dict) -> str:
    """Name the adjusted-p column that defined significance, per fdr_strategy.

    Legacy results without an explicit strategy used global pooled BH.
    """
    strategy = ((pb or {}).get("multiple_testing", {}) or {}).get(
        "fdr_strategy", "global"
    )
    return "padj_local (per-cluster BH)" if strategy == "per_cluster" else "padj_global"


def _lfc_shrinkage_clause(pb: dict) -> str:
    """Disclose apeGLM LFC shrinkage when it was requested (C4)."""
    shrink = (pb or {}).get("lfc_shrinkage") or {}
    if not shrink.get("requested"):
        return ""
    return (
        "Reported log2 fold changes are apeGLM-shrunken estimates (pydeseq2 "
        "lfc_shrink); the effect-size threshold is applied to the shrunken "
        "value while p-values are unchanged, and the unshrunken MLE is kept as "
        "log2fc_raw. "
    )


def _fmt_stat(value) -> str:
    """Compact numeric display that preserves very small nonzero values."""
    if not isinstance(value, (int, float)):
        return str(value)
    if value == 0:
        return "0"
    value = float(value)
    if abs(value) < 1e-3:
        return f"{value:.2e}"
    return f"{value:.4g}"


def _group_label(groupby: str | None, n: int | None = None) -> str:
    """Human-readable label for an obs grouping column."""
    if not groupby:
        return "groups"
    if groupby in {"cell_type", "celltype", "cell_type_celltypist"}:
        singular = "cell type"
    elif groupby == "leiden":
        singular = "Leiden cluster"
    else:
        singular = f"{groupby} group"
    return singular if n == 1 else f"{singular}s"


def _label_cell_type(value) -> str:
    if isinstance(value, dict):
        return (value.get("cell_type")
                or value.get("celltypist_label")
                or value.get("label")
                or "")
    return str(value) if value else ""


def _annotation_state(findings: dict) -> dict:
    ct_block = findings.get("cell_types") or {}
    labels = [
        _label_cell_type(v)
        for v in (ct_block.get("cell_types", {}) or {}).values()
    ]
    labels = [x for x in labels if x]
    invalid = {"annotation_failed", "failed", "unknown", "nan", "none"}
    valid = [
        x for x in labels
        if x.strip().lower() not in invalid
        and not x.strip().lower().startswith("unresolved cluster")
    ]
    source = "unknown"
    for v in (ct_block.get("cell_types", {}) or {}).values():
        if isinstance(v, dict) and v.get("annotation_source"):
            source = str(v.get("annotation_source"))
            break
    label_col = ct_block.get("label_col")
    return {
        "has_valid": bool(valid),
        "labels": valid,
        "n_unique": len(set(valid)),
        "source": source,
        "label_col": label_col,
        "is_marker_fallback": source in {
            "marker_fallback", "unresolved_marker_fallback",
        },
    }


def _top_de_blocks(pb: dict, limit: int = 5) -> list[tuple[str, str, dict]]:
    rows = []
    for group, info in (pb.get("per_group", {}) or {}).items():
        for comp_key, comp in (info.get("per_comparison", {}) or {}).items():
            if comp.get("status") == "success":
                rows.append((str(group), str(comp_key), comp))
    rows.sort(key=lambda row: row[2].get("n_significant", 0), reverse=True)
    return rows[:limit]


def _top_pathway_blocks(pwp: dict, limit: int = 3) -> list[tuple[str, dict]]:
    blocks = list((pwp.get("per_cluster", {}) or {}).items())
    blocks.sort(key=lambda kv: kv[1].get("n_significant", 0), reverse=True)
    return blocks[:limit]


def _gene_name(rec: dict) -> str:
    return str(rec.get("symbol") or rec.get("gene") or rec.get("name") or "?")


def _gene_brief(rec: dict) -> str:
    gene = _gene_name(rec)
    lfc = rec.get("log2fc", rec.get("log2FoldChange"))
    padj = rec.get("padj_global", rec.get("padj", rec.get("padj_local")))
    details = []
    if isinstance(lfc, (int, float)):
        details.append(f"log2FC={lfc:+.2f}")
    if isinstance(padj, (int, float)):
        details.append(f"FDR={_fmt_stat(padj)}")
    return f"{gene} ({', '.join(details)})" if details else gene


def _top_directional_genes(comp: dict, direction: str,
                           limit: int = 3) -> list[str]:
    records = comp.get("top_genes") or comp.get("all_sig") or []
    if direction == "up":
        rows = [r for r in records if isinstance(r, dict)
                and r.get("log2fc", 0) > 0]
        rows.sort(key=lambda r: r.get("log2fc", 0), reverse=True)
    else:
        rows = [r for r in records if isinstance(r, dict)
                and r.get("log2fc", 0) < 0]
        rows.sort(key=lambda r: r.get("log2fc", 0))
    return [_gene_brief(r) for r in rows[:limit]]


def _find_pathway_block(pwp: dict, group: str, comp_key: str) -> tuple[str, dict]:
    per_cluster = pwp.get("per_cluster", {}) or {}
    if not per_cluster:
        return "", {}
    group_s = str(group)
    comp_s = str(comp_key)
    candidates = [
        f"{group_s}::{comp_s}",
        f"{group_s}__{comp_s}",
        f"{group_s} {comp_s}",
        f"{group_s}_{comp_s}",
    ]
    for key in candidates:
        if key in per_cluster:
            return key, per_cluster[key]
    for key, block in per_cluster.items():
        key_s = str(key)
        if group_s in key_s and comp_s in key_s:
            return key_s, block
    return "", {}


def _top_pathway_terms(block: dict, limit: int = 3) -> list[str]:
    terms = []
    for db_name, rows in (block.get("results", {}) or {}).items():
        for row in rows or []:
            if not isinstance(row, dict):
                continue
            term = row.get("term") or row.get("Term")
            if not term:
                continue
            padj = _term_value(
                row, "adjusted_p", "Adjusted P-value", "adj_p", "padj",
                default=None,
            )
            label = f"{db_name}: {term}"
            if isinstance(padj, (int, float)):
                label += f" (FDR={_fmt_stat(padj)})"
            terms.append(label)
            if len(terms) >= limit:
                return terms
    return terms


