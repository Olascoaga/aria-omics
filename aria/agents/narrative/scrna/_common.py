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
    "_top_pathway_terms", "_term_value", "_comp_fdr_strategy",
    "_primary_de_counts", "_record_padj", "_finite_number",
]

import math


def _term_value(term: dict, *keys, default=""):
    """Return the first present key's value from a term dict. Leaf helper shared
    by the pathway selectors here and the supplementary-table export."""
    for key in keys:
        if key in term:
            return term.get(key)
    return default


def _finite_number(value) -> float | None:
    """Return ``value`` as a finite float, else ``None`` (bools excluded)."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def _comp_fdr_strategy(comp: dict, pb: dict | None = None) -> str | None:
    """Primary FDR family of a pseudobulk comparison.

    Returns ``"per_cluster"`` / ``"global"`` when the script declared it (on the
    comparison or in ``pb["multiple_testing"]``), else ``None`` (legacy results,
    which used global pooled BH).
    """
    strategy = (comp or {}).get("fdr_strategy")
    if not strategy and pb:
        strategy = ((pb.get("multiple_testing") or {}).get("fdr_strategy"))
    return str(strategy) if strategy else None


def _primary_de_counts(comp: dict, pb: dict | None = None) -> dict:
    """Significance counts under the PRIMARY FDR family (Wave 4A).

    ``top_genes`` / ``all_sig`` and the per-group bar chart follow the primary
    family declared by ``fdr_strategy`` (``n_significant``/``n_up``/``n_down``
    in rna_pseudobulk_de). Tables and prose must report the same family so they
    never print e.g. ``0 up`` next to a list of up-regulated top genes. Legacy
    results without a strategy keep the historical global-first lookup.
    """
    comp = comp or {}
    strategy = _comp_fdr_strategy(comp, pb)

    def _pick(*keys):
        for key in keys:
            if comp.get(key) is not None:
                return comp.get(key)
        return 0

    if strategy == "per_cluster":
        n_sig = _pick("n_significant", "n_significant_local")
        n_up = _pick("n_up", "n_up_local")
        n_down = _pick("n_down", "n_down_local")
        secondary = comp.get("n_significant_global")
        label, secondary_label = "per-cluster FDR", "global"
    elif strategy == "global":
        n_sig = _pick("n_significant", "n_significant_global")
        n_up = _pick("n_up", "n_up_global")
        n_down = _pick("n_down", "n_down_global")
        secondary = comp.get("n_significant_local")
        label, secondary_label = "global-FDR", "local"
    else:
        n_sig = _pick("n_significant_global", "n_significant")
        n_up = _pick("n_up_global", "n_up")
        n_down = _pick("n_down_global", "n_down")
        secondary = comp.get("n_significant_local", n_sig)
        label, secondary_label = "global-FDR", "local"
    return {
        "strategy": strategy or "global",
        "n_significant": n_sig,
        "n_up": n_up,
        "n_down": n_down,
        "n_secondary": secondary,
        "label": label,
        "secondary_label": secondary_label,
    }


def _record_padj(rec: dict, strategy: str | None = None):
    """Adjusted p of a gene record under the primary FDR family.

    ``rec["padj"]`` is the primary column in rna_pseudobulk_de; legacy records
    without a declared strategy keep the global-first lookup.
    """
    if strategy == "per_cluster":
        keys = ("padj", "padj_local", "padj_global")
    elif strategy == "global":
        keys = ("padj_global", "padj", "padj_local")
    else:
        keys = ("padj_global", "padj", "padj_local")
    for key in keys:
        if rec.get(key) is not None:
            return rec.get(key)
    return None


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


def _gene_brief(rec: dict, strategy: str | None = None) -> str:
    gene = _gene_name(rec)
    lfc = rec.get("log2fc", rec.get("log2FoldChange"))
    padj = _record_padj(rec, strategy)
    details = []
    if _finite_number(lfc) is not None:
        details.append(f"log2FC={lfc:+.2f}")
    if _finite_number(padj) is not None:
        details.append(f"FDR={_fmt_stat(padj)}")
    return f"{gene} ({', '.join(details)})" if details else gene


def _top_directional_genes(comp: dict, direction: str,
                           limit: int = 3) -> list[str]:
    strategy = _comp_fdr_strategy(comp)
    records = comp.get("top_genes") or comp.get("all_sig") or []
    if direction == "up":
        rows = [r for r in records if isinstance(r, dict)
                and r.get("log2fc", 0) > 0]
        rows.sort(key=lambda r: r.get("log2fc", 0), reverse=True)
    else:
        rows = [r for r in records if isinstance(r, dict)
                and r.get("log2fc", 0) < 0]
        rows.sort(key=lambda r: r.get("log2fc", 0))
    return [_gene_brief(r, strategy) for r in rows[:limit]]


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
    # Wave 4A: no unanchored substring fallback — "T" must not match
    # "CD8_T::..." and cluster "1" must not match "10::...". Only accept keys
    # whose delimiter-split tokens equal the group and comparison exactly.
    for key, block in per_cluster.items():
        key_s = str(key)
        for sep in ("::", "__", " "):
            if sep in key_s:
                head, _, tail = key_s.partition(sep)
                if head == group_s and tail == comp_s:
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
                row, "padj", "adjusted_p", "Adjusted P-value", "adj_p",
                default=None,
            )
            label = f"{db_name}: {term}"
            if _finite_number(padj) is not None:
                label += f" (FDR={_fmt_stat(padj)})"
            terms.append(label)
            if len(terms) >= limit:
                return terms
    return terms


