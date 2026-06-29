"""Bulk RNA evidence adapter (ADR-057 S5).

Pure and deterministic — no LLM. Flattens the AUDITED bulk RNA findings (per
contrast: significant gene-level differential expression + ORA/GSEA pathway
terms) into ``EvidenceSignal`` items: ``entity`` is a real gene symbol or pathway
term, ``audited_node_ref`` is the run-ledger node that produced it
(``ledger://bulk/...``), and ``caveats_inherited`` carries the confounds that node
already flags (low replication, unmodeled batch). Honest-empty when there are no
bulk findings, and — when a ``run_ledger`` is supplied — it emits a layer's
signals only if that ledger node actually ran (no fabrication).
"""

from __future__ import annotations

from typing import Any

from aria.agents.narrative.run_ledger import _node_index

from ..types import EvidenceSignal

_DE_NODE = "ledger://bulk/differential_expression"
_ORA_NODE = "ledger://bulk/pathway_enrichment"
_MAX_GENES_PER_CONTRAST = 30
_MAX_TERMS_PER_CONTRAST = 10


def _bulk_findings(agent_results: dict) -> dict:
    bk = (agent_results or {}).get("bulk_rna_agent", {})
    findings = bk.get("findings", bk) if isinstance(bk, dict) else {}
    return findings if isinstance(findings, dict) else {}


def _node_ran(run_ledger: dict | None, node_id: str) -> bool:
    if run_ledger is None:
        return True  # caller asserts the inputs are audited
    node = _node_index(run_ledger).get(node_id)
    return bool(node) and node.get("status") == "ran"


def _num(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _direction(value: Any) -> str:
    v = _num(value)
    if v is None:
        return "na"
    return "up" if v > 0 else "down" if v < 0 else "na"


def _design_caveats(findings: dict, exp_ctx: dict | None) -> list[str]:
    caveats: list[str] = []
    if findings.get("low_power_warning") or findings.get("low_power"):
        caveats.append("low replication")
    design = ((exp_ctx or {}).get("design", {}) or {})
    covariates = [str(c).lower() for c in (design.get("covariates", []) or [])]
    if design.get("has_batch") and not any("batch" in c for c in covariates):
        caveats.append("batch effect not modeled")
    return caveats


def _contrast_label(contrast: dict, index: int) -> str:
    """A stable label for the contrast a signal was measured in (its context)."""
    for key in ("name", "contrast", "comparison", "label", "id"):
        value = contrast.get(key)
        if value:
            return str(value)
    return f"contrast_{index}"


def _ora_terms(contrast: dict) -> list[str]:
    pw = contrast.get("pathways")
    items: list[Any] = []
    if isinstance(pw, dict):
        for value in pw.values():
            if isinstance(value, list):
                items.extend(value)
    elif isinstance(pw, list):
        items = pw
    terms: list[str] = []
    for item in items:
        if isinstance(item, dict):
            name = item.get("term") or item.get("name") or item.get("pathway")
        else:
            name = item
        if name:
            terms.append(str(name))
        if len(terms) >= _MAX_TERMS_PER_CONTRAST:
            break
    return terms


def bulk_rna_evidence(
    agent_results: dict,
    run_ledger: dict | None = None,
    exp_ctx: dict | None = None,
) -> list[EvidenceSignal]:
    """Flatten audited bulk RNA findings into EvidenceSignal items."""
    findings = _bulk_findings(agent_results)
    if not findings:
        return []

    de_ran = _node_ran(run_ledger, _DE_NODE)
    ora_ran = _node_ran(run_ledger, _ORA_NODE)
    base_caveats = _design_caveats(findings, exp_ctx)

    signals: list[EvidenceSignal] = []
    seen: set[tuple[str, str, str]] = set()
    for c_index, contrast in enumerate(findings.get("contrasts") or []):
        if not isinstance(contrast, dict):
            continue
        label = _contrast_label(contrast, c_index)
        caveats = list(base_caveats)
        if (
            contrast.get("low_power") or contrast.get("low_power_warning")
        ) and "low replication" not in caveats:
            caveats.append("low replication")

        if de_ran:
            n_genes = 0
            for gene in contrast.get("top_genes") or []:
                if not isinstance(gene, dict):
                    continue
                symbol = gene.get("symbol") or gene.get("gene")
                lfc = gene.get("log2fc")
                if not symbol or lfc is None:
                    continue
                key = ("gene", str(symbol), label)
                if key in seen:
                    continue
                seen.add(key)
                signals.append(
                    EvidenceSignal(
                        entity=str(symbol),
                        entity_kind="gene",
                        modality="bulk_RNA",
                        measure="log2fc",
                        audited_node_ref=_DE_NODE,
                        value=_num(lfc),
                        direction=_direction(lfc),
                        caveats_inherited=list(caveats),
                        context=label,
                    )
                )
                n_genes += 1
                if n_genes >= _MAX_GENES_PER_CONTRAST:
                    break

        if ora_ran:
            for term in _ora_terms(contrast):
                key = ("pathway", term, label)
                if key in seen:
                    continue
                seen.add(key)
                signals.append(
                    EvidenceSignal(
                        entity=term,
                        entity_kind="pathway",
                        modality="bulk_RNA",
                        measure="ora",
                        audited_node_ref=_ORA_NODE,
                        direction="na",
                        caveats_inherited=list(caveats),
                        context=label,
                    )
                )
    return signals
