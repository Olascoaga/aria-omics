"""Shared chromatin adapter helpers (ADR-057 S7-S8).

Bulk ATAC and scATAC both interpret accessibility through the same audited
TF-motif enrichment (``motifs.per_group[*].top_motifs``). This module factors the
common pieces so the two adapters stay behavior-identical and DRY.
"""

from __future__ import annotations

from typing import Any

from aria.agents.narrative.run_ledger import _node_index

from ..types import EvidenceSignal


def node_ran(run_ledger: dict | None, node_id: str) -> bool:
    if run_ledger is None:
        return True
    node = _node_index(run_ledger).get(node_id)
    return bool(node) and node.get("status") == "ran"


def num(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def motif_signals(
    motifs: dict,
    *,
    modality: str,
    node: str,
    base_caveats: list[str],
    max_motifs: int,
    seen: set[str] | None = None,
) -> list[EvidenceSignal]:
    """Flatten audited ``top_motifs`` per group into tf_motif EvidenceSignals.

    Keeps only motifs enriched at the run's padj threshold; the enriched peak
    group is recorded in ``measure`` and the log2 enrichment in ``value``. Dedups
    by TF name. Shared verbatim by the bulk ATAC (S7) and scATAC (S8) adapters.
    """
    out: list[EvidenceSignal] = []
    if not isinstance(motifs, dict):
        return out
    padj_max = num(motifs.get("padj_max"))
    seen = seen if seen is not None else set()
    for group, gdata in (motifs.get("per_group") or {}).items():
        if not isinstance(gdata, dict):
            continue
        for motif in gdata.get("top_motifs") or []:
            if not isinstance(motif, dict):
                continue
            name = motif.get("name")
            if not name:
                continue
            padj = num(motif.get("padj"))
            if padj is not None and padj_max is not None and padj >= padj_max:
                continue
            key = str(name).strip().lower()
            if not key or key in seen:
                continue
            seen.add(key)
            out.append(
                EvidenceSignal(
                    entity=str(name),
                    entity_kind="tf_motif",
                    modality=modality,  # type: ignore[arg-type]
                    measure=f"motif_enrich({group})",
                    audited_node_ref=node,
                    value=num(motif.get("log2_enrichment")),
                    direction="na",
                    caveats_inherited=list(base_caveats),
                )
            )
            if len(out) >= max_motifs:
                return out
    return out
