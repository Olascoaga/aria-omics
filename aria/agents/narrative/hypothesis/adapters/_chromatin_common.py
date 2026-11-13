"""Shared chromatin adapter helpers (ADR-057 S7-S8; ADR-078 footprinting).

Bulk ATAC and scATAC both interpret accessibility through the same audited
TF-motif enrichment (``motifs.per_group[*].top_motifs``) and, since ADR-078, the
same replicate-inferential TOBIAS differential footprinting. This module factors
the common pieces so the two adapters stay behavior-identical and DRY.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

from aria.agents.narrative.run_ledger import _node_index

from ..caveats import (
    BATCH,
    FOOTPRINT_ASSOCIATIVE,
    LOW_REPLICATION,
    PARALOG_MOTIF_AMBIGUITY,
)
from ..types import EvidenceSignal

FOOTPRINT_NODE = "ledger://chromatin/differential_tf_footprinting"


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
    group is recorded in ``measure``/``context`` and the log2 enrichment in
    ``value``. Dedups by (TF name, group): the SAME TF enriched in two clusters
    is two distinct, independent signals — not collapsed (H4). Shared verbatim by
    the bulk ATAC (S7) and scATAC (S8) adapters.
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
            name_key = str(name).strip().lower()
            key = (name_key, str(group))
            if not name_key or key in seen:
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
                    context=str(group),
                )
            )
            if len(out) >= max_motifs:
                return out
    return out


def _footprint_tf_symbols(name: Any) -> tuple[list[str], bool]:
    """Split a JASPAR/TOBIAS motif name into clean TF symbols.

    ``"FOSL2::JUND(var.2)"`` -> ``(["FOSL2", "JUND"], True)``; the boolean flags a
    ``::`` dimer, whose regulator identity is ambiguous between its subunits.
    """
    raw = str(name or "").split("(")[0].strip()
    if not raw:
        return [], False
    is_dimer = "::" in raw
    parts = raw.split("::") if is_dimer else [raw]
    return [p.strip() for p in parts if p.strip()], is_dimer


def _footprint_rows_from_table(path: Any, alpha: float) -> list[dict]:
    """All BH-significant rows of the persisted replicate-inference TSV."""
    if not path or not Path(str(path)).is_file():
        return []
    rows: list[dict] = []
    try:
        with open(str(path), encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                padj = num(row.get("padj"))
                if padj is None or not math.isfinite(padj) or padj > alpha:
                    continue
                rows.append({"tf": row.get("tf") or row.get("motif_id"),
                             "change": num(row.get("change")), "padj": padj})
    except (OSError, csv.Error, UnicodeDecodeError):
        return []
    return rows


def footprint_signals(
    footprinting: dict,
    *,
    modality: str,
    node: str = FOOTPRINT_NODE,
    extra_caveats: list[str] | None = None,
    max_tfs: int = 20,
) -> list[EvidenceSignal]:
    """Flatten BH-significant replicate footprint results into tf_motif signals.

    STRICT INFERENTIAL GATE (ADR-078): only a run whose replicate inference
    succeeded emits signals, and only motifs with ``padj <= alpha`` across motifs.
    Descriptive / pooled / batch-confounded footprinting never becomes evidence.
    Every signal inherits ``FOOTPRINT_ASSOCIATIVE`` (occupancy is not regulation);
    a ``::`` dimer is split into its subunits, each carrying
    ``PARALOG_MOTIF_AMBIGUITY``.
    """
    if not (isinstance(footprinting, dict) and footprinting.get("ran") is True):
        return []
    summary = footprinting.get("differential_summary") or {}
    if not isinstance(summary, dict):
        return []
    inference = summary.get("inference") or {}
    if not isinstance(inference, dict) or inference.get("status") != "success":
        return []

    ga = str(footprinting.get("group_a") or "group_a")
    gb = str(footprinting.get("group_b") or "group_b")
    alpha = num(inference.get("alpha")) or num(summary.get("alpha")) or 0.05
    context = f"{ga}_vs_{gb}"

    caveats = [FOOTPRINT_ASSOCIATIVE]
    if inference.get("low_power_warning") and LOW_REPLICATION not in caveats:
        caveats.append(LOW_REPLICATION)
    if ((inference.get("batch_adjusted")
         or "batch" in str(inference.get("test", "")).lower())
            and BATCH not in caveats):
        caveats.append(BATCH)
    for code in extra_caveats or []:
        if code not in caveats:
            caveats.append(code)

    rows = _footprint_rows_from_table(
        footprinting.get("results_table") or summary.get("results_table"), alpha)
    if not rows:
        rows = [r for r in ((summary.get(f"top_toward_{ga}") or [])
                            + (summary.get(f"top_toward_{gb}") or []))
                if isinstance(r, dict)]
    candidates = []
    for row in rows:
        padj = num(row.get("padj"))
        if padj is None or not math.isfinite(padj) or padj > alpha:
            continue
        change = num(row.get("change"))
        candidates.append((padj, -abs(change or 0.0), row.get("tf"), change))
    candidates.sort(key=lambda item: (item[0], item[1]))

    out: list[EvidenceSignal] = []
    seen: set[tuple[str, str]] = set()
    for _padj, _neg_abs, name, change in candidates:
        symbols, is_dimer = _footprint_tf_symbols(name)
        row_caveats = list(caveats)
        if is_dimer and PARALOG_MOTIF_AMBIGUITY not in row_caveats:
            row_caveats.append(PARALOG_MOTIF_AMBIGUITY)
        for symbol in symbols:
            key = (symbol.lower(), context)
            if key in seen:
                continue
            seen.add(key)
            effect = change or 0.0
            out.append(
                EvidenceSignal(
                    entity=symbol,
                    entity_kind="tf_motif",
                    modality=modality,  # type: ignore[arg-type]
                    measure=f"footprint_score({ga}_vs_{gb})",
                    audited_node_ref=node,
                    value=change,
                    direction=("up" if effect > 0
                               else "down" if effect < 0 else "na"),
                    caveats_inherited=list(row_caveats),
                    context=context,
                )
            )
            if len(out) >= max_tfs:
                return out
    return out
