"""Bulk ATAC evidence adapter (ADR-057 S7).

Pure and deterministic — no LLM. The hypothesis-relevant, inline, NAMED entities
of bulk ATAC are the transcription factors whose motifs are over-represented in a
condition's differentially accessible peaks. This adapter flattens the audited
TF-motif enrichment into ``EvidenceSignal`` items: ``entity`` is a real TF name,
``measure`` records the peak group it was enriched in, ``value`` is the log2
enrichment, and ``audited_node_ref`` is ``ledger://chromatin/motif_enrichment``
(bulk ATAC routes through the ``chromatin`` ledger group).

Every signal inherits the two permanent bulk ATAC caveats: an enriched motif is
an associative database match, NOT evidence of TF binding or activity; and bulk
ATAC replication is typically low, so effects are directional rather than
FDR-calibrated when the DA layer flagged low power. Raw DA peak coordinates and
the CSV-backed peak-annotation / ORA gene layers are intentionally not emitted as
entities here (coordinates are not groundable names; gene/pathway names live only
in the artifact CSVs). BH-significant replicate TF footprints (ADR-078) are
emitted as tf_motif signals carrying ``footprint_associative``. A layer is emitted
only when its ledger node actually ran.
"""

from __future__ import annotations

from aria.agents.narrative.run_ledger import _chromatin_findings

from ..caveats import LOW_REPLICATION, MOTIF_NOT_BINDING
from ..types import EvidenceSignal
from ._chromatin_common import (
    FOOTPRINT_NODE,
    footprint_signals,
    motif_signals,
    node_ran,
)

_MOTIF_NODE = "ledger://chromatin/motif_enrichment"
_MAX_MOTIFS = 40

# H16: structured caveat codes (human gloss lives in caveats.py).
_ASSOCIATIVE_CAVEAT = MOTIF_NOT_BINDING
_LOW_REPLICATION_CAVEAT = LOW_REPLICATION


def _low_power(findings: dict) -> bool:
    da = findings.get("differential_accessibility") or {}
    comparisons = da.get("comparisons") or []
    return any(
        isinstance(c, dict) and c.get("low_power_warning") for c in comparisons
    )


def _bulk_chromatin_modality(findings: dict) -> str:
    for key in ("preprocessing", "qc", "peaks", "peak_counts",
                "differential_accessibility", "motifs", "footprinting"):
        sub = findings.get(key)
        if isinstance(sub, dict):
            dt = sub.get("data_type")
            if dt in {"bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"}:
                return str(dt)
    return "bulk_ATAC"


def bulk_atac_evidence(
    agent_results: dict,
    run_ledger: dict | None = None,
    exp_ctx: dict | None = None,
) -> list[EvidenceSignal]:
    """Flatten audited bulk chromatin (`bulk_ATAC`, `ChIP`, `CUT_AND_RUN`,
    `CUT_AND_TAG`) TF-motif enrichment and footprints into EvidenceSignal items."""
    findings = _chromatin_findings(agent_results)
    if not findings:
        return []

    modality = _bulk_chromatin_modality(findings)
    signals: list[EvidenceSignal] = []
    low_power = _low_power(findings)

    # 1. TF-motif enrichment over the DA/DB peak sets.
    motifs = findings.get("motifs") or {}
    if motifs.get("ran") and node_ran(run_ledger, _MOTIF_NODE):
        caveats = [_ASSOCIATIVE_CAVEAT]
        if low_power:
            caveats.append(_LOW_REPLICATION_CAVEAT)
        signals += motif_signals(
            motifs,
            modality=modality,
            node=_MOTIF_NODE,
            base_caveats=caveats,
            max_motifs=_MAX_MOTIFS,
        )

    # 2. ADR-078: replicate-inferential differential TF footprinting. Independent
    # of the motif layer (a run without motif enrichment still grounds its
    # BH-significant footprints); descriptive footprinting emits nothing.
    footprinting = findings.get("footprinting") or {}
    if footprinting.get("ran") and node_ran(run_ledger, FOOTPRINT_NODE):
        signals += footprint_signals(
            footprinting,
            modality=modality,
            node=FOOTPRINT_NODE,
            extra_caveats=[LOW_REPLICATION] if low_power else [],
        )

    return signals
