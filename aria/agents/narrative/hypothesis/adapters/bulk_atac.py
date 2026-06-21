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
in the artifact CSVs). A layer is emitted only when its ledger node actually ran.
"""

from __future__ import annotations

from aria.agents.narrative.run_ledger import _chromatin_findings

from ..types import EvidenceSignal
from ._chromatin_common import motif_signals, node_ran

_MOTIF_NODE = "ledger://chromatin/motif_enrichment"
_MAX_MOTIFS = 40

_ASSOCIATIVE_CAVEAT = (
    "enriched motif is an associative database match, not evidence of TF "
    "binding or activity"
)
_LOW_REPLICATION_CAVEAT = (
    "low replication (bulk ATAC n is typically small; directional, not "
    "FDR-calibrated)"
)


def _low_power(findings: dict) -> bool:
    da = findings.get("differential_accessibility") or {}
    comparisons = da.get("comparisons") or []
    return any(
        isinstance(c, dict) and c.get("low_power_warning") for c in comparisons
    )


def bulk_atac_evidence(
    agent_results: dict,
    run_ledger: dict | None = None,
    exp_ctx: dict | None = None,
) -> list[EvidenceSignal]:
    """Flatten audited bulk ATAC TF-motif enrichment into EvidenceSignal items."""
    findings = _chromatin_findings(agent_results)
    if not findings:
        return []

    motifs = findings.get("motifs") or {}
    if not (motifs.get("ran") and node_ran(run_ledger, _MOTIF_NODE)):
        return []

    caveats = [_ASSOCIATIVE_CAVEAT]
    if _low_power(findings):
        caveats.append(_LOW_REPLICATION_CAVEAT)

    return motif_signals(
        motifs,
        modality="bulk_ATAC",
        node=_MOTIF_NODE,
        base_caveats=caveats,
        max_motifs=_MAX_MOTIFS,
    )
