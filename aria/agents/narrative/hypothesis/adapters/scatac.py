"""scATAC evidence adapter (ADR-057 S8) — the last single-modality adapter.

Pure and deterministic — no LLM. Flattens the audited scATAC findings into
``EvidenceSignal`` items across its named, inline layers:

  - per-cluster TF-motif enrichment -> tf_motif signals (the most robust scATAC
    interpretive layer; shared with bulk ATAC via ``_chromatin_common``).
  - regulatory gene-activity scores  -> gene signals, caveated as moderate and
    NOT RNA expression.
  - peak-to-gene links               -> gene signals, caveated as associative
    correlation, not regulatory mechanism.

Per-cluster Wilcoxon DA markers are peak coordinates (not groundable names, like
bulk ATAC) and are intentionally not emitted as entities. Each layer is emitted
only when its run-ledger node ran; every signal carries the permanent scientific
caveats so the devils_advocate gate forces the LLM to own them.
"""

from __future__ import annotations

from aria.agents.narrative.run_ledger import _chromatin_findings

from ..caveats import (
    GENE_ACTIVITY_PROXY,
    MOTIF_NOT_BINDING,
    PEAK2GENE_ASSOCIATIVE,
)
from ..types import EvidenceSignal
from ._chromatin_common import motif_signals, node_ran, num

_MOTIF_NODE = "ledger://chromatin/motif_enrichment"
_REG_NODE = "ledger://chromatin/regulatory_layers"
_MAX_MOTIFS = 40
_MAX_GENES = 20

# H16: caveats are structured codes (the gate's source of truth); the human gloss
# is rendered downstream from aria/agents/narrative/hypothesis/caveats.py.
_MOTIF_CAVEAT = MOTIF_NOT_BINDING
_GENE_ACTIVITY_CAVEAT = GENE_ACTIVITY_PROXY
_PEAK2GENE_CAVEAT = PEAK2GENE_ASSOCIATIVE


def _corr_direction(value) -> str:
    s = str(value or "").lower()
    if "pos" in s or "up" in s:
        return "up"
    if "neg" in s or "down" in s:
        return "down"
    return "na"


def scatac_evidence(
    agent_results: dict,
    run_ledger: dict | None = None,
    exp_ctx: dict | None = None,
) -> list[EvidenceSignal]:
    """Flatten audited scATAC findings into EvidenceSignal items."""
    findings = _chromatin_findings(agent_results)
    if not findings:
        return []

    signals: list[EvidenceSignal] = []

    # 1. per-cluster TF motif enrichment (most robust)
    motifs = findings.get("motifs") or {}
    if motifs.get("ran") and node_ran(run_ledger, _MOTIF_NODE):
        signals += motif_signals(
            motifs,
            modality="scATAC",
            node=_MOTIF_NODE,
            base_caveats=[_MOTIF_CAVEAT],
            max_motifs=_MAX_MOTIFS,
        )

    # 2. regulatory layers: gene-activity + peak-to-gene named genes
    regulatory = findings.get("regulatory") or {}
    if regulatory.get("ran") and node_ran(run_ledger, _REG_NODE):
        seen: set[tuple[str, str]] = set()

        def _add_gene(gene, measure, *, value=None, direction="na", caveat):
            ent = str(gene or "").strip()
            if not ent:
                return
            key = (measure, ent.lower())
            if key in seen:
                return
            seen.add(key)
            signals.append(
                EvidenceSignal(
                    entity=ent,
                    entity_kind="gene",
                    modality="scATAC",
                    measure=measure,
                    audited_node_ref=_REG_NODE,
                    value=value,
                    direction=direction,
                    caveats_inherited=[caveat],
                )
            )

        gene_scores = regulatory.get("gene_scores") or {}
        if gene_scores.get("ran"):
            n = 0
            for row in gene_scores.get("top_genes") or []:
                if not isinstance(row, dict):
                    continue
                _add_gene(
                    row.get("gene"), "gene_activity",
                    value=num(row.get("mean_gene_activity_score")),
                    caveat=_GENE_ACTIVITY_CAVEAT,
                )
                n += 1
                if n >= _MAX_GENES:
                    break

        peak_to_gene = regulatory.get("peak_to_gene") or {}
        if peak_to_gene.get("ran"):
            n = 0
            for row in peak_to_gene.get("top_links") or []:
                if not isinstance(row, dict):
                    continue
                _add_gene(
                    row.get("gene"), "peak2gene",
                    value=num(row.get("correlation")),
                    direction=_corr_direction(row.get("direction")),
                    caveat=_PEAK2GENE_CAVEAT,
                )
                n += 1
                if n >= _MAX_GENES:
                    break

    return signals
