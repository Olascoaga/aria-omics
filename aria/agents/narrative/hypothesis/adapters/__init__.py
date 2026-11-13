"""Per-modality evidence adapters for the SPECULATIVE tier (ADR-057 S5-S8).

Each adapter is a pure, deterministic function that flattens one modality's
AUDITED artifacts into the modality-agnostic ``EvidenceSignal`` vocabulary the
HypothesisAgent reasons over. The grounding verifier lives one level up; the
adapter's job is to emit only real entities, each tagged with the run-ledger node
that produced it and the confounds that node already carries. V48 integration
becomes a 5th adapter (S10) emitting the same ``EvidenceSignal``.
"""

from __future__ import annotations

from ..caveats import CROSS_MODAL_UNPAIRED_COHORT, HIC_LOOP_ASSOCIATIVE
from .bulk_atac import bulk_atac_evidence
from .bulk_rna import bulk_rna_evidence
from .cross_modal import (
    CROSS_MODAL_NODE,
    SAMPLE_LINK_MEASURE,
    TF_CONCORDANCE_MEASURE,
    chromatin_tf_lineage_key,
    cross_modal_evidence,
    tag_chromatin_tf_lineage,
)
from .hic import HIC_LOOP_MEASURE, HIC_LOOPS_NODE, hic_evidence
from .integration import PAIRED_CIS_REGULATORY_KEY, integration_evidence
from .scatac import scatac_evidence
from .scrna import scrna_evidence

__all__ = [
    "bulk_rna_evidence",
    "scrna_evidence",
    "bulk_atac_evidence",
    "scatac_evidence",
    "hic_evidence",
    "HIC_LOOPS_NODE",
    "HIC_LOOP_MEASURE",
    "HIC_LOOP_ASSOCIATIVE",
    "integration_evidence",
    "PAIRED_CIS_REGULATORY_KEY",
    "cross_modal_evidence",
    "tag_chromatin_tf_lineage",
    "chromatin_tf_lineage_key",
    "CROSS_MODAL_NODE",
    "SAMPLE_LINK_MEASURE",
    "TF_CONCORDANCE_MEASURE",
    "CROSS_MODAL_UNPAIRED_COHORT",
]
