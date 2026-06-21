"""Per-modality evidence adapters for the SPECULATIVE tier (ADR-057 S5-S8).

Each adapter is a pure, deterministic function that flattens one modality's
AUDITED artifacts into the modality-agnostic ``EvidenceSignal`` vocabulary the
HypothesisAgent reasons over. The grounding verifier lives one level up; the
adapter's job is to emit only real entities, each tagged with the run-ledger node
that produced it and the confounds that node already carries. V48 integration
becomes a 5th adapter (S10) emitting the same ``EvidenceSignal``.
"""

from __future__ import annotations

from .bulk_atac import bulk_atac_evidence
from .bulk_rna import bulk_rna_evidence
from .scrna import scrna_evidence

__all__ = ["bulk_rna_evidence", "scrna_evidence", "bulk_atac_evidence"]
