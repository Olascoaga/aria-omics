"""SPECULATIVE hypothesis layer (ADR-057).

Creativity in the reasoning, zero invention in the facts. This package holds the
modality-agnostic Evidence Interface (``EvidenceSignal``), the hypothesis schema
(``Hypothesis`` / ``DiscriminatingExperiment``), and the grounding verifier that
rejects any hypothesis naming an entity absent from the audited evidence.

S1 scope: types + grounding verifier. The falsifiability gate (S2),
devils_advocate (S3), ledger quarantine (S4), and per-modality evidence adapters
(S5-S8) land in later slices. See
``memory/audit/ARIA_PLAN_HypothesisAgent_2026-06-21.md``.
"""

from __future__ import annotations

from .adapters import (
    bulk_atac_evidence,
    bulk_rna_evidence,
    scatac_evidence,
    scrna_evidence,
)
from .devils_advocate import (
    check_devils_advocate,
    declared_confounds,
    visible_confounds,
)
from .proposer_llm import (
    LLMProposer,
    build_proposer_prompt,
    parse_hypotheses,
)
from .ranking import rank_hypotheses
from .report_section import (
    build_speculative_section,
    gather_evidence,
    render_speculative_section_html,
)
from .gates import GateResult, check_falsifiability, check_language
from .grounding import (
    GroundingResult,
    build_evidence_index,
    verify_hypothesis_grounding,
)
from .quarantine import (
    SPECULATIVE_TIER,
    SpeculativePromotionError,
    assert_no_speculative_promotion,
    is_hypothesis_node,
    quarantine_hypotheses,
    quarantine_node_id,
)
from .types import DiscriminatingExperiment, EvidenceSignal, Hypothesis

__all__ = [
    "EvidenceSignal",
    "Hypothesis",
    "DiscriminatingExperiment",
    "GroundingResult",
    "build_evidence_index",
    "verify_hypothesis_grounding",
    "GateResult",
    "check_falsifiability",
    "check_language",
    "check_devils_advocate",
    "visible_confounds",
    "declared_confounds",
    "SPECULATIVE_TIER",
    "SpeculativePromotionError",
    "assert_no_speculative_promotion",
    "is_hypothesis_node",
    "quarantine_hypotheses",
    "quarantine_node_id",
    "bulk_rna_evidence",
    "scrna_evidence",
    "bulk_atac_evidence",
    "scatac_evidence",
    "LLMProposer",
    "build_proposer_prompt",
    "parse_hypotheses",
    "rank_hypotheses",
    "build_speculative_section",
    "gather_evidence",
    "render_speculative_section_html",
]
