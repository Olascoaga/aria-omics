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

from .devils_advocate import (
    check_devils_advocate,
    declared_confounds,
    visible_confounds,
)
from .gates import GateResult, check_falsifiability, check_language
from .grounding import (
    GroundingResult,
    build_evidence_index,
    verify_hypothesis_grounding,
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
]
