"""
Narrative kernel primitives.

This package defines the structured contract used by NarrativeAgent plugins.
The legacy report helpers still own modality-specific prose and figures; the
kernel wraps those outputs in auditable blocks before rendering.
"""

from aria.agents.narrative.types import Caveat, EvidenceItem, NarrativeBlock

__all__ = ["Caveat", "EvidenceItem", "NarrativeBlock"]
