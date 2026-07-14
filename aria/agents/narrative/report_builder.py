"""Compatibility facade for the concern-oriented report builder.

The implementation lives in :mod:`aria.agents.narrative.reporting`. Existing
imports keep using this module so ``NarrativeAgent`` and external test harnesses
retain the same contract.
"""

from aria.agents.narrative.reporting import (
    ReportBuilderMixin,
    _executive_summary_numbers,
    _normalize_exec_summary_number,
)

__all__ = [
    "ReportBuilderMixin",
    "_executive_summary_numbers",
    "_normalize_exec_summary_number",
]
