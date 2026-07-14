"""Composed report-builder mixin used by :class:`NarrativeAgent`."""

from aria.agents.narrative.reporting.executive import ExecutiveSummaryMixin
from aria.agents.narrative.reporting.methodology import MethodologyMixin
from aria.agents.narrative.reporting.render import ReportRenderMixin
from aria.agents.narrative.reporting.sections import ReportSectionsMixin


class ReportBuilderMixin(
    ReportRenderMixin,
    ExecutiveSummaryMixin,
    MethodologyMixin,
    ReportSectionsMixin,
):
    """Compose report-building concerns while preserving the ``self`` API."""


__all__ = ["ReportBuilderMixin"]
