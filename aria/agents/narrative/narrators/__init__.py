"""Built-in modality narrators."""

from aria.agents.narrative.narrators.scrna import ScrnaNarrator
from aria.agents.narrative.narrators.bulk_rna import BulkRnaNarrator
from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
from aria.agents.narrative.narrators.hic import HiCNarrator
from aria.agents.narrative.narrators.integration import IntegrationNarrator

__all__ = [
    "ScrnaNarrator",
    "BulkRnaNarrator",
    "ChromatinNarrator",
    "HiCNarrator",
    "IntegrationNarrator",
]
