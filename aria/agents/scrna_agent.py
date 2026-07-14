"""ARIA scRNAAgent — compatibility facade.

The implementation moved to the cohesive ``aria.agents.scrna`` subpackage in A7
(a mixin split of the former 2.3k-line class into QC/clustering, annotation,
DE/pathway, pseudobulk and trajectory/cell-communication concerns). This module
is kept as the stable import path: consumers (orchestrator dynamic import,
tests) still do ``from aria.agents.scrna_agent import scRNAAgent``.

The agent remains subprocess-only (never imports scanpy) and its behavior is
unchanged; it is pinned by ``tests/test_scrna_agent_contract.py``. New code
should import from ``aria.agents.scrna`` directly.
"""
from __future__ import annotations

from aria.agents.scrna import SCRNA_SYSTEM, scRNAAgent

__all__ = ["scRNAAgent", "SCRNA_SYSTEM"]
