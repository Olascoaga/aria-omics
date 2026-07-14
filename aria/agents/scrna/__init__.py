"""scRNAAgent subpackage (A7 extraction of the scrna_agent.py monolith).

The single 2.3k-line ``scRNAAgent`` class is split into concern mixins composed
by ``agent.scRNAAgent``. The public surface (the class + SCRNA_SYSTEM prompt) is
re-exported here and again from the ``aria/agents/scrna_agent.py`` facade so
existing consumers keep importing the same path. Behavior pinned by
``tests/test_scrna_agent_contract.py``.
"""
from __future__ import annotations

from aria.agents.scrna._base import SCRNA_SYSTEM
from aria.agents.scrna.agent import scRNAAgent

__all__ = ["scRNAAgent", "SCRNA_SYSTEM"]
