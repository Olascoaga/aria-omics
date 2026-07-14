"""
ARIA NarrativeAgent — scRNA / pseudobulk extension (compatibility facade)
-------------------------------------------------------------------------
The implementation moved to the cohesive ``aria.agents.narrative.scrna``
subpackage in A7 (formatting/text/tables/figures split out of a 2.3k-line
monolith). This module is kept as a stable import path: existing consumers
(``narrative_agent``, ``report_builder``, ``devils_advocate``, the scRNA
narrator, tests) still do ``from aria.agents import _narrative_scrna`` and read
``_narrative_scrna.<name>``.

The functions remain pure (no LLM, no bus) and unchanged; behavior is pinned by
``tests/test_narrative_scrna_contract.py``. New code should import from
``aria.agents.narrative.scrna`` directly.
"""
from __future__ import annotations

from aria.agents.narrative.scrna import *  # noqa: F401,F403
from aria.agents.narrative.scrna import __all__  # noqa: F401
