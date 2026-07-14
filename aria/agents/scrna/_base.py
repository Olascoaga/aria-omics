"""Shared base for the scRNAAgent subpackage: imports, logger, system prompt.

Extracted from aria/agents/scrna_agent.py (A7). Each concern mixin imports `*`
from here so every moved method keeps its original module-level references.
"""
from __future__ import annotations

import hashlib
import json
import logging
import re
import uuid
from pathlib import Path

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence
from aria.llm.provider import LLMProvider, TaskTier
from aria.llm.prompt_boundary import (
    PromptDataField,
    build_untrusted_prompt,
    system_with_untrusted_boundary,
)
from aria.llm.parameter_advisor import ParameterAdvisor
from aria.memory.memory import ARIAMemory

log = logging.getLogger("aria.scrna")


SCRNA_SYSTEM = """
You are ARIA's scRNAAgent — a specialist in single-cell RNA-seq analysis.

Your expertise:
- scRNA-seq QC: doublet detection, mitochondrial filtering, MAD thresholds
- Normalization: scran, log1p; feature selection: HVGs
- Dimensionality reduction: PCA, UMAP, tSNE
- Clustering: Leiden algorithm, resolution selection
- Cell type annotation: marker-based, known cell type databases
- Differential expression: Wilcoxon per cluster, pseudo-bulk per condition
- Trajectory analysis: PAGA, DPT pseudotime, RNA velocity (scVelo)
- Cell-cell communication: LIANA, CellChat, NicheNet
- Batch correction: Harmony, scVI

Critical knowledge:
- Marker gene extraction: use names[cluster] NOT names[0]
- Single cells are not replicates: use pseudo-bulk for condition comparisons
- Mitochondrial % cutoffs must be context-aware (stressed cells have high MT%)
- Leiden resolution is data-dependent: always use ParameterAdvisor
- Harmony works on PCA embeddings: recompute neighbors after correction
- PAGA needs a root cell for meaningful pseudotime direction
""".strip()


__all__ = [
    "hashlib", "json", "logging", "re", "uuid", "Path", "BaseAgent",
    "Confidence", "LLMProvider", "TaskTier", "PromptDataField",
    "build_untrusted_prompt", "system_with_untrusted_boundary",
    "ParameterAdvisor", "ARIAMemory", "log", "SCRNA_SYSTEM",
]
