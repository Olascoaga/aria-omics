"""S5 (pre-integration audit): Docs Drift Guard.

Codex finding #2: public chromatin docs drifted behind the real product (bulk ATAC
called "scaffolded", scATAC called "alpha"). This fence makes the docs verifiable
against the single source of truth (MODALITY_VALIDATION in the orchestrator):

  1. validation_status.md carries a machine-readable "Authoritative modality tiers"
     table that must EQUAL MODALITY_VALIDATION_LEVELS exactly (code <-> doc drift).
  2. the chromatin roadmap must not reassert the specific obsolete states this slice
     corrected (scATAC "alpha runs"; bulk ATAC as a measured-QC/MACS3-only or
     scaffolded-count-matrix slice).

Parses files only — runs in any env.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_VALIDATION_STATUS = _ROOT / "docs" / "validation_status.md"
_CHROMATIN_ROADMAP = _ROOT / "docs" / "workflows" / "chromatin_roadmap.md"

_TABLE_RE = re.compile(
    r"<!-- MODALITY_TIERS_TABLE_START -->(.*?)<!-- MODALITY_TIERS_TABLE_END -->",
    re.DOTALL,
)


def _doc_tiers() -> dict[str, str]:
    text = _VALIDATION_STATUS.read_text(encoding="utf-8")
    m = _TABLE_RE.search(text)
    assert m, "validation_status.md is missing the MODALITY_TIERS_TABLE markers"
    tiers: dict[str, str] = {}
    for line in m.group(1).splitlines():
        line = line.strip()
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) != 2:
            continue
        modality, tier = cells
        if modality in ("Modality", "---") or tier in ("Tier", "---"):
            continue
        tiers[modality] = tier
    return tiers


def test_doc_tier_table_matches_single_source():
    from aria.agents.orchestrator_agent import MODALITY_VALIDATION_LEVELS

    doc = _doc_tiers()
    src = dict(MODALITY_VALIDATION_LEVELS)
    # The doc table is the public projection of the dispatchable modalities; it must
    # match the single source exactly (no missing, no extra, no wrong tier).
    assert doc == src, (
        "validation_status.md modality tiers drifted from MODALITY_VALIDATION:\n"
        f"  doc only:  { {k: v for k, v in doc.items() if doc.get(k) != src.get(k)} }\n"
        f"  code only: { {k: v for k, v in src.items() if src.get(k) != doc.get(k)} }"
    )


def test_roadmap_does_not_reassert_corrected_obsolete_states():
    text = _CHROMATIN_ROADMAP.read_text(encoding="utf-8").lower()
    # scATAC is beta (ADR-048), not alpha: the old "available for reviewed alpha
    # runs" phrasing must not return.
    assert "alpha runs" not in text, (
        "chromatin_roadmap.md still calls scATAC an 'alpha' run (it is beta)"
    )
    # bulk ATAC is the full workflow, not a measured-QC/MACS3-only or scaffolded
    # count-matrix slice.
    assert "scaffolded peak-by-sample count" not in text, (
        "chromatin_roadmap.md still calls the bulk ATAC count matrix 'scaffolded'"
    )
    assert "beta slice for measured qc + macs3" not in text, (
        "chromatin_roadmap.md still frames bulk ATAC as a QC+MACS3-only slice"
    )


def test_roadmap_points_to_single_source():
    text = _CHROMATIN_ROADMAP.read_text(encoding="utf-8")
    assert "MODALITY_VALIDATION" in text, (
        "chromatin_roadmap.md should point readers to the single-source tier table"
    )
