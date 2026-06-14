"""scATAC / chromatin figure generation (W0.1, scATAC P0 pre-preprint plan).

Mirrors ``aria.agents._narrative_scrna.generate_figures``: render publication
figures in the chromatin conda stack via ``EnvironmentManager`` and write their
paths into ``findings["figures"]``, which ``ChromatinNarrator.figures()`` surfaces.

No fabrication: a figure is produced ONLY when its underlying data exists (a
clustered ``.h5ad`` with an embedding) and an env manager is available. With
``env_manager=None`` or no ``.h5ad`` the figure set stays empty (honest absence),
exactly like the modality's QC ``None`` discipline.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from aria.agents.narrative.narrators.chromatin import _CHROMATIN_MODALITY_KEYS

log = logging.getLogger(__name__)


def live_findings(agent_result: dict) -> Optional[dict]:
    """Return the LIVE findings dict that ``unwrap_chromatin_findings`` surfaces.

    ``unwrap_chromatin_findings`` rebuilds a fresh dict when the per-modality
    wrapper is present, so a caller that wants a mutation (adding ``figures``) to be
    visible on the narrator's re-unwrap must mutate the SAME nested dict unwrap
    reads from. This returns that object.
    """
    if not isinstance(agent_result, dict):
        return None
    findings = agent_result.get("findings", agent_result)
    if not isinstance(findings, dict):
        return None
    for key, value in findings.items():
        if key in _CHROMATIN_MODALITY_KEYS and isinstance(value, dict):
            nested = value.get("findings")
            return nested if isinstance(nested, dict) else value
    return findings


def _umap_color_keys(findings: dict) -> list:
    """Obs columns to color the UMAP by, in priority order. ``rna_figure_umap``
    silently skips columns that are absent, so listing fallbacks is safe."""
    da = findings.get("differential_accessibility") or {}
    lsi = findings.get("lsi") or findings.get("lsi_clustering") or {}
    groupby = da.get("groupby") or lsi.get("cluster_key") or "leiden"
    keys: list = []
    for candidate in (groupby, "leiden", "sample", "sample_id", "batch",
                      "condition"):
        if candidate and candidate not in keys:
            keys.append(candidate)
    return keys


def generate_figures(findings: dict, h5ad_path: Optional[str], output_dir,
                     env_manager=None) -> dict:
    """Generate scATAC figures into ``findings['figures']`` (mutates + returns it).

    W0.1 renders the **UMAP** (colored by cluster + sample/condition when present)
    from the clustered ``.h5ad`` via the generic ``rna_figure_umap.py`` run in the
    ``chromatin`` stack (the script is modality-agnostic: it plots an existing
    embedding and skips missing obs columns). Later W0.x items add fragment-size /
    QC violins / marker heatmap / DA volcano-MA / motif dotplot here.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figs = findings.setdefault("figures", {})

    if h5ad_path and env_manager is not None:
        color_keys = _umap_color_keys(findings)
        try:
            res = env_manager.run_in_stack(
                stack="chromatin",
                script_path="aria/scripts/rna_figure_umap.py",
                params={
                    "h5ad_path":  str(h5ad_path),
                    "color_by":   color_keys,
                    "output_dir": str(output_dir),
                },
            )
            if res.get("status") == "success":
                for key, path in (res.get("figures") or {}).items():
                    figs[f"umap_{key}"] = path
            else:
                log.warning(
                    "chromatin UMAP figure failed: %s — %s",
                    res.get("error_type"), str(res.get("details", ""))[:200],
                )
        except Exception as e:  # subprocess crash must not abort the report
            log.warning("chromatin UMAP figure subprocess crashed: %s", e)

    return findings
