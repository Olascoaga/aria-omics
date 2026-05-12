"""
ARIA scRNA — inject condition obs column from DesignAgent group mapping.

Used by scrna_agent._inject_condition_obs to bridge between the
DesignAgent output (sample_id → group_name mapping) and the obs schema
that rna_pseudobulk_de expects (one obs column per cell with the
condition label).

The cell-level "sample" column is auto-discovered in this priority:
    obs.sample_id  →  obs.batch  →  obs['orig.ident']  →  guess from
    obs_name prefix when names are formatted as "<sample>_<barcode>".

Subprocess interface (aria.scripts._base.run_script):
    python rna_inject_condition.py <in.json> <out.json>

Input params:
    data_path:    str  — h5ad to read
    groups:       dict — {group_name: [sample_ids]} from DesignAgent
    factor:       str  — name of the obs column to write (e.g. 'age_group')
    batch_col:    str  (optional) — propagated for downstream reference
    output_path:  str  — where to write the modified h5ad

Output:
    {status, output_path, condition_col, replicate_col, matched_cells,
     unmatched_cells, samples_seen, reason}
"""

from __future__ import annotations

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from aria.scripts._base import run_script


def _build_sample_to_group(groups: dict) -> dict:
    """Invert {group: [samples]} → {sample: group}."""
    mapping: dict = {}
    for grp, samples in (groups or {}).items():
        for s in samples or []:
            mapping[str(s)] = str(grp)
    return mapping


def _pick_sample_col(adata) -> tuple:
    """
    Return (sample_col, sample_values_array). The sample column is the
    one whose values are donor / replicate IDs that we can map to
    DesignAgent groups. Tries known columns; falls back to parsing the
    obs_names prefix.
    """
    candidates = ("sample_id", "batch", "orig.ident", "sample", "donor")
    for c in candidates:
        if c in adata.obs.columns:
            return c, adata.obs[c].astype(str).values
    # Fallback: <sample>_<barcode> prefix
    import numpy as np
    prefixes = [str(n).rsplit("_", 1)[0] for n in adata.obs_names]
    return "_obs_name_prefix", np.asarray(prefixes)


def inject(params: dict) -> dict:
    import anndata as ad
    import numpy as np
    import pandas as pd
    from pathlib import Path

    data_path   = params["data_path"]
    groups      = params.get("groups", {}) or {}
    factor      = params.get("factor", "condition")
    batch_col   = params.get("batch_col")
    output_path = params["output_path"]

    if not groups:
        return {"status": "skipped", "reason": "no_groups_in_design"}

    s2g = _build_sample_to_group(groups)
    if not s2g:
        return {"status": "skipped", "reason": "groups_empty"}

    adata = ad.read_h5ad(data_path)
    sample_col, sample_values = _pick_sample_col(adata)

    # Map each cell's sample → group. Cells from samples not in the
    # design mapping become NaN (we drop them before writing so
    # pseudobulk doesn't see ambiguous labels).
    mapped = pd.Series(sample_values).map(s2g)
    matched = int(mapped.notna().sum())
    unmatched = int(mapped.isna().sum())
    if matched == 0:
        return {
            "status": "skipped",
            "reason": (
                f"no cells matched any design sample. "
                f"sample_col='{sample_col}', "
                f"samples seen: {sorted(set(sample_values))[:5]}…, "
                f"design samples: {sorted(s2g.keys())[:5]}…"
            ),
            "samples_seen":    sorted(set(map(str, sample_values)))[:50],
            "design_samples":  sorted(s2g.keys()),
        }

    # Keep only matched cells (drop unmatched silently — they'd just
    # contribute to NaN factor levels in DESeq2 design).
    mask = mapped.notna().values
    if unmatched:
        adata = adata[mask].copy()
        mapped = mapped[mask].reset_index(drop=True)

    adata.obs[factor] = pd.Categorical(mapped.astype(str))

    # Surface the chosen replicate column so the caller can pass it to
    # rna_pseudobulk_de without re-deriving.
    replicate_col = (
        "sample_id" if "sample_id" in adata.obs.columns
        else sample_col
    )
    # If we synthesised the prefix, also persist it as sample_id so
    # pseudobulk has a stable obs column to group on.
    if sample_col == "_obs_name_prefix":
        adata.obs["sample_id"] = pd.Categorical(
            [str(n).rsplit("_", 1)[0] for n in adata.obs_names]
        )
        replicate_col = "sample_id"

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

    return {
        "status":           "success",
        "output_path":      output_path,
        "condition_col":    factor,
        "replicate_col":    replicate_col,
        "batch_col":        batch_col,
        "matched_cells":    matched,
        "unmatched_cells":  unmatched,
        "n_groups":         len(set(mapped)),
        "samples_seen":     sorted(set(map(str, sample_values)))[:50],
    }


if __name__ == "__main__":
    run_script(inject)
