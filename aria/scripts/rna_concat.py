"""
ARIA RNA Concatenation Script
------------------------------
Concatenates a list of post-QC AnnData objects into a single .h5ad ready
for batch correction (Harmony) and downstream clustering.

Design notes:
  - Inner join on var_names so every cell has values for the same gene
    set (required by Harmony / PCA — outer joins fill with NaN and break
    sparse matmuls).
  - obs["sample_id"] and obs["batch"] are populated from the manifest
    (overwriting whatever rna_qc wrote, so the caller controls the labels).
  - obs_names are prefixed with sample_id to guarantee global uniqueness.
  - adata.raw is preserved when present on every input.

Input params:
    samples:   [ {path: str, sample_id: str, ...arbitrary metadata}, ...]
               Each path must point to an .h5ad (the output of rna_qc.py).
               Metadata keys other than path/sample_id are propagated to
               obs as flat columns (useful for condition / age / etc.).
    output_dir:  str (optional) — where to write concatenated.h5ad
                                   (default: parent of first input)
    join:        "inner" | "outer" (default: "inner")

Output:
    {
      "status":          "success" | "error",
      "n_samples":       int,
      "n_cells_total":   int,
      "n_genes_shared":  int,
      "per_sample":      [{sample_id, n_cells, n_genes}, ...],
      "output_path":     str,
      "batch_col":       "batch",
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def _cache_params(params: dict) -> dict:
    samples = []
    for sample in params.get("samples") or []:
        samples.append({
            k: str(v)
            for k, v in sorted(sample.items())
        })
    return {
        "samples": samples,
        "join": str(params.get("join", "inner")),
    }


def _cache_matches(cached: dict, expected: dict) -> bool:
    return cached.get("cache_version") == 2 and cached.get("cache_params") == expected


def rna_concat(params: dict) -> dict:
    import anndata as ad
    import scanpy as sc
    import json
    from pathlib import Path

    samples    = params.get("samples") or []
    output_dir = params.get("output_dir")
    join       = params.get("join", "inner")
    cache_params = _cache_params(params)

    if len(samples) < 2:
        return {
            "status":     "error",
            "error_type": "TooFewSamples",
            "details":    (f"rna_concat needs at least 2 samples; "
                           f"received {len(samples)}."),
        }

    out_dir = Path(output_dir) if output_dir else Path(samples[0]["path"]).parent
    out_path = out_dir / "concatenated.h5ad"
    summary_path = out_path.with_suffix(".summary.json")
    input_paths = [Path(s.get("path", "")) for s in samples if s.get("path")]
    if (out_path.exists() and input_paths and
            all(p.exists() for p in input_paths) and
            out_path.stat().st_mtime >= max(p.stat().st_mtime for p in input_paths)):
        try:
            if summary_path.exists():
                with open(summary_path) as f:
                    cached = json.load(f)
                if not _cache_matches(cached, cache_params):
                    raise ValueError("stale concat cache")
                cached["resumed"] = True
                cached["warnings"] = (
                    cached.get("warnings", []) +
                    ["[resume] rna_concat skipped (valid concatenated.h5ad exists)"]
                )
                return cached
            # Legacy outputs without signed summaries are rerun to avoid stale
            # sample metadata or join settings leaking into downstream design.
            pass
        except Exception:
            pass

    # ── Load each input and stamp provenance ─────────────────────────────
    adatas      = []
    per_sample  = []
    extra_keys: set = set()
    for s in samples:
        path = s.get("path")
        sid  = s.get("sample_id")
        if not path or not sid:
            return {
                "status":     "error",
                "error_type": "InvalidManifest",
                "details":    (f"Each sample entry needs both 'path' and "
                               f"'sample_id'. Got: {s}"),
            }
        if not Path(path).exists():
            return {
                "status":     "error",
                "error_type": "FileNotFound",
                "details":    f"Missing input: {path}",
            }
        a = sc.read_h5ad(path)
        a.var_names_make_unique()
        # Overwrite (or set) provenance — caller is authoritative.
        a.obs["sample_id"] = str(sid)
        a.obs["batch"]     = str(sid)
        # Propagate any other metadata keys as flat obs columns. This is
        # how condition / age / sex labels reach Harmony / DE downstream.
        for k, v in s.items():
            if k in ("path", "sample_id"):
                continue
            a.obs[k] = v
            extra_keys.add(k)
        # Prefix obs_names to keep cell barcodes globally unique after concat.
        a.obs_names = [f"{sid}_{b}" for b in a.obs_names]
        adatas.append(a)
        per_sample.append({
            "sample_id": str(sid),
            "n_cells":   int(a.n_obs),
            "n_genes":   int(a.n_vars),
        })

    # ── Concatenate ──────────────────────────────────────────────────────
    # outer_join would fill missing genes with 0 / NaN; inner is the safe
    # default for downstream Harmony + PCA which expect a complete matrix.
    combined = ad.concat(
        adatas,
        axis=0, join=join,
        merge="same",          # keep var entries only when identical across inputs
        label=None,            # we already set obs["batch"] / obs["sample_id"]
        index_unique=None,     # we prefixed obs_names manually above
    )

    # Re-establish .raw from one of the inputs so DE / pathway scripts can
    # still address non-HVG genes (rna_integration also preserves raw, so
    # there's a chain of custody from QC through integration to DE).
    if adatas[0].raw is not None:
        # adata.raw needs the same var/obs alignment — easiest is to set
        # raw to a copy of the combined object before any HVG subsetting.
        combined.raw = combined.copy()

    # ── Write ────────────────────────────────────────────────────────────
    out_dir.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(str(out_path))

    result = {
        "status":          "success",
        "n_samples":       len(samples),
        "n_cells_total":   int(combined.n_obs),
        "n_genes_shared":  int(combined.n_vars),
        "per_sample":      per_sample,
        "extra_obs_keys":  sorted(extra_keys),
        "output_path":     str(out_path),
        "batch_col":       "batch",
        "join":            join,
        "cache_version":   2,
        "cache_params":    cache_params,
    }
    try:
        with open(summary_path, "w") as f:
            json.dump(result, f, indent=2)
    except Exception:
        pass
    return result


if __name__ == "__main__":
    run_script(rna_concat)
