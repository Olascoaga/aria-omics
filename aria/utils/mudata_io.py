"""
ARIA MuData (.h5mu) I/O helpers
-------------------------------
Reads same-cell paired RNA+ATAC MuData files for the v4.6 chromatin entry path
(C8, audit 2026-05-29). This is a REAL reader, not a placeholder: when the
MuData tooling is absent it returns a structured blocker rather than fabricating
modality metadata (ADR-002). The previous state of the world had no `.h5mu`
detection at all, so the planned `hc11_paired.h5mu` validation input was
undetectable and would have fallen into an unvalidated path.
"""

from __future__ import annotations

from pathlib import Path


# obs-column name fragments that commonly carry per-cell ATAC fragment counts
# and mitochondrial fractions. Matched case-insensitively as substrings; ARIA
# only reports these metrics when a real column is present (no fabrication).
_FRAG_COL_HINTS = (
    "atac_fragments", "n_fragments", "nfrags", "ncount_atac",
    "passed_filters", "fragments",
)
_MITO_COL_HINTS = ("pct_reads_mito", "pct_mito", "mito_frac", "percent_mito")

_ATAC_KEYS = ("atac", "peaks", "chromatin", "accessibility")
_RNA_KEYS = ("rna", "gex", "expression")


def _import_mudata():
    """Return the `mudata` module, or None if neither mudata nor muon is
    importable in the current environment."""
    try:
        import mudata
        return mudata
    except ImportError:
        try:
            import muon  # noqa: F401  (muon pulls mudata as a dependency)
            import mudata
            return mudata
        except ImportError:
            return None


def _validate_path(path: str) -> dict | None:
    p = Path(path)
    if not p.exists():
        return {"status": "error", "error_type": "FileNotFound",
                "details": f"h5mu path does not exist: {path}"}
    if p.suffix.lower() != ".h5mu":
        return {"status": "error", "error_type": "NotAnH5mu",
                "details": f"expected a .h5mu (MuData) file, got: {p.name}"}
    return None


def _match_modality(mod_names, keys) -> str | None:
    for name in mod_names:
        low = str(name).lower()
        if any(k in low for k in keys):
            return str(name)
    return None


def read_h5mu_summary(path: str, backed: bool = True) -> dict:
    """
    Summarize a `.h5mu` without fabricating anything: modality names, per-modality
    (n_obs, n_vars), shared obs columns, and which modality looks like ATAC / RNA.

    Structured errors (no exceptions leak to the caller):
      - FileNotFound         : path missing
      - NotAnH5mu            : wrong suffix
      - MuDataToolingMissing : neither mudata nor muon importable
      - H5muReadFailed       : the reader raised
    """
    bad = _validate_path(path)
    if bad:
        return bad

    md = _import_mudata()
    if md is None:
        return {"status": "error", "error_type": "MuDataToolingMissing",
                "details": ("reading .h5mu requires the 'mudata' (or 'muon') "
                            "package from the v4.6 chromatin stack; it is not "
                            "importable in the current environment.")}

    try:
        try:
            mdata = md.read_h5mu(str(path), backed="r" if backed else None)
        except TypeError:
            mdata = md.read_h5mu(str(path))
    except Exception as exc:
        return {"status": "error", "error_type": "H5muReadFailed",
                "details": f"{type(exc).__name__}: {str(exc)[:200]}"}

    modalities = {}
    for name, mod in mdata.mod.items():
        modalities[str(name)] = {
            "n_obs": int(mod.n_obs),
            "n_vars": int(mod.n_vars),
            "obs_columns": [str(c) for c in mod.obs.columns][:50],
        }

    atac_mod = _match_modality(mdata.mod.keys(), _ATAC_KEYS)
    rna_mod = _match_modality(mdata.mod.keys(), _RNA_KEYS)

    return {
        "status": "success",
        "path": str(path),
        "n_modalities": len(modalities),
        "modalities": modalities,
        "n_obs": int(mdata.n_obs),
        "shared_obs_columns": [str(c) for c in mdata.obs.columns][:50],
        "atac_modality": atac_mod,
        "rna_modality": rna_mod,
        "is_paired_rna_atac": bool(atac_mod and rna_mod),
        "tooling": getattr(md, "__name__", "mudata"),
    }


def read_h5mu_atac(path: str):
    """
    Load the ATAC modality of a `.h5mu` as an AnnData for QC, plus real per-cell
    metrics derived from its obs (fragment counts / mito fraction) WHEN those
    columns exist — otherwise they are reported as not-computed, never invented.

    Returns either:
      {"status": "success", "modality": str, "adata": AnnData,
       "n_obs": int, "n_vars": int,
       "fragment_count_col": str|None, "mito_col": str|None,
       "median_fragments_per_cell": float|None, "median_mito_fraction": float|None}
    or a structured error dict (same error_types as read_h5mu_summary, plus
    NoAtacModality).
    """
    bad = _validate_path(path)
    if bad:
        return bad

    md = _import_mudata()
    if md is None:
        return {"status": "error", "error_type": "MuDataToolingMissing",
                "details": ("reading .h5mu requires the 'mudata' (or 'muon') "
                            "package from the v4.6 chromatin stack; it is not "
                            "importable in the current environment.")}

    try:
        mdata = md.read_h5mu(str(path))
    except Exception as exc:
        return {"status": "error", "error_type": "H5muReadFailed",
                "details": f"{type(exc).__name__}: {str(exc)[:200]}"}

    atac_name = _match_modality(mdata.mod.keys(), _ATAC_KEYS)
    if atac_name is None:
        return {"status": "error", "error_type": "NoAtacModality",
                "details": (f"no ATAC-like modality among {list(mdata.mod.keys())}; "
                            f"expected a name containing one of {_ATAC_KEYS}.")}

    adata = mdata.mod[atac_name]
    obs = adata.obs

    def _find_col(hints):
        for col in obs.columns:
            low = str(col).lower()
            if any(h in low for h in hints):
                return str(col)
        return None

    frag_col = _find_col(_FRAG_COL_HINTS)
    mito_col = _find_col(_MITO_COL_HINTS)

    median_frags = None
    median_mito = None
    try:
        import numpy as np
        if frag_col is not None:
            vals = np.asarray(obs[frag_col].values, dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size:
                median_frags = float(np.median(vals))
        if mito_col is not None:
            vals = np.asarray(obs[mito_col].values, dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size:
                median_mito = float(np.median(vals))
    except Exception:
        pass

    return {
        "status": "success",
        "modality": atac_name,
        "adata": adata,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "fragment_count_col": frag_col,
        "mito_col": mito_col,
        "median_fragments_per_cell": median_frags,
        "median_mito_fraction": median_mito,
    }
