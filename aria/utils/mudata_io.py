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

from collections import Counter
from pathlib import Path
import re
from typing import Any


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
    names = [str(n) for n in mod_names]
    for name in names:
        if name.lower() in keys:
            return name
    for name in names:
        low = name.lower()
        if any(k in low for k in keys):
            return name
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


def _inspect_count_matrix(matrix: Any) -> dict[str, Any]:
    """Inspect a count matrix without densifying CSR/CSC sparse storage."""
    import numpy as np
    try:
        from scipy import sparse as sp
        is_sparse = sp.issparse(matrix)
    except ImportError:
        is_sparse = hasattr(matrix, "data") and hasattr(matrix, "nnz")

    if matrix is None:
        return {
            "present": False,
            "is_raw_counts": False,
            "is_finite": False,
            "has_negative": False,
            "is_integer_valued": False,
            "is_binary": False,
            "n_nonzero": 0,
            "max_value": None,
            "reason": "matrix_is_none",
        }

    if is_sparse:
        vals = np.asarray(matrix.data)
        nnz = int(vals.size)
    else:
        arr = np.asarray(matrix)
        vals = arr.ravel()
        nnz = int(np.count_nonzero(vals))

    if vals.size == 0:
        return {
            "present": True,
            "is_raw_counts": True,
            "is_finite": True,
            "has_negative": False,
            "is_integer_valued": True,
            "is_binary": False,
            "n_nonzero": 0,
            "max_value": 0.0,
            "reason": "empty_or_all_zero",
        }

    if np.issubdtype(vals.dtype, np.bool_):
        max_val = float(np.max(vals.astype(np.int8)))
        return {
            "present": True,
            "is_raw_counts": True,
            "is_finite": True,
            "has_negative": False,
            "is_integer_valued": True,
            "is_binary": bool(nnz > 0),
            "n_nonzero": nnz,
            "max_value": max_val,
            "reason": None,
        }

    if np.issubdtype(vals.dtype, np.integer) or np.issubdtype(vals.dtype, np.unsignedinteger):
        has_neg = bool(np.any(vals < 0))
        max_val = float(np.max(vals))
        return {
            "present": True,
            "is_raw_counts": not has_neg,
            "is_finite": True,
            "has_negative": has_neg,
            "is_integer_valued": True,
            "is_binary": bool(not has_neg and nnz > 0 and max_val <= 1.0),
            "n_nonzero": nnz,
            "max_value": max_val,
            "reason": "negative_values" if has_neg else None,
        }

    vals_f = vals.astype(np.float64, copy=False)
    is_finite = bool(np.all(np.isfinite(vals_f)))
    if not is_finite:
        return {
            "present": True,
            "is_raw_counts": False,
            "is_finite": False,
            "has_negative": False,
            "is_integer_valued": False,
            "is_binary": False,
            "n_nonzero": nnz,
            "max_value": None,
            "reason": "non_finite_values",
        }

    has_neg = bool(np.any(vals_f < 0))
    is_int = bool(np.all(np.floor(vals_f) == vals_f))
    max_val = float(np.max(vals_f))
    is_raw = bool(not has_neg and is_int)
    reason = None
    if has_neg:
        reason = "negative_values"
    elif not is_int:
        reason = "fractional_values"
    return {
        "present": True,
        "is_raw_counts": is_raw,
        "is_finite": True,
        "has_negative": has_neg,
        "is_integer_valued": is_int,
        "is_binary": bool(is_raw and nnz > 0 and max_val <= 1.0),
        "n_nonzero": nnz,
        "max_value": max_val,
        "reason": reason,
    }


def _select_raw_counts_matrix(
    adata: Any,
    *,
    modality_label: str,
    require_raw_counts: bool = True,
) -> dict[str, Any]:
    """Select verified raw integer counts from `X`, `layers`, or `raw.X`."""
    candidates: list[tuple[str, Any]] = [("X", adata.X)]
    layers = getattr(adata, "layers", None)
    if layers is not None:
        for layer_key in ("counts", "raw_counts"):
            if layer_key in layers:
                candidates.append((f"layers['{layer_key}']", layers[layer_key]))
    raw_obj = getattr(adata, "raw", None)
    if raw_obj is not None and getattr(raw_obj, "shape", None) == adata.shape:
        candidates.append(("raw.X", raw_obj.X))

    inspections: dict[str, dict[str, Any]] = {}
    for source, mat in candidates:
        info = _inspect_count_matrix(mat)
        inspections[source] = info
        if info["present"] and info["is_raw_counts"]:
            return {
                "ok": True,
                "source": source,
                "matrix": mat,
                "inspection": info,
                "inspections": inspections,
            }

    x_info = inspections.get("X", {})
    if not require_raw_counts and x_info.get("present") and x_info.get("is_finite"):
        return {
            "ok": True,
            "source": "X",
            "matrix": adata.X,
            "inspection": x_info,
            "inspections": inspections,
        }

    reasons = ", ".join(
        f"{src}={info.get('reason') or 'invalid'}"
        for src, info in inspections.items()
    )
    return {
        "ok": False,
        "error_type": "NonIntegerCounts",
        "details": (
            f"Modality '{modality_label}' does not expose non-negative integer "
            f"raw counts ({reasons}). Supply raw integer counts in `.X` or "
            f"`.layers['counts']`."
        ),
        "inspections": inspections,
    }


def read_h5mu_paired(
    path: str,
    *,
    require_raw_counts: bool = True,
) -> dict[str, Any]:
    """Load and validate a paired RNA+ATAC `.h5mu` for V48 Slice I1.

    Enforces:
      1. Both RNA and ATAC modalities exist (`NoRnaModality`, `NoAtacModality`).
      2. Neither modality contains duplicate `obs_names` (`DuplicateBarcodes`).
      3. Shared barcodes are computed in deterministic `rna.obs_names` order,
         never via unordered `set` iteration (`NoSharedBarcodes` when empty).
      4. Raw non-negative integer count matrices are verified without densifying
         sparse storage (`NonIntegerCounts` when `require_raw_counts=True`).
    """
    bad = _validate_path(path)
    if bad:
        return bad

    md = _import_mudata()
    if md is None:
        return {
            "status": "error",
            "error_type": "MuDataToolingMissing",
            "details": (
                "reading .h5mu requires the 'mudata' (or 'muon') package from "
                "the chromatin stack (aria-chromatin-env); it is not importable "
                "in the current environment."
            ),
        }

    try:
        mdata = md.read_h5mu(str(path))
    except Exception as exc:
        return {
            "status": "error",
            "error_type": "H5muReadFailed",
            "details": f"{type(exc).__name__}: {str(exc)[:200]}",
        }

    mod_keys = list(mdata.mod.keys())
    rna_name = _match_modality(mod_keys, _RNA_KEYS)
    if rna_name is None:
        return {
            "status": "error",
            "error_type": "NoRnaModality",
            "details": (
                f"no RNA-like modality among {mod_keys}; expected a name "
                f"matching one of {_RNA_KEYS}."
            ),
        }

    atac_name = _match_modality(
        [k for k in mod_keys if k != rna_name],
        _ATAC_KEYS,
    )
    if atac_name is None:
        return {
            "status": "error",
            "error_type": "NoAtacModality",
            "details": (
                f"no ATAC-like modality among {mod_keys}; expected a name "
                f"matching one of {_ATAC_KEYS}."
            ),
        }

    rna_adata = mdata.mod[rna_name]
    atac_adata = mdata.mod[atac_name]

    rna_obs = [str(bc) for bc in rna_adata.obs_names]
    atac_obs = [str(bc) for bc in atac_adata.obs_names]

    rna_dups = len(rna_obs) - len(set(rna_obs))
    atac_dups = len(atac_obs) - len(set(atac_obs))
    if rna_dups > 0 or atac_dups > 0:
        return {
            "status": "error",
            "error_type": "DuplicateBarcodes",
            "details": (
                f"Paired .h5mu contains duplicate cell barcodes "
                f"(rna_duplicates={rna_dups}, atac_duplicates={atac_dups}). "
                f"Barcodes must be unique per modality."
            ),
        }

    atac_set = set(atac_obs)
    rna_set = set(rna_obs)
    # Deterministic canonical order: preserve RNA obs_names sequence.
    shared_barcodes = [bc for bc in rna_obs if bc in atac_set]
    if not shared_barcodes:
        return {
            "status": "error",
            "error_type": "NoSharedBarcodes",
            "details": (
                f"RNA modality '{rna_name}' (n={len(rna_obs)}, sample="
                f"{rna_obs[:3]}) and ATAC modality '{atac_name}' "
                f"(n={len(atac_obs)}, sample={atac_obs[:3]}) share zero cell "
                f"barcodes. ARIA refuses silent barcode normalization across "
                f"modalities."
            ),
        }

    rna_mat_info = _select_raw_counts_matrix(
        rna_adata,
        modality_label=rna_name,
        require_raw_counts=require_raw_counts,
    )
    if not rna_mat_info["ok"]:
        return {
            "status": "error",
            "error_type": rna_mat_info["error_type"],
            "details": rna_mat_info["details"],
        }

    atac_mat_info = _select_raw_counts_matrix(
        atac_adata,
        modality_label=atac_name,
        require_raw_counts=require_raw_counts,
    )
    if not atac_mat_info["ok"]:
        return {
            "status": "error",
            "error_type": atac_mat_info["error_type"],
            "details": atac_mat_info["details"],
        }

    def _find_col(obs_df, hints):
        for col in obs_df.columns:
            low = str(col).lower()
            if any(h in low for h in hints):
                return str(col)
        return None

    frag_col = _find_col(atac_adata.obs, _FRAG_COL_HINTS)
    mito_col = _find_col(atac_adata.obs, _MITO_COL_HINTS)
    median_frags = None
    median_mito = None
    try:
        import numpy as np
        if frag_col is not None:
            vals = np.asarray(atac_adata.obs[frag_col].values, dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size:
                median_frags = float(np.median(vals))
        if mito_col is not None:
            vals = np.asarray(atac_adata.obs[mito_col].values, dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size:
                median_mito = float(np.median(vals))
    except Exception:
        pass

    n_rna = len(rna_obs)
    n_atac = len(atac_obs)
    n_shared = len(shared_barcodes)

    return {
        "status": "success",
        "path": str(path),
        "mdata": mdata,
        "rna_modality": rna_name,
        "atac_modality": atac_name,
        "rna_adata": rna_adata,
        "atac_adata": atac_adata,
        "rna_counts_source": rna_mat_info["source"],
        "atac_counts_source": atac_mat_info["source"],
        "rna_counts_matrix": rna_mat_info["matrix"],
        "atac_counts_matrix": atac_mat_info["matrix"],
        "rna_is_raw_counts": bool(rna_mat_info["inspection"]["is_raw_counts"]),
        "atac_is_raw_counts": bool(atac_mat_info["inspection"]["is_raw_counts"]),
        "atac_is_binary": bool(atac_mat_info["inspection"]["is_binary"]),
        "n_rna_cells": n_rna,
        "n_atac_cells": n_atac,
        "n_shared_cells": n_shared,
        "n_rna_only": int(len(rna_set - atac_set)),
        "n_atac_only": int(len(atac_set - rna_set)),
        "shared_fraction_rna": float(n_shared / n_rna) if n_rna else 0.0,
        "shared_fraction_atac": float(n_shared / n_atac) if n_atac else 0.0,
        "exact_same_order": bool(rna_obs == atac_obs),
        "shared_barcodes": shared_barcodes,
        "atac_fragment_count_col": frag_col,
        "atac_mito_col": mito_col,
        "median_atac_fragments_per_cell": median_frags,
        "median_atac_mito_fraction": median_mito,
    }


_GEM_SUFFIX_RE = re.compile(r"-\d+$")


def _decode_h5_scalar(val: Any) -> str | None:
    if val is None:
        return None
    if isinstance(val, (bytes, bytearray)):
        s = val.decode("utf-8", errors="replace")
        return s if s != "" else None
    try:
        import numpy as np
        if isinstance(val, np.bytes_):
            s = val.decode("utf-8", errors="replace")
            return s if s != "" else None
        if isinstance(val, (float, np.floating)) and np.isnan(val):
            return None
    except ImportError:
        pass
    s = str(val)
    return s if s != "" else None


def _decode_h5_str_list(raw_arr: Any) -> list[str | None]:
    return [_decode_h5_scalar(v) for v in raw_arr]


def _read_h5_obs_series(obs_obj: Any, col: str) -> tuple[bool, list[str | None] | str]:
    """Read a single `/obs/<col>` column from an `.h5ad` via `h5py` only."""
    import h5py
    import numpy as np

    if isinstance(obs_obj, h5py.Dataset):
        # Legacy AnnData (<0.7) compound structured array for `/obs`.
        names = obs_obj.dtype.names or ()
        if col not in names:
            return False, f"Column '{col}' not found in compound /obs ({list(names)[:30]})."
        raw_vals = obs_obj[col][:]
        return True, _decode_h5_str_list(raw_vals)

    if col not in obs_obj:
        avail = [k for k in obs_obj.keys() if k != "__categories"]
        return False, f"Column '{col}' not found in /obs (available: {avail[:30]})."

    item = obs_obj[col]
    if isinstance(item, h5py.Group):
        if "categories" in item and "codes" in item:
            cats = _decode_h5_str_list(item["categories"][:])
            codes = np.asarray(item["codes"][:], dtype=np.int64)
            n_cats = len(cats)
            out: list[str | None] = [
                cats[int(c)] if 0 <= int(c) < n_cats else None
                for c in codes
            ]
            return True, out
        if "values" in item and "mask" in item:
            vals = _decode_h5_str_list(item["values"][:])
            mask = np.asarray(item["mask"][:], dtype=bool)
            out = [
                None if bool(m) else v
                for v, m in zip(vals, mask)
            ]
            return True, out
        return False, f"Unsupported HDF5 group encoding for /obs/{col}: keys={list(item.keys())}."

    if isinstance(item, h5py.Dataset):
        arr = item[:]
        # Legacy AnnData categorical stored via `/obs/__categories/<col>`.
        if (
            np.issubdtype(arr.dtype, np.integer)
            and "__categories" in obs_obj
            and col in obs_obj["__categories"]
        ):
            cats = _decode_h5_str_list(obs_obj["__categories"][col][:])
            codes = np.asarray(arr, dtype=np.int64)
            n_cats = len(cats)
            out = [
                cats[int(c)] if 0 <= int(c) < n_cats else None
                for c in codes
            ]
            return True, out
        return True, _decode_h5_str_list(arr)

    return False, f"Unsupported HDF5 object type for /obs/{col}: {type(item).__name__}."


def read_h5ad_obs_columns(
    h5ad_path: str,
    columns: list[str] | tuple[str, ...],
) -> dict[str, Any]:
    """Read `obs_names` and selected `/obs` columns from `.h5ad` without loading `X`."""
    p = Path(h5ad_path)
    if not p.exists():
        return {
            "status": "error",
            "error_type": "ReferenceNotFound",
            "details": f"Reference .h5ad file does not exist: {h5ad_path}",
        }
    if p.suffix.lower() != ".h5ad":
        return {
            "status": "error",
            "error_type": "InvalidReferenceFormat",
            "details": f"Expected a .h5ad reference file, got: {p.name}",
        }

    try:
        import h5py
    except ImportError as exc:
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": f"Reading .h5ad /obs metadata requires h5py: {exc}",
        }

    try:
        with h5py.File(str(p), "r") as f:
            if "obs" not in f:
                return {
                    "status": "error",
                    "error_type": "InvalidReferenceH5ad",
                    "details": f"Reference .h5ad has no '/obs' group: {h5ad_path}",
                }
            obs_obj = f["obs"]
            idx_attr = obs_obj.attrs.get("_index", "_index")
            if isinstance(idx_attr, (bytes, bytearray)):
                idx_key = idx_attr.decode("utf-8", errors="replace")
            else:
                idx_key = str(idx_attr)

            if isinstance(obs_obj, h5py.Dataset):
                names = obs_obj.dtype.names or ()
                if idx_key not in names:
                    idx_key = "index" if "index" in names else (names[0] if names else "")
                if not idx_key:
                    return {
                        "status": "error",
                        "error_type": "InvalidReferenceH5ad",
                        "details": "Compound /obs dataset has no index field.",
                    }
                raw_index = obs_obj[idx_key][:]
            else:
                if idx_key not in obs_obj:
                    for fallback in ("_index", "index"):
                        if fallback in obs_obj:
                            idx_key = fallback
                            break
                if idx_key not in obs_obj:
                    return {
                        "status": "error",
                        "error_type": "InvalidReferenceH5ad",
                        "details": (
                            f"Could not locate obs index '{idx_key}' in "
                            f"reference .h5ad: {h5ad_path}"
                        ),
                    }
                raw_index = obs_obj[idx_key][:]

            obs_names = [str(_decode_h5_scalar(v) or "") for v in raw_index]
            col_data: dict[str, list[str | None]] = {}
            for col in columns:
                ok, series_or_err = _read_h5_obs_series(obs_obj, col)
                if not ok:
                    return {
                        "status": "error",
                        "error_type": "ReferenceColumnMissing",
                        "details": str(series_or_err),
                    }
                assert isinstance(series_or_err, list)
                col_data[col] = series_or_err
    except Exception as exc:
        return {
            "status": "error",
            "error_type": "ReferenceReadFailed",
            "details": f"{type(exc).__name__}: {str(exc)[:200]}",
        }

    return {
        "status": "success",
        "obs_names": obs_names,
        "columns": col_data,
        "n_obs": len(obs_names),
    }


def join_obs_labels_from_h5ad(
    query_barcodes: list[str],
    h5ad_path: str,
    *,
    label_col: str = "subclass",
    sample_id: str | None = None,
    reference_sample_col: str | None = None,
    unassigned_label: str = "unassigned",
) -> dict[str, Any]:
    """Deterministic, bijective `/obs` barcode metadata join from a reference `.h5ad`.

    Mechanism: `obs_barcode_join` (exact or explicit `sample_id`-scoped 1-to-1
    barcode correspondence; never an algorithmic label projection, never guesses
    across ambiguous 1-to-N barcode collisions).
    """
    cols_to_read = [label_col]
    if reference_sample_col and reference_sample_col not in cols_to_read:
        cols_to_read.append(reference_sample_col)

    loaded = read_h5ad_obs_columns(h5ad_path, cols_to_read)
    if loaded.get("status") != "success":
        return loaded

    ref_barcodes: list[str] = loaded["obs_names"]
    ref_labels: list[str | None] = loaded["columns"][label_col]

    sid = str(sample_id).strip() if sample_id is not None else ""
    if reference_sample_col:
        if not sid:
            return {
                "status": "error",
                "error_type": "InvalidSampleScope",
                "details": (
                    f"reference_sample_col='{reference_sample_col}' requires a "
                    f"non-empty sample_id parameter."
                ),
            }
        ref_samples = loaded["columns"][reference_sample_col]
        scoped_pairs = [
            (bc, lbl)
            for bc, lbl, smp in zip(ref_barcodes, ref_labels, ref_samples)
            if smp == sid
        ]
        if not scoped_pairs:
            # Allow case-insensitive sample_id match ONLY when it maps to a
            # single unique sample level in reference_sample_col.
            sid_low = sid.lower()
            ci_levels = sorted({
                str(smp) for smp in ref_samples
                if smp is not None and str(smp).lower() == sid_low
            })
            if len(ci_levels) == 1:
                matched_level = ci_levels[0]
                scoped_pairs = [
                    (bc, lbl)
                    for bc, lbl, smp in zip(ref_barcodes, ref_labels, ref_samples)
                    if smp == matched_level
                ]
            elif len(ci_levels) > 1:
                return {
                    "status": "error",
                    "error_type": "AmbiguousSampleScope",
                    "details": (
                        f"sample_id={sid!r} matched multiple case-distinct "
                        f"levels in obs['{reference_sample_col}']: {ci_levels}."
                    ),
                }
        if not scoped_pairs:
            return {
                "status": "error",
                "error_type": "SampleIdNotFoundInReference",
                "details": (
                    f"No rows in reference .h5ad had "
                    f"obs['{reference_sample_col}'] == {sid!r}."
                ),
            }
        ref_barcodes = [p[0] for p in scoped_pairs]
        ref_labels = [p[1] for p in scoped_pairs]

    q_list = [str(b) for b in query_barcodes]
    if len(q_list) != len(set(q_list)):
        return {
            "status": "error",
            "error_type": "DuplicateQueryBarcodes",
            "details": "Query barcodes contain duplicates; 1-to-1 join refused.",
        }

    def _strip_gem(bc: str) -> str:
        return _GEM_SUFFIX_RE.sub("", bc)

    def _extract_after_sample_prefix(bc: str, sample_tok: str, *, ci: bool = False) -> str | None:
        n = len(sample_tok)
        if len(bc) <= n + 1:
            return None
        sep = bc[n]
        if sep not in ("_", "#"):
            return None
        head = bc[:n]
        if (head.lower() == sample_tok.lower()) if ci else (head == sample_tok):
            return bc[n + 1:]
        return None

    strategies: list[tuple[str, Any, Any]] = [
        ("exact", lambda q: q, lambda r: r),
    ]
    if sid:
        strategies.extend([
            ("sample_prefix_underscore", lambda q, s=sid: f"{s}_{q}", lambda r: r),
            ("query_has_sample_prefix_underscore", lambda q: q, lambda r, s=sid: f"{s}_{r}"),
            ("sample_prefix_hash", lambda q, s=sid: f"{s}#{q}", lambda r: r),
            ("query_has_sample_prefix_hash", lambda q: q, lambda r, s=sid: f"{s}#{r}"),
            ("sample_suffix_dash", lambda q, s=sid: f"{q}-{s}", lambda r: r),
            ("sample_suffix_underscore", lambda q, s=sid: f"{q}_{s}", lambda r: r),
            (
                "sample_prefix_ci",
                lambda q: q,
                lambda r, s=sid: _extract_after_sample_prefix(r, s, ci=True),
            ),
            (
                "sample_prefix_stripped_gem",
                lambda q: _strip_gem(q),
                lambda r, s=sid: (
                    _strip_gem(rest)
                    if (rest := _extract_after_sample_prefix(r, s, ci=True)) is not None
                    else None
                ),
            ),
        ])
        if reference_sample_col:
            _trail_re = re.compile(r"(?:^|[_#])([ACGTNacgtn]{6,}(?:-\d+)?)$")

            def _trailing_bc(bc: str) -> str | None:
                m = _trail_re.search(bc)
                return m.group(1) if m else None

            strategies.extend([
                (
                    "trailing_barcode_within_sample_scope",
                    lambda q: _trailing_bc(q) or q,
                    _trailing_bc,
                ),
                (
                    "stripped_gem_within_sample_scope",
                    lambda q: _strip_gem(_trailing_bc(q) or q),
                    lambda r: (
                        _strip_gem(tb)
                        if (tb := _trailing_bc(r)) is not None
                        else _strip_gem(r)
                    ),
                ),
            ])

    chosen_strategy: str | None = None
    chosen_ref_map: dict[str, str | None] = {}
    chosen_q_keys: list[str] = []

    for strat_name, q_fn, r_fn in strategies:
        q_keys = [q_fn(q) for q in q_list]
        if len(q_keys) != len(set(q_keys)):
            # Query barcode transform collapsed distinct cells (e.g. AAAC-1 & AAAC-2).
            q_key_set = set(q_keys)
            r_keys_probe = {rk for rk in (r_fn(r) for r in ref_barcodes) if rk is not None}
            if q_key_set & r_keys_probe:
                return {
                    "status": "error",
                    "error_type": "AmbiguousBarcodeJoin",
                    "details": (
                        f"Barcode strategy '{strat_name}' introduces 1-to-N "
                        f"collisions across query cells (e.g. multiple GEM "
                        f"wells sharing a core sequence). Join refused."
                    ),
                }
            continue

        ref_buckets: dict[str, list[str | None]] = {}
        for r_bc, r_lbl in zip(ref_barcodes, ref_labels):
            rk = r_fn(r_bc)
            if rk is None:
                continue
            ref_buckets.setdefault(rk, []).append(r_lbl)

        q_key_set = set(q_keys)
        matched_keys = [k for k in q_keys if k in ref_buckets]
        if not matched_keys:
            continue

        # Enforce global 1-to-1 uniqueness across candidate reference keys in
        # the matched namespace (any duplicate key in reference fails closed).
        colliding_ref = [
            rk for rk, vals in ref_buckets.items()
            if len(vals) > 1 and (rk in q_key_set or strat_name != "exact")
        ]
        if colliding_ref:
            return {
                "status": "error",
                "error_type": "AmbiguousBarcodeJoin",
                "details": (
                    f"Barcode strategy '{strat_name}' encountered duplicate "
                    f"reference keys (n_collisions={len(colliding_ref)}, "
                    f"sample={colliding_ref[:3]}). ARIA refuses ambiguous "
                    f"1-to-N metadata joins."
                ),
            }

        chosen_strategy = strat_name
        chosen_ref_map = {rk: vals[0] for rk, vals in ref_buckets.items() if len(vals) == 1}
        chosen_q_keys = q_keys
        break

    if chosen_strategy is None:
        return {
            "status": "success",
            "mechanism": "obs_barcode_join",
            "match_strategy": None,
            "n_query": len(q_list),
            "n_reference_in_scope": len(ref_barcodes),
            "n_matched": 0,
            "n_unmatched": len(q_list),
            "match_fraction": 0.0,
            "labels": [unassigned_label] * len(q_list),
            "matched_mask": [False] * len(q_list),
            "label_counts": {unassigned_label: len(q_list)} if q_list else {},
            "sample_query_barcodes": q_list[:3],
            "sample_reference_barcodes": ref_barcodes[:3],
        }

    labels_out: list[str] = []
    matched_mask: list[bool] = []
    n_matched = 0
    for qk in chosen_q_keys:
        if qk in chosen_ref_map:
            raw_lbl = chosen_ref_map[qk]
            if raw_lbl is not None and str(raw_lbl).strip() != "":
                labels_out.append(str(raw_lbl))
                matched_mask.append(True)
                n_matched += 1
            else:
                # Matched barcode had NA/empty label in reference -> unassigned.
                labels_out.append(unassigned_label)
                matched_mask.append(False)
        else:
            labels_out.append(unassigned_label)
            matched_mask.append(False)

    n_query = len(q_list)
    counts = dict(sorted(Counter(labels_out).items(), key=lambda kv: (-kv[1], kv[0])))

    return {
        "status": "success",
        "mechanism": "obs_barcode_join",
        "match_strategy": chosen_strategy,
        "n_query": n_query,
        "n_reference_in_scope": len(ref_barcodes),
        "n_matched": n_matched,
        "n_unmatched": n_query - n_matched,
        "match_fraction": float(n_matched / n_query) if n_query else 0.0,
        "labels": labels_out,
        "matched_mask": matched_mask,
        "label_counts": counts,
    }
