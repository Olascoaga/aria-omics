"""Small helpers for validity-based h5ad resume logic."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def h5ad_is_readable(path: str | Path) -> bool:
    """Return True only when an h5ad can be opened by anndata."""
    try:
        adata = read_h5ad(path, backed="r")
        if hasattr(adata, "file") and adata.file is not None:
            adata.file.close()
        return True
    except Exception:
        return False


def read_h5ad(path: str | Path, **kwargs: Any):
    """Read an h5ad with a clearer error for truncated or unreadable files."""
    try:
        import anndata as ad
        return ad.read_h5ad(path, **kwargs)
    except Exception as exc:
        raise ValueError(f"Unreadable h5ad file: {path}: {exc}") from exc
