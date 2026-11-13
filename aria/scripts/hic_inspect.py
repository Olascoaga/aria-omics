"""ARIA Hi-C File Inspector (ADR-082).

validation_level: production (stack="hic").

Reads Hi-C / Micro-C metadata WITHOUT loading dense contact matrices into RAM.
Supports ``.cool``, ``.mcool``, ``.hic``, ``.pairs`` / ``.pairs.gz``, and
aligned chimeric ``.bam`` files.
"""

from __future__ import annotations

import gzip
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from aria.scripts._base import run_script


def hic_inspect(params: dict) -> dict:
    files = [str(f) for f in (params.get("files") or []) if f]
    genome = str(params.get("genome") or "hg38")

    valid_files = [f for f in files if Path(f).exists()]
    if not valid_files:
        return {
            "status": "error",
            "error_type": "FileNotFound",
            "details": "No valid Hi-C files found.",
        }

    all_resolutions: set[int] = set()
    all_formats: set[str] = set()
    all_chroms: list[str] = []
    seen_chroms: set[str] = set()
    per_file: list[dict] = []
    needs_preprocessing = False

    for fpath in valid_files:
        info = _inspect_file(fpath)
        per_file.append({"file": fpath, **info})
        all_resolutions.update(int(r) for r in (info.get("resolutions") or []))
        all_formats.add(str(info.get("format") or "unknown"))
        if info.get("needs_preprocessing"):
            needs_preprocessing = True
        for ch in info.get("chromosomes") or []:
            if ch not in seen_chroms:
                seen_chroms.add(ch)
                all_chroms.append(ch)

    sorted_res = sorted(all_resolutions)
    from aria.utils.hic_constants import RAM_ESTIMATES_GB

    estimates = {r: RAM_ESTIMATES_GB.get(r, 1.0) for r in sorted_res}
    return {
        "status": "success",
        "data_type": "HiC",
        "analysis": "hic_inspect",
        "validation_level": "production",
        "genome": genome,
        "files": valid_files,
        "formats": sorted(all_formats),
        "available_resolutions": sorted_res,
        "chromosomes": all_chroms,
        "estimated_sizes_gb": estimates,
        "needs_preprocessing": bool(needs_preprocessing),
        "per_file": per_file,
    }


def _inspect_file(fpath: str) -> dict:
    name = Path(fpath).name.lower()
    if name.endswith(".mcool"):
        return _inspect_mcool(fpath)
    if name.endswith(".cool"):
        return _inspect_cool(fpath)
    if name.endswith(".hic"):
        return _inspect_hic(fpath)
    if name.endswith((".pairs", ".pairs.gz")):
        return _inspect_pairs(fpath)
    if name.endswith((".bam", ".sam")):
        return {
            "format": "bam",
            "resolutions": [],
            "chromosomes": [],
            "needs_preprocessing": True,
        }
    return {"format": "unknown", "resolutions": [], "chromosomes": []}


def _inspect_mcool(fpath: str) -> dict:
    try:
        import cooler

        resolutions = cooler.fileops.list_coolers(fpath)
        res_values: list[int] = []
        for r in resolutions:
            tail = r.split("/")[-1]
            if tail.isdigit():
                res_values.append(int(tail))
        chroms: list[str] = []
        if res_values:
            finest = f"{fpath}::resolutions/{min(res_values)}"
            clr = cooler.Cooler(finest)
            chroms = list(clr.chromnames)
        return {
            "format": "mcool",
            "resolutions": sorted(res_values),
            "chromosomes": chroms,
            "needs_preprocessing": False,
        }
    except Exception as exc:
        return {
            "format": "mcool",
            "resolutions": [],
            "chromosomes": [],
            "error": str(exc)[:200],
        }


def _inspect_cool(fpath: str) -> dict:
    try:
        import cooler

        clr = cooler.Cooler(fpath)
        return {
            "format": "cool",
            "resolutions": [int(clr.binsize)],
            "chromosomes": list(clr.chromnames),
            "needs_preprocessing": False,
        }
    except Exception as exc:
        return {
            "format": "cool",
            "resolutions": [],
            "chromosomes": [],
            "error": str(exc)[:200],
        }


def _inspect_hic(fpath: str) -> dict:
    try:
        import hicstraw

        hic = hicstraw.HiCFile(fpath)
        return {
            "format": "hic",
            "resolutions": sorted(int(r) for r in hic.getResolutions()),
            "chromosomes": [
                str(c.name) for c in hic.getChromosomes()
                if str(c.name) not in {"All", "ALL", "all"}
            ],
            "needs_preprocessing": True,
        }
    except Exception as exc:
        return {
            "format": "hic",
            "resolutions": [],
            "chromosomes": [],
            "error": str(exc)[:200],
        }


def _inspect_pairs(fpath: str) -> dict:
    chroms: list[str] = []
    opener = gzip.open if str(fpath).lower().endswith(".gz") else open
    try:
        with opener(fpath, "rt", encoding="utf-8", errors="replace") as fh:
            for idx, line in enumerate(fh):
                if idx > 500:
                    break
                if line.startswith("#chromsize:"):
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        chroms.append(parts[1])
                elif not line.startswith("#"):
                    break
    except Exception:
        pass
    return {
        "format": "pairs",
        "resolutions": [],
        "chromosomes": chroms,
        "needs_preprocessing": True,
    }


if __name__ == "__main__":
    run_script(hic_inspect)
