"""
ARIA Hi-C File Inspector
-------------------------
Reads Hi-C file metadata WITHOUT loading any matrix data.
Fast pre-flight check for available resolutions, chromosomes,
and file format. Executed before any analysis.

Input params:
    files:  list — .cool, .mcool, or .hic files
    genome: str

Output:
    {
      "status":               "success",
      "files":                list,
      "formats":              [str],
      "available_resolutions": [int],  — sorted ascending
      "chromosomes":          [str],
      "estimated_sizes_gb":   {resolution: gb}
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def hic_inspect(params: dict) -> dict:
    from pathlib import Path

    files  = params.get("files", [])
    genome = params.get("genome", "hg38")

    valid_files = [f for f in files if Path(f).exists()]
    if not valid_files:
        return {
            "status":     "error",
            "error_type": "FileNotFound",
            "details":    "No valid Hi-C files found.",
        }

    all_resolutions = set()
    all_formats     = set()
    all_chroms      = set()

    for fpath in valid_files:
        info = _inspect_file(fpath)
        all_resolutions.update(info.get("resolutions", []))
        all_formats.add(info.get("format", "unknown"))
        all_chroms.update(info.get("chromosomes", []))

    sorted_res = sorted(all_resolutions)

    # Estimate RAM for each resolution (human genome)
    from aria.agents.genome_arch_agent import RAM_ESTIMATES_GB
    estimates = {r: RAM_ESTIMATES_GB.get(r, 100.0) for r in sorted_res}

    return {
        "status":                "success",
        "files":                 valid_files,
        "formats":               list(all_formats),
        "available_resolutions": sorted_res,
        "chromosomes":           sorted(all_chroms)[:25],
        "estimated_sizes_gb":    estimates,
    }


def _inspect_file(fpath: str) -> dict:
    """Inspect a single Hi-C file — metadata only."""
    from pathlib import Path
    ext = "".join(Path(fpath).suffixes)

    if ".mcool" in ext:
        return _inspect_mcool(fpath)
    elif ".cool" in ext:
        return _inspect_cool(fpath)
    elif ext == ".hic":
        return _inspect_hic(fpath)
    else:
        return {"format": "unknown", "resolutions": [], "chromosomes": []}


def _inspect_mcool(fpath: str) -> dict:
    try:
        import cooler
        resolutions = cooler.fileops.list_coolers(fpath)
        # resolutions are paths like "/resolutions/1000000"
        res_values  = []
        for r in resolutions:
            try:
                val = int(r.split("/")[-1])
                res_values.append(val)
            except ValueError:
                pass
        # Get chromosomes from finest resolution
        if res_values:
            finest = f"{fpath}::resolutions/{min(res_values)}"
            clr    = cooler.Cooler(finest)
            chroms = list(clr.chromnames)
        else:
            chroms = []
        return {
            "format":      "mcool",
            "resolutions": sorted(res_values),
            "chromosomes": chroms,
        }
    except Exception as e:
        return {"format": "mcool", "resolutions": [], "chromosomes": [],
                "error": str(e)[:100]}


def _inspect_cool(fpath: str) -> dict:
    try:
        import cooler
        clr = cooler.Cooler(fpath)
        return {
            "format":      "cool",
            "resolutions": [clr.binsize],
            "chromosomes": list(clr.chromnames),
        }
    except Exception as e:
        return {"format": "cool", "resolutions": [], "chromosomes": [],
                "error": str(e)[:100]}


def _inspect_hic(fpath: str) -> dict:
    try:
        import hicstraw
        hic = hicstraw.HiCFile(fpath)
        return {
            "format":      "hic",
            "resolutions": hic.getResolutions(),
            "chromosomes": [c.name for c in hic.getChromosomes()
                            if c.name != "All"],
        }
    except Exception as e:
        return {"format": "hic", "resolutions": [], "chromosomes": [],
                "error": str(e)[:100]}


if __name__ == "__main__":
    run_script(hic_inspect)
