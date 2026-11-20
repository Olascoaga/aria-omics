"""
ARIA Hi-C constants
-------------------
Shared constants for Hi-C / Micro-C resolution sizing. Lives in `utils/` so
both the orchestrating agent (`aria.agents.genome_arch_agent`) and the
subprocess-only inspection script (`aria.scripts.hic_inspect`) can import
without violating the scripts-cannot-import-agents contract.
"""

from __future__ import annotations


# Approximate RAM requirements for human genome (hg38) full-genome analysis.
# Scale down ~10x for mouse, ~50x for drosophila (see GENOME_SCALE in
# genome_arch_agent).
RAM_ESTIMATES_GB: dict[int, float] = {
    1_000_000: 0.05,    # 1Mb
    500_000:   0.2,     # 500kb
    100_000:   1.0,     # 100kb
    40_000:    8.0,     # 40kb
    10_000:    100.0,   # 10kb
    5_000:     400.0,   # 5kb
}


def estimate_ram_gb(resolution_bp: int) -> float:
    """Estimate human full-genome dense matrix RAM (GB), monotonic O(1/res^2)."""
    r = int(resolution_bp)
    if r in RAM_ESTIMATES_GB:
        return RAM_ESTIMATES_GB[r]
    if r <= 0:
        return 1000.0
    return round(100.0 * ((10_000.0 / float(r)) ** 2), 4)


def list_mcool_resolutions(path: str) -> dict[int, str]:
    """Return ``{binsize_bp: internal_cooler_group}`` for a ``.mcool`` file.

    Supports both modern ``/resolutions/<bp>`` layouts and legacy zoom-level
    ``.mcool`` paths (``/0``, ``/1``, ...) by querying ``cooler.Cooler.binsize``.
    """
    import cooler

    groups = cooler.fileops.list_coolers(str(path))
    out: dict[int, str] = {}
    fallback_groups: list[str] = []
    for grp in groups:
        grp_s = str(grp)
        tail = grp_s.rstrip("/").split("/")[-1]
        if "/resolutions/" in grp_s and tail.isdigit():
            res = int(tail)
            if res > 0:
                out.setdefault(res, grp_s)
                continue
        fallback_groups.append(grp_s)
    for grp_s in fallback_groups:
        try:
            clr = cooler.Cooler(f"{path}::{grp_s}")
            if clr.binsize is not None and int(clr.binsize) > 0:
                out.setdefault(int(clr.binsize), grp_s)
        except Exception:
            continue
    return dict(sorted(out.items()))
