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
