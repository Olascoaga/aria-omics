"""Regression tests for D7 layer-contract fixes in Hi-C."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


def test_d7_011_hic_adapter_falls_through_to_qc_n_samples_when_loops_omits_it(tmp_path: Path):
    """D7-011: hic_topology._compute_loops does not write 'n_samples' in the
    loops dict (sample count lives in qc['n_samples']). hypothesis/adapters/hic.py
    must not default loops.get('n_samples', 1) and falsely tag n>=2 Hi-C cohorts
    with LOW_REPLICATION."""
    hic_py = Path("/home/medusa/anaconda3/envs/aria-hic-env/bin/python")
    if not hic_py.is_file():
        pytest.skip("aria-hic-env not available")
    pytest.importorskip("litellm")
    pytest.importorskip("pydantic")
    from aria.agents.narrative.hypothesis.adapters.hic import _low_replication, hic_evidence
    from aria.agents.narrative.hypothesis.caveats import LOW_REPLICATION
    from aria.utils.environment_manager import env_manager

    cool_path = tmp_path / "sample.cool"
    build_code = f"""
import cooler, numpy as np, pandas as pd
res = 25000
n_bins = 40
bins = pd.DataFrame({{
    'chrom': ['chr1'] * n_bins,
    'start': [i * res for i in range(n_bins)],
    'end': [(i + 1) * res for i in range(n_bins)],
    'weight': np.ones(n_bins, dtype=float),
}})
rows = [(i, j, max(1, int(100 / (1 + j - i)))) for i in range(n_bins) for j in range(i, min(n_bins, i + 15))]
pixels = pd.DataFrame(rows, columns=['bin1_id', 'bin2_id', 'count'])
cooler.create_cooler({str(cool_path)!r}, bins, pixels, dtypes={{'count': int}})
"""
    subprocess.run([str(hic_py), "-c", build_code], check=True)
    loops_res = env_manager.run_in_stack(
        "hic",
        "aria/scripts/hic_topology.py",
        {
            "files": [str(cool_path)],
            "analysis": "loops",
            "resolution": 25000,
            "fdr_threshold": 0.1,
            "output_dir": str(tmp_path / "loops_out"),
        },
    )
    assert loops_res["status"] == "success", loops_res
    assert "n_samples" not in loops_res

    loops_with_hit = {
        **loops_res,
        "n_loops": 1,
        "top_loops": [{
            "chrom1": "chr1", "start1": 25000, "end1": 50000,
            "chrom2": "chr1", "start2": 200000, "end2": 225000,
            "count": 80, "qval": 0.001,
            "anchor1_genes": "MYC", "loop_class": "enhancer_promoter",
        }],
    }
    findings_n4 = {"qc": {"status": "success", "n_samples": 4}, "loops": loops_with_hit}
    assert _low_replication(findings_n4) is False
    sigs_n4 = hic_evidence({"genome_arch_agent": {"status": "done", "findings": findings_n4}})
    assert sigs_n4 and LOW_REPLICATION not in sigs_n4[0].caveats_inherited

    # Still fail-closed (LOW_REPLICATION=True) when qc.n_samples < 2 or unknown:
    findings_n1 = {"qc": {"status": "success", "n_samples": 1}, "loops": loops_with_hit}
    assert _low_replication(findings_n1) is True
    assert _low_replication({"loops": loops_with_hit}) is True
