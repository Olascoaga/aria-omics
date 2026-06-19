"""V47 bulk ATAC opening: beta QC + peak-calling slice.

Bulk ATAC is dispatchable only behind the readiness acknowledgement gate. The
validated first slice is measured QC plus MACS3 peak calling; DA remains a
structured scaffold skip until the replicate-aware count-matrix lane closes.
"""

from pathlib import Path


class _FakeBulkAtacEnv:
    def __init__(self):
        self.calls = []

    def run_in_stack(self, *, stack, script_path, params):
        self.calls.append((stack, Path(script_path).name, params))
        if script_path.endswith("chromatin_qc.py"):
            return {
                "status": "success",
                "data_type": "bulk_ATAC",
                "n_samples": 1,
                "mito_fraction": 0.02,
                "dup_rate": 0.11,
                "frip": None,
                "tss_enrichment": None,
                "metrics_not_computed": [
                    "frip (requires called peaks - run chromatin_peaks first)",
                    "tss_enrichment (computed post peak calling)",
                ],
                "qc_complete": False,
                "pass_qc": None,
                "warnings": [],
            }
        if script_path.endswith("chromatin_peaks.py"):
            return {
                "status": "success",
                "data_type": "bulk_ATAC",
                "n_peaks": 12345,
                "peaks_path": "/tmp/bulk_atac_peaks.narrowPeak",
                "consensus_peaks_path": None,
                "frip": 0.31,
                "genome": params.get("genome"),
                "macs3_cmd": "macs3 callpeak ...",
                "warnings": [],
            }
        if script_path.endswith("chromatin_peak_counts.py"):
            return {
                "status": "success",
                "data_type": "bulk_ATAC",
                "validation_level": "scaffold",
                "analysis": "peak_count_matrix",
                "counting_method": "bedtools coverage -counts",
                "counts_matrix_path": "/tmp/bulk_atac_peak_counts.tsv",
                "sample_metadata_path": "/tmp/bulk_atac_samples.tsv",
                "n_peaks": 12345,
                "n_samples": len(params.get("files") or []),
                "sample_ids": ["sample"],
                "warnings": [],
            }
        raise AssertionError(script_path)


def _agent(env):
    from aria.agents.chromatin_agent import ChromatinAgent

    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = env
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_finding = lambda *args, **kwargs: None
    return agent


def test_bulk_atac_agent_runs_qc_and_peak_calling_but_skips_da(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    env = _FakeBulkAtacEnv()
    agent = _agent(env)

    result = agent._run_bulk_atac(
        "exp",
        {"genome": "hg38"},
        {"comparison": "treated_vs_control"},
        [str(bam)],
    )

    assert result["status"] == "done"
    findings = result["findings"]
    assert findings["qc"]["data_type"] == "bulk_ATAC"
    assert findings["qc"]["validation_level"] == "beta"
    assert findings["peaks"]["n_peaks"] == 12345
    assert findings["peaks"]["validation_level"] == "beta"
    assert findings["peak_counts"]["status"] == "success"
    assert findings["peak_counts"]["validation_level"] == "scaffold"
    assert findings["peak_counts"]["analysis"] == "peak_count_matrix"
    assert findings["differential_accessibility"]["status"] == "skipped"
    assert findings["differential_accessibility"]["reason"] == "bulk_atac_da_not_validated"
    assert findings["differential_accessibility"]["validation_level"] == "scaffold"

    assert [(stack, script) for stack, script, _ in env.calls] == [
        ("chromatin", "chromatin_qc.py"),
        ("chromatin", "chromatin_peaks.py"),
        ("chromatin", "chromatin_peak_counts.py"),
    ]
    peaks_params = env.calls[1][2]
    assert peaks_params["data_type"] == "bulk_ATAC"
    assert peaks_params["macs3_params"]["format"] == "BAMPE"
    count_params = env.calls[2][2]
    assert count_params["peaks_path"] == "/tmp/bulk_atac_peaks.narrowPeak"


def test_chromatin_peaks_does_not_fabricate_frip(monkeypatch):
    import subprocess
    from aria.scripts.chromatin_peaks import _compute_frip

    def boom(*args, **kwargs):
        raise FileNotFoundError("tool missing")

    monkeypatch.setattr(subprocess, "run", boom)
    assert _compute_frip("/tmp/x.bam", "/tmp/x.narrowPeak") is None


def test_bulk_atac_narrator_surfaces_qc_and_peak_calling_under_strict_gate():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.render_blocks import render_blocks

    findings = {
        "bulk_ATAC": {
            "status": "done",
            "findings": {
                "qc": {
                    "status": "success",
                    "data_type": "bulk_ATAC",
                    "n_samples": 2,
                    "mito_fraction": 0.03,
                    "dup_rate": 0.12,
                    "metrics_not_computed": [],
                    "qc_complete": False,
                    "pass_qc": None,
                },
                "peaks": {
                    "status": "success",
                    "data_type": "bulk_ATAC",
                    "n_peaks": 23456,
                    "genome": "hg38",
                    "peaks_path": "/tmp/bulk_atac_peaks.narrowPeak",
                    "frip": None,
                    "warnings": [],
                },
                "peak_counts": {
                    "status": "success",
                    "data_type": "bulk_ATAC",
                    "counting_method": "bedtools coverage -counts",
                    "counts_matrix_path": "/tmp/bulk_atac_peak_counts.tsv",
                    "sample_metadata_path": "/tmp/bulk_atac_samples.tsv",
                    "n_peaks": 23456,
                    "n_samples": 2,
                    "warnings": [],
                },
            },
        }
    }

    blocks = ChromatinNarrator().collect(
        "chromatin_agent", {"status": "done", "findings": findings}, {})
    ids = {block.id for block in blocks}
    assert "chromatin.qc.bulk_ATAC" in ids
    assert "chromatin.peak_calling.bulk_ATAC" in ids
    assert "chromatin.peak_count_matrix.bulk_ATAC" in ids
    html = render_blocks(blocks, strict=True)
    assert "Bulk ATAC QC measured 2 sample" in html
    assert "MACS3 peak calling for bulk ATAC identified 23456 peaks" in html
    assert "Bulk ATAC peak counting produced a matrix with 23456 peaks across 2 samples" in html
