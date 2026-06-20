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
        if script_path.endswith("chromatin_bulk_diffacc.py"):
            return {
                "status": "success",
                "ran": True,
                "data_type": "bulk_ATAC",
                "validation_level": "beta",
                "analysis": "differential_accessibility",
                "method": "replicate-level DESeq2 over bulk ATAC peak counts",
                "n_comparisons_success": 1,
                "n_replicate_samples": 6,
                "n_peaks_tested": 12345,
                "n_sig_total": 42,
                "padj_max": 0.05,
                "lfc_min": 0.5,
                "output_csv": "/tmp/bulk_atac_da_summary.csv",
                "replicate_counts_path": "/tmp/bulk_atac_replicate_counts.tsv",
                "comparisons": [{
                    "test": "treated", "reference": "control",
                    "status": "success", "n_sig": 42,
                    "n_up": 20, "n_down": 22,
                }],
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


def test_bulk_atac_agent_runs_qc_peaks_counts_and_da(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    env = _FakeBulkAtacEnv()
    agent = _agent(env)

    result = agent._run_bulk_atac(
        "exp",
        {
            "genome": "hg38",
            "sample_metadata": {
                "sample": {"condition": "treated", "replicate": "r1"},
            },
            "comparisons": [["treated", "control"]],
        },
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
    assert findings["differential_accessibility"]["status"] == "success"
    assert findings["differential_accessibility"]["ran"] is True
    assert findings["differential_accessibility"]["validation_level"] == "beta"
    assert findings["differential_accessibility"]["n_sig_total"] == 42

    # DA dispatches to the rna stack: pydeseq2 lives in aria-rna-env, not the
    # chromatin env. The QC/peaks/counts steps stay in the chromatin stack.
    assert [(stack, script) for stack, script, _ in env.calls] == [
        ("chromatin", "chromatin_qc.py"),
        ("chromatin", "chromatin_peaks.py"),
        ("chromatin", "chromatin_peak_counts.py"),
        ("rna", "chromatin_bulk_diffacc.py"),
    ]
    peaks_params = env.calls[1][2]
    assert peaks_params["data_type"] == "bulk_ATAC"
    assert peaks_params["macs3_params"]["format"] == "BAMPE"
    count_params = env.calls[2][2]
    assert count_params["peaks_path"] == "/tmp/bulk_atac_peaks.narrowPeak"
    da_params = env.calls[3][2]
    assert da_params["counts_matrix_path"] == "/tmp/bulk_atac_peak_counts.tsv"
    assert da_params["comparisons"] == [["treated", "control"]]
    # No replicate override supplied -> the script's production floor applies.
    assert "min_replicates_per_condition" not in da_params


def test_bulk_atac_da_runs_in_rna_stack_with_replicate_override(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    env = _FakeBulkAtacEnv()
    agent = _agent(env)

    agent._run_bulk_atac(
        "exp",
        {
            "genome": "hg38",
            "sample_metadata": {
                "sample": {"condition": "treated", "replicate": "r1"},
            },
            "comparisons": [["treated", "control"]],
            "min_replicates_per_condition": 2,
        },
        {"comparison": "treated_vs_control"},
        [str(bam)],
    )

    da_call = [c for c in env.calls if c[1] == "chromatin_bulk_diffacc.py"][0]
    assert da_call[0] == "rna"
    assert da_call[2]["min_replicates_per_condition"] == 2


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
                "differential_accessibility": {
                    "status": "success",
                    "ran": True,
                    "data_type": "bulk_ATAC",
                    "method": "replicate-level DESeq2 over bulk ATAC peak counts",
                    "n_comparisons_success": 1,
                    "n_replicate_samples": 6,
                    "n_peaks_tested": 23456,
                    "n_sig_total": 42,
                    "padj_max": 0.05,
                    "lfc_min": 0.5,
                    "output_csv": "/tmp/bulk_atac_da_summary.csv",
                    "replicate_counts_path": "/tmp/bulk_atac_replicate_counts.tsv",
                    "comparisons": [{
                        "test": "treated", "reference": "control",
                        "status": "success", "n_sig": 42,
                    }],
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
    assert "chromatin.differential_accessibility.bulk_atac" in ids
    html = render_blocks(blocks, strict=True)
    assert "Bulk ATAC QC measured 2 sample" in html
    assert "MACS3 peak calling for bulk ATAC identified 23456 peaks" in html
    assert "Bulk ATAC peak counting produced a matrix with 23456 peaks across 2 samples" in html
    assert "Bulk ATAC DESeq2 differential accessibility ran 1 comparison" in html


def test_bulk_atac_methods_report_caller_counting_and_model_params():
    """The methods section reports real aligner-input/caller/counting/model
    parameters for bulk ATAC and never falls back to the scATAC Wilcoxon
    marker language."""
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    findings = {
        "bulk_ATAC": {
            "status": "done",
            "findings": {
                "qc": {
                    "status": "success", "data_type": "bulk_ATAC",
                    "n_samples": 4, "mito_fraction": 0.03,
                },
                "peaks": {
                    "status": "success", "data_type": "bulk_ATAC",
                    "n_peaks": 124544, "genome": "hg38",
                    "consensus_peaks_path": "/tmp/consensus.narrowPeak",
                    "macs3_cmd": ("macs3 callpeak -t a.bam b.bam -f BAMPE "
                                  "-g hs --nomodel --keep-dup all -q 0.05"),
                },
                "peak_counts": {
                    "status": "success", "data_type": "bulk_ATAC",
                    "counting_method": "bedtools coverage -sorted -counts",
                    "n_peaks": 124544, "n_samples": 4,
                },
                "differential_accessibility": {
                    "status": "success", "ran": True, "data_type": "bulk_ATAC",
                    "method": "replicate-level DESeq2 over bulk ATAC peak counts",
                    "n_replicate_samples": 4,
                    "condition_col": "cell_line", "replicate_col": "replicate",
                    "covariates_adjusted": [],
                    "padj_max": 0.05, "lfc_min": 0.5,
                    "min_replicates_per_condition": 2,
                    "comparisons": [{
                        "test": "K562", "reference": "GM12878",
                        "status": "success",
                        "fitted_design_formula": "~ cell_line",
                        "low_power_warning": True,
                    }],
                },
            },
        }
    }
    agent_result = {"status": "done", "findings": findings}
    methods = ChromatinNarrator().methods("chromatin_agent", agent_result, {})
    blob = " ".join(methods)

    # caller params (incl. the exact reproducible MACS3 command)
    assert "MACS3" in blob and "macs3 callpeak" in blob and "hg38" in blob
    assert "consensus peak universe" in blob
    # counting params
    assert "bedtools coverage -sorted -counts" in blob
    assert "peak-by-sample count matrix" in blob
    # model params
    assert "replicate-level DESeq2" in blob
    assert "cell_line" in blob and "replicate" in blob
    assert "~ cell_line" in blob
    assert "|log2FC| >= 0.5" in blob and "adjusted p <= 0.05" in blob
    assert "At least 2 replicates per condition" in blob
    assert "low-power warning" in blob
    # honesty: ARIA does not realign, and no scATAC marker language leaks
    assert "does not realign" in blob
    assert "Wilcoxon" not in blob
    assert "per-cluster" not in blob
