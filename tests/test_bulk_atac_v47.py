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
        if script_path.endswith("chromatin_motifs.py"):
            return {
                "status": "success",
                "ran": True,
                "reason": None,
                "method": params.get("method", "hypergeometric"),
                "genome_fasta": "/tmp/hg38.fa",
                "motif_source": {"collection": "JASPAR2024_CORE_vertebrates",
                                 "release": "2024", "n_motifs": 100},
                "n_groups": len(params.get("regions") or {}),
                "per_group": {g: {"n_regions": len(v), "n_enriched": 3}
                              for g, v in (params.get("regions") or {}).items()},
                "output_csv": "/tmp/bulk_atac_motifs/motif_enrichment.csv",
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
                "motifs": {
                    "status": "success", "ran": True, "data_type": "bulk_ATAC",
                    "method": "hypergeometric",
                    "motif_source": {
                        "collection": "JASPAR2024_CORE_vertebrates",
                        "release": "2024"},
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
    # motif interpretation params (associative, both directions, offline)
    assert "TF motif enrichment" in blob
    assert "both conditions, no one-sided pruning" in blob
    assert "JASPAR2024_CORE_vertebrates" in blob
    assert "associative database match" in blob
    # honesty: ARIA does not realign, and no scATAC marker language leaks
    assert "does not realign" in blob
    assert "Wilcoxon" not in blob
    assert "per-cluster" not in blob


def _write_da_full_csv(path, rows):
    """rows: list of (peak, log2fc, significant_bool[, padj])."""
    import csv
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["peak", "baseMean", "log2FoldChange", "pvalue",
                         "padj", "significant"])
        for row in rows:
            peak, lfc, sig = row[:3]
            padj = row[3] if len(row) > 3 else 0.01
            writer.writerow([peak, 100.0, lfc, 0.001, padj, str(bool(sig))])


def test_bulk_da_motif_regions_splits_both_directions(tmp_path):
    from aria.agents.chromatin_agent import _bulk_da_motif_regions

    csv_path = tmp_path / "bulk_atac_da_K562_vs_GM12878_full.csv"
    _write_da_full_csv(csv_path, [
        ("chr1:100-200", 2.0, True),     # up in K562 (test)
        ("chr1:300-400", -1.5, True),    # up in GM12878 (reference)
        ("chr2:500-600", 0.1, False),    # not significant -> background only
        ("chr2:700-800", -3.0, True),    # up in GM12878
    ])
    comparisons = [{
        "test": "K562", "reference": "GM12878", "status": "success",
        "full_results_csv": str(csv_path),
    }]

    regions, background, warnings = _bulk_da_motif_regions(comparisons)

    assert set(regions) == {"K562_vs_GM12878::up_in_K562",
                            "K562_vs_GM12878::up_in_GM12878"}
    assert regions["K562_vs_GM12878::up_in_K562"] == ["chr1:100-200"]
    assert regions["K562_vs_GM12878::up_in_GM12878"] == [
        "chr2:700-800", "chr1:300-400"]
    # background is the FULL tested-peak universe incl. the non-significant peak
    assert set(background) == {"chr1:100-200", "chr1:300-400",
                               "chr2:500-600", "chr2:700-800"}
    assert warnings == []


def test_bulk_da_motif_regions_honest_when_no_csv():
    from aria.agents.chromatin_agent import _bulk_da_motif_regions

    regions, background, warnings = _bulk_da_motif_regions([
        {"test": "a", "reference": "b", "status": "success",
         "full_results_csv": "/nonexistent/path.csv"},
        {"test": "c", "reference": "d", "status": "failed"},
    ])
    assert regions == {}
    assert background == []
    assert any("not found" in w for w in warnings)


def test_bulk_da_motif_regions_caps_by_significance_not_csv_order(tmp_path):
    from aria.agents.chromatin_agent import _bulk_da_motif_regions

    csv_path = tmp_path / "bulk_atac_da_full.csv"
    rows = [
        ("chr1:1-2", 8.0, True, 0.049),
        ("chr1:3-4", 9.0, True, 0.048),
    ]
    rows.extend(
        (f"chr2:{100 + i * 10}-{105 + i * 10}", 1.1, True, 1e-8 + i * 1e-10)
        for i in range(5000)
    )
    _write_da_full_csv(csv_path, rows)

    regions, background, warnings = _bulk_da_motif_regions([{
        "test": "GroupA", "reference": "GroupB", "status": "success",
        "full_results_csv": str(csv_path),
    }])

    picked = regions["GroupA_vs_GroupB::up_in_GroupA"]
    assert len(picked) == 5000
    assert "chr1:1-2" not in picked
    assert "chr1:3-4" not in picked
    assert picked[0] == "chr2:100-105"
    assert len(background) == 5002
    assert any("ranking by padj" in w for w in warnings)


def test_bulk_atac_motif_enrichment_dispatched_after_da(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    csv_path = tmp_path / "da_full.csv"
    _write_da_full_csv(csv_path, [
        ("chr1:100-200", 2.0, True),
        ("chr1:300-400", -2.0, True),
    ])

    class _Env(_FakeBulkAtacEnv):
        def run_in_stack(self, *, stack, script_path, params):
            if script_path.endswith("chromatin_bulk_diffacc.py"):
                self.calls.append((stack, Path(script_path).name, params))
                return {
                    "status": "success", "ran": True, "data_type": "bulk_ATAC",
                    "validation_level": "beta",
                    "analysis": "differential_accessibility",
                    "method": "replicate-level DESeq2 over bulk ATAC peak counts",
                    "n_comparisons_success": 1, "n_replicate_samples": 4,
                    "n_peaks_tested": 2, "n_sig_total": 2,
                    "padj_max": 0.05, "lfc_min": 0.5,
                    "comparisons": [{
                        "test": "K562", "reference": "GM12878",
                        "status": "success", "full_results_csv": str(csv_path),
                    }],
                    "warnings": [],
                }
            return super().run_in_stack(
                stack=stack, script_path=script_path, params=params)

    env = _Env()
    agent = _agent(env)
    result = agent._run_bulk_atac(
        "exp",
        {"genome": "hg38", "comparisons": [["K562", "GM12878"]],
         "min_replicates_per_condition": 2},
        {"comparison": "K562_vs_GM12878"},
        [str(bam)],
    )

    motifs = result["findings"].get("motifs")
    assert motifs is not None and motifs["ran"] is True
    assert motifs["data_type"] == "bulk_ATAC"
    assert motifs["validation_level"] == "beta"

    motif_call = [c for c in env.calls if c[1] == "chromatin_motifs.py"][0]
    # motif enrichment runs on the chromatin stack (snapatac2 + genome FASTA)
    assert motif_call[0] == "chromatin"
    mp = motif_call[2]
    assert mp["assembly"] == "hg38"
    assert set(mp["regions"]) == {"K562_vs_GM12878::up_in_K562",
                                  "K562_vs_GM12878::up_in_GM12878"}
    assert "chr1:100-200" in mp["background"]
    assert (mp["foreground_truncation_strategy"]
            == "rank_by_padj_then_abs_log2fc_before_cap")
    # no data_path/da_csv: bulk ATAC feeds explicit regions + background
    assert "data_path" not in mp and "da_csv" not in mp


def test_bulk_atac_motif_skipped_when_no_da_peaks(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    # Shared fake DA result carries comparisons WITHOUT full_results_csv ->
    # no DA peaks to interpret -> motif step is honestly skipped (no dispatch).
    env = _FakeBulkAtacEnv()
    agent = _agent(env)
    result = agent._run_bulk_atac(
        "exp",
        {"genome": "hg38", "comparisons": [["treated", "control"]]},
        {"comparison": "treated_vs_control"},
        [str(bam)],
    )
    assert "motifs" not in result["findings"]
    assert not any(c[1] == "chromatin_motifs.py" for c in env.calls)
