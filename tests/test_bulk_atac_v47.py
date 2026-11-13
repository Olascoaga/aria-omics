"""Bulk ATAC lane: production spine (ADR-077).

FASTQ -> BAM -> QC -> MACS3 -> peak-count matrix -> replicate-gated DESeq2 DA ->
motifs / annotation / ORA are stamped ``production``; condition-level TOBIAS
footprinting stays ``beta``. These tests pin the dispatch order, the design-
contract trigger, pre-computed matrix intake, the C3 ledger node for raw
preprocessing, strict W-CLAIM support, and the dual PNG+SVG bulk figures.
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
        if script_path.endswith("chromatin_peak_annotation.py"):
            return {
                "status": "success",
                "ran": True,
                "data_type": "bulk_ATAC",
                "validation_level": "beta",
                "analysis": "peak_annotation",
                "method": "nearest-TSS genomic annotation",
                "gtf": "annotation.gtf",
                "promoter_upstream": 2000,
                "promoter_downstream": 2000,
                "feature_distribution_overall": {
                    "Promoter": 12, "Intronic": 20, "Distal Intergenic": 10},
                "comparisons": [{"test": "treated", "reference": "control",
                                 "n_annotated": 42,
                                 "annotation_csv": "/tmp/ann.csv"}],
                "warnings": [],
            }
        if script_path.endswith("chromatin_peak_ora.py"):
            return {
                "status": "success",
                "ran": True,
                "data_type": "bulk_ATAC",
                "validation_level": "beta",
                "analysis": "peak_ora",
                "method": "functional ORA (local hypergeometric)",
                "organism": "Homo sapiens",
                "gtf": "annotation.gtf",
                "comparisons": [{"test": "treated", "reference": "control",
                                 "ora_method": "local_hypergeometric",
                                 "pathways_summary": {"GO_BP": 5},
                                 "figures": {}}],
                "figures": {},
                "warnings": [],
            }
        if script_path.endswith("chromatin_footprint_tobias.py"):
            return {
                "status": "success",
                "ran": False,
                "reason": "TOBIAS not installed",
                "method": "tobias",
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
    assert findings["qc"]["validation_level"] == "production"
    assert findings["peaks"]["n_peaks"] == 12345
    assert findings["peaks"]["validation_level"] == "production"
    assert findings["peak_counts"]["status"] == "success"
    assert findings["peak_counts"]["validation_level"] == "production"
    assert findings["peak_counts"]["analysis"] == "peak_count_matrix"
    assert findings["differential_accessibility"]["status"] == "success"
    assert findings["differential_accessibility"]["ran"] is True
    assert (findings["differential_accessibility"]["validation_level"]
            == "production")
    assert findings["peak_annotation"]["validation_level"] == "production"
    assert findings["peak_ora"]["validation_level"] == "production"
    # TOBIAS footprinting promoted to production in ADR-078.
    assert findings["footprinting"]["validation_level"] == "production"
    # Post-peak-calling FRiP closes the QC gate (TSSe still unmeasured here).
    assert findings["qc"]["frip"] == 0.31
    assert not any(str(m).startswith("frip")
                   for m in findings["qc"]["metrics_not_computed"])
    assert findings["qc"]["pass_qc"] is None
    assert findings["differential_accessibility"]["n_sig_total"] == 42

    # DA dispatches to the rna stack: pydeseq2 lives in aria-rna-env, not the
    # chromatin env. The QC/peaks/counts steps stay in the chromatin stack.
    assert [(stack, script) for stack, script, _ in env.calls] == [
        ("chromatin", "chromatin_qc.py"),
        ("chromatin", "chromatin_peaks.py"),
        ("chromatin", "chromatin_peak_counts.py"),
        ("rna", "chromatin_bulk_diffacc.py"),
        # B2/B3: genomic peak annotation + functional ORA run in the rna stack
        # after DA. (Motif enrichment self-skips here: the mock DA comparisons
        # carry no per-peak data, so _bulk_da_motif_regions yields no regions.)
        ("rna", "chromatin_peak_annotation.py"),
        ("rna", "chromatin_peak_ora.py"),
        ("tobias", "chromatin_footprint_tobias.py"),
    ]
    assert findings["peak_annotation"]["analysis"] == "peak_annotation"
    assert findings["peak_ora"]["analysis"] == "peak_ora"
    assert findings["footprinting"]["analysis"] == "differential_tf_footprinting"
    assert findings["footprinting"]["ran"] is False
    peaks_params = env.calls[1][2]
    assert peaks_params["data_type"] == "bulk_ATAC"
    assert peaks_params["macs3_params"]["format"] == "BAMPE"
    count_params = env.calls[2][2]
    assert count_params["peaks_path"] == "/tmp/bulk_atac_peaks.narrowPeak"
    da_params = env.calls[3][2]
    assert da_params["counts_matrix_path"] == "/tmp/bulk_atac_peak_counts.tsv"
    assert da_params["comparisons"] == [["treated", "control"]]
    footprint_params = env.calls[-1][2]
    assert footprint_params["mode"] == "bulk"
    assert footprint_params["group_a"] == "treated"
    assert footprint_params["group_b"] == "control"
    assert footprint_params["condition_bams"] == {"treated": [str(bam)]}
    assert footprint_params["peaks_bed"] == "/tmp/bulk_atac_peaks.narrowPeak"
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
                    "covariates_dropped": [{
                        "covariate": "batch",
                        "reason": (
                            "not constant within biological replicate after "
                            "technical-sample aggregation"),
                    }],
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
    assert "not adjusted for after replicate aggregation: batch" in blob
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
    assert motifs["validation_level"] == "production"

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
    # F7: the caller declares the REAL modality so the motif engine does not
    # resolve thresholds under the hardcoded scATAC policy.
    assert mp["modality"] == "bulk_ATAC"
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


# ══════════════════════════════════════════════════════════════════════════════
# ADR-077 Track 2: analytical closure of the bulk ATAC lane
# ══════════════════════════════════════════════════════════════════════════════

import json

import pytest


def _scripts(env):
    return [script for _stack, script, _params in env.calls]


@pytest.mark.parametrize("exp_extra", [
    {"design": {"main_factor": "condition",
                "contrasts": [{"test": "treated", "reference": "control"}]}},
    {"comparisons": [{"test": "treated", "reference": "control"}]},
])
def test_design_contract_comparison_triggers_counts_and_da_without_intent(
        tmp_path, exp_extra):
    """A confirmed design contrast (or exp_ctx comparisons) must reach DA even
    when the parsed intent carries no free-text comparison."""
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    env = _FakeBulkAtacEnv()
    agent = _agent(env)

    result = agent._run_bulk_atac(
        "exp", {"genome": "hg38", **exp_extra}, {}, [str(bam)])

    findings = result["findings"]
    assert "chromatin_peak_counts.py" in _scripts(env)
    assert "chromatin_bulk_diffacc.py" in _scripts(env)
    assert findings["differential_accessibility"]["ran"] is True
    da_params = next(p for _s, n, p in env.calls
                     if n == "chromatin_bulk_diffacc.py")
    assert da_params["comparisons"] == [
        {"test": "treated", "reference": "control"}]


def test_no_comparison_anywhere_skips_counts_and_da(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    env = _FakeBulkAtacEnv()
    result = _agent(env)._run_bulk_atac("exp", {"genome": "hg38"}, {},
                                         [str(bam)])
    assert _scripts(env) == ["chromatin_qc.py", "chromatin_peaks.py"]
    assert "differential_accessibility" not in result["findings"]


def test_design_columns_propagate_to_footprint_replicate_design(tmp_path):
    bams = []
    for name in ("a1", "a2", "b1", "b2"):
        path = tmp_path / f"{name}.bam"
        path.write_bytes(b"")
        bams.append(str(path))
    env = _FakeBulkAtacEnv()
    _agent(env)._run_bulk_atac(
        "exp",
        {
            "genome": "hg38",
            "design": {"main_factor": "cell_line", "replicate_col": "donor"},
            "comparisons": [["treated", "control"]],
            "sample_ids": ["a1", "a2", "b1", "b2"],
            "sample_metadata": {
                "a1": {"cell_line": "treated", "donor": "d1"},
                "a2": {"cell_line": "treated", "donor": "d2"},
                "b1": {"cell_line": "control", "donor": "d1"},
                "b2": {"cell_line": "control", "donor": "d2"},
            },
        },
        {},
        bams,
    )
    da_params = next(p for _s, n, p in env.calls
                     if n == "chromatin_bulk_diffacc.py")
    assert da_params["condition_col"] == "cell_line"
    assert da_params["replicate_col"] == "donor"
    fp = next(p for _s, n, p in env.calls
              if n == "chromatin_footprint_tobias.py")
    assert fp["condition_bams"] == {"treated": bams[:2], "control": bams[2:]}
    assert fp["replicate_bams"] == {
        "treated": {"d1": [bams[0]], "d2": [bams[1]]},
        "control": {"d1": [bams[2]], "d2": [bams[3]]},
    }


def test_qc_gate_closes_with_peak_frip_and_assets_are_forwarded(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")

    class _Env(_FakeBulkAtacEnv):
        def run_in_stack(self, *, stack, script_path, params):
            out = super().run_in_stack(stack=stack, script_path=script_path,
                                       params=params)
            if script_path.endswith("chromatin_qc.py"):
                out = {**out, "tss_enrichment": 7.2,
                       "metrics_not_computed": ["frip"]}
            if script_path.endswith("chromatin_peaks.py"):
                out = {**out, "frip_per_sample": {"sample": 0.31}}
            return out

    env = _Env()
    result = _agent(env)._run_bulk_atac(
        "exp",
        {"genome": "hg38", "gtf": "/ref/genes.gtf",
         "blacklist_bed": "/ref/blacklist.bed",
         "reference_fasta": "/ref/genome.fa"},
        {},
        [str(bam)],
    )
    qc = result["findings"]["qc"]
    assert qc["frip"] == 0.31 and qc["frip_per_sample"] == {"sample": 0.31}
    assert qc["metrics_not_computed"] == []
    assert qc["qc_complete"] is True and qc["pass_qc"] is True
    qc_params = env.calls[0][2]
    assert qc_params["gtf"] == "/ref/genes.gtf"
    assert qc_params["reference_fasta"] == "/ref/genome.fa"
    peaks_params = env.calls[1][2]
    assert peaks_params["blacklist_bed"] == "/ref/blacklist.bed"
    assert peaks_params["reference_fasta"] == "/ref/genome.fa"


def test_minimal_qc_params_without_optional_assets(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    env = _FakeBulkAtacEnv()
    _agent(env)._run_bulk_atac("exp", {"genome": "hg38"}, {}, [str(bam)])
    assert set(env.calls[0][2]) == {"data_type", "files", "genome"}
    assert "blacklist_bed" not in env.calls[1][2]


@pytest.mark.parametrize("suffix,sep", [(".tsv", "\t"), (".csv", ",")])
def test_precomputed_peak_count_matrix_intake_skips_bam_stages(
        tmp_path, suffix, sep):
    matrix = tmp_path / f"consensus_peak_counts{suffix}"
    samples = ["c1", "c2", "c3", "t1", "t2", "t3"]
    rows = [sep.join(["peak_id", *samples])]
    for i in range(5):
        rows.append(sep.join([f"chr1:{i * 100}-{i * 100 + 50}",
                              *[str(10 + i + j) for j in range(6)]]))
    matrix.write_text("\n".join(rows) + "\n")
    sample_metadata = {
        s: {"condition": "control" if s.startswith("c") else "treated",
            "replicate": f"r{s[1]}"}
        for s in samples
    }
    env = _FakeBulkAtacEnv()
    result = _agent(env)._run_bulk_atac(
        "exp",
        {"genome": "hg38", "sample_metadata": sample_metadata,
         "comparisons": [["treated", "control"]]},
        {},
        [str(matrix)],
    )

    assert result["status"] == "done"
    scripts = _scripts(env)
    for skipped in ("chromatin_qc.py", "chromatin_peaks.py",
                    "chromatin_peak_counts.py", "chromatin_footprint_tobias.py"):
        assert skipped not in scripts
    assert scripts[0] == "chromatin_bulk_diffacc.py"
    assert "chromatin_peak_annotation.py" in scripts
    assert "chromatin_peak_ora.py" in scripts

    pc = result["findings"]["peak_counts"]
    assert pc["status"] == "success" and pc["ran"] is True
    assert pc["source"] == "precomputed_matrix"
    assert pc["data_type"] == "bulk_ATAC"
    assert pc["validation_level"] == "production"
    assert pc["n_peaks"] == 5 and pc["n_samples"] == 6
    counts_path = Path(pc["counts_matrix_path"])
    assert counts_path.suffix == ".tsv"
    assert counts_path.read_text().splitlines()[0].split("\t")[0] == "peak_id"
    meta_lines = Path(pc["sample_metadata_path"]).read_text().splitlines()
    assert meta_lines[0].split("\t") == ["sample_id", "condition", "replicate"]
    assert len(meta_lines) == 7

    da_params = env.calls[0][2]
    assert da_params["counts_matrix_path"] == str(counts_path)
    assert da_params["sample_metadata_path"] == pc["sample_metadata_path"]
    assert (result["findings"]["differential_accessibility"]["validation_level"]
            == "production")


def test_precomputed_matrix_via_exp_ctx_path_and_missing_metadata(tmp_path):
    matrix = tmp_path / "counts.tsv"
    matrix.write_text("peak_id\ta\tb\nchr1:1-10\t1\t2\n")
    env = _FakeBulkAtacEnv()
    result = _agent(env)._run_bulk_atac(
        "exp", {"counts_matrix_path": str(matrix),
                "comparisons": [["a", "b"]]}, {}, [])
    assert env.calls == []
    da = result["findings"]["differential_accessibility"]
    assert da["status"] == "skipped"
    assert "sample_metadata_missing" in da["reason"]


def test_fragments_and_peak_tables_are_not_count_matrices():
    from aria.agents.chromatin_agent import ChromatinAgent

    is_matrix = ChromatinAgent._is_peak_count_matrix_path
    assert is_matrix("counts.tsv") and is_matrix("counts.csv")
    assert not is_matrix("atac_fragments.tsv.gz")
    assert not is_matrix("sample.fragments.tsv")
    assert not is_matrix("peaks.bed")
    assert not is_matrix("x_peaks.narrowPeak")
    assert ChromatinAgent._detect_precomputed_counts_matrix(
        {"counts_matrix_path": "c.tsv"}, ["s.bam"]) is None


class _FakeRawFastqBulkAtacEnv(_FakeBulkAtacEnv):
    def __init__(self, bam_dir):
        super().__init__()
        self.bam_dir = Path(bam_dir)

    def run_in_stack(self, *, stack, script_path, params):
        if script_path.endswith("rna_fastq_qc.py"):
            self.calls.append((stack, Path(script_path).name, params))
            return {"status": "success", "n_samples": 2,
                    "samples": [{"name": s, "paired": True,
                                 "r1_trimmed": f"{s}_R1.fq.gz",
                                 "r2_trimmed": f"{s}_R2.fq.gz"}
                                for s in ("ctl", "trt")],
                    "multiqc_report": None}
        if script_path.endswith("atac_align.py"):
            self.calls.append((stack, Path(script_path).name, params))
            return {
                "status": "success", "validation_level": "production",
                "aligner": "bwa-mem2", "read_layout": "paired",
                "filter": {"mapq": 30, "samtools_f": 2, "samtools_F": 1804,
                           "remove_mito": True},
                "bam_files": [
                    {"name": s, "sample_id": s, "status": "success",
                     "paired": True, "mapped_reads_post_filter": 1000 + i,
                     "bam": str(self.bam_dir / f"{s}.filtered.bam")}
                    for i, s in enumerate(("ctl", "trt"))
                ],
            }
        out = super().run_in_stack(stack=stack, script_path=script_path,
                                   params=params)
        if script_path.endswith("chromatin_qc.py"):
            out = {**out, "n_samples": 2, "mito_fraction": 0.03,
                   "dup_rate": 0.12}
        return out


def _raw_fastq_bulk_findings(tmp_path):
    fq_dir = tmp_path / "fastq"
    fq_dir.mkdir()
    fastqs = []
    for s in ("ctl", "trt"):
        for r in ("R1", "R2"):
            path = fq_dir / f"{s}_{r}.fastq.gz"
            path.write_bytes(b"")
            fastqs.append(str(path))
    bam_dir = tmp_path / "aligned"
    bam_dir.mkdir()
    env = _FakeRawFastqBulkAtacEnv(bam_dir)
    manifest = {"samples": [{"sample_id": "ctl"}, {"sample_id": "trt"}]}
    result = _agent(env)._run_bulk_atac(
        "exp",
        {"genome": "hg38", "genome_fasta": str(tmp_path / "genome.fa"),
         "bulk_fastq_manifest": manifest, "sample_ids": ["ctl", "trt"],
         "comparisons": [["treated", "control"]]},
        {},
        fastqs,
    )
    return env, result, fastqs, manifest


def test_raw_fastq_preprocessing_params_and_production_stamp(tmp_path):
    env, result, fastqs, manifest = _raw_fastq_bulk_findings(tmp_path)
    assert result["status"] == "done"
    qc_params = env.calls[0][2]
    assert env.calls[0][:2] == ("rnaseq", "rna_fastq_qc.py")
    assert qc_params["fastq_files"] == fastqs
    assert qc_params["fastq_dir"] == str(Path(fastqs[0]).parent)
    assert qc_params["bulk_fastq_manifest"] == manifest
    assert qc_params["sample_ids"] == ["ctl", "trt"]
    pre = result["findings"]["preprocessing"]
    assert pre["data_type"] == "bulk_ATAC"
    assert pre["validation_level"] == "production"
    assert pre["read_layout"] == "paired"
    assert pre["fastq_files"] == fastqs


def test_raw_fastq_ledger_node_and_strict_wclaim_report(tmp_path):
    from aria.agents.narrative.evidence_verifier import (
        verify_block_claim_support,
    )
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.run_ledger import build_run_ledger
    from aria.agents.narrative.validators import validate_blocks

    _env, result, fastqs, _manifest = _raw_fastq_bulk_findings(tmp_path)
    agent_results = {"chromatin_agent": {
        "status": "done", "findings": {"bulk_ATAC": result}}}

    ledger = build_run_ledger({}, agent_results)
    by_node = {e["node_id"]: e for e in ledger["entries"]}
    node = by_node["ledger://chromatin/bulk_atac_preprocessing"]
    assert node["status"] == "executed"
    assert node["divergence"] is False
    prov = node["provenance"]
    assert prov["script"] == "aria/scripts/atac_align.py"
    assert prov["environment_stack"] == "atacseq"
    assert prov["inputs"] == fastqs
    assert [Path(p).name for p in prov["artifacts"]] == [
        "ctl.filtered.bam", "trt.filtered.bam"]
    assert prov["validation_level"] == "production"
    # The scATAC node is untouched (not run) on a bulk lane.
    assert by_node["ledger://chromatin/scatac_preprocessing"]["status"] == \
        "not_run"

    blocks = validate_blocks(ChromatinNarrator().collect(
        "chromatin_agent", agent_results["chromatin_agent"], {}))
    by_id = {b.id: b for b in blocks}
    pre_block = by_id["chromatin.bulk_atac_preprocessing"]
    assert pre_block.analysis == "bulk_atac_preprocessing"
    assert pre_block.metadata["validation_level"] == "production"
    assert "2 sample(s) (paired layout)" in pre_block.claim
    # W-CLAIM: lowercase prose (no tool/format acronyms) after the capital.
    assert pre_block.claim[1:] == pre_block.claim[1:].lower()
    for block in blocks:
        manifest = verify_block_claim_support(block, strict=True)
        assert manifest["status"] == "supported", block.id
    assert (by_id["chromatin.peak_count_matrix.bulk_ATAC"]
            .metadata["validation_level"] == "production")
    assert (by_id["chromatin.differential_accessibility.bulk_atac"]
            .metadata["validation_level"] == "production")
    caveat_text = " ".join(
        c.text for c in by_id["chromatin.peak_count_matrix.bulk_ATAC"].caveats)
    assert "scaffolded" not in caveat_text
    assert "negative-binomial" in caveat_text
    assert "counts-per-million" in caveat_text

    pytest.importorskip("litellm")
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    agent.memory = type("M", (), {"db_path": ":memory:"})()
    report = agent._render_html_report(
        experiment_id="exp_bulk_atac_c3",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={"summary": "Bulk ATAC raw FASTQ differential accessibility"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={"high": [], "medium": [], "low": [],
                          "insufficient": []},
        methods="bulk ATAC bwa-mem2 + MACS3 + pyDESeq2",
        decisions=[],
        agent_results=agent_results,
        report_dir=tmp_path / "report",
    )
    assert report.exists()
    methodology = json.loads((report.parent / "methodology.json").read_text())
    assert methodology["claim_compilation"]["n_withheld"] == 0
    assert (methodology["run_ledger"]["claim_ledger_verification"]
            ["n_violations"] == 0)
    claim_analyses = {c["analysis"] for c in methodology["claims"]}
    assert "bulk_atac_preprocessing" in claim_analyses
    assert "differential_accessibility" in claim_analyses


def test_failed_alignment_is_not_an_executed_ledger_node():
    from aria.agents.narrative.run_ledger import build_run_ledger

    findings = {"bulk_ATAC": {"status": "failed", "findings": {
        "preprocessing": {"status": "failed", "ran": False,
                          "data_type": "bulk_ATAC", "stage": "align",
                          "analysis": "bulk_atac_alignment"}}}}
    ledger = build_run_ledger(
        {"design_intelligence": {"recommended": ["bwa-mem2 alignment"]}},
        {"chromatin_agent": {"status": "done", "findings": findings}})
    node = next(e for e in ledger["entries"]
                if e["node_id"] == "ledger://chromatin/bulk_atac_preprocessing")
    assert node["status"] == "error"
    assert node["divergence"] is True
    assert node["provenance"]["artifacts"] == []


def test_narrator_tables_surface_counts_cpm_and_per_comparison_csv():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    findings = {"bulk_ATAC": {"status": "done", "findings": {
        "peak_counts": {"status": "success", "data_type": "bulk_ATAC",
                        "counts_matrix_path": "counts.tsv",
                        "counts_cpm_path": "cpm.tsv",
                        "n_peaks": 3, "n_samples": 2},
        "differential_accessibility": {
            "status": "success", "ran": True, "data_type": "bulk_ATAC",
            "output_csv": "summary.csv",
            "comparisons": [
                {"test": "t", "reference": "c", "status": "success",
                 "full_results_csv": "t_vs_c_full.csv"},
                {"test": "u", "reference": "c", "status": "success",
                 "output_csv": "u_vs_c.csv"},
            ]},
    }}}
    tables = ChromatinNarrator().tables(
        "chromatin_agent", {"status": "done", "findings": findings})
    paths = {t["path"] for t in tables}
    assert {"counts.tsv", "cpm.tsv", "summary.csv", "t_vs_c_full.csv",
            "u_vs_c.csv"} <= paths
    block = ChromatinNarrator()._peak_counts_block(
        findings["bulk_ATAC"]["findings"]["peak_counts"])
    assert any(e.path == "cpm.tsv" for e in block.evidence)


def _write_raw_counts(path, n_peaks=60):
    import numpy as np

    rng = np.random.default_rng(7)
    samples = ["ctl_r1", "ctl_r2", "ctl_r3", "trt_r1", "trt_r2", "trt_r3"]
    lines = ["\t".join(["peak_id", *samples])]
    for i in range(n_peaks):
        base = rng.integers(20, 200)
        shift = 4 if i < 10 else 1
        vals = [int(rng.poisson(base)) for _ in range(3)] + [
            int(rng.poisson(base * shift)) for _ in range(3)]
        lines.append("\t".join([f"chr1:{i * 500}-{i * 500 + 300}",
                                *map(str, vals)]))
    path.write_text("\n".join(lines) + "\n")
    return samples


def test_chromatin_figure_bulk_writes_dual_png_svg(tmp_path):
    from aria.scripts.chromatin_figure_bulk import chromatin_figure_bulk

    counts = tmp_path / "bulk_atac_peak_counts.tsv"
    samples = _write_raw_counts(counts)
    out = chromatin_figure_bulk({
        "output_dir": str(tmp_path / "figs"),
        "counts_path": str(counts),
        "sample_metadata": {s: {"condition": s.split("_")[0]}
                            for s in samples},
        "qc": {"frip_per_sample": {s: 0.25 + 0.01 * i
                                   for i, s in enumerate(samples)},
               "tss_enrichment_per_sample": {s: 6.0 for s in samples},
               "dup_rate": 0.11},
    })
    assert out["status"] == "success"
    assert out["validation_level"] == "production"
    assert set(out["figures"]) == {"bulk_atac_pca", "bulk_atac_qc_summary"}
    for key, png in out["figures"].items():
        png = Path(png)
        assert png.suffix == ".png" and png.stat().st_size > 0, key
        svg = png.with_suffix(".svg")
        assert svg.exists() and svg.stat().st_size > 0, key
        assert svg.read_text()[:200].lstrip().startswith("<?xml")


def test_chromatin_figure_bulk_omits_degenerate_panels(tmp_path):
    from aria.scripts.chromatin_figure_bulk import chromatin_figure_bulk

    one = tmp_path / "one.tsv"
    one.write_text("peak_id\ts1\nchr1:1-2\t5\nchr1:3-4\t6\n")
    out = chromatin_figure_bulk({
        "output_dir": str(tmp_path / "figs"), "counts_path": str(one),
        "qc": {"frip": None, "tss_enrichment": None, "dup_rate": None},
    })
    assert out["status"] == "success"
    assert out["figures"] == {}
    assert set(out["skipped"]) == {"bulk_atac_pca", "bulk_atac_qc_summary"}
    assert not list((tmp_path / "figs").glob("*.png"))


def test_generate_figures_merges_bulk_panels_in_process_and_via_env(tmp_path):
    from aria.agents import _narrative_chromatin as nc

    counts = tmp_path / "bulk_atac_peak_counts.tsv"
    samples = _write_raw_counts(counts)
    meta = tmp_path / "samples.tsv"
    meta.write_text("sample_id\tcondition\n" + "".join(
        f"{s}\t{s.split('_')[0]}\n" for s in samples))
    findings = {
        "qc": {"status": "success", "data_type": "bulk_ATAC",
               "frip": 0.3, "tss_enrichment": 8.0, "dup_rate": 0.1},
        "peaks": {"status": "success", "data_type": "bulk_ATAC",
                  "frip": 0.3, "n_peaks": 60},
        "peak_counts": {"status": "success", "data_type": "bulk_ATAC",
                        "counts_matrix_path": str(counts),
                        "sample_metadata_path": str(meta)},
    }
    nc.generate_figures(findings, None, tmp_path / "figs_inproc",
                        env_manager=None)
    assert {"bulk_atac_pca", "bulk_atac_qc_summary"} <= set(findings["figures"])

    class _Env:
        def __init__(self):
            self.calls = []

        def run_in_stack(self, *, stack, script_path, params):
            self.calls.append((stack, script_path, params))
            if script_path.endswith("chromatin_figure_bulk.py"):
                from aria.scripts.chromatin_figure_bulk import (
                    render_bulk_figures,
                )
                return render_bulk_figures(params)
            return {"status": "error", "error_type": "Unexpected"}

    env = _Env()
    findings2 = {k: dict(v) for k, v in findings.items() if k != "figures"}
    nc.generate_figures(findings2, None, tmp_path / "figs_env",
                        env_manager=env)
    call = next(c for c in env.calls
                if c[1] == "aria/scripts/chromatin_figure_bulk.py")
    assert call[0] == "rna"
    assert call[2]["counts_path"] == str(counts)
    assert call[2]["sample_metadata"] == str(meta)
    assert {"bulk_atac_pca", "bulk_atac_qc_summary"} <= set(
        findings2["figures"])

    # scATAC findings never trigger the bulk figure script.
    env2 = _Env()
    nc.generate_figures({"qc": {"status": "success", "data_type": "scATAC"}},
                        None, tmp_path / "figs_sc", env_manager=env2)
    assert not any(c[1].endswith("chromatin_figure_bulk.py")
                   for c in env2.calls)
