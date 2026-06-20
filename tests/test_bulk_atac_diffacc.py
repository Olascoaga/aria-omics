"""V47 bulk ATAC differential accessibility guards."""

from pathlib import Path


def test_bulk_atac_diffacc_contract_registered():
    from aria.utils.script_contracts import contract_for_script

    contract = contract_for_script("aria/scripts/chromatin_bulk_diffacc.py")
    assert contract is not None
    assert contract.validation_level == "beta"
    assert {f.name for f in contract.success_outputs} >= {
        "ran", "data_type", "analysis",
    }


def test_bulk_atac_diffacc_requires_explicit_comparison(tmp_path):
    from aria.scripts.chromatin_bulk_diffacc import chromatin_bulk_diffacc

    counts = tmp_path / "counts.tsv"
    counts.write_text("peak_id\tA\tB\nchr1:1-10\t1\t2\n")
    meta = tmp_path / "samples.tsv"
    meta.write_text(
        "sample_id\tcondition\treplicate\n"
        "A\tctrl\tr1\nB\ttreated\tr2\n"
    )

    res = chromatin_bulk_diffacc({
        "data_type": "bulk_ATAC",
        "counts_matrix_path": str(counts),
        "sample_metadata_path": str(meta),
    })

    assert res["status"] == "success"
    assert res["ran"] is False
    assert "no explicit comparison" in res["reason"]


def test_prepare_replicate_matrix_discloses_unusable_covariates():
    import pandas as pd
    from aria.scripts.chromatin_bulk_diffacc import _prepare_replicate_matrix

    counts = pd.DataFrame({
        "c1a": [1, 2], "c1b": [3, 4], "c2": [5, 6],
        "t1": [7, 8], "t2": [9, 10],
    }, index=["chr1:1-10", "chr1:20-30"])
    metadata = pd.DataFrame({
        "condition": ["ctrl", "ctrl", "ctrl", "treated", "treated"],
        "replicate": ["r1", "r1", "r2", "r1", "r2"],
        "batch": ["b1", "b2", "b1", "b1", "b1"],
        "lane": ["L1", "L1", "L2", "L1", "L2"],
    }, index=["c1a", "c1b", "c2", "t1", "t2"])

    rep_counts, rep_meta, usable, dropped = _prepare_replicate_matrix(
        counts, metadata, "condition", "replicate",
        ["batch", "lane", "missing"],
    )

    assert list(rep_counts.columns) == [
        "ctrl::r1", "ctrl::r2", "treated::r1", "treated::r2"]
    assert rep_counts.loc["chr1:1-10", "ctrl::r1"] == 4
    assert rep_meta.loc["ctrl::r1", "lane"] == "L1"
    assert usable == ["lane"]
    by_cov = {d["covariate"]: d for d in dropped}
    assert by_cov["missing"]["reason"] == "not present in the sample metadata"
    assert by_cov["batch"]["affected_replicates"] == ["ctrl::r1"]
    assert "not constant" in by_cov["batch"]["reason"]


def test_bulk_atac_diffacc_runs_shared_deseq2_and_gates_convergence(
        tmp_path, monkeypatch):
    import pandas as pd

    from aria.scripts.chromatin_bulk_diffacc import chromatin_bulk_diffacc
    import aria.scripts.rna_bulk_de as rna_bulk_de

    counts = tmp_path / "counts.tsv"
    counts.write_text(
        "peak_id\tc1\tc2\tc3\tt1\tt2\tt3\n"
        "chr1:1-10\t10\t11\t12\t80\t81\t82\n"
        "chr1:20-30\t30\t29\t31\t5\t4\t6\n"
        "chr2:1-9\t7\t7\t7\t7\t7\t7\n"
    )
    meta = tmp_path / "samples.tsv"
    meta.write_text(
        "sample_id\tcondition\treplicate\n"
        "c1\tctrl\tr1\nc2\tctrl\tr2\nc3\tctrl\tr3\n"
        "t1\ttreated\tr1\nt2\ttreated\tr2\nt3\ttreated\tr3\n"
    )
    captured = {}

    def fake_deseq2(count_df, meta_df, design_factor, numerator, denominator,
                    padj_thr, lfc_thr, **kwargs):
        captured["count_columns"] = list(count_df.columns)
        captured["metadata_index"] = list(meta_df.index)
        captured["design_factor"] = design_factor
        captured["comparison"] = (numerator, denominator)
        captured["expose_convergence"] = kwargs.get("expose_convergence")
        results = pd.DataFrame({
            "log2FoldChange": [2.1, -1.7, 0.2],
            "padj": [0.001, 0.01, 0.9],
            "lfc_converged": [True, False, True],
            "baseMean": [40.0, 20.0, 7.0],
        }, index=["chr1:1-10", "chr1:20-30", "chr2:1-9"])
        return {
            "status": "success",
            "results": results,
            "n_replicates": {"test": 3, "ref": 3},
            "low_power_warning": False,
            "lfc_shrinkage": {"requested": True, "applied": True},
            "lfc_threshold_test": {"applied": True},
            "fitted_design_formula": "~ condition",
            "power_estimate_at_lfc_min": {"power": 0.8},
        }, []

    monkeypatch.setattr(rna_bulk_de, "_run_deseq2", fake_deseq2)

    res = chromatin_bulk_diffacc({
        "data_type": "bulk_ATAC",
        "counts_matrix_path": str(counts),
        "sample_metadata_path": str(meta),
        "comparisons": [["treated", "ctrl"]],
        "output_dir": str(tmp_path / "da"),
    })

    assert res["status"] == "success"
    assert res["ran"] is True
    assert res["validation_level"] == "beta"
    assert captured["count_columns"] == [
        "ctrl::r1", "ctrl::r2", "ctrl::r3",
        "treated::r1", "treated::r2", "treated::r3",
    ]
    assert captured["design_factor"] == "condition"
    assert captured["comparison"] == ("treated", "ctrl")
    assert captured["expose_convergence"] is True
    comp = res["comparisons"][0]
    assert comp["status"] == "success"
    assert comp["n_sig"] == 1
    assert comp["n_up"] == 1
    assert comp["n_down"] == 0
    assert comp["n_nonconverged_excluded"] == 1
    assert Path(comp["full_results_csv"]).exists()
    assert Path(res["output_csv"]).exists()


def test_bulk_atac_diffacc_reports_dropped_covariates(tmp_path, monkeypatch):
    import pandas as pd

    from aria.scripts.chromatin_bulk_diffacc import chromatin_bulk_diffacc
    import aria.scripts.rna_bulk_de as rna_bulk_de

    counts = tmp_path / "counts.tsv"
    counts.write_text(
        "peak_id\tc1a\tc1b\tc2\tc3\tt1\tt2\tt3\n"
        "chr1:1-10\t10\t11\t12\t13\t80\t81\t82\n"
        "chr1:20-30\t30\t29\t31\t32\t5\t4\t6\n"
    )
    meta = tmp_path / "samples.tsv"
    meta.write_text(
        "sample_id\tcondition\treplicate\tbatch\tlane\n"
        "c1a\tctrl\tr1\tb1\tL1\n"
        "c1b\tctrl\tr1\tb2\tL1\n"
        "c2\tctrl\tr2\tb1\tL2\n"
        "c3\tctrl\tr3\tb1\tL3\n"
        "t1\ttreated\tr1\tb1\tL1\n"
        "t2\ttreated\tr2\tb1\tL2\n"
        "t3\ttreated\tr3\tb1\tL3\n"
    )
    captured = {}

    def fake_deseq2(count_df, meta_df, design_factor, numerator, denominator,
                    padj_thr, lfc_thr, **kwargs):
        captured["covariates"] = kwargs.get("covariates")
        results = pd.DataFrame({
            "log2FoldChange": [2.1, -1.7],
            "padj": [0.001, 0.8],
            "lfc_converged": [True, True],
            "baseMean": [40.0, 20.0],
        }, index=["chr1:1-10", "chr1:20-30"])
        return {
            "status": "success",
            "results": results,
            "n_replicates": {"test": 3, "ref": 3},
            "fitted_design_formula": "~ lane + condition",
            "covariates_adjusted": ["lane"],
            "covariates_dropped": [],
        }, []

    monkeypatch.setattr(rna_bulk_de, "_run_deseq2", fake_deseq2)

    res = chromatin_bulk_diffacc({
        "data_type": "bulk_ATAC",
        "counts_matrix_path": str(counts),
        "sample_metadata_path": str(meta),
        "comparisons": [["treated", "ctrl"]],
        "covariates": ["batch", "lane", "missing"],
        "output_dir": str(tmp_path / "da"),
    })

    assert res["status"] == "success" and res["ran"] is True
    assert captured["covariates"] == ["lane"]
    assert res["covariates_adjusted"] == ["lane"]
    dropped = {d["covariate"]: d for d in res["covariates_dropped"]}
    assert set(dropped) == {"batch", "missing"}
    assert dropped["batch"]["affected_replicates"] == ["ctrl::r1"]
    assert any("Covariate 'batch' was not adjusted" in w
               for w in res["warnings"])
    assert res["comparisons"][0]["covariates_dropped"]


def test_bulk_atac_diffacc_blocks_dropped_required_covariate(tmp_path):
    from aria.scripts.chromatin_bulk_diffacc import chromatin_bulk_diffacc

    counts = tmp_path / "counts.tsv"
    counts.write_text(
        "peak_id\tc1a\tc1b\tt1\tt2\n"
        "chr1:1-10\t10\t11\t80\t81\n"
    )
    meta = tmp_path / "samples.tsv"
    meta.write_text(
        "sample_id\tcondition\treplicate\tbatch\n"
        "c1a\tctrl\tr1\tb1\n"
        "c1b\tctrl\tr1\tb2\n"
        "t1\ttreated\tr1\tb1\n"
        "t2\ttreated\tr2\tb1\n"
    )

    res = chromatin_bulk_diffacc({
        "data_type": "bulk_ATAC",
        "counts_matrix_path": str(counts),
        "sample_metadata_path": str(meta),
        "comparisons": [["treated", "ctrl"]],
        "covariates": ["batch"],
        "required_covariates": ["batch"],
    })

    assert res["status"] == "success"
    assert res["ran"] is False
    assert "required covariate" in res["reason"]
    assert res["covariates_dropped"][0]["covariate"] == "batch"


def test_bulk_atac_diffacc_skips_when_metadata_lacks_replicates(tmp_path):
    from aria.scripts.chromatin_bulk_diffacc import chromatin_bulk_diffacc

    counts = tmp_path / "counts.tsv"
    counts.write_text("peak_id\tA\tB\nchr1:1-10\t1\t2\n")
    meta = tmp_path / "samples.tsv"
    meta.write_text("sample_id\tcondition\nA\tctrl\nB\ttreated\n")

    res = chromatin_bulk_diffacc({
        "data_type": "bulk_ATAC",
        "counts_matrix_path": str(counts),
        "sample_metadata_path": str(meta),
        "comparisons": [["treated", "ctrl"]],
    })

    assert res["status"] == "success"
    assert res["ran"] is False
    assert "replicate" in res["reason"]
