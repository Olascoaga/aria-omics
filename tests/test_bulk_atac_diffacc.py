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
