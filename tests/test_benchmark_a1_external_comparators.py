import json
from pathlib import Path


def test_a1_external_comparator_contract_registered():
    from aria.utils.script_contracts import SCRIPT_CONTRACTS

    key = "aria/scripts/benchmark_a1_external_comparators.py"
    assert key in SCRIPT_CONTRACTS
    assert SCRIPT_CONTRACTS[key].validation_level == "beta"


def test_a1_external_r_runner_uses_portable_base_io():
    script = Path("aria/scripts/benchmark_a1_external_comparators.R").read_text()

    assert 'dest = "output_dir"' in script
    assert 'dest = "output_json"' in script
    assert "library(data.table)" not in script
    assert "read.delim" in script
    assert "write.table" in script


def test_a1_external_comparator_runner_reports_missing_rscript(tmp_path):
    from aria.scripts.benchmark_a1_external_comparators import (
        benchmark_a1_external_comparators,
    )

    result = benchmark_a1_external_comparators({
        "output_dir": str(tmp_path),
        "quick": True,
        "rscript_bin": "definitely-not-rscript",
    })

    assert result["status"] == "error"
    assert result["error_type"] == "MissingDependency"
    assert Path(result["artifacts"]["counts_tsv"]).exists()
    assert Path(result["artifacts"]["truth_tsv"]).exists()


def test_a1_external_comparator_runner_scores_mocked_r_output(tmp_path, monkeypatch):
    import pandas as pd
    from aria.scripts import benchmark_a1_external_comparators as mod

    def fake_which(_name):
        return "/usr/bin/Rscript"

    class Result:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run(cmd, *args, **kwargs):
        r_out = Path(cmd[cmd.index("--output-dir") + 1])
        r_json = Path(cmd[cmd.index("--output-json") + 1])
        counts_path = Path(cmd[cmd.index("--counts") + 1])
        truth = pd.read_csv(counts_path.parent / "a1_truth.tsv", sep="\t")
        genes = truth["gene"].astype(str).tolist()
        is_de = truth["is_true_de"].astype(bool).tolist()
        true_lfc = truth["true_log2fc"].astype(float).tolist()
        # Plant known-truth R comparator table: exact true_log2fc, padj=0.001 on DE, 0.8 on null.
        table = pd.DataFrame({
            "gene": genes,
            "log2FoldChange": true_lfc,
            "pvalue": [1e-5 if d else 0.7 for d in is_de],
            "padj": [1e-3 if d else 0.8 for d in is_de],
        })
        n_called_true = sum(1 for d, lfc in zip(is_de, true_lfc) if d and abs(lfc) >= 0.5)
        methods = {}
        for key in ("deseq2", "edgeR_QLF", "limma_voom"):
            path = r_out / f"{key}.tsv"
            table.to_csv(path, sep="\t", index=False)
            methods[key] = {
                "status": "success",
                "method": key,
                "results_tsv": str(path),
                "n_tested": len(genes),
                "n_called": n_called_true,
            }
        r_json.write_text(json.dumps({
            "status": "success",
            "methods": methods,
            "versions": {"R": "mock"},
        }), encoding="utf-8")
        return Result()

    monkeypatch.setattr(mod.shutil, "which", fake_which)
    monkeypatch.setattr(mod.subprocess, "run", fake_run)

    result = mod.benchmark_a1_external_comparators({
        "output_dir": str(tmp_path),
        "quick": True,
    })

    assert result["status"] == "success"
    assert set(result["methods"]) == {"deseq2", "edgeR_QLF", "limma_voom"}
    # T10 (tri-audit 2026-06-14): the A1 relativizer keeps machine-absolute paths
    # out of the public manifest, so artifacts["manifest_json"] is RELATIVE (a
    # basename when output_dir is outside the repo, as a tmp_path is). The earlier
    # `Path(manifest_json).exists()` wrongly assumed it resolved from the CWD.
    # Verify the public path is non-absolute and the file exists where it was
    # actually written (under output_dir). This is the A1 generator E2E exercise.
    manifest_json = result["artifacts"]["manifest_json"]
    assert not manifest_json.startswith("/")
    assert (tmp_path / Path(manifest_json).name).exists()
    # D3-021: assert exact analytical scoring metrics on known-truth R comparator table.
    scores = result["methods"]["deseq2"]["scores"]
    sig_sc = scores["significant_call_concordance"]
    lfc_sc = scores["lfc_concordance"]
    rank_sc = scores["ranking_concordance"]
    assert sig_sc["n_called"] == sig_sc["n_true_de"]
    assert sig_sc["recall"] == 1.0
    assert sig_sc["precision"] == 1.0
    assert sig_sc["empirical_fdr"] == 0.0
    assert lfc_sc["pearson_true_de"] > 0.99
    assert lfc_sc["spearman_true_de"] > 0.99
    assert rank_sc["top_k_jaccard"] == 1.0
