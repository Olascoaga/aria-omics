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
        truth = pd.read_csv(Path(cmd[cmd.index("--counts") + 1]), sep="\t")
        genes = truth["gene"].astype(str).tolist()
        table = pd.DataFrame({
            "gene": genes,
            "log2FoldChange": [idx / max(len(genes), 1) for idx, _ in enumerate(genes)],
            "pvalue": [1.0] * len(genes),
            "padj": [1.0] * len(genes),
        })
        methods = {}
        for key in ("deseq2", "edgeR_QLF", "limma_voom"):
            path = r_out / f"{key}.tsv"
            table.to_csv(path, sep="\t", index=False)
            methods[key] = {
                "status": "success",
                "method": key,
                "results_tsv": str(path),
                "n_tested": len(genes),
                "n_called": 0,
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
    assert Path(result["artifacts"]["manifest_json"]).exists()
    assert result["methods"]["deseq2"]["scores"]["significant_call_concordance"]["n_called"] == 0
