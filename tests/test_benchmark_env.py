from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCH_ENV = ROOT / "envs" / "aria-bench-env.yml"


def test_aria_bench_env_spec_exists_and_is_separate_from_runtime_envs():
    text = BENCH_ENV.read_text(encoding="utf-8")

    assert "name: aria-bench-env" in text
    assert "stack=\"benchmark\"" in text
    assert "aria-rna-env" in text
    assert "aria-chromatin-env" not in text


def test_aria_bench_env_declares_frozen_reference_comparator_stack():
    text = BENCH_ENV.read_text(encoding="utf-8").lower()

    for dep in (
        "r-base",
        "bioconductor-deseq2",
        "bioconductor-edger",
        "bioconductor-limma",
        "bioconductor-apeglm",
        "bioconductor-muscat",
        "bioconductor-mast",
        "bioconductor-speckle",
        "bioconductor-singler",
        "bioconductor-fgsea",
        "bioconductor-gsva",
    ):
        assert dep in text
