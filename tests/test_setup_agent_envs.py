from pathlib import Path


def test_setup_agent_uses_dedicated_ingestion_env_for_scrna_fastq():
    from aria.agents.setup_agent import ARIA_ENVS, SetupAgent

    agent = SetupAgent.__new__(SetupAgent)
    envs = agent._needed_envs(
        {"scRNA": ["sample_R1.fastq.gz", "sample_R2.fastq.gz"]},
        has_fastq=True,
    )

    assert "aria-ingestion-env" in envs
    assert "aria-rnaseq-env" not in envs
    assert "aria-rna-env" in envs
    assert ARIA_ENVS["aria-ingestion-env"]["stack"] == "ingestion"
    assert (Path("envs") / "aria-ingestion-env.yml").exists()


def test_setup_agent_keeps_bulk_fastq_on_rnaseq_env():
    from aria.agents.setup_agent import SetupAgent

    agent = SetupAgent.__new__(SetupAgent)
    envs = agent._needed_envs(
        {"bulk_RNA_raw": ["sample_R1.fastq.gz", "sample_R2.fastq.gz"]},
        has_fastq=True,
    )

    assert "aria-rnaseq-env" in envs
    assert "aria-ingestion-env" not in envs


def test_setup_agent_does_not_require_tobias_for_chip_or_cutandrun():
    from aria.agents.setup_agent import SetupAgent

    agent = SetupAgent.__new__(SetupAgent)
    for mod in ("ChIP", "CUT_AND_RUN", "CUT_AND_TAG"):
        envs = agent._needed_envs({mod: ["sample.bam"]}, has_fastq=False)
        assert "aria-chromatin-env" in envs
        assert "aria-rna-env" in envs
        assert "aria-tobias-env" not in envs

    atac_envs = agent._needed_envs({"bulk_ATAC": ["sample.bam"]}, has_fastq=False)
    assert "aria-tobias-env" in atac_envs
