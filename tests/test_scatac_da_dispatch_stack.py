"""scATAC differential accessibility must dispatch to the rna stack.

`chromatin_diffacc.py` runs per-cluster wilcoxon (scanpy) + pseudobulk DESeq2
(pydeseq2). pydeseq2 0.5.4 requires numpy>=2, but aria-chromatin-env is pinned
numpy<2 (snapatac2/episcanpy/muon). aria-rna-env already ships scanpy + anndata +
pydeseq2 on numpy 2, so DA runs there with the same pydeseq2 the rest of ARIA's
DE/DA core uses. QC / LSI / motif / regulatory stay in the chromatin stack.
"""

from pathlib import Path


class _FakeScatacEnv:
    def __init__(self):
        self.calls = []

    def run_in_stack(self, *, stack, script_path, params):
        name = Path(script_path).name
        self.calls.append((stack, name, params))
        if name == "chromatin_qc.py":
            return {"status": "success", "data_type": "scATAC"}
        if name == "chromatin_lsi_clustering.py":
            return {"status": "success",
                    "output_path": str(Path(params["output_dir"]) / "clustered.h5ad")}
        if name == "chromatin_diffacc.py":
            return {"status": "success",
                    "per_cluster": {"output_csv": "/tmp/da.csv"}}
        if name in ("chromatin_motifs.py", "chromatin_regulatory.py"):
            return {"status": "success"}
        raise AssertionError(script_path)


def _agent(env):
    from aria.agents.chromatin_agent import ChromatinAgent

    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = env
    agent.publish_status = lambda *a, **k: None
    agent.publish_finding = lambda *a, **k: None
    agent.publish_escalation = lambda *a, **k: None
    agent._publish_qc_finding = lambda *a, **k: None
    return agent


def test_scatac_da_dispatches_to_rna_stack(tmp_path):
    env = _FakeScatacEnv()
    agent = _agent(env)
    h5mu = tmp_path / "sample.h5mu"
    h5mu.write_bytes(b"")

    result = agent._run_scatac_matrix(
        "exp",
        {
            "genome": "hg38",
            # genome_fasta present -> the motif-genome escalation branch is skipped.
            "genome_fasta": str(tmp_path / "genome.fa"),
            "condition_col": "condition",
            "replicate_col": "replicate",
            "comparisons": [["old", "young"]],
        },
        {},
        str(h5mu),
    )

    assert result["status"] == "done"
    by_script = {name: stack for stack, name, _ in env.calls}
    assert by_script["chromatin_diffacc.py"] == "rna"
    # the rest of the scATAC pipeline stays in the chromatin stack
    assert by_script["chromatin_qc.py"] == "chromatin"
    assert by_script["chromatin_lsi_clustering.py"] == "chromatin"
    assert by_script["chromatin_motifs.py"] == "chromatin"
    assert by_script["chromatin_regulatory.py"] == "chromatin"
    # confirmed cross-condition design is forwarded to the DA script
    da_params = next(p for s, n, p in env.calls if n == "chromatin_diffacc.py")
    assert da_params["comparisons"] == [["old", "young"]]
