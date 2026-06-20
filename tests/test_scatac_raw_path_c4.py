"""C4 Codex audit guards for the raw scATAC fragments/BAM path.

Raw scATAC inputs do not yet produce the cell x peak matrix required by the
validated matrix pipeline. The live agent must disclose that limitation instead
of computing dead LSI params or pretending LSI/DA/regulatory layers ran.
"""

from pathlib import Path


class _RawScatacEnv:
    def __init__(self):
        self.calls = []

    def run_in_stack(self, *, stack, script_path, params):
        name = Path(script_path).name
        self.calls.append((stack, name, params))
        if name == "chromatin_qc.py":
            return {"status": "success", "data_type": "scATAC", "warnings": []}
        if name == "chromatin_peaks.py":
            return {
                "status": "success",
                "data_type": "scATAC",
                "n_peaks": 10,
                "peaks_path": "/tmp/scatac_peaks.narrowPeak",
            }
        if name == "chromatin_motifs.py":
            return {"status": "success", "ran": True, "data_type": "scATAC"}
        raise AssertionError(f"unexpected raw scATAC dispatch: {name}")


def _agent(env):
    from aria.agents.chromatin_agent import ChromatinAgent

    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = env
    agent.publish_status = lambda *a, **k: None
    agent.publish_finding = lambda *a, **k: None
    agent.publish_escalation = lambda *a, **k: None
    agent._publish_qc_finding = lambda *a, **k: None
    agent._publish_peaks_finding = lambda *a, **k: None
    return agent


def test_raw_scatac_path_discloses_matrix_pipeline_skip_without_dead_lsi():
    env = _RawScatacEnv()
    agent = _agent(env)

    res = agent._run_scatac(
        "exp",
        {"genome": "hg38", "organism": "Homo sapiens"},
        {"user_question": "compare accessibility"},
        ["/tmp/fragments.tsv.gz"],
    )

    assert res["status"] == "done"
    findings = res["findings"]
    assert "lsi_params" not in findings
    assert "lsi" not in findings
    assert "differential_accessibility" not in findings
    assert "regulatory" not in findings
    matrix = findings["matrix_pipeline"]
    assert matrix["status"] == "skipped"
    assert matrix["ran"] is False
    assert matrix["reason"] == "raw_scatac_requires_peak_matrix_h5mu"
    assert matrix["required_input"] == ".h5mu peak matrix"
    assert set(matrix["skipped_steps"]) == {
        "lsi_clustering",
        "differential_accessibility",
        "regulatory_layers",
    }

    dispatched = {name for _stack, name, _params in env.calls}
    assert dispatched == {"chromatin_qc.py", "chromatin_peaks.py"}


def test_raw_scatac_matrix_skip_becomes_narrative_limitation():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    result = {
        "status": "done",
        "findings": {
            "matrix_pipeline": {
                "status": "skipped",
                "ran": False,
                "data_type": "scATAC",
                "analysis": "matrix_pipeline",
                "validation_level": "beta",
                "reason": "raw_scatac_requires_peak_matrix_h5mu",
                "message": "LSI/DA/regulatory layers require a processed .h5mu.",
                "required_input": ".h5mu peak matrix",
                "completed_steps": ["qc", "peak_calling"],
                "skipped_steps": [
                    "lsi_clustering",
                    "differential_accessibility",
                    "regulatory_layers",
                ],
            }
        },
    }

    blocks = ChromatinNarrator().collect("chromatin_agent", result, {})

    block = next(b for b in blocks
                 if b.id == "chromatin.scatac_matrix_pipeline.skipped")
    assert block.status == "limitation"
    assert block.confidence == "insufficient"
    assert block.metrics["reason"] == "raw_scatac_requires_peak_matrix_h5mu"
    assert any(".h5mu" in c.text for c in block.caveats)
