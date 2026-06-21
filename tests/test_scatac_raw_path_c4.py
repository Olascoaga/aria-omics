"""C4 Codex audit guards for the raw scATAC fragments/BAM path.

The B2b fragments->matrix bridge lets raw scATAC fragments reach the full matrix
pipeline (LSI/DA/regulatory). When the bridge CAN build a cell x peak matrix, the
agent dispatches the full pipeline; when it CANNOT (snapatac2/genome unavailable,
unreadable fragments), it must honestly disclose the skip instead of computing
dead LSI params or pretending the layers ran.
"""

from pathlib import Path


class _RawScatacEnv:
    """Raw fragments path where the fragments->matrix bridge cannot build a
    matrix (e.g. snapatac2 unavailable) -> honest matrix_pipeline skip."""

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
        if name == "chromatin_fragments_to_matrix.py":
            return {"status": "skipped", "ran": False,
                    "analysis": "fragments_to_peak_matrix",
                    "reason": "snapatac2_unavailable"}
        if name == "chromatin_motifs.py":
            return {"status": "success", "ran": True, "data_type": "scATAC"}
        raise AssertionError(f"unexpected raw scATAC dispatch: {name}")


class _BridgeOkEnv:
    """Raw fragments path where the bridge builds a matrix -> full pipeline."""

    def __init__(self):
        self.calls = []

    def run_in_stack(self, *, stack, script_path, params):
        name = Path(script_path).name
        self.calls.append((stack, name, params))
        if name == "chromatin_qc.py":
            return {"status": "success", "data_type": "scATAC", "warnings": []}
        if name == "chromatin_peaks.py":
            return {"status": "success", "data_type": "scATAC", "n_peaks": 10,
                    "peaks_path": "/tmp/scatac_peaks.narrowPeak"}
        if name == "chromatin_fragments_to_matrix.py":
            return {"status": "success", "ran": True,
                    "output_path": "/tmp/aria_chromatin/fragments_peak_matrix.h5ad",
                    "n_cells": 1200, "n_peaks": 50000, "peak_mode": "provided"}
        if name == "chromatin_lsi_clustering.py":
            return {"status": "success",
                    "output_path": "/tmp/aria_chromatin/lsi_clustered.h5ad"}
        if name == "chromatin_diffacc.py":
            return {"status": "success", "per_cluster": {"output_csv": "/tmp/da.csv"}}
        if name in ("chromatin_motifs.py", "chromatin_regulatory.py"):
            return {"status": "success"}
        if name == "chromatin_footprint_tobias.py":
            return {"status": "success", "ran": False, "reason": "no_tobias"}
        raise AssertionError(f"unexpected dispatch: {name}")


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
    # The skip reason now propagates from the bridge attempt (honest cause).
    assert matrix["reason"] == "snapatac2_unavailable"
    assert set(matrix["skipped_steps"]) == {
        "lsi_clustering",
        "differential_accessibility",
        "regulatory_layers",
    }
    # The bridge was attempted and recorded.
    assert findings["fragments_to_matrix"]["status"] == "skipped"
    dispatched = {name for _stack, name, _params in env.calls}
    assert dispatched == {"chromatin_qc.py", "chromatin_peaks.py",
                          "chromatin_fragments_to_matrix.py"}


def test_raw_scatac_bridge_success_runs_full_matrix_pipeline():
    env = _BridgeOkEnv()
    agent = _agent(env)

    res = agent._run_scatac(
        "exp",
        {"genome": "hg38", "organism": "Homo sapiens",
         "condition_col": "condition", "replicate_col": "replicate",
         "comparisons": [["old", "young"]]},
        {"user_question": "compare accessibility"},
        ["/tmp/fragments.tsv.gz"],
    )

    assert res["status"] == "done"
    findings = res["findings"]
    # The bridge built the matrix and the full pipeline ran from raw fragments.
    assert findings["fragments_to_matrix"]["status"] == "success"
    assert findings["lsi"]["status"] == "success"
    assert findings["differential_accessibility"]["status"] == "success"
    assert "matrix_pipeline" not in findings  # no honest-skip placeholder
    dispatched = {name for _stack, name, _params in env.calls}
    assert "chromatin_lsi_clustering.py" in dispatched
    assert "chromatin_diffacc.py" in dispatched
    # The fragments-side QC was NOT recomputed on the derived matrix.
    assert sum(1 for _s, n, _p in env.calls if n == "chromatin_qc.py") == 1


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
