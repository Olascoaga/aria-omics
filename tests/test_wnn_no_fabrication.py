"""P0-8 regression: integration_wnn must not fabricate scientific output.

The WNN script (a v4.7 scaffold, dispatch-gated by INTEGRATION_VALIDATION) shipped
three fabrications: `_load_atac` returned an empty `AnnData()` as the "peak
matrix", `_mock_wnn` returned hardcoded cluster counts / modality weights as a
"success", and a `0.6/0.4` weight fallback was invented when muon did not expose
weights. None of this is real science. Per ADR-002/ADR-025 (P0-8), the
unimplemented peak-matrix step is now an explicit structural blocker and the
fabricated success paths are gone — without building real WNN.
"""

import pytest

from aria.scripts import integration_wnn as wnn


def test_load_atac_is_an_explicit_blocker_not_a_fake_matrix(tmp_path):
    frag = tmp_path / "fragments.tsv.gz"
    frag.write_bytes(b"")
    # The real peak-matrix construction is not implemented: raise, never return
    # a fabricated empty/placeholder AnnData.
    with pytest.raises(NotImplementedError):
        wnn._load_atac([str(frag)], "hg38", "Homo sapiens")


def test_mock_wnn_fabrication_helper_is_removed():
    assert not hasattr(wnn, "_mock_wnn"), \
        "_mock_wnn fabricates a success result and must not exist (ADR-002)"


def test_integration_wnn_never_fabricates_a_success_under_mock_flag(tmp_path):
    rna = tmp_path / "rna.h5ad"
    rna.write_bytes(b"")
    frag = tmp_path / "fragments.tsv.gz"
    frag.write_bytes(b"")
    result = wnn.integration_wnn({
        "rna_files": [str(rna)],
        "atac_files": [str(frag)],
        "mock": True,          # legacy mock opt-in must no longer fabricate
        "output_dir": str(tmp_path / "out"),
    })
    assert result["status"] == "error"
    assert result["error_type"] in {"MissingDependency", "NotImplemented"}
    # No fabricated modality weights / cluster counts leak through.
    assert "mean_rna_weight" not in result
    assert "n_joint_clusters" not in result


def test_module_has_no_hardcoded_fabricated_weights():
    import inspect
    src = inspect.getsource(wnn)
    for needle in ("0.62", "0.38", "= 0.6\n", "= 0.4\n"):
        assert needle not in src, f"fabricated weight literal {needle!r} still present"
