"""S2 (pre-integration audit): atomic artifact writes in aria/scripts/_base.py.

Everything downstream assumes resume-by-file-validity: an artifact on disk is
complete. A plain json.dump(f) breaks that — a crash mid-serialization leaves a
truncated file that corrupts resume. _atomic_write_json must guarantee the final
artifact is always either absent/previous or fully written, never partial.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from aria.scripts._base import _atomic_write_json, _json_serializer


def test_atomic_write_produces_valid_complete_json(tmp_path):
    out = tmp_path / "result.json"
    payload = {"status": "success", "n_peaks": 12345, "vals": [1, 2, 3]}
    _atomic_write_json(out, payload)
    assert json.loads(out.read_text()) == payload


def test_serialization_failure_leaves_no_partial_file(tmp_path):
    out = tmp_path / "result.json"
    # An un-serializable object (no default handler) must raise BEFORE any file
    # is created — not leave a truncated artifact.
    with pytest.raises(TypeError):
        _atomic_write_json(out, {"bad": object()})
    assert not out.exists()
    # And no temp leftovers.
    assert list(tmp_path.glob("*.tmp")) == []


def test_failed_write_preserves_existing_artifact(tmp_path):
    """The resume-critical guarantee: a failed re-write keeps the previous valid
    artifact intact (a plain json.dump would have truncated it)."""
    out = tmp_path / "result.json"
    good = {"status": "success", "generation": 1}
    _atomic_write_json(out, good)

    # Re-write with an un-serializable payload -> must fail without touching out.
    with pytest.raises(TypeError):
        _atomic_write_json(out, {"status": "success", "bad": object()})

    assert json.loads(out.read_text()) == good  # original intact
    assert list(tmp_path.glob("*.tmp")) == []   # no leftovers


def test_replace_failure_cleans_tmp_and_raises(tmp_path, monkeypatch):
    """If the atomic os.replace itself fails, the temp file is cleaned up and the
    error propagates (no silent partial state)."""
    out = tmp_path / "result.json"

    def boom(src, dst):
        raise OSError("simulated replace failure")

    monkeypatch.setattr(os, "replace", boom)
    with pytest.raises(OSError):
        _atomic_write_json(out, {"status": "success"})
    assert not out.exists()
    assert list(tmp_path.glob("*.tmp")) == []


def test_uses_json_serializer_for_numpy(tmp_path):
    np = pytest.importorskip("numpy")
    out = tmp_path / "result.json"
    _atomic_write_json(out, {"x": np.int64(7), "y": np.float64(1.5)},
                       default=_json_serializer)
    data = json.loads(out.read_text())
    assert data == {"x": 7, "y": 1.5}


def test_creates_parent_dir(tmp_path):
    out = tmp_path / "nested" / "deep" / "result.json"
    _atomic_write_json(out, {"status": "success"})
    assert out.exists()


def test_run_script_writes_atomically(tmp_path, monkeypatch):
    """End-to-end: run_script routes its success write through the atomic helper
    (valid JSON, no temp leftovers)."""
    import sys
    from aria.scripts import _base

    inp = tmp_path / "in.json"
    outp = tmp_path / "out.json"
    inp.write_text(json.dumps({"k": "v"}))
    monkeypatch.setattr(sys, "argv", ["script.py", str(inp), str(outp)])

    _base.run_script(lambda params: {"echo": params["k"]})

    data = json.loads(outp.read_text())
    assert data == {"echo": "v", "status": "success"}
    assert list(tmp_path.glob("*.tmp")) == []
