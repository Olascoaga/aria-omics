"""F6 (preprint audit 2026-06-19): every committed bulk ATAC benchmark artifact must
carry reproducibility provenance.

The two V47 artifacts shipped with NO provenance (no git_commit/aria_version/tool
versions/seed/run date). This guard fails on any JSON in
``docs/benchmark_results/bulk_atac/`` that lacks the required provenance keys, so the
regression cannot return. It also pins the stamper's two non-negotiable behaviours:
the scientific body is preserved verbatim, and stamping is idempotent.
"""
import datetime as dt
import json
from pathlib import Path

from aria.benchmarks.bulk_atac_v47_provenance import (
    SOURCE_COMMITS,
    stamp_artifacts,
    stamp_body,
)

ROOT = Path(__file__).resolve().parents[1]
BENCH_DIR = ROOT / "docs" / "benchmark_results" / "bulk_atac"


def test_every_bulk_atac_json_has_provenance():
    """The committed surface: each artifact has aria_version + a provenance block
    with git_commit, source_commit, tool_versions, and a run date."""
    json_files = sorted(BENCH_DIR.rglob("*.json"))
    assert json_files, "no bulk ATAC benchmark artifacts found"
    for json_file in json_files:
        data = json.loads(json_file.read_text(encoding="utf-8"))
        rel = json_file.relative_to(ROOT)
        assert "aria_version" in data, f"{rel}: missing aria_version"
        prov = data.get("provenance")
        assert isinstance(prov, dict), f"{rel}: missing provenance block"
        for key in ("git_commit", "aria_version", "source_commit",
                    "tool_versions", "generated_utc"):
            assert prov.get(key), f"{rel}: provenance missing {key}"


def test_stamping_preserves_scientific_body():
    """Stamping adds only the owned provenance keys; the measured science is untouched."""
    sample = {"benchmark": "x", "result": {"n_significant": 68722}, "markers": [1, 2, 3]}
    out = stamp_body(dict(sample), "v47_k562_gm12878_da_validation.json", repo_root=ROOT)
    for key, value in sample.items():
        assert out[key] == value
    assert out["provenance"]["source_commit"] == SOURCE_COMMITS[
        "v47_k562_gm12878_da_validation.json"]
    assert out["provenance"]["tool_versions"]["pydeseq2"] == "0.5.4"
    # honest: deterministic pipeline -> no fabricated seed
    assert out["provenance"]["seed"] is None


def test_stamping_is_idempotent():
    """Re-stamping an already-stamped body does not nest or duplicate provenance."""
    fixed = dt.datetime(2026, 6, 20, tzinfo=dt.timezone.utc)
    body = {"benchmark": "x"}
    once = stamp_body(dict(body), "v47_k562_gm12878_da_validation.json",
                      repo_root=ROOT, now=fixed)
    twice = stamp_body(dict(once), "v47_k562_gm12878_da_validation.json",
                       repo_root=ROOT, now=fixed)
    assert once == twice
    assert "provenance" not in once["provenance"]


def test_stamp_artifacts_round_trips_real_files(tmp_path):
    """The runner reads the real committed bodies, stamps a copy, and reports them all."""
    work = tmp_path / "bulk_atac"
    work.mkdir()
    for name in SOURCE_COMMITS:
        src = BENCH_DIR / name
        if src.exists():
            (work / name).write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    manifest = stamp_artifacts(bench_dir=work, repo_root=ROOT)
    assert manifest["status"] in {"ok", "partial"}
    for entry in manifest["stamped"]:
        data = json.loads((work / entry["artifact"]).read_text(encoding="utf-8"))
        assert data["provenance"]["git_commit"]
        assert data["provenance"]["source_commit"] == entry["source_commit"]
