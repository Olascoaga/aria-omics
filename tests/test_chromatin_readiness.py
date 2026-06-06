"""Stage 3 v4.6 readiness gate (audit 2026-05-29):

- B9: chromatin_qc must not fabricate metrics. TSS enrichment, FRiP, and the
  per-assay FRiP constants are gone; barcode counts are real; not-computed
  metrics are reported as None, never as placeholders.
- C8: the `.h5mu` MuData entry path is detected by DataAudit and read by a real
  reader that blocks cleanly when the tooling is absent.
- B12: EnvironmentManager no longer reads a (never-written) env_aliases.json.
"""

import inspect
import os
from pathlib import Path

import pytest


# ── B9: no fabricated chromatin QC metrics ───────────────────────────────────

def test_tss_and_frip_are_not_fabricated(tmp_path):
    from aria.scripts.chromatin_qc import (
        _compute_tss_enrichment, _estimate_frip, _estimate_frip_bulk,
    )
    f = tmp_path / "fragments.tsv"
    f.write_text("chr1\t100\t250\tAAA\t1\n")
    # TSS used to be base + hash(frag_file) % 30; FRiP a constant 0.35; the bulk
    # FRiP a per-assay constant. All must now be not-computed.
    assert _compute_tss_enrichment(str(f), "hg38") is None
    assert _estimate_frip(str(f)) is None
    assert _estimate_frip_bulk("bulk_ATAC") is None
    assert _estimate_frip_bulk("ChIP") is None


def test_fragment_sizes_counts_real_barcodes(tmp_path):
    from aria.scripts.chromatin_qc import _compute_fragment_sizes
    lines = [
        "chr1\t100\t250\tAAA\t1",
        "chr1\t100\t260\tBBB\t1",
        "chr2\t100\t180\tAAA\t1",
        "chr2\t100\t900\tCCC\t1",
    ]
    f = tmp_path / "fragments.tsv"
    f.write_text("\n".join(lines) + "\n")
    out = _compute_fragment_sizes(str(f))
    assert out["n_barcodes"] == 3          # AAA, BBB, CCC — not the old len(set()) == 0
    assert out["n_fragments"] == 4
    assert out["scan_truncated"] is False
    assert out["median_size"] > 0


# ── C8: .h5mu detection + real reader ─────────────────────────────────────────

def _audit_agent():
    from aria.agents.data_audit_agent import DataAuditAgent
    return DataAuditAgent.__new__(DataAuditAgent)


def test_data_audit_detects_and_classifies_h5mu(tmp_path):
    agent = _audit_agent()
    p = tmp_path / "hc11_paired.h5mu"
    p.write_bytes(b"\x00")
    # scanned by extension
    assert p in agent._scan_directory(tmp_path)
    # and classified as the scATAC/chromatin entry, not "unknown"
    classified = agent._classify_files([p])
    assert str(p) in classified.get("scATAC", [])
    assert str(p) not in classified.get("unknown", [])


def test_read_h5mu_summary_file_not_found(tmp_path):
    from aria.utils.mudata_io import read_h5mu_summary
    res = read_h5mu_summary(str(tmp_path / "missing.h5mu"))
    assert res["status"] == "error"
    assert res["error_type"] == "FileNotFound"


def test_read_h5mu_summary_blocks_without_tooling(tmp_path, monkeypatch):
    from aria.utils import mudata_io
    p = tmp_path / "x.h5mu"
    p.write_bytes(b"\x00")
    monkeypatch.setattr(mudata_io, "_import_mudata", lambda: None)
    res = mudata_io.read_h5mu_summary(str(p))
    assert res["status"] == "error"
    assert res["error_type"] == "MuDataToolingMissing"


def test_chromatin_qc_h5mu_routing_blocks_without_tooling(tmp_path, monkeypatch):
    from aria.scripts import chromatin_qc as cq
    from aria.utils import mudata_io
    p = tmp_path / "hc11_paired.h5mu"
    p.write_bytes(b"\x00")
    monkeypatch.setattr(mudata_io, "_import_mudata", lambda: None)
    res = cq.chromatin_qc({"data_type": "scATAC", "files": [str(p)]})
    assert res["status"] == "error"
    # honest structured blocker, not a fabricated QC result
    assert res["error_type"] == "MissingDependency"


def test_read_h5mu_real_roundtrip(tmp_path):
    md = pytest.importorskip("mudata")
    import anndata as ad
    import numpy as np
    from aria.utils.mudata_io import read_h5mu_summary, read_h5mu_atac

    rna = ad.AnnData(np.ones((10, 5), dtype="float32"))
    atac = ad.AnnData(np.ones((10, 8), dtype="float32"))
    atac.obs["atac_fragments"] = np.arange(10, dtype=float) + 100.0
    mdata = md.MuData({"rna": rna, "atac": atac})
    p = tmp_path / "paired.h5mu"
    mdata.write(str(p))

    summ = read_h5mu_summary(str(p))
    assert summ["status"] == "success"
    assert summ["is_paired_rna_atac"] is True
    assert summ["atac_modality"] == "atac" and summ["rna_modality"] == "rna"

    atac_res = read_h5mu_atac(str(p))
    assert atac_res["status"] == "success"
    assert atac_res["n_obs"] == 10 and atac_res["n_vars"] == 8
    assert atac_res["fragment_count_col"] == "atac_fragments"
    assert atac_res["median_fragments_per_cell"] is not None


def test_chromatin_qc_h5mu_real_is_honest(tmp_path):
    md = pytest.importorskip("mudata")
    import anndata as ad
    import numpy as np
    from aria.scripts.chromatin_qc import chromatin_qc

    rna = ad.AnnData(np.ones((12, 5), dtype="float32"))
    atac = ad.AnnData(np.ones((12, 9), dtype="float32"))
    atac.obs["atac_fragments"] = np.arange(12, dtype=float) + 50.0
    mdata = md.MuData({"rna": rna, "atac": atac})
    p = tmp_path / "hc11_paired.h5mu"
    mdata.write(str(p))

    res = chromatin_qc({"data_type": "scATAC", "files": [str(p)]})
    assert res["status"] == "success"
    assert res["input_kind"] == "h5mu"
    assert res["n_cells"] == 12 and res["n_peaks"] == 9
    # not-computed metrics stay None, not fabricated
    assert res["frip"] is None and res["tss_enrichment"] is None
    assert res["qc_complete"] is False
    assert res["pass_qc"] is None
    assert any("tss_enrichment" in m for m in res["metrics_not_computed"])


def test_chromatin_qc_real_hc11_h5mu_validation_input_is_honest():
    """Dataset-gated v4.6 entry-path validation.

    This does not promote scATAC to production. It fixes the expected behavior
    for the real `.h5mu` validation input: ARIA can read the ATAC modality and
    report matrix dimensions, while FRiP/TSS remain not-computed until the v4.6
    chromatin stack adds peak calling and reference-backed TSS enrichment.
    """
    pytest.importorskip("mudata")
    from aria.scripts.chromatin_qc import chromatin_qc

    default_path = (
        "/home/medusa/Samael/Erosion/data_inputs/"
        "muon_processed/hc11_paired.h5mu"
    )
    p = Path(os.environ.get("ARIA_V46_H5MU_VALIDATION_INPUT", default_path))
    if not p.exists():
        pytest.skip(f"v4.6 .h5mu validation input not present: {p}")

    res = chromatin_qc({
        "data_type": "scATAC",
        "files": [str(p)],
        "genome": "mm10",
        "organism": "Mus musculus",
    })

    assert res["status"] == "success"
    assert res["input_kind"] == "h5mu"
    assert res["atac_modality"] == "atac"
    assert res["n_cells"] == 3143
    assert res["n_peaks"] == 60990
    assert res["frip"] is None
    assert res["tss_enrichment"] is None
    assert res["qc_complete"] is False
    assert res["pass_qc"] is None
    missing = "; ".join(res["metrics_not_computed"])
    assert "tss_enrichment" in missing
    assert "frip" in missing


# ── B12: no dead env-alias path ──────────────────────────────────────────────

def test_resolve_env_does_not_read_an_alias_file():
    from aria.utils import environment_manager as em
    src = inspect.getsource(em.EnvironmentManager._resolve_env)
    # The dead alias-resolution branch read+parsed ~/.aria/env_aliases.json,
    # which no code ever writes. Resolution must no longer load any file (B12).
    assert "read_text" not in src
    assert "json.loads" not in src


def test_resolve_env_ignores_any_alias_file(tmp_path, monkeypatch):
    from aria.utils import environment_manager as em
    mgr = em.EnvironmentManager.__new__(em.EnvironmentManager)
    mgr._conda_ok = True
    mgr._env_cache = {"rna": False}          # rna env not installed
    mgr._env_cache_key = "process"           # served from cache, no subprocess
    aria_dir = tmp_path / ".aria"
    aria_dir.mkdir(parents=True)
    (aria_dir / "env_aliases.json").write_text('{"aria-rna-env": "hijacked-env"}')
    monkeypatch.setenv("HOME", str(tmp_path))
    resolved = mgr._resolve_env("rna")
    assert resolved == mgr.FALLBACK_ENV
    assert resolved != "hijacked-env"
