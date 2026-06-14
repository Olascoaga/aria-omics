"""W0.1 (scATAC P0 preprint plan, 2026-06-14): chromatin figure-generation scaffold.

The ChromatinNarrator previously returned no figures (`figures()` -> []). W0.1 adds
`_narrative_chromatin.generate_figures` (mirroring the scRNA figure pipeline: render
in the conda stack via EnvironmentManager, attach PNG paths to findings["figures"])
and surfaces them through the narrator. No fabrication: a figure appears only when
its data + env exist.
"""

from __future__ import annotations

from pathlib import Path

import pytest


def _scatac_envelope(inner: dict) -> dict:
    """A ChromatinAgent envelope with the per-modality wrapper unwrap expects."""
    return {"findings": {"scATAC": {"findings": inner}}}


def test_narrator_figures_surface_findings_figures():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    env = _scatac_envelope({"figures": {"umap_leiden": "/tmp/umap_leiden.png"}})
    figs = ChromatinNarrator().figures("chromatin_agent", env)
    assert {"umap_leiden"} == {f["id"] for f in figs}
    f = figs[0]
    assert f["path"] == "/tmp/umap_leiden.png"
    assert "caption" in f


def test_narrator_figures_empty_when_no_figures():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    env = _scatac_envelope({"qc": {"status": "success"}})
    assert ChromatinNarrator().figures("chromatin_agent", env) == []


def test_live_findings_returns_mutable_nested_dict():
    from aria.agents import _narrative_chromatin
    from aria.agents.narrative.narrators.chromatin import unwrap_chromatin_findings

    env = _scatac_envelope({"qc": {"status": "success"}})
    live = _narrative_chromatin.live_findings(env)
    live["figures"] = {"umap_leiden": "/tmp/x.png"}
    # The mutation must be visible through unwrap (what the narrator reads).
    assert unwrap_chromatin_findings(env).get("figures") == {"umap_leiden": "/tmp/x.png"}


def test_generate_figures_skips_without_env_or_h5ad(tmp_path):
    from aria.agents import _narrative_chromatin

    findings = {"lsi": {"output_path": "/x.h5ad"}}
    out = _narrative_chromatin.generate_figures(
        findings, "/x.h5ad", tmp_path / "figures", env_manager=None)
    assert out["figures"] == {}

    findings2 = {"lsi": {"output_path": None}}
    out2 = _narrative_chromatin.generate_figures(
        findings2, None, tmp_path / "figures2", env_manager=object())
    assert out2["figures"] == {}


class _FakeEM:
    def __init__(self):
        self.calls = []

    def run_in_stack(self, stack, script_path, params):
        self.calls.append({"stack": stack, "script_path": script_path,
                           "params": params})
        out = Path(params["output_dir"]) / "leiden.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(b"\x89PNG")
        return {"status": "success", "figures": {"leiden": str(out)}}


def test_generate_figures_renders_umap_in_chromatin_stack(tmp_path):
    from aria.agents import _narrative_chromatin

    em = _FakeEM()
    findings = {"lsi": {"output_path": "/clustered.h5ad"},
                "differential_accessibility": {"groupby": "leiden"}}
    out = _narrative_chromatin.generate_figures(
        findings, "/clustered.h5ad", tmp_path / "figures", env_manager=em)

    assert out["figures"]["umap_leiden"].endswith("leiden.png")
    assert len(em.calls) == 1
    call = em.calls[0]
    assert call["stack"] == "chromatin"
    assert call["script_path"] == "aria/scripts/rna_figure_umap.py"
    assert "leiden" in call["params"]["color_by"]
    assert call["params"]["h5ad_path"] == "/clustered.h5ad"


def test_report_builder_wires_chromatin_figures_to_narrator(tmp_path, monkeypatch):
    import aria.utils.environment_manager as em_mod
    from aria.agents.narrative.report_builder import ReportBuilderMixin
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    monkeypatch.setattr(em_mod, "env_manager", _FakeEM(), raising=False)

    env = _scatac_envelope({
        "lsi": {"output_path": str(tmp_path / "clustered.h5ad")},
        "differential_accessibility": {"groupby": "leiden"},
    })
    agent_results = {"chromatin_agent": env}

    rb = ReportBuilderMixin.__new__(ReportBuilderMixin)
    rb._generate_chromatin_figures(agent_results, tmp_path)

    figs = ChromatinNarrator().figures("chromatin_agent", env)
    assert any(f["id"].startswith("umap_") for f in figs)
