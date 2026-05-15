from __future__ import annotations

import os
import json
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_package_compiles():
    result = subprocess.run(
        [sys.executable, "-m", "compileall", "-q", "aria"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr or result.stdout


def test_aria_env_loader_reads_export_file(tmp_path, monkeypatch):
    from aria.utils import env_loader

    env_file = tmp_path / ".env"
    env_file.write_text(
        "\n".join([
            "# private ARIA keys",
            "export ANTHROPIC_API_KEY='anthropic-test'",
            'OPENAI_API_KEY="openai-test"',
            "GEMINI_API_KEY=gemini-test",
        ]),
        encoding="utf-8",
    )

    for key in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "GEMINI_API_KEY"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setattr(env_loader, "_LOADED_PATHS", set())

    loaded = env_loader.load_aria_env(env_file)

    assert loaded == {
        "ANTHROPIC_API_KEY": "anthropic-test",
        "OPENAI_API_KEY": "openai-test",
        "GEMINI_API_KEY": "gemini-test",
    }
    assert os.environ["ANTHROPIC_API_KEY"] == "anthropic-test"
    assert os.environ["OPENAI_API_KEY"] == "openai-test"
    assert os.environ["GEMINI_API_KEY"] == "gemini-test"


def test_aria_env_loader_preserves_existing_env(tmp_path, monkeypatch):
    from aria.utils import env_loader

    env_file = tmp_path / ".env"
    env_file.write_text("export ANTHROPIC_API_KEY=file-value\n",
                        encoding="utf-8")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "terminal-value")
    monkeypatch.setattr(env_loader, "_LOADED_PATHS", set())

    loaded = env_loader.load_aria_env(env_file)

    assert loaded == {}
    assert os.environ["ANTHROPIC_API_KEY"] == "terminal-value"


def test_llm_provider_loads_aria_env_file(tmp_path, monkeypatch):
    from aria.llm.provider import LLMProvider
    from aria.utils import env_loader

    env_file = tmp_path / ".env"
    env_file.write_text("export ANTHROPIC_API_KEY=provider-test\n",
                        encoding="utf-8")
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    monkeypatch.setenv("ARIA_ENV_FILE", str(env_file))
    monkeypatch.setattr(env_loader, "_LOADED_PATHS", set())

    provider = LLMProvider()

    assert provider.api_keys["anthropic"] == "provider-test"


def test_bulk_rna_legacy_script_passes():
    env = os.environ.copy()
    env["ARIA_ALLOW_MOCKS"] = "1"
    result = subprocess.run(
        [sys.executable, "tests/test_bulk_rna.py"],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_h5ad_obs_design_inference(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    from aria.agents.data_audit_agent import DataAuditAgent

    obs = pd.DataFrame(
        {
            "age_group": ["young"] * 12 + ["young"] * 12
                         + ["old"] * 12 + ["old"] * 12,
            "donor_id": ["y1"] * 12 + ["y2"] * 12
                        + ["o1"] * 12 + ["o2"] * 12,
            "subclass": ["OPC", "Microglia"] * 24,
            "Gender": ["F"] * 24 + ["M"] * 24,
        },
        index=[f"cell_{i}" for i in range(48)],
    )
    adata = ad.AnnData(
        X=np.ones((48, 4), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(4)]),
    )
    h5ad_path = tmp_path / "hippo.h5ad"
    adata.write_h5ad(h5ad_path)

    agent = DataAuditAgent.__new__(DataAuditAgent)
    inferred = agent._infer_h5ad_design(
        [str(h5ad_path)], "compare young versus old oligodendrocytes"
    )

    assert inferred["source"] == "h5ad_obs"
    assert inferred["condition_col"] == "age_group"
    assert inferred["replicate_col"] == "donor_id"
    assert inferred["groupby_col"] == "subclass"
    assert inferred["groups"] == {
        "old": ["o1", "o2"],
        "young": ["y1", "y2"],
    }
    assert inferred["covariates"] == ["Gender"]
    assert inferred["pseudobulk"]["from_obs"] is True

    summary = agent._build_checkpoint_summary(
        {"scRNA": [str(h5ad_path)]}, "unknown", "unknown", [],
        inferred_design=inferred,
    )
    assert "covariates: Gender" in summary


def test_data_audit_infers_hg38_from_homo_sapiens_question(tmp_path):
    from aria.agents.data_audit_agent import DataAuditAgent

    agent = DataAuditAgent.__new__(DataAuditAgent)
    genome, organism = agent._infer_genome_organism(
        [], tmp_path, "Compare aging signatures in Homo sapiens hippocampus"
    )

    assert genome == "hg38"
    assert organism == "Homo sapiens"


def test_data_audit_infers_human_h5ad_from_gene_symbols(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    from aria.agents.data_audit_agent import DataAuditAgent

    genes = [
        "SAMD11", "ISG15", "TMEM88B", "PRDM16", "MEGF6", "C1QA",
        "C1QB", "C1QC", "NCMAP", "C1orf141",
    ]
    adata = ad.AnnData(
        X=np.ones((4, len(genes)), dtype=np.float32),
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(4)]),
        var=pd.DataFrame(index=genes),
    )
    h5ad_path = tmp_path / "hippocampus_RNA.h5ad"
    adata.write_h5ad(h5ad_path)

    genome, organism = DataAuditAgent._infer_h5ad_genome_organism(
        [str(h5ad_path)]
    )

    assert genome == "hg38"
    assert organism == "Homo sapiens"


def test_tui_drains_stale_pasted_lines(monkeypatch):
    from aria import tui

    class FakeStdin:
        def __init__(self, lines):
            self.lines = list(lines)

        def isatty(self):
            return True

        def readline(self):
            return self.lines.pop(0) if self.lines else ""

    fake_stdin = FakeStdin([
        "second line\n",
        "\n",
        "third line\n",
    ])

    def fake_select(reads, writes, errors, timeout):
        return (reads, writes, errors) if fake_stdin.lines else ([], [], [])

    monkeypatch.setattr(sys, "stdin", fake_stdin)
    monkeypatch.setattr(tui.select, "select", fake_select)

    assert tui._read_pasted_stdin_lines() == [
        "second line",
        "",
        "third line",
    ]


def test_tui_multiline_question_uses_end_sentinel(monkeypatch):
    from aria import tui

    answers = iter([
        "first line",
        "",
        "third line",
        "END",
    ])
    monkeypatch.setattr(tui.console, "input", lambda prompt: next(answers))

    assert tui._read_multiline_question() == "first line\n\nthird line"


def test_data_audit_ignores_stale_aria_h5ad_intermediates():
    from aria.agents.data_audit_agent import DataAuditAgent

    agent = DataAuditAgent.__new__(DataAuditAgent)
    classified = {
        "scRNA": [
            "/data/GSE278576_hippocampus_RNA.h5ad",
            "/data/qc_filtered.h5ad",
            "/data/clustered.h5ad",
        ],
        "unknown": ["/data/notes.txt"],
    }

    filtered, ignored = agent._filter_aria_intermediate_outputs(classified)

    assert filtered["scRNA"] == ["/data/GSE278576_hippocampus_RNA.h5ad"]
    assert ignored == ["/data/qc_filtered.h5ad", "/data/clustered.h5ad"]


def test_data_audit_keeps_intermediate_if_it_is_the_only_scrna_input():
    from aria.agents.data_audit_agent import DataAuditAgent

    agent = DataAuditAgent.__new__(DataAuditAgent)
    classified = {"scRNA": ["/data/qc_filtered.h5ad"]}

    filtered, ignored = agent._filter_aria_intermediate_outputs(classified)

    assert filtered["scRNA"] == ["/data/qc_filtered.h5ad"]
    assert ignored == []


def test_scrna_pseudobulk_uses_h5ad_obs_design(tmp_path):
    from aria.agents.scrna_agent import scRNAAgent

    h5ad_path = tmp_path / "annotated.h5ad"
    h5ad_path.write_text("placeholder")
    calls = []

    class FakeEnv:
        def run_in_stack(self, **kwargs):
            calls.append(kwargs)
            script = kwargs["script_path"]
            if script.endswith("rna_diff_abundance.py"):
                # T1.1 wire-up: agent runs DA before pseudobulk. Return a
                # benign 'no significant shift' result so composition
                # covariate defaults to OFF for this test.
                return {
                    "status": "success",
                    "method": "poisson_offset_glm",
                    "groupby": kwargs["params"]["groupby"],
                    "condition_col": kwargs["params"]["condition_col"],
                    "replicate_col": kwargs["params"]["replicate_col"],
                    "covariates": kwargs["params"].get("covariates", []),
                    "significance_alpha": 0.10,
                    "any_significant": False,
                    "per_comparison": {},
                    "n_replicates_per_group": {},
                    "warnings": [],
                }
            if script.endswith("rna_pseudobulk_de.py"):
                params = kwargs["params"]
                return {
                    "status": "success",
                    "groupby": params["groupby"],
                    "condition_col": params["condition_col"],
                    "replicate_col": params["replicate_col"],
                    "covariates": params["covariates"],
                    "thresholds": {},
                    "n_groups": 0,
                    "per_group": {},
                }
            raise AssertionError(f"unexpected script: {script}")

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.env = FakeEnv()
    agent.publish_finding = lambda *args, **kwargs: None
    agent._inject_condition_obs = lambda *args, **kwargs: (_ for _ in ()).throw(
        AssertionError("condition injection should be skipped")
    )

    result = agent._run_pseudobulk(
        "exp",
        str(h5ad_path),
        {
            "organism": "Homo sapiens",
            "design": {
                "groups": {"old": ["o1", "o2"], "young": ["y1", "y2"]},
                "main_factor": "age_group",
                "pseudobulk": {
                    "from_obs": True,
                    "condition_col": "age_group",
                    "replicate_col": "donor_id",
                    "groupby_col": "subclass",
                    "covariates": ["Gender"],
                    "comparisons": [["old", "young"]],
                },
            },
        },
        {"summary": "compare young vs old"},
        {},
    )

    assert result["status"] == "success"
    # Last call is the pseudobulk; first call is the DA pre-pass.
    pb_call = next(
        c for c in calls if c["script_path"].endswith("rna_pseudobulk_de.py")
    )
    params = pb_call["params"]
    assert params["data_path"] == str(h5ad_path)
    assert params["condition_col"] == "age_group"
    assert params["replicate_col"] == "donor_id"
    assert params["groupby"] == "subclass"
    assert params["covariates"] == ["Gender"]
    # With DA returning any_significant=False, composition covariate is OFF.
    assert params["composition_covariate"] is False


def test_rna_qc_uses_existing_h5ad_obs_metrics_for_processed_input(tmp_path, monkeypatch):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "numba-cache"))

    from aria.scripts.rna_qc import rna_qc

    n_cells = 80
    obs = pd.DataFrame(
        {
            "nFeature_RNA": [900] * n_cells,
            "nCount_RNA": [2200] * n_cells,
            "percent.mt": [0.2] * n_cells,
            "age_group": ["20-39"] * 40 + ["80-100"] * 40,
            "orig.ident": ["d1"] * 20 + ["d2"] * 20 + ["d3"] * 20 + ["d4"] * 20,
            "subclass": ["Oligo", "OPC"] * 40,
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    adata = ad.AnnData(
        X=np.random.default_rng(0).normal(size=(n_cells, 20)),
        obs=obs,
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(20)]),
    )
    input_path = tmp_path / "processed.h5ad"
    adata.write_h5ad(input_path)

    result = rna_qc(
        {
            "data_path": str(input_path),
            "organism": "Homo sapiens",
            "output_dir": str(tmp_path),
        }
    )

    assert result["status"] == "success"
    assert result["n_cells_after"] == n_cells
    assert result["scrublet"]["ran"] is False
    assert any("existing h5ad obs QC metrics" in w for w in result["warnings"])


def test_rna_qc_empty_h5ad_returns_structured_error(tmp_path, monkeypatch):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "numba-cache"))

    from aria.scripts.rna_qc import rna_qc

    empty = ad.AnnData(
        X=np.empty((0, 0), dtype=np.float32),
        obs=pd.DataFrame(index=[]),
        var=pd.DataFrame(index=[]),
    )
    input_path = tmp_path / "qc_filtered.h5ad"
    empty.write_h5ad(input_path)

    result = rna_qc({
        "data_path": str(input_path),
        "organism": "Homo sapiens",
        "output_dir": str(tmp_path),
    })

    assert result["status"] == "error"
    assert result["error_type"] == "EmptyAnnData"
    assert result["n_cells_after"] == 0


def test_clustering_cache_requires_matching_parameters():
    from aria.scripts.rna_clustering import _cache_matches, _cache_params

    params = {
        "data_path": "/tmp/qc.h5ad",
        "resolution": 0.4,
        "n_neighbors": 15,
        "max_cells": 100_000,
        "seed": 0,
    }
    expected = _cache_params(params)

    assert _cache_matches(
        {"cache_version": 4, "cache_params": expected},
        expected,
    )
    stale = {**expected, "resolution": 0.8}
    assert not _cache_matches(
        {"cache_version": 4, "cache_params": stale},
        expected,
    )
    assert not _cache_matches(
        {"cache_version": 3, "cache_params": expected},
        expected,
    )


def test_integration_cache_requires_matching_parameters():
    from aria.scripts.rna_integration import _cache_matches, _cache_params

    params = {
        "data_path": "/tmp/qc.h5ad",
        "batch_col": "donor_id",
        "max_cells": 250_000,
        "seed": 0,
    }
    expected = _cache_params(params)

    assert _cache_matches(
        {"cache_version": 2, "cache_params": expected},
        expected,
    )
    stale = {**expected, "batch_col": "sample_id"}
    assert not _cache_matches(
        {"cache_version": 2, "cache_params": stale},
        expected,
    )
    assert not _cache_matches(
        {"cache_version": 1, "cache_params": expected},
        expected,
    )


def test_qc_cache_requires_matching_parameters():
    from aria.scripts.rna_qc import _cache_matches, _cache_params

    params = {
        "data_path": "/tmp/raw.h5ad",
        "organism": "Homo sapiens",
        "min_genes": 200,
        "run_scrublet": False,
        "biological_context": {"user_question": "compare stressed cells"},
    }
    expected = _cache_params(params)

    assert _cache_matches(
        {"cache_version": 2, "cache_params": expected},
        expected,
    )
    stale = {**expected, "min_genes": 500}
    assert not _cache_matches(
        {"cache_version": 2, "cache_params": stale},
        expected,
    )


def test_concat_cache_requires_matching_manifest():
    from aria.scripts.rna_concat import _cache_matches, _cache_params

    params = {
        "samples": [
            {"path": "/tmp/a.h5ad", "sample_id": "a", "condition": "ctrl"},
            {"path": "/tmp/b.h5ad", "sample_id": "b", "condition": "treated"},
        ],
        "join": "inner",
    }
    expected = _cache_params(params)

    assert _cache_matches(
        {"cache_version": 2, "cache_params": expected},
        expected,
    )
    stale = _cache_params({
        **params,
        "samples": [
            {"path": "/tmp/a.h5ad", "sample_id": "a", "condition": "ctrl"},
            {"path": "/tmp/b.h5ad", "sample_id": "b", "condition": "ctrl"},
        ],
    })
    assert not _cache_matches(
        {"cache_version": 2, "cache_params": stale},
        expected,
    )


def test_marker_fallback_annotation_is_explicit_and_conservative():
    from aria.agents.scrna_agent import scRNAAgent

    labels = scRNAAgent._marker_based_annotation({
        "0": ["MBP", "PLP1", "MOG", "ACTB"],
        "1": ["GENE_A", "GENE_B"],
    })

    assert labels["0"]["cell_type"] == "Oligodendrocyte"
    assert labels["0"]["annotation_source"] == "marker_fallback"
    assert labels["0"]["confidence"] == "medium"
    assert labels["1"]["cell_type"] == "Unresolved cluster 1"
    assert labels["1"]["confidence"] == "low"


def test_sample_id_strips_intermediate_stems_to_avoid_recursive_naming():
    from aria.agents.scrna_agent import scRNAAgent

    # Real 10x stems pass through cleanly.
    assert (
        scRNAAgent._sample_id_from_path(
            "/tmp/GSE278576_hc11_raw_feature_bc_matrix.h5"
        )
        == "GSE278576_hc11"
    )

    # `qc_filtered_<name>` prefix is stripped to avoid double-prefixing.
    assert (
        scRNAAgent._sample_id_from_path(
            "/tmp/scrna_multi/qc_filtered_GSE278576_hc11.h5ad"
        )
        == "GSE278576_hc11"
    )

    # ARIA intermediate stems (annotated, clustered, qc_filtered, ...) fall
    # back to a hash-based id so we never produce qc_filtered_annotated.h5ad
    # or qc_filtered_qc_filtered.h5ad on accidental re-feed.
    for stem in ("annotated", "clustered", "qc_filtered", "trajectory",
                 "concatenated", "integrated", "with_condition"):
        sid = scRNAAgent._sample_id_from_path(f"/tmp/scrna_multi/{stem}.h5ad")
        assert sid.startswith("sample_"), sid
        assert stem not in sid


def test_predefined_celltype_col_picks_subclass_from_obs(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    from aria.agents.scrna_agent import scRNAAgent

    obs = pd.DataFrame(
        {"subclass": ["OPC", "OPC", "Astro", "Astro", "Micro", "Micro"]},
        index=[f"c_{i}" for i in range(6)],
    )
    adata = ad.AnnData(
        X=np.ones((6, 3), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=[f"g_{i}" for i in range(3)]),
    )
    h5ad_path = tmp_path / "annotated.h5ad"
    adata.write_h5ad(h5ad_path)

    agent = scRNAAgent.__new__(scRNAAgent)
    exp_ctx = {
        "design": {"pseudobulk": {"groupby_col": "subclass"}},
    }
    picked = agent._predefined_celltype_col(str(h5ad_path), exp_ctx)
    assert picked == "subclass"

    # Should reject when the obs column is missing in the working h5ad.
    exp_ctx_missing = {
        "design": {"pseudobulk": {"groupby_col": "celltype_missing"}},
    }
    assert agent._predefined_celltype_col(str(h5ad_path), exp_ctx_missing) is None

    # Should reject "leiden" as a non-informative fallback marker.
    exp_ctx_leiden = {"design": {"pseudobulk": {"groupby_col": "leiden"}}}
    assert agent._predefined_celltype_col(str(h5ad_path), exp_ctx_leiden) is None


def test_scrna_infers_and_materializes_cell_focus(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    from aria.agents.scrna_agent import scRNAAgent

    obs = pd.DataFrame(
        {"subclass": ["OPC"] * 3 + ["Oligo"] * 4 + ["Astro"] * 5},
        index=[f"c_{i}" for i in range(12)],
    )
    adata = ad.AnnData(
        X=np.ones((12, 3), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=[f"g_{i}" for i in range(3)]),
    )
    h5ad_path = tmp_path / "hippo.h5ad"
    adata.write_h5ad(h5ad_path)

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.publish_finding = lambda *args, **kwargs: None
    exp_ctx = {
        "user_question": "focus on oligodendrocyte trajectories from OPCs",
        "design": {"pseudobulk": {"groupby_col": "subclass"}},
    }
    focus = agent._prepare_focused_h5ads(
        "exp-focus",
        [str(h5ad_path)],
        exp_ctx,
        {"summary": "OPC to oligodendrocyte lineage"},
    )

    assert focus["status"] == "success"
    assert focus["values"] == ["OPC", "Oligo"]
    assert focus["n_cells_before"] == 12
    assert focus["n_cells_after"] == 7
    focused = ad.read_h5ad(focus["files"][0])
    assert set(focused.obs["subclass"].astype(str)) == {"OPC", "Oligo"}
    assert focused.n_obs == 7


def test_scrna_cell_focus_skips_when_all_groups_requested(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    from aria.agents.scrna_agent import scRNAAgent

    obs = pd.DataFrame(
        {"subclass": ["OPC", "Oligo", "Astro"] * 2},
        index=[f"c_{i}" for i in range(6)],
    )
    adata = ad.AnnData(
        X=np.ones((6, 2), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    h5ad_path = tmp_path / "all_groups.h5ad"
    adata.write_h5ad(h5ad_path)

    focus = scRNAAgent._infer_cell_focus_values(
        [str(h5ad_path)],
        "subclass",
        {"user_question": "compare OPC Oligo and Astro"},
        {"summary": ""},
    )

    assert focus == set()


def test_scrna_cell_focus_does_not_expand_generic_neurons(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    from aria.agents.scrna_agent import scRNAAgent

    obs = pd.DataFrame(
        {"subclass": ["OPC", "Oligo", "CA1", "DG"] * 2},
        index=[f"c_{i}" for i in range(8)],
    )
    adata = ad.AnnData(
        X=np.ones((8, 2), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    h5ad_path = tmp_path / "hippo_neurons.h5ad"
    adata.write_h5ad(h5ad_path)

    focus = scRNAAgent._infer_cell_focus_values(
        [str(h5ad_path)],
        "subclass",
        {"user_question": "hippocampus neurons and oligodendrocytes"},
        {"summary": "oligodendrocyte aging"},
    )

    assert focus == set()


def test_scrna_cell_focus_uses_only_explicit_focus_clause(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    from aria.agents.scrna_agent import scRNAAgent

    obs = pd.DataFrame(
        {"subclass": ["Astro", "Oligo", "OPC"] * 2},
        index=[f"c_{i}" for i in range(6)],
    )
    adata = ad.AnnData(
        X=np.ones((6, 2), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    h5ad_path = tmp_path / "astro_focus.h5ad"
    adata.write_h5ad(h5ad_path)

    focus = scRNAAgent._infer_cell_focus_values(
        [str(h5ad_path)],
        "subclass",
        {
            "user_question": (
                "Focus only on astrocytes. Do not run oligodendrocyte "
                "trajectory analysis."
            )
        },
        {"summary": "astrocyte aging and oligodendrocyte trajectory not needed"},
    )

    assert focus == {"Astro"}


def test_design_intelligence_scrna_microglia_feasibility(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    from aria.agents.design_intelligence import DesignIntelligence

    obs = pd.DataFrame(
        {
            "subclass": ["Microglia"] * 6 + ["Oligo"] * 6,
            "age_group": ["20-39"] * 3 + ["80-100"] * 3
                         + ["20-39"] * 3 + ["80-100"] * 3,
            "orig.ident": ["y1", "y2", "y3", "o1", "o2", "o3"] * 2,
            "Gender": ["F", "M", "F", "M", "F", "M"] * 2,
            "nCount_RNA": [1000] * 12,
        },
        index=[f"c_{i}" for i in range(12)],
    )
    adata = ad.AnnData(
        X=np.ones((12, 3), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["APOE", "TREM2", "C1QA"]),
    )
    h5ad_path = tmp_path / "microglia.h5ad"
    adata.write_h5ad(h5ad_path)

    exp_ctx = {
        "modalities": {"scRNA": [str(h5ad_path)]},
        "user_question": "focus only on microglia aging signatures",
        "design": {
            "groups": {"20-39": ["y1", "y2", "y3"],
                       "80-100": ["o1", "o2", "o3"]},
            "pseudobulk": {
                "condition_col": "age_group",
                "replicate_col": "orig.ident",
                "groupby_col": "subclass",
                "covariates": ["Gender"],
            },
        },
    }

    di = DesignIntelligence().evaluate(
        exp_ctx, {"summary": "microglia aging and ligand receptor signaling"}
    )

    assert di["status"] == "success"
    assert any("pseudobulk" in item.lower() for item in di["recommended"])
    assert any("Microglia" in item for item in di["recommended"])
    assert any("LIANA" in item for item in di["unsupported"])
    assert any("RNA velocity" in item for item in di["unsupported"])


def test_orchestrator_plan_summary_includes_design_intelligence():
    from aria.agents.orchestrator_agent import OrchestratorAgent

    orch = OrchestratorAgent.__new__(OrchestratorAgent)
    summary = orch._format_plan_summary({
        "steps": [{"order": 1, "agent": "scrna_agent",
                   "analysis": "Analyze scRNA"}],
        "estimated_complexity": "low",
        "rationale": "test",
        "design_intelligence": {
            "recommended": ["Microglia pseudobulk DE"],
            "optional": ["Pathway enrichment"],
            "unsupported": ["RNA velocity is not supported"],
            "warnings": [],
        },
    })

    assert "Design Intelligence" in summary
    assert "Microglia pseudobulk DE" in summary
    assert "RNA velocity is not supported" in summary
    assert "recommended analyses" in summary
    assert "optional supported analyses" in summary


def test_scrna_annotation_from_obs_skips_celltypist(tmp_path):
    from aria.agents.scrna_agent import scRNAAgent

    h5ad_path = tmp_path / "annotated.h5ad"
    h5ad_path.write_text("placeholder")

    class FakeEnv:
        def run_in_stack(self, **kwargs):
            raise AssertionError(
                f"no subprocess should run; got {kwargs.get('script_path')}"
            )

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.env = FakeEnv()
    agent.publish_finding = lambda *args, **kwargs: None

    annotation = agent._annotation_from_obs(
        "exp",
        str(h5ad_path),
        cell_type_col="subclass",
        cluster_sizes={"OPC": 100, "Astro": 80, "Micro": 50},
        top_markers={"OPC": ["PDGFRA"], "Astro": ["AQP4"], "Micro": ["P2RY12"]},
    )

    assert annotation["label_col"] == "subclass"
    assert annotation["annotation_source"] == "input_obs"
    assert annotation["celltypist"]["ran"] is False
    assert set(annotation["cell_types"]) == {"OPC", "Astro", "Micro"}
    for entry in annotation["cell_types"].values():
        assert entry["annotation_source"] == "input_obs"


def test_scrna_pseudobulk_skips_when_annotation_unrecoverable(tmp_path):
    from aria.agents.scrna_agent import scRNAAgent

    h5ad_path = tmp_path / "clustered.h5ad"
    h5ad_path.write_text("placeholder")

    class FakeEnv:
        def run_in_stack(self, **kwargs):
            raise AssertionError("pseudobulk dispatch should be skipped")

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.env = FakeEnv()
    agent.publish_finding = lambda *args, **kwargs: None

    annotation = {
        "annotation_source": "unrecoverable_matrix",
        "label_col": None,
        "celltypist": {"ran": False, "error_type": "InvalidExpressionMatrix"},
    }
    result = agent._run_pseudobulk(
        "exp",
        str(h5ad_path),
        {
            "organism": "Homo sapiens",
            "design": {
                "groups": {"young": ["y1", "y2"], "old": ["o1", "o2"]},
                "main_factor": "age_group",
                "pseudobulk": {
                    "from_obs": True,
                    "condition_col": "age_group",
                    "replicate_col": "donor_id",
                    "groupby_col": None,
                    "covariates": [],
                    "comparisons": [["old", "young"]],
                },
            },
        },
        {"summary": "young vs old"},
        annotation,
    )
    assert result["status"] == "skipped"
    assert "annotation" in result["reason"]


def test_clustering_skip_leiden_when_cluster_col_provided():
    from aria.scripts.rna_clustering import _cache_params

    expected = _cache_params({
        "data_path": "/tmp/qc.h5ad",
        "resolution": 0.4,
        "cluster_col": "subclass",
    })
    assert expected["cluster_col"] == "subclass"

    expected_default = _cache_params({
        "data_path": "/tmp/qc.h5ad",
        "resolution": 0.4,
    })
    assert expected_default["cluster_col"] == ""


def test_clustering_predefined_groupby_keeps_all_cells(tmp_path, monkeypatch):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "numba-cache"))

    from aria.scripts.rna_clustering import rna_clustering

    n_cells = 12
    adata = ad.AnnData(
        X=np.random.default_rng(1).normal(size=(n_cells, 5)),
        obs=pd.DataFrame(
            {"subclass": ["OPC"] * 4 + ["Astro"] * 4 + ["Micro"] * 4},
            index=[f"cell_{i}" for i in range(n_cells)],
        ),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(5)]),
    )
    input_path = tmp_path / "processed.h5ad"
    adata.write_h5ad(input_path)

    result = rna_clustering({
        "data_path": str(input_path),
        "resolution": 0.5,
        "cluster_col": "subclass",
        "max_cells": 5,
    })

    assert result["status"] == "success"
    assert result["predef_clusters"] is True
    assert result["sketch_used"] is False
    assert result["n_cells_used"] == n_cells
    assert result["groupby"] == "subclass"
    assert result["cluster_sizes"] == {"OPC": 4, "Astro": 4, "Micro": 4}


def test_scrna_narrative_reports_predefined_groupby_not_leiden():
    from aria.agents import _narrative_scrna

    findings = {
        "clustering": {
            "n_clusters": 18,
            "groupby": "subclass",
            "predef_clusters": True,
        },
        "clustering_decision": {
            "groupby": "subclass",
            "predef_clusters": True,
        },
        "pseudobulk_de": {
            "groupby": "subclass",
            "n_groups": 18,
            "thresholds": {"padj_max": 0.05, "lfc_min": 0.5},
            "per_group": {},
        },
        "trajectory": {
            "status": "done",
            "paga": {
                "n_connections": 153,
                "n_strong": 40,
                "strong_threshold": 0.05,
                "max_connectivity": 1.0,
            },
            "pseudotime": {"computed": False},
            "velocity": {"computed": False},
        },
    }

    summary = _narrative_scrna.summarize_scrna_text(findings)
    methods = _narrative_scrna.build_scrna_methods(findings)

    assert "obs['subclass']" in summary
    assert "Leiden clustering was skipped" in summary
    assert "subclass groups" in summary
    assert "subclasss" not in summary
    assert "active lineage transitions" not in summary
    assert "not proof of active differentiation" in summary
    assert "skipped Leiden clustering" in methods


def test_scrna_narrative_surfaces_de_per_cluster_failure():
    from aria.agents import _narrative_scrna

    findings = {
        "clustering": {
            "n_clusters": 18,
            "groupby": "subclass",
            "predef_clusters": True,
        },
        "differential_expression": {
            "status": "error",
            "error_type": "Timeout",
            "details": "Execution exceeded 3600s limit.",
        },
        "pseudobulk_de": {
            "groupby": "subclass",
            "n_groups": 18,
            "thresholds": {"padj_max": 0.05, "lfc_min": 0.5},
            "per_group": {},
        },
    }

    summary = _narrative_scrna.summarize_scrna_text(findings)
    assert "Per-cluster marker discovery" in summary
    assert "Timeout" in summary
    assert "Execution exceeded 3600s limit" in summary
    # The downstream pseudobulk path is preserved as the valid signal.
    assert "pseudobulk" in summary.lower()


def test_cellcomm_table_labels_rank_metric_not_generic_score():
    from aria.agents import _narrative_scrna

    findings = {
        "cell_communication": {
            "status": "done",
            "method": "liana_rank_aggregate (specificity_rank)",
            "top_interactions": [{
                "source": "Microglia",
                "target": "CA1",
                "ligand": "APOE",
                "receptor": "TREM2",
                "score": 0.0,
                "rank": 1,
                "rank_metric": "specificity_rank",
                "cellphone_pval": 0.0,
            }],
        }
    }

    rows = _narrative_scrna.extract_cellcomm_table(findings)
    assert "#1" in rows
    assert "specificity_rank" in rows
    assert "<td><code>0.0</code></td>" not in rows
    assert "<td>—</td>" in rows


def test_integrated_interpretation_sanitizes_long_prompt():
    from aria.agents import _narrative_scrna

    findings = {
        "qc": {"n_cells_after": 1000},
        "cell_types": {"annotation_source": "input_obs"},
        "pseudobulk_de": {"per_group": {}},
    }
    intent = {
        "summary": (
            "Compare aging signatures in hippocampus snRNA-seq. Use age_group "
            "20-39 as young.\n\nRun all applicable single-cell analyses. "
            "Do not run ATAC."
        )
    }

    text = _narrative_scrna.build_scrna_integrated_interpretation(
        findings, intent
    )

    assert "Compare aging signatures in hippocampus snRNA-seq" in text
    assert "Run all applicable" not in text
    assert "Do not run ATAC" not in text


def test_narrative_html_escapes_methods_and_uses_package_version(tmp_path):
    from aria import __version__
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    report = agent._render_html_report(
        experiment_id="exp_test",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={"summary": "test"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={
            "high": [], "medium": [], "low": [], "insufficient": []
        },
        methods="Low-count genes (<10 counts in fewer than 25% of samples).",
        decisions=[],
        agent_results={},
        report_dir=tmp_path / "report",
    )
    html = report.read_text(encoding="utf-8")
    assert f"Generated by ARIA v{__version__}" in html
    assert "&lt;10 counts" in html
    assert "Generated by ARIA v0.2" not in html
    assert "Provenance" in html
    assert (tmp_path / "report" / "methodology.json").exists()


def test_bulk_interpretation_guard_downgrades_causal_language():
    from aria.agents.narrative_agent import NarrativeAgent

    guarded = NarrativeAgent._guard_bulk_interpretation(
        "GeneA acts as a broad upstream regulator. "
        "These data position pathway regulators as hierarchical "
        "gatekeepers. GeneB enforces cell identity."
    )

    assert "acts as a broad upstream regulator" not in guarded
    assert "hierarchical gatekeepers" not in guarded
    assert "GeneB enforces cell identity" not in guarded
    assert "not direct evidence of transcriptional causality" in guarded


def test_findings_table_formats_raw_error_dicts():
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent.__new__(NarrativeAgent)
    html = agent._build_findings_table({
        "high": [],
        "medium": [{
            "summary": {
                "error": "Integration requires at least 2 modalities."
            },
            "agent": "integration_agent",
        }],
        "low": [],
        "insufficient": [],
    })

    assert "Integration was skipped" in html
    assert "{&#x27;error&#x27;" not in html


def test_apply_cluster_labels_writes_real_obs_column(tmp_path):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    from aria.scripts.rna_apply_cluster_labels import rna_apply_cluster_labels

    adata = ad.AnnData(
        X=np.ones((4, 3), dtype=np.float32),
        obs=pd.DataFrame(
            {"leiden": ["0", "0", "1", "1"]},
            index=[f"cell_{i}" for i in range(4)],
        ),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(3)]),
    )
    input_path = tmp_path / "clustered.h5ad"
    adata.write_h5ad(input_path)

    result = rna_apply_cluster_labels({
        "data_path": str(input_path),
        "labels": {
            "0": {"cell_type": "Oligodendrocyte"},
            "1": {"cell_type": "Microglia"},
        },
        "label_col": "cell_type_marker",
        "output_dir": str(tmp_path / "labels"),
    })

    assert result["status"] == "success"
    assert result["label_col"] == "cell_type_marker"
    labeled = ad.read_h5ad(result["output_path"])
    assert list(labeled.obs["cell_type_marker"]) == [
        "Oligodendrocyte",
        "Oligodendrocyte",
        "Microglia",
        "Microglia",
    ]


# ── T1.3: low-power warning at n=2 ────────────────────────────────────────────

def test_pseudobulk_low_power_warning_at_n2(tmp_path):
    """T1.3 — pseudobulk DE with exactly 2 replicates per condition must
    surface low_power_warning on every analyzable block. The DE itself
    still runs when the caller explicitly passes
    min_replicates_per_condition=2; the warning is what changes."""
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")

    rng = np.random.default_rng(0)
    n_cells_per_block = 30
    rows = []
    for cond in ("ctrl", "treat"):
        for rep_i in range(1, 3):  # n=2 replicates per condition
            for ct in ("A", "B"):
                for _ in range(n_cells_per_block):
                    rows.append({
                        "condition":  cond,
                        "donor":      f"{cond}_{rep_i}",
                        "cell_type":  ct,
                    })
    obs = pd.DataFrame(rows)
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    n_genes = 60
    X = rng.negative_binomial(20, 0.5, (len(obs), n_genes)).astype(np.float32)
    mask = (obs["condition"] == "treat") & (obs["cell_type"] == "A")
    X[mask.to_numpy(), :10] *= 5
    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=[f"GENE_{i:04d}" for i in range(n_genes)]),
    )
    input_path = tmp_path / "pb.h5ad"
    adata.write_h5ad(input_path)

    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de
    result = rna_pseudobulk_de({
        "data_path":                     str(input_path),
        "groupby":                       "cell_type",
        "condition_col":                 "condition",
        "replicate_col":                 "donor",
        "comparisons":                   [["treat", "ctrl"]],
        "output_dir":                    str(tmp_path / "out"),
        "min_replicates_per_condition":  2,
        "min_cells_per_pseudosample":    5,
        "padj_max":                      0.5,
        "lfc_min":                       0.0,
    })

    assert result["status"] == "success", result
    per_group = result["per_group"]
    assert per_group, "Expected per_group entries"

    any_success = False
    for group, group_info in per_group.items():
        for comp_key, comp in (group_info.get("per_comparison") or {}).items():
            if comp.get("status") != "success":
                continue
            any_success = True
            assert comp.get("low_power_warning") is True, (
                f"{group} {comp_key} should be flagged low_power; got {comp}"
            )
            assert "n=" in (comp.get("low_power_reason") or "")
            assert comp["n_replicates"] == {"test": 2, "ref": 2}

    assert any_success, "Expected at least one successful per_comparison block"


def test_pseudobulk_default_min_replicates_is_three(tmp_path):
    """T1.3 — the script's default min_replicates_per_condition is now 3.
    A run that omits the param and has only 2 replicates per condition
    must skip the comparison with a clear reason."""
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")

    rng = np.random.default_rng(1)
    rows = []
    for cond in ("ctrl", "treat"):
        for rep_i in range(1, 3):
            for _ in range(30):
                rows.append({"condition": cond,
                             "donor":     f"{cond}_{rep_i}",
                             "cell_type": "A"})
    obs = pd.DataFrame(rows)
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    X = rng.negative_binomial(20, 0.5, (len(obs), 40)).astype(np.float32)
    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=[f"GENE_{i:04d}" for i in range(40)]),
    )
    input_path = tmp_path / "pb.h5ad"
    adata.write_h5ad(input_path)

    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de
    result = rna_pseudobulk_de({
        "data_path":     str(input_path),
        "groupby":       "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons":   [["treat", "ctrl"]],
        "output_dir":    str(tmp_path / "out"),
        # min_replicates_per_condition NOT set — default must be 3
    })
    assert result["status"] == "success"
    # With default min_replicates=3 and only n=2 per condition, the group
    # is skipped wholesale before per_comparison even runs (the script
    # gates on len(pseudosamples) < 2 * min_replicates_per_condition).
    skipped_any = False
    for group_info in result["per_group"].values():
        group_status = group_info.get("status")
        if group_status == "skipped":
            skipped_any = True
            assert "pseudosamples" in group_info.get("reason", "") \
                or "replicates" in group_info.get("reason", "")
            continue
        for comp in (group_info.get("per_comparison") or {}).values():
            if comp.get("status") == "skipped":
                skipped_any = True
                assert "not enough replicates" in comp.get("reason", "")
    assert skipped_any, (
        "Default min_replicates=3 must reject n=2 designs somewhere: "
        f"{result['per_group']}"
    )


def test_design_intelligence_downgrades_pseudobulk_at_n2():
    """T1.3 — Design Intelligence must move pseudobulk DE from
    'recommended' to 'optional' (with a caveat) when any group has
    fewer than three biological replicates."""
    from aria.agents.design_intelligence import DesignIntelligence

    di = DesignIntelligence()

    exp_context_n3 = {
        "design": {
            "main_factor": "condition",
            "groups": {"ctrl": ["c1", "c2", "c3"], "treat": ["t1", "t2", "t3"]},
            "pseudobulk": {
                "condition_col": "condition",
                "replicate_col": "donor",
                "groupby_col":   "cell_type",
            },
        },
        "modalities": {"scRNA": []},
    }
    exp_context_n2 = {
        "design": {
            "main_factor": "condition",
            "groups": {"ctrl": ["c1", "c2"], "treat": ["t1", "t2"]},
            "pseudobulk": {
                "condition_col": "condition",
                "replicate_col": "donor",
                "groupby_col":   "cell_type",
            },
        },
        "modalities": {"scRNA": []},
    }

    profile_n3 = di._scrna_profile(exp_context_n3, {})
    profile_n2 = di._scrna_profile(exp_context_n2, {})

    rec_n3 = " ".join(profile_n3["recommended"])
    assert "pseudobulk DE" in rec_n3, rec_n3

    rec_n2 = " ".join(profile_n2["recommended"])
    opt_n2 = " ".join(profile_n2["optional"])
    warn_n2 = " ".join(profile_n2["warnings"])
    assert "pseudobulk DE" not in rec_n2, (
        "Pseudobulk DE must not be 'recommended' at n=2: " + rec_n2
    )
    assert "low-power caveat" in opt_n2, opt_n2
    assert "n=2" in warn_n2 or "replicates" in warn_n2, warn_n2


def test_design_intelligence_downgrades_bulk_at_n2():
    """T1.3 — same downgrade rule applies to bulk RNA DE."""
    from aria.agents.design_intelligence import DesignIntelligence

    di = DesignIntelligence()
    ctx = {"design": {"groups": {"ctrl": ["c1", "c2"], "treat": ["t1", "t2"]}}}
    profile = di._bulk_rna_profile("bulk_RNA", ctx, {})

    rec = " ".join(profile["recommended"])
    opt = " ".join(profile["optional"])
    warn = " ".join(profile["warnings"])
    assert "DESeq2 differential expression with biological replicates" not in rec
    assert "low-power caveat" in opt, opt
    assert "n=2" in warn or "replicates" in warn, warn


def test_scrna_pseudobulk_gate_blocks_n2_in_recommended_only():
    """T1.3 — CP2 recommended-only must not run n=2 pseudobulk DE.
    Design Intelligence marks it optional, so it only runs when optional
    supported analyses were selected."""
    from aria.agents.scrna_agent import scRNAAgent

    agent = scRNAAgent.__new__(scRNAAgent)
    intent = {"summary": "compare treated vs control single-cell RNA-seq"}
    exp_ctx = {
        "run_optional_supported": False,
        "design": {
            "groups": {"control": ["c1", "c2"], "treated": ["t1", "t2"]},
            "pseudobulk": {
                "condition_col": "condition",
                "replicate_col": "donor",
                "groupby_col": "cell_type",
            },
        },
    }

    assert agent._needs_pseudobulk(intent, exp_ctx) is False


def test_scrna_pseudobulk_gate_allows_n2_when_optional_selected():
    """T1.3 — n=2 pseudobulk remains supported with warnings when CP2
    selected recommended + optional supported analyses."""
    from aria.agents.scrna_agent import scRNAAgent

    agent = scRNAAgent.__new__(scRNAAgent)
    intent = {"summary": "compare treated vs control single-cell RNA-seq"}
    exp_ctx = {
        "run_optional_supported": True,
        "design": {
            "groups": {"control": ["c1", "c2"], "treated": ["t1", "t2"]},
            "pseudobulk": {
                "condition_col": "condition",
                "replicate_col": "donor",
                "groupby_col": "cell_type",
            },
        },
    }

    assert agent._needs_pseudobulk(intent, exp_ctx) is True


# ── T1.1: differential abundance + composition correction ─────────────────────

def test_diff_abundance_detects_2x_shift(tmp_path):
    """T1.1 — when a cell type's abundance doubles in one condition, the
    Poisson-offset GLM (or Fisher-exact fallback) must flag it as
    significant and assign direction='up'."""
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    rng = np.random.default_rng(42)
    rows = []
    # 4 replicates per condition × 3 cell types. Cell type "B" doubles in
    # condition "treat": 50 -> 100 cells. Types A and C stay flat.
    base_counts = {"A": 100, "B": 50, "C": 80}
    for cond in ("ctrl", "treat"):
        for rep_i in range(1, 5):
            for ct in ("A", "B", "C"):
                n = base_counts[ct]
                if ct == "B" and cond == "treat":
                    n = base_counts[ct] * 2
                # Small per-replicate jitter so the GLM has variance to fit.
                jitter = int(rng.normal(0, max(1, n * 0.05)))
                for _ in range(max(1, n + jitter)):
                    rows.append({
                        "condition": cond,
                        "donor":     f"{cond}_{rep_i}",
                        "cell_type": ct,
                    })
    obs = pd.DataFrame(rows)
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    X = np.zeros((len(obs), 5), dtype=np.float32)
    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=[f"GENE_{i}" for i in range(5)]),
    )
    h5ad_path = tmp_path / "da.h5ad"
    adata.write_h5ad(h5ad_path)

    from aria.scripts.rna_diff_abundance import rna_diff_abundance
    result = rna_diff_abundance({
        "data_path":     str(h5ad_path),
        "groupby":       "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons":   [["treat", "ctrl"]],
        "output_dir":    str(tmp_path / "out"),
    })

    assert result["status"] == "success", result
    assert result["any_significant"] is True
    comp = result["per_comparison"]["treat_vs_ctrl"]
    assert comp["status"] == "success"
    # B doubled in absolute terms in treat -> significantly UP in proportion.
    b_row = next(r for r in comp["per_cell_type"] if r["name"] == "B")
    assert b_row["significant"] is True, b_row
    assert b_row["direction"] == "up", b_row
    assert b_row["log2_fold_change"] > 0.5, b_row
    # A and C are flat in ABSOLUTE counts, but because B doubled the total
    # number of cells grew, so A and C drop in PROPORTION. That is the
    # whole point of composition analysis and the GLM correctly flags it
    # as a downward shift. The narrative makes this visible.
    a_row = next(r for r in comp["per_cell_type"] if r["name"] == "A")
    c_row = next(r for r in comp["per_cell_type"] if r["name"] == "C")
    assert a_row["log2_fold_change"] < 0, a_row
    assert c_row["log2_fold_change"] < 0, c_row


def test_diff_abundance_no_signal_returns_none_significant(tmp_path):
    """T1.1 — flat data must yield no significant shifts and the agent's
    decision logic relies on this to keep composition_covariate OFF."""
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    rng = np.random.default_rng(0)
    rows = []
    for cond in ("ctrl", "treat"):
        for rep_i in range(1, 5):
            for ct in ("A", "B"):
                # ~80 cells of each type in each replicate, identical between
                # conditions modulo small jitter.
                n = 80 + int(rng.normal(0, 3))
                for _ in range(max(1, n)):
                    rows.append({
                        "condition": cond,
                        "donor":     f"{cond}_{rep_i}",
                        "cell_type": ct,
                    })
    obs = pd.DataFrame(rows)
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    adata = ad.AnnData(
        X=np.zeros((len(obs), 3), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=[f"GENE_{i}" for i in range(3)]),
    )
    h5ad_path = tmp_path / "da_flat.h5ad"
    adata.write_h5ad(h5ad_path)

    from aria.scripts.rna_diff_abundance import rna_diff_abundance
    result = rna_diff_abundance({
        "data_path":     str(h5ad_path),
        "groupby":       "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons":   [["treat", "ctrl"]],
        "output_dir":    str(tmp_path / "out"),
    })
    assert result["status"] == "success"
    assert result["any_significant"] is False
    comp = result["per_comparison"]["treat_vs_ctrl"]
    assert comp["n_significant"] == 0


def test_pseudobulk_de_passes_composition_covariate_when_da_significant(tmp_path):
    """T1.1 — scrna_agent._run_pseudobulk must pass
    composition_covariate=True to rna_pseudobulk_de whenever the DA
    pre-pass returns any_significant=True."""
    from aria.agents.scrna_agent import scRNAAgent

    h5ad_path = tmp_path / "annotated.h5ad"
    h5ad_path.write_text("placeholder")
    calls = []

    class FakeEnv:
        def run_in_stack(self, **kwargs):
            calls.append(kwargs)
            script = kwargs["script_path"]
            if script.endswith("rna_diff_abundance.py"):
                return {
                    "status": "success",
                    "method": "poisson_offset_glm",
                    "groupby": kwargs["params"]["groupby"],
                    "condition_col": kwargs["params"]["condition_col"],
                    "replicate_col": kwargs["params"]["replicate_col"],
                    "covariates": kwargs["params"].get("covariates", []),
                    "significance_alpha": 0.10,
                    "any_significant": True,
                    "per_comparison": {
                        "old_vs_young": {
                            "status": "success",
                            "n_significant": 1,
                            "n_replicates": {"test": 2, "ref": 2},
                            "per_cell_type": [{
                                "name": "Microglia",
                                "n_test": 200, "n_ref": 100,
                                "prop_test": 0.2, "prop_ref": 0.1,
                                "log2_fold_change": 1.0,
                                "pval": 0.001, "padj": 0.003,
                                "direction": "up", "significant": True,
                            }],
                        }
                    },
                    "n_replicates_per_group": {"old": 2, "young": 2},
                    "warnings": [],
                }
            if script.endswith("rna_pseudobulk_de.py"):
                return {
                    "status": "success",
                    "groupby": kwargs["params"]["groupby"],
                    "condition_col": kwargs["params"]["condition_col"],
                    "replicate_col": kwargs["params"]["replicate_col"],
                    "covariates": kwargs["params"]["covariates"],
                    "thresholds": {},
                    "n_groups": 0,
                    "per_group": {},
                }
            raise AssertionError(f"unexpected script: {script}")

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.env = FakeEnv()
    agent.publish_finding = lambda *args, **kwargs: None
    agent._inject_condition_obs = lambda *args, **kwargs: (_ for _ in ()).throw(
        AssertionError("condition injection should be skipped")
    )

    result = agent._run_pseudobulk(
        "exp",
        str(h5ad_path),
        {
            "organism": "Homo sapiens",
            "design": {
                "groups": {"old": ["o1", "o2"], "young": ["y1", "y2"]},
                "main_factor": "age_group",
                "pseudobulk": {
                    "from_obs": True,
                    "condition_col": "age_group",
                    "replicate_col": "donor_id",
                    "groupby_col": "subclass",
                    "comparisons": [["old", "young"]],
                },
            },
        },
        {"summary": "compare old vs young"},
        {},
    )

    assert result["status"] == "success"
    pb_call = next(
        c for c in calls if c["script_path"].endswith("rna_pseudobulk_de.py")
    )
    assert pb_call["params"]["composition_covariate"] is True, (
        "DA flagged a significant shift; pseudobulk MUST receive "
        "composition_covariate=True"
    )
    assert pb_call["params"]["min_replicates_per_condition"] == 3, (
        "Recommended-only scRNA pseudobulk must keep the publication-ready "
        "n>=3 replicate floor; n=2 is optional-supported only."
    )
    # Findings must surface the differential_abundance block back to the agent.
    assert result.get("differential_abundance") is not None
    assert result["differential_abundance"]["any_significant"] is True


def test_pseudobulk_composition_covariate_flag_propagates_into_blocks(tmp_path):
    """T1.1 — when called with composition_covariate=True and the
    covariate actually varies across pseudosamples, each successful
    per_comparison block must carry corrected_for_composition=True."""
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")

    rng = np.random.default_rng(11)
    # 3 donors per condition, 2 cell types. We give donors very different
    # totals so the composition log-ratio varies across pseudosamples.
    donor_totals = {
        "ctrl_1": (40, 200), "ctrl_2": (30, 250), "ctrl_3": (50, 220),
        "treat_1": (90, 220), "treat_2": (110, 200), "treat_3": (100, 180),
    }
    rows = []
    for donor, (n_a, n_b) in donor_totals.items():
        cond = "treat" if donor.startswith("treat") else "ctrl"
        for _ in range(n_a):
            rows.append({"condition": cond, "donor": donor, "cell_type": "A"})
        for _ in range(n_b):
            rows.append({"condition": cond, "donor": donor, "cell_type": "B"})
    obs = pd.DataFrame(rows)
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    n_genes = 40
    X = rng.negative_binomial(20, 0.5, (len(obs), n_genes)).astype(np.float32)
    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=[f"GENE_{i:04d}" for i in range(n_genes)]),
    )
    h5ad_path = tmp_path / "pb_comp.h5ad"
    adata.write_h5ad(h5ad_path)

    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de
    result = rna_pseudobulk_de({
        "data_path":     str(h5ad_path),
        "groupby":       "cell_type",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons":   [["treat", "ctrl"]],
        "output_dir":    str(tmp_path / "out"),
        "composition_covariate": True,
        "min_cells_per_pseudosample": 5,
        "padj_max":      0.5,
        "lfc_min":       0.0,
    })

    assert result["status"] == "success"
    assert result["composition_covariate_requested"] is True
    mt = result.get("multiple_testing") or {}
    assert mt.get("local_method") == "BH"
    assert mt.get("global_method") == "BH"
    assert mt.get("n_tests_global", 0) > 0
    found_corrected = False
    for group_info in result["per_group"].values():
        for comp in (group_info.get("per_comparison") or {}).values():
            if comp.get("status") != "success":
                continue
            assert comp["corrected_for_composition"] is True, comp
            assert isinstance(comp.get("power_estimate_at_lfc_min"), float)
            assert 0.0 <= comp["power_estimate_at_lfc_min"] <= 1.0
            assert "n_significant_local" in comp
            assert "n_significant_global" in comp
            for rec in comp.get("top_genes", []) + comp.get("all_sig", []):
                assert "padj_local" in rec
                assert "padj_global" in rec
                assert rec["padj"] == rec["padj_global"]
            found_corrected = True
    assert found_corrected, "Expected at least one composition-corrected block"


def test_global_fdr_is_more_conservative_than_local_family():
    """T1.2 — global BH across all block-gene tests can only reduce the
    number of calls for a fixed alpha in this constructed family."""
    import numpy as np
    from aria.scripts.rna_pseudobulk_de import _bh_correct, _global_bh

    blocks = []
    for block_i in range(5):
        pvals = np.ones(1000, dtype=float)
        if block_i == 0:
            pvals[0] = 1e-7
        blocks.append(pvals)

    alpha = 2e-4
    n_local = sum(int((_bh_correct(p) < alpha).sum()) for p in blocks)
    global_padj = _global_bh(np.concatenate(blocks))
    n_global = int((global_padj < alpha).sum())

    assert n_global <= n_local
    assert n_local == 1
    assert n_global == 0


def test_pathway_ora_receives_explicit_background(tmp_path):
    """T1.4 — scRNA ORA receives the dataset-expressed gene universe, not
    only the significant DE list."""
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    from aria.agents.scrna_agent import scRNAAgent

    adata = ad.AnnData(
        X=np.array([[1, 0, 2], [0, 0, 3]], dtype=np.float32),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    h5ad_path = tmp_path / "expr.h5ad"
    adata.write_h5ad(h5ad_path)
    calls = []

    class FakeEnv:
        def run_in_stack(self, **kwargs):
            calls.append(kwargs)
            return {
                "status": "success",
                "organism": "Homo sapiens",
                "databases": {},
                "background_size": len(kwargs["params"]["background_genes"]),
                "background_source": "dataset_expressed_genes",
                "per_cluster": {},
                "output_csv": None,
            }

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.env = FakeEnv()
    agent.publish_finding = lambda *args, **kwargs: None
    de_result = {
        "de_genes_by_cluster": {
            "0": [
                {"gene": "G1", "log2fc": 1.0, "padj": 0.01},
                {"gene": "G2", "log2fc": 1.0, "padj": 0.01},
                {"gene": "G3", "log2fc": 1.0, "padj": 0.01},
            ]
        }
    }

    result = agent._run_pathway_per_cluster(
        "exp", de_result, {"organism": "Homo sapiens"}, data_path=str(h5ad_path)
    )

    assert result["status"] == "success"
    params = calls[0]["params"]
    assert set(params["background_genes"]) == {"G1", "G3"}


def test_power_estimate_monotone_in_n():
    """T1.5 — power should increase with replicate count when dispersion
    and effect size are held fixed."""
    from aria.utils.power_estimation import pseudobulk_power_estimate

    common = {
        "dispersion_estimate": 0.12,
        "target_log2fc": 0.5,
        "alpha": 0.05,
        "mean_expression": 100.0,
    }
    p2 = pseudobulk_power_estimate(2, **common)
    p3 = pseudobulk_power_estimate(3, **common)
    p10 = pseudobulk_power_estimate(10, **common)

    assert p10 > p3 > p2


def test_provenance_block_contains_required_fields():
    from aria.utils.provenance import collect_provenance

    prov = collect_provenance()
    for key in {
        "aria_version", "git_sha", "git_dirty", "python_version",
        "platform", "conda_env", "timestamp_utc",
    }:
        assert key in prov


def test_hash_params_order_invariant():
    from aria.utils.provenance import hash_params

    assert hash_params({"a": 1, "b": {"c": 2}}) == \
        hash_params({"b": {"c": 2}, "a": 1})


def test_hash_file_stable_under_chunk_size(tmp_path):
    from aria.utils.provenance import hash_file

    p = tmp_path / "blob.bin"
    p.write_bytes(b"abc123" * 1000)
    assert hash_file(p, chunk_size=7) == hash_file(p, chunk_size=8192)


def test_lockfile_embed_renders_warning_when_missing():
    from aria.agents.narrative_agent import NarrativeAgent

    html = NarrativeAgent._build_lockfile_section()
    assert "conda lockfiles" in html.lower() or ".linux-64.lock" in html


def test_methodology_json_emitted_with_required_keys(tmp_path):
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    agent.memory = type("M", (), {"db_path": ":memory:"})()
    report = agent._render_html_report(
        experiment_id="exp_test",
        exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
        intent={"summary": "test"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={
            "high": [], "medium": [], "low": [], "insufficient": []
        },
        methods="methods",
        decisions=[],
        agent_results={},
        report_dir=tmp_path / "report2",
    )
    methodology = report.parent / "methodology.json"
    assert methodology.exists()
    data = json.loads(methodology.read_text())
    for key in {
        "provenance", "design", "design_intelligence", "thresholds",
        "seeds", "tools", "decisions",
    }:
        assert key in data


def test_reproducible_mode_writes_memory_snapshot(tmp_path):
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    db = tmp_path / "memory.db"
    db.write_bytes(b"sqlite")
    agent.memory = type("M", (), {"db_path": str(db)})()
    report = agent._render_html_report(
        experiment_id="exp_test",
        exp_ctx={
            "organism": "Homo sapiens",
            "genome": "hg38",
            "reproducible_mode": True,
        },
        intent={"summary": "test"},
        executive_summary="ok",
        findings_sections={"conflicts": "none"},
        grouped_findings={
            "high": [], "medium": [], "low": [], "insufficient": []
        },
        methods="methods",
        decisions=[],
        agent_results={},
        report_dir=tmp_path / "report3",
    )
    html = report.read_text()
    assert "&lt;timestamp redacted for byte-identity&gt;" in html
    assert (report.parent / "memory_snapshot.sqlite").exists()


def test_reproducible_mode_byte_identity(tmp_path):
    from aria.agents.narrative_agent import NarrativeAgent

    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.reports_dir = tmp_path
    agent.memory = type("M", (), {"db_path": ":memory:"})()
    common = {
        "experiment_id": "exp_test",
        "exp_ctx": {
            "organism": "Homo sapiens",
            "genome": "hg38",
            "reproducible_mode": True,
            "provenance": {
                "aria_version": "test",
                "git_sha": "abc",
                "git_dirty": False,
                "python_version": "py",
                "platform": "linux",
                "conda_env": "",
                "timestamp_utc": "will be redacted",
            },
        },
        "intent": {"summary": "test"},
        "executive_summary": "ok",
        "findings_sections": {"conflicts": "none"},
        "grouped_findings": {
            "high": [], "medium": [], "low": [], "insufficient": []
        },
        "methods": "methods",
        "decisions": [],
        "agent_results": {},
    }
    r1 = agent._render_html_report(report_dir=tmp_path / "r1", **common)
    r2 = agent._render_html_report(report_dir=tmp_path / "r2", **common)

    assert r1.read_bytes() == r2.read_bytes()
