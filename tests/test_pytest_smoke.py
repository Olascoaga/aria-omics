from __future__ import annotations

import os
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
    assert inferred["pseudobulk"]["from_obs"] is True


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
    params = calls[0]["params"]
    assert params["data_path"] == str(h5ad_path)
    assert params["condition_col"] == "age_group"
    assert params["replicate_col"] == "donor_id"
    assert params["groupby"] == "subclass"
    assert params["covariates"] == ["Gender"]


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
