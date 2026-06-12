"""B-PB1: production scRNA pseudobulk must use raw QC counts, not lognorm raw."""

from __future__ import annotations


def test_scrna_agent_dispatches_raw_qc_counts_to_pseudobulk(tmp_path, monkeypatch):
    import sys
    import types

    monkeypatch.setitem(
        sys.modules,
        "litellm",
        types.SimpleNamespace(completion=lambda *args, **kwargs: None),
    )

    from aria.agents.scrna_agent import scRNAAgent

    annotated_h5ad = tmp_path / "annotated.h5ad"
    qc_h5ad = tmp_path / "qc_filtered.h5ad"
    annotated_h5ad.write_text("placeholder", encoding="utf-8")
    qc_h5ad.write_text("placeholder", encoding="utf-8")
    calls = []

    class FakeEnv:
        def run_in_stack(self, **kwargs):
            calls.append(kwargs)
            script = kwargs["script_path"]
            if script.endswith("rna_diff_abundance.py"):
                return {"status": "success", "any_significant": False}
            if script.endswith("rna_pseudobulk_de.py"):
                return {
                    "status": "success",
                    "groupby": kwargs["params"]["groupby"],
                    "condition_col": kwargs["params"]["condition_col"],
                    "replicate_col": kwargs["params"]["replicate_col"],
                    "covariates": kwargs["params"].get("covariates", []),
                    "thresholds": {},
                    "n_groups": 0,
                    "per_group": {},
                    "count_source": "raw_counts",
                    "count_source_data_path": kwargs["params"].get("counts_data_path"),
                    "lognorm_recovered": False,
                }
            raise AssertionError(f"unexpected script: {script}")

    agent = scRNAAgent.__new__(scRNAAgent)
    agent.env = FakeEnv()
    agent._workspace = lambda *a, **k: tmp_path
    agent._log_decision = lambda *a, **k: None
    agent.publish_finding = lambda *a, **k: None

    result = agent._run_pseudobulk(
        "exp",
        str(annotated_h5ad),
        {
            "design": {
                "groups": {"treated": ["t1", "t2", "t3"], "ctrl": ["c1", "c2", "c3"]},
                "main_factor": "condition",
                "pseudobulk": {
                    "from_obs": True,
                    "condition_col": "condition",
                    "replicate_col": "donor",
                    "groupby_col": "cell_type",
                    "comparisons": [["treated", "ctrl"]],
                },
            },
        },
        {},
        {},
        raw_counts_h5ad=str(qc_h5ad),
    )

    assert result["status"] == "success"
    assert result["pseudobulk_de"]["count_source"] == "raw_counts"
    assert result["pseudobulk_de"]["lognorm_recovered"] is False
    assert result["pseudobulk_de"]["count_source_data_path"] == str(qc_h5ad)
    pb_call = next(
        c for c in calls if c["script_path"].endswith("rna_pseudobulk_de.py")
    )
    assert pb_call["params"]["data_path"] == str(annotated_h5ad)
    assert pb_call["params"]["counts_data_path"] == str(qc_h5ad)


def test_pseudobulk_uses_counts_data_path_instead_of_lognorm_raw(tmp_path, monkeypatch):
    pytest = __import__("pytest")
    pytest.importorskip("scipy.sparse")
    np = __import__("pytest").importorskip("numpy")
    pd = pytest.importorskip("pandas")
    import aria.utils.safe_h5ad as safe_h5ad

    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de

    obs_names = [f"cell{i}" for i in range(12)]
    var_names = [f"Gene{i}" for i in range(6)]
    raw_counts = np.array(
        [[5 + (i % 3), 2, 0, 1, 3, 4] for i in range(len(obs_names))],
        dtype=np.int64,
    )
    obs = pd.DataFrame(
        {
            "cell_type": ["TypeA"] * 12,
            "condition": ["treated"] * 6 + ["ctrl"] * 6,
            "donor": [f"t{i // 2}" for i in range(6)] + [f"c{i // 2}" for i in range(6)],
            "total_counts": raw_counts.sum(axis=1),
        },
        index=obs_names,
    )
    counts_path = tmp_path / "qc_filtered.h5ad"
    counts_path.write_text("placeholder", encoding="utf-8")

    lib = raw_counts.sum(axis=1, keepdims=True)
    lognorm = np.log1p(raw_counts / lib * 10000.0)
    labelled_path = tmp_path / "annotated.h5ad"
    labelled_path.write_text("placeholder", encoding="utf-8")

    class _Raw:
        def __init__(self, X, var_names):
            self.X = X
            self.var_names = pd.Index(var_names)

    class _Adata:
        def __init__(self, X, obs, var_names, raw=None):
            self.X = X
            self.obs = obs
            self.var_names = pd.Index(var_names)
            self.obs_names = obs.index
            self.raw = raw
            self.n_obs = len(obs)

        def __getitem__(self, item):
            rows, cols = item
            if not isinstance(rows, slice):
                row_index = self.obs.index[rows]
                obs = self.obs.loc[row_index].copy()
                X = self.X[rows]
            else:
                obs = self.obs.copy()
                X = self.X
            if not isinstance(cols, slice):
                X = X[:, cols]
                var_names_out = self.var_names[cols]
            else:
                var_names_out = self.var_names
            return _Adata(X, obs, list(var_names_out), self.raw)

        def copy(self):
            return _Adata(self.X.copy(), self.obs.copy(), list(self.var_names), self.raw)

    labelled = _Adata(
        lognorm[:, :3],
        obs.copy(),
        var_names[:3],
        raw=_Raw(lognorm, var_names),
    )
    counts = _Adata(raw_counts, obs.copy(), var_names)

    def _fake_read_h5ad(path):
        if str(path) == str(labelled_path):
            return labelled
        if str(path) == str(counts_path):
            return counts
        raise AssertionError(f"unexpected h5ad path: {path}")

    monkeypatch.setattr(safe_h5ad, "read_h5ad", _fake_read_h5ad)
    result = rna_pseudobulk_de(
        {
            "data_path": str(labelled_path),
            "counts_data_path": str(counts_path),
            "groupby": "cell_type",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["treated", "ctrl"]],
            "min_cells_per_pseudosample": 1,
            # Keep this light and pydeseq2-free: count-source resolution happens
            # before the per-comparison replicate gate.
            "min_replicates_per_condition": 99,
            "output_dir": str(tmp_path / "out"),
        }
    )

    assert result["status"] == "success"
    assert result["count_source"] == "raw_counts"
    assert result["lognorm_recovered"] is False
    assert result["count_source_data_path"] == str(counts_path)
