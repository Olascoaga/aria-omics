"""A-CMT1/A-CLUST1 guards for marker comments and duplicate computation."""

from __future__ import annotations


def test_marker_comments_do_not_call_lognorm_raw_counts():
    from pathlib import Path

    clustering_src = Path("aria/scripts/rna_clustering.py").read_text()
    de_src = Path("aria/scripts/rna_de_per_cluster.py").read_text()

    assert "use raw counts if available" not in clustering_src.lower()
    assert "use raw counts if present" not in de_src.lower()
    assert "log-normalized" in clustering_src.lower()
    assert "log-normalized" in de_src.lower()


def test_de_per_cluster_reuses_existing_rank_genes_groups(tmp_path, monkeypatch):
    import sys
    import types

    import numpy as np

    import aria.utils.safe_h5ad as safe_h5ad
    from aria.scripts import rna_de_per_cluster as mod

    class _Series:
        def __init__(self, values):
            self._values = np.asarray(values)

        def nunique(self):
            return len(set(map(str, self._values)))

        def unique(self):
            return np.unique(self._values)

        def astype(self, _kind):
            return _Series([str(v) for v in self._values])

        def __eq__(self, other):
            return self._values == other

    class _Obs:
        columns = ["leiden"]

        def __init__(self):
            self._cols = {"leiden": _Series(["0"] * 4 + ["1"] * 4)}

        def __getitem__(self, key):
            return self._cols[key]

    class _Raw:
        var_names = ["marker_a", "marker_b"]

        def __init__(self):
            self.X = np.ones((8, 2), dtype=float)

    dtype = [("0", "O"), ("1", "O")]

    class _Adata:
        def __init__(self):
            self.obs = _Obs()
            self.raw = _Raw()
            self.X = self.raw.X
            self.var_names = self.raw.var_names
            self.uns = {
                "rank_genes_groups": {
                    "names": np.array(
                        [("marker_a", "marker_b"), ("marker_b", "marker_a")],
                        dtype=dtype,
                    ),
                    "logfoldchanges": np.array([(2.0, 2.0), (2.0, 2.0)], dtype=dtype),
                    "pvals_adj": np.array([(0.001, 0.001), (0.001, 0.001)], dtype=dtype),
                }
            }

    def _rank_genes_groups(*_args, **_kwargs):
        raise AssertionError("existing clustering markers should be reused")

    fake_scanpy = types.SimpleNamespace(
        tl=types.SimpleNamespace(rank_genes_groups=_rank_genes_groups)
    )

    class _FakeDataFrame:
        def __init__(self, rows):
            self.rows = rows

        def to_csv(self, path, index=False):
            return None

    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setitem(sys.modules, "pandas", types.SimpleNamespace(DataFrame=_FakeDataFrame))
    monkeypatch.setattr(safe_h5ad, "read_h5ad", lambda _path: _Adata())

    input_path = tmp_path / "clustered.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    result = mod.rna_de_per_cluster(
        {
            "data_path": str(input_path),
            "groupby": "leiden",
            "output_dir": str(tmp_path),
        }
    )

    assert result["status"] == "success"
    assert result["marker_source"] == "rank_genes_groups"
    assert result["n_significant_total"] == 4
