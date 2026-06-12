"""A-MARK1: per-cluster marker calls must honor pct_in filtering."""

from __future__ import annotations


def test_de_per_cluster_filters_low_pct_in_markers(tmp_path, monkeypatch):
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
            self._cols = {"leiden": _Series(["0"] * 10 + ["1"] * 10)}

        def __getitem__(self, key):
            return self._cols[key]

    class _Raw:
        var_names = ["rare_marker", "common_marker"]

        def __init__(self):
            x = np.zeros((20, 2), dtype=float)
            x[0, 0] = 1.0          # rare_marker pct_in = 0.10 in cluster 0
            x[:5, 1] = 1.0         # common_marker pct_in = 0.50 in cluster 0
            self.X = x

    class _Adata:
        def __init__(self):
            self.obs = _Obs()
            self.raw = _Raw()
            self.X = self.raw.X
            self.var_names = self.raw.var_names
            self.uns = {}

    adata = _Adata()

    def _rank_genes_groups(adata_arg, **_kwargs):
        dtype = [("0", "O"), ("1", "O")]
        adata_arg.uns["aria_de"] = {
            "names": np.array(
                [("rare_marker", "common_marker"),
                 ("common_marker", "rare_marker")],
                dtype=dtype,
            ),
            "logfoldchanges": np.array(
                [(3.0, 3.0), (3.0, 3.0)],
                dtype=dtype,
            ),
            "pvals_adj": np.array(
                [(0.001, 0.001), (0.001, 0.001)],
                dtype=dtype,
            ),
        }

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
    monkeypatch.setattr(safe_h5ad, "read_h5ad", lambda _path: adata)

    input_path = tmp_path / "clustered.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    result = mod.rna_de_per_cluster(
        {
            "data_path": str(input_path),
            "groupby": "leiden",
            "padj_max": 0.05,
            "lfc_min": 0.5,
            "min_pct_in": 0.2,
            "output_dir": str(tmp_path),
        }
    )

    assert result["status"] == "success"
    assert result["min_pct_in"] == 0.2
    cluster0 = result["de_genes_by_cluster"]["0"]
    assert [r["gene"] for r in cluster0] == ["common_marker"]
    assert cluster0[0]["pct_in"] == 0.5
    assert result["n_sig_by_cluster"]["0"] == 1
