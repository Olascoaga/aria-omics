"""B-TRAJ1 guards for DPT root selection."""

from __future__ import annotations


# ── T13 (tri-audit 2026-06-14): progenitor root tokens must match WHOLE WORDS, so
# "stem" does not fire on "brainstem". ──────────────────────────────────────────

def test_progenitor_label_matches_whole_words_only():
    from aria.scripts.rna_trajectory import _is_progenitor_label

    # Genuine progenitor-like populations (whole word, incl. plural).
    assert _is_progenitor_label("Neural stem cells")
    assert _is_progenitor_label("Radial glia progenitors")
    assert _is_progenitor_label("Oligodendrocyte precursor cells")
    assert _is_progenitor_label("stem")
    # Substrings must NOT match.
    assert not _is_progenitor_label("Brainstem neurons")
    assert not _is_progenitor_label("Astrocytes")
    assert not _is_progenitor_label("Microglia")


# ── T14 (tri-audit 2026-06-14): the PAGA "strong" connectivity threshold is 0.05
# (the code), not the 0.005 the stale comment claimed. Lock the operative value. ──

def test_paga_strong_connectivity_threshold_is_0_05():
    from aria.scripts.rna_trajectory import _PAGA_STRONG_CONNECTIVITY

    assert _PAGA_STRONG_CONNECTIVITY == 0.05
    # A 0.02 edge is weak (below 0.05); were the threshold wrongly lowered to the
    # commented 0.005 it would be counted as strong.
    assert not (0.02 >= _PAGA_STRONG_CONNECTIVITY)
    assert 0.2 >= _PAGA_STRONG_CONNECTIVITY


def test_trajectory_skips_dpt_without_biological_root(tmp_path, monkeypatch):
    import sys
    import types

    import numpy as np

    import aria.utils.safe_h5ad as safe_h5ad
    from aria.scripts import rna_trajectory as mod

    class _Series:
        def __init__(self, values):
            self.values = np.asarray(values)

        def unique(self):
            return np.unique(self.values)

    class _Obs:
        columns = ["leiden", "n_genes_by_counts"]

        def __init__(self):
            self._cols = {
                "leiden": _Series(["0", "0", "1", "1"]),
                "n_genes_by_counts": _Series([20, 2, 25, 30]),
            }

        def __getitem__(self, key):
            return self._cols[key]

    class _Adata:
        def __init__(self):
            self.obs = _Obs()
            self.uns = {}
            self.obsm = {"X_pca": np.ones((4, 2), dtype=float)}
            self.layers = {}

        def write_h5ad(self, path):
            self.output_path = path

    def _paga(adata, groups):
        adata.uns["paga"] = {
            "connectivities": np.array([[0.0, 0.2], [0.2, 0.0]], dtype=float)
        }

    def _dpt(_adata):
        raise AssertionError("DPT must not run without an explicit or progenitor root")

    fake_scanpy = types.SimpleNamespace(
        pp=types.SimpleNamespace(neighbors=lambda *_args, **_kwargs: None),
        tl=types.SimpleNamespace(
            paga=_paga,
            umap=lambda *_args, **_kwargs: None,
            diffmap=lambda *_args, **_kwargs: None,
            dpt=_dpt,
        ),
    )
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setattr(safe_h5ad, "read_h5ad", lambda _path: _Adata())

    input_path = tmp_path / "clustered.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")

    result = mod.rna_trajectory(
        {
            "data_path": str(input_path),
            "cell_type_col": "cell_type",
            "output_dir": str(tmp_path),
        }
    )

    assert result["status"] == "success"
    assert result["paga"]["n_connections"] == 1
    assert result["pseudotime"] == {
        "computed": False,
        "reason": "root_unresolved",
        "root_used": None,
        "root_selection_policy": "explicit_or_progenitor_label_required",
    }
