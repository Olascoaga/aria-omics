"""B-QC1/B-QC2 + N-QC1 guards for scRNA QC threshold policy."""

from __future__ import annotations


def test_mad_bounds_log_space_is_less_aggressive_on_skewed_upper_tail():
    """N-QC1: total_counts / n_genes_by_counts are right-skewed (log-normal),
    so an upper MAD bound on RAW values clips the legitimate high tail. The
    sc-best-practices approach computes MAD in log space."""
    from aria.scripts.rna_qc import _mad_bounds

    counts = [1000, 1200, 1500, 1800, 2000, 2500, 3000, 4000, 6000, 10000]

    _, hi_raw = _mad_bounds(counts)
    _, hi_log = _mad_bounds(counts, log_space=True)

    # Log-space upper bound is less aggressive on the skewed right tail.
    assert hi_log > hi_raw
    # A legitimately high-count cell clipped by raw MAD survives in log space.
    assert 6000 > hi_raw       # raw MAD would discard it
    assert 6000 <= hi_log      # log-space MAD keeps it


def test_mad_bounds_log_space_keeps_lower_bound_nonnegative():
    from aria.scripts.rna_qc import _mad_bounds

    counts = [500, 800, 1000, 1200, 1500, 2000, 3000, 5000]
    lo, hi = _mad_bounds(counts, log_space=True)
    assert lo >= 0.0
    assert hi > lo


def test_mad_bounds_zero_mad_does_not_collapse_to_median():
    """T6: when most cells share an identical value MAD=0, so the bounds collapse
    to the median (zero width) and clip legitimate minimal variation. A minimum-
    width fallback must keep a near-median cell while still removing a gross
    outlier."""
    from aria.scripts.rna_qc import _mad_bounds

    # 10 cells at 2000 + 1 at 2050 + 1 gross outlier -> MAD over abs-devs is 0.
    counts = [2000] * 10 + [2050, 200000]
    lo, hi = _mad_bounds(counts, log_space=True)

    # Zero-width bug would set hi == 2000 and discard the 2050 cell.
    assert hi >= 2050           # legitimate minimal variation survives
    assert hi < 200000          # the gross outlier is still removed
    assert 0.0 <= lo < 2000     # lower bound stays below the bulk, non-negative


def test_mad_bounds_floor_does_not_widen_legitimate_nonzero_mad():
    """The minimum-width floor must only ever rescue a degenerate MAD=0; it must
    not widen a healthy, non-degenerate spread."""
    from aria.scripts.rna_qc import _mad_bounds

    counts = [1000, 1200, 1500, 1800, 2000, 2500, 3000, 4000, 6000, 10000]
    lo_floor, hi_floor = _mad_bounds(counts, log_space=True)
    # Same call with the floor disabled must match (MAD is well above the floor).
    lo_raw, hi_raw = _mad_bounds(counts, log_space=True, min_mad_rel=0.0)
    assert (lo_floor, hi_floor) == (lo_raw, hi_raw)


def _write_processed_qc_h5ad(tmp_path, name, *, mt_values, counts, genes):
    pytest = __import__("pytest")
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    n_cells = len(mt_values)
    obs = pd.DataFrame(
        {
            "nFeature_RNA": genes,
            "nCount_RNA": counts,
            "percent.mt": mt_values,
        },
        index=[f"{name}_cell_{i}" for i in range(n_cells)],
    )
    adata = ad.AnnData(
        X=np.ones((n_cells, 8), dtype=float),
        obs=obs,
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(8)]),
    )
    path = tmp_path / f"{name}.h5ad"
    adata.write_h5ad(path)
    return path


def test_mt_threshold_is_not_driven_by_user_question(tmp_path, monkeypatch):
    pytest = __import__("pytest")
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "numba-cache"))
    monkeypatch.setenv("NUMBA_DISABLE_JIT", "1")
    pytest.importorskip("scanpy")

    from aria.scripts.rna_qc import rna_qc

    neutral_path = _write_processed_qc_h5ad(
        tmp_path,
        "neutral",
        mt_values=[10.0] * 6 + [30.0] * 6,
        counts=[2000] * 12,
        genes=[600] * 12,
    )
    stress_path = _write_processed_qc_h5ad(
        tmp_path,
        "stress",
        mt_values=[10.0] * 6 + [30.0] * 6,
        counts=[2000] * 12,
        genes=[600] * 12,
    )

    neutral = rna_qc(
        {
            "data_path": str(neutral_path),
            "organism": "Homo sapiens",
            "biological_context": {"user_question": "compare baseline groups"},
            "output_dir": str(tmp_path / "neutral_out"),
        }
    )
    stress = rna_qc(
        {
            "data_path": str(stress_path),
            "organism": "Homo sapiens",
            "biological_context": {"user_question": "senescence stress response"},
            "output_dir": str(tmp_path / "stress_out"),
        }
    )

    assert neutral["status"] == "success"
    assert stress["status"] == "success"
    assert neutral["mt_threshold_used"] == stress["mt_threshold_used"]
    assert neutral["n_cells_after"] == stress["n_cells_after"]
    assert neutral["stress_context"] is False
    assert stress["stress_context"] is False


def test_qc_applies_bilateral_mad_for_counts_and_genes(tmp_path, monkeypatch):
    pytest = __import__("pytest")
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "numba-cache"))
    monkeypatch.setenv("NUMBA_DISABLE_JIT", "1")
    pytest.importorskip("scanpy")
    ad = pytest.importorskip("anndata")

    from aria.scripts.rna_qc import rna_qc

    input_path = _write_processed_qc_h5ad(
        tmp_path,
        "outlier",
        mt_values=[5.0] * 11,
        counts=[2000] * 10 + [200000],
        genes=[600] * 10 + [9000],
    )

    result = rna_qc(
        {
            "data_path": str(input_path),
            "organism": "Homo sapiens",
            "output_dir": str(tmp_path / "outlier_out"),
        }
    )

    assert result["status"] == "success"
    assert result["n_cells_before"] == 11
    assert result["n_cells_after"] == 10
    # T6: with a degenerate MAD=0 (10 identical + 1 outlier) the upper bound no
    # longer collapses to the exact median; the minimum-width fallback widens it
    # above the bulk while still removing the gross outlier.
    assert result["max_counts_used"] > 2000.0
    assert result["max_counts_used"] < 200000.0
    assert result["max_genes_used"] > 600.0
    assert result["max_genes_used"] < 9000.0

    filtered = ad.read_h5ad(result["output_path"])
    assert "outlier_cell_10" not in set(filtered.obs_names)
