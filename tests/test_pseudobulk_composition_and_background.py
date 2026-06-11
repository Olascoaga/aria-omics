"""A3 (audit 2026-06-11): two silent-degradation fixes in rna_pseudobulk_de.

1. Composition covariate: when collinearity with the contrast cannot be MEASURED
   (no variance in the composition or the condition vector within a block),
   `_abs_corr` returns None and the OLD code fell through to ADD the covariate.
   Conservative fix: omit it and record `collinearity_unmeasurable`.
2. Expressed-gene background: a failure in the mask computation fell back to all
   genes SILENTLY. Fix: return an explicit degradation flag.

Failing-first: the helpers do not exist yet (ImportError).
"""

import numpy as np

from aria.scripts.rna_pseudobulk_de import (
    _composition_covariate_decision,
    _expressed_background,
)

_MAX = 0.8


def test_composition_omitted_when_collinearity_unmeasurable():
    # zero variance in the composition vector -> correlation is undefined
    comp = np.array([1.0, 1.0, 1.0, 1.0])
    cond = np.array([0.0, 0.0, 1.0, 1.0])
    use, reason = _composition_covariate_decision(comp, cond, _MAX)
    assert use is False
    assert reason == "collinearity_unmeasurable"


def test_composition_omitted_when_collinear():
    comp = np.array([0.1, 0.15, 0.95, 1.0])
    cond = np.array([0.0, 0.0, 1.0, 1.0])  # near-perfect correlation
    use, reason = _composition_covariate_decision(comp, cond, _MAX)
    assert use is False
    assert "collinear_with_condition" in reason


def test_composition_used_when_independent():
    comp = np.array([0.5, 0.1, 0.6, 0.2])
    cond = np.array([0.0, 0.0, 1.0, 1.0])  # |corr| ~ 0.24 < 0.8
    use, reason = _composition_covariate_decision(comp, cond, _MAX)
    assert use is True
    assert reason is None


def test_expressed_background_no_degradation_on_valid_counts():
    counts = np.array([[0, 1], [2, 0], [0, 0]])  # gene0 and gene1 both expressed
    genes, degraded = _expressed_background(counts, ["g0", "g1"])
    assert degraded is False
    assert set(genes) == {"g0", "g1"}


def test_expressed_background_degraded_falls_back_to_all_genes():
    class _Boom:
        def __gt__(self, other):
            raise RuntimeError("cannot compare")

    genes, degraded = _expressed_background(_Boom(), ["g0", "g1"])
    assert degraded is True
    assert genes == ["g0", "g1"]
