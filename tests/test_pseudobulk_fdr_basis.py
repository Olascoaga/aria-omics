"""B2 (audit 2026-06-11): the pseudobulk local vs global BH families do NOT share
a base hypothesis set.

DESeq2 independent filtering drops low-count genes from the LOCAL padj (sets
padj_local = NaN, excluded from the local BH denominator), while the GLOBAL pool
keeps every gene with a finite p-value. The robustness_multiverse compared the
two as if they were two corrections of one identical test set. The fix DISCLOSES
and QUANTIFIES the gap instead of hiding it.

Failing-first: the disclosure helper does not exist yet (ImportError).
"""

import numpy as np
import pandas as pd

from aria.scripts.rna_pseudobulk_de import _fdr_filtering_basis


def _block(pvals, padj_local):
    return {"results": pd.DataFrame({"pvalue": pvals, "padj_local": padj_local})}


def test_basis_flags_independent_filtered_genes():
    # 5 genes have a finite p-value; one (gene B in block 1) was independent-
    # filtered locally (padj_local NaN) yet still sits in the global pool.
    blocks = [
        _block([0.001, 0.02, 0.5], [0.01, np.nan, 0.6]),
        _block([0.4, 0.001], [0.5, 0.01]),
    ]
    basis = _fdr_filtering_basis(blocks)
    assert basis["n_global_pool_tests"] == 5
    assert basis["n_independent_filtered_into_global"] == 1
    assert basis["same_base_hypothesis_set"] is False
    # the disclosure must name independent filtering, not claim one shared table
    assert "independent filtering" in basis["note"].lower()


def test_basis_same_when_no_filtering():
    blocks = [_block([0.001, 0.02], [0.01, 0.05])]
    basis = _fdr_filtering_basis(blocks)
    assert basis["n_independent_filtered_into_global"] == 0
    assert basis["same_base_hypothesis_set"] is True


def test_basis_ignores_nonfinite_pvalues():
    # a NaN p-value is not part of the global pool at all
    blocks = [_block([0.001, np.nan], [0.01, np.nan])]
    basis = _fdr_filtering_basis(blocks)
    assert basis["n_global_pool_tests"] == 1
    assert basis["n_independent_filtered_into_global"] == 0
