"""Synthetic differential-expression benchmark with a known ground truth.

X6 (senior audit 2026-05-28): end-to-end tests prove flow, not biological
accuracy. This module simulates a single-cell-like dataset from a negative
binomial model where the set of truly differential genes is KNOWN, runs ARIA's
real pseudobulk DE path against it, and reports recovery metrics (recall on
true-DE genes and empirical FDR on null genes). A dependency change that breaks
DE precision then fails a test instead of shipping silently.

Per ADR-011 the dataset uses only neutral synthetic labels (GENE_####,
COND_A/COND_B, donor d##, cell type ctype0); it is a generated ground-truth
benchmark, not a named golden dataset.

The simulator is dependency-light (numpy/pandas/anndata). pydeseq2 is only
needed to *run* the benchmark, imported lazily by the pseudobulk script.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class SyntheticDEDataset:
    """A simulated dataset plus its ground truth."""
    adata: Any                              # AnnData (cells x genes), integer counts in X
    de_genes: dict[str, str]                # gene -> "up" | "down" (truth in COND_B vs COND_A)
    null_genes: list[str]                   # genes with no true effect
    params: dict[str, Any] = field(default_factory=dict)

    @property
    def n_de(self) -> int:
        return len(self.de_genes)


@dataclass
class DEBenchmarkResult:
    status: str                             # "pass" | "fail" | "error"
    recall: float                           # fraction of true-DE genes recovered
    empirical_fdr: float                    # fraction of calls that are truly null
    n_true_de: int
    n_called: int
    n_true_positive: int
    n_false_positive: int
    tolerances: dict[str, float]
    messages: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "recall": round(self.recall, 4),
            "empirical_fdr": round(self.empirical_fdr, 4),
            "n_true_de": self.n_true_de,
            "n_called": self.n_called,
            "n_true_positive": self.n_true_positive,
            "n_false_positive": self.n_false_positive,
            "tolerances": self.tolerances,
            "messages": self.messages,
        }


def simulate_pseudobulk_dataset(
    *,
    n_genes: int = 1500,
    n_de: int = 150,
    donors_per_condition: int = 6,
    cells_per_donor: int = 80,
    dispersion: float = 0.2,
    min_abs_log2fc: float = 1.0,
    max_abs_log2fc: float = 2.0,
    seed: int = 7,
) -> SyntheticDEDataset:
    """Simulate an unpaired two-condition single-cell dataset.

    Negative binomial counts (Gamma-Poisson). Donors are nested in condition
    (distinct donors per group), so pseudobulk aggregation by donor is the
    biological replication unit. ``n_de`` genes get a fold change in COND_B; the
    rest are null. Deterministic for a given seed.
    """
    import numpy as np
    import pandas as pd
    import anndata as ad

    rng = np.random.default_rng(seed)

    gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    # Base mean expression per gene (lognormal spread, floored so DESeq2 has signal).
    base_mean = np.clip(np.exp(rng.normal(1.4, 1.0, size=n_genes)), 1.0, None)

    # Pick true-DE genes and assign a per-gene fold change / direction.
    de_idx = rng.choice(n_genes, size=n_de, replace=False)
    de_mask = np.zeros(n_genes, dtype=bool)
    de_mask[de_idx] = True
    log2fc = np.zeros(n_genes, dtype=float)
    signs = rng.choice([-1.0, 1.0], size=n_de)
    mags = rng.uniform(min_abs_log2fc, max_abs_log2fc, size=n_de)
    log2fc[de_idx] = signs * mags
    cond_b_fc = np.power(2.0, log2fc)   # multiplicative factor applied in COND_B

    conditions = ["COND_A", "COND_B"]
    rows_counts = []
    obs_records = []
    cell_counter = 0
    for ci, cond in enumerate(conditions):
        fc = cond_b_fc if cond == "COND_B" else np.ones(n_genes)
        for d in range(donors_per_condition):
            donor = f"d{ci * donors_per_condition + d:02d}"
            # Per-donor size factor adds realistic between-replicate variation.
            donor_sf = float(np.exp(rng.normal(0.0, 0.2)))
            gene_mean = base_mean * fc * donor_sf            # (n_genes,)
            for _ in range(cells_per_donor):
                # Gamma-Poisson NB sampling, vectorized over genes.
                shape = 1.0 / dispersion
                gamma = rng.gamma(shape=shape, scale=gene_mean * dispersion)
                counts = rng.poisson(gamma)
                rows_counts.append(counts.astype(np.int64))
                obs_records.append({
                    "condition": cond,
                    "donor": donor,
                    "ctype": "ctype0",
                })
                cell_counter += 1

    X = np.vstack(rows_counts)
    obs = pd.DataFrame(obs_records, index=[f"cell_{i}" for i in range(cell_counter)])
    var = pd.DataFrame(index=gene_names)
    adata = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var)

    de_genes = {
        gene_names[i]: ("up" if log2fc[i] > 0 else "down")
        for i in de_idx
    }
    null_genes = [gene_names[i] for i in range(n_genes) if not de_mask[i]]

    return SyntheticDEDataset(
        adata=adata,
        de_genes=de_genes,
        null_genes=null_genes,
        params={
            "n_genes": n_genes, "n_de": n_de,
            "donors_per_condition": donors_per_condition,
            "cells_per_donor": cells_per_donor,
            "dispersion": dispersion,
            "min_abs_log2fc": min_abs_log2fc,
            "max_abs_log2fc": max_abs_log2fc,
            "seed": seed,
        },
    )


@dataclass
class SyntheticBulkDEDataset:
    """A simulated bulk RNA-seq count matrix plus its ground truth."""
    counts: Any                             # DataFrame (genes x samples), integer counts
    metadata: Any                           # DataFrame (samples x design), with "condition"
    de_genes: dict[str, str]                # gene -> "up" | "down" (truth in COND_B vs COND_A)
    null_genes: list[str]                   # genes with no true effect
    params: dict[str, Any] = field(default_factory=dict)

    @property
    def n_de(self) -> int:
        return len(self.de_genes)


def simulate_bulk_dataset(
    *,
    n_genes: int = 2000,
    n_de: int = 200,
    replicates_per_condition: int = 5,
    dispersion: float = 0.2,
    min_abs_log2fc: float = 1.0,
    max_abs_log2fc: float = 2.5,
    seed: int = 11,
) -> SyntheticBulkDEDataset:
    """Simulate a two-condition bulk RNA-seq count matrix (genes x samples).

    Negative binomial counts (Gamma-Poisson) at bulk sequencing depth, with one
    sample per biological replicate (the replication unit for bulk DESeq2).
    ``n_de`` genes get a fold change in COND_B; the rest are null. Deterministic
    for a given seed. Neutral labels only (ADR-011).
    """
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)

    gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    # Bulk depth: a higher base mean than the single-cell simulator.
    base_mean = np.clip(np.exp(rng.normal(3.0, 1.0, size=n_genes)), 5.0, None)

    de_idx = rng.choice(n_genes, size=n_de, replace=False)
    de_mask = np.zeros(n_genes, dtype=bool)
    de_mask[de_idx] = True
    log2fc = np.zeros(n_genes, dtype=float)
    signs = rng.choice([-1.0, 1.0], size=n_de)
    mags = rng.uniform(min_abs_log2fc, max_abs_log2fc, size=n_de)
    log2fc[de_idx] = signs * mags
    cond_b_fc = np.power(2.0, log2fc)

    sample_cols, sample_names, conds = [], [], []
    for cond in ("COND_A", "COND_B"):
        fc = cond_b_fc if cond == "COND_B" else np.ones(n_genes)
        for r in range(replicates_per_condition):
            sample = f"{cond}_r{r}"
            # Per-sample size factor for realistic between-replicate variation.
            sf = float(np.exp(rng.normal(0.0, 0.15)))
            gene_mean = base_mean * fc * sf
            shape = 1.0 / dispersion
            gamma = rng.gamma(shape=shape, scale=gene_mean * dispersion)
            sample_cols.append(rng.poisson(gamma).astype(np.int64))
            sample_names.append(sample)
            conds.append(cond)

    counts = pd.DataFrame(
        np.array(sample_cols).T, index=gene_names, columns=sample_names)
    metadata = pd.DataFrame({"condition": conds}, index=sample_names)

    de_genes = {gene_names[i]: ("up" if log2fc[i] > 0 else "down") for i in de_idx}
    null_genes = [gene_names[i] for i in range(n_genes) if not de_mask[i]]

    return SyntheticBulkDEDataset(
        counts=counts,
        metadata=metadata,
        de_genes=de_genes,
        null_genes=null_genes,
        params={
            "n_genes": n_genes, "n_de": n_de,
            "replicates_per_condition": replicates_per_condition,
            "dispersion": dispersion,
            "min_abs_log2fc": min_abs_log2fc,
            "max_abs_log2fc": max_abs_log2fc,
            "seed": seed,
        },
    )


def run_bulk_de_benchmark(
    dataset: SyntheticBulkDEDataset | None = None,
    *,
    min_recall: float = 0.5,
    max_empirical_fdr: float = 0.2,
    padj_max: float = 0.05,
    lfc_min: float = 0.5,
    lfc_shrink: bool = True,
    **sim_kwargs,
) -> DEBenchmarkResult:
    """Run ARIA's real bulk DE (`_run_deseq2`, incl. apeGLM) on the synthetic
    matrix and score recovery — the bulk analog of the pseudobulk benchmark.

    Requires pydeseq2. Returns a :class:`DEBenchmarkResult`; ``status`` is "pass"
    only when recall and empirical FDR are within tolerance.
    """
    from aria.scripts.rna_bulk_de import _run_deseq2

    if dataset is None:
        dataset = simulate_bulk_dataset(**sim_kwargs)

    tolerances = {"min_recall": min_recall, "max_empirical_fdr": max_empirical_fdr}
    result, _warnings = _run_deseq2(
        dataset.counts, dataset.metadata, "condition", "COND_B", "COND_A",
        padj_thr=padj_max, lfc_thr=lfc_min, lfc_shrink=lfc_shrink,
    )

    if result.get("status") != "success":
        return DEBenchmarkResult(
            status="error", recall=0.0, empirical_fdr=1.0,
            n_true_de=dataset.n_de, n_called=0,
            n_true_positive=0, n_false_positive=0,
            tolerances=tolerances,
            messages=[f"bulk DE did not succeed: {result.get('error_type')} "
                      f"{result.get('details', '')}"],
        )

    called = set(result.get("sig_genes", []) or [])
    true_de = set(dataset.de_genes)
    null = set(dataset.null_genes)
    tp = called & true_de
    fp = called & null
    recall = len(tp) / max(len(true_de), 1)
    empirical_fdr = len(fp) / max(len(called), 1)

    messages = [
        f"recall={recall:.3f} (>= {min_recall}); "
        f"empirical_fdr={empirical_fdr:.3f} (<= {max_empirical_fdr}); "
        f"called={len(called)}, true_de={len(true_de)}",
    ]
    status = "pass" if (recall >= min_recall and empirical_fdr <= max_empirical_fdr) else "fail"

    return DEBenchmarkResult(
        status=status,
        recall=recall,
        empirical_fdr=empirical_fdr,
        n_true_de=len(true_de),
        n_called=len(called),
        n_true_positive=len(tp),
        n_false_positive=len(fp),
        tolerances=tolerances,
        messages=messages,
    )


def run_pseudobulk_de_benchmark(
    dataset: SyntheticDEDataset | None = None,
    *,
    workdir: str | None = None,
    min_recall: float = 0.5,
    max_empirical_fdr: float = 0.2,
    padj_max: float = 0.05,
    lfc_min: float = 0.5,
    **sim_kwargs,
) -> DEBenchmarkResult:
    """Run ARIA's real pseudobulk DE on the synthetic data and score recovery.

    Requires pydeseq2 (imported by the pseudobulk script). Returns a
    :class:`DEBenchmarkResult`; ``status`` is "pass" only when recall and
    empirical FDR are within tolerance.
    """
    import tempfile
    from pathlib import Path
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de

    if dataset is None:
        dataset = simulate_pseudobulk_dataset(**sim_kwargs)

    tmp = workdir or tempfile.mkdtemp(prefix="aria_x6_")
    h5ad_path = str(Path(tmp) / "synthetic.h5ad")
    dataset.adata.write_h5ad(h5ad_path)

    result = rna_pseudobulk_de({
        "data_path": h5ad_path,
        "groupby": "ctype",
        "condition_col": "condition",
        "replicate_col": "donor",
        "comparisons": [["COND_B", "COND_A"]],
        "use_raw": False,
        "min_replicates_per_condition": 3,
        "padj_max": padj_max,
        "lfc_min": lfc_min,
        "output_dir": tmp,
    })

    if result.get("status") != "success":
        return DEBenchmarkResult(
            status="error", recall=0.0, empirical_fdr=1.0,
            n_true_de=dataset.n_de, n_called=0,
            n_true_positive=0, n_false_positive=0,
            tolerances={"min_recall": min_recall, "max_empirical_fdr": max_empirical_fdr},
            messages=[f"pseudobulk DE did not succeed: {result.get('error_type')} "
                      f"{result.get('details', '')}"],
        )

    block = (
        result.get("per_group", {})
        .get("ctype0", {})
        .get("per_comparison", {})
        .get("COND_B_vs_COND_A", {})
    )
    called = {rec["gene"] for rec in block.get("all_sig", []) or []}

    true_de = set(dataset.de_genes)
    null = set(dataset.null_genes)
    tp = called & true_de
    fp = called & null
    recall = len(tp) / max(len(true_de), 1)
    empirical_fdr = len(fp) / max(len(called), 1)

    messages = [
        f"recall={recall:.3f} (>= {min_recall}); "
        f"empirical_fdr={empirical_fdr:.3f} (<= {max_empirical_fdr}); "
        f"called={len(called)}, true_de={len(true_de)}",
    ]
    status = "pass" if (recall >= min_recall and empirical_fdr <= max_empirical_fdr) else "fail"

    return DEBenchmarkResult(
        status=status,
        recall=recall,
        empirical_fdr=empirical_fdr,
        n_true_de=len(true_de),
        n_called=len(called),
        n_true_positive=len(tp),
        n_false_positive=len(fp),
        tolerances={"min_recall": min_recall, "max_empirical_fdr": max_empirical_fdr},
        messages=messages,
    )
