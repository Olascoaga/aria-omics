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


@dataclass
class NegativeControlResult:
    """W-CALIB negative control: false-positive rate under a permuted null.

    The recovery benchmarks (recall + empirical FDR on a dataset WITH signal)
    answer "does ARIA find true effects?". This answers the complementary
    calibration question — "does ARIA stay quiet when there is NO effect?" — by
    permuting the condition labels (destroying any real association) and checking
    the real DE path does not over-call. Every gene called under the permuted
    null is, by construction, a false positive; a well-calibrated method keeps
    the false-positive rate at or below the nominal alpha.
    """
    status: str                             # "pass" | "fail" | "error"
    false_positive_rate: float              # mean fraction of tested genes called over permutations
    max_false_positive_rate: float          # worst single permutation
    nominal_alpha: float                    # the FDR/alpha the method was run at
    n_permutations: int
    n_tested: int                           # genes tested per permutation (last run)
    calls_per_permutation: list[int]        # significant count in each permutation
    tolerances: dict[str, float]
    messages: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "false_positive_rate": round(self.false_positive_rate, 4),
            "max_false_positive_rate": round(self.max_false_positive_rate, 4),
            "nominal_alpha": self.nominal_alpha,
            "n_permutations": self.n_permutations,
            "n_tested": self.n_tested,
            "calls_per_permutation": self.calls_per_permutation,
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


# ── W-CALIB: label-permutation negative controls (empirical type-I / FDR) ─────
#
# A recovery benchmark proves ARIA finds true effects. A negative control proves
# it does not invent them: under permuted (null) labels the false-positive rate
# must stay near the nominal alpha. Together they make the calibration statement
# the v4.6 gate asks for — empirical FDR ≈ nominal on bulk and pseudobulk, not
# just recall.


def run_bulk_de_negative_control(
    dataset: SyntheticBulkDEDataset | None = None,
    *,
    n_permutations: int = 5,
    nominal_alpha: float = 0.05,
    max_false_positive_rate: float = 0.05,
    lfc_min: float = 0.5,
    lfc_shrink: bool = True,
    seed: int = 11,
    **sim_kwargs,
) -> NegativeControlResult:
    """Permute bulk condition labels and confirm ARIA's real bulk DE stays quiet.

    The simulated matrix carries real signal, but the condition labels are
    shuffled across samples before each DE run, so there is NO true association
    and every significant call is a false positive. A calibrated method keeps the
    mean false-positive rate (significant / tested) at or below
    ``max_false_positive_rate``. Requires pydeseq2.
    """
    import numpy as np
    from aria.scripts.rna_bulk_de import _run_deseq2

    if dataset is None:
        dataset = simulate_bulk_dataset(seed=seed, **sim_kwargs)

    tolerances = {
        "nominal_alpha": nominal_alpha,
        "max_false_positive_rate": max_false_positive_rate,
    }
    rng = np.random.default_rng(seed)
    labels = list(dataset.metadata["condition"].values)
    n_tested = 0
    calls: list[int] = []
    rates: list[float] = []
    errors: list[str] = []

    for _ in range(n_permutations):
        permuted = list(labels)
        rng.shuffle(permuted)
        meta = dataset.metadata.copy()
        meta["condition"] = permuted
        result, _w = _run_deseq2(
            dataset.counts, meta, "condition", "COND_B", "COND_A",
            padj_thr=nominal_alpha, lfc_thr=lfc_min, lfc_shrink=lfc_shrink,
        )
        if result.get("status") != "success":
            errors.append(f"{result.get('error_type')} {result.get('details', '')}")
            continue
        # Genes actually tested = rows of the DESeq2 results table (post-filter);
        # fall back to the full matrix if the script shape is unavailable.
        res_df = result.get("results")
        tested = int(len(res_df)) if res_df is not None else 0
        if not tested:
            tested = int(dataset.counts.shape[0])
        n_called = len(result.get("sig_genes", []) or [])
        n_tested = tested
        calls.append(n_called)
        rates.append(n_called / max(tested, 1))

    if not rates:
        return NegativeControlResult(
            status="error", false_positive_rate=1.0, max_false_positive_rate=1.0,
            nominal_alpha=nominal_alpha, n_permutations=n_permutations,
            n_tested=0, calls_per_permutation=calls, tolerances=tolerances,
            messages=["bulk DE did not succeed under any permutation: "
                      + "; ".join(errors)],
        )

    mean_fpr = float(sum(rates) / len(rates))
    worst = float(max(rates))
    status = "pass" if mean_fpr <= max_false_positive_rate else "fail"
    messages = [
        f"label-permutation null: mean false-positive rate={mean_fpr:.4f} "
        f"(<= {max_false_positive_rate}; nominal alpha={nominal_alpha}); "
        f"worst={worst:.4f}; calls/perm={calls}; n_tested={n_tested}",
    ]
    return NegativeControlResult(
        status=status, false_positive_rate=mean_fpr, max_false_positive_rate=worst,
        nominal_alpha=nominal_alpha, n_permutations=len(rates), n_tested=n_tested,
        calls_per_permutation=calls, tolerances=tolerances, messages=messages,
    )


def run_pseudobulk_de_negative_control(
    dataset: SyntheticDEDataset | None = None,
    *,
    workdir: str | None = None,
    n_permutations: int = 4,
    nominal_alpha: float = 0.05,
    max_false_positive_rate: float = 0.05,
    lfc_min: float = 0.5,
    seed: int = 11,
    **sim_kwargs,
) -> NegativeControlResult:
    """Permute the donor→condition map and confirm pseudobulk DE stays quiet.

    Donors are the replication unit (nested in condition), so the permutation
    reassigns whole donors to conditions at random with balanced group sizes —
    a valid label-permutation null that preserves the pseudobulk structure. Under
    it the real pseudobulk DE path should call essentially nothing; the mean
    false-positive rate must stay at or below ``max_false_positive_rate``.
    Requires pydeseq2.
    """
    import tempfile
    from pathlib import Path
    import numpy as np
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de

    if dataset is None:
        dataset = simulate_pseudobulk_dataset(seed=seed, **sim_kwargs)

    tolerances = {
        "nominal_alpha": nominal_alpha,
        "max_false_positive_rate": max_false_positive_rate,
    }
    rng = np.random.default_rng(seed)
    obs = dataset.adata.obs
    donors = list(dict.fromkeys(obs["donor"].tolist()))   # stable unique order
    n_a = int((obs.drop_duplicates("donor")["condition"] == "COND_A").sum())
    tmp_root = workdir or tempfile.mkdtemp(prefix="aria_wcalib_neg_")

    calls: list[int] = []
    rates: list[float] = []
    errors: list[str] = []
    n_tested = 0

    for i in range(n_permutations):
        shuffled = list(donors)
        rng.shuffle(shuffled)
        # First n_a donors -> COND_A, the rest -> COND_B (balanced like the truth).
        donor_to_cond = {d: ("COND_A" if j < n_a else "COND_B")
                         for j, d in enumerate(shuffled)}
        adata = dataset.adata.copy()
        adata.obs = adata.obs.copy()
        adata.obs["condition"] = [donor_to_cond[d] for d in adata.obs["donor"]]

        run_dir = str(Path(tmp_root) / f"perm{i}")
        Path(run_dir).mkdir(parents=True, exist_ok=True)
        h5ad_path = str(Path(run_dir) / "permuted.h5ad")
        adata.write_h5ad(h5ad_path)

        result = rna_pseudobulk_de({
            "data_path": h5ad_path,
            "groupby": "ctype",
            "condition_col": "condition",
            "replicate_col": "donor",
            "comparisons": [["COND_B", "COND_A"]],
            "use_raw": False,
            "min_replicates_per_condition": 3,
            "padj_max": nominal_alpha,
            "lfc_min": lfc_min,
            "output_dir": run_dir,
        })
        if result.get("status") != "success":
            errors.append(f"{result.get('error_type')} {result.get('details', '')}")
            continue
        block = (
            result.get("per_group", {}).get("ctype0", {})
            .get("per_comparison", {}).get("COND_B_vs_COND_A", {})
        )
        n_called = len(block.get("all_sig", []) or [])
        tested = int(block.get("n_tested") or block.get("n_genes_tested") or 0)
        if not tested:
            tested = int(dataset.adata.shape[1])
        n_tested = tested
        calls.append(n_called)
        rates.append(n_called / max(tested, 1))

    if not rates:
        return NegativeControlResult(
            status="error", false_positive_rate=1.0, max_false_positive_rate=1.0,
            nominal_alpha=nominal_alpha, n_permutations=n_permutations,
            n_tested=0, calls_per_permutation=calls, tolerances=tolerances,
            messages=["pseudobulk DE did not succeed under any permutation: "
                      + "; ".join(errors)],
        )

    mean_fpr = float(sum(rates) / len(rates))
    worst = float(max(rates))
    status = "pass" if mean_fpr <= max_false_positive_rate else "fail"
    messages = [
        f"donor-permutation null: mean false-positive rate={mean_fpr:.4f} "
        f"(<= {max_false_positive_rate}; nominal alpha={nominal_alpha}); "
        f"worst={worst:.4f}; calls/perm={calls}; n_tested={n_tested}",
    ]
    return NegativeControlResult(
        status=status, false_positive_rate=mean_fpr, max_false_positive_rate=worst,
        nominal_alpha=nominal_alpha, n_permutations=len(rates), n_tested=n_tested,
        calls_per_permutation=calls, tolerances=tolerances, messages=messages,
    )


def run_calibration_suite(
    *,
    seed: int = 11,
    quick: bool = False,
) -> dict[str, Any]:
    """Run the full W-CALIB suite and assemble a structured calibration manifest.

    Combines recovery (recall + empirical FDR) with the label-permutation
    negative control (false-positive rate ≈ nominal) for both the bulk and the
    pseudobulk DE paths, and returns a single dict ready to embed in a report's
    provenance / methodology. Requires pydeseq2; the caller must gate on it.

    ``quick`` shrinks the simulations for the doctor command; the pytest gate and
    CI use the full configs. This never fabricates: it runs the REAL DE code and
    records exactly what it measured, including an overall ``status``.
    """
    # The bulk path uses a conservative lfcThreshold-in-Wald test, so a too-small
    # matrix is underpowered for recall; the full config is the genuinely powered
    # calibration gate. ``quick`` is a faster smoke for the doctor command (it
    # uses looser recovery tolerances at the call site, below).
    if quick:
        bulk_kw = dict(n_genes=600, n_de=60, replicates_per_condition=6)
        pb_kw = dict(n_genes=600, n_de=60, donors_per_condition=5, cells_per_donor=40)
        n_perms = 3
        bulk_min_recall = pb_min_recall = 0.4
    else:
        bulk_kw = dict(n_genes=1000, n_de=120, replicates_per_condition=6)
        pb_kw = dict(n_genes=1200, n_de=120, donors_per_condition=6, cells_per_donor=80)
        n_perms = 3
        bulk_min_recall = pb_min_recall = 0.5

    bulk_recovery = run_bulk_de_benchmark(seed=seed, min_recall=bulk_min_recall, **bulk_kw)
    bulk_neg = run_bulk_de_negative_control(seed=seed, n_permutations=n_perms, **bulk_kw)
    pb_recovery = run_pseudobulk_de_benchmark(seed=seed, min_recall=pb_min_recall, **pb_kw)
    pb_neg = run_pseudobulk_de_negative_control(seed=seed, n_permutations=n_perms, **pb_kw)

    paths = {
        "bulk": {
            "recovery": bulk_recovery.as_dict(),
            "negative_control": bulk_neg.as_dict(),
        },
        "pseudobulk": {
            "recovery": pb_recovery.as_dict(),
            "negative_control": pb_neg.as_dict(),
        },
    }
    all_results = [bulk_recovery, bulk_neg, pb_recovery, pb_neg]
    if all(r.status == "pass" for r in all_results):
        status = "pass"
    elif any(r.status == "error" for r in all_results):
        status = "error"
    else:
        status = "fail"

    return {
        "status": status,
        "measured": True,
        "seed": seed,
        "quick": quick,
        "paths": paths,
        "summary": {
            "bulk_recall": bulk_recovery.recall,
            "bulk_empirical_fdr": bulk_recovery.empirical_fdr,
            "bulk_null_fpr": bulk_neg.false_positive_rate,
            "pseudobulk_recall": pb_recovery.recall,
            "pseudobulk_empirical_fdr": pb_recovery.empirical_fdr,
            "pseudobulk_null_fpr": pb_neg.false_positive_rate,
        },
    }
