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
    true_log2fc: dict[str, float] = field(default_factory=dict)
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
        true_log2fc={gene_names[i]: float(log2fc[i]) for i in range(n_genes)},
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


def simulate_pseudoreplication_null_dataset(
    *,
    n_genes: int = 600,
    donors_per_condition: int = 4,
    cells_per_donor: int = 80,
    dispersion: float = 0.2,
    donor_gene_sd: float = 0.9,
    seed: int = 31,
) -> SyntheticDEDataset:
    """Simulate a no-condition-effect scRNA dataset with donor heterogeneity.

    There is no planted COND_B-vs-COND_A effect. Donors carry stable per-gene
    expression offsets, so a cell-level test that treats cells as independent
    replicates is anti-conservative, while donor-aware pseudobulk keeps the
    inferential unit at the biological replicate.
    """
    import numpy as np
    import pandas as pd
    import anndata as ad

    rng = np.random.default_rng(seed)
    gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    base_mean = np.clip(np.exp(rng.normal(1.5, 1.0, size=n_genes)), 1.0, None)

    rows_counts = []
    obs_records = []
    cell_counter = 0
    for ci, cond in enumerate(("COND_A", "COND_B")):
        for d in range(donors_per_condition):
            donor = f"d{ci * donors_per_condition + d:02d}"
            donor_gene_log2 = rng.normal(0.0, donor_gene_sd, size=n_genes)
            donor_gene_fc = np.power(2.0, donor_gene_log2)
            donor_sf = float(np.exp(rng.normal(0.0, 0.12)))
            gene_mean = base_mean * donor_gene_fc * donor_sf
            for _ in range(cells_per_donor):
                cell_sf = float(np.exp(rng.normal(0.0, 0.08)))
                shape = 1.0 / dispersion
                gamma = rng.gamma(shape=shape, scale=gene_mean * cell_sf * dispersion)
                rows_counts.append(rng.poisson(gamma).astype(np.int64))
                obs_records.append({
                    "condition": cond,
                    "donor": donor,
                    "ctype": "ctype0",
                })
                cell_counter += 1

    X = np.vstack(rows_counts)
    obs = pd.DataFrame(obs_records, index=[f"cell_{i}" for i in range(cell_counter)])
    var = pd.DataFrame(index=gene_names)
    return SyntheticDEDataset(
        adata=ad.AnnData(X=X.astype(np.float32), obs=obs, var=var),
        de_genes={},
        null_genes=gene_names,
        true_log2fc={gene: 0.0 for gene in gene_names},
        params={
            "n_genes": n_genes,
            "n_de": 0,
            "donors_per_condition": donors_per_condition,
            "cells_per_donor": cells_per_donor,
            "dispersion": dispersion,
            "donor_gene_sd": donor_gene_sd,
            "seed": seed,
            "null_model": "donor_gene_heterogeneity_no_condition_effect",
        },
    )


def _bh_adjust(pvalues: Any) -> Any:
    try:
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(pvalues, method="fdr_bh")
        return padj
    except Exception:
        from aria.utils.stats import bh_correct
        return bh_correct(pvalues)


def _naive_cell_level_null_calls(
    dataset: SyntheticDEDataset,
    *,
    nominal_alpha: float = 0.05,
    min_abs_log2fc: float = 0.0,
) -> dict[str, Any]:
    """Cell-level Welch tests that intentionally ignore donor structure."""
    import numpy as np

    X = np.asarray(dataset.adata.X, dtype=float)
    cond = np.asarray(dataset.adata.obs["condition"].values)
    a = np.log1p(X[cond == "COND_A"])
    b = np.log1p(X[cond == "COND_B"])
    try:
        from scipy import stats
        _stat, pvals = stats.ttest_ind(b, a, axis=0, equal_var=False, nan_policy="omit")
    except Exception:
        # Normal approximation fallback; the benchmark lane normally has SciPy.
        import math
        mean_diff = b.mean(axis=0) - a.mean(axis=0)
        se = np.sqrt(b.var(axis=0, ddof=1) / max(b.shape[0], 1)
                     + a.var(axis=0, ddof=1) / max(a.shape[0], 1))
        z = np.divide(mean_diff, se, out=np.zeros_like(mean_diff), where=se > 0)
        pvals = np.array([math.erfc(abs(float(v)) / math.sqrt(2.0)) for v in z])
    pvals = np.asarray(pvals, dtype=float)
    pvals[~np.isfinite(pvals)] = 1.0
    padj = np.asarray(_bh_adjust(pvals), dtype=float)

    mean_a = np.asarray(X[cond == "COND_A"].mean(axis=0), dtype=float)
    mean_b = np.asarray(X[cond == "COND_B"].mean(axis=0), dtype=float)
    log2fc = np.log2((mean_b + 1.0) / (mean_a + 1.0))
    called_mask = (padj < nominal_alpha) & (np.abs(log2fc) >= min_abs_log2fc)
    genes = list(dataset.adata.var_names)
    called = [genes[i] for i, flag in enumerate(called_mask) if bool(flag)]
    return {
        "method": "naive_cell_level_welch_treats_cells_as_replicates",
        "status": "success",
        "nominal_alpha": nominal_alpha,
        "min_abs_log2fc": min_abs_log2fc,
        "n_cells_condition_a": int(a.shape[0]),
        "n_cells_condition_b": int(b.shape[0]),
        "n_tested": int(len(genes)),
        "n_called": int(len(called)),
        "false_positive_rate": round(len(called) / max(len(genes), 1), 4),
        "called_genes_top": called[:25],
    }


def _run_pseudobulk_null_calls(
    dataset: SyntheticDEDataset,
    *,
    workdir: str | None = None,
    nominal_alpha: float = 0.05,
    lfc_min: float = 0.0,
) -> dict[str, Any]:
    """Run ARIA pseudobulk DE on a known-null dataset and count false calls."""
    import tempfile
    from pathlib import Path
    from aria.scripts.rna_pseudobulk_de import rna_pseudobulk_de

    tmp = workdir or tempfile.mkdtemp(prefix="aria_a2_null_")
    h5ad_path = str(Path(tmp) / "pseudoreplication_null.h5ad")
    dataset.adata.write_h5ad(h5ad_path)
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
        "fdr_strategy": "per_cluster",
        "output_dir": tmp,
    })
    if result.get("status") != "success":
        return {
            "method": "aria_donor_aware_pseudobulk",
            "status": "error",
            "error_type": result.get("error_type"),
            "details": result.get("details", ""),
            "nominal_alpha": nominal_alpha,
            "n_tested": int(dataset.adata.shape[1]),
            "n_called": 0,
            "false_positive_rate": 1.0,
        }
    block = (
        result.get("per_group", {})
        .get("ctype0", {})
        .get("per_comparison", {})
        .get("COND_B_vs_COND_A", {})
    )
    called = [rec.get("gene") for rec in block.get("all_sig", []) or []]
    n_tested = int(block.get("n_tested") or block.get("n_genes_tested") or dataset.adata.shape[1])
    return {
        "method": "aria_donor_aware_pseudobulk",
        "status": "success",
        "nominal_alpha": nominal_alpha,
        "lfc_min": lfc_min,
        "n_biological_replicates_condition_a": int(
            dataset.adata.obs.drop_duplicates("donor")
            .query("condition == 'COND_A'")
            .shape[0]
        ),
        "n_biological_replicates_condition_b": int(
            dataset.adata.obs.drop_duplicates("donor")
            .query("condition == 'COND_B'")
            .shape[0]
        ),
        "n_tested": n_tested,
        "n_called": int(len(called)),
        "false_positive_rate": round(len(called) / max(n_tested, 1), 4),
        "called_genes_top": [g for g in called if g][:25],
    }


def write_pseudobulk_a2_figure(manifest: dict[str, Any], path: str) -> str:
    """Write a dependency-free SVG summary for preliminary Fig. 2."""
    from pathlib import Path

    axes = manifest.get("axes", {})
    recovery = axes.get("donor_aware_recovery", {})
    guard = axes.get("pseudoreplication_guard", {})
    vals = [
        ("Pseudobulk recall", float(recovery.get("recall", 0.0)), "#2F6F73"),
        ("Pseudobulk empirical FDR", float(recovery.get("empirical_fdr", 1.0)), "#7C3AED"),
        ("Pseudobulk null FPR", float(guard.get("pseudobulk_false_positive_rate", 1.0)), "#2F6F73"),
        ("Naive cell-level null FPR", float(guard.get("naive_cell_level_false_positive_rate", 0.0)), "#A23B3B"),
    ]
    width, height = 790, 410
    left, top = 230, 74
    bar_w, bar_h, gap = 410, 40, 27
    rows = []
    for i, (label, value, fill) in enumerate(vals):
        y = top + i * (bar_h + gap)
        v = max(0.0, min(1.0, value))
        rows.append(
            f'<text x="28" y="{y + 26}" font-size="16" fill="#1f2933">{label}</text>'
            f'<rect x="{left}" y="{y}" width="{bar_w}" height="{bar_h}" fill="#e8edf0"/>'
            f'<rect x="{left}" y="{y}" width="{bar_w * v:.1f}" height="{bar_h}" fill="{fill}"/>'
            f'<text x="{left + bar_w + 18}" y="{y + 26}" font-size="16" fill="#1f2933">{value:.3f}</text>'
        )
    status = str(manifest.get("status", "unknown")).upper()
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="#ffffff"/>'
        '<text x="28" y="36" font-size="22" font-weight="700" fill="#111827">'
        'Fig. 2 preliminary: A2 donor-aware pseudobulk</text>'
        f'<text x="660" y="36" font-size="18" font-weight="700" fill="#2F6F73">{status}</text>'
        '<text x="28" y="382" font-size="13" fill="#4b5563">'
        'Synthetic donor-null anti-pattern; Kang + muscat external lane remains pending.</text>'
        + "".join(rows)
        + "</svg>"
    )
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg, encoding="utf-8")
    return str(out)


def run_pseudobulk_a2_benchmark(
    *,
    seed: int = 23,
    quick: bool = False,
    output_dir: str | None = None,
    manifest_name: str = "a2_pseudobulk_v4.5.5.json",
    figure_name: str = "fig2_a2_pseudobulk_v4.5.5.svg",
    artifact_version: str = "v4.5.5",
) -> dict[str, Any]:
    """Execute Benchmark A2's preliminary donor-aware pseudobulk lane."""
    import json
    import subprocess
    from pathlib import Path

    recovery_kw = (
        dict(n_genes=700, n_de=70, donors_per_condition=5, cells_per_donor=45)
        if quick else
        dict(n_genes=1200, n_de=120, donors_per_condition=6, cells_per_donor=80)
    )
    null_kw = (
        dict(n_genes=400, donors_per_condition=4, cells_per_donor=45, donor_gene_sd=0.9)
        if quick else
        dict(n_genes=600, donors_per_condition=4, cells_per_donor=80, donor_gene_sd=0.9)
    )
    recovery = run_pseudobulk_de_benchmark(seed=seed, min_recall=0.5, **recovery_kw)
    null_ds = simulate_pseudoreplication_null_dataset(seed=seed + 8, **null_kw)
    pb_null = _run_pseudobulk_null_calls(null_ds, nominal_alpha=0.05, lfc_min=0.0)
    naive_null = _naive_cell_level_null_calls(null_ds, nominal_alpha=0.05)

    pb_fpr = float(pb_null.get("false_positive_rate", 1.0))
    naive_fpr = float(naive_null.get("false_positive_rate", 0.0))
    inflation_ratio = naive_fpr / max(pb_fpr, 1.0 / max(int(pb_null.get("n_tested", 1)), 1))
    axes = {
        "donor_aware_recovery": recovery.as_dict(),
        "pseudoreplication_guard": {
            "pseudobulk_false_positive_rate": round(pb_fpr, 4),
            "naive_cell_level_false_positive_rate": round(naive_fpr, 4),
            "inflation_ratio_vs_pseudobulk_floor": round(float(inflation_ratio), 4),
            "pseudobulk": pb_null,
            "naive_cell_level": naive_null,
        },
    }
    axis_pass = {
        "donor_aware_recovery": recovery.status == "pass",
        "pseudobulk_null_control": (
            pb_null.get("status") == "success" and pb_fpr <= 0.10
        ),
        "cell_level_antipattern_detected": naive_fpr >= max(0.20, pb_fpr * 3.0),
    }
    try:
        git_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
        git_status = subprocess.check_output(
            ["git", "status", "--porcelain"], text=True
        ).strip()
    except Exception:
        git_commit = "unknown"
        git_status = "unknown"
    manifest = {
        "status": "pass" if all(axis_pass.values()) else "fail",
        "benchmark": "A2_pseudobulk_scrna",
        "benchmark_version": "v1",
        "artifact_version": artifact_version,
        "scope": "preliminary_synthetic_donor_aware",
        "method_under_test": "ARIA donor-aware scRNA pseudobulk DE",
        "external_reference_lane": {
            "status": "pending",
            "reason": "Kang + muscat requires local benchmark data and aria-bench-env",
        },
        "datasets": {
            "recovery": recovery_kw | {"seed": seed},
            "pseudoreplication_null": null_ds.params,
        },
        "tolerances": {
            "min_recall": 0.5,
            "max_recovery_empirical_fdr": 0.2,
            "max_pseudobulk_null_fpr": 0.10,
            "min_naive_cell_level_null_fpr": 0.20,
            "min_naive_vs_pseudobulk_inflation": 3.0,
        },
        "axis_pass": axis_pass,
        "axes": axes,
        "provenance": {
            "git_commit": git_commit,
            "git_dirty": bool(git_status),
            "runtime_note": (
                "artifact_version is fixed for the v4.5 benchmark lane; local "
                "working tree may contain unrelated v4.6 changes."
            ),
        },
        "messages": [
            "A2 preliminary validates that ARIA keeps donor/sample as the "
            "inferential unit and exposes the cell-level pseudoreplication "
            "anti-pattern on a controlled donor-null simulation.",
            "Kang et al. + muscat remains the external A2 lane once local data "
            "and aria-bench-env are available.",
        ],
    }
    if output_dir:
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        manifest_path = outdir / manifest_name
        figure_path = outdir / figure_name
        write_pseudobulk_a2_figure(manifest, str(figure_path))
        manifest["artifacts"] = {
            "manifest_json": str(manifest_path),
            "figure_svg": str(figure_path),
        }
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


@dataclass
class SyntheticBulkDEDataset:
    """A simulated bulk RNA-seq count matrix plus its ground truth."""
    counts: Any                             # DataFrame (genes x samples), integer counts
    metadata: Any                           # DataFrame (samples x design), with "condition"
    de_genes: dict[str, str]                # gene -> "up" | "down" (truth in COND_B vs COND_A)
    null_genes: list[str]                   # genes with no true effect
    true_log2fc: dict[str, float] = field(default_factory=dict)
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
        true_log2fc={gene_names[i]: float(log2fc[i]) for i in range(n_genes)},
        params={
            "n_genes": n_genes, "n_de": n_de,
            "replicates_per_condition": replicates_per_condition,
            "dispersion": dispersion,
            "min_abs_log2fc": min_abs_log2fc,
            "max_abs_log2fc": max_abs_log2fc,
            "seed": seed,
        },
    )


def _rank_biased_overlap(left: list[str], right: list[str], p: float = 0.9) -> float:
    """Finite-list rank-biased overlap for top-list concordance."""
    if not left or not right:
        return 0.0
    depth = min(len(left), len(right))
    seen_left: set[str] = set()
    seen_right: set[str] = set()
    score = 0.0
    for d in range(1, depth + 1):
        seen_left.add(left[d - 1])
        seen_right.add(right[d - 1])
        agreement = len(seen_left & seen_right) / d
        score += (1.0 - p) * (p ** (d - 1)) * agreement
    return float(score)


def _finite_float(value: Any, default: float = 0.0) -> float:
    try:
        import math
        val = float(value)
        return val if math.isfinite(val) else default
    except Exception:
        return default


def score_bulk_de_a1(
    dataset: SyntheticBulkDEDataset,
    de_result: dict[str, Any],
    negative_control: NegativeControlResult,
    *,
    padj_max: float = 0.05,
    lfc_min: float = 0.5,
    top_k: int | None = None,
    min_recall: float = 0.5,
    max_empirical_fdr: float = 0.2,
    max_null_fpr: float = 0.05,
    min_lfc_spearman: float = 0.7,
    min_top_k_jaccard: float = 0.35,
) -> dict[str, Any]:
    """Score Benchmark A1's four bulk-DE axes from a real ARIA DESeq2 result."""
    import pandas as pd

    tolerances = {
        "padj_max": padj_max,
        "lfc_min": lfc_min,
        "min_recall": min_recall,
        "max_empirical_fdr": max_empirical_fdr,
        "max_null_false_positive_rate": max_null_fpr,
        "min_lfc_spearman": min_lfc_spearman,
        "min_top_k_jaccard": min_top_k_jaccard,
    }
    if de_result.get("status") != "success":
        return {
            "status": "error",
            "tolerances": tolerances,
            "messages": [
                f"bulk DE did not succeed: {de_result.get('error_type')} "
                f"{de_result.get('details', '')}".strip()
            ],
        }

    results_df = de_result.get("results")
    if results_df is None or not hasattr(results_df, "copy"):
        return {
            "status": "error",
            "tolerances": tolerances,
            "messages": ["bulk DE result did not include a results DataFrame"],
        }

    df = results_df.copy()
    df.index = [str(idx) for idx in df.index]
    true_de = set(dataset.de_genes)
    null = set(dataset.null_genes)
    called = set(de_result.get("sig_genes", []) or [])
    tp = called & true_de
    fp = called & null

    recall = len(tp) / max(len(true_de), 1)
    empirical_fdr = len(fp) / max(len(called), 1)
    precision = len(tp) / max(len(called), 1)

    truth = pd.Series(dataset.true_log2fc, name="true_log2fc", dtype=float)
    est = df.get("log2FoldChange")
    if est is None:
        est = pd.Series(dtype=float)
    lfc_frame = pd.DataFrame({
        "estimated": pd.to_numeric(est, errors="coerce"),
        "true": truth,
    }).dropna()
    lfc_de_frame = lfc_frame[lfc_frame["true"].abs() > 0]
    if lfc_frame.empty:
        pearson = spearman = pearson_de = spearman_de = 0.0
    else:
        pearson = _finite_float(lfc_frame["estimated"].corr(lfc_frame["true"], method="pearson"))
        spearman = _finite_float(lfc_frame["estimated"].corr(lfc_frame["true"], method="spearman"))
        if lfc_de_frame.empty:
            pearson_de = spearman_de = 0.0
        else:
            pearson_de = _finite_float(
                lfc_de_frame["estimated"].corr(lfc_de_frame["true"], method="pearson")
            )
            spearman_de = _finite_float(
                lfc_de_frame["estimated"].corr(lfc_de_frame["true"], method="spearman")
            )

    if top_k is None:
        top_k = min(100, max(1, dataset.n_de))
    ranked_est = (
        df.assign(
            _padj=pd.to_numeric(df.get("padj"), errors="coerce").fillna(1.0),
            _abs_lfc=pd.to_numeric(df.get("log2FoldChange"), errors="coerce").abs().fillna(0.0),
        )
        .sort_values(["_padj", "_abs_lfc"], ascending=[True, False])
        .index.astype(str)
        .tolist()
    )
    ranked_truth = (
        pd.Series(dataset.true_log2fc, dtype=float)
        .abs()
        .sort_values(ascending=False)
        .index.astype(str)
        .tolist()
    )
    est_top = set(ranked_est[:top_k])
    truth_top = set(ranked_truth[:top_k])
    top_intersection = est_top & truth_top
    top_k_jaccard = len(top_intersection) / max(len(est_top | truth_top), 1)
    top_k_truth_recall = len(top_intersection) / max(len(truth_top), 1)
    rbo = _rank_biased_overlap(ranked_est[:top_k], ranked_truth[:top_k], p=0.9)

    axes = {
        "fdr_calibration": {
            "status": negative_control.status,
            "nominal_alpha": negative_control.nominal_alpha,
            "false_positive_rate": round(negative_control.false_positive_rate, 4),
            "max_false_positive_rate": round(negative_control.max_false_positive_rate, 4),
            "calls_per_permutation": negative_control.calls_per_permutation,
            "n_permutations": negative_control.n_permutations,
            "n_tested": negative_control.n_tested,
        },
        "lfc_concordance": {
            "pearson_all_tested": round(pearson, 4),
            "spearman_all_tested": round(spearman, 4),
            "pearson_true_de": round(pearson_de, 4),
            "spearman_true_de": round(spearman_de, 4),
            "n_genes_scored": int(len(lfc_frame)),
            "n_true_de_scored": int(len(lfc_de_frame)),
        },
        "ranking_concordance": {
            "top_k": int(top_k),
            "top_k_jaccard": round(top_k_jaccard, 4),
            "top_k_truth_recall": round(top_k_truth_recall, 4),
            "rank_biased_overlap_p0.9": round(rbo, 4),
        },
        "significant_call_concordance": {
            "recall": round(recall, 4),
            "precision": round(precision, 4),
            "empirical_fdr": round(empirical_fdr, 4),
            "n_true_de": len(true_de),
            "n_called": len(called),
            "n_true_positive": len(tp),
            "n_false_positive": len(fp),
        },
    }
    axis_pass = {
        "fdr_calibration": (
            negative_control.status == "pass"
            and negative_control.false_positive_rate <= max_null_fpr
        ),
        "lfc_concordance": spearman_de >= min_lfc_spearman,
        "ranking_concordance": top_k_jaccard >= min_top_k_jaccard,
        "significant_call_concordance": (
            recall >= min_recall and empirical_fdr <= max_empirical_fdr
        ),
    }
    status = "pass" if all(axis_pass.values()) else "fail"
    return {
        "status": status,
        "benchmark": "A1_bulk_de",
        "benchmark_version": "v1",
        "scope": "preliminary_synthetic_truth",
        "method_under_test": "ARIA bulk RNA DESeq2 path (_run_deseq2, apeGLM-enabled)",
        "comparison": {"numerator": "COND_B", "denominator": "COND_A"},
        "dataset": dataset.params,
        "thresholds": {"padj": padj_max, "lfc_threshold": lfc_min},
        "tolerances": tolerances,
        "axis_pass": axis_pass,
        "axes": axes,
        "messages": [
            "A1 preliminary validates ARIA's bulk DE path against synthetic truth; "
            "external DESeq2/edgeR/limma comparators remain assigned to aria-bench-env."
        ],
    }


def sweep_bulk_de_lfc_threshold(
    dataset: SyntheticBulkDEDataset,
    *,
    thresholds: tuple[float, ...] = (0.0, 0.25, 0.5, 1.0),
    padj_max: float = 0.05,
    policy_threshold: float = 0.5,
) -> dict[str, Any]:
    """Quantify the recall/precision frontier of ARIA's Wald ``lfcThreshold``.

    Runs the real bulk DE path (``_run_deseq2``, apeGLM-enabled) on a single
    synthetic dataset at several Wald lfcThreshold values. ``lfc_thr=0.0``
    reproduces the standard DESeq2 null (H0: LFC = 0); higher thresholds test
    H0: ``|LFC| <= thr`` and deliberately trade recall for precision.

    This isolates how much of any recall gap versus standard DESeq2/edgeR/limma
    is ARIA's *effect-size policy* (a user-controlled choice) rather than an
    engine difference: the ``lfc_threshold == 0`` point is the matched-null
    DESeq2-equivalence reference, computed here (not hardcoded), so a report can
    show the frontier instead of a single conservative recall number.
    """
    from aria.scripts.rna_bulk_de import _run_deseq2

    true_de = set(map(str, dataset.de_genes))
    null = set(map(str, dataset.null_genes))
    n_true = len(true_de)

    points: list[dict[str, Any]] = []
    for thr in thresholds:
        thr = float(thr)
        res, _w = _run_deseq2(
            dataset.counts, dataset.metadata, "condition", "COND_B", "COND_A",
            padj_thr=padj_max, lfc_thr=thr, lfc_shrink=True,
        )
        if res.get("status") != "success":
            points.append({
                "lfc_threshold": thr,
                "status": "error",
                "error_type": res.get("error_type"),
                "details": res.get("details", ""),
            })
            continue
        called = set(map(str, res.get("sig_genes", []) or []))
        tp = called & true_de
        fp = called & null
        points.append({
            "lfc_threshold": thr,
            "status": "success",
            "n_called": len(called),
            "n_true_positive": len(tp),
            "n_false_positive": len(fp),
            "recall": round(len(tp) / max(n_true, 1), 4),
            "precision": round(len(tp) / max(len(called), 1), 4),
            "empirical_fdr": round(len(fp) / max(len(called), 1), 4),
            "is_matched_null": thr == 0.0,
            "is_policy_default": thr == float(policy_threshold),
        })
    return {
        "description": (
            "Wald lfcThreshold recall/precision frontier on one synthetic truth. "
            "lfc_threshold=0 is the matched-null DESeq2 equivalence reference; "
            "higher thresholds are ARIA's deliberate effect-size policy."
        ),
        "padj_max": padj_max,
        "policy_threshold": float(policy_threshold),
        "n_true_de": n_true,
        "points": points,
    }


def write_bulk_de_a1_figure(manifest: dict[str, Any], path: str) -> str:
    """Write a dependency-free SVG summary for preliminary Fig. 1."""
    from pathlib import Path

    axes = manifest.get("axes", {})
    vals = [
        ("FDR calibration", 1.0 - float(axes.get("fdr_calibration", {}).get("false_positive_rate", 1.0))),
        ("LFC Spearman", float(axes.get("lfc_concordance", {}).get("spearman_true_de", 0.0))),
        ("Top-k Jaccard", float(axes.get("ranking_concordance", {}).get("top_k_jaccard", 0.0))),
        ("Call recall", float(axes.get("significant_call_concordance", {}).get("recall", 0.0))),
    ]
    width, height = 760, 390
    left, top = 190, 70
    bar_w, bar_h, gap = 430, 42, 26
    rows = []
    for i, (label, value) in enumerate(vals):
        y = top + i * (bar_h + gap)
        v = max(0.0, min(1.0, value))
        fill = "#2F6F73" if v >= 0.5 else "#A23B3B"
        rows.append(
            f'<text x="28" y="{y + 27}" font-size="17" fill="#1f2933">{label}</text>'
            f'<rect x="{left}" y="{y}" width="{bar_w}" height="{bar_h}" fill="#e8edf0"/>'
            f'<rect x="{left}" y="{y}" width="{bar_w * v:.1f}" height="{bar_h}" fill="{fill}"/>'
            f'<text x="{left + bar_w + 18}" y="{y + 27}" font-size="17" fill="#1f2933">{v:.2f}</text>'
        )
    status = manifest.get("status", "unknown").upper()
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="#ffffff"/>'
        '<text x="28" y="34" font-size="22" font-weight="700" fill="#111827">'
        'Fig. 1 preliminary: A1 bulk DE validation</text>'
        f'<text x="630" y="34" font-size="18" font-weight="700" fill="#2F6F73">{status}</text>'
        '<text x="28" y="366" font-size="13" fill="#4b5563">'
        'Synthetic truth; external R comparators are separate aria-bench-env work.</text>'
        + "".join(rows)
        + "</svg>"
    )
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg, encoding="utf-8")
    return str(out)


def run_bulk_de_a1_benchmark(
    *,
    seed: int = 11,
    quick: bool = False,
    output_dir: str | None = None,
    manifest_name: str = "a1_bulk_de_manifest.json",
    figure_name: str = "fig1_a1_bulk_de.svg",
) -> dict[str, Any]:
    """Execute Benchmark A1's preliminary synthetic bulk-DE lane."""
    import json
    from pathlib import Path
    from aria.scripts.rna_bulk_de import _run_deseq2
    from aria.version import __version__, collect_version_metadata

    sim_kwargs = (
        dict(n_genes=600, n_de=60, replicates_per_condition=6)
        if quick else
        dict(n_genes=1000, n_de=120, replicates_per_condition=6)
    )
    n_perms = 2 if quick else 3
    dataset = simulate_bulk_dataset(seed=seed, **sim_kwargs)
    de_result, warnings = _run_deseq2(
        dataset.counts, dataset.metadata, "condition", "COND_B", "COND_A",
        padj_thr=0.05, lfc_thr=0.5, lfc_shrink=True,
    )
    negative = run_bulk_de_negative_control(
        dataset=dataset, seed=seed, n_permutations=n_perms,
        nominal_alpha=0.05, max_false_positive_rate=0.05, lfc_min=0.5,
    )
    manifest = score_bulk_de_a1(dataset, de_result, negative)
    frontier = sweep_bulk_de_lfc_threshold(
        dataset,
        thresholds=(0.0, 0.5) if quick else (0.0, 0.25, 0.5, 1.0),
        padj_max=0.05,
        policy_threshold=0.5,
    )
    manifest["lfc_threshold_frontier"] = frontier
    matched = next(
        (p for p in frontier["points"]
         if p.get("is_matched_null") and p.get("status") == "success"),
        None,
    )
    if matched is not None:
        manifest.setdefault("messages", []).append(
            "lfc_threshold_frontier: at lfc_threshold=0 (matched DESeq2 null) "
            f"ARIA recall={matched['recall']} / empirical_fdr={matched['empirical_fdr']} "
            f"(n_called={matched['n_called']}); the scored axes use ARIA's default "
            "lfcThreshold=0.5 effect-size policy. The recall difference is policy, "
            "not engine."
        )
    manifest.update({
        "aria_version": __version__,
        "provenance": collect_version_metadata(),
        "warnings": warnings,
    })
    if output_dir:
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        manifest_path = outdir / manifest_name
        figure_path = outdir / figure_name
        write_bulk_de_a1_figure(manifest, str(figure_path))
        manifest["artifacts"] = {
            "manifest_json": str(manifest_path),
            "figure_svg": str(figure_path),
        }
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


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


# ── W-CALIB: spike-in dose-response (effect-size calibration ladder) ─────────
#
# Recovery + negative controls answer "finds true effects?" and "quiet under the
# null?". Spike-ins add the dose-response question an ERCC ladder answers: across
# a ladder of KNOWN |log2FC| levels, does detection rise monotonically from ~0 at
# level 0 (true nulls) to high at the strongest level, are the level-0 spike-ins
# kept below alpha, and are the estimated effect sizes close to the truth?


@dataclass
class SpikeInDataset:
    """A bulk matrix with spike-in genes at a ladder of known fold-changes."""
    counts: Any                              # DataFrame (genes x samples)
    metadata: Any                            # DataFrame (samples x design)
    spike_true_log2fc: dict[str, float]      # spike gene -> signed true log2fc
    spike_level: dict[str, float]            # spike gene -> |true log2fc| level
    null_genes: list[str]                    # background null genes
    levels: list[float]
    params: dict[str, Any] = field(default_factory=dict)


@dataclass
class SpikeInResult:
    status: str                              # "pass" | "fail" | "error"
    levels: list[float]
    detection_rate_by_level: dict[str, float]  # str(level) -> fraction significant
    lfc_mae: float                           # mean |estimated - true| log2fc on spike-ins
    null_spike_fpr: float                    # detection rate at level 0 (true nulls)
    nominal_alpha: float
    n_spike_per_level: int
    tolerances: dict[str, float]
    messages: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "levels": self.levels,
            "detection_rate_by_level": {
                k: round(v, 4) for k, v in self.detection_rate_by_level.items()
            },
            "lfc_mae": round(self.lfc_mae, 4),
            "null_spike_fpr": round(self.null_spike_fpr, 4),
            "nominal_alpha": self.nominal_alpha,
            "n_spike_per_level": self.n_spike_per_level,
            "tolerances": self.tolerances,
            "messages": self.messages,
        }


def simulate_spike_in_bulk_dataset(
    *,
    levels: tuple[float, ...] = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0),
    genes_per_level: int = 15,
    n_background: int = 1500,
    replicates_per_condition: int = 6,
    dispersion: float = 0.2,
    seed: int = 17,
) -> SpikeInDataset:
    """Simulate a bulk matrix with a ladder of known-fold-change spike-in genes.

    ``genes_per_level`` spike-in genes are planted at each ``|log2FC|`` level (sign
    random, except level 0 which is a true null). The remaining ``n_background``
    genes are null. Negative-binomial counts at bulk depth; deterministic for a
    seed. Neutral labels only (ADR-011).
    """
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)

    n_spike = genes_per_level * len(levels)
    n_genes = n_background + n_spike
    gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    base_mean = np.clip(np.exp(rng.normal(3.0, 1.0, size=n_genes)), 5.0, None)

    # The last n_spike genes are the spike-ins, ordered by level.
    spike_true_log2fc: dict[str, float] = {}
    spike_level: dict[str, float] = {}
    log2fc = np.zeros(n_genes, dtype=float)
    si = n_background
    for lvl in levels:
        for _ in range(genes_per_level):
            sign = 1.0 if lvl == 0.0 else float(rng.choice([-1.0, 1.0]))
            val = sign * float(lvl)
            log2fc[si] = val
            spike_true_log2fc[gene_names[si]] = val
            spike_level[gene_names[si]] = float(lvl)
            si += 1
    cond_b_fc = np.power(2.0, log2fc)

    sample_cols, sample_names, conds = [], [], []
    for cond in ("COND_A", "COND_B"):
        fc = cond_b_fc if cond == "COND_B" else np.ones(n_genes)
        for r in range(replicates_per_condition):
            sf = float(np.exp(rng.normal(0.0, 0.15)))
            gene_mean = base_mean * fc * sf
            shape = 1.0 / dispersion
            gamma = rng.gamma(shape=shape, scale=gene_mean * dispersion)
            sample_cols.append(rng.poisson(gamma).astype(np.int64))
            sample_names.append(f"{cond}_r{r}")
            conds.append(cond)

    counts = pd.DataFrame(
        np.array(sample_cols).T, index=gene_names, columns=sample_names)
    metadata = pd.DataFrame({"condition": conds}, index=sample_names)
    null_genes = [gene_names[i] for i in range(n_background)]

    return SpikeInDataset(
        counts=counts,
        metadata=metadata,
        spike_true_log2fc=spike_true_log2fc,
        spike_level=spike_level,
        null_genes=null_genes,
        levels=list(levels),
        params={
            "levels": list(levels), "genes_per_level": genes_per_level,
            "n_background": n_background,
            "replicates_per_condition": replicates_per_condition,
            "dispersion": dispersion, "seed": seed,
        },
    )


def run_bulk_de_spike_in(
    dataset: SpikeInDataset | None = None,
    *,
    nominal_alpha: float = 0.05,
    max_null_fpr: float = 0.1,
    min_top_detection: float = 0.5,
    max_lfc_mae: float = 1.0,
    lfc_min: float = 0.5,
    lfc_shrink: bool = True,
    seed: int = 17,
    **sim_kwargs,
) -> SpikeInResult:
    """Run ARIA's real bulk DE on a spike-in ladder and score effect-size calibration.

    Measures per-level detection rate (a dose-response curve), the level-0 false
    positive rate (true-null spike-ins), and the mean absolute error between the
    apeGLM-shrunken estimated log2FC and the known truth on the spike-ins.
    ``status`` is "pass" only when the null spike-ins stay below ``max_null_fpr``,
    the strongest level clears ``min_top_detection``, and the effect-size MAE is
    within ``max_lfc_mae``. Requires pydeseq2.
    """
    from aria.scripts.rna_bulk_de import _run_deseq2

    if dataset is None:
        dataset = simulate_spike_in_bulk_dataset(seed=seed, **sim_kwargs)

    tolerances = {
        "max_null_fpr": max_null_fpr,
        "min_top_detection": min_top_detection,
        "max_lfc_mae": max_lfc_mae,
        "nominal_alpha": nominal_alpha,
    }
    levels = dataset.levels
    n_spike_per_level = sum(1 for v in dataset.spike_level.values()
                            if v == levels[0]) if levels else 0

    result, _warnings = _run_deseq2(
        dataset.counts, dataset.metadata, "condition", "COND_B", "COND_A",
        padj_thr=nominal_alpha, lfc_thr=lfc_min, lfc_shrink=lfc_shrink,
    )
    if result.get("status") != "success":
        return SpikeInResult(
            status="error", levels=levels, detection_rate_by_level={},
            lfc_mae=float("nan"), null_spike_fpr=1.0, nominal_alpha=nominal_alpha,
            n_spike_per_level=n_spike_per_level, tolerances=tolerances,
            messages=[f"bulk DE did not succeed: {result.get('error_type')} "
                      f"{result.get('details', '')}"],
        )

    called = set(result.get("sig_genes", []) or [])
    results_df = result.get("results")
    est_log2fc = {}
    if results_df is not None:
        est_log2fc = {str(g): float(v)
                      for g, v in results_df["log2FoldChange"].items()}

    # Detection rate per level.
    detection_rate_by_level: dict[str, float] = {}
    by_level: dict[float, list[str]] = {}
    for gene, lvl in dataset.spike_level.items():
        by_level.setdefault(lvl, []).append(gene)
    for lvl in levels:
        genes = by_level.get(lvl, [])
        n_called = sum(1 for g in genes if g in called)
        detection_rate_by_level[str(lvl)] = n_called / max(len(genes), 1)

    null_spike_fpr = detection_rate_by_level.get(str(levels[0]), 0.0) if levels else 0.0

    # Effect-size accuracy (shrunken estimate vs truth) over spike-ins present.
    errs = [abs(est_log2fc[g] - true)
            for g, true in dataset.spike_true_log2fc.items() if g in est_log2fc]
    lfc_mae = float(sum(errs) / len(errs)) if errs else float("nan")

    top_detection = detection_rate_by_level.get(str(levels[-1]), 0.0) if levels else 0.0
    ok = (
        null_spike_fpr <= max_null_fpr
        and top_detection >= min_top_detection
        and lfc_mae <= max_lfc_mae
    )
    messages = [
        f"null_spike_fpr={null_spike_fpr:.3f} (<= {max_null_fpr}); "
        f"top_detection={top_detection:.3f} (>= {min_top_detection}); "
        f"lfc_mae={lfc_mae:.3f} (<= {max_lfc_mae})",
    ]
    return SpikeInResult(
        status="pass" if ok else "fail",
        levels=levels,
        detection_rate_by_level=detection_rate_by_level,
        lfc_mae=lfc_mae,
        null_spike_fpr=null_spike_fpr,
        nominal_alpha=nominal_alpha,
        n_spike_per_level=n_spike_per_level,
        tolerances=tolerances,
        messages=messages,
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
        spike_kw = dict(levels=(0.0, 1.0, 2.0, 3.0), genes_per_level=10,
                        n_background=500, replicates_per_condition=6)
        n_perms = 3
        bulk_min_recall = pb_min_recall = 0.4
    else:
        bulk_kw = dict(n_genes=1000, n_de=120, replicates_per_condition=6)
        pb_kw = dict(n_genes=1200, n_de=120, donors_per_condition=6, cells_per_donor=80)
        spike_kw = dict(levels=(0.0, 0.5, 1.0, 1.5, 2.0, 3.0), genes_per_level=15,
                        n_background=1000, replicates_per_condition=6)
        n_perms = 3
        bulk_min_recall = pb_min_recall = 0.5

    bulk_recovery = run_bulk_de_benchmark(seed=seed, min_recall=bulk_min_recall, **bulk_kw)
    bulk_neg = run_bulk_de_negative_control(seed=seed, n_permutations=n_perms, **bulk_kw)
    pb_recovery = run_pseudobulk_de_benchmark(seed=seed, min_recall=pb_min_recall, **pb_kw)
    pb_neg = run_pseudobulk_de_negative_control(seed=seed, n_permutations=n_perms, **pb_kw)
    spike = run_bulk_de_spike_in(seed=seed, **spike_kw)

    paths = {
        "bulk": {
            "recovery": bulk_recovery.as_dict(),
            "negative_control": bulk_neg.as_dict(),
            "spike_in": spike.as_dict(),
        },
        "pseudobulk": {
            "recovery": pb_recovery.as_dict(),
            "negative_control": pb_neg.as_dict(),
        },
    }
    all_results = [bulk_recovery, bulk_neg, pb_recovery, pb_neg, spike]
    if all(r.status == "pass" for r in all_results):
        status = "pass"
    elif any(r.status == "error" for r in all_results):
        status = "error"
    else:
        status = "fail"

    top_level = str(spike.levels[-1]) if spike.levels else ""
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
            "bulk_spike_null_fpr": spike.null_spike_fpr,
            "bulk_spike_top_detection": spike.detection_rate_by_level.get(top_level, 0.0),
            "bulk_spike_lfc_mae": spike.lfc_mae,
        },
    }
