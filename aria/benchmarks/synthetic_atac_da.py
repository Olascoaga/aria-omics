"""Typed scATAC differential-accessibility (DA) calibration hook — SCAFFOLD (P3-2).

P3-2 (pre-4.6 polish): leave the calibration suite with a TYPED SLOT for scATAC
differential accessibility against a simulated ground truth, so the v4.6
chromatin DA implementation is born with a numerical regression — recall on
truly-differential peaks and empirical FDR on null peaks — instead of only flow
tests, exactly as bulk/pseudobulk DE already have (X6 / W-CALIB).

This is a SCAFFOLD. The ground-truth simulator and the recovery scoring are real
and usable now, but ARIA has NO validated scATAC DA backend yet. So
``run_atac_da_benchmark`` REQUIRES the caller to inject the DA function (``da_fn``);
with none it raises ``NotImplementedError`` rather than fabricating calibration
metrics (same no-fabrication contract as the integration scaffolds). When v4.6
lands a real DA caller, wire it as ``da_fn`` and this becomes a live regression —
no other change needed.

Dependency-light (numpy / pandas / anndata). Neutral synthetic labels only
(ADR-011): PEAK_#####, COND_A / COND_B, cell c#####. It is a generated
ground-truth benchmark, not a named golden dataset.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

from aria.benchmarks.synthetic_de import DEBenchmarkResult


@dataclass
class SyntheticATACDADataset:
    """A simulated scATAC accessibility matrix plus its ground truth."""
    adata: Any                              # AnnData (cells x peaks), integer accessibility counts in X
    da_peaks: dict[str, str]                # peak -> "up" | "down" (truth in COND_B vs COND_A)
    null_peaks: list[str]                   # peaks with no true effect
    params: dict[str, Any] = field(default_factory=dict)

    @property
    def n_da(self) -> int:
        return len(self.da_peaks)


# The contract v4.6 fills: given the dataset, return the set/list of peaks the
# real scATAC DA path calls significant (in COND_B vs COND_A).
ATACDACaller = Callable[[SyntheticATACDADataset], "set[str] | list[str]"]


def simulate_atac_da_dataset(
    *,
    n_peaks: int = 3000,
    n_da: int = 300,
    n_cells_per_condition: int = 400,
    base_rate: float = 0.06,
    min_abs_log2fc: float = 1.0,
    max_abs_log2fc: float = 2.5,
    seed: int = 11,
) -> SyntheticATACDADataset:
    """Simulate a two-condition scATAC accessibility matrix (cells x peaks).

    scATAC accessibility is sparse and near-binary, so each (peak, cell) is drawn
    from a Poisson with a low per-peak base accessibility rate; ``n_da`` peaks get
    a fold change in their rate in COND_B and the rest are null. Deterministic for
    a given seed. Neutral labels only (ADR-011).
    """
    import numpy as np
    import pandas as pd
    import anndata as ad

    rng = np.random.default_rng(seed)

    peak_names = [f"PEAK_{i:05d}" for i in range(n_peaks)]
    # Per-peak base accessibility rate (mean counts per cell), heavy-tailed-ish
    # but bounded so peaks are sparse like real scATAC.
    base = np.clip(rng.lognormal(mean=np.log(base_rate), sigma=0.6, size=n_peaks),
                   1e-3, 1.0)

    da_idx = rng.choice(n_peaks, size=n_da, replace=False)
    da_mask = np.zeros(n_peaks, dtype=bool)
    da_mask[da_idx] = True
    log2fc = np.zeros(n_peaks, dtype=float)
    signs = rng.choice([-1.0, 1.0], size=n_da)
    mags = rng.uniform(min_abs_log2fc, max_abs_log2fc, size=n_da)
    log2fc[da_idx] = signs * mags
    cond_b_fc = np.power(2.0, log2fc)

    rows, conds, cell_names = [], [], []
    for cond in ("COND_A", "COND_B"):
        rate_vec = base * (cond_b_fc if cond == "COND_B" else np.ones(n_peaks))
        for c in range(n_cells_per_condition):
            # Per-cell depth factor for realistic library-size variation.
            depth = float(np.exp(rng.normal(0.0, 0.25)))
            rows.append(rng.poisson(rate_vec * depth).astype(np.int64))
            conds.append(cond)
            cell_names.append(f"{cond}_c{c:05d}")

    X = np.asarray(rows, dtype=np.int64)
    obs = pd.DataFrame({"condition": conds}, index=cell_names)
    var = pd.DataFrame(index=peak_names)
    adata = ad.AnnData(X=X, obs=obs, var=var)

    da_peaks = {peak_names[i]: ("up" if log2fc[i] > 0 else "down") for i in da_idx}
    null_peaks = [peak_names[i] for i in range(n_peaks) if not da_mask[i]]

    return SyntheticATACDADataset(
        adata=adata,
        da_peaks=da_peaks,
        null_peaks=null_peaks,
        params={
            "n_peaks": n_peaks, "n_da": n_da,
            "n_cells_per_condition": n_cells_per_condition,
            "base_rate": base_rate,
            "min_abs_log2fc": min_abs_log2fc,
            "max_abs_log2fc": max_abs_log2fc,
            "seed": seed,
        },
    )


def run_atac_da_benchmark(
    dataset: SyntheticATACDADataset | None = None,
    *,
    da_fn: ATACDACaller | None = None,
    min_recall: float = 0.5,
    max_empirical_fdr: float = 0.2,
    **sim_kwargs,
) -> DEBenchmarkResult:
    """Score a scATAC DA caller against the simulated ground truth.

    ``da_fn`` is the real differential-accessibility function (injected by v4.6):
    given the dataset it returns the peaks it calls significant. Recovery is then
    scored exactly like the DE benchmarks (recall on true-DA peaks, empirical FDR
    on null peaks). With no ``da_fn`` this raises ``NotImplementedError`` — ARIA
    has no validated scATAC DA backend yet and refuses to fabricate calibration
    metrics. This is the typed slot v4.6 fills.
    """
    if dataset is None:
        dataset = simulate_atac_da_dataset(**sim_kwargs)

    tolerances = {"min_recall": min_recall, "max_empirical_fdr": max_empirical_fdr}

    if da_fn is None:
        raise NotImplementedError(
            "scATAC differential accessibility is a v4.6 scaffold: no validated DA "
            "backend is wired yet. Inject the real DA caller as `da_fn` once v4.6 "
            "implements it; ARIA refuses to fabricate calibration metrics."
        )

    called = set(da_fn(dataset) or [])
    true_da = set(dataset.da_peaks)
    null = set(dataset.null_peaks)
    tp = called & true_da
    fp = called & null
    recall = len(tp) / max(len(true_da), 1)
    empirical_fdr = len(fp) / max(len(called), 1)

    status = "pass" if (recall >= min_recall and empirical_fdr <= max_empirical_fdr) else "fail"
    messages = [
        f"recall={recall:.3f} (>= {min_recall}); "
        f"empirical_fdr={empirical_fdr:.3f} (<= {max_empirical_fdr}); "
        f"called={len(called)}, true_da={len(true_da)}",
    ]

    return DEBenchmarkResult(
        status=status,
        recall=recall,
        empirical_fdr=empirical_fdr,
        n_true_de=len(true_da),
        n_called=len(called),
        n_true_positive=len(tp),
        n_false_positive=len(fp),
        tolerances=tolerances,
        messages=messages,
    )
