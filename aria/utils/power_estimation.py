"""Approximate RNA-seq power estimates for ARIA reports.

These helpers intentionally use a closed-form Wald approximation rather than
simulation. DESeq2 models counts with a negative binomial likelihood
(Love et al. 2014); for report-level planning language we approximate the
standard error of a log fold-change from group size, mean expression, and
dispersion. The result is a transparent design diagnostic, not a substitute
for a dedicated power analysis.
"""

from __future__ import annotations

from math import isfinite, log
from statistics import NormalDist


def _as_pair(n_per_group) -> tuple[float, float]:
    def number(value):
        return float("nan") if isinstance(value, bool) else float(value)

    if isinstance(n_per_group, dict):
        vals = list(n_per_group.values())
    elif isinstance(n_per_group, (list, tuple)):
        vals = list(n_per_group)
    else:
        vals = [n_per_group]
    if not vals:
        return 0.0, 0.0
    if len(vals) >= 2:
        return number(vals[0]), number(vals[1])
    n = number(vals[0])
    return n, n


def bulk_power_estimate(n_per_group,
                        mean_expression: float = 100.0,
                        dispersion: float = 0.1,
                        target_log2fc: float = 0.5,
                        alpha: float = 0.05) -> float:
    """Approximate two-sided Wald power for a bulk RNA-seq contrast."""
    n1, n2 = _as_pair(n_per_group)
    if not isfinite(n1) or not isfinite(n2) or n1 <= 1 or n2 <= 1:
        return 0.0
    if any(isinstance(v, bool) for v in (mean_expression, dispersion, target_log2fc, alpha)):
        return 0.0
    mean_expression = float(1.0 if mean_expression is None else mean_expression)
    dispersion = float(0.0 if dispersion is None else dispersion)
    mean_expression = max(mean_expression, 1.0) if isfinite(mean_expression) else 1.0
    dispersion = max(dispersion, 1e-8) if isfinite(dispersion) else 0.1
    target = abs(float(0.0 if target_log2fc is None else target_log2fc))
    alpha = float(0.05 if alpha is None else alpha)
    if not isfinite(target) or not isfinite(alpha) or alpha <= 0:
        return 0.0
    alpha = min(max(alpha, 1e-12), 0.999999)

    variance_log = (1.0 / n1 + 1.0 / n2) * (
        (1.0 / mean_expression) + dispersion
    )
    se_log2 = (variance_log ** 0.5) / log(2.0)
    if not isfinite(se_log2) or se_log2 <= 0:
        return 0.0

    nd = NormalDist()
    z_alpha = nd.inv_cdf(1.0 - alpha / 2.0)
    z_effect = target / se_log2
    power = nd.cdf(-z_alpha - z_effect) + (1.0 - nd.cdf(z_alpha - z_effect))
    if not isfinite(power):
        return 0.0
    return round(float(min(max(power, 0.0), 1.0)), 4)


def pseudobulk_power_estimate(n_per_group,
                              dispersion_estimate: float,
                              target_log2fc: float,
                              alpha: float,
                              mean_expression: float = 100.0) -> float:
    """Approximate two-sided Wald power for one scRNA pseudobulk block."""
    return bulk_power_estimate(
        n_per_group=n_per_group,
        mean_expression=mean_expression,
        dispersion=dispersion_estimate,
        target_log2fc=target_log2fc,
        alpha=alpha,
    )

