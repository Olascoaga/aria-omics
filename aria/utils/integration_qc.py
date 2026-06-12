"""Integration quality red-flags (audit item X8).

Harmony (and any batch integration) can *overcorrect*: it mixes batches so well
that genuine biological structure is destroyed. ARIA used to report the batch
silhouette as a passive number; this turns the danger signals into explicit,
actionable findings.

Convention (matches `rna_integration.py`): the silhouette is computed on BATCH
labels, where LOWER (more negative) = better mixing. The cluster silhouette
(from `rna_clustering.py`) is computed on biological clusters, where HIGHER =
cleaner structure. Overcorrection = batches well-mixed but cluster structure
collapsed.

Pure / dependency-light so it is fully unit-testable.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class QCIssue:
    severity: str           # "warning" | "blocking"
    check: str
    message: str
    recommendation: str

    def as_dict(self) -> dict[str, str]:
        return {
            "severity": self.severity,
            "check": self.check,
            "message": self.message,
            "recommendation": self.recommendation,
        }


def assess_integration_quality(
    silhouette_before: float | None,
    silhouette_after: float | None,
    cluster_silhouette: float | None = None,
    *,
    batch_high: float = 0.10,
    cluster_low: float = 0.0,
) -> dict[str, Any]:
    """Flag under-mixing and overcorrection from integration silhouettes.

    Args:
        silhouette_before/after: BATCH-label silhouette before/after integration
            (lower = better mixing).
        cluster_silhouette: biological-cluster silhouette on the integrated
            embedding (higher = cleaner structure), if available.
        batch_high: batch silhouette above which batches are still separated
            (under-correction).
        cluster_low: cluster silhouette at/below which biological structure is
            considered lost (overcorrection signal).

    Returns: {"status": "clean"|"warnings", "issues": [QCIssue.as_dict(), ...],
              "metrics": {...}}.
    """
    issues: list[QCIssue] = []

    sb = _as_float(silhouette_before)
    sa = _as_float(silhouette_after)
    cs = _as_float(cluster_silhouette)

    # Integration made batch separation WORSE (mixing regressed).
    if sb is not None and sa is not None and sa > sb + 1e-9:
        issues.append(QCIssue(
            "warning",
            "integration_worsened_mixing",
            (f"Batch silhouette rose after integration "
             f"({sb:.3f} -> {sa:.3f}); lower means better mixing, so "
             f"integration did not improve batch overlap."),
            "Check the batch column and whether integration was needed; "
            "consider an alternative method or no integration.",
        ))

    # Batch effect still dominates after integration (under-correction).
    if sa is not None and sa > batch_high:
        issues.append(QCIssue(
            "warning",
            "residual_batch_effect",
            (f"Batch silhouette after integration is {sa:.3f} (> {batch_high}); "
             f"batches remain separated, so downstream structure may still be "
             f"batch-driven."),
            "Inspect the embedding colored by batch; consider stronger "
            "integration or adding batch as a covariate in DE.",
        ))

    # Biological cluster structure collapsed (overcorrection).
    if cs is not None and cs <= cluster_low:
        sev = "warning"
        msg = (f"Cluster silhouette on the integrated embedding is {cs:.3f} "
               f"(<= {cluster_low}); cells are not closer to their own cluster "
               f"than to others, a sign of overcorrection or absent structure.")
        rec = ("Verify integration did not erase biology: compare marker "
               "separation before/after, lower integration strength, or skip "
               "integration if batches are biologically distinct.")
        # N-INT1: strong mixing AND collapsed structure is the classic
        # overcorrection signature — escalate from a soft warning to blocking so
        # the report does not bury a destroyed-biology run among minor warnings.
        if sb is not None and sa is not None and sa < sb - 0.05:
            sev = "blocking"
            msg = ("Batches were strongly mixed but " + msg[0].lower() + msg[1:]
                   + " This pattern is the classic overcorrection signature.")
        issues.append(QCIssue(sev, "possible_overcorrection", msg, rec))

    status = "warnings" if issues else "clean"
    return {
        "status": status,
        "issues": [i.as_dict() for i in issues],
        "metrics": {
            "silhouette_before": sb,
            "silhouette_after": sa,
            "cluster_silhouette": cs,
        },
    }


def _as_float(x) -> float | None:
    if x is None:
        return None
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    if f != f:  # NaN
        return None
    return f
