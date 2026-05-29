"""Annotation-coherence red-flags for reused obs labels (audit item X9).

When ARIA reuses pre-existing `obs` annotations (e.g. `obs['subclass']`) and
skips Leiden/CellTypist, it trusts the labels blindly. A label from another
species or atlas, or an over-merged/noise group, would then propagate into DE
and the report unchecked.

This module does NOT use a hardcoded marker -> cell-type dictionary (that would
violate ADR-011 and could not generalize). It is purely data-driven: it asks
whether each reused label has a DISTINCT transcriptional signature (its own top
marker genes). Labels with no/weak distinct markers are flagged as low
confidence. When markers could not be computed at all (e.g. a Seurat-scaled X
with no usable log-norm layer), it flags the annotation as UNVERIFIED rather
than silently trusting it.

It cannot detect a transcriptionally-distinct-but-misnamed cluster without an
external reference — that limitation is stated explicitly, not hidden.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class QCIssue:
    severity: str
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


def assess_annotation_coherence(
    top_markers: dict[str, list] | None,
    cluster_sizes: dict[str, int] | None = None,
    *,
    reused: bool = True,
    markers_verified: bool = True,
    min_markers: int = 3,
) -> dict[str, Any]:
    """Assess whether reused annotation labels are transcriptionally distinct.

    Args:
        top_markers: {label: [marker_gene, ...]} from a rank-genes step.
        cluster_sizes: {label: n_cells} (for context in messages).
        reused: True when labels came from existing obs (skipped Leiden).
        markers_verified: False when markers could not be computed on a
            suitable expression layer (e.g. scaled-only X).
        min_markers: a label with fewer distinct top markers than this is
            considered to lack a clear signature.

    Returns: {"status": "clean"|"warnings"|"unverified",
              "per_label": {...}, "issues": [...]}.
    """
    top_markers = top_markers or {}
    cluster_sizes = cluster_sizes or {}
    issues: list[QCIssue] = []

    n_markers_per_label = {
        str(label): len([m for m in (markers or []) if str(m).strip()
                         and str(m).lower() != "nan"])
        for label, markers in top_markers.items()
    }

    # Could not verify at all (reuse path with no computed markers).
    if not markers_verified or (n_markers_per_label and
                                all(v == 0 for v in n_markers_per_label.values())):
        if reused:
            issues.append(QCIssue(
                "warning",
                "annotation_unverified",
                ("Reused obs annotations were not marker-verified (no suitable "
                 "log-normalized expression was available to recompute marker "
                 "genes). Labels are trusted as provided; ARIA cannot confirm "
                 "each label reflects a distinct transcriptional state."),
                "If the labels come from another atlas/species, validate a few "
                "canonical markers manually, or let ARIA re-cluster/annotate.",
            ))
        return {
            "status": "unverified" if reused else "clean",
            "per_label": {k: {"n_markers": v, "distinct": None}
                          for k, v in n_markers_per_label.items()},
            "issues": [i.as_dict() for i in issues],
        }

    # Marker-verified: flag labels lacking a distinct signature.
    per_label: dict[str, Any] = {}
    weak_labels: list[str] = []
    for label, n in n_markers_per_label.items():
        distinct = n >= min_markers
        per_label[label] = {"n_markers": n, "distinct": distinct}
        if not distinct:
            weak_labels.append(label)

    if weak_labels:
        sizes = ", ".join(
            f"{lbl} (n={cluster_sizes.get(lbl, '?')}, markers={n_markers_per_label[lbl]})"
            for lbl in weak_labels
        )
        issues.append(QCIssue(
            "warning",
            "labels_without_distinct_markers",
            (f"{len(weak_labels)} annotation label(s) lack a distinct marker "
             f"signature (< {min_markers} significant markers): {sizes}. "
             f"They may be noise, doublets, or over-merged groups rather than "
             f"coherent cell types."),
            "Inspect these labels; consider merging, excluding, or re-annotating "
            "them before drawing cell-type-specific conclusions.",
        ))

    return {
        "status": "warnings" if issues else "clean",
        "per_label": per_label,
        "issues": [i.as_dict() for i in issues],
    }
