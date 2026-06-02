"""Ambient-RNA contamination red-flags for scRNA (audit item P1-4).

Ambient ("soup") RNA — cell-free transcripts captured alongside each cell —
leaks the most abundant genes into every droplet, so a highly expressed gene
from one population shows up as a "marker" everywhere. ARIA does not run
SoupX/decontX by default (it is an opt-in step), so this assessor at least
DETECTS the signature and warns; it never corrects.

It is data-driven (ADR-011: no hardcoded gene or cell-type lists). It reuses the
per-cluster `top_markers` ARIA already computes during clustering and measures
how many clusters list each gene as a top marker. A gene that is "top" in a
large fraction of clusters is non-specific (ambient/soup-like or housekeeping).
When many top-marker slots are filled by such ubiquitous genes, ARIA flags
possible ambient contamination and recommends the optional decontamination
step.

Limitation (stated, not hidden): ubiquitous markers can also reflect genuine
shared biology (e.g. broadly expressed lineage genes) — this is a flag to
review, not proof of contamination. A rigorous estimate needs the raw
(unfiltered) droplet matrix or a decontamination model, which is the opt-in
correction step, not this detector.
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


def _clean_genes(markers) -> list[str]:
    out = []
    for m in (markers or []):
        s = str(m).strip()
        if s and s.lower() != "nan":
            out.append(s)
    return out


def assess_ambient_contamination(
    top_markers: dict[str, list] | None,
    *,
    ubiquity_fraction: float = 0.6,
    max_report: int = 15,
) -> dict[str, Any]:
    """Flag ambient-like contamination from cross-cluster top-marker ubiquity.

    Args:
        top_markers: {cluster_label: [top_marker_gene, ...]} from rank-genes.
        ubiquity_fraction: a gene listed as a top marker in at least this
            fraction of clusters is considered non-specific (ubiquitous).
        max_report: cap on the number of ubiquitous genes reported.

    Returns: {"status": "clean"|"warnings"|"unverified", "issues": [...],
              "metrics": {"n_clusters", "ubiquitous_genes", "ubiquity_fraction",
                          "ubiquitous_marker_share"}}.
    """
    top_markers = top_markers or {}
    per_cluster = {str(k): _clean_genes(v) for k, v in top_markers.items()}
    clusters_with_markers = {k: v for k, v in per_cluster.items() if v}
    n_clusters = len(clusters_with_markers)

    metrics: dict[str, Any] = {
        "n_clusters": n_clusters,
        "ubiquitous_genes": [],
        "ubiquity_fraction": float(ubiquity_fraction),
        "ubiquitous_marker_share": 0.0,
    }

    # Need at least two annotated clusters with markers to talk about
    # cross-cluster ubiquity at all.
    if n_clusters < 2:
        return {"status": "unverified", "issues": [], "metrics": metrics}

    # Count, for each gene, how many clusters list it as a top marker.
    gene_cluster_count: dict[str, int] = {}
    total_slots = 0
    for genes in clusters_with_markers.values():
        total_slots += len(genes)
        for g in set(genes):       # one vote per cluster
            gene_cluster_count[g] = gene_cluster_count.get(g, 0) + 1

    threshold = max(2, int(round(ubiquity_fraction * n_clusters)))
    ubiquitous = sorted(
        (g for g, c in gene_cluster_count.items() if c >= threshold),
        key=lambda g: (-gene_cluster_count[g], g),
    )

    # Share of top-marker slots occupied by ubiquitous genes (severity proxy).
    ubi_slots = sum(gene_cluster_count[g] for g in ubiquitous)
    share = round(ubi_slots / total_slots, 4) if total_slots else 0.0

    metrics["ubiquitous_genes"] = ubiquitous[:max_report]
    metrics["ubiquitous_marker_share"] = share

    issues: list[QCIssue] = []
    if ubiquitous:
        shown = ", ".join(ubiquitous[:max_report])
        more = "" if len(ubiquitous) <= max_report else \
               f" (+{len(ubiquitous) - max_report} more)"
        issues.append(QCIssue(
            "warning",
            "possible_ambient_contamination",
            (f"{len(ubiquitous)} gene(s) appear as top markers in >= "
             f"{threshold}/{n_clusters} clusters ({shown}{more}), filling "
             f"{share * 100:.0f}% of top-marker slots. Such ubiquitous "
             f"'markers' are a signature of ambient (soup) RNA leaking abundant "
             f"transcripts into every cell, though they can also reflect shared "
             f"biology — review before trusting cross-cluster overlap."),
            ("Consider the optional ambient-RNA decontamination step "
             "(SoupX/decontX) before pseudobulk/DE, or confirm these genes are "
             "genuinely co-expressed; do not over-interpret markers shared "
             "across unrelated cell types."),
        ))

    return {
        "status": "warnings" if issues else "clean",
        "issues": [i.as_dict() for i in issues],
        "metrics": metrics,
    }
