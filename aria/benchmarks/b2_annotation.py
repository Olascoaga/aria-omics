"""B2 multi-annotator kit: labeling-sheet I/O + inter-annotator agreement +
adjudication for the claim/narrative governance gold standard.

The preliminary B2 (``governance_b2``) hand-labels claims inline. The preprint B2
needs an INDEPENDENT human gold: ≥2 annotators label the same claims blind against
a published rubric, agreement is quantified (Cohen's κ for two raters, Krippendorff's
α for ≥2 with missing data), and disagreements are adjudicated into final labels.

This module is pure (csv + math only); it neither runs ARIA governance nor decides
labels. The failure taxonomy is the single source ``B2_LABELS``; sheets and loaders
validate against it so a typo'd label is caught, not silently scored.
"""
from __future__ import annotations

import csv
import io
from typing import Iterable, Mapping, Sequence

# Single-source failure taxonomy (must match docs/.../B2_RUBRIC.md and governance_b2).
B2_LABELS = (
    "clean",
    "unsupported",
    "overclaim",
    "fabricated",
    "causal_inflation",
    "missing_caveat",
)

_SHEET_FIELDS = ("claim_id", "modality", "analysis", "claim_text",
                 "evidence_summary", "label", "notes")


def export_labeling_sheet(corpus_rows: Sequence[Mapping[str, object]]) -> str:
    """Serialize claims to a blind labeling CSV (one row per claim, empty label).

    ``corpus_rows`` items provide claim_id/modality/analysis/claim_text/
    evidence_summary. ``label`` and ``notes`` are left blank for the annotator. The
    gold/reference label is intentionally NOT exported (blind annotation).
    """
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=_SHEET_FIELDS)
    writer.writeheader()
    for r in corpus_rows:
        writer.writerow({
            "claim_id": r["claim_id"],
            "modality": r.get("modality", ""),
            "analysis": r.get("analysis", ""),
            "claim_text": r.get("claim_text", ""),
            "evidence_summary": r.get("evidence_summary", ""),
            "label": "",
            "notes": "",
        })
    return buf.getvalue()


def load_annotations(csv_text: str) -> dict[str, str]:
    """Load a filled labeling sheet -> {claim_id: label}. Validates the taxonomy.

    Blank labels are skipped (not yet labelled). An unknown label raises, so a
    mislabelled sheet is caught before it pollutes agreement/scoring.
    """
    out: dict[str, str] = {}
    reader = csv.DictReader(io.StringIO(csv_text))
    for row in reader:
        cid = (row.get("claim_id") or "").strip()
        label = (row.get("label") or "").strip()
        if not cid or not label:
            continue
        if label not in B2_LABELS:
            raise ValueError(
                f"claim {cid!r}: unknown label {label!r}; "
                f"allowed: {', '.join(B2_LABELS)}")
        out[cid] = label
    return out


def _shared_ids(a: Mapping[str, str], b: Mapping[str, str]) -> list[str]:
    return sorted(set(a) & set(b))


def cohen_kappa(a: Mapping[str, str], b: Mapping[str, str]) -> dict[str, object]:
    """Cohen's κ between two annotators over their commonly-labelled claims."""
    ids = _shared_ids(a, b)
    n = len(ids)
    if n == 0:
        return {"kappa": None, "n": 0, "po": None, "pe": None,
                "reason": "no commonly labelled claims"}
    la = [a[i] for i in ids]
    lb = [b[i] for i in ids]
    po = sum(1 for x, y in zip(la, lb) if x == y) / n
    labels = set(la) | set(lb)
    pe = 0.0
    for lab in labels:
        pa = sum(1 for x in la if x == lab) / n
        pb = sum(1 for y in lb if y == lab) / n
        pe += pa * pb
    kappa = 1.0 if pe == 1.0 else (po - pe) / (1.0 - pe)
    return {"kappa": round(kappa, 4), "n": n, "po": round(po, 4), "pe": round(pe, 4)}


def krippendorff_alpha(
    annotations: Sequence[Mapping[str, str]],
) -> dict[str, object]:
    """Krippendorff's α (nominal) over ≥2 annotators, tolerant of missing labels.

    α = 1 - Do/De, with nominal disagreement (mismatched label pairs) aggregated
    over all units that have ≥2 ratings. α=1 perfect, 0 chance, <0 worse-than-chance.
    """
    units: dict[str, list[str]] = {}
    for ann in annotations:
        for cid, label in ann.items():
            units.setdefault(cid, []).append(label)
    # keep only units with >=2 ratings
    units = {u: v for u, v in units.items() if len(v) >= 2}
    if not units:
        return {"alpha": None, "n_units": 0, "reason": "no unit has >=2 ratings"}

    # value-by-value coincidence over pairable ratings
    label_counts: dict[str, int] = {}
    total_pairable = 0
    observed_disagreement = 0.0
    for ratings in units.values():
        m = len(ratings)
        total_pairable += m
        for lab in ratings:
            label_counts[lab] = label_counts.get(lab, 0) + 1
        # within-unit disagreement: fraction of mismatched ordered pairs / (m-1)
        mismatch = sum(1 for i in range(m) for j in range(m)
                       if i != j and ratings[i] != ratings[j])
        observed_disagreement += mismatch / (m - 1)
    do = observed_disagreement / total_pairable
    n = total_pairable
    # expected nominal disagreement from the marginal label distribution
    de = 1.0 - sum((c / n) ** 2 for c in label_counts.values())
    de = de * n / (n - 1) if n > 1 else de
    alpha = 1.0 if de == 0 else 1.0 - do / de
    return {"alpha": round(alpha, 4), "n_units": len(units),
            "do": round(do, 4), "de": round(de, 4)}


def adjudicate(
    annotations: Sequence[Mapping[str, str]],
    overrides: Mapping[str, str] | None = None,
) -> dict[str, dict[str, object]]:
    """Resolve final gold labels per claim.

    Unanimous (or majority) annotator labels stand; ties or conflicts are flagged
    ``needs_adjudication`` unless an ``overrides`` (adjudicator) label resolves them.
    Returns ``{claim_id: {gold, agreement, votes, adjudicated}}``.
    """
    overrides = overrides or {}
    votes: dict[str, dict[str, int]] = {}
    for ann in annotations:
        for cid, label in ann.items():
            votes.setdefault(cid, {})[label] = votes.setdefault(cid, {}).get(label, 0) + 1
    out: dict[str, dict[str, object]] = {}
    for cid, tally in votes.items():
        top = max(tally.values())
        winners = sorted(k for k, v in tally.items() if v == top)
        unanimous = len(tally) == 1
        if cid in overrides:
            gold, agreement, adj = overrides[cid], ("adjudicated"), True
            if gold not in B2_LABELS:
                raise ValueError(f"adjudicator label {gold!r} not in taxonomy")
        elif len(winners) == 1 and (unanimous or top * 2 > sum(tally.values())):
            gold, agreement, adj = winners[0], (
                "unanimous" if unanimous else "majority"), False
        else:
            gold, agreement, adj = None, "needs_adjudication", False
        out[cid] = {"gold": gold, "agreement": agreement,
                    "votes": dict(tally), "adjudicated": adj}
    return out
