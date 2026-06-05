"""Deterministic cross-analysis pattern detection for the biological synthesis.

Pure set/sign math over the structured bulk DE + pathway results — NO LLM, no
biology hardcoded. Detects, with full traceability to counts:

- within-contrast convergence: a contrast produced both DE and pathway enrichment;
- cross-contrast convergence: two contrasts sharing a reference share DE genes,
  and (when up/down sets are available) how many move the SAME direction;
- cross-contrast divergence: contrast-specific genes and pathway terms;
- reliability: power / low-power / outlier flags that bound interpretation.

Every number a downstream claim states must come from one of these counts, so the
strict evidence verifier can confirm the discussion against the evidence card.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Any


@dataclass
class WithinContrastPattern:
    name: str
    n_de: int
    n_pathway_terms: int
    converges: bool                 # DE present AND pathway enrichment present

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name, "n_de": self.n_de,
            "n_pathway_terms": self.n_pathway_terms, "converges": self.converges,
        }


@dataclass
class CrossContrastPattern:
    contrast_a: str
    contrast_b: str
    shared_reference: str | None    # common denominator, if any (e.g. "WT")
    n_shared_genes: int
    n_specific_a: int
    n_specific_b: int
    n_direction_concordant: int     # shared genes moving the same way (or -1 if unknown)
    n_direction_discordant: int
    direction_known: bool
    shared_pathway_terms: list[str]
    specific_terms_a: list[str]
    specific_terms_b: list[str]
    n_terms_a: int
    n_terms_b: int

    @property
    def n_shared_terms(self) -> int:
        return len(self.shared_pathway_terms)

    def as_dict(self) -> dict[str, Any]:
        return {
            "contrast_a": self.contrast_a, "contrast_b": self.contrast_b,
            "shared_reference": self.shared_reference,
            "n_shared_genes": self.n_shared_genes,
            "n_specific_a": self.n_specific_a, "n_specific_b": self.n_specific_b,
            "n_direction_concordant": self.n_direction_concordant,
            "n_direction_discordant": self.n_direction_discordant,
            "direction_known": self.direction_known,
            "n_shared_terms": self.n_shared_terms,
            "shared_pathway_terms": self.shared_pathway_terms[:20],
            "specific_terms_a": self.specific_terms_a[:20],
            "specific_terms_b": self.specific_terms_b[:20],
            "n_terms_a": self.n_terms_a, "n_terms_b": self.n_terms_b,
        }


def _gene_ids(contrast: dict) -> set[str]:
    ids = contrast.get("all_sig_gene_ids") or contrast.get("all_sig_genes") or []
    return {str(g) for g in ids}


def _direction_sets(contrast: dict) -> tuple[set[str], set[str], bool]:
    up = contrast.get("up_gene_ids")
    down = contrast.get("down_gene_ids")
    if up is None and down is None:
        return set(), set(), False
    return ({str(g) for g in (up or [])},
            {str(g) for g in (down or [])}, True)


def _pathway_terms(contrast: dict) -> list[str]:
    """Ordered, de-duplicated enriched terms (ORA output order = significance)."""
    terms: list[str] = []
    seen: set[str] = set()
    for rows in (contrast.get("pathways") or {}).values():
        for row in rows or []:
            if isinstance(row, dict):
                t = row.get("term") or row.get("Term")
                if t and str(t) not in seen:
                    seen.add(str(t))
                    terms.append(str(t))
    return terms


def _successful(contrasts: list[dict]) -> list[dict]:
    return [c for c in (contrasts or [])
            if isinstance(c, dict) and c.get("status", "success") == "success"
            and int(c.get("n_significant", 0) or 0) > 0]


def detect_bulk_patterns(contrasts: list[dict]) -> dict[str, Any]:
    """Detect within- and cross-contrast patterns over bulk DE + pathway results.

    ``contrasts`` is ``agent_results['bulk_rna_agent']['findings']['contrasts']``.
    Returns a serializable manifest of patterns; empty sections when the data does
    not support them (e.g. no pathways => no convergence claim).
    """
    valid = _successful(contrasts)

    within: list[WithinContrastPattern] = []
    for c in valid:
        terms = _pathway_terms(c)
        within.append(WithinContrastPattern(
            name=str(c.get("name", "contrast")),
            n_de=int(c.get("n_significant", 0) or 0),
            n_pathway_terms=len(terms),
            converges=bool(terms),
        ))

    cross: list[CrossContrastPattern] = []
    for a, b in combinations(valid, 2):
        ga, gb = _gene_ids(a), _gene_ids(b)
        if not ga or not gb:
            continue
        shared = ga & gb
        up_a, down_a, known_a = _direction_sets(a)
        up_b, down_b, known_b = _direction_sets(b)
        direction_known = known_a and known_b
        concordant = discordant = -1
        if direction_known:
            concordant = len((up_a & up_b) | (down_a & down_b))
            discordant = len((up_a & down_b) | (down_a & up_b))
        ta, tb = _pathway_terms(a), _pathway_terms(b)
        tb_set, ta_set = set(tb), set(ta)
        den_a = str(a.get("denominator") or "")
        den_b = str(b.get("denominator") or "")
        cross.append(CrossContrastPattern(
            contrast_a=str(a.get("name", "A")),
            contrast_b=str(b.get("name", "B")),
            shared_reference=den_a if den_a and den_a == den_b else None,
            n_shared_genes=len(shared),
            n_specific_a=len(ga - gb),
            n_specific_b=len(gb - ga),
            n_direction_concordant=concordant,
            n_direction_discordant=discordant,
            direction_known=direction_known,
            # Ordered by each contrast's ORA significance (most enriched first).
            shared_pathway_terms=[t for t in ta if t in tb_set],
            specific_terms_a=[t for t in ta if t not in tb_set],
            specific_terms_b=[t for t in tb if t not in ta_set],
            n_terms_a=len(ta), n_terms_b=len(tb),
        ))

    return {
        "modalities_present": ["bulk RNA-seq"] if valid else [],
        "n_contrasts": len(valid),
        "within_contrast": [w.as_dict() for w in within],
        "cross_contrast": [x.as_dict() for x in cross],
        "reliability": _reliability(valid),
    }


def _reliability(contrasts: list[dict]) -> dict[str, Any]:
    powers = [c.get("power_estimate_at_lfc_min") for c in contrasts
              if isinstance(c.get("power_estimate_at_lfc_min"), (int, float))]
    low_power = [str(c.get("name")) for c in contrasts
                 if c.get("low_power_warning")]
    return {
        "min_power": round(min(powers), 3) if powers else None,
        "max_power": round(max(powers), 3) if powers else None,
        "low_power_contrasts": low_power,
    }
