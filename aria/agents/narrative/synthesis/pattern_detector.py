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
    top_shared_genes: list[str]     # shared genes, ranked by significance (symbols)
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
            "top_shared_genes": self.top_shared_genes[:20],
            "n_terms_a": self.n_terms_a, "n_terms_b": self.n_terms_b,
        }


def _gene_id_list(contrast: dict) -> list[str]:
    """Significant gene IDs, ranked by significance (the script sorts by padj)."""
    ids = contrast.get("all_sig_gene_ids") or contrast.get("all_sig_genes") or []
    return [str(g) for g in ids]


def _gene_ids(contrast: dict) -> set[str]:
    return set(_gene_id_list(contrast))


def _id_to_symbol(contrast: dict) -> dict[str, str]:
    ids = contrast.get("all_sig_gene_ids") or []
    syms = contrast.get("all_sig_genes") or []
    return {str(i): str(s) for i, s in zip(ids, syms)}


def _looks_like_raw_id(symbol: str) -> bool:
    s = str(symbol or "")
    return s.startswith(("ENSG", "ENSMUS", "ENST", "ENS")) or s.isdigit() or not s


def _top_shared_symbols(a: dict, b_ids: set[str], k: int = 10) -> list[str]:
    """Shared genes (ranked by contrast a's significance), named by symbol.

    Skips genes that only have a raw Ensembl/numeric id so the named list uses
    community-friendly symbols a reviewer expects.
    """
    sym_map = _id_to_symbol(a)
    out: list[str] = []
    for gid in _gene_id_list(a):
        if gid not in b_ids:
            continue
        sym = sym_map.get(gid, gid)
        if not _looks_like_raw_id(sym) and sym not in out:
            out.append(sym)
        if len(out) >= k:
            break
    return out


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
            top_shared_genes=_top_shared_symbols(a, gb),
            n_terms_a=len(ta), n_terms_b=len(tb),
        ))

    return {
        "modalities_present": ["bulk RNA-seq"] if valid else [],
        "n_contrasts": len(valid),
        "within_contrast": [w.as_dict() for w in within],
        "cross_contrast": [x.as_dict() for x in cross],
        "reliability": _reliability(valid),
    }


def _scrna_findings(agent_result: dict) -> dict[str, Any]:
    """Return the structured scRNA findings envelope, tolerating legacy shapes."""
    if not isinstance(agent_result, dict):
        return {}
    findings = agent_result.get("findings", agent_result)
    if isinstance(findings, dict) and "scRNA" in findings:
        nested = findings.get("scRNA") or {}
        if isinstance(nested, dict):
            return nested.get("findings", nested) or {}
    return findings if isinstance(findings, dict) else {}


def _scrna_pathway_terms(block: dict) -> list[str]:
    terms: list[str] = []
    seen: set[str] = set()
    for rows in (block.get("results", {}) or {}).values():
        for row in rows or []:
            if isinstance(row, dict):
                term = row.get("term") or row.get("Term")
                if term and str(term) not in seen:
                    seen.add(str(term))
                    terms.append(str(term))
    return terms


def _scrna_successful_pseudobulk(findings: dict) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    pb = (findings or {}).get("pseudobulk_de") or {}
    for group, info in (pb.get("per_group", {}) or {}).items():
        for comparison, comp in (info.get("per_comparison", {}) or {}).items():
            if comp.get("status", "success") != "success":
                continue
            n_sig = int(comp.get("n_significant",
                                 comp.get("n_significant_global", 0)) or 0)
            if n_sig <= 0:
                continue
            out.append({
                "group": str(group),
                "comparison": str(comparison),
                "block_key": f"{group}::{comparison}",
                "n_de": n_sig,
                "n_up": int(comp.get("n_up", comp.get("n_up_global", 0)) or 0),
                "n_down": int(comp.get("n_down",
                                        comp.get("n_down_global", 0)) or 0),
                "composition_corrected": bool(
                    comp.get("corrected_for_composition")
                ),
                "power": comp.get("power_estimate_at_lfc_min"),
                "low_power": bool(comp.get("low_power_warning")),
                "top_genes": [
                    str(g.get("gene"))
                    for g in (comp.get("top_genes") or [])
                    if isinstance(g, dict) and g.get("gene")
                ][:8],
            })
    return sorted(out, key=lambda x: x["n_de"], reverse=True)


def detect_scrna_patterns(agent_result: dict) -> dict[str, Any]:
    """Detect scRNA synthesis patterns from measured structured outputs.

    The detector does not infer cell identities, mechanisms, or gene programs from
    names. It only summarizes support already produced by validated scRNA steps:
    donor-level pseudobulk DE, per-block ORA, abundance shifts, LIANA, and
    trajectory context.
    """
    findings = _scrna_findings(agent_result)
    if not findings:
        return {"modalities_present": [], "n_cells": None}

    qc = findings.get("qc") or {}
    ann = findings.get("cell_types") or {}
    has_labels = bool((ann.get("cell_types") or {}))
    pwp = findings.get("pseudobulk_pathways") or {}
    pathway_blocks = pwp.get("per_cluster", {}) or {}

    pb_blocks = _scrna_successful_pseudobulk(findings)
    for block in pb_blocks:
        pw = pathway_blocks.get(block["block_key"]) or {}
        terms = _scrna_pathway_terms(pw)
        block["n_pathway_terms"] = int(pw.get("n_significant", len(terms)) or 0)
        block["top_pathway_terms"] = terms[:5]

    da = findings.get("differential_abundance") or {}
    shifted: list[str] = []
    n_da_tests = 0
    for comp in (da.get("per_comparison", {}) or {}).values():
        rows = comp.get("per_cell_type", []) or []
        n_da_tests += len(rows)
        for row in rows:
            if row.get("significant") and row.get("name"):
                shifted.append(str(row["name"]))

    ccc = findings.get("cell_communication") or {}
    traj = findings.get("trajectory") or {}
    paga = traj.get("paga", {}) or {}
    pt = traj.get("pseudotime", {}) or {}

    powers = [b["power"] for b in pb_blocks
              if isinstance(b.get("power"), (int, float))]
    return {
        "modalities_present": ["scRNA-seq"],
        "n_cells": qc.get("n_cells_after"),
        "resolution": "cell-type" if has_labels else "cluster",
        "n_pseudobulk_blocks": len(pb_blocks),
        "strongest_pseudobulk": pb_blocks[0] if pb_blocks else None,
        "n_pathway_supported_blocks": sum(
            1 for b in pb_blocks if b.get("n_pathway_terms", 0) > 0
        ),
        "n_abundance_tests": n_da_tests,
        "n_abundance_shifts": len(shifted),
        "abundance_shift_labels": shifted[:8],
        "cellcomm": {
            "ran": ccc.get("status") in {"done", "success"},
            "n_interactions": int(ccc.get("n_interactions", 0) or 0),
            "n_cell_types": int(ccc.get("n_cell_types", 0) or 0),
        },
        "trajectory": {
            "ran": traj.get("status") in {"done", "success"},
            "n_connections": int(paga.get("n_connections", 0) or 0),
            "n_strong": int(paga.get("n_strong", 0) or 0),
            "dpt_computed": bool(pt.get("computed")),
        },
        "reliability": {
            "min_power": round(min(powers), 3) if powers else None,
            "max_power": round(max(powers), 3) if powers else None,
            "low_power_blocks": [
                f"{b['group']}::{b['comparison']}" for b in pb_blocks
                if b.get("low_power")
            ],
            "lognorm_recovered": bool(
                (findings.get("pseudobulk_de") or {}).get("lognorm_recovered")
            ),
        },
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
