"""Cross-modal footprint <-> RNA concordance (scATAC P4.3 governance tie-in).

Pure synthesis over structured results (ADR-028 Slice 2, cross-modal): for each top
differential-binding TF from TOBIAS BINDetect (footprint occupancy higher in cell-type
group A or B), check whether the TF's OWN gene is concordantly higher in that group's
RNA. A concordant pair (footprint up in A AND the TF gene's RNA up in A) is ASSOCIATIVE
cross-modal support, never causal regulation. No LLM, no hardcoded biology — the TF set
comes from the run, the RNA means come from the paired data.
"""
from __future__ import annotations

import math
from typing import Any


def _finite_mean(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        num = float(value)
    except (TypeError, ValueError):
        return None
    return num if math.isfinite(num) else None


def _tf_genes(tf_name: str) -> list[str]:
    """Resolve a JASPAR/TOBIAS motif name to its candidate gene symbol(s): split a
    dimer (``FOSL2::JUND``) into components and uppercase. The TF gene RNA is checked
    for any component (a dimer is concordant if either subunit's RNA agrees)."""
    raw = str(tf_name).split("(")[0].strip()
    parts = raw.split("::") if "::" in raw else [raw]
    return [p.strip().upper() for p in parts if p.strip()]


def _rna_favored_group(genes: list[str], rna_group_means: dict[str, dict[str, float]],
                       group_a: str, group_b: str,
                       target_group: str | None = None) -> tuple[str | None, dict | None]:
    """Return the group with the higher RNA mean across the TF's gene(s) (and the
    component means used), or ``(None, None)`` when no component gene has finite
    non-tied RNA means (ties return ``None`` / neutral instead of biasing group_a)."""
    ci_means: dict[str, dict[str, Any]] = {}
    for k, v in (rna_group_means or {}).items():
        if isinstance(v, dict):
            ku = str(k).strip().upper()
            if ku and (ku not in ci_means or str(k) == ku):
                ci_means[ku] = v

    non_tied: list[tuple[str, float, float, str]] = []
    tied: list[tuple[str, float, float]] = []
    for g in genes:
        gu = str(g).strip().upper()
        means = ci_means.get(gu)
        if not means or group_a not in means or group_b not in means:
            continue
        ma = _finite_mean(means[group_a])
        mb = _finite_mean(means[group_b])
        if ma is None or mb is None:
            continue
        if ma == mb:
            tied.append((gu, ma, mb))
            continue
        grp = group_a if ma > mb else group_b
        non_tied.append((gu, ma, mb, grp))

    if not non_tied:
        if tied:
            g0, ma0, mb0 = tied[0]
            return None, {"gene": g0, group_a: round(ma0, 4),
                          group_b: round(mb0, 4)}
        return None, None

    pool = non_tied
    if target_group is not None:
        supporting = [c for c in non_tied if c[3] == target_group]
        if supporting:
            pool = supporting
    best = max(pool, key=lambda c: abs(c[1] - c[2]))
    return best[3], {"gene": best[0], group_a: round(best[1], 4),
                     group_b: round(best[2], 4)}


def footprint_rna_concordance(differential_summary: dict,
                              rna_group_means: dict[str, dict[str, float]],
                              group_a: str, group_b: str,
                              top_n: int = 10) -> dict[str, Any]:
    """Annotate the top differential-binding TFs (toward each group) with whether the
    TF gene's RNA is concordantly higher in that same group. Returns per-direction
    annotated lists plus concordance counts. A TF whose gene RNA is unavailable is
    marked ``rna_checked=False`` and excluded from the rate (honest, not assumed)."""
    out: dict[str, Any] = {"group_a": group_a, "group_b": group_b}
    n_checked = n_concordant = 0
    for favored in (group_a, group_b):
        annotated = []
        for tf in (differential_summary.get(f"top_toward_{favored}") or [])[:top_n]:
            genes = _tf_genes(tf.get("tf", ""))
            rna_group, means = _rna_favored_group(genes, rna_group_means,
                                                  group_a, group_b,
                                                  target_group=favored)
            checked = rna_group is not None
            concordant = checked and rna_group == favored
            if checked:
                n_checked += 1
                n_concordant += int(concordant)
            annotated.append({**tf, "rna_checked": checked,
                              "rna_favored_group": rna_group,
                              "rna_concordant": concordant if checked else None,
                              "rna_means": means})
        out[f"toward_{favored}"] = annotated
    out["n_checked"] = n_checked
    out["n_concordant"] = n_concordant
    out["concordance_rate"] = (round(n_concordant / n_checked, 4)
                               if n_checked else None)
    return out
