"""V48 Slice I6 — External concordance & determinism scoring for candidate eRegulons.

Pure NumPy/SciPy scoring layer comparing ARIA's deterministic paired candidate
eRegulons (`aria/scripts/integration_eregulons.py`, `aria/utils/eregulon_core.py`)
against an external enhancer-GRN reference execution (`scenicplus==1.0.0`) as well
as quantifying SCENIC+'s own multi-seed stochastic GBM ceiling (`seed_a` vs `seed_b`).

Evaluates four complementary concordance tiers (`ADR-064`, `Concilio d9842f`):
  1. Shared testable TF & `(TF, sign_class)` module discovery (precision, recall,
     F1, Jaccard, exact analytical hypergeometric expected Jaccard, fold enrichment,
     and Fisher exact one-sided p-value).
  2. `(TF, gene, sign_class)` target-edge overlap (global universe + conditional
     per-shared-module target-gene Jaccard distribution).
  3. `(TF, peak, gene)` regulatory triplet overlap under both exact canonical
     peak coordinates (`chr:start-end`) and `>=50%` reciprocal genomic overlap.
  4. Per-cell regulon AUCell activity Spearman correlation across shared cells,
     calibrated against a seeded size-matched random-gene-set AUCell baseline.

Honesty contract (`ADR-002`): returns `status="not_run"` with an explicit reason
when external reference outputs are absent; never fabricates concordance values.
"""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import stats

from aria.utils.peak2gene_core import canonical_contig, parse_peak_coord

_SIGN_MAP: dict[str, str] = {
    "pos_tf_pos_p2g": "pos_tf_pos_p2g",
    "+/+": "pos_tf_pos_p2g",
    "+_+": "pos_tf_pos_p2g",
    "neg_tf_pos_p2g": "neg_tf_pos_p2g",
    "-/+": "neg_tf_pos_p2g",
    "-_+": "neg_tf_pos_p2g",
    "pos_tf_neg_p2g": "pos_tf_neg_p2g",
    "+/-": "pos_tf_neg_p2g",
    "+_-": "pos_tf_neg_p2g",
    "neg_tf_neg_p2g": "neg_tf_neg_p2g",
    "-/-": "neg_tf_neg_p2g",
    "-_-": "neg_tf_neg_p2g",
}


def canonical_sign_class(raw: str | None) -> str:
    """Normalize ARIA or SCENIC+ sign labels into ARIA's non-causal sign classes."""
    if not raw:
        return "pos_tf_pos_p2g"
    key = str(raw).strip().lower()
    if key in _SIGN_MAP:
        return _SIGN_MAP[key]
    if "neg" in key and "tf" in key and "neg" in key.split("tf", 1)[-1]:
        return "neg_tf_neg_p2g"
    if key.startswith("neg") or key.startswith("-"):
        return "neg_tf_pos_p2g"
    return "pos_tf_pos_p2g"


def canonical_peak_key(peak: str) -> str:
    """Canonicalize peak strings (`chr1-100-200` / `1:100-200`) to `chr1:100-200`."""
    parsed = parse_peak_coord(str(peak))
    if parsed is None:
        return str(peak).strip()
    chrom, s, e = parsed
    c = canonical_contig(chrom) or chrom
    return f"{c}:{int(s)}-{int(e)}"


def _rankdata_average(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64).ravel()
    n = arr.size
    if n == 0:
        return np.empty(0, dtype=np.float64)
    order = np.argsort(arr, kind="mergesort")
    sorted_vals = arr[order]
    ranks = np.empty(n, dtype=np.float64)
    i = 0
    while i < n:
        j = i + 1
        while j < n and sorted_vals[j] == sorted_vals[i]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + (j - 1)) + 1.0
        i = j
    return ranks


def spearman_corr(a: Sequence[float] | np.ndarray, b: Sequence[float] | np.ndarray) -> float | None:
    """Pure-NumPy Spearman rank correlation; returns None if constant or <2 finite pairs."""
    xa = np.asarray(a, dtype=np.float64).ravel()
    xb = np.asarray(b, dtype=np.float64).ravel()
    if xa.size != xb.size or xa.size < 2:
        return None
    mask = np.isfinite(xa) & np.isfinite(xb)
    if int(mask.sum()) < 2:
        return None
    ra = _rankdata_average(xa[mask])
    rb = _rankdata_average(xb[mask])
    sa = float(np.std(ra))
    sb = float(np.std(rb))
    if sa <= 1e-15 or sb <= 1e-15:
        return None
    val = float(np.corrcoef(ra, rb)[0, 1])
    return float(np.clip(val, -1.0, 1.0)) if np.isfinite(val) else None


def set_overlap_summary(
    set_a: set[Any],
    set_b: set[Any],
    *,
    universe_size: int | None = None,
) -> dict[str, Any]:
    """Compute overlap, precision, recall, F1, Jaccard, and analytical null expectation."""
    n_a = len(set_a)
    n_b = len(set_b)
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    prec = float(inter) / float(n_a) if n_a > 0 else (1.0 if n_b == 0 else 0.0)
    rec = float(inter) / float(n_b) if n_b > 0 else (1.0 if n_a == 0 else 0.0)
    f1 = (2.0 * prec * rec / (prec + rec)) if (prec + rec) > 0.0 else 0.0
    jacc = float(inter) / float(union) if union > 0 else (1.0 if (n_a == 0 and n_b == 0) else 0.0)

    out: dict[str, Any] = {
        "n_aria": int(n_a),
        "n_reference": int(n_b),
        "n_intersection": int(inter),
        "n_union": int(union),
        "precision_vs_ref": round(prec, 6),
        "recall_vs_ref": round(rec, 6),
        "f1": round(f1, 6),
        "jaccard": round(jacc, 6),
    }

    if universe_size is not None and int(universe_size) >= union and int(universe_size) > 0:
        u = int(universe_size)
        exp_inter = (float(n_a) * float(n_b)) / float(u)
        denom_exp = float(n_a + n_b) - exp_inter
        exp_jacc = (exp_inter / denom_exp) if denom_exp > 0.0 else 0.0
        fe_inter = (float(inter) / exp_inter) if exp_inter > 0.0 else None
        # Fisher / hypergeometric one-sided enrichment tail P(X >= inter)
        if inter > 0 and n_a > 0 and n_b > 0:
            pval = float(stats.hypergeom.sf(inter - 1, u, n_b, n_a))
        else:
            pval = 1.0
        out["universe_size"] = u
        out["expected_intersection_null"] = round(exp_inter, 6)
        out["expected_jaccard_null"] = round(exp_jacc, 6)
        out["fold_enrichment_over_null"] = (
            round(fe_inter, 4) if fe_inter is not None and math.isfinite(fe_inter) else None
        )
        out["hypergeom_pvalue"] = float(pval)
    return out


def _reciprocal_triplet_overlap_count(
    triplets_a: set[tuple[str, str, str]],
    triplets_b: set[tuple[str, str, str]],
    *,
    min_reciprocal_frac: float = 0.50,
) -> dict[str, Any]:
    """Count triplets in `A` matched by `>=min_reciprocal_frac` genomic peak overlap in `B`."""
    if not triplets_a or not triplets_b:
        return {
            "n_aria_matched_reciprocal_50pct": 0,
            "n_ref_matched_reciprocal_50pct": 0,
            "reciprocal_50pct_jaccard_approx": 0.0,
        }

    def _index(trips: set[tuple[str, str, str]]) -> dict[tuple[str, str, str], list[tuple[int, int, str]]]:
        idx: dict[tuple[str, str, str], list[tuple[int, int, str]]] = {}
        for tf, pk, gn in trips:
            parsed = parse_peak_coord(pk)
            if parsed is None:
                continue
            c = canonical_contig(parsed[0]) or parsed[0]
            idx.setdefault((tf, gn, c), []).append((int(parsed[1]), int(parsed[2]), pk))
        return idx

    idx_b = _index(triplets_b)
    matched_a = 0
    matched_b_keys: set[tuple[str, str, str]] = set()
    for tf, pk_a, gn in triplets_a:
        pa = parse_peak_coord(pk_a)
        if pa is None:
            if (tf, pk_a, gn) in triplets_b:
                matched_a += 1
                matched_b_keys.add((tf, pk_a, gn))
            continue
        ca = canonical_contig(pa[0]) or pa[0]
        sa, ea = int(pa[1]), int(pa[2])
        len_a = max(1, ea - sa)
        cands = idx_b.get((tf, gn, ca), [])
        hit_a = False
        for sb, eb, pk_b in cands:
            ov = min(ea, eb) - max(sa, sb)
            if ov <= 0:
                continue
            len_b = max(1, eb - sb)
            if (float(ov) / float(len_a) >= min_reciprocal_frac) and (
                float(ov) / float(len_b) >= min_reciprocal_frac
            ):
                hit_a = True
                matched_b_keys.add((tf, pk_b, gn))
        if hit_a:
            matched_a += 1

    matched_b = len(matched_b_keys)
    eff_inter = min(matched_a, matched_b)
    eff_union = len(triplets_a) + len(triplets_b) - eff_inter
    jacc = float(eff_inter) / float(eff_union) if eff_union > 0 else 0.0
    return {
        "n_aria_matched_reciprocal_50pct": int(matched_a),
        "n_ref_matched_reciprocal_50pct": int(matched_b),
        "reciprocal_50pct_jaccard_approx": round(jacc, 6),
    }


def _normalize_regulon_records(
    triplets: Sequence[Mapping[str, Any]],
    *,
    tf_universe: set[str] | None = None,
    gene_universe: set[str] | None = None,
) -> dict[str, Any]:
    """Extract canonical TFs, `(TF, sign)`, `(TF, gene, sign)` edges, and `(TF, peak, gene)` triplets."""
    tfs: set[str] = set()
    modules: dict[tuple[str, str], set[str]] = {}
    edges: set[tuple[str, str, str]] = set()
    canon_triplets: set[tuple[str, str, str]] = set()
    signed_triplets: set[tuple[str, str, str, str]] = set()

    for row in triplets:
        tf = str(row.get("tf") or row.get("TF") or "").strip()
        gene = str(row.get("gene") or row.get("Gene") or "").strip()
        peak = str(row.get("peak") or row.get("Region") or "").strip()
        if not tf or not gene:
            continue
        # Enforce identical self-loop exclusion on both ARIA and SCENIC+ (`g != TF`)
        if tf.upper() == gene.upper():
            continue
        if tf_universe is not None and tf not in tf_universe:
            continue
        if gene_universe is not None and gene not in gene_universe:
            continue
        sc = canonical_sign_class(
            row.get("sign_class") or row.get("sign") or row.get("direction")
        )
        tfs.add(tf)
        modules.setdefault((tf, sc), set()).add(gene)
        edges.add((tf, gene, sc))
        if peak:
            cpk = canonical_peak_key(peak)
            canon_triplets.add((tf, cpk, gene))
            signed_triplets.add((tf, cpk, gene, sc))

    return {
        "tfs": tfs,
        "modules": modules,
        "edges": edges,
        "triplets": canon_triplets,
        "signed_triplets": signed_triplets,
    }


def score_eregulon_concordance(
    aria_triplets: Sequence[Mapping[str, Any]] | None,
    ref_triplets: Sequence[Mapping[str, Any]] | None,
    *,
    shared_tf_universe: Sequence[str] | None = None,
    shared_gene_universe: Sequence[str] | None = None,
    aria_cell_aucell: Mapping[str, Sequence[float] | np.ndarray] | None = None,
    ref_cell_aucell: Mapping[str, Sequence[float] | np.ndarray] | None = None,
    random_aucell_baseline_spearmans: Sequence[float] | None = None,
    restrict_to_sign_classes: Sequence[str] | None = ("pos_tf_pos_p2g", "neg_tf_pos_p2g"),
) -> dict[str, Any]:
    """Score external concordance between ARIA candidate eRegulons and SCENIC+."""
    if aria_triplets is None or ref_triplets is None:
        return {
            "status": "not_run",
            "reason": "External reference triplets or ARIA triplets were not provided.",
        }

    tf_univ = set(str(x) for x in shared_tf_universe) if shared_tf_universe is not None else None
    gn_univ = set(str(x) for x in shared_gene_universe) if shared_gene_universe is not None else None
    allowed_signs = (
        {canonical_sign_class(s) for s in restrict_to_sign_classes}
        if restrict_to_sign_classes is not None
        else None
    )

    def _filter_signs(rows: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
        if allowed_signs is None:
            return list(rows)
        return [
            r
            for r in rows
            if canonical_sign_class(
                r.get("sign_class") or r.get("sign") or r.get("direction")
            )
            in allowed_signs
        ]

    a_norm = _normalize_regulon_records(
        _filter_signs(aria_triplets), tf_universe=tf_univ, gene_universe=gn_univ
    )
    r_norm = _normalize_regulon_records(
        _filter_signs(ref_triplets), tf_universe=tf_univ, gene_universe=gn_univ
    )

    n_tf_u = len(tf_univ) if tf_univ is not None else None
    n_gn_u = len(gn_univ) if gn_univ is not None else None
    n_sign_u = len(allowed_signs) if allowed_signs is not None else 4
    mod_u = (n_tf_u * n_sign_u) if n_tf_u is not None else None
    edge_u = (
        (n_tf_u * max(1, n_gn_u - 1) * n_sign_u)
        if (n_tf_u is not None and n_gn_u is not None)
        else None
    )

    tf_summary = set_overlap_summary(a_norm["tfs"], r_norm["tfs"], universe_size=n_tf_u)
    mod_summary = set_overlap_summary(
        set(a_norm["modules"].keys()),
        set(r_norm["modules"].keys()),
        universe_size=mod_u,
    )
    edge_summary = set_overlap_summary(
        a_norm["edges"], r_norm["edges"], universe_size=edge_u
    )
    triplet_exact = set_overlap_summary(a_norm["triplets"], r_norm["triplets"])
    triplet_recip = _reciprocal_triplet_overlap_count(
        a_norm["triplets"], r_norm["triplets"], min_reciprocal_frac=0.50
    )

    # Conditional target-gene overlap within shared `(TF, sign_class)` modules
    shared_mod_keys = sorted(set(a_norm["modules"].keys()) & set(r_norm["modules"].keys()))
    per_mod_jacc: list[float] = []
    per_mod_rows: list[dict[str, Any]] = []
    cond_edges_a: set[tuple[str, str, str]] = set()
    cond_edges_r: set[tuple[str, str, str]] = set()
    for tf, sc in shared_mod_keys:
        ga = a_norm["modules"][(tf, sc)]
        gr = r_norm["modules"][(tf, sc)]
        ov = set_overlap_summary(
            ga, gr, universe_size=(max(1, n_gn_u - 1) if n_gn_u is not None else None)
        )
        per_mod_jacc.append(float(ov["jaccard"]))
        per_mod_rows.append(
            {
                "tf": tf,
                "sign_class": sc,
                "n_aria_genes": len(ga),
                "n_ref_genes": len(gr),
                "n_shared_genes": ov["n_intersection"],
                "jaccard": ov["jaccard"],
                "fold_enrichment_over_null": ov.get("fold_enrichment_over_null"),
            }
        )
        for g in ga:
            cond_edges_a.add((tf, g, sc))
        for g in gr:
            cond_edges_r.add((tf, g, sc))

    cond_edge_u = (
        len(shared_mod_keys) * max(1, n_gn_u - 1)
        if (n_gn_u is not None and shared_mod_keys)
        else None
    )
    cond_edge_summary = set_overlap_summary(
        cond_edges_a, cond_edges_r, universe_size=cond_edge_u
    )
    per_mod_rows.sort(key=lambda d: (-d["jaccard"], -d["n_shared_genes"], d["tf"]))

    # Per-sign breakdown
    by_sign: dict[str, Any] = {}
    sign_iter = sorted(allowed_signs) if allowed_signs else sorted(set(_SIGN_MAP.values()))
    for sc in sign_iter:
        a_sc_tfs = {tf for (tf, s) in a_norm["modules"] if s == sc}
        r_sc_tfs = {tf for (tf, s) in r_norm["modules"] if s == sc}
        by_sign[sc] = set_overlap_summary(a_sc_tfs, r_sc_tfs, universe_size=n_tf_u)

    # Per-cell AUCell Spearman correlation on matched `(TF, sign_class)` keys
    aucell_summary: dict[str, Any] = {"status": "not_computed"}
    if aria_cell_aucell is not None and ref_cell_aucell is not None:
        shared_act_keys = sorted(set(aria_cell_aucell.keys()) & set(ref_cell_aucell.keys()))
        corrs: list[float] = []
        per_reg_corrs: list[dict[str, Any]] = []
        for k in shared_act_keys:
            rho = spearman_corr(aria_cell_aucell[k], ref_cell_aucell[k])
            if rho is not None:
                corrs.append(rho)
                per_reg_corrs.append({"module_key": k, "spearman_rho": round(rho, 6)})
        per_reg_corrs.sort(key=lambda d: -d["spearman_rho"])
        if corrs:
            arr_c = np.asarray(corrs, dtype=np.float64)
            q25, q75 = np.percentile(arr_c, [25.0, 75.0])
            aucell_summary = {
                "status": "computed",
                "n_shared_modules_scored": len(corrs),
                "median_spearman": round(float(np.median(arr_c)), 6),
                "mean_spearman": round(float(np.mean(arr_c)), 6),
                "iqr_spearman": round(float(q75 - q25), 6),
                "q25_spearman": round(float(q25), 6),
                "q75_spearman": round(float(q75), 6),
                "frac_spearman_gt_0_5": round(float(np.mean(arr_c >= 0.5)), 6),
                "frac_spearman_gt_0_7": round(float(np.mean(arr_c >= 0.7)), 6),
                "top_matched_modules": per_reg_corrs[:15],
            }
            if random_aucell_baseline_spearmans:
                b_arr = np.asarray(
                    [x for x in random_aucell_baseline_spearmans if x is not None and np.isfinite(x)],
                    dtype=np.float64,
                )
                if b_arr.size > 0:
                    b_med = float(np.median(b_arr))
                    aucell_summary["random_geneset_baseline_median_spearman"] = round(b_med, 6)
                    aucell_summary["median_spearman_excess_over_random"] = round(
                        float(np.median(arr_c)) - b_med, 6
                    )
        else:
            aucell_summary = {
                "status": "empty_overlap",
                "n_shared_modules_scored": 0,
            }

    j_arr = np.asarray(per_mod_jacc, dtype=np.float64)
    return {
        "status": "computed",
        "tf_concordance": tf_summary,
        "module_sign_concordance": mod_summary,
        "by_sign_class": by_sign,
        "edge_concordance_global": edge_summary,
        "edge_concordance_shared_modules": cond_edge_summary,
        "per_shared_module_gene_jaccard": {
            "n_shared_modules": len(per_mod_jacc),
            "median_jaccard": round(float(np.median(j_arr)), 6) if j_arr.size else None,
            "mean_jaccard": round(float(np.mean(j_arr)), 6) if j_arr.size else None,
            "max_jaccard": round(float(np.max(j_arr)), 6) if j_arr.size else None,
            "top_shared_modules": per_mod_rows[:20],
        },
        "triplet_concordance": {
            **triplet_exact,
            **triplet_recip,
        },
        "cell_activity_concordance": aucell_summary,
    }


def score_eregulon_determinism(
    eregulons_a: Sequence[Mapping[str, Any]],
    triplets_a: Sequence[Mapping[str, Any]],
    eregulons_b: Sequence[Mapping[str, Any]],
    triplets_b: Sequence[Mapping[str, Any]],
    *,
    rna_aucell_a: np.ndarray | None = None,
    rna_aucell_b: np.ndarray | None = None,
) -> dict[str, Any]:
    """Compare two ARIA runs (`seed=0` vs `seed=0` bit-exact, or `seed=0` vs `seed=1`)."""
    mods_a = {
        (str(r["tf"]), canonical_sign_class(r.get("sign_class"))): float(r.get("regulon_padj", 1.0))
        for r in eregulons_a
    }
    mods_b = {
        (str(r["tf"]), canonical_sign_class(r.get("sign_class"))): float(r.get("regulon_padj", 1.0))
        for r in eregulons_b
    }
    trips_a = {
        (
            str(t["tf"]),
            canonical_peak_key(str(t["peak"])),
            str(t["gene"]),
            canonical_sign_class(t.get("sign_class")),
        )
        for t in triplets_a
    }
    trips_b = {
        (
            str(t["tf"]),
            canonical_peak_key(str(t["peak"])),
            str(t["gene"]),
            canonical_sign_class(t.get("sign_class")),
        )
        for t in triplets_b
    }

    mod_ov = set_overlap_summary(set(mods_a.keys()), set(mods_b.keys()))
    trip_ov = set_overlap_summary(trips_a, trips_b)

    shared_mods = sorted(set(mods_a.keys()) & set(mods_b.keys()))
    max_padj_diff = (
        max(abs(mods_a[k] - mods_b[k]) for k in shared_mods) if shared_mods else 0.0
    )

    aucell_max_diff: float | None = None
    if (
        rna_aucell_a is not None
        and rna_aucell_b is not None
        and rna_aucell_a.shape == rna_aucell_b.shape
    ):
        aucell_max_diff = (
            float(np.max(np.abs(rna_aucell_a - rna_aucell_b)))
            if rna_aucell_a.size > 0
            else 0.0
        )

    bit_exact = (
        set(mods_a.keys()) == set(mods_b.keys())
        and trips_a == trips_b
        and max_padj_diff == 0.0
        and (aucell_max_diff is None or aucell_max_diff == 0.0)
    )
    return {
        "bit_exact": bool(bit_exact),
        "module_jaccard": mod_ov["jaccard"],
        "triplet_jaccard": trip_ov["jaccard"],
        "n_modules_a": len(mods_a),
        "n_modules_b": len(mods_b),
        "n_triplets_a": len(trips_a),
        "n_triplets_b": len(trips_b),
        "max_shared_padj_abs_diff": float(max_padj_diff),
        "max_rna_aucell_abs_diff": aucell_max_diff,
    }
