"""Benchmark A1 reference-data lane: SEQC/MAQC TaqMan validation.

The synthetic A1 lane (``synthetic_de.py``) proves ARIA's bulk DE recovers a
*simulated* truth. This module validates the same DE path against an *external*
reference: the MAQC/SEQC samples (A = UHRR, B = HBRR, optional titration
mixtures C/D) with TaqMan qPCR log2(A/B) as the validated ground truth for ~1000
genes (MAQC-I, Nature Biotechnology 2006; SEQC, Nature Biotechnology 2014).

It is data-gated and never fabricates: the scorer runs only on a real staged
bundle, and ``run_seqc_maqc_a1_benchmark`` skips honestly (``status="skipped"``)
when the bundle is absent. ARIA computes nothing here it cannot measure — TaqMan
truth and the count matrix come entirely from the staged reference bundle.

Bundle schema (``ARIA_SEQC_MAQC_BUNDLE`` / ``--bundle``):
- ``counts.tsv``  : column ``gene`` + one column per sample (raw integer counts).
- ``samples.tsv`` : columns ``sample``, ``group`` (group in {A, B, C, D}).
- ``taqman.tsv``  : columns ``gene``, ``log2_ab`` (TaqMan log2(A/B) truth).
- ``manifest.json``: provenance (source, sha256, dates) written by the fetcher.

Standard SEQC/MAQC DE-evaluation axes reported here:
1. LFC concordance: Spearman/Pearson of ARIA's RNA-seq log2(A/B) vs TaqMan.
2. TaqMan-DE detection: AUC of ARIA's significance ranking, with TaqMan-DE
   genes (|log2| >= threshold) as positives and TaqMan-null genes as negatives,
   plus recall/empirical-FDR at the padj gate.
3. Titration monotonicity: fraction of TaqMan-DE genes whose group means are
   monotone across A -> C -> D -> B (requires the C/D mixtures; else not
   computed).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

__all__ = [
    "load_seqc_maqc_bundle",
    "score_seqc_maqc_a1",
    "run_seqc_maqc_a1_benchmark",
    "write_seqc_maqc_a1_figure",
    "run_seqc_maqc_multisite",
    "write_seqc_multisite_figure",
    "score_ercc_dose_response",
    "run_ercc_dose_response",
    "write_ercc_figure",
]


def _finite_float(value: Any) -> float:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return 0.0
    if f != f or f in (float("inf"), float("-inf")):
        return 0.0
    return f


def _auc(pos_scores: list[float], neg_scores: list[float]) -> float | None:
    """Tie-aware ROC AUC = P(score(positive) > score(negative)) via the
    Mann-Whitney U statistic. Returns None when either class is empty."""
    import numpy as np

    pos = np.asarray(pos_scores, dtype=float)
    neg = np.asarray(neg_scores, dtype=float)
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    if pos.size == 0 or neg.size == 0:
        return None
    allv = np.concatenate([pos, neg])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(allv.size, dtype=float)
    ranks[order] = np.arange(1, allv.size + 1, dtype=float)
    # Average ranks for ties so the AUC is unbiased under score ties.
    _, inv, counts = np.unique(allv, return_inverse=True, return_counts=True)
    cum = np.cumsum(counts)
    start = cum - counts
    avg = (start + cum + 1) / 2.0
    ranks = avg[inv]
    rank_sum_pos = ranks[: pos.size].sum()
    u = rank_sum_pos - pos.size * (pos.size + 1) / 2.0
    return float(u / (pos.size * neg.size))


def load_seqc_maqc_bundle(bundle_dir: str | Path) -> dict[str, Any] | None:
    """Load a staged SEQC/MAQC reference bundle. Returns None if absent so the
    runner can skip honestly; raises only when a present bundle is malformed."""
    import json
    import pandas as pd

    bundle = Path(bundle_dir)
    counts_p = bundle / "counts.tsv"
    samples_p = bundle / "samples.tsv"
    taqman_p = bundle / "taqman.tsv"
    if not (counts_p.exists() and samples_p.exists() and taqman_p.exists()):
        return None

    counts = pd.read_csv(counts_p, sep="\t")
    if "gene" not in counts.columns:
        raise ValueError("counts.tsv must have a 'gene' column")
    counts = counts.set_index("gene")
    counts.index = [str(g) for g in counts.index]
    counts = counts[~counts.index.duplicated(keep="first")]

    samples = pd.read_csv(samples_p, sep="\t")
    for col in ("sample", "group"):
        if col not in samples.columns:
            raise ValueError(f"samples.tsv must have a '{col}' column")
    samples["sample"] = samples["sample"].astype(str)
    samples["group"] = samples["group"].astype(str).str.upper()

    taqman = pd.read_csv(taqman_p, sep="\t")
    for col in ("gene", "log2_ab"):
        if col not in taqman.columns:
            raise ValueError(f"taqman.tsv must have a '{col}' column")
    taqman["gene"] = taqman["gene"].astype(str)
    taqman = taqman[~taqman["gene"].duplicated(keep="first")]
    taqman_log2 = {
        row.gene: _finite_float(row.log2_ab)
        for row in taqman.itertuples(index=False)
        if str(row.log2_ab) not in ("", "nan", "NA")
    }

    manifest = {}
    manifest_p = bundle / "manifest.json"
    if manifest_p.exists():
        try:
            manifest = json.loads(manifest_p.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            manifest = {}

    # Optional ERCC spike-in dose-response files (absent in older bundles).
    ercc_counts = ercc_truth = None
    ercc_counts_p = bundle / "ercc_counts.tsv"
    ercc_truth_p = bundle / "ercc_truth.tsv"
    if ercc_counts_p.exists() and ercc_truth_p.exists():
        ec = pd.read_csv(ercc_counts_p, sep="\t")
        if "ercc_id" in ec.columns:
            ec = ec.set_index("ercc_id")
            ec.index = [str(i) for i in ec.index]
            ercc_counts = ec[~ec.index.duplicated(keep="first")]
        et = pd.read_csv(ercc_truth_p, sep="\t")
        if "ercc_id" in et.columns:
            et["ercc_id"] = et["ercc_id"].astype(str)
            ercc_truth = et[~et["ercc_id"].duplicated(keep="first")].set_index("ercc_id")

    return {
        "counts": counts,
        "samples": samples,
        "taqman_log2_ab": taqman_log2,
        "ercc_counts": ercc_counts,
        "ercc_truth": ercc_truth,
        "manifest": manifest,
        "bundle_dir": str(bundle),
    }


def _group_means_cpm(counts, samples) -> dict[str, Any]:
    """Library-size-normalized (CPM) mean expression per group."""
    import numpy as np
    import pandas as pd

    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    means: dict[str, Any] = {}
    for grp, sub in samples.groupby("group"):
        cols = [s for s in sub["sample"].tolist() if s in cpm.columns]
        if cols:
            means[grp] = cpm[cols].mean(axis=1)
    return means


def score_seqc_maqc_a1(
    bundle: dict[str, Any],
    *,
    numerator: str = "A",
    denominator: str = "B",
    padj_max: float = 0.05,
    lfc_thr: float = 0.0,
    taqman_de_threshold: float = 1.0,
    taqman_null_threshold: float = 0.2,
    signal_floor: float = 0.5,
    min_lfc_pearson: float = 0.7,
    min_auc: float = 0.8,
) -> dict[str, Any]:
    """Run ARIA's real bulk DE on A vs B and score it against TaqMan truth.

    ``lfc_thr=0.0`` is used by default so the reference concordance is measured
    against ARIA's standard DESeq2-equivalent null (matched to TaqMan's plain
    log2(A/B)); the effect-size policy frontier is the synthetic A1 lane's job.
    """
    import numpy as np
    import pandas as pd
    from aria.scripts.rna_bulk_de import _run_deseq2

    counts = bundle["counts"]
    samples = bundle["samples"]
    taqman = bundle["taqman_log2_ab"]

    grp = samples.set_index("sample")["group"]
    keep = [s for s in counts.columns if s in grp.index and grp[s] in (numerator, denominator)]
    ab_counts = counts[keep]
    meta = pd.DataFrame({"group": [grp[s] for s in keep]}, index=keep)
    n_num = int((meta["group"] == numerator).sum())
    n_den = int((meta["group"] == denominator).sum())
    if n_num < 2 or n_den < 2:
        return {
            "status": "error",
            "messages": [
                f"need >=2 replicates per group; got {numerator}={n_num}, "
                f"{denominator}={n_den}"
            ],
        }

    de_result, warnings = _run_deseq2(
        ab_counts, meta, "group", numerator, denominator,
        padj_thr=padj_max, lfc_thr=lfc_thr, lfc_shrink=True,
    )
    if de_result.get("status") != "success":
        return {
            "status": "error",
            "messages": [
                f"bulk DE did not succeed: {de_result.get('error_type')} "
                f"{de_result.get('details', '')}".strip()
            ],
            "warnings": warnings,
        }

    df = de_result.get("results")
    df = df.copy()
    df.index = [str(i) for i in df.index]
    est_lfc = pd.to_numeric(df.get("log2FoldChange"), errors="coerce")
    pvalue = pd.to_numeric(df.get("pvalue"), errors="coerce")
    padj = pd.to_numeric(df.get("padj"), errors="coerce")

    # 1. LFC concordance vs TaqMan on the overlap.
    overlap = [g for g in df.index if g in taqman]
    conc = pd.DataFrame({
        "aria": est_lfc.reindex(overlap).values,
        "taqman": [taqman[g] for g in overlap],
    }, index=overlap).replace([np.inf, -np.inf], np.nan).dropna()
    if conc.empty:
        spearman = pearson = spearman_signal = 0.0
        n_signal = 0
    else:
        spearman = _finite_float(conc["aria"].corr(conc["taqman"], method="spearman"))
        pearson = _finite_float(conc["aria"].corr(conc["taqman"], method="pearson"))
        # Rank concordance restricted to genes TaqMan resolves as real signal:
        # the panel-wide Spearman is swamped by null ties when most genes are
        # non-DE, so report a signal-gene Spearman that is not.
        signal = conc[conc["taqman"].abs() >= signal_floor]
        n_signal = int(len(signal))
        spearman_signal = (
            _finite_float(signal["aria"].corr(signal["taqman"], method="spearman"))
            if n_signal >= 3 else 0.0
        )

    # 2. TaqMan-DE detection AUC + recall/empirical-FDR at the padj gate.
    score = (-np.log10(pvalue.clip(lower=1e-300))).reindex(overlap)
    pos = [g for g in overlap if abs(taqman[g]) >= taqman_de_threshold]
    neg = [g for g in overlap if abs(taqman[g]) <= taqman_null_threshold]
    auc = _auc(
        [float(score.get(g, np.nan)) for g in pos],
        [float(score.get(g, np.nan)) for g in neg],
    )
    called = set(map(str, de_result.get("sig_genes", []) or []))
    pos_set, neg_set = set(pos), set(neg)
    tp = called & pos_set
    fp = called & neg_set
    recall = len(tp) / max(len(pos_set), 1)
    emp_fdr_ref = len(fp) / max(len(called & (pos_set | neg_set)), 1)

    # 3. Titration monotonicity across A -> C -> D -> B (needs the mixtures).
    means = _group_means_cpm(counts, samples)
    have_titration = all(g in means for g in (numerator, "C", "D", denominator))
    if have_titration and pos:
        order = [numerator, "C", "D", denominator]
        seq = pd.DataFrame({g: means[g] for g in order})
        mono = 0
        for g in pos:
            if g not in seq.index:
                continue
            v = seq.loc[g].values.astype(float)
            if np.all(np.isfinite(v)) and (
                np.all(np.diff(v) >= -1e-9) or np.all(np.diff(v) <= 1e-9)
            ):
                mono += 1
        titration = {
            "status": "computed",
            "order": order,
            "n_evaluated": len(pos),
            "n_monotone": int(mono),
            "fraction_monotone": round(mono / max(len(pos), 1), 4),
        }
    else:
        titration = {
            "status": "not_computed",
            "reason": "titration mixtures C/D absent from the bundle"
            if not have_titration else "no TaqMan-DE positives",
        }

    axes = {
        "lfc_concordance": {
            "pearson": round(pearson, 4),
            "spearman": round(spearman, 4),
            "spearman_signal": round(spearman_signal, 4),
            "signal_floor": signal_floor,
            "n_genes_scored": int(len(conc)),
            "n_signal_genes": n_signal,
        },
        "taqman_de_detection": {
            "auc": None if auc is None else round(auc, 4),
            "n_positive": len(pos_set),
            "n_negative": len(neg_set),
            "recall_at_padj": round(recall, 4),
            "empirical_fdr_ref": round(emp_fdr_ref, 4),
            "n_called_total": len(called),
            "taqman_de_threshold": taqman_de_threshold,
            "taqman_null_threshold": taqman_null_threshold,
        },
        "titration_monotonicity": titration,
    }
    axis_pass = {
        "lfc_concordance": pearson >= min_lfc_pearson,
        "taqman_de_detection": auc is not None and auc >= min_auc,
    }
    status = "pass" if all(axis_pass.values()) else "fail"
    return {
        "status": status,
        "benchmark": "A1_seqc_maqc_reference",
        "benchmark_version": "v1",
        "scope": "external_reference_taqman_truth",
        "method_under_test": "ARIA bulk RNA DESeq2 path (_run_deseq2, apeGLM-enabled)",
        "comparison": {"numerator": numerator, "denominator": denominator},
        "samples": {
            "n_numerator": n_num,
            "n_denominator": n_den,
            "groups_present": sorted(set(grp.values)),
        },
        "thresholds": {"padj": padj_max, "lfc_threshold": lfc_thr},
        "tolerances": {"min_lfc_pearson": min_lfc_pearson, "min_auc": min_auc},
        "axis_pass": axis_pass,
        "axes": axes,
        "source": bundle.get("manifest", {}),
        "warnings": warnings,
        "messages": [
            "A1 reference lane validates ARIA's bulk DE against MAQC/SEQC TaqMan "
            "truth; TaqMan log2(A/B) and the count matrix come from the staged "
            "bundle, not from ARIA."
        ],
    }


def write_seqc_maqc_a1_figure(manifest: dict[str, Any], path: str) -> str:
    """Dependency-free SVG for the reference A1 figure."""
    axes = manifest.get("axes", {})
    det = axes.get("taqman_de_detection", {})
    titr = axes.get("titration_monotonicity", {})
    vals = [
        ("LFC Pearson vs TaqMan", _finite_float(axes.get("lfc_concordance", {}).get("pearson", 0.0))),
        ("TaqMan-DE AUC", _finite_float(det.get("auc", 0.0))),
        ("Titration monotone", _finite_float(titr.get("fraction_monotone", 0.0))
         if titr.get("status") == "computed" else 0.0),
    ]
    width = 760
    left, top = 290, 70
    bar_w, bar_h, gap = 330, 42, 26
    rows = []
    for i, (label, value) in enumerate(vals):
        y = top + i * (bar_h + gap)
        v = max(0.0, min(1.0, value))
        fill = "#2F6F73" if v >= 0.5 else "#A23B3B"
        shown = (
            "n/a" if (vals[i][0] == "Titration monotone" and titr.get("status") != "computed")
            else f"{v:.2f}"
        )
        rows.append(
            f'<text x="28" y="{y + 27}" font-size="16" fill="#1f2933">{label}</text>'
            f'<rect x="{left}" y="{y}" width="{bar_w}" height="{bar_h}" fill="#e8edf0"/>'
            f'<rect x="{left}" y="{y}" width="{bar_w * v:.1f}" height="{bar_h}" fill="{fill}"/>'
            f'<text x="{left + bar_w + 16}" y="{y + 27}" font-size="15" fill="#1f2933">{shown}</text>'
        )
    status = manifest.get("status", "unknown").upper()
    height = top + len(vals) * (bar_h + gap) + 40
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="#ffffff"/>'
        '<text x="28" y="34" font-size="21" font-weight="700" fill="#111827">'
        'Fig. 1 reference: A1 bulk DE vs MAQC/SEQC TaqMan</text>'
        f'<text x="650" y="34" font-size="17" font-weight="700" fill="#2F6F73">{status}</text>'
        f'<text x="28" y="{height - 14}" font-size="12" fill="#4b5563">'
        'External TaqMan qPCR truth; counts + truth from the staged reference bundle.</text>'
        + "".join(rows)
        + "</svg>"
    )
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg, encoding="utf-8")
    return str(out)


def run_seqc_maqc_a1_benchmark(
    bundle_dir: str | Path,
    *,
    output_dir: str | None = None,
    manifest_name: str = "a1_seqc_maqc_v4.5.5.json",
    figure_name: str = "fig1_a1_seqc_maqc_v4.5.5.svg",
    **score_kwargs: Any,
) -> dict[str, Any]:
    """Execute the A1 SEQC/MAQC reference lane, or skip honestly if no bundle."""
    import json

    bundle = load_seqc_maqc_bundle(bundle_dir)
    if bundle is None:
        return {
            "status": "skipped",
            "benchmark": "A1_seqc_maqc_reference",
            "reason": "no reference bundle staged",
            "expected_bundle": {
                "dir": str(bundle_dir),
                "files": ["counts.tsv", "samples.tsv", "taqman.tsv"],
                "bootstrap": "scripts/fetch_seqc_maqc_reference.py",
            },
            "messages": [
                "A1 SEQC/MAQC reference lane skipped: stage the bundle (see "
                "scripts/fetch_seqc_maqc_reference.py) then re-run. Nothing is "
                "fabricated."
            ],
        }

    manifest = score_seqc_maqc_a1(bundle, **score_kwargs)
    if output_dir:
        from aria.version import __version__, collect_version_metadata
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        manifest["aria_version"] = __version__
        manifest["provenance"] = collect_version_metadata()
        manifest_path = outdir / manifest_name
        figure_path = outdir / figure_name
        if manifest.get("status") in ("pass", "fail"):
            write_seqc_maqc_a1_figure(manifest, str(figure_path))
            manifest["artifacts"] = {
                "manifest_json": str(manifest_path),
                "figure_svg": str(figure_path),
            }
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
        )
    return manifest


def _site_de_lfc(bundle, numerator, denominator, padj_max, lfc_thr):
    """Run ARIA's real bulk DE for one site; return (de_result, est_lfc, n_num,
    n_den). Shared by the multi-site lane; the single-site scorer is unchanged."""
    import pandas as pd
    from aria.scripts.rna_bulk_de import _run_deseq2

    counts = bundle["counts"]
    grp = bundle["samples"].set_index("sample")["group"]
    keep = [s for s in counts.columns if s in grp.index and grp[s] in (numerator, denominator)]
    meta = pd.DataFrame({"group": [grp[s] for s in keep]}, index=keep)
    n_num = int((meta["group"] == numerator).sum())
    n_den = int((meta["group"] == denominator).sum())
    if n_num < 2 or n_den < 2:
        return None, None, n_num, n_den
    de_result, _w = _run_deseq2(
        counts[keep], meta, "group", numerator, denominator,
        padj_thr=padj_max, lfc_thr=lfc_thr, lfc_shrink=True,
    )
    if de_result.get("status") != "success":
        return de_result, None, n_num, n_den
    df = de_result["results"].copy()
    df.index = [str(i) for i in df.index]
    est_lfc = pd.to_numeric(df.get("log2FoldChange"), errors="coerce")
    est_lfc = est_lfc[est_lfc.notna()]
    return de_result, est_lfc, n_num, n_den


def _taqman_summary(de_result, est_lfc, taqman, *, signal_floor, taqman_de_threshold,
                    taqman_null_threshold):
    """Compact per-site TaqMan concordance (Pearson + signal Spearman + AUC)."""
    import numpy as np
    import pandas as pd

    df = de_result["results"].copy()
    df.index = [str(i) for i in df.index]
    pvalue = pd.to_numeric(df.get("pvalue"), errors="coerce")
    overlap = [g for g in est_lfc.index if g in taqman]
    conc = pd.DataFrame({
        "aria": est_lfc.reindex(overlap).values,
        "taqman": [taqman[g] for g in overlap],
    }, index=overlap).replace([np.inf, -np.inf], np.nan).dropna()
    if conc.empty:
        return {"pearson": 0.0, "spearman_signal": 0.0, "auc": None, "n_overlap": 0}
    pearson = _finite_float(conc["aria"].corr(conc["taqman"], method="pearson"))
    signal = conc[conc["taqman"].abs() >= signal_floor]
    spearman_signal = (
        _finite_float(signal["aria"].corr(signal["taqman"], method="spearman"))
        if len(signal) >= 3 else 0.0
    )
    score = (-np.log10(pvalue.clip(lower=1e-300))).reindex(overlap)
    pos = [g for g in overlap if abs(taqman[g]) >= taqman_de_threshold]
    neg = [g for g in overlap if abs(taqman[g]) <= taqman_null_threshold]
    auc = _auc(
        [float(score.get(g, np.nan)) for g in pos],
        [float(score.get(g, np.nan)) for g in neg],
    )
    return {
        "pearson": round(pearson, 4),
        "spearman_signal": round(spearman_signal, 4),
        "auc": None if auc is None else round(auc, 4),
        "n_overlap": int(len(conc)),
    }


def run_seqc_maqc_multisite(
    site_bundles: dict[str, str | Path],
    *,
    numerator: str = "A",
    denominator: str = "B",
    padj_max: float = 0.05,
    lfc_thr: float = 0.0,
    signal_floor: float = 0.5,
    taqman_de_threshold: float = 1.0,
    taqman_null_threshold: float = 0.2,
    min_cross_site_pearson: float = 0.9,
    output_dir: str | None = None,
    manifest_name: str = "a1_seqc_multisite_v4.5.5.json",
    figure_name: str = "fig1_a1_seqc_multisite_v4.5.5.svg",
) -> dict[str, Any]:
    """Cross-site SEQC reproducibility: run ARIA's bulk DE (A vs B) at each site
    and report the pairwise log2FC concordance between sites plus each site's
    TaqMan concordance. Cross-site log2FC correlation is the SEQC reproducibility
    metric; high values mean the DE result does not depend on the sequencing
    site. Sites with no staged bundle are skipped honestly."""
    import json
    import numpy as np
    import pandas as pd

    per_site: dict[str, Any] = {}
    lfc_by_site: dict[str, Any] = {}
    for site, bdir in site_bundles.items():
        bundle = load_seqc_maqc_bundle(bdir)
        if bundle is None:
            per_site[site] = {"status": "skipped", "reason": "no bundle staged"}
            continue
        de_result, est_lfc, n_num, n_den = _site_de_lfc(
            bundle, numerator, denominator, padj_max, lfc_thr
        )
        if est_lfc is None:
            per_site[site] = {
                "status": "error",
                "reason": (de_result or {}).get("error_type", "insufficient_replicates"),
                "n_numerator": n_num, "n_denominator": n_den,
            }
            continue
        summary = _taqman_summary(
            de_result, est_lfc, bundle["taqman_log2_ab"],
            signal_floor=signal_floor, taqman_de_threshold=taqman_de_threshold,
            taqman_null_threshold=taqman_null_threshold,
        )
        lfc_by_site[site] = est_lfc
        per_site[site] = {
            "status": "success",
            "n_numerator": n_num, "n_denominator": n_den,
            "n_genes_tested": int(len(est_lfc)),
            "taqman": summary,
        }

    # Cross-site pairwise log2FC concordance on commonly-tested genes.
    sites = sorted(lfc_by_site)
    pearson_matrix: dict[str, dict[str, float | None]] = {}
    spearman_matrix: dict[str, dict[str, float | None]] = {}
    offdiag_pearson: list[float] = []
    pair_n: list[int] = []
    for s1 in sites:
        pearson_matrix[s1] = {}
        spearman_matrix[s1] = {}
        for s2 in sites:
            common = lfc_by_site[s1].index.intersection(lfc_by_site[s2].index)
            joint = pd.DataFrame({
                "a": lfc_by_site[s1].reindex(common).values,
                "b": lfc_by_site[s2].reindex(common).values,
            }).replace([np.inf, -np.inf], np.nan).dropna()
            if len(joint) < 3:
                pearson_matrix[s1][s2] = spearman_matrix[s1][s2] = None
                continue
            p = round(_finite_float(joint["a"].corr(joint["b"], method="pearson")), 4)
            sp = round(_finite_float(joint["a"].corr(joint["b"], method="spearman")), 4)
            pearson_matrix[s1][s2] = p
            spearman_matrix[s1][s2] = sp
            if s1 < s2:
                offdiag_pearson.append(p)
                pair_n.append(int(len(joint)))

    cross_site = {
        "sites": sites,
        "n_sites": len(sites),
        "pearson_matrix": pearson_matrix,
        "spearman_matrix": spearman_matrix,
        "mean_offdiagonal_pearson": round(float(np.mean(offdiag_pearson)), 4) if offdiag_pearson else None,
        "min_offdiagonal_pearson": round(float(np.min(offdiag_pearson)), 4) if offdiag_pearson else None,
        "n_pairs": len(offdiag_pearson),
        "median_pair_genes": int(np.median(pair_n)) if pair_n else 0,
    }
    cross_ok = (
        cross_site["min_offdiagonal_pearson"] is not None
        and cross_site["min_offdiagonal_pearson"] >= min_cross_site_pearson
    )
    taqman_ok = all(
        v.get("taqman", {}).get("pearson", 0.0) >= 0.7
        for v in per_site.values() if v.get("status") == "success"
    )
    n_ok_sites = sum(1 for v in per_site.values() if v.get("status") == "success")
    status = "pass" if (cross_ok and taqman_ok and n_ok_sites >= 2) else (
        "fail" if n_ok_sites >= 2 else "skipped"
    )

    manifest = {
        "status": status,
        "benchmark": "A1_seqc_multisite_reproducibility",
        "benchmark_version": "v1",
        "scope": "external_reference_cross_site",
        "method_under_test": "ARIA bulk RNA DESeq2 path (_run_deseq2, apeGLM-enabled)",
        "comparison": {"numerator": numerator, "denominator": denominator},
        "tolerances": {"min_cross_site_pearson": min_cross_site_pearson},
        "per_site": per_site,
        "cross_site": cross_site,
        "messages": [
            "Cross-site SEQC reproducibility: ARIA's A-vs-B log2FC is correlated "
            "between sequencing sites; counts + TaqMan truth come from the staged "
            "per-site bundles, not from ARIA."
        ],
    }
    if output_dir and n_ok_sites >= 2:
        from aria.version import __version__, collect_version_metadata
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        manifest["aria_version"] = __version__
        manifest["provenance"] = collect_version_metadata()
        figure_path = outdir / figure_name
        manifest_path = outdir / manifest_name
        write_seqc_multisite_figure(manifest, str(figure_path))
        manifest["artifacts"] = {
            "manifest_json": str(manifest_path),
            "figure_svg": str(figure_path),
        }
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
        )
    return manifest


def write_seqc_multisite_figure(manifest: dict[str, Any], path: str) -> str:
    """Dependency-free SVG: per-site TaqMan Pearson + the cross-site mean."""
    per_site = manifest.get("per_site", {})
    cross = manifest.get("cross_site", {})
    rows_data = [
        (site, _finite_float(v.get("taqman", {}).get("pearson", 0.0)))
        for site, v in sorted(per_site.items())
        if v.get("status") == "success"
    ]
    mean_cross = cross.get("mean_offdiagonal_pearson")
    if mean_cross is not None:
        rows_data.append((f"cross-site mean ({cross.get('n_pairs', 0)} pairs)",
                          _finite_float(mean_cross)))

    width = 760
    left, top = 320, 70
    bar_w, bar_h, gap = 300, 38, 22
    rows = []
    for i, (label, value) in enumerate(rows_data):
        y = top + i * (bar_h + gap)
        v = max(0.0, min(1.0, value))
        fill = "#2F6F73" if "cross-site" in label else "#3b6ea2"
        if v < 0.5:
            fill = "#A23B3B"
        rows.append(
            f'<text x="28" y="{y + 25}" font-size="15" fill="#1f2933">{label}</text>'
            f'<rect x="{left}" y="{y}" width="{bar_w}" height="{bar_h}" fill="#e8edf0"/>'
            f'<rect x="{left}" y="{y}" width="{bar_w * v:.1f}" height="{bar_h}" fill="{fill}"/>'
            f'<text x="{left + bar_w + 14}" y="{y + 25}" font-size="14" fill="#1f2933">{v:.3f}</text>'
        )
    status = manifest.get("status", "unknown").upper()
    height = top + len(rows_data) * (bar_h + gap) + 36
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="#ffffff"/>'
        '<text x="28" y="34" font-size="20" font-weight="700" fill="#111827">'
        'Fig. 1 reference: SEQC cross-site reproducibility (A vs B)</text>'
        f'<text x="660" y="34" font-size="16" font-weight="700" fill="#2F6F73">{status}</text>'
        f'<text x="28" y="{height - 12}" font-size="12" fill="#4b5563">'
        'Per-site log2FC vs TaqMan (blue) and mean pairwise cross-site log2FC concordance (teal).</text>'
        + "".join(rows)
        + "</svg>"
    )
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg, encoding="utf-8")
    return str(out)


def _ols_slope(x, y) -> float | None:
    """Least-squares slope of y on x (no intercept assumption)."""
    import numpy as np

    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 3 or np.allclose(x.var(), 0.0):
        return None
    return float(np.polyfit(x, y, 1)[0])


def score_ercc_dose_response(
    bundle: dict[str, Any],
    *,
    numerator: str = "A",
    denominator: str = "B",
    min_count_sum: int = 10,
    min_fc_pearson: float = 0.5,
    min_dynamic_range_pearson: float = 0.9,
) -> dict[str, Any]:
    """Score ARIA's recovery of the ERCC spike-in dose-response from the staged
    bundle: (1) fold-change recovery — measured log2(A/B) of each ERCC vs the
    known Mix1/Mix2 log2 ratio, per subgroup; (2) dynamic-range linearity —
    measured CPM vs the known input concentration across orders of magnitude.

    ERCC counts (per sample) and the ERCC concentration/fold-change truth come
    entirely from the bundle. Mix 1 is assumed spiked into the numerator group
    (SEQC convention); the fitted slope sign reports whether that holds."""
    import numpy as np
    import pandas as pd

    ercc = bundle.get("ercc_counts")
    truth = bundle.get("ercc_truth")
    if ercc is None or truth is None:
        return {"status": "skipped", "reason": "no ERCC files in bundle"}

    gene_counts = bundle["counts"]
    grp = bundle["samples"].set_index("sample")["group"]
    a_cols = [s for s in ercc.columns if s in grp.index and grp[s] == numerator]
    b_cols = [s for s in ercc.columns if s in grp.index and grp[s] == denominator]
    if len(a_cols) < 2 or len(b_cols) < 2:
        return {"status": "error", "messages": ["need >=2 reps per group for ERCC"]}

    # CPM normalize ERCC by the gene-matrix library size (the main library).
    cpm = pd.DataFrame(index=ercc.index)
    for s in a_cols + b_cols:
        lib = float(gene_counts[s].sum()) or np.nan
        cpm[s] = ercc[s] / lib * 1e6
    mean_a = cpm[a_cols].mean(axis=1)
    mean_b = cpm[b_cols].mean(axis=1)
    total = ercc[a_cols + b_cols].sum(axis=1)

    common = [e for e in ercc.index if e in truth.index]
    detected = [
        e for e in common
        if total.get(e, 0) >= min_count_sum
        and mean_a.get(e, 0) > 0 and mean_b.get(e, 0) > 0
    ]

    # 1. Fold-change recovery vs the known Mix1/Mix2 log2 ratio.
    measured_fc = {e: float(np.log2(mean_a[e] / mean_b[e])) for e in detected}
    expected_fc = {e: float(truth.loc[e, "log2_mix1_mix2"]) for e in detected}
    mvec = [measured_fc[e] for e in detected]
    evec = [expected_fc[e] for e in detected]
    fc_pearson = (
        _finite_float(pd.Series(mvec).corr(pd.Series(evec))) if len(detected) >= 3 else 0.0
    )
    slope = _ols_slope(evec, mvec)
    by_subgroup = {}
    sub = truth["subgroup"].astype(str)
    for g in sorted(set(sub.loc[detected])):
        ids = [e for e in detected if sub.get(e) == g]
        by_subgroup[g] = {
            "n": len(ids),
            "expected_log2": round(float(np.mean([expected_fc[e] for e in ids])), 3),
            "measured_log2_mean": round(float(np.mean([measured_fc[e] for e in ids])), 3),
            "measured_log2_sd": round(float(np.std([measured_fc[e] for e in ids])), 3),
        }

    # 2. Dynamic-range linearity: measured CPM vs known input concentration.
    # Each (sample, ERCC) uses its mix's nominal concentration.
    log_meas, log_nom = [], []
    for e in detected:
        c1 = float(truth.loc[e, "conc_mix1"])
        c2 = float(truth.loc[e, "conc_mix2"])
        for s in a_cols:
            if cpm.loc[e, s] > 0 and c1 > 0:
                log_meas.append(np.log2(cpm.loc[e, s])); log_nom.append(np.log2(c1))
        for s in b_cols:
            if cpm.loc[e, s] > 0 and c2 > 0:
                log_meas.append(np.log2(cpm.loc[e, s])); log_nom.append(np.log2(c2))
    if len(log_nom) >= 3:
        dr_pearson = _finite_float(pd.Series(log_meas).corr(pd.Series(log_nom)))
        dr_slope = _ols_slope(log_nom, log_meas)
        dynamic_range_log10 = round(
            float((np.max(log_nom) - np.min(log_nom)) / np.log2(10)), 2
        )
    else:
        dr_pearson = 0.0
        dr_slope = None
        dynamic_range_log10 = 0.0

    axes = {
        "fold_change_recovery": {
            "pearson": round(fc_pearson, 4),
            "slope_measured_vs_expected": None if slope is None else round(slope, 3),
            "n_ercc_detected": len(detected),
            "n_ercc_total": int(len(ercc.index)),
            "by_subgroup": by_subgroup,
        },
        "dynamic_range": {
            "pearson_log_cpm_vs_log_conc": round(dr_pearson, 4),
            "slope": None if dr_slope is None else round(dr_slope, 3),
            "dynamic_range_orders_of_magnitude": dynamic_range_log10,
            "n_points": len(log_nom),
        },
    }
    axis_pass = {
        "fold_change_recovery": abs(fc_pearson) >= min_fc_pearson,
        "dynamic_range": dr_pearson >= min_dynamic_range_pearson,
    }
    status = "pass" if all(axis_pass.values()) else "fail"
    return {
        "status": status,
        "benchmark": "A1_ercc_dose_response",
        "benchmark_version": "v1",
        "scope": "external_reference_ercc",
        "method_under_test": "ARIA CPM normalization on the SEQC gene library",
        "comparison": {"numerator": numerator, "denominator": denominator,
                       "mix_assumption": f"{numerator}=Mix1, {denominator}=Mix2"},
        "tolerances": {"min_fc_pearson": min_fc_pearson,
                       "min_dynamic_range_pearson": min_dynamic_range_pearson},
        "axis_pass": axis_pass,
        "axes": axes,
        "source": bundle.get("manifest", {}),
        "messages": [
            "ERCC dose-response: measured log2(A/B) vs known Mix1/Mix2 ratio "
            "(per subgroup) and measured CPM vs known input concentration. ERCC "
            "counts + concentration truth come from the staged bundle."
        ],
    }


def write_ercc_figure(manifest: dict[str, Any], path: str) -> str:
    """Dependency-free SVG: per-subgroup expected vs measured ERCC log2FC."""
    axes = manifest.get("axes", {})
    fc = axes.get("fold_change_recovery", {})
    dr = axes.get("dynamic_range", {})
    by = fc.get("by_subgroup", {})
    rows_data = [
        (f"subgroup {g} (exp {v['expected_log2']:+.2f})", v["measured_log2_mean"],
         v["expected_log2"])
        for g, v in sorted(by.items())
    ]
    width = 760
    left, mid, top = 300, 380, 80
    row_h, gap = 36, 20
    span = 3.2  # log2 axis half-range for the centered bars
    rows = []
    for i, (label, measured, expected) in enumerate(rows_data):
        y = top + i * (row_h + gap)
        mx = mid + (max(-span, min(span, measured)) / span) * (width - mid - 90)
        ex = mid + (max(-span, min(span, expected)) / span) * (width - mid - 90)
        rows.append(
            f'<text x="20" y="{y + 23}" font-size="14" fill="#1f2933">{label}</text>'
            f'<line x1="{ex:.1f}" y1="{y}" x2="{ex:.1f}" y2="{y + row_h}" '
            f'stroke="#A23B3B" stroke-width="3"/>'
            f'<circle cx="{mx:.1f}" cy="{y + row_h / 2:.1f}" r="7" fill="#2F6F73"/>'
            f'<text x="{mx + 12:.1f}" y="{y + 23}" font-size="13" fill="#1f2933">'
            f'{measured:+.2f}</text>'
        )
    status = manifest.get("status", "unknown").upper()
    height = top + len(rows_data) * (row_h + gap) + 56
    yaxis = mid
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="#ffffff"/>'
        '<text x="20" y="34" font-size="20" font-weight="700" fill="#111827">'
        'Fig. 1 reference: ERCC dose-response (A vs B)</text>'
        f'<text x="640" y="34" font-size="16" font-weight="700" fill="#2F6F73">{status}</text>'
        f'<line x1="{yaxis}" y1="{top - 10}" x2="{yaxis}" y2="{height - 50}" '
        'stroke="#cbd5e1" stroke-width="1"/>'
        f'<text x="20" y="58" font-size="12" fill="#4b5563">'
        f'Red bar = expected log2(Mix1/Mix2); teal dot = measured. FC Pearson '
        f'{fc.get("pearson", 0):.3f}, slope {fc.get("slope_measured_vs_expected")}; '
        f'dynamic-range Pearson {dr.get("pearson_log_cpm_vs_log_conc", 0):.3f} over '
        f'{dr.get("dynamic_range_orders_of_magnitude", 0)} log10.</text>'
        + "".join(rows)
        + "</svg>"
    )
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg, encoding="utf-8")
    return str(out)


def run_ercc_dose_response(
    bundle_dir: str | Path,
    *,
    output_dir: str | None = None,
    manifest_name: str = "a1_ercc_dose_response_v4.5.5.json",
    figure_name: str = "fig1_a1_ercc_dose_response_v4.5.5.svg",
    **score_kwargs: Any,
) -> dict[str, Any]:
    """Execute the ERCC dose-response lane, or skip honestly if no ERCC bundle."""
    import json

    bundle = load_seqc_maqc_bundle(bundle_dir)
    if bundle is None:
        return {"status": "skipped", "reason": "no bundle staged",
                "benchmark": "A1_ercc_dose_response"}
    manifest = score_ercc_dose_response(bundle, **score_kwargs)
    if output_dir and manifest.get("status") in ("pass", "fail"):
        from aria.version import __version__, collect_version_metadata
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        manifest["aria_version"] = __version__
        manifest["provenance"] = collect_version_metadata()
        figure_path = outdir / figure_name
        manifest_path = outdir / manifest_name
        write_ercc_figure(manifest, str(figure_path))
        manifest["artifacts"] = {"manifest_json": str(manifest_path),
                                 "figure_svg": str(figure_path)}
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
        )
    return manifest
