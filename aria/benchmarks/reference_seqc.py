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

    return {
        "counts": counts,
        "samples": samples,
        "taqman_log2_ab": taqman_log2,
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
