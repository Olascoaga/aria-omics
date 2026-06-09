"""Benchmark A2 external reference lane: ARIA vs muscat on the Kang dataset.

A2 (Fig 2) is donor-aware pseudobulk scRNA DE. The preliminary lane proved it on
synthetic donor data. This module compares ARIA's pseudobulk DE against the
reference method *muscat* (Crowell et al., Nat Commun 2020) on the real *Kang
et al. 2018* PBMC dataset (8 donors, control vs IFN-beta stimulation) — the
canonical multi-sample scRNA-DE benchmark.

Apples-to-apples like the A1 external lane: muscat aggregates cells to per-cluster
pseudobulk (genes x donor-condition samples) and runs its DE; ARIA's real bulk DE
core (`_run_deseq2`) runs on the *same* exported pseudobulk matrices. The
comparison is therefore of the DE statistics on identical donor-level inputs, not
of two different aggregations. ARIA fabricates nothing: the Kang counts and the
muscat reference results come entirely from the staged R export.

Export schema (written by aria/scripts/benchmark_a2_external_muscat.R):
- ``samples.tsv``       : columns ``sample`` (= sample_id), ``group`` (ctrl/stim).
- ``pb_<cluster>.tsv``  : column ``gene`` + one column per sample (pseudobulk sums).
- ``muscat_<cluster>.tsv``: columns ``gene``, ``logFC``, ``p_val``, ``p_adj``.
- ``clusters.json``     : clusters, contrast, dataset provenance.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

__all__ = [
    "score_aria_vs_muscat",
    "run_kang_muscat_benchmark",
    "write_kang_muscat_figure",
]


def _finite(value: Any) -> float:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return 0.0
    return f if f == f and f not in (float("inf"), float("-inf")) else 0.0


def _jaccard(a: set, b: set) -> float:
    return len(a & b) / max(len(a | b), 1)


def _score_cluster(pb_path, muscat_path, samples, *, numerator, denominator,
                   padj_max, top_k):
    import numpy as np
    import pandas as pd
    from aria.scripts.rna_bulk_de import _run_deseq2

    pb = pd.read_csv(pb_path, sep="\t")
    if "gene" not in pb.columns:
        return {"status": "error", "reason": "pb matrix lacks 'gene'"}
    pb = pb.set_index("gene")
    pb.index = [str(g) for g in pb.index]
    pb = pb[~pb.index.duplicated(keep="first")]

    grp = samples.set_index("sample")["group"]
    keep = [s for s in pb.columns if s in grp.index and grp[s] in (numerator, denominator)]
    meta = pd.DataFrame({"group": [grp[s] for s in keep]}, index=keep)
    n_num = int((meta["group"] == numerator).sum())
    n_den = int((meta["group"] == denominator).sum())
    if n_num < 2 or n_den < 2:
        return {"status": "skipped", "reason": f"too few samples ({n_num}/{n_den})"}

    de, _w = _run_deseq2(
        pb[keep].round().astype(int), meta, "group", numerator, denominator,
        padj_thr=padj_max, lfc_thr=0.0, lfc_shrink=True,
    )
    if de.get("status") != "success":
        return {"status": "error", "reason": de.get("error_type", "aria_de_failed")}

    a = de["results"].copy()
    a.index = [str(i) for i in a.index]
    aria_lfc = pd.to_numeric(a.get("log2FoldChange"), errors="coerce")
    aria_padj = pd.to_numeric(a.get("padj"), errors="coerce")
    aria_sig = set(map(str, de.get("sig_genes", []) or []))

    m = pd.read_csv(muscat_path, sep="\t")
    m["gene"] = m["gene"].astype(str)
    m = m.set_index("gene")
    mus_lfc = pd.to_numeric(m.get("logFC"), errors="coerce")
    mus_padj = pd.to_numeric(m.get("p_adj"), errors="coerce")
    mus_sig = set(mus_padj[mus_padj < padj_max].index.astype(str))

    common = [g for g in a.index if g in m.index]
    conc = pd.DataFrame({
        "aria": aria_lfc.reindex(common).values,
        "muscat": mus_lfc.reindex(common).values,
    }, index=common).replace([np.inf, -np.inf], np.nan).dropna()
    spearman = _finite(conc["aria"].corr(conc["muscat"], method="spearman")) if len(conc) >= 3 else 0.0
    pearson = _finite(conc["aria"].corr(conc["muscat"], method="pearson")) if len(conc) >= 3 else 0.0
    # Rank concordance restricted to genes muscat calls significant: the
    # panel-wide Spearman is swamped by the non-DE majority (noise vs noise),
    # so a signal-gene Spearman reflects real agreement on the DE genes.
    sig_common = [g for g in conc.index if g in mus_sig]
    spearman_sig = (
        _finite(conc.loc[sig_common, "aria"].corr(conc.loc[sig_common, "muscat"],
                                                  method="spearman"))
        if len(sig_common) >= 3 else 0.0
    )

    # Ranking concordance on the top-k by significance.
    aria_rank = aria_padj.reindex(common).fillna(1.0).sort_values().index.tolist()
    mus_rank = mus_padj.reindex(common).fillna(1.0).sort_values().index.tolist()
    k = min(top_k, len(common))
    top_jacc = _jaccard(set(aria_rank[:k]), set(mus_rank[:k])) if k else 0.0

    inter = aria_sig & mus_sig
    # Direction agreement on shared significant genes.
    shared = [g for g in inter if g in conc.index]
    dir_agree = (
        sum(1 for g in shared if np.sign(conc.loc[g, "aria"]) == np.sign(conc.loc[g, "muscat"]))
        / max(len(shared), 1)
    ) if shared else 0.0

    return {
        "status": "success",
        "n_samples": {"numerator": n_num, "denominator": n_den},
        "n_genes_common": int(len(conc)),
        "lfc_pearson": round(pearson, 4),
        "lfc_spearman": round(spearman, 4),
        "lfc_spearman_sig": round(spearman_sig, 4),
        "top_k_jaccard": round(top_jacc, 4),
        "n_aria_sig": len(aria_sig),
        "n_muscat_sig": len(mus_sig),
        "n_shared_sig": len(inter),
        "sig_jaccard": round(_jaccard(aria_sig, mus_sig), 4),
        "shared_sig_direction_agreement": round(dir_agree, 4),
    }


def score_aria_vs_muscat(
    export_dir: str | Path,
    *,
    numerator: str = "stim",
    denominator: str = "ctrl",
    padj_max: float = 0.05,
    top_k: int = 100,
    min_mean_lfc_pearson: float = 0.7,
    min_mean_sig_jaccard: float = 0.3,
) -> dict[str, Any]:
    """Score ARIA's pseudobulk DE against muscat per cluster on the Kang export."""
    import numpy as np
    import pandas as pd

    export = Path(export_dir)
    clusters_p = export / "clusters.json"
    samples_p = export / "samples.tsv"
    if not (clusters_p.exists() and samples_p.exists()):
        return {"status": "skipped", "reason": "no muscat export staged",
                "benchmark": "A2_kang_muscat"}

    meta = json.loads(clusters_p.read_text(encoding="utf-8"))
    samples = pd.read_csv(samples_p, sep="\t")
    samples["sample"] = samples["sample"].astype(str)
    samples["group"] = samples["group"].astype(str)

    per_cluster: dict[str, Any] = {}
    for cl in meta.get("clusters", []):
        pb_p = export / f"pb_{cl}.tsv"
        mus_p = export / f"muscat_{cl}.tsv"
        if not (pb_p.exists() and mus_p.exists()):
            per_cluster[cl] = {"status": "skipped", "reason": "missing cluster files"}
            continue
        per_cluster[cl] = _score_cluster(
            pb_p, mus_p, samples,
            numerator=numerator, denominator=denominator,
            padj_max=padj_max, top_k=top_k,
        )

    ok = [v for v in per_cluster.values() if v.get("status") == "success"]
    if not ok:
        return {"status": "error", "benchmark": "A2_kang_muscat",
                "reason": "no cluster scored", "per_cluster": per_cluster}

    summary = {
        "n_clusters_scored": len(ok),
        "mean_lfc_pearson": round(float(np.mean([v["lfc_pearson"] for v in ok])), 4),
        "min_lfc_pearson": round(float(np.min([v["lfc_pearson"] for v in ok])), 4),
        "mean_lfc_spearman": round(float(np.mean([v["lfc_spearman"] for v in ok])), 4),
        "mean_lfc_spearman_sig": round(float(np.mean([v["lfc_spearman_sig"] for v in ok])), 4),
        "mean_sig_jaccard": round(float(np.mean([v["sig_jaccard"] for v in ok])), 4),
        "mean_top_k_jaccard": round(float(np.mean([v["top_k_jaccard"] for v in ok])), 4),
        "mean_shared_sig_direction_agreement": round(
            float(np.mean([v["shared_sig_direction_agreement"] for v in ok])), 4),
        "total_shared_sig": int(sum(v["n_shared_sig"] for v in ok)),
    }
    axis_pass = {
        "lfc_concordance": summary["mean_lfc_pearson"] >= min_mean_lfc_pearson,
        "sig_overlap": summary["mean_sig_jaccard"] >= min_mean_sig_jaccard,
    }
    status = "pass" if all(axis_pass.values()) else "fail"
    return {
        "status": status,
        "benchmark": "A2_kang_muscat",
        "benchmark_version": "v1",
        "scope": "external_reference_kang_muscat",
        "method_under_test": "ARIA pseudobulk DE (_run_deseq2) vs muscat pbDS",
        "comparison": {"numerator": numerator, "denominator": denominator},
        "dataset": meta.get("dataset", "muscat::example_sce (Kang 2018)"),
        "tolerances": {"min_mean_lfc_pearson": min_mean_lfc_pearson,
                       "min_mean_sig_jaccard": min_mean_sig_jaccard},
        "axis_pass": axis_pass,
        "summary": summary,
        "per_cluster": per_cluster,
        "messages": [
            "A2 external: ARIA's pseudobulk DE vs muscat on identical per-cluster "
            "pseudobulk from the real Kang dataset. Counts + muscat reference come "
            "from the staged export, not from ARIA."
        ],
    }


def write_kang_muscat_figure(manifest: dict[str, Any], path: str) -> str:
    """Dependency-free SVG: per-cluster ARIA-vs-muscat LFC Spearman + sig Jaccard."""
    per = manifest.get("per_cluster", {})
    rows_data = [
        (cl, _finite(v.get("lfc_pearson", 0.0)), _finite(v.get("sig_jaccard", 0.0)))
        for cl, v in sorted(per.items()) if v.get("status") == "success"
    ]
    width = 780
    left, bar_w, bar_h, gap, top = 230, 230, 26, 16, 78
    rows = []
    for i, (label, sp, jc) in enumerate(rows_data):
        y = top + i * (bar_h + gap)
        rows.append(
            f'<text x="16" y="{y + 18}" font-size="13" fill="#1f2933">{label[:26]}</text>'
            f'<rect x="{left}" y="{y}" width="{bar_w}" height="{bar_h}" fill="#e8edf0"/>'
            f'<rect x="{left}" y="{y}" width="{bar_w * max(0, min(1, sp)):.1f}" height="{bar_h}" fill="#2F6F73"/>'
            f'<text x="{left + bar_w + 8}" y="{y + 18}" font-size="12" fill="#1f2933">r={sp:.2f}</text>'
            f'<rect x="{left + bar_w + 70}" y="{y}" width="{bar_w}" height="{bar_h}" fill="#e8edf0"/>'
            f'<rect x="{left + bar_w + 70}" y="{y}" width="{bar_w * max(0, min(1, jc)):.1f}" height="{bar_h}" fill="#3b6ea2"/>'
            f'<text x="{left + 2 * bar_w + 78}" y="{y + 18}" font-size="12" fill="#1f2933">J={jc:.2f}</text>'
        )
    s = manifest.get("summary", {})
    status = manifest.get("status", "unknown").upper()
    height = top + len(rows_data) * (bar_h + gap) + 44
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="#ffffff"/>'
        '<text x="16" y="32" font-size="19" font-weight="700" fill="#111827">'
        'Fig. 2 reference: ARIA vs muscat pseudobulk DE (Kang)</text>'
        f'<text x="660" y="32" font-size="16" font-weight="700" fill="#2F6F73">{status}</text>'
        f'<text x="16" y="54" font-size="12" fill="#4b5563">'
        f'Teal = log2FC Pearson (ARIA vs muscat); blue = significant-gene Jaccard. '
        f'Mean r={s.get("mean_lfc_pearson", 0)}, mean J={s.get("mean_sig_jaccard", 0)}, '
        f'{s.get("total_shared_sig", 0)} shared significant.</text>'
        + "".join(rows)
        + "</svg>"
    )
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg, encoding="utf-8")
    return str(out)


def run_kang_muscat_benchmark(
    export_dir: str | Path,
    *,
    output_dir: str | None = None,
    manifest_name: str = "a2_kang_muscat_v4.5.5.json",
    figure_name: str = "fig2_a2_kang_muscat_v4.5.5.svg",
    **score_kwargs: Any,
) -> dict[str, Any]:
    """Score the staged muscat/Kang export, or skip honestly if absent."""
    manifest = score_aria_vs_muscat(export_dir, **score_kwargs)
    if output_dir and manifest.get("status") in ("pass", "fail"):
        from aria.version import __version__, collect_version_metadata
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        manifest["aria_version"] = __version__
        manifest["provenance"] = collect_version_metadata()
        figure_path = outdir / figure_name
        manifest_path = outdir / manifest_name
        write_kang_muscat_figure(manifest, str(figure_path))
        manifest["artifacts"] = {"manifest_json": str(manifest_path),
                                 "figure_svg": str(figure_path)}
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
        )
    return manifest
