"""scATAC P3b — SnapATAC2 reference driver (PRODUCES the outputs P3a scores).

This is the drop-in companion to the P3a concordance harness
(``benchmark_scatac_concordance.py``). P3a CONSUMES a reference tool's outputs;
this driver PRODUCES them by running SnapATAC2 inside ``aria-bench-atac-env`` and
emitting the same JSON schema P3a reads:

    {"clusters": [<per-cell label>, ...],   # ALIGNED to ARIA's cell order
     "da_peaks": ["chr1:100-200", ...],     # reference DA-peak ids/coords
     "motifs":   ["MOTIF_A", ...]}          # reference motif ranking

The reference pipeline (spectral LSI -> leiden -> marker peaks -> motif
enrichment) runs over the SAME peak-matrix ``.h5ad`` ARIA clustered, so the cell
order matches by construction and cluster ARI/NMI compares the same cells. The
SnapATAC2 stack is installed in a dedicated env (``aria-bench-atac-env``), kept
separate from the R/Bioconductor ``aria-bench-env`` so SnapATAC2's ``numpy>=2``
does not disturb the validated RNA comparator lock.

No-fabrication contract (ADR-002): with SnapATAC2 absent or the input missing the
driver returns a structured honest ``not_run`` instead of inventing reference
outputs. The motif axis degrades to ``motifs=None`` + a reason (never a fake
ranking) when a local genome FASTA or the motif database is unavailable; ARIA's
motif lane uses JASPAR2024 while SnapATAC2's bundled database is CIS-BP, so the
exact-name motif concordance is a conservative lower bound (TF-name level) and is
labelled as such in ``notes``.

scATAC stays ``alpha`` + ``requires_ack``. This driver is benchmark scaffolding,
never a dispatch path for a user run and never a modality promotion.
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any


_DEFAULT_GENOME_FASTA = Path(os.path.expanduser("~/.aria/genomes/hg38/genome.fa"))


def _ranked_motifs_from_enrichment(enrichment: dict[str, Any]) -> list[str]:
    """Aggregate SnapATAC2 per-group motif enrichment into one ranked TF list.

    ``snap.tl.motif_enrichment`` returns ``{group: polars.DataFrame}``. We keep,
    per motif, its strongest signal across groups (smallest adjusted p-value, then
    largest log2 fold-change), then rank most-enriched first. TF/gene names are
    upper-cased so a CIS-BP entry can match a JASPAR2024 TF name when they agree.
    """
    best: dict[str, tuple[float, float]] = {}
    for df in enrichment.values():
        cols = {c.lower(): c for c in df.columns}
        name_col = cols.get("name") or cols.get("id") or list(df.columns)[0]
        padj_col = (cols.get("adjusted p-value") or cols.get("adjusted_pvalue")
                    or cols.get("adjusted p value") or cols.get("p-value")
                    or cols.get("pvalue"))
        fc_col = (cols.get("log2(fold change)") or cols.get("log2_fold_change")
                  or cols.get("fold change") or cols.get("log2fc"))
        for row in df.iter_rows(named=True):
            raw = row.get(name_col)
            if raw is None:
                continue
            name = str(raw).split("_")[0].split("+")[0].strip().upper()
            if not name:
                continue
            padj = float(row.get(padj_col, 1.0)) if padj_col else 1.0
            fc = float(row.get(fc_col, 0.0)) if fc_col else 0.0
            prev = best.get(name)
            if prev is None or (padj, -fc) < (prev[0], -prev[1]):
                best[name] = (padj, fc)
    return [name for name, _ in sorted(best.items(), key=lambda kv: (kv[1][0], -kv[1][1]))]


def _run_reference_pipeline(
    input_path: Path, params: dict[str, Any]
) -> dict[str, Any]:
    """Run SnapATAC2 on ARIA's peak matrix and return the P3a output schema.

    Cell order is preserved end-to-end (no cell filtering; the matrix is already
    ARIA-QC'd and clustered), so ``clusters`` aligns to ARIA's ``clusters`` by
    construction. ``da_peaks`` is the union of per-cluster marker peaks; ``motifs``
    is a best-effort ranked TF list (honest ``None`` when genome/motif assets are
    unavailable).
    """
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp
    import snapatac2 as snap

    resolution = float(params.get("resolution", 1.0))
    random_state = int(params.get("random_state", 0))
    n_features = int(params.get("n_features", 50000))
    marker_pvalue = float(params.get("marker_pvalue", 0.01))

    src = ad.read_h5ad(str(input_path))
    X = src.X
    X = (X if sp.issparse(X) else sp.csr_matrix(X)).tocsr().astype(np.float32)
    obs_names = [str(x) for x in src.obs_names]
    var_names = [str(x) for x in src.var_names]

    work_dir = params.get("work_dir") or tempfile.mkdtemp(prefix="p3b_snap_")
    backing = Path(work_dir) / "snapatac2_reference.h5ad"
    if backing.exists():
        backing.unlink()
    adata = snap.AnnData(X=X, filename=str(backing))
    adata.obs_names = obs_names
    adata.var_names = var_names

    snap.pp.select_features(adata, n_features=n_features, verbose=False)
    snap.tl.spectral(adata, n_comps=30, random_state=random_state)
    snap.pp.knn(adata, n_neighbors=50, random_state=random_state)
    snap.tl.leiden(adata, resolution=resolution, random_state=random_state)
    clusters = [str(x) for x in adata.obs["leiden"]]
    if len(clusters) != len(obs_names):
        raise RuntimeError(
            "SnapATAC2 changed the cell count/order; clusters would not align "
            f"to ARIA ({len(clusters)} vs {len(obs_names)} cells).")

    markers = snap.tl.marker_regions(adata, groupby="leiden", pvalue=marker_pvalue)
    da_peaks = sorted({p for peaks in markers.values() for p in peaks})

    motifs: list[str] | None = None
    motif_note: str
    genome_fasta = params.get("genome_fasta")
    fasta_path = Path(genome_fasta) if genome_fasta else _DEFAULT_GENOME_FASTA
    if params.get("enable_motifs", True) is False:
        motif_note = "motif enrichment disabled via params (enable_motifs=False)."
    elif not fasta_path.exists():
        motif_note = (f"motif enrichment skipped: genome FASTA not found at "
                      f"{fasta_path} (no auto-download; honest None per ADR-002).")
    elif not markers:
        motif_note = "motif enrichment skipped: no marker regions to test."
    else:
        try:
            motif_db = snap.datasets.cis_bp(unique=True)
            enrichment = snap.tl.motif_enrichment(
                motifs=motif_db, regions=markers, genome_fasta=fasta_path)
            motifs = _ranked_motifs_from_enrichment(enrichment)
            motif_note = (
                "SnapATAC2 motif enrichment uses CIS-BP; ARIA uses JASPAR2024. "
                "Names are upper-cased TF symbols so concordance is a conservative "
                "lower bound at the TF-name level (database namespaces differ).")
        except Exception as exc:  # noqa: BLE001 - honest degradation, never fake
            motifs = None
            motif_note = (f"motif enrichment failed ({type(exc).__name__}: {exc}); "
                          "reported as honest None rather than a fabricated ranking.")

    return {
        "clusters": clusters,
        "da_peaks": da_peaks,
        "motifs": motifs,
        "notes": {
            "n_cells": len(clusters),
            "n_clusters": len(set(clusters)),
            "n_marker_groups": len(markers),
            "n_da_peaks": len(da_peaks),
            "resolution": resolution,
            "random_state": random_state,
            "n_features": n_features,
            "marker_pvalue": marker_pvalue,
            "snapatac2_version": snap.__version__,
            "motif_note": motif_note,
        },
    }


def benchmark_scatac_reference_snapatac2(params: dict) -> dict:
    tool = "SnapATAC2"
    input_path_str = params.get("input_path")
    if not input_path_str or not Path(input_path_str).exists():
        return {
            "status": "not_run",
            "error_type": "MissingDependency",
            "tool": tool,
            "reason": "input_path missing or unreadable (expected a .h5ad / "
                      ".h5mu peak matrix or fragments file matching the ARIA run).",
        }

    try:
        import snapatac2  # noqa: F401
    except ImportError:
        return {
            "status": "not_run",
            "error_type": "MissingDependency",
            "tool": tool,
            "reason": (
                "SnapATAC2 is not installed in this environment. The P3b reference "
                "driver runs in aria-bench-atac-env, whose install (SnapATAC2) + "
                "public datasets are infra-gated on Samael; ARIA does not fabricate "
                "reference outputs."),
        }

    try:
        outputs = _run_reference_pipeline(Path(input_path_str), params)
    except NotImplementedError as exc:
        return {
            "status": "not_run",
            "error_type": "NotImplemented",
            "tool": tool,
            "reason": str(exc),
        }

    result: dict[str, Any] = {
        "status": "success",
        "tool": tool,
        "clusters": outputs.get("clusters"),
        "da_peaks": outputs.get("da_peaks"),
        "motifs": outputs.get("motifs"),
        "notes": outputs.get("notes"),
    }
    output_json = params.get("output_json")
    if output_json:
        out = Path(output_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
        result["artifacts"] = {"output_json": str(out)}
    return result


if __name__ == "__main__":
    from aria.scripts._base import run_script
    run_script(benchmark_scatac_reference_snapatac2)
