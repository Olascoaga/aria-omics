"""
ARIA scATAC fragments -> cell x peak matrix bridge
---------------------------------------------------
Builds the validated cell x peak AnnData that the scATAC matrix pipeline
(chromatin_lsi_clustering -> differential accessibility -> motif/regulatory)
requires, starting from a raw single-cell fragments file (the canonical scATAC
artifact: a barcoded, position-sorted `fragments.tsv.gz`).

This closes the C4 gap: before this, raw scATAC fragments only supported QC +
MACS3 peak calling + motif enrichment; LSI/DA/regulatory needed a pre-made
`.h5mu`. With this bridge, a fragments file (e.g. one downloaded by the GEO
connector, or produced by chromap from FASTQ) reaches the full pipeline.

Executed inside aria-chromatin-env (snapatac2) via EnvironmentManager.

Two peak modes:
  - provided  : a peak BED/narrowPeak is given -> make_peak_matrix directly
                (deterministic, no preliminary clustering).
  - de_novo   : no peaks given -> import -> tile matrix -> spectral+leiden
                (preliminary, only to group cells for peak calling) -> snapatac2
                MACS3 per group -> merge_peaks -> make_peak_matrix.

Input params:
  fragments_file: str   — barcoded fragments.tsv(.gz)
  genome:         str   — assembly (hg38/mm10/...); resolves chrom sizes + (de
                          novo) the MACS3 effective genome size
  peak_file:      str   — (optional) BED/narrowPeak peak set -> "provided" mode
  library_manifest: list — optional E5 library/sample/donor metadata keyed by
                           barcode prefixes in a merged fragments file
  output_dir:     str
  min_fragments:  int   (default 1000) — per-cell fragment floor
  resolution:     float (default 1.0)  — preliminary leiden resolution (de novo)
  macs3_qvalue:   float (default 0.05) — de novo peak-calling q-value

Output:
  {
    "status": "success",
    "output_path": str,        — cell x peak .h5ad for the LSI pipeline
    "n_cells": int,
    "n_peaks": int,
    "peak_mode": "provided" | "de_novo",
    "peak_source": str,
    "validation_level": "beta",
    "warnings": [str]
  }
  On a missing dependency / unknown assembly / no usable fragments it returns a
  structured not-run result (never a fabricated matrix).
"""

from __future__ import annotations
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def chromatin_fragments_to_matrix(params: dict) -> dict:
    fragments_file = params.get("fragments_file")
    genome = params.get("genome", "hg38")
    peak_file = params.get("peak_file")
    library_manifest = params.get("library_manifest") or []
    output_dir = Path(params.get("output_dir")
                      or (Path(fragments_file).parent / "aria_chromatin"
                          if fragments_file else "aria_chromatin"))
    min_fragments = int(params.get("min_fragments", 1000))
    resolution = float(params.get("resolution", 1.0))
    macs3_qvalue = float(params.get("macs3_qvalue", 0.05))
    warnings: list[str] = []

    def not_run(reason: str, **extra) -> dict:
        return {"status": "skipped", "ran": False,
                "analysis": "fragments_to_peak_matrix",
                "validation_level": "beta", "reason": reason,
                "warnings": warnings, **extra}

    if not fragments_file or not Path(fragments_file).exists():
        return not_run("fragments_file_missing",
                       message=f"fragments file not found: {fragments_file}")

    # Resolve the assembly -> snapatac2 genome object (chrom sizes + annotation),
    # reusing ARIA's governed mapping (same as QC/TSSe).
    from aria.utils import genomes
    attr = genomes.snapatac2_attr(genome)
    if not attr:
        return not_run("unknown_assembly",
                       message=(f"no snapatac2 genome for assembly "
                                f"'{genome or 'unknown'}'"))

    from aria.utils import privacy
    if not privacy.egress_allowed():
        return not_run("air_gapped_genome",
                       message=(f"snapatac2 genome '{attr}' is a governed "
                                f"auto-fetch and egress is disabled"))

    try:
        import snapatac2 as snap
    except ImportError as e:
        return not_run("snapatac2_unavailable", message=str(e))

    gobj = getattr(snap.genome, attr, None)
    if gobj is None:
        return not_run("snapatac2_genome_missing",
                       message=f"snapatac2 has no genome '{attr}'")

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "fragments_peak_matrix.h5ad"

    # Resume-by-file-validity: a prior valid matrix with matching peak intent.
    required_obs = ({"library_id", "sample_id", "donor_id"}
                    if library_manifest else set())
    if _matrix_valid(out_path, required_obs=required_obs):
        warnings.append("[resume] reused existing fragments_peak_matrix.h5ad")
        meta = _matrix_meta(out_path)
        return {"status": "success", "output_path": str(out_path),
                "validation_level": "beta", "resumed": True,
                "peak_mode": "provided" if peak_file else "de_novo",
                "warnings": warnings, **meta}

    import_fn = (getattr(snap.pp, "import_fragments", None)
                 or getattr(snap.pp, "import_data", None))
    if import_fn is None:
        return not_run("snapatac2_no_importer",
                       message="snapatac2 exposes no import_fragments/import_data")

    try:
        adata = import_fn(fragments_file, chrom_sizes=gobj,
                          min_num_fragments=min_fragments,
                          sorted_by_barcode=False, file=None)
    except Exception as e:
        return not_run("fragment_import_failed", message=str(e))

    n_cells = int(adata.n_obs)
    if n_cells == 0:
        return not_run("no_cells_pass_fragment_floor",
                       message=(f"no barcodes had >= {min_fragments} fragments"))

    n_library_annotated = _annotate_library_metadata(
        adata, library_manifest
    ) if library_manifest else 0
    if library_manifest and n_library_annotated != n_cells:
        return not_run(
            "library_metadata_correspondence_failed",
            message=(
                f"typed library manifest annotated {n_library_annotated}/{n_cells} "
                "fragment barcodes"
            ),
        )

    peak_mode = "provided" if peak_file and Path(peak_file).exists() else "de_novo"
    try:
        if peak_mode == "provided":
            peak_matrix = snap.pp.make_peak_matrix(
                adata, peak_file=str(peak_file), inplace=False)
            peak_source = str(peak_file)
        else:
            peak_source, peak_matrix = _de_novo_peaks(
                snap, adata, gobj, resolution, macs3_qvalue, warnings)
            if peak_matrix is None:
                return not_run("de_novo_peak_calling_failed",
                               message=peak_source)
    except Exception as e:
        return not_run("peak_matrix_failed", message=f"{peak_mode}: {e}")

    # Carry the per-cell QC the importer computed (n_fragment, tsse if present)
    # so the matrix pipeline keeps cell-level provenance.
    try:
        peak_matrix.write_h5ad(out_path)
    except Exception as e:
        return not_run("write_failed", message=str(e))

    n_peaks = int(peak_matrix.n_vars)
    return {
        "status": "success",
        "output_path": str(out_path),
        "n_cells": int(peak_matrix.n_obs),
        "n_peaks": n_peaks,
        "peak_mode": peak_mode,
        "peak_source": peak_source,
        "min_fragments": min_fragments,
        "genome": genome,
        "validation_level": "beta",
        "n_libraries": len(library_manifest) if library_manifest else None,
        "n_library_annotated": n_library_annotated,
        "warnings": warnings,
    }


def _de_novo_peaks(snap, adata, gobj, resolution, qvalue, warnings):
    """Preliminary clustering only to group cells for MACS3, then consensus peaks.

    The clustering here is NOT the reported clustering — it exists solely so
    MACS3 can call peaks per cell group. The returned cell x peak matrix is then
    re-clustered by the validated LSI pipeline downstream.
    """
    snap.pp.add_tile_matrix(adata)
    snap.pp.select_features(adata)
    snap.tl.spectral(adata)
    snap.pp.knn(adata)
    snap.tl.leiden(adata, resolution=resolution)
    n_groups = int(adata.obs["leiden"].nunique())
    warnings.append(f"de novo: {n_groups} preliminary cell groups for MACS3 "
                    f"(re-clustered downstream by the LSI pipeline)")
    snap.tl.macs3(adata, groupby="leiden", qvalue=qvalue)
    merged = snap.tl.merge_peaks(adata.uns["macs3"], gobj)
    peaks = merged["Peaks"] if "Peaks" in merged else merged
    peak_matrix = snap.pp.make_peak_matrix(adata, use_rep=peaks, inplace=False)
    return f"de_novo_macs3 (q<{qvalue}, {n_groups} groups)", peak_matrix


def _annotate_library_metadata(adata, libraries: list[dict]) -> int:
    """Map prefixed fragment barcodes back to typed library metadata."""
    if not libraries:
        return 0
    rows = {}
    for library in libraries:
        library_id = str(library.get("library_id") or "")
        prefix = str(library.get("barcode_prefix") or f"{library_id}#")
        metadata = dict(library.get("metadata") or {})
        metadata.update({
            "library_id": library_id,
            "sample_id": str(library.get("sample_id") or ""),
            "donor_id": str(
                library.get("donor_id") or metadata.get("donor_id") or ""
            ),
        })
        rows[prefix] = metadata

    annotations = []
    matched = 0
    for barcode in map(str, adata.obs.index):
        record = next(
            (metadata for prefix, metadata in rows.items()
             if barcode.startswith(prefix)),
            None,
        )
        annotations.append(record or {})
        matched += int(record is not None)

    columns = sorted({key for row in annotations for key in row})
    for column in columns:
        adata.obs[column] = [row.get(column) for row in annotations]
    return matched


def _matrix_valid(path: Path, required_obs: set[str] | None = None) -> bool:
    try:
        if not path.exists() or path.stat().st_size < 1024:
            return False
        import anndata as ad
        a = ad.read_h5ad(path, backed="r")
        ok = a.n_obs > 0 and a.n_vars > 0
        if required_obs:
            ok = ok and required_obs <= set(a.obs.columns)
        a.file.close()
        return ok
    except Exception:
        return False


def _matrix_meta(path: Path) -> dict:
    try:
        import anndata as ad
        a = ad.read_h5ad(path, backed="r")
        meta = {"n_cells": int(a.n_obs), "n_peaks": int(a.n_vars)}
        a.file.close()
        return meta
    except Exception:
        return {}


if __name__ == "__main__":
    run_script(chromatin_fragments_to_matrix)
