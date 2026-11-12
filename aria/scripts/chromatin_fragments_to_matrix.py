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
    "validation_level": "production",
    "warnings": [str]
  }
  On a missing dependency / unknown assembly / no usable fragments it returns a
  structured not-run result (never a fabricated matrix).
"""

from __future__ import annotations
import gzip
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def _resolve_local_chrom_sizes(params: dict) -> tuple[dict[str, int] | None, str | None, list[str]]:
    """Resolve only explicitly supplied local assembly assets, with no network IO."""
    warnings = []

    def positive_sizes(items):
        sizes = {}
        for chrom, length in items:
            value = int(length)
            if isinstance(length, bool) or value <= 0 or float(length) != value or not str(chrom).strip():
                raise ValueError("contig names must be non-empty and lengths must be positive integers")
            if str(chrom) in sizes:
                raise ValueError(f"duplicate contig: {chrom}")
            sizes[str(chrom)] = value
        if not sizes:
            raise ValueError("chromosome sizes are empty")
        return sizes

    def parse_table(path):
        with open(path) as handle:
            return positive_sizes(line.split()[:2] for line in handle
                                  if line.strip() and not line.startswith("#"))

    sizes = params.get("chrom_sizes")
    if isinstance(sizes, dict):
        try:
            return positive_sizes(sizes.items()), "chrom_sizes", warnings
        except (ValueError, TypeError, OverflowError) as exc:
            warnings.append(f"invalid chrom_sizes: {exc}")
    for key in ("chrom_sizes", "genome_fai"):
        path = params.get(key)
        if isinstance(path, (str, Path)) and path:
            try:
                return parse_table(path), str(path), warnings
            except (OSError, ValueError, TypeError, OverflowError) as exc:
                warnings.append(f"cannot read {key}: {exc}")
    fasta = params.get("genome_fasta")
    if fasta:
        fasta = Path(fasta)
        if Path(str(fasta) + ".fai").is_file():
            try:
                return parse_table(str(fasta) + ".fai"), str(fasta) + ".fai", warnings
            except (OSError, ValueError, TypeError, OverflowError) as exc:
                warnings.append(f"cannot read FASTA index: {exc}")
        try:
            opener = gzip.open if fasta.suffix == ".gz" else open
            lengths, contig, length = [], None, 0
            with opener(fasta, "rt") as handle:
                for line in handle:
                    if line.startswith(">"):
                        if contig is not None:
                            lengths.append((contig, length))
                        header = line[1:].split()
                        if not header:
                            raise ValueError("FASTA header has no contig name")
                        contig, length = header[0], 0
                    elif line.strip():
                        if contig is None:
                            raise ValueError("FASTA sequence has no header")
                        length += len("".join(line.split()))
            if contig is not None:
                lengths.append((contig, length))
            return positive_sizes(lengths), str(fasta), warnings
        except (OSError, ValueError, TypeError, OverflowError) as exc:
            warnings.append(f"cannot read genome_fasta: {exc}")
    return None, None, warnings


def chromatin_fragments_to_matrix(params: dict) -> dict:
    from aria.utils.stage_manifest import build_stage_manifest, stage_is_current, write_stage_manifest

    fragments_file = params.get("fragments_file")
    fragments_files = params.get("fragments_files")
    if fragments_files is None:
        fragments_files = [fragments_file] if fragments_file else []
    genome = params.get("genome", "hg38")
    peak_file = params.get("peak_file")
    library_manifest = params.get("library_manifest") or []
    if isinstance(library_manifest, dict):
        library_manifest = library_manifest.get("libraries") or []
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
                "validation_level": "production", "reason": reason,
                "warnings": warnings, **extra}

    if (not isinstance(fragments_files, list) or not fragments_files
            or any(not isinstance(path, (str, Path)) or not Path(path).is_file()
                   for path in fragments_files)):
        return not_run("fragments_file_missing", message="one or more fragment files are missing")
    fragments_files = [str(path) for path in fragments_files]
    fragments_file = fragments_files[0]
    if peak_file and not Path(peak_file).is_file():
        return not_run("peak_file_missing", message=f"peak file not found: {peak_file}")

    gobj, local_source, local_warnings = _resolve_local_chrom_sizes(params)
    warnings.extend(local_warnings)
    attr = None
    if gobj is None:
        from aria.utils import genomes, privacy
        attr = genomes.snapatac2_attr(genome)
        if not attr:
            return not_run("unknown_assembly", message=f"no snapatac2 genome for assembly '{genome or 'unknown'}'")
        if not privacy.egress_allowed():
            return not_run("air_gapped_genome",
                           message=f"snapatac2 genome '{attr}' is a governed auto-fetch and egress is disabled")
    try:
        import snapatac2 as snap
    except ImportError as exc:
        return not_run("snapatac2_unavailable", message=str(exc))
    if gobj is None:
        gobj = getattr(snap.genome, attr, None)
        if gobj is None:
            return not_run("snapatac2_genome_missing", message=f"snapatac2 has no genome '{attr}'")

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "fragments_peak_matrix.h5ad"
    samples = []
    if len(fragments_files) >= 2:
        from aria.utils.scatac_cohort_union import resolve_sample_ids
        try:
            samples = resolve_sample_ids(fragments_files, params.get("sample_ids"))
        except ValueError as exc:
            return not_run("invalid_sample_ids", message=str(exc))
        if "library_manifest" not in params:
            library_manifest = []
            for sample in samples:
                meta = dict((params.get("sample_metadata") or {}).get(sample) or {})
                library_manifest.append({
                    "library_id": sample, "sample_id": sample,
                    "donor_id": meta.get("donor_id") or meta.get("replicate") or sample,
                    "barcode_prefix": f"{sample}#", "metadata": meta,
                })

    manifest_path = output_dir / "chromatin_fragments_to_matrix.stage_manifest.json"
    inputs = [(f"fragments_{index}", path) for index, path in enumerate(fragments_files)]
    for key in ("peak_file", "genome_fai", "genome_fasta", "chrom_sizes"):
        value = params.get(key)
        if isinstance(value, (str, Path)) and value and Path(value).is_file():
            inputs.append((key, str(value)))
    if local_source and local_source != "chrom_sizes" and Path(local_source).is_file():
        inputs.append(("resolved_chrom_sizes", local_source))
    stage_params = {"genome": genome, "min_fragments": min_fragments,
                    "resolution": resolution, "macs3_qvalue": macs3_qvalue,
                    "library_manifest": library_manifest, "sample_metadata": params.get("sample_metadata"),
                    "sample_ids": samples, "chrom_sizes": gobj if isinstance(gobj, dict) else None}
    tool_version = getattr(snap, "__version__", "")
    required_obs = ({"library_id", "sample_id", "donor_id"} if library_manifest else set())
    if (_matrix_valid(out_path, required_obs=required_obs) and stage_is_current(
            manifest_path, inputs=inputs, params=stage_params, tool_version=tool_version)[0]):
        warnings.append("[resume] reused existing fragments_peak_matrix.h5ad")
        return {"status": "success", "output_path": str(out_path),
                "validation_level": "production", "resumed": True,
                "peak_mode": "provided" if peak_file else "de_novo",
                "warnings": warnings, **_matrix_meta(out_path)}

    if len(fragments_files) >= 2:
        from aria.scripts.chromatin_scatac_align import _merge_library_fragments, _bgzip_tabix
        combined = output_dir / "combined_cohort.fragments.tsv"
        records = [{"library_id": sample, "fragments_file": path}
                   for sample, path in zip(samples, fragments_files)]
        union = _merge_library_fragments(records, combined)
        if union.get("status") != "success":
            return not_run("fragment_union_failed", message=union.get("message"))
        compressed = output_dir / "combined_cohort.fragments.tsv.gz"
        result = _bgzip_tabix(combined, compressed, warnings)
        fragments_file = str(compressed if result == "ok" else combined)
        if result != "ok":
            warnings.append(f"bgzip/tabix unavailable ({result}); leaving plain cohort TSV")

    import_fn = (getattr(snap.pp, "import_fragments", None) or getattr(snap.pp, "import_data", None))
    if import_fn is None:
        return not_run("snapatac2_no_importer", message="snapatac2 exposes no import_fragments/import_data")
    try:
        adata = import_fn(fragments_file, chrom_sizes=gobj, min_num_fragments=min_fragments,
                          sorted_by_barcode=False, file=None)
    except Exception as exc:
        return not_run("fragment_import_failed", message=str(exc))
    n_cells = int(adata.n_obs)
    if n_cells == 0:
        return not_run("no_cells_pass_fragment_floor", message=f"no barcodes had >= {min_fragments} fragments")
    n_library_annotated = _annotate_library_metadata(adata, library_manifest) if library_manifest else 0
    if library_manifest and n_library_annotated != n_cells:
        return not_run("library_metadata_correspondence_failed",
                       message=f"typed library manifest annotated {n_library_annotated}/{n_cells} fragment barcodes")

    peak_mode = "provided" if peak_file else "de_novo"
    try:
        if peak_mode == "provided":
            peak_matrix = snap.pp.make_peak_matrix(adata, peak_file=str(peak_file), inplace=False)
            peak_source = str(peak_file)
        else:
            peak_source, peak_matrix = _de_novo_peaks(snap, adata, gobj, resolution, macs3_qvalue, warnings)
            if peak_matrix is None:
                return not_run("de_novo_peak_calling_failed", message=peak_source)
        # Some snapatac2 versions return only obs_names on the peak matrix.
        # Align by barcode so metadata and importer QC survive reordering.
        if not peak_matrix.obs_names.isin(adata.obs_names).all():
            return not_run("peak_matrix_failed", message="peak matrix contains unknown fragment barcodes")
        for column in adata.obs.columns:
            if column not in peak_matrix.obs.columns:
                peak_matrix.obs[column] = adata.obs[column].reindex(peak_matrix.obs_names)
    except Exception as exc:
        return not_run("peak_matrix_failed", message=f"{peak_mode}: {exc}")
    try:
        peak_matrix.write_h5ad(out_path)
        write_stage_manifest(manifest_path, build_stage_manifest(
            stage="chromatin_fragments_to_matrix", inputs=inputs, params=stage_params, tool_version=tool_version))
    except Exception as exc:
        return not_run("write_failed", message=str(exc))
    return {"status": "success", "output_path": str(out_path),
            "n_cells": int(peak_matrix.n_obs), "n_peaks": int(peak_matrix.n_vars),
            "peak_mode": peak_mode, "peak_source": peak_source, "min_fragments": min_fragments,
            "genome": genome, "validation_level": "production",
            "n_libraries": len(library_manifest) if library_manifest else None,
            "n_library_annotated": n_library_annotated, "warnings": warnings}


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
