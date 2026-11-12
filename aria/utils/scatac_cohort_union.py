"""Sparse scATAC cohort assembly using interval-overlap consensus peaks (ADR-043)."""

from __future__ import annotations

from bisect import bisect_right
from pathlib import Path
import re


def resolve_sample_ids(paths: list[str], sample_ids: list[str] | None = None) -> list[str]:
    """Resolve safe, deterministic, unique sample prefixes in input order."""
    explicit = sample_ids is not None and len(sample_ids) == len(paths)
    resolved, used = [], set()
    for index, path in enumerate(paths):
        if explicit:
            candidate = str(sample_ids[index]).strip()
            if not candidate or "#" in candidate:
                raise ValueError("sample_ids must be non-empty and must not contain '#'")
        else:
            candidate = Path(path).name
            for suffix in (".fragments.tsv.gz", ".fragments.tsv", ".tsv.gz", ".bed.gz",
                           ".h5ad", ".h5mu", ".tsv", ".bed"):
                if candidate.endswith(suffix):
                    candidate = candidate[:-len(suffix)]
                    break
            candidate = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate).strip("_.-") or "sample"
        base, number = candidate, 2
        while candidate in used:
            candidate = f"{base}_{number}"
            number += 1
        used.add(candidate)
        resolved.append(candidate)
    return resolved


def _peak_intervals(adata) -> list[tuple[str, int, int]]:
    from aria.utils.peak2gene_core import parse_peak_coord

    coordinates = []
    has_columns = {"chrom", "start", "end"} <= set(adata.var.columns)
    for index, name in enumerate(adata.var_names):
        coord = parse_peak_coord(str(name))
        if coord is None and has_columns:
            row = adata.var.iloc[index]
            try:
                start, end = int(row["start"]), int(row["end"])
                if float(row["start"]) != start or float(row["end"]) != end:
                    raise ValueError("interval coordinates must be integers")
                coord = (str(row["chrom"]), start, end)
            except (TypeError, ValueError, OverflowError):
                coord = None
        if coord is None or not coord[0] or coord[1] < 0 or coord[2] <= coord[1]:
            raise ValueError(f"feature {name!r} is not a valid genomic interval")
        coordinates.append(coord)
    if not coordinates:
        raise ValueError("input has no genomic peak features")
    return coordinates


def _raw_counts(adata):
    """Select finite, nonnegative integer counts; reject normalized-only inputs."""
    import numpy as np
    from scipy import sparse

    for matrix in (adata.layers.get("counts"), adata.X):
        if matrix is None:
            continue
        matrix = sparse.csr_matrix(matrix)
        values = matrix.data
        if (np.all(np.isfinite(values)) and np.all(values >= 0)
                and np.all(values == np.floor(values))
                and (not values.size or values.max() < np.iinfo(np.int64).max)):
            return matrix.astype(np.int64)
    raise ValueError("input requires raw, nonnegative integer counts in layers['counts'] or X")


def merge_scatac_matrices(
    paths: list[str], output_path: str, *, sample_ids: list[str] | None = None,
    sample_metadata: dict[str, dict] | None = None,
) -> dict:
    """Project each donor's counts onto a sparse genomic consensus and concatenate.

    Half-open intervals merge only on positive overlap. Adjacent intervals remain
    separate. Each source feature contributes once to its containing consensus.
    """
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from aria.utils.safe_h5ad import read_h5ad

    def error(error_type, details):
        return {"status": "error", "error_type": error_type, "details": str(details)}

    if not isinstance(paths, list) or len(paths) < 2:
        return error("InvalidInputs", "at least two input matrices are required")
    if Path(output_path).suffix.lower() != ".h5ad":
        return error("InvalidOutputPath", "cohort output must have a .h5ad suffix")
    if Path(output_path).resolve() in {Path(path).resolve() for path in paths}:
        return error("InvalidOutputPath", "output must differ from every input")
    try:
        samples = resolve_sample_ids(paths, sample_ids)
        metadata = sample_metadata or {}
        datasets, all_intervals = [], []
        for path in paths:
            if not Path(path).is_file():
                return error("MissingInput", path)
            suffix = Path(path).suffix.lower()
            if suffix == ".h5mu":
                import mudata
                mdata = mudata.read_h5mu(path)
                if "atac" not in mdata.mod:
                    return error("MissingATACModality", path)
                adata = mdata.mod["atac"]
            elif suffix == ".h5ad":
                adata = read_h5ad(path)
            else:
                return error("UnsupportedInput", path)
            try:
                intervals = _peak_intervals(adata)
            except ValueError as exc:
                return error("NonGenomicPeakFeatures", f"{path}: {exc}")
            try:
                counts = _raw_counts(adata)
            except ValueError as exc:
                return error("NonIntegerCounts", f"{path}: {exc}")
            datasets.append((adata, counts, intervals))
            all_intervals.extend(intervals)

        consensus = []
        for chrom, start, end in sorted(all_intervals):
            if consensus and consensus[-1][0] == chrom and start < consensus[-1][2]:
                previous = consensus[-1]
                consensus[-1] = (chrom, previous[1], max(previous[2], end))
            else:
                consensus.append((chrom, start, end))
        chrom_index = {}
        for index, (chrom, start, end) in enumerate(consensus):
            starts, targets = chrom_index.setdefault(chrom, ([], []))
            starts.append(start)
            targets.append(index)
        var = pd.DataFrame(consensus, columns=["chrom", "start", "end"],
                           index=[f"{chrom}:{start}-{end}" for chrom, start, end in consensus])
        projected = []
        for sample, (adata, counts, intervals) in zip(samples, datasets):
            columns = []
            for chrom, start, end in intervals:
                starts, targets = chrom_index[chrom]
                target = targets[bisect_right(starts, start) - 1]
                if not (consensus[target][1] <= start < end <= consensus[target][2]):
                    raise ValueError("source interval is not contained in its consensus")
                columns.append(target)
            projection = sparse.csr_matrix(
                (np.ones(len(intervals), dtype=np.int64), (np.arange(len(intervals)), columns)),
                shape=(len(intervals), len(consensus)),
            )
            obs = adata.obs.copy()
            prefix = f"{sample}#"
            obs.index = pd.Index([str(cell) if str(cell).startswith(prefix) else prefix + str(cell)
                                  for cell in obs.index])
            if not obs.index.is_unique:
                return error("DuplicateCellBarcodes", sample)
            meta = dict(metadata.get(sample) or {})
            fallback_donor = meta.get("donor_id") or meta.get("replicate") or sample
            if "donor_id" in obs:
                obs["donor_id"] = obs["donor_id"].astype(object).where(
                    obs["donor_id"].notna() & (obs["donor_id"].astype(str) != ""), fallback_donor)
            else:
                obs["donor_id"] = fallback_donor
            for key, value in meta.items():
                if key not in {"sample_id", "donor_id"}:
                    obs[key] = value
            obs["sample_id"] = sample
            projected.append(ad.AnnData(X=(counts @ projection).tocsr(), obs=obs, var=var.copy()))
        merged = ad.concat(projected, join="outer", merge="same")
        merged.var = var.copy()
        merged.layers["counts"] = merged.X.copy()
        merged.uns["peak_provenance"] = "cohort_interval_union_consensus"
        merged.uns["samples"] = samples
        destination = Path(output_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        merged.write_h5ad(destination)
        return {"status": "success", "output_path": str(destination), "n_inputs": len(paths),
                "n_cells": int(merged.n_obs), "n_consensus_peaks": int(merged.n_vars),
                "peak_provenance": "cohort_interval_union_consensus", "samples": samples}
    except Exception as exc:
        return error(type(exc).__name__, exc)
