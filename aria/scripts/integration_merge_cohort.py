"""Sparse paired multi-donor cohort assembly with ADR-043 interval union semantics."""

from __future__ import annotations

import os
import sys
from bisect import bisect_right
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script
from aria.utils.scatac_cohort_union import resolve_sample_ids, _peak_intervals
from aria.utils.mudata_io import (
    MultiomeInputError,
    read_h5mu_paired,
    paired_input_files,
    _write_paired_atomic,
    _stage_cache,
    _stage_complete,
)


def integration_merge_cohort(params: dict[str, Any]) -> dict[str, Any]:
    try:
        import anndata as ad
        import mudata
        import numpy as np
        import pandas as pd
        from scipy import sparse

        def fail(kind, details):
            raise MultiomeInputError(kind, str(details))

        paths = params.get("inputs")
        if not isinstance(paths, list) or len(paths) < 2:
            fail("InvalidInputs", "At least two paired inputs are required")
        if params.get("sample_ids") is not None and len(params["sample_ids"]) != len(
            paths
        ):
            fail("InvalidSampleIds", "sample_ids must match inputs in length")
        samples = resolve_sample_ids(paths, params.get("sample_ids"))
        metadata = params.get("sample_metadata") or {}
        join = params.get("rna_gene_join", "inner")
        if join not in {"inner", "outer"}:
            fail("InvalidParameter", "rna_gene_join must be inner or outer")
        min_fraction = float(params.get("min_shared_barcode_fraction", 0.50))
        if not 0 <= min_fraction <= 1:
            fail("InvalidParameter", "min_shared_barcode_fraction must be in [0, 1]")
        directory = Path(params.get("output_dir") or "/tmp/aria_cohort")
        output = str(params.get("output_path") or directory / "cohort_merged.h5mu")
        if Path(output).resolve() in {Path(p).resolve() for p in paths}:
            fail("InvalidOutputPath", "Output must differ from every input")
        inputs = [
            entry
            for i, path in enumerate(paths)
            for entry in paired_input_files(f"sample:{i}", path)
        ]
        normalized = dict(
            params,
            sample_ids=samples,
            rna_gene_join=join,
            min_shared_barcode_fraction=min_fraction,
        )
        cached, state = _stage_cache(
            "integration_merge_cohort",
            directory,
            inputs,
            normalized,
            params.get("resume", True),
        )
        if cached:
            return cached
        datasets, all_intervals, per_sample = [], [], {}
        for sample, path in zip(samples, paths):
            paired = read_h5mu_paired(path, require_raw_counts=True)
            if paired["status"] != "success":
                return paired
            if (
                min(paired["shared_fraction_rna"], paired["shared_fraction_atac"])
                < min_fraction
            ):
                fail(
                    "LowBarcodeOverlap",
                    f"{sample}: shared barcode fraction below {min_fraction}",
                )
            shared = paired["shared_barcodes"]
            rna, atac = paired["rna_adata"], paired["atac_adata"]
            if not rna.var_names.is_unique:
                fail("DuplicateGeneNames", sample)
            rrows, arows = rna.obs_names.get_indexer(
                shared
            ), atac.obs_names.get_indexer(shared)
            rcounts = sparse.csr_matrix(paired["rna_counts_matrix"])[rrows]
            acounts = sparse.csr_matrix(paired["atac_counts_matrix"])[arows]
            if any(
                m.data.size and m.data.max() >= np.iinfo(np.int64).max
                for m in (rcounts, acounts)
            ):
                fail("NonIntegerCounts", "Counts exceed int64 storage")
            rcounts, acounts = rcounts.astype(np.int64), acounts.astype(np.int64)
            try:
                intervals = _peak_intervals(atac)
            except ValueError as exc:
                fail("NonGenomicPeakFeatures", exc)
            # Preserve global and modality metadata; reject conflicting donor identities.
            obs = paired["mdata"].obs.reindex(shared).copy()
            for source in (rna.obs.reindex(shared), atac.obs.reindex(shared)):
                for col in source:
                    if col == "donor_id" and col in obs:
                        conflict = (
                            obs[col].notna()
                            & source[col].notna()
                            & (obs[col].astype(str) != source[col].astype(str))
                        )
                        if conflict.any():
                            fail("ConflictingDonorMetadata", sample)
                    obs[col] = (
                        obs[col].combine_first(source[col])
                        if col in obs
                        else source[col]
                    )
            meta = dict(metadata.get(sample) or metadata.get(Path(path).name) or {})
            fallback = meta.get("donor_id") or sample
            if "donor_id" in obs:
                obs["donor_id"] = (
                    obs["donor_id"]
                    .astype(object)
                    .where(
                        obs["donor_id"].notna() & (obs["donor_id"].astype(str) != ""),
                        fallback,
                    )
                )
            else:
                obs["donor_id"] = fallback
            for key, value in meta.items():
                if key not in {"donor_id", "sample_id"}:
                    obs[key] = value
            obs["sample_id"] = sample
            prefix = sample + "#"
            obs.index = pd.Index(
                [b if b.startswith(prefix) else prefix + b for b in shared]
            )
            if not obs.index.is_unique:
                fail("DuplicateBarcodes", sample)
            datasets.append((rna.var.copy(), rcounts, acounts, intervals, obs))
            all_intervals.extend(intervals)
            per_sample[sample] = len(shared)
        genes = list(datasets[0][0].index)
        if join == "inner":
            common = set(genes)
            for var, *_ in datasets[1:]:
                common.intersection_update(var.index)
            genes = [g for g in genes if g in common]
            if len(genes) < 3:
                fail(
                    "EmptyGeneIntersection",
                    f"Only {len(genes)} shared RNA genes; at least 3 required",
                )
        else:
            seen = set(genes)
            for var, *_ in datasets[1:]:
                for gene in var.index:
                    if gene not in seen:
                        genes.append(gene)
                        seen.add(gene)
        # Keep annotations where all available values agree, including union-only genes.
        gene_var = pd.DataFrame(index=genes)
        columns = list(
            dict.fromkeys(col for var, *_ in datasets for col in var.columns)
        )
        for col in columns:
            values, consistent = {}, True
            for var, *_ in datasets:
                if col not in var:
                    continue
                for gene, value in var[col].items():
                    if gene not in gene_var.index or pd.isna(value):
                        continue
                    if gene in values and values[gene] != value:
                        consistent = False
                        break
                    values[gene] = value
                if not consistent:
                    break
            if consistent:
                gene_var[col] = pd.Series(values).reindex(genes)
        consensus = []
        for chrom, start, end in sorted(all_intervals):
            if consensus and consensus[-1][0] == chrom and start < consensus[-1][2]:
                previous = consensus[-1]
                consensus[-1] = (chrom, previous[1], max(previous[2], end))
            else:
                consensus.append((chrom, start, end))
        chrom_index = {}
        for i, (chrom, start, _) in enumerate(consensus):
            starts, targets = chrom_index.setdefault(chrom, ([], []))
            starts.append(start)
            targets.append(i)
        peak_var = pd.DataFrame(
            consensus,
            columns=["chrom", "start", "end"],
            index=[f"{c}:{s}-{e}" for c, s, e in consensus],
        )
        gene_index = {g: i for i, g in enumerate(genes)}
        rparts, aparts, observations = [], [], []
        for var, rcounts, acounts, intervals, obs in datasets:
            source, target = [], []
            for i, g in enumerate(var.index):
                if g in gene_index:
                    source.append(i)
                    target.append(gene_index[g])
            gprojection = sparse.csr_matrix(
                (np.ones(len(source), dtype=np.int64), (source, target)),
                shape=(len(var), len(genes)),
            )
            columns = []
            for chrom, start, end in intervals:
                starts, targets = chrom_index[chrom]
                target = targets[bisect_right(starts, start) - 1]
                if not consensus[target][1] <= start < end <= consensus[target][2]:
                    fail(
                        "InvalidPeakProjection",
                        "Source peak is not contained in consensus",
                    )
                columns.append(target)
            projection = sparse.csr_matrix(
                (
                    np.ones(len(intervals), dtype=np.int64),
                    (np.arange(len(intervals)), columns),
                ),
                shape=(len(intervals), len(consensus)),
            )
            rparts.append((rcounts @ gprojection).tocsr())
            aparts.append((acounts @ projection).tocsr())
            observations.append(obs)
        obs = pd.concat(observations)
        if not obs.index.is_unique:
            fail("DuplicateBarcodes", "Sample prefixes collide across inputs")
        rna = ad.AnnData(
            sparse.vstack(rparts, format="csr"), obs=obs.copy(), var=gene_var
        )
        atac = ad.AnnData(
            sparse.vstack(aparts, format="csr"), obs=obs.copy(), var=peak_var
        )
        for mod in (rna, atac):
            mod.layers["counts"] = mod.X.copy()
        out = mudata.MuData({"rna": rna, "atac": atac})
        out.obs = obs.copy()
        out.uns["aria_cohort_merge"] = {
            "inputs": list(map(str, paths)),
            "sample_ids": samples,
            "rna_gene_join": join,
            "peak_provenance": "cohort_interval_union_consensus",
            "per_sample_cells": per_sample,
        }
        _write_paired_atomic(out, output)
        return _stage_complete(
            {
                "status": "success",
                "output_path": output,
                "n_samples": len(samples),
                "sample_ids": samples,
                "n_cells_total": len(obs),
                "per_sample_cells": per_sample,
                "n_genes": len(genes),
                "rna_gene_join": join,
                "n_consensus_peaks": len(consensus),
                "validation_level": "production",
                "warnings": [],
            },
            state,
        )
    except Exception as exc:
        return {
            "status": "error",
            "error_type": getattr(
                exc,
                "error_type",
                (
                    "MissingDependency"
                    if isinstance(exc, ImportError)
                    else type(exc).__name__
                ),
            ),
            "details": str(exc),
        }


if __name__ == "__main__":
    run_script(integration_merge_cohort)
