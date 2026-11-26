"""Bulk RNA-seq count/metadata loading, correspondence enforcement, metadata
inference, technical-replicate aggregation and comparison resolution (A7 split).

Re-exported from aria.scripts.rna_bulk_de."""
from __future__ import annotations
import os
import re
import warnings

from pathlib import Path

from aria.utils.count_classifier import classify_matrix, validate_raw_count_matrix
from aria.scripts.rna_bulk.deseq2 import _run_deseq2
from aria.scripts.rna_bulk.gtf_io import _strip_gene_version

VALID_COUNT_UNION_POLICIES = ("error", "intersect", "union_zero_fill")


def _resolve_unique_sample_names(paths: list[str]) -> tuple[list[str] | None, str | None]:
    """Derive unique per-sample column names, disambiguating colliding stems by parent dir."""
    resolved_paths = [Path(p) for p in paths]
    stems = [p.stem for p in resolved_paths]
    if len(set(stems)) == len(stems):
        return stems, None

    stem_counts: dict[str, int] = {}
    for s in stems:
        stem_counts[s] = stem_counts.get(s, 0) + 1

    names: list[str] = []
    for p, stem in zip(resolved_paths, stems, strict=True):
        if stem_counts[stem] == 1:
            names.append(stem)
        else:
            parent = p.parent.name.strip()
            names.append(f"{parent}_{stem}" if parent else stem)

    if len(set(names)) != len(names):
        # Try one additional ancestor directory for nested collisions (a/sub/counts.tsv vs b/sub/counts.tsv)
        name_counts: dict[str, int] = {}
        for n in names:
            name_counts[n] = name_counts.get(n, 0) + 1
        refined: list[str] = []
        for p, stem, cand in zip(resolved_paths, stems, names, strict=True):
            if name_counts[cand] == 1:
                refined.append(cand)
            else:
                grand = p.parent.parent.name.strip()
                parent = p.parent.name.strip()
                parts = [x for x in (grand, parent, stem) if x]
                refined.append("_".join(parts) if parts else cand)
        names = refined

    if len(set(names)) != len(names):
        dups = sorted({n for n in names if names.count(n) > 1})
        return None, (
            f"Multiple per-sample count files collide on sample name(s) {dups}. "
            "Provide distinct filenames or parent directories for each sample."
        )
    return names, None


def _load_counts(
    files: list,
    allow_nonraw: bool = False,
    count_union_policy: str = "error",
) -> tuple:
    """
    Load counts matrix from various formats.
    Returns (DataFrame, warnings_list, count_meta).
    genes × samples orientation enforced.

    Raw-count guard (audit 2026-05-29, B10 / P-RAWCLASS): DESeq2 requires raw
    integer counts. The loaded matrix is classified before rounding so that
    TPM/CPM/FPKM/log-normalized/scaled inputs are NOT silently coerced into
    pseudo-counts. Non-raw matrices are hard-refused unless ``allow_nonraw`` is
    set, in which case they are coerced at explicit low confidence. ``count_meta``
    always carries ``count_source`` / ``kind`` (or ``refused`` details).

    Multi-file merge policy (D4-001): Ensembl version suffixes are stripped
    consistently before comparing per-file gene sets. If gene sets still differ,
    ``count_union_policy`` governs the behavior:
      - ``"error"`` (default): fail closed with mismatch size and examples;
      - ``"intersect"``: retain only the shared gene intersection;
      - ``"union_zero_fill"``: take the union and fill missing genes with 0
        (only for per-sample files that legitimately omit zero-count genes).
    """
    import pandas as pd
    warnings = []

    policy = str(count_union_policy or "error").strip()
    if policy not in VALID_COUNT_UNION_POLICIES:
        err = (
            f"Invalid count_union_policy={count_union_policy!r}. Valid options "
            "for merging per-sample count files are: 'error' (default; fail "
            "closed if per-sample gene sets differ), 'intersect' (keep only "
            "genes present in every file), or 'union_zero_fill' (outer union "
            "filling missing genes with 0, appropriate only when per-sample "
            "files legitimately omit zero-count genes)."
        )
        return None, [err], {
            "refused": True,
            "error_type": "InvalidCountUnionPolicy",
            "details": err,
        }

    valid = [f for f in files if Path(f).exists()]
    if not valid:
        return None, ["No valid count files found."], {}

    # Detect format
    count_files = [f for f in valid
                   if any(f.endswith(x) for x in
                          [".tsv", ".csv", ".txt", ".counts",
                           ".featureCounts", ".htseq"])]
    if not count_files:
        # Try all valid files
        count_files = valid

    count_union_meta: dict | None = None

    # Try to load as single matrix or merge multiple per-sample files
    if len(count_files) == 1:
        sep = "," if str(count_files[0]).lower().endswith(".csv") else "\t"
        try:
            counts = pd.read_csv(count_files[0], sep=sep, index_col=0,
                                  comment="#")
            # Drop non-count columns (featureCounts format)
            drop_cols = [c for c in counts.columns
                         if c in ("Chr","Start","End","Strand","Length")]
            counts = counts.drop(columns=drop_cols, errors="ignore")
            # Keep only numeric columns
            counts = counts.select_dtypes(include="number")
            if counts.empty:
                return None, ["Count matrix has no numeric columns."], {}
        except Exception as e:
            return None, [f"Failed to load {count_files[0]}: {e}"], {}
    else:
        # Multiple per-sample files — merge by gene ID (fail closed if any sample fails)
        sample_names, name_err = _resolve_unique_sample_names(count_files)
        if name_err is not None or sample_names is None:
            return None, [name_err or "Duplicate per-sample count filenames."], {
                "refused": True,
                "error_type": "DuplicateSampleNames",
                "details": name_err or "Duplicate per-sample count filenames.",
            }

        frames = []
        per_file_sets: list[tuple[str, set[str]]] = []
        for f, sample_name in zip(count_files, sample_names, strict=True):
            sep = "," if str(f).lower().endswith(".csv") else "\t"
            try:
                df = pd.read_csv(f, sep=sep, index_col=0, comment="#",
                                  header=None)
                num_df = df.select_dtypes(include="number")
                if num_df.empty and len(df) >= 2:
                    # Standard per-sample files often carry a header row (e.g. gene_id\tcount).
                    df_hdr = pd.read_csv(
                        f, sep=sep, index_col=0, comment="#", header=0
                    )
                    drop_cols = [
                        c for c in df_hdr.columns
                        if c in ("Chr", "Start", "End", "Strand", "Length")
                    ]
                    df_hdr = df_hdr.drop(columns=drop_cols, errors="ignore")
                    num_df = df_hdr.select_dtypes(include="number")
                if num_df.empty:
                    return None, [f"Count file {f} has no numeric columns."], {}
                if num_df.shape[1] != 1:
                    err = (
                        f"Per-sample count file {f} contains {num_df.shape[1]} "
                        "numeric columns (expected 1). Pass a single combined "
                        "matrix file or one count column per per-sample file."
                    )
                    return None, [err], {
                        "refused": True,
                        "error_type": "InvalidPerSampleCountFile",
                        "details": err,
                    }
                num_df = num_df.copy()
                num_df.columns = [sample_name]
                num_df.index = pd.Index(
                    [_strip_gene_version(str(g)) for g in num_df.index]
                )
                if not num_df.index.is_unique:
                    dups = sorted(
                        set(num_df.index[num_df.index.duplicated()].tolist())
                    )
                    err = (
                        f"Count file {f} contains {len(dups)} duplicate gene "
                        f"ID(s) after stripping Ensembl version suffixes "
                        f"(e.g. {dups[:10]}). Resolve duplicate gene entries "
                        "before merging."
                    )
                    return None, [err], {
                        "refused": True,
                        "error_type": "DuplicateGeneIds",
                        "details": err,
                    }
                per_file_sets.append((sample_name, set(num_df.index.tolist())))
                frames.append(num_df)
            except (OSError, ValueError, UnicodeDecodeError,
                    pd.errors.ParserError, pd.errors.EmptyDataError) as e:
                return None, [f"Failed to load count file {f}: {e}"], {}

        union_genes = set().union(*(gset for _, gset in per_file_sets))
        inter_genes = set.intersection(*(gset for _, gset in per_file_sets))
        sym_diff = union_genes - inter_genes
        per_file_missing = {
            label: len(union_genes - gset) for label, gset in per_file_sets
        }
        count_union_meta = {
            "policy": policy,
            "n_genes_union": len(union_genes),
            "n_genes_intersection": len(inter_genes),
            "n_genes_mismatch": len(sym_diff),
            "per_file_missing": per_file_missing,
        }

        if sym_diff:
            examples = sorted(sym_diff)[:10]
            if policy == "error":
                err = (
                    f"Per-sample count files have mismatched gene sets after "
                    f"stripping Ensembl version suffixes: {len(sym_diff)} "
                    f"gene(s) differ across files (union={len(union_genes)}, "
                    f"intersection={len(inter_genes)}, examples={examples}, "
                    f"per_file_missing={per_file_missing}). A gene-set "
                    "mismatch is a data-integrity issue that must be resolved "
                    "explicitly (ARIA never silently intersects or zero-fills "
                    "missing genes). Pass count_union_policy explicitly: "
                    "'error' (default: fail closed on any gene-set mismatch), "
                    "'intersect' (keep only genes present in every file), or "
                    "'union_zero_fill' (outer union filling missing genes "
                    "with 0; use only if per-sample files legitimately omit "
                    "zero-count genes, such as some GEO supplementary files)."
                )
                return None, [err], {
                    "refused": True,
                    "error_type": "CountGeneSetMismatch",
                    "details": err,
                    "count_union": count_union_meta,
                }
            if policy == "intersect":
                if not inter_genes:
                    err = (
                        "Per-sample count files have 0 shared genes under "
                        f"count_union_policy='intersect' "
                        f"(union={len(union_genes)}, examples={examples})."
                    )
                    return None, [err], {
                        "refused": True,
                        "error_type": "EmptyGeneIntersection",
                        "details": err,
                        "count_union": count_union_meta,
                    }
                counts = pd.concat(frames, axis=1, join="inner")
                warnings.append(
                    f"Per-sample count files had {len(sym_diff)} mismatched "
                    f"gene(s) (union={len(union_genes)}, "
                    f"intersection={len(inter_genes)}); restricted to the "
                    "shared gene intersection under "
                    f"count_union_policy='intersect' "
                    f"(per_file_missing={per_file_missing})."
                )
            else:  # policy == "union_zero_fill"
                counts = pd.concat(frames, axis=1, join="outer").fillna(0)
                warnings.append(
                    f"Per-sample count files had {len(sym_diff)} mismatched "
                    f"gene(s) (union={len(union_genes)}, "
                    f"intersection={len(inter_genes)}); filled missing genes "
                    "with 0 under count_union_policy='union_zero_fill' "
                    f"(per_file_missing={per_file_missing})."
                )
        else:
            counts = pd.concat(frames, axis=1, join="inner")

    # Ensure genes × samples (more genes than samples in typical experiments)
    if counts.shape[1] > counts.shape[0]:
        warnings.append(
            f"Transposing matrix: detected {counts.shape[1]} rows × "
            f"{counts.shape[0]} cols — expected genes as rows."
        )
        counts = counts.T

    # Raw-count guard (B10 / P-RAWCLASS): classify BEFORE rounding so a
    # normalized matrix is not silently turned into pseudo-counts.
    source_hint = ";".join(str(f) for f in count_files)
    info = classify_matrix(
        counts.values,
        gene_ids=list(counts.index),
        source_hint=source_hint,
    )
    if info["is_raw_counts"]:
        count_source = "raw_counts"
    elif allow_nonraw:
        count_source = "coerced_nonraw"
        warnings.append(
            f"Count matrix does not look like raw counts "
            f"(kind={info['kind']}, score={info.get('raw_count_score', 0):.2f}, "
            f"max={info['max']:.2f}); coercing to "
            f"integers because allow_nonraw_counts=True. DESeq2 results are "
            f"LOW CONFIDENCE — supply raw counts for a valid negative-binomial "
            f"fit."
        )
    else:
        return None, warnings, {
            "refused":    True,
            "error_type": "NonRawCounts",
            "kind":       info["kind"],
            "raw_count_score": info.get("raw_count_score"),
            "confidence": info.get("confidence"),
            "sub_scores": info.get("sub_scores", {}),
            "score_basis": info.get("score_basis", {}),
            "details": (
                f"Count matrix does not look like raw counts "
                f"(kind={info['kind']}, score={info.get('raw_count_score', 0):.2f}, "
                f"max={info['max']:.2f}, "
                f"min={info['min']:.2f}). DESeq2 requires raw integer counts; "
                f"TPM/CPM/FPKM/log-normalized/scaled inputs are invalid. Supply "
                f"a raw-count matrix, or set allow_nonraw_counts=True to coerce "
                f"at low confidence."
            ),
        }

    # B6: full vectorized validation before rounding the WHOLE matrix. classify_matrix
    # decided from a <=200-row sample; a fractional/negative/NaN/inf value in an
    # unsampled row would otherwise be silently rounded (or crash `.astype(int)` on
    # NaN). A matrix accepted as raw MUST be all finite non-negative integers
    # everywhere; a coerced (allow_nonraw) matrix may carry fractional/negative
    # values the user opted to round, but NaN/inf still cannot be coerced.
    full_check = validate_raw_count_matrix(counts)
    if not full_check["valid"]:
        block = count_source == "raw_counts" or full_check["n_nonfinite"] > 0
        if block:
            examples = full_check.get("examples") or []
            example_txt = (f" e.g. {examples}" if examples else "")
            return None, warnings, {
                "refused":    True,
                "error_type": "InvalidRawCounts",
                "kind":       info["kind"],
                "raw_count_score": info.get("raw_count_score"),
                "confidence": info.get("confidence"),
                "sub_scores": info.get("sub_scores", {}),
                "score_basis": info.get("score_basis", {}),
                "full_validation": {k: full_check[k] for k in (
                    "n_nonfinite", "n_negative", "n_noninteger", "n_offending")},
                "details": (
                    f"Count matrix failed full raw-count validation: "
                    f"{full_check['reason']}{example_txt}. classify_matrix scores a "
                    f"200-row sample, but DESeq2 needs every value to be a finite "
                    f"non-negative integer; a raw matrix cannot contain a "
                    f"non-integer/negative value and no matrix can carry NaN/inf "
                    f"into integer rounding. Supply a clean raw-count matrix."
                ),
            }

    # Round to integers (required by DESeq2)
    counts = counts.round().astype(int)
    counts.columns = counts.columns.astype(str)

    meta_out = {
        "count_source": count_source,
        "kind": info["kind"],
        "raw_count_score": info.get("raw_count_score"),
        "confidence": info.get("confidence"),
        "sub_scores": info.get("sub_scores", {}),
        "score_basis": info.get("score_basis", {}),
    }
    if count_union_meta is not None:
        meta_out["count_union"] = count_union_meta
    return counts, warnings, meta_out


def _enforce_metadata_correspondence(counts, metadata, excluded_samples=None):
    """B5: every count-matrix column must have an aligned metadata row.

    Partial metadata must never silently reduce the analysis. When an explicit
    metadata TSV covers only a subset of the columns, ``_load_or_infer_metadata``
    returns ``meta.loc[common]`` and ``_run_deseq2`` later subsets
    ``counts[meta_sub.index]`` — dropping the unmatched columns from the fit with
    no error and no disclosure. This gate closes that path.

    Count columns without a metadata row fail closed with the exact list, UNLESS
    they are named in ``excluded_samples`` — an explicit, audited exclusion. Only
    genuine orphan columns may be excluded; naming a sample that has metadata, or
    a name that is not a count column, is a misuse and also fails closed.

    Returns ``(counts_kept, disclosure)``. Raises ``ValueError`` on any
    unauthorized or invalid correspondence gap.
    """
    count_cols = [str(c) for c in counts.columns]
    if len(set(count_cols)) != len(count_cols):
        dup_cols = sorted(
            {c for c in count_cols if count_cols.count(c) > 1}
        )
        raise ValueError(
            "Count matrix contains duplicate sample column names: "
            f"{dup_cols}. Every count column must have a unique sample ID."
        )
    meta_ids = [str(i) for i in metadata.index]
    if len(set(meta_ids)) != len(meta_ids):
        dup_meta = sorted(
            {i for i in meta_ids if meta_ids.count(i) > 1}
        )
        raise ValueError(
            "Metadata table contains duplicate sample IDs: "
            f"{dup_meta}. Collapse technical lanes or fix duplicate rows "
            "before differential expression."
        )
    have_meta = set(meta_ids)
    orphans = [c for c in count_cols if c not in have_meta]

    # De-duplicate the requested exclusions while preserving order.
    seen: set[str] = set()
    excluded: list[str] = []
    for sample in (excluded_samples or []):
        name = str(sample)
        if name not in seen:
            seen.add(name)
            excluded.append(name)

    orphan_set = set(orphans)
    invalid = [name for name in excluded if name not in orphan_set]
    if invalid:
        raise ValueError(
            "excluded_samples may only name count-matrix columns that lack "
            f"metadata; these are not unmatched count columns: {invalid}. "
            f"Unmatched columns are: {orphans if orphans else 'none'}."
        )

    unauthorized = [c for c in orphans if c not in seen]
    if unauthorized:
        raise ValueError(
            "count-matrix columns have no aligned metadata row: "
            f"{unauthorized}. Provide metadata for every sample (total "
            "correspondence), or list them explicitly in excluded_samples to "
            "drop them from the analysis with an audit record. DE will not run "
            "on a silently reduced sample set."
        )

    disclosure = {
        "excluded_samples": list(excluded),
        "n_excluded":       len(excluded),
        "n_count_columns":  len(count_cols),
        "n_analyzed":       len(count_cols) - len(excluded),
    }
    if excluded:
        kept = [c for c in counts.columns if str(c) not in seen]
        counts = counts[kept]
    return counts, disclosure


def _metadata_inference_allowed(params: dict) -> bool:
    """Explicit legacy/dev escape hatch for column-name metadata inference."""
    opt_in = (
        params.get("allow_inferred_metadata")
        or params.get("allow_metadata_inference")
        or params.get("legacy_metadata_inference")
    )
    if opt_in:
        return True
    return os.environ.get(
        "ARIA_ALLOW_BULK_METADATA_INFERENCE", ""
    ).strip().lower() in {"1", "true", "yes", "on"}


def _load_or_infer_metadata(counts, metadata_file: str,
                              design_factor: str,
                              allow_inference: bool = False) -> tuple:
    """
    Load explicit metadata or infer groups from column names.

    Supported naming patterns:
      ctrl_1, ctrl_2, treat_1, treat_2     → condition = {ctrl, treat}
      WT_rep1, KO_rep1                      → condition = {WT, KO}
      sample_A_1, sample_B_1               → condition = {A, B}
    """
    import pandas as pd
    import re
    warnings = []

    counts.columns = counts.columns.astype(str)
    if not counts.columns.is_unique:
        dup_cols = sorted(set(counts.columns[counts.columns.duplicated()].tolist()))
        raise ValueError(
            "Count matrix contains duplicate sample column names: "
            f"{dup_cols}. Every sample column must have a unique ID."
        )

    # Try loading explicit metadata file
    if metadata_file and Path(metadata_file).exists():
        meta = None
        try:
            meta = pd.read_csv(metadata_file, sep="\t", index_col=0)
            meta.index = meta.index.astype(str)
        except Exception as e:
            fallback = (
                "Attempting automatic detection."
                if allow_inference
                else "Automatic detection is disabled in production."
            )
            warnings.append(
                f"Metadata file load failed: {e}. "
                f"{fallback}"
            )
        if meta is not None:
            if not meta.index.is_unique:
                dup_ids = sorted(set(meta.index[meta.index.duplicated()].tolist()))
                raise ValueError(
                    f"Metadata file '{metadata_file}' contains duplicate "
                    f"sample IDs: {dup_ids}. Each sample must appear at most "
                    "once before differential expression."
                )
            # Align to count matrix samples
            common = [s for s in counts.columns if s in meta.index]
            if len(common) < 2:
                fallback = (
                    "Falling back to automatic detection."
                    if allow_inference
                    else "Automatic detection is disabled in production."
                )
                warnings.append(
                    f"Metadata file has {len(common)} matching samples. "
                    f"{fallback}"
                )
            else:
                return meta.loc[common], warnings
    elif metadata_file:
        warnings.append(f"Metadata file not found: {metadata_file}.")

    if not allow_inference:
        return None, warnings + [
            "Production bulk RNA DE requires a valid metadata TSV aligned to "
            "the count matrix; column-name metadata inference is disabled."
        ]

    # Automatic group detection from column names
    samples = list(counts.columns)
    groups  = _infer_groups(samples)

    if groups is None or len(set(groups.values())) < 2:
        return None, warnings + [
            "Could not infer experimental groups from sample names. "
            "Please provide a metadata TSV file with columns: "
            "sample, condition (and optionally: batch, replicate)."
        ]

    meta = pd.DataFrame({
        "sample":      samples,
        design_factor: [groups[s] for s in samples],
    }, index=samples)

    warnings.append(
        f"Experimental groups inferred from sample names: "
        f"{dict(set((v,sum(1 for x in groups.values() if x==v)) for v in set(groups.values())))}. "
        f"If this is incorrect, provide a metadata file."
    )

    return meta, warnings


def _aggregate_technical_replicates(
    counts,
    metadata,
    *,
    design_factor: str,
    unit_col: str,
    covariates: list | None = None,
) -> tuple:
    """Sum raw-count libraries into condition-scoped biological units.

    Every design value must be invariant within a unit.  This is deliberately
    fail-closed: a unit spanning conditions/covariates or an unmapped library
    cannot be repaired by guessing.
    """
    import pandas as pd
    from aria.utils.design_matrix import validate_design_matrix

    if unit_col not in metadata.columns:
        raise ValueError(
            f"technical replicate column '{unit_col}' is absent from metadata"
        )
    missing_metadata = [
        sample for sample in counts.columns if sample not in metadata.index
    ]
    if missing_metadata:
        raise ValueError(
            "technical replicate metadata is missing count libraries: "
            f"{missing_metadata}"
        )
    meta = metadata.loc[list(counts.columns)].copy()
    if (
        meta[unit_col].isna().any()
        or meta[unit_col].astype(str).str.strip().eq("").any()
    ):
        raise ValueError("every technical library must map to a biological unit")

    design_cols = [design_factor, *[c for c in (covariates or []) if c]]
    absent = [column for column in design_cols if column not in meta.columns]
    if absent:
        raise ValueError(
            "technical replicate metadata is missing design columns: "
            f"{absent}"
        )

    unit_rows: list[dict] = []
    aggregated: dict[str, object] = {}
    members: dict[str, list[str]] = {}
    for unit, rows in meta.groupby(unit_col, sort=False, observed=True):
        unit_id = str(unit)
        sample_ids = [str(sample) for sample in rows.index]
        row = {}
        for column in design_cols:
            values = rows[column].dropna().astype(str).unique().tolist()
            if len(values) != 1:
                raise ValueError(
                    f"biological unit '{unit_id}' spans multiple {column} values: "
                    f"{values}"
                )
            row[column] = values[0]
        aggregated[unit_id] = counts[sample_ids].sum(axis=1)
        unit_rows.append({unit_col: unit_id, **row})
        members[unit_id] = sample_ids

    unit_counts = pd.DataFrame(aggregated, index=counts.index)
    unit_metadata = pd.DataFrame(unit_rows).set_index(unit_col)
    design_check = validate_design_matrix(
        unit_metadata,
        condition_col=design_factor,
        covariates=list(covariates or []),
        min_replicates_per_condition=2,
    )
    provenance = {
        "ran": True,
        "method": "sum_raw_counts_by_biological_unit",
        "unit_column": unit_col,
        "n_input_libraries": int(counts.shape[1]),
        "n_biological_units": int(unit_counts.shape[1]),
        "replicates_per_condition": {
            str(group): int(count)
            for group, count in unit_metadata[design_factor].value_counts().items()
        },
        "residual_degrees_of_freedom": design_check.get(
            "residual_degrees_of_freedom", 0
        ),
        "design_rank": design_check.get("rank", 0),
        "members": members,
    }
    return unit_counts, unit_metadata, provenance


def _infer_groups(samples: list) -> dict | None:
    """
    Detect condition groups from sample names using regex patterns.
    Returns {sample_name: group_label} or None if detection fails.
    """
    import re

    # Pattern 1: condition_replicate (ctrl_1, treat_1, ctrl_2, treat_2)
    p1 = re.compile(r'^([A-Za-z][A-Za-z0-9]+)[_\-](\d+)$')
    matches = {s: m.group(1) for s in samples
               if (m := p1.match(s))}
    if len(matches) == len(samples) and len(set(matches.values())) >= 2:
        return matches

    # Pattern 2: prefix_suffix with letters (WT_rep1, KO_rep1)
    p2 = re.compile(r'^([A-Za-z][A-Za-z0-9]+)[_\-]([Rr]ep\d+|[A-Za-z]\d*)$')
    matches = {s: m.group(1) for s in samples
               if (m := p2.match(s))}
    if len(matches) == len(samples) and len(set(matches.values())) >= 2:
        return matches

    # Pattern 3: split by last underscore, use prefix as group
    if all("_" in s for s in samples):
        groups = {s: "_".join(s.split("_")[:-1]) for s in samples}
        if len(set(groups.values())) >= 2:
            return groups

    # Pattern 4: alphabetical prefix before first digit
    p4 = re.compile(r'^([A-Za-z][A-Za-z0-9\-]*?)(\d.*)$')
    matches = {s: m.group(1).rstrip("_-") for s in samples
               if (m := p4.match(s))}
    if len(matches) == len(samples) and len(set(matches.values())) >= 2:
        return matches

    return None


def _resolve_comparison(metadata, design_factor: str,
                          comparison: dict) -> tuple:
    """
    Resolve which groups to compare.
    P0-5: never infer or substitute numerator/reference levels.
    """
    warnings = []
    groups   = sorted(metadata[design_factor].unique())

    if len(groups) < 2:
        return comparison, [
            f"Only one group found in '{design_factor}': {groups}. "
            f"Cannot run differential expression."
        ]

    num = comparison.get("numerator",   "")
    den = comparison.get("denominator", "")

    if not num or not den:
        warnings.append(
            "Explicit numerator and denominator are required; no comparison "
            "was inferred from group names."
        )
        return {"numerator": "", "denominator": ""}, warnings

    if num not in groups:
        warnings.append(f"Numerator '{num}' not in {groups}.")
    if den not in groups:
        warnings.append(f"Denominator '{den}' not in {groups}.")

    return {"numerator": num, "denominator": den}, warnings


# ── Sample QC ─────────────────────────────────────────────────────────────────

# ══════════════════════════════════════════════════════════════════════════════
# METHODOLOGY LAYER — explicit normalization & dimensionality reduction
# ══════════════════════════════════════════════════════════════════════════════
#
# Decisions baked into this module (all justified in methods section of report):
#
#   1. DESeq2 DE test        → raw counts (DESeq2 normalizes internally)
#   2. PCA / MDS             → VST + top N variable protein_coding genes
#   3. Heatmap (padj top)    → log2(counts+1) + row z-score
#   4. Heatmap (|log2FC| top)→ log2(counts+1) + row z-score (NEW in v3.8)
#   5. TPM (supplementary)   → gene-length × library-size normalized (NEW)
#
# Why VST over log2(raw+1) + StandardScaler for PCA:
#   - VST is the DESeq2 authors' recommendation (Love, Huber, Anders 2014).
#   - Produces homoscedastic values — equal variance across expression levels.
#   - StandardScaler z-scores per gene over-weight low-variance genes and
#     compress high-variance (biologically meaningful) genes.
#   - log2(raw+1) doesn't correct for library size differences across samples.
#
# Why top N variable protein_coding only for DR:
#   - Pseudogenes/rRNAs dominate raw variance without biological meaning.
#   - Low-variance genes are noise; top 2000 captures informative signal.
#   - Matches the standard DESeq2 vignette workflow.


