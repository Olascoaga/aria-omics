"""Verified same-cell multiome assembly from ARC, paired matrices, or raw FASTQs."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from bisect import bisect_right
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script
from aria.utils.mudata_io import (
    MultiomeInputError,
    _open_text,
    _read_10x_matrix,
    _tenx_modality,
    _select_raw_counts_matrix,
    _write_paired_atomic,
    _stage_cache,
    _stage_complete,
    paired_input_files,
    read_h5mu_paired,
    validate_paired_multiome,
)


def _fail(kind, details):
    raise MultiomeInputError(kind, str(details))


def _load_matrix(path, feature_type):
    import anndata as ad

    if not Path(path).exists():
        _fail("FileNotFound", path)
    if Path(path).suffix.lower() == ".h5ad":
        return ad.read_h5ad(path)
    return _tenx_modality(_read_10x_matrix(path), feature_type)


def _raw(adata, label):
    import numpy as np
    from scipy import sparse

    selected = _select_raw_counts_matrix(
        adata, modality_label=label, require_raw_counts=True
    )
    if not selected["ok"]:
        _fail(selected["error_type"], selected["details"])
    matrix = sparse.csr_matrix(selected["matrix"])
    if matrix.data.size and matrix.data.max() >= np.iinfo(np.int64).max:
        _fail("NonIntegerCounts", "Counts exceed int64 storage")
    adata.X = matrix.astype(np.int64)
    adata.layers["counts"] = adata.X.copy()
    return adata


def _fragment_rows(path):
    """Stream validated, half-open fragments and their observed multiplicity."""
    with _open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                _fail("InvalidFragments", "Fragments need chrom, start, end, barcode")
            try:
                start, end = int(parts[1]), int(parts[2])
                count = int(parts[4]) if len(parts) > 4 else 1
            except ValueError:
                _fail("InvalidFragments", "Noninteger fragment coordinates/count")
            if start < 0 or end <= start or count < 1:
                _fail("InvalidFragments", "Invalid fragment interval/count")
            yield parts[0], start, end, parts[3], count


def _quantify_fragments(fragments_path, peaks_path):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    intervals = []
    with _open_text(peaks_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            row = line.split()
            chrom, start, end = row[0], int(row[1]), int(row[2])
            if start < 0 or end <= start:
                _fail("NonGenomicPeakFeatures", line)
            intervals.append((chrom, start, end))
    if not intervals or len(set(intervals)) != len(intervals):
        _fail("NonGenomicPeakFeatures", "BED must contain unique, valid peaks")
    # Prefix maxima permit early termination even for nested/overlapping peaks.
    index = {}
    for col, (chrom, start, end) in enumerate(intervals):
        index.setdefault(chrom, []).append((start, end, col))
    for chrom, values in index.items():
        values.sort()
        starts, maxima, running = [], [], -1
        for start, end, _ in values:
            starts.append(start)
            running = max(running, end)
            maxima.append(running)
        index[chrom] = (values, starts, maxima)
    cells, counts = {}, {}
    for chrom, start, end, barcode, count in _fragment_rows(fragments_path):
        row = cells.setdefault(barcode, len(cells))
        if chrom not in index:
            continue
        values, starts, maxima = index[chrom]
        pos = bisect_right(starts, end - 1) - 1
        while pos >= 0 and maxima[pos] > start:
            pstart, pend, col = values[pos]
            if pstart < end and pend > start:
                key = (row, col)
                counts[key] = counts.get(key, 0) + count
            pos -= 1
    keys = list(counts)
    matrix = sparse.csr_matrix(
        (
            np.array([counts[k] for k in keys], dtype=np.int64),
            ([k[0] for k in keys], [k[1] for k in keys]),
        ),
        shape=(len(cells), len(intervals)),
        dtype=np.int64,
    )
    return ad.AnnData(
        matrix,
        obs=pd.DataFrame(index=list(cells)),
        var=pd.DataFrame(index=[f"{c}:{s}-{e}" for c, s, e in intervals]),
    )


def _split_suffix(barcode):
    return (barcode[:-2], "-1") if barcode.endswith("-1") else (barcode, "")


def _revcomp(barcode):
    base, suffix = _split_suffix(barcode)
    return base.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1] + suffix


def _check_bijection(rna, atac):
    if (
        len(rna) != len(atac)
        or len(set(rna)) != len(rna)
        or len(set(atac)) != len(atac)
    ):
        _fail(
            "NonBijectiveBarcodeMap",
            "Barcode maps must have unique RNA and ATAC entries and equal lengths",
        )


def _translate_barcodes(rna, atac, params):
    import csv

    namespace = params.get("barcode_namespace", "shared_barcodes")
    rnames, anames = list(map(str, rna.obs_names)), list(map(str, atac.obs_names))
    if len(set(rnames)) != len(rnames) or len(set(anames)) != len(anames):
        _fail("DuplicateBarcodes", "Input cell barcodes must be unique")
    translation = {"mode": namespace, "orientation": "direct", "n_translated": 0}
    if namespace == "shared_barcodes":
        # Preserve exact matches; only translate a complete, unique suffix bijection.
        if set(rnames) != set(anames) and params.get("normalize_10x_suffix", True):
            rn, an = [
                [_split_suffix(b)[0] for b in names] for names in (rnames, anames)
            ]
            if (
                len(set(rn)) == len(rn)
                and len(set(an)) == len(an)
                and set(rn) == set(an)
            ):
                mapping = dict(zip(rn, rnames))
                atac.obs_names = [mapping[b] for b in an]
                translation["n_translated"] = sum(
                    a != b for a, b in zip(anames, atac.obs_names)
                )
        return translation
    if namespace == "barcode_map_tsv":
        with _open_text(params["barcode_map_path"]) as handle:
            lines = [
                line for line in handle if line.strip() and not line.startswith("#")
            ]
        rows = list(
            csv.reader(lines, delimiter="\t" if lines and "\t" in lines[0] else ",")
        )
        if rows and rows[0] == ["rna_barcode", "atac_barcode"]:
            rows = rows[1:]
        if not rows or any(len(row) != 2 or not all(row) for row in rows):
            _fail(
                "NonBijectiveBarcodeMap",
                "Expected two columns: rna_barcode, atac_barcode",
            )
        rmap, amap = ([row[i] for row in rows] for i in (0, 1))
        _check_bijection(rmap, amap)
        mapping = dict(zip(amap, rmap))
    elif namespace == "10x_arc_whitelist_pair":

        def whitelist(key):
            with _open_text(params[key]) as handle:
                entries = [line.strip() for line in handle]
            if not entries or any(not entry for entry in entries):
                _fail(
                    "NonBijectiveBarcodeMap",
                    "Whitelists must contain one nonempty barcode per line",
                )
            return entries

        rmap, amap = whitelist("gex_whitelist_path"), whitelist("atac_whitelist_path")
        _check_bijection(rmap, amap)
        # Whitelists normally omit GEM suffixes. Strip only -1, preserving it on output.
        rb, ab = [[_split_suffix(b)[0] for b in values] for values in (rmap, amap)]
        _check_bijection(rb, ab)
        direct = dict(zip(ab, rmap))
        reverse = {_revcomp(b): a for b, a in direct.items()}

        def translated(table):
            out = []
            for barcode in anames:
                base, suffix = _split_suffix(barcode)
                target = table.get(base)
                if target is None:
                    out.append(None)
                else:
                    target_base, target_suffix = _split_suffix(target)
                    out.append(target_base + (suffix or target_suffix))
            return out

        candidates = {
            "direct": translated(direct),
            "reverse_complement": translated(reverse),
        }
        orientation = params.get("atac_barcode_orientation") or "auto"
        if orientation not in {"auto", *candidates}:
            _fail(
                "InvalidParameter",
                "atac_barcode_orientation must be direct, reverse_complement, or auto",
            )
        if orientation == "auto":
            rset = set(rnames)
            orientation = max(
                candidates, key=lambda key: sum(b in rset for b in candidates[key])
            )
        translation["orientation"] = orientation
        mapping = {
            original: target
            for original, target in zip(anames, candidates[orientation])
            if target is not None
        }
    else:
        _fail("InvalidParameter", f"Unknown barcode_namespace: {namespace}")
    # Unmapped cells cannot establish pairing; do not let identity fallbacks create false matches.
    keep = [i for i, b in enumerate(anames) if b in mapping]
    translated = [mapping[anames[i]] for i in keep]
    if len(set(translated)) != len(translated):
        _fail("NonBijectiveBarcodeMap", "Translation produces duplicate ATAC barcodes")
    translation["n_translated"] = len(keep)
    return translation, keep, translated


def _tool(name, env, explicit=None):
    import shutil

    if explicit:
        candidate = shutil.which(str(explicit))
        return candidate
    prefixes = [Path(sys.prefix)]
    if os.environ.get("CONDA_PREFIX"):
        prefixes.insert(0, Path(os.environ["CONDA_PREFIX"]))
    candidates = [
        candidate
        for prefix in prefixes
        for candidate in (
            prefix / "bin" / name,
            prefix.parent / env / "bin" / name,
            prefix / "envs" / env / "bin" / name,
        )
    ]
    return next(
        (str(p) for p in candidates if p.is_file() and os.access(p, os.X_OK)),
        shutil.which(name),
    )


def _raw_fastqs(params, output_dir):
    """Execute local tools or consume explicitly delegated real count/fragment outputs."""
    import subprocess

    platform = params.get("platform")
    allowed = {
        "10x_multiome_arc_v1",
        "share_seq",
        "paired_tag",
        "issaac_seq",
        "dogma_seq",
        "custom_paired",
    }
    if platform not in allowed:
        _fail(
            "InvalidParameter", f"Raw assembly requires platform in {sorted(allowed)}"
        )
    rna = params.get("rna_counts_out")
    fragments = params.get("atac_fragments_out")
    if rna and not Path(rna).exists():
        _fail("FileNotFound", rna)
    if fragments and not Path(fragments).is_file():
        _fail("FileNotFound", fragments)
    kb = None if rna else _tool("kb", "aria-ingestion-env", params.get("kb_bin"))
    chromap = (
        None
        if fragments
        else _tool("chromap", "aria-atacseq-env", params.get("chromap_bin"))
    )
    if (not rna and not kb) or (not fragments and not chromap):
        _fail(
            "MissingDependency",
            "Raw multiome quantification requires kb and chromap (or delegated outputs)",
        )
    required = []
    if not rna:
        required += [
            params.get("kb_index") or params.get("index_path"),
            params.get("t2g_path"),
        ]
    if not fragments:
        required += [params.get("genome_fasta"), params.get("chromap_index")]
    if any(not value or not Path(value).is_file() for value in required):
        _fail(
            "MissingReference",
            "Provide existing kb_index/t2g_path and genome_fasta/chromap_index as needed",
        )
    for key, needed in (("rna_fastqs", not rna), ("atac_fastqs", not fragments)):
        fastqs = params.get(key)
        if needed and (
            not isinstance(fastqs, list)
            or len(fastqs) < 2
            or any(not Path(f).is_file() for f in fastqs)
        ):
            _fail("MissingInput", f"{key} requires existing paired FASTQ files")
    output_dir.mkdir(parents=True, exist_ok=True)

    def execute(command):
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            timeout=int(params.get("timeout_seconds", 86400)),
        )
        if result.returncode:
            _fail(
                "QuantificationFailed",
                f"{Path(command[0]).name}: {result.stderr[-2000:]}",
            )

    if not rna:
        from aria.utils.raw_ingestion import build_kb_count_command

        technology = params.get("technology") or (
            "10xv3" if platform == "10x_multiome_arc_v1" else None
        )
        if not technology:
            _fail(
                "MissingReference",
                "Non-ARC raw GEX requires an explicit kb technology/read geometry",
            )
        import tempfile

        # Each invocation gets an empty tool directory; interrupted outputs cannot be reused.
        kb_dir = Path(tempfile.mkdtemp(prefix="kb_counts.", dir=output_dir))
        command = build_kb_count_command(
            fastq_files=params["rna_fastqs"],
            index_path=params.get("kb_index") or params["index_path"],
            t2g_path=params["t2g_path"],
            output_dir=str(kb_dir),
            chemistry=technology,
            technology=technology,
            threads=int(params.get("threads", 8)),
        )
        command[0] = kb
        execute(command)
        import anndata as ad
        import pandas as pd
        from scipy import sparse
        from scipy.io import mmread

        matrix_dir = next(
            (
                kb_dir / name
                for name in ("counts_filtered", "counts_unfiltered", "")
                if (kb_dir / name / "cells_x_genes.mtx").is_file()
            ),
            None,
        )
        if matrix_dir is None:
            _fail("QuantificationFailed", "kb did not produce a cells_x_genes matrix")
        with open(matrix_dir / "cells_x_genes.barcodes.txt") as handle:
            barcodes = [line.strip() for line in handle]
        with open(matrix_dir / "cells_x_genes.genes.txt") as handle:
            genes = [line.strip().split("\t")[0] for line in handle]
        rna_adata = ad.AnnData(
            sparse.csr_matrix(mmread(matrix_dir / "cells_x_genes.mtx")),
            obs=pd.DataFrame(index=barcodes),
            var=pd.DataFrame(index=genes),
        )
        rna_adata.var["gene_ids"] = genes
        rna_adata.var_names_make_unique()
    else:
        rna_adata = _load_matrix(rna, "Gene Expression")
    if not fragments:
        fragments = str(output_dir / "chromap.fragments.tsv")
        Path(fragments).unlink(missing_ok=True)
        fastqs = params["atac_fastqs"]
        if len(fastqs) != 3:
            _fail("MissingInput", "ATAC requires [R1, R2, barcode_FASTQ]")
        command = [
            chromap,
            "--preset",
            "atac",
            "-x",
            params["chromap_index"],
            "-r",
            params["genome_fasta"],
            "-1",
            fastqs[0],
            "-2",
            fastqs[1],
            "-b",
            fastqs[2],
            "-o",
            fragments,
            "-t",
            str(params.get("threads", 8)),
        ]
        if params.get("read_format"):
            command += ["--read-format", params["read_format"]]
        elif platform != "10x_multiome_arc_v1":
            _fail(
                "MissingReference",
                "Non-ARC raw ATAC requires explicit chromap read_format",
            )
        if params.get("atac_whitelist_path"):
            command += ["--barcode-whitelist", params["atac_whitelist_path"]]
        execute(command)
        if not Path(fragments).is_file():
            _fail("QuantificationFailed", "chromap did not produce fragments")
    peaks = params.get("peaks_path")
    if not peaks:
        from aria.scripts.chromatin_peaks import chromatin_peaks

        called = chromatin_peaks(
            {
                "data_type": "scATAC",
                "files": [fragments],
                "output_dir": str(output_dir / "peaks"),
                "genome": params.get("genome", "hg38"),
                "macs3_params": {"format": "BEDPE"},
            }
        )
        if called.get("status") != "success":
            _fail(
                called.get("error_type", "MissingDependency"),
                called.get("details", "Peak calling failed"),
            )
        peaks = called.get("consensus_peaks_path") or called["peaks_path"]
    return rna_adata, _quantify_fragments(fragments, peaks)


def integration_assemble_multiome(params: dict[str, Any]) -> dict[str, Any]:
    try:
        import mudata

        mode = params.get("mode")
        if mode in (None, "", "auto"):
            mode = (
                "cellranger_arc"
                if params.get("arc_path")
                else (
                    "raw_multiome_fastq"
                    if params.get("rna_fastqs") or params.get("rna_counts_out")
                    else "paired_matrices"
                )
            )
        output_dir = Path(params.get("output_dir") or "/tmp/aria_multiome")
        output = str(
            params.get("output_path") or output_dir / "assembled_multiome.h5mu"
        )
        inputs = []
        for key in ("arc_path", "rna_path", "atac_path", "rna_counts_out"):
            if params.get(key):
                inputs += paired_input_files(key, params[key])
        for key in (
            "fragments_path",
            "peaks_path",
            "barcode_map_path",
            "gex_whitelist_path",
            "atac_whitelist_path",
            "atac_fragments_out",
            "kb_index",
            "index_path",
            "t2g_path",
            "genome_fasta",
            "chromap_index",
        ):
            if params.get(key):
                inputs.append((key, str(params[key])))
        for key in ("rna_fastqs", "atac_fastqs"):
            inputs += [
                (f"{key}:{i}", str(p)) for i, p in enumerate(params.get(key) or [])
            ]
        if Path(output).resolve() in {Path(p).resolve() for _, p in inputs}:
            _fail("InvalidOutputPath", "Output must differ from inputs")
        cached, state = _stage_cache(
            "integration_assemble_multiome",
            output_dir,
            inputs,
            dict(params, mode=mode),
            params.get("resume", True),
        )
        if cached:
            return cached
        translation = {
            "mode": "shared_barcodes",
            "orientation": "direct",
            "n_translated": 0,
        }
        if mode == "cellranger_arc":
            if not params.get("arc_path"):
                _fail("MissingInput", "cellranger_arc requires arc_path")
            paired = read_h5mu_paired(params["arc_path"])
            if paired["status"] != "success":
                return paired
            rna, atac = paired["rna_adata"], paired["atac_adata"]
        elif mode in {"paired_matrices", "raw_multiome_fastq"}:
            if mode == "raw_multiome_fastq":
                rna, atac = _raw_fastqs(params, output_dir)
            else:
                if not params.get("rna_path") or not (
                    params.get("atac_path")
                    or (params.get("fragments_path") and params.get("peaks_path"))
                ):
                    _fail(
                        "MissingInput",
                        "paired_matrices requires rna_path and atac_path or fragments_path + peaks_path",
                    )
                rna = _load_matrix(params["rna_path"], "Gene Expression")
                atac = (
                    _load_matrix(params["atac_path"], "Peaks")
                    if params.get("atac_path")
                    else _quantify_fragments(
                        params["fragments_path"], params["peaks_path"]
                    )
                )
            n_atac_original = atac.n_obs
            translated = _translate_barcodes(rna, atac, params)
            if isinstance(translated, tuple):
                translation, rows, names = translated
                atac = atac[rows].copy()
                atac.obs_names = names
            else:
                translation = translated
            paired = validate_paired_multiome(mudata.MuData({"rna": rna, "atac": atac}))
            if paired["status"] != "success":
                return paired
            paired["n_atac_cells"] = n_atac_original
            paired["shared_fraction_atac"] = (
                paired["n_shared_cells"] / n_atac_original if n_atac_original else 0.0
            )
        else:
            _fail("InvalidParameter", f"Unknown intake mode: {mode}")
        minimum = int(params.get("min_shared_cells", 10))
        fraction = float(params.get("min_shared_barcode_fraction", 0.25))
        if minimum < 1 or not 0 <= fraction <= 1:
            _fail("InvalidParameter", "Invalid shared cell thresholds")
        if paired["n_shared_cells"] < minimum:
            _fail(
                "InsufficientSharedCells",
                f"Need {minimum} shared cells; found {paired['n_shared_cells']}",
            )
        if (
            min(paired["shared_fraction_rna"], paired["shared_fraction_atac"])
            < fraction
        ):
            _fail("LowBarcodeOverlap", f"Shared fraction below {fraction}")
        shared = paired["shared_barcodes"]
        rna = _raw(rna[shared].copy(), "rna")
        atac = _raw(atac[shared].copy(), "atac")
        from aria.utils.scatac_cohort_union import _peak_intervals

        atac.var_names = [f"{c}:{s}-{e}" for c, s, e in _peak_intervals(atac)]
        if not rna.var_names.is_unique or not atac.var_names.is_unique:
            _fail("DuplicateFeatureNames", "Feature names must be unique")
        if params.get("fragments_path"):
            totals, mito = {}, {}
            for chrom, _, _, barcode, count in _fragment_rows(params["fragments_path"]):
                totals[barcode] = totals.get(barcode, 0) + count
                if chrom.lower() in {"chrm", "chrmt", "mt", "m"}:
                    mito[barcode] = mito.get(barcode, 0) + count
            atac.obs["atac_fragments"] = [totals.get(b, float("nan")) for b in shared]
            atac.obs["atac_mito_frac"] = [
                mito.get(b, 0) / totals[b] if b in totals else float("nan")
                for b in shared
            ]
        mdata = mudata.MuData({"rna": rna, "atac": atac})
        _write_paired_atomic(mdata, output)
        result = {
            "status": "success",
            "output_path": output,
            "mode": mode,
            "platform": params.get("platform") or "unspecified",
            **{
                key: paired[key]
                for key in (
                    "n_rna_cells",
                    "n_atac_cells",
                    "n_shared_cells",
                    "shared_fraction_rna",
                    "shared_fraction_atac",
                )
            },
            "n_genes": rna.n_vars,
            "n_peaks": atac.n_vars,
            "barcode_translation": translation,
            "validation_level": "production",
            "warnings": [],
        }
        return _stage_complete(result, state)
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
    run_script(integration_assemble_multiome)
