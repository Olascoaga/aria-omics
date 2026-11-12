"""
ARIA Multiplexed Well-Barcoded Bulk RNA-seq Quantification via STARsolo (ADR-075)
---------------------------------------------------------------------------------
Quantifies pooled plate-based bulk RNA-seq libraries (BRB-seq, DRUG-seq,
Prime-seq, QuantSeq-Pool, CEL-Seq2 bulk plate, multiplexed STRT) directly from
raw paired FASTQs into a gene × sample UMI-deduplicated count matrix using
``STAR --soloType CB_UMI_Simple`` (``1MM_Directional_UMItools``).

Executed inside ``aria-rnaseq-env`` via ``EnvironmentManager``.
Content-addressed resume is enforced via ``aria/utils/stage_manifest.py`` (A4).
"""

from __future__ import annotations

import gzip
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script
from aria.scripts.rna_align import (
    _INDEX_MANIFEST,
    _build_star_index,
    _index_exists,
    _parse_star_log,
    _sa_index_nbases,
    _star_version,
)
from aria.utils.rna_fastq_protocol import validate_pooled_plate_manifest
from aria.utils.stage_manifest import (
    build_stage_manifest,
    stage_is_current,
    write_stage_manifest,
)


def _open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def _resolve_pooled_fastqs(
    fastq_files: list[str],
    manifest: dict[str, Any],
) -> tuple[list[str], list[str], str | None]:
    """Return ``(r1_files, r2_files, error_message)`` sorted deterministically."""
    r1_explicit = [str(p) for p in (manifest.get("r1_fastqs") or [])]
    r2_explicit = [str(p) for p in (manifest.get("r2_fastqs") or [])]
    if r1_explicit and r2_explicit:
        if len(r1_explicit) != len(r2_explicit):
            return [], [], (
                f"pooled_plate_manifest has unequal r1_fastqs ({len(r1_explicit)}) "
                f"and r2_fastqs ({len(r2_explicit)})."
            )
        return r1_explicit, r2_explicit, None

    r1: list[str] = []
    r2: list[str] = []
    for raw in sorted(str(p) for p in (fastq_files or [])):
        name = Path(raw).name
        stem = re.sub(r"\.(?:fastq|fq)(?:\.gz)?$", "", name, flags=re.IGNORECASE)
        if re.search(r"_(?:R1|1)(?:_\d+)?$", stem, re.IGNORECASE):
            r1.append(raw)
        elif re.search(r"_(?:R2|2)(?:_\d+)?$", stem, re.IGNORECASE):
            r2.append(raw)
    if not r1 or not r2:
        return [], [], (
            "Pooled plate bulk RNA-seq requires paired R1 and R2 FASTQ files "
            "(one carrying well barcode + UMI and the other carrying cDNA)."
        )
    if len(r1) != len(r2):
        return [], [], (
            f"Unequal R1 ({len(r1)}) and R2 ({len(r2)}) FASTQ lane counts for "
            "pooled_plate_bulk."
        )
    return r1, r2, None


def _parse_starsolo_raw_matrix(
    raw_dir: Path,
    barcode_to_sample: dict[str, str],
    output_dir: Path,
    warnings: list[str],
) -> tuple[Path, Path | None, int, dict[str, int]]:
    """Convert STARsolo ``Solo.out/Gene/raw`` MatrixMarket into ``counts_matrix.tsv``."""
    import pandas as pd

    bc_path = next(
        (raw_dir / n for n in ("barcodes.tsv", "barcodes.tsv.gz") if (raw_dir / n).exists()),
        None,
    )
    feat_path = next(
        (raw_dir / n for n in ("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")
         if (raw_dir / n).exists()),
        None,
    )
    mtx_path = next(
        (raw_dir / n for n in ("matrix.mtx", "matrix.mtx.gz") if (raw_dir / n).exists()),
        None,
    )
    if not (bc_path and feat_path and mtx_path):
        raise FileNotFoundError(
            f"STARsolo raw Gene output incomplete in {raw_dir}: expected "
            "barcodes.tsv, features.tsv, and matrix.mtx."
        )

    with _open_text(bc_path) as fh:
        barcodes = [ln.strip().split("\t", 1)[0].upper() for ln in fh if ln.strip()]

    gene_ids: list[str] = []
    symbol_map: dict[str, str] = {}
    with _open_text(feat_path) as fh:
        for ln in fh:
            if not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            gid = parts[0].strip().split(".")[0]
            sym = parts[1].strip() if len(parts) > 1 and parts[1].strip() else gid
            gene_ids.append(gid)
            symbol_map[gid] = sym

    bc_to_col_idx = {bc: idx for idx, bc in enumerate(barcodes)}
    missing_bcs = [bc for bc in barcode_to_sample if bc not in bc_to_col_idx]
    if missing_bcs:
        raise ValueError(
            f"Declared barcodes missing from STARsolo output: {missing_bcs}"
        )

    ordered_items = list(barcode_to_sample.items())
    sample_names = [sample_id for _, sample_id in ordered_items]
    target_col_map = {
        bc_to_col_idx[bc] + 1: sample_pos
        for sample_pos, (bc, _) in enumerate(ordered_items)
    }

    n_features = len(gene_ids)
    n_samples = len(sample_names)
    dense = [[0] * n_samples for _ in range(n_features)]

    with _open_text(mtx_path) as fh:
        header = fh.readline()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"Invalid MatrixMarket header in {mtx_path}")
        dims_seen = False
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith("%"):
                continue
            if not dims_seen:
                dims_seen = True
                continue
            row_s, col_s, val_s = s.split()[:3]
            col_1based = int(col_s)
            sample_pos = target_col_map.get(col_1based)
            if sample_pos is None:
                continue
            row_0based = int(row_s) - 1
            dense[row_0based][sample_pos] += int(round(float(val_s)))

    df = pd.DataFrame(dense, index=gene_ids, columns=sample_names)
    df.index.name = "gene_id"
    # Collapse duplicate stripped Ensembl IDs if present
    if df.index.has_duplicates:
        df = df.groupby(level=0, sort=False).sum()

    per_sample_totals = {
        col: int(df[col].sum()) for col in sample_names
    }
    zero_samples = [s for s, tot in per_sample_totals.items() if tot == 0]
    if zero_samples:
        warnings.append(
            f"Pooled plate sample(s) have 0 assigned UMI counts: {zero_samples}. "
            "Verify well barcode orientation and whitelist."
        )

    n_before = len(df)
    df = df[df.sum(axis=1) > 0]
    n_after = len(df)
    if n_before - n_after > 0:
        warnings.append(
            f"{n_before - n_after} genes with zero UMI counts across all "
            "barcoded samples removed."
        )

    out_matrix = output_dir / "counts_matrix.tsv"
    df.to_csv(out_matrix, sep="\t")

    sym_path = output_dir / "counts_with_symbols.tsv"
    df_sym = df.copy()
    df_sym.insert(
        0,
        "gene_symbol",
        [symbol_map.get(gid, gid) for gid in df_sym.index],
    )
    df_sym.to_csv(sym_path, sep="\t")
    return out_matrix, sym_path, len(df), per_sample_totals


def rna_starsolo_bulk(params: dict[str, Any]) -> dict[str, Any]:
    fastq_files = [str(p) for p in (params.get("fastq_files") or [])]
    genome_dir = Path(str(params["genome_dir"]))
    genome_fasta = str(params.get("genome_fasta") or "")
    gtf_file = str(params.get("gtf_file") or "")
    threads = int(params.get("threads", 8))
    preset_hint = params.get("preset_hint")
    warnings: list[str] = []

    manifest, manifest_errors = validate_pooled_plate_manifest(
        params.get("pooled_plate_manifest"),
        preset_hint=preset_hint,
    )
    if manifest_errors or manifest is None:
        return {
            "status": "error",
            "error_type": "InvalidPooledPlateManifest",
            "details": "; ".join(manifest_errors or ["Missing pooled_plate_manifest."]),
            "warnings": warnings,
        }

    r1_files, r2_files, fq_err = _resolve_pooled_fastqs(fastq_files, manifest)
    if fq_err:
        return {
            "status": "error",
            "error_type": "InvalidPooledPlateFastqs",
            "details": fq_err,
            "warnings": warnings,
        }
    for fq in [*r1_files, *r2_files]:
        if not Path(fq).exists():
            return {
                "status": "error",
                "error_type": "MissingFastqInput",
                "details": f"FASTQ input does not exist: {fq}",
                "warnings": warnings,
            }

    output_dir = Path(
        params.get("output_dir")
        or str(Path(r1_files[0]).parent.parent / "aria_processing" / "starsolo_bulk")
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Ensure STAR genome index exists and is current (A4) ────────────────
    sa_nbases = _sa_index_nbases(genome_fasta)
    index_manifest = str(genome_dir / _INDEX_MANIFEST)
    index_inputs = [("fasta", genome_fasta or None), ("gtf", gtf_file or None)]
    index_params = {"sa_index_nbases": sa_nbases}
    star_version = _star_version()
    need_build = not _index_exists(genome_dir)
    if not need_build and genome_fasta and gtf_file:
        idx_ok, idx_reason = stage_is_current(
            index_manifest,
            inputs=index_inputs,
            params=index_params,
            tool_version=star_version,
        )
        if not idx_ok:
            warnings.append(f"[A4] rebuilding STAR index (stale: {idx_reason})")
            need_build = True
    if need_build:
        if not (genome_fasta and gtf_file):
            return {
                "status": "error",
                "error_type": "MissingGenomeIndex",
                "details": (
                    f"STAR index not found at {genome_dir}; supply genome_fasta "
                    "and gtf_file or prebuild the STAR index."
                ),
            }
        build_res = _build_star_index(
            genome_dir,
            genome_fasta,
            gtf_file,
            threads,
            warnings,
            sa_index_nbases=sa_nbases,
        )
        if build_res != "ok":
            return {
                "status": "error",
                "error_type": "IndexBuildFailed",
                "details": build_res,
                "warnings": warnings,
            }
        write_stage_manifest(
            index_manifest,
            build_stage_manifest(
                stage="star_index",
                inputs=index_inputs,
                params=index_params,
                tool_version=star_version,
            ),
        )

    barcode_to_sample: dict[str, str] = manifest["barcode_to_sample"]
    sample_names = list(barcode_to_sample.values())

    # Materialize or validate barcode whitelist
    if manifest.get("barcode_whitelist"):
        wl_path = Path(manifest["barcode_whitelist"])
        if not wl_path.is_file():
            return {
                "status": "error",
                "error_type": "MissingBarcodeWhitelist",
                "details": f"barcode_whitelist does not exist: {wl_path}",
            }
        wl_set = {
            ln.strip().split()[0].upper()
            for ln in wl_path.read_text(encoding="utf-8").splitlines()
            if ln.strip()
        }
        missing_wl = [bc for bc in barcode_to_sample if bc not in wl_set]
        if missing_wl:
            return {
                "status": "error",
                "error_type": "BarcodeWhitelistMismatch",
                "details": (
                    "Barcodes in barcode_to_sample are absent from "
                    f"barcode_whitelist: {missing_wl}"
                ),
            }
    else:
        wl_path = output_dir / "well_barcode_whitelist.txt"
        wl_path.write_text(
            "\n".join(sorted(barcode_to_sample.keys())) + "\n",
            encoding="utf-8",
        )

    stage_manifest_path = output_dir / "starsolo_bulk.stage.json"
    metrics_cache_path = output_dir / "starsolo_bulk.metrics.json"
    matrix_file = output_dir / "counts_matrix.tsv"
    symbols_file = output_dir / "counts_with_symbols.tsv"
    index_ref = (
        index_manifest if os.path.exists(index_manifest) else str(genome_dir / "SAindex")
    )
    stage_inputs = [
        *[(f"r1_{i}", p) for i, p in enumerate(r1_files, 1)],
        *[(f"r2_{i}", p) for i, p in enumerate(r2_files, 1)],
        ("whitelist", str(wl_path)),
        ("gtf", gtf_file or None),
        ("index", index_ref),
    ]
    stage_params = {
        "cb_start": manifest["cb_start"],
        "cb_len": manifest["cb_len"],
        "umi_start": manifest["umi_start"],
        "umi_len": manifest["umi_len"],
        "cdna_read": manifest["cdna_read"],
        "strand": manifest["strand"],
        "barcode_to_sample": barcode_to_sample,
    }

    current, _ = stage_is_current(
        str(stage_manifest_path),
        inputs=stage_inputs,
        params=stage_params,
        tool_version=star_version,
    )
    if current and matrix_file.exists() and metrics_cache_path.exists():
        cached = json.loads(metrics_cache_path.read_text(encoding="utf-8"))
        warnings.append(
            f"[resume] STARsolo pooled quantification skipped "
            f"({cached.get('n_genes', 0)} genes × {len(sample_names)} samples)"
        )
        cached["warnings"] = warnings
        cached["resumed"] = True
        return cached

    # STARsolo requires cDNA read(s) FIRST and Barcode+UMI read(s) SECOND.
    if manifest["cdna_read"] == "R2":
        cdna_arg = ",".join(r2_files)
        bc_arg = ",".join(r1_files)
    else:
        cdna_arg = ",".join(r1_files)
        bc_arg = ",".join(r2_files)

    is_gz = str(r1_files[0]).lower().endswith(".gz")
    solo_prefix = str(output_dir / "starsolo_")
    cmd = [
        "STAR",
        "--runThreadN", str(threads),
        "--genomeDir", str(genome_dir),
        "--readFilesIn", cdna_arg, bc_arg,
    ]
    if is_gz:
        cmd += ["--readFilesCommand", "zcat"]
    if gtf_file and Path(gtf_file).is_file():
        cmd += ["--sjdbGTFfile", gtf_file]
    cmd += [
        "--outFileNamePrefix", solo_prefix,
        "--outSAMtype", "None",
        "--soloType", "CB_UMI_Simple",
        "--soloCBwhitelist", str(wl_path),
        "--soloCBstart", str(manifest["cb_start"]),
        "--soloCBlen", str(manifest["cb_len"]),
        "--soloUMIstart", str(manifest["umi_start"]),
        "--soloUMIlen", str(manifest["umi_len"]),
        "--soloBarcodeReadLength", "0",
        "--soloStrand", str(manifest["strand"]),
        "--soloFeatures", "Gene",
        "--soloUMIdedup", "1MM_Directional_UMItools",
    ]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except FileNotFoundError:
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": "STAR is required in aria-rnaseq-env for STARsolo quantification.",
        }
    except subprocess.TimeoutExpired:
        return {
            "status": "error",
            "error_type": "STARsoloTimeout",
            "details": "STARsolo timed out (>2h).",
        }

    if proc.returncode != 0:
        return {
            "status": "error",
            "error_type": "STARsoloFailed",
            "details": (proc.stderr or proc.stdout or "STARsolo failed")[-400:],
            "warnings": warnings,
        }

    raw_gene_dir = Path(f"{solo_prefix}Solo.out") / "Gene" / "raw"
    try:
        out_matrix, out_sym, n_genes, per_sample_totals = _parse_starsolo_raw_matrix(
            raw_gene_dir,
            barcode_to_sample=barcode_to_sample,
            output_dir=output_dir,
            warnings=warnings,
        )
    except Exception as exc:
        return {
            "status": "error",
            "error_type": "STARsoloMatrixParseFailed",
            "details": str(exc),
            "warnings": warnings,
        }

    if n_genes == 0:
        return {
            "status": "error",
            "error_type": "EmptySTARsoloMatrix",
            "details": (
                "STARsolo produced zero non-zero genes across all declared "
                "well barcodes. Check barcode coordinates, strand, and genome."
            ),
            "warnings": warnings,
        }

    star_stats = _parse_star_log(Path(f"{solo_prefix}Log.final.out"))
    strand_used = 1 if manifest["strand"] == "Forward" else 2
    payload = {
        "status": "success",
        "counts_matrix": str(out_matrix),
        "counts_with_symbols": str(out_sym) if out_sym and out_sym.exists() else None,
        "n_genes": n_genes,
        "n_samples": len(sample_names),
        "sample_names": sample_names,
        "strand_used": strand_used,
        "read_layout": "pooled_plate_barcoded",
        "library_protocol": "pooled_plate_bulk",
        "umi_deduplicated": True,
        "pct_unique": star_stats.get("pct_unique", 0.0),
        "n_reads": star_stats.get("n_reads", 0),
        "per_sample_umi_totals": per_sample_totals,
        "output_dir": str(output_dir),
        "warnings": warnings,
    }
    metrics_cache_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_stage_manifest(
        str(stage_manifest_path),
        build_stage_manifest(
            stage="starsolo_bulk",
            inputs=stage_inputs,
            params=stage_params,
            tool_version=star_version,
        ),
    )
    return payload


if __name__ == "__main__":
    run_script(rna_starsolo_bulk)
