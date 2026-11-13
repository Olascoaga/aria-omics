"""
ARIA Paired Candidate eRegulons Script (V48 Slice I6 — Tier 2b Regulatory Layer)
--------------------------------------------------------------------------------
Assembles candidate enhancer-driven regulons (`TF -> cis peak -> target gene`)
from paired same-cell RNA + ATAC data (`.h5mu` preferred via `read_h5mu_paired`,
or deterministically barcode-aligned `.h5ad` / 10x raw count matrices). Executed
inside `aria-integration-env` by `EnvironmentManager`.

Governance & Statistical Contract (`ADR-002`, `ADR-005`, `ADR-011`, `ADR-064`):
  1. Offline, checksum-verified reference assets only (`load_local_motif_collection`,
     local GTF, and local `.fai`-indexed genome FASTA). Missing assets return
     structured errors/skips and NEVER fabricate regulons or make network calls.
  2. Reuses `aria.utils.peak2gene_core` (`ADR-068`) + `aria.utils.eregulon_core`:
     fixed-width `500 bp` peak-center motif scanning, GC-stratified CMH cistrome
     enrichment, `(GC, accessibility)` trans-chromosomal peak-to-gene FDR, and
     `(mean_expr, cells_detected)` 2D bin-stratified module FDR (`bh_correct`
     applied across all tested `(motif_id, tf, sign_class)` modules prior to
     best-motif selection per `(tf, sign_class)`).
  3. Strictly associative terminology: sign classes (`pos_tf_pos_p2g`,
     `neg_tf_pos_p2g`, etc.), paralog motif-sharing ambiguity flags, and AUCell
     activity scores describe paired multi-omic co-variation only and NEVER
     establish direct transcriptional activation, repression, or causality
     without orthogonal perturbation or chromatin-conformation validation.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
import sys
import tempfile
from typing import Any
import uuid

import numpy as np
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.utils.mudata_io import (
    _stage_cache,
    _stage_complete,
    paired_input_files,
    MultiomeInputError,
)
from aria.scripts._base import run_script
from aria.utils.eregulon_core import (
    ContigMismatchError,
    compute_candidate_eregulons,
)
from aria.utils.motifs import DEFAULT_COLLECTION, load_local_motif_collection
from aria.utils.mudata_io import _select_raw_counts_matrix, read_h5mu_paired
from aria.utils.peak2gene_core import (
    parse_gtf_tss_detailed,
    resolve_gtf_path,
    resolve_indexed_fasta_path,
)
from aria.utils.reference_integrity import (
    public_integrity_result,
    reference_is_usable,
    verify_reference_file,
)
from aria.utils.wnn_diagnostics import sanitize_for_hdf5_uns

_EREGULON_CSV_COLUMNS = (
    "regulon_id",
    "tf",
    "sign_class",
    "best_motif_id",
    "motif_source",
    "tf_symbol_match_mode",
    "n_motifs_tested_for_tf",
    "n_candidate_genes",
    "n_target_genes",
    "n_enhancer_peaks",
    "n_triplets",
    "module_zscore",
    "module_pvalue",
    "regulon_padj",
    "regulon_padj_scope",
    "gsea_leading_edge_es",
    "mean_r_tf_gene",
    "mean_r_tf_gene_partial_depth",
    "depth_confounded",
    "cistrome_cmh_z",
    "cistrome_cmh_padj",
    "cistrome_log2_fe",
    "motif_family_ambiguity",
    "ambiguous_paralog_tfs",
    "mean_rna_aucell",
    "mean_atac_aucell",
    "rna_atac_activity_corr",
    "target_genes",
)

_TRIPLET_CSV_COLUMNS = (
    "regulon_id",
    "tf",
    "motif_id",
    "motif_source",
    "sign_class",
    "peak",
    "gene",
    "r_tf_gene",
    "r_tf_gene_partial_depth",
    "r_tf_peak",
    "r_tf_peak_partial_depth",
    "r_peak_gene",
    "padj_peak_gene",
    "distance_bp",
    "region_class",
    "regulon_padj",
)

_ASSOCIATIVE_CAVEAT = (
    "Candidate eRegulons combine PWM sequence matches, cis peak-to-gene "
    "co-variation, and TF-to-target expression correlation across paired cells; "
    "they represent associative regulatory candidates (subject to paralog motif "
    "sharing and indirect co-expression) and do not establish causal regulation "
    "without independent perturbation or binding assays."
)


def _write_csv_atomic(
    path: Path,
    columns: tuple[str, ...],
    rows: list[dict[str, Any]],
) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    os.close(fd)
    tmp_p = Path(tmp_name)
    try:
        with open(tmp_p, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=list(columns), extrasaction="ignore"
            )
            writer.writeheader()
            for r in rows:
                writer.writerow(r)
        os.replace(str(tmp_p), str(path))
    finally:
        if tmp_p.exists():
            try:
                tmp_p.unlink()
            except OSError:
                pass
    return str(path)


def _format_eregulon_csv_rows(eregulons: list[dict[str, Any]]) -> list[dict[str, Any]]:
    formatted: list[dict[str, Any]] = []
    for reg in eregulons:
        row = dict(reg)
        amb = reg.get("ambiguous_paralog_tfs") or []
        row["ambiguous_paralog_tfs"] = ";".join(
            f"{a['tf']}({a['peak_jaccard']:.2f})" for a in amb
        )
        row["target_genes"] = ";".join(str(g) for g in (reg.get("target_genes") or []))
        formatted.append(row)
    return formatted


def _load_precomputed_p2g_csv(csv_path: Path) -> dict[str, Any]:
    links: list[dict[str, Any]] = []
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            links.append(
                {
                    "gene": str(r["gene"]),
                    "peak": str(r["peak"]),
                    "correlation": float(r["correlation"]),
                    "padj": (
                        float(r["padj"])
                        if r.get("padj") not in (None, "")
                        else None
                    ),
                    "distance_bp": int(r["distance_bp"]),
                    "region_class": str(r.get("region_class") or "cis"),
                }
            )
    return {
        "links": links,
        "n_tested": len(links),
        "n_fdr_significant": len(links),
        "background_match": "precomputed_csv",
    }


def _load_rna_adata(rna_path: str):
    import scanpy as sc

    p = Path(rna_path)
    if p.suffix == ".h5ad":
        return sc.read_h5ad(str(p))
    if p.suffix == ".h5":
        return sc.read_10x_h5(str(p))
    if p.is_dir():
        return sc.read_10x_mtx(str(p), var_names="gene_symbols", cache=False)
    return None


def _load_atac_adata(atac_files: list[str]):
    import anndata as ad
    import scanpy as sc

    for f in atac_files:
        p = Path(f)
        if p.suffix == ".h5ad":
            return ad.read_h5ad(str(p)), None
        if p.suffix == ".h5":
            return sc.read_10x_h5(str(p), gex_only=False), None
        if p.is_dir():
            return sc.read_10x_mtx(str(p), cache=False), None
        low = p.name.lower()
        if low.endswith((".tsv.gz", ".fragments.tsv.gz", ".bed.gz", ".bam")):
            return None, (
                f"Raw chromatin fragment/alignment input ({p.name}) must be "
                "quantified into a peak-by-cell count matrix (.h5mu or .h5ad) "
                "before candidate eRegulon assembly."
            )
    return None, "Candidate eRegulon assembly requires a valid ATAC peak-by-cell matrix."


def integration_eregulons(params: dict) -> dict[str, Any]:
    h5mu_raw = params.get("h5mu_path") or params.get("data_path")
    rna_files = list(params.get("rna_files") or [])
    atac_files = list(params.get("atac_files") or [])
    genome = str(params.get("genome") or "hg38")
    organism = str(params.get("organism") or "Homo sapiens")
    output_dir = Path(str(params.get("output_dir") or "/tmp/aria_eregulons"))
    seed = int(params.get("seed", 0) or 0)

    if params.get("distance_bp") is not None and params.get("distance_bp") != "":
        distance_bp = int(params["distance_bp"])
        distance_kb = round(float(distance_bp) / 1000.0, 3)
    else:
        distance_kb_val = float(params.get("distance_kb", 500))
        distance_bp = int(round(distance_kb_val * 1000.0))
        distance_kb = (
            int(distance_kb_val)
            if float(distance_kb_val).is_integer()
            else round(distance_kb_val, 3)
        )

    min_distance_bp = int(params.get("min_distance_bp", 0) or 0)
    promoter_bp = int(params.get("promoter_bp", 2000) or 2000)
    p2g_min_corr = float(params.get("p2g_min_corr", 0.05))
    p2g_max_padj = float(params.get("p2g_max_padj", 0.05))
    min_shared_cells = int(params.get("min_shared_cells", 50) or 50)
    min_cells_detected = int(params.get("min_cells_detected", 10) or 10)
    scan_window_bp = int(params.get("scan_window_bp", 500) or 500)
    motif_pvalue = float(params.get("motif_pvalue", 5e-5))
    scanner_backend = str(params.get("scanner_backend") or "auto")
    min_tf_gene_corr = float(params.get("min_tf_gene_corr", 0.05))
    min_target_genes = int(params.get("min_target_genes", 5) or 5)
    max_regulon_padj = float(params.get("max_regulon_padj", 0.05))
    n_null_permutations = int(params.get("n_null_permutations", 500))
    aucell_top_frac = float(params.get("aucell_top_frac", 0.05))
    compute_cell_activity = bool(params.get("compute_cell_activity", True))
    standard_chroms_only = bool(params.get("standard_chroms_only", True))
    sign_classes = tuple(
        params.get("sign_classes") or ("pos_tf_pos_p2g", "neg_tf_pos_p2g")
    )

    warnings: list[str] = []
    output_dir.mkdir(parents=True, exist_ok=True)

    if not (0.0 < motif_pvalue < 1.0):
        return {
            "status": "error",
            "error_type": "InvalidParameter",
            "details": f"motif_pvalue must lie strictly in (0, 1), got {motif_pvalue}",
        }

    # Auto-detect paired .h5mu in rna_files/atac_files if h5mu_path was omitted
    resolved_h5mu: str | None = str(h5mu_raw).strip() if h5mu_raw else None
    if not resolved_h5mu:
        for cand in [*rna_files, *atac_files]:
            if str(cand).lower().endswith(".h5mu") and Path(str(cand)).exists():
                resolved_h5mu = str(cand)
                break

    if resolved_h5mu and not Path(resolved_h5mu).exists():
        return {
            "status": "error",
            "error_type": "MissingH5mu",
            "details": f"Paired .h5mu file does not exist: {resolved_h5mu}",
        }

    if not resolved_h5mu:
        valid_rna = [str(f) for f in rna_files if Path(str(f)).exists()]
        valid_atac = [str(f) for f in atac_files if Path(str(f)).exists()]
        if not valid_rna:
            return {
                "status": "error",
                "error_type": "MissingRNA",
                "details": "No valid RNA files or paired .h5mu provided for eRegulons.",
            }
        if not valid_atac:
            return {
                "status": "error",
                "error_type": "MissingATAC",
                "details": "No valid ATAC files provided for eRegulons.",
            }
    else:
        valid_rna = []
        valid_atac = []

    # ── Resolve offline MEME motif collection (W-PRIV / ADR-002) ──
    motif_file = params.get("motif_file")
    motif_collection = str(params.get("motif_collection") or DEFAULT_COLLECTION)
    if motif_file:
        mf_path = Path(str(motif_file))
        if not mf_path.is_file():
            return {
                "status": "error",
                "error_type": "MissingMotifCollection",
                "details": f"motif_file not found: {motif_file}",
            }
        integrity = verify_reference_file(str(mf_path))
        if not reference_is_usable(integrity):
            return {
                "status": "error",
                "error_type": "InvalidMotifChecksum",
                "details": "motif_file checksum validation failed.",
                "reference_integrity": public_integrity_result(integrity),
            }
        meme_path = str(mf_path)
        motif_source_meta: dict[str, Any] = {
            "collection": mf_path.stem,
            "source": "explicit_file",
            "release": "custom",
            "integrity": public_integrity_result(integrity),
        }
    else:
        loaded = load_local_motif_collection(motif_collection)
        if loaded is None:
            return {
                "status": "error",
                "error_type": "MissingMotifCollection",
                "details": (
                    f"Motif collection '{motif_collection}' is not staged locally "
                    "under ARIA_MOTIF_DIR; stage it via scripts/fetch_motifs.py "
                    "(offline-only, never fabricated)."
                ),
            }
        meme_path, motif_source_meta = loaded

    # ── Resolve local GTF and `.fai`-indexed genome FASTA ──
    explicit_gtf = params.get("gtf_path") or params.get("gene_annotation_gtf")
    try:
        gtf_resolved = resolve_gtf_path(
            str(explicit_gtf) if explicit_gtf else None,
            genome=genome,
        )
    except FileNotFoundError as exc:
        return {
            "status": "error",
            "error_type": "MissingGTF",
            "details": str(exc),
        }
    if not gtf_resolved:
        return {
            "status": "error",
            "error_type": "MissingGTF",
            "details": (
                f"Candidate eRegulon assembly requires a local GTF annotation "
                f"for genome='{genome}'."
            ),
        }

    explicit_fasta = params.get("genome_fasta") or params.get("fasta_path")
    if explicit_fasta and not Path(str(explicit_fasta)).is_file():
        return {
            "status": "error",
            "error_type": "MissingIndexedFASTA",
            "details": f"genome_fasta not found: {explicit_fasta}",
        }
    fasta_resolved = resolve_indexed_fasta_path(
        str(explicit_fasta) if explicit_fasta else None,
        genome=genome,
    )
    if not fasta_resolved:
        return {
            "status": "error",
            "error_type": "MissingIndexedFASTA",
            "details": (
                f"Candidate eRegulon motif scanning requires a local `.fai`-indexed "
                f"reference FASTA for genome='{genome}'."
            ),
        }

    try:
        stage_inputs = (
            paired_input_files("h5mu_input", resolved_h5mu)
            if resolved_h5mu
            else [
                item
                for role, paths in (("rna", valid_rna), ("atac", valid_atac))
                for i, path in enumerate(paths)
                for item in paired_input_files(f"{role}:{i}", path)
            ]
        )
    except MultiomeInputError as exc:
        return {"status": "error", "error_type": exc.error_type, "details": str(exc)}
    for role, path in (("gtf", gtf_resolved), ("genome_fasta", fasta_resolved),
                       ("fasta_index", str(fasta_resolved) + ".fai" if fasta_resolved else None),
                       ("peaks", params.get("peaks_path")), ("motifs", meme_path), ("precomputed_p2g", params.get("precomputed_p2g_csv"))):
        if path and Path(str(path)).is_file():
            stage_inputs.append((role, str(path)))
    stage_params = dict(params)
    stage_params.update(
        distance_bp=distance_bp,
        min_distance_bp=min_distance_bp,
        promoter_bp=promoter_bp,
        p2g_min_corr=p2g_min_corr,
        p2g_max_padj=p2g_max_padj,
        min_shared_cells=min_shared_cells,
        min_cells_detected=min_cells_detected,
        scan_window_bp=scan_window_bp,
        motif_pvalue=motif_pvalue,
        scanner_backend=scanner_backend,
        min_tf_gene_corr=min_tf_gene_corr,
        min_target_genes=min_target_genes,
        max_regulon_padj=max_regulon_padj,
        n_null_permutations=n_null_permutations,
        aucell_top_frac=aucell_top_frac,
        compute_cell_activity=compute_cell_activity,
        standard_chroms_only=standard_chroms_only,
        sign_classes=sign_classes,
        seed=seed,
        genome=genome,
        organism=organism,
    )
    cached, stage_state = _stage_cache(
        "integration_eregulons",
        output_dir,
        stage_inputs,
        stage_params,
        params.get("resume", True),
        required=("output_path", "triplets_path"),
    )
    if cached:
        return cached

    # ── Ingest paired RNA + ATAC raw count matrices in canonical barcode order ──
    mdata = None
    if resolved_h5mu:
        paired = read_h5mu_paired(resolved_h5mu, require_raw_counts=True)
        if paired.get("status") != "success":
            return {
                "status": "error",
                "error_type": str(paired.get("error_type") or "InvalidPairedH5mu"),
                "details": str(
                    paired.get("details") or "Paired .h5mu validation failed."
                ),
                "warnings": warnings,
            }
        mdata = paired["mdata"]
        rna_adata = paired["rna_adata"]
        atac_adata = paired["atac_adata"]
        barcodes = list(paired["shared_barcodes"])
        if len(barcodes) < min_shared_cells:
            return {
                "status": "error",
                "error_type": "InsufficientSharedCells",
                "details": (
                    f"Only {len(barcodes)} shared cells in paired .h5mu "
                    f"(minimum required: {min_shared_cells})."
                ),
            }
        r_idx_map = {str(b): i for i, b in enumerate(rna_adata.obs_names)}
        a_idx_map = {str(b): i for i, b in enumerate(atac_adata.obs_names)}
        r_order = [r_idx_map[b] for b in barcodes]
        a_order = [a_idx_map[b] for b in barcodes]
        rna_mat = paired["rna_counts_matrix"][r_order, :]
        atac_mat = paired["atac_counts_matrix"][a_order, :]
        rna_gene_names = [str(g) for g in rna_adata.var_names]
        atac_peak_names = [str(p) for p in atac_adata.var_names]
        input_mode = "paired_h5mu"
    else:
        rna_adata = _load_rna_adata(valid_rna[0])
        if rna_adata is None:
            return {
                "status": "error",
                "error_type": "InvalidRNAInput",
                "details": f"Could not load RNA matrix from {valid_rna[0]}",
            }
        atac_adata, atac_err = _load_atac_adata(valid_atac)
        if atac_adata is None:
            return {
                "status": "error",
                "error_type": "ATACMatrixNotLoadable",
                "details": str(atac_err),
            }
        r_idx_map = {str(b): i for i, b in enumerate(rna_adata.obs_names)}
        a_idx_map = {str(b): i for i, b in enumerate(atac_adata.obs_names)}
        atac_set = set(a_idx_map)
        seen_bc: set[str] = set()
        barcodes = []
        for bc in (str(x) for x in rna_adata.obs_names):
            if bc in atac_set and bc not in seen_bc:
                barcodes.append(bc)
                seen_bc.add(bc)
        if len(barcodes) < min_shared_cells:
            return {
                "status": "error",
                "error_type": "InsufficientSharedCells",
                "details": (
                    f"Only {len(barcodes)} shared barcodes between RNA and ATAC "
                    f"(minimum required: {min_shared_cells})."
                ),
            }
        rna_sub = rna_adata[[r_idx_map[b] for b in barcodes], :].copy()
        atac_sub = atac_adata[[a_idx_map[b] for b in barcodes], :].copy()
        rna_chk = _select_raw_counts_matrix(
            rna_sub, modality_label="rna", require_raw_counts=True
        )
        if rna_chk.get("matrix") is None or not rna_chk.get("inspection", {}).get("is_raw_counts", False):
            return {
                "status": "error",
                "error_type": "NonRawCountsRNA",
                "details": str(
                    rna_chk.get("details")
                    or f"Unverified RNA count matrix ({rna_chk.get('source')})"
                ),
            }
        atac_chk = _select_raw_counts_matrix(
            atac_sub, modality_label="atac", require_raw_counts=True
        )
        if atac_chk.get("matrix") is None or not atac_chk.get("inspection", {}).get("is_raw_counts", False):
            return {
                "status": "error",
                "error_type": "NonRawCountsATAC",
                "details": str(
                    atac_chk.get("details")
                    or f"Unverified ATAC count matrix ({atac_chk.get('source')})"
                ),
            }
        rna_mat = rna_chk["matrix"]
        atac_mat = atac_chk["matrix"]
        rna_gene_names = [str(g) for g in rna_sub.var_names]
        atac_peak_names = [str(p) for p in atac_sub.var_names]
        input_mode = "paired_matrices"

    gene_tss_records, gtf_stats = parse_gtf_tss_detailed(
        gtf_resolved, wanted=set(rna_gene_names)
    )
    if not gene_tss_records:
        return {
            "status": "error",
            "error_type": "EmptyGTFAnnotation",
            "details": f"No gene TSS records parsed from GTF: {gtf_resolved}",
            "gtf_stats": gtf_stats,
        }

    precomputed_p2g: dict[str, Any] | None = None
    p2g_csv = params.get("precomputed_p2g_csv")
    if p2g_csv and Path(str(p2g_csv)).is_file():
        precomputed_p2g = _load_precomputed_p2g_csv(Path(str(p2g_csv)))

    try:
        core_res = compute_candidate_eregulons(
            rna_matrix=rna_mat,
            rna_gene_names=rna_gene_names,
            atac_matrix=atac_mat,
            atac_peak_names=atac_peak_names,
            gene_tss_records=gene_tss_records,
            meme_path=meme_path,
            fasta_path=fasta_resolved,
            distance_bp=distance_bp,
            min_distance_bp=min_distance_bp,
            promoter_bp=promoter_bp,
            p2g_min_corr=p2g_min_corr,
            p2g_max_padj=p2g_max_padj,
            min_cells_detected=min_cells_detected,
            scan_window_bp=scan_window_bp,
            motif_pvalue=motif_pvalue,
            scanner_backend=scanner_backend,
            min_tf_gene_corr=min_tf_gene_corr,
            min_target_genes=min_target_genes,
            max_regulon_padj=max_regulon_padj,
            n_null_permutations=n_null_permutations,
            sign_classes=sign_classes,
            compute_cell_activity=compute_cell_activity,
            aucell_top_frac=aucell_top_frac,
            standard_chroms_only=standard_chroms_only,
            seed=seed,
            precomputed_p2g=precomputed_p2g,
        )
    except ContigMismatchError as exc:
        return {
            "status": "error",
            "error_type": "ContigMismatch",
            "details": str(exc),
        }

    eregulons: list[dict[str, Any]] = core_res["eregulons"]
    triplets: list[dict[str, Any]] = core_res["triplets"]

    ereg_csv = _write_csv_atomic(
        output_dir / "candidate_eregulons.csv",
        _EREGULON_CSV_COLUMNS,
        _format_eregulon_csv_rows(eregulons),
    )
    trip_csv = _write_csv_atomic(
        output_dir / "eregulon_triplets.csv",
        _TRIPLET_CSV_COLUMNS,
        triplets,
    )

    n_pos_pos = sum(1 for r in eregulons if r["sign_class"] == "pos_tf_pos_p2g")
    n_neg_pos = sum(1 for r in eregulons if r["sign_class"] == "neg_tf_pos_p2g")
    n_other_sign = len(eregulons) - n_pos_pos - n_neg_pos
    n_paralog_ambig = sum(1 for r in eregulons if r["motif_family_ambiguity"])
    n_depth_confounded = sum(1 for r in eregulons if r["depth_confounded"])

    if n_paralog_ambig > 0:
        warnings.append(
            f"{n_paralog_ambig} candidate eRegulon(s) share cis-enhancer peaks/motifs "
            "with paralogous TFs (`motif_family_ambiguity=True`); inspect "
            "`ambiguous_paralog_tfs` before attributing TF specificity."
        )
    if n_depth_confounded > 0:
        warnings.append(
            f"{n_depth_confounded} candidate eRegulon(s) showed substantial attenuation "
            "after conditioning on library sequencing depth (`depth_confounded=True`)."
        )

    # Optional atomic persistence to `.h5mu` (`update_h5mu`, `write_h5mu`, or `output_h5mu`)
    h5mu_output_path: str | None = None
    update_in_place = bool(params.get("update_h5mu", False))
    write_h5mu_flag = (
        update_in_place
        or bool(params.get("write_h5mu", False))
        or bool(params.get("output_h5mu"))
    )
    if write_h5mu_flag and mdata is not None:
        default_target = (
            Path(str(resolved_h5mu))
            if (update_in_place and resolved_h5mu)
            else (output_dir / "eregulons_paired.h5mu")
        )
        target_h5mu = Path(str(params.get("output_h5mu") or default_target))
        target_h5mu.parent.mkdir(parents=True, exist_ok=True)

        # Align AUCell matrices `(len(barcodes), n_eregulons)` to `mdata.obs_names`
        rna_auc = core_res["rna_aucell"]
        atac_auc = core_res["atac_aucell"]
        if rna_auc.shape[1] > 0:
            bc_to_row = {b: i for i, b in enumerate(barcodes)}
            m_obs = [str(b) for b in mdata.obs_names]
            full_rna_auc = np.full(
                (mdata.n_obs, rna_auc.shape[1]), np.nan, dtype=np.float32
            )
            full_atac_auc = np.full(
                (mdata.n_obs, atac_auc.shape[1]), np.nan, dtype=np.float32
            )
            for m_i, b in enumerate(m_obs):
                src_i = bc_to_row.get(b)
                if src_i is not None:
                    full_rna_auc[m_i] = rna_auc[src_i]
                    full_atac_auc[m_i] = atac_auc[src_i]
            mdata.obsm["X_eregulon_rna"] = full_rna_auc
            mdata.obsm["X_eregulon_atac"] = full_atac_auc

        mdata.uns["aria_eregulons"] = sanitize_for_hdf5_uns(
            {
                "slice": "V48_I6_candidate_eregulons",
                "seed": seed,
                "distance_bp": distance_bp,
                "scan_window_bp": core_res["scan_window_bp"],
                "motif_pvalue": motif_pvalue,
                "scanner_backend": core_res["scanner_backend"],
                "regulon_ids": [r["regulon_id"] for r in eregulons],
                "n_eregulons": len(eregulons),
                "n_triplets": len(triplets),
                "n_unique_tfs": core_res["n_unique_tfs"],
                "n_paralog_ambiguous": n_paralog_ambig,
                "n_depth_confounded": n_depth_confounded,
                "associative_only": True,
                "top_eregulons": eregulons[:25],
            }
        )
        tmp_h5mu = target_h5mu.parent / (
            f".{target_h5mu.stem}.tmp.{os.getpid()}.{uuid.uuid4().hex[:8]}.h5mu"
        )
        try:
            mdata.write_h5mu(str(tmp_h5mu))
            os.replace(str(tmp_h5mu), str(target_h5mu))
            h5mu_output_path = str(target_h5mu)
        finally:
            if tmp_h5mu.exists():
                try:
                    tmp_h5mu.unlink()
                except OSError:
                    pass

    result = {
        "status": "success",
        "validation_level": "production",
        "input_mode": input_mode,
        "organism": organism,
        "genome": genome,
        "seed": seed,
        "n_shared_cells": len(barcodes),
        "distance_bp": distance_bp,
        "distance_kb": distance_kb,
        "scan_window_bp": core_res["scan_window_bp"],
        "motif_pvalue": motif_pvalue,
        "scanner_backend": core_res["scanner_backend"],
        "motif_source": {
            "collection": motif_source_meta.get("collection"),
            "release": motif_source_meta.get("release"),
            "source": motif_source_meta.get("source"),
            "n_motifs": core_res["n_motifs_total"],
            "integrity": motif_source_meta.get("integrity"),
        },
        "p2g_summary": core_res["p2g_summary"],
        "n_motif_tf_mappings": core_res["n_motif_tf_mappings"],
        "n_unmapped_motifs": core_res["n_unmapped_motifs"],
        "unmapped_motifs_sample": core_res["unmapped_motifs_sample"],
        "n_candidate_modules_tested": core_res["n_candidate_modules_tested"],
        "n_eregulons": len(eregulons),
        "n_triplets": len(triplets),
        "n_unique_tfs": core_res["n_unique_tfs"],
        "n_unique_enhancers": core_res["n_unique_enhancers"],
        "n_unique_target_genes": core_res["n_unique_target_genes"],
        "n_pos_tf_pos_p2g": n_pos_pos,
        "n_neg_tf_pos_p2g": n_neg_pos,
        "n_other_sign_classes": n_other_sign,
        "n_paralog_ambiguous": n_paralog_ambig,
        "n_depth_confounded": n_depth_confounded,
        "top_eregulons": eregulons[:30],
        "output_path": ereg_csv,
        "triplets_path": trip_csv,
        "h5mu_output_path": h5mu_output_path,
        "associative_only": True,
        "perturbation_validated": False,
        "chromatin_conformation_validated": False,
        "caveat": _ASSOCIATIVE_CAVEAT,
        "warnings": warnings,
    }
    return _stage_complete(result, stage_state)


if __name__ == "__main__":
    run_script(integration_eregulons)
