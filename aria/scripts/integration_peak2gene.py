"""
ARIA Paired Peak-to-Gene Regulatory Linking Script (V48 Slice I5 — Tier 2b)
---------------------------------------------------------------------------
Computes cis peak-to-gene accessibility-expression associations across same-cell
paired multiome data (`.h5mu` preferred via `read_h5mu_paired`, or deterministically
barcode-aligned RNA + ATAC `.h5ad` / 10x matrices). Executed inside
`aria-integration-env` by `EnvironmentManager`.

Statistical & Governance Contract (ADR-002, ADR-005, ADR-011, ADR-050, ADR-064):
  1. Reuses the unified `aria.utils.peak2gene_core` engine shared with
     `aria/scripts/chromatin_regulatory.py` (`_peak_to_gene`), eliminating
     duplicate statistical implementations across lanes.
  2. Verified raw-count ingestion (`read_h5mu_paired` / `_select_raw_counts_matrix`)
     and canonical barcode ordering (never unordered `set` iteration, never blind
     double `log1p` normalization, never arbitrary 5,000-gene truncation).
  3. `Signac::LinkPeaks`-faithful normalization (`log1p(CP10k)` RNA vs raw-count
     ATAC peaks, or optional `atac_transform="tfidf"`), feature prevalence gate
     (`min_cells_detected=10`), and deterministic `(GC, accessibility)` or
     `(accessibility, length)` binned trans-chromosomal background null (`zscore`,
     two-sided `p_value`, one-sided `p_value_signac`, and shared
     `aria.utils.stats.bh_correct` `padj` across all tested cis pairs prior to
     effect-size thresholding).
  4. Strictly associative interpretation boundary: correlation signs, z-scores,
     and TSS distances describe paired cis co-variation only and NEVER establish
     enhancer, silencer, or poised regulatory mechanism without independent
     orthogonal validation (e.g. CRISPRi/a, reporter assay, Hi-C/ChIA-PET).
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
import sys
import tempfile
from typing import Any
import uuid

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script
from aria.utils.mudata_io import _select_raw_counts_matrix, read_h5mu_paired
from aria.utils.peak2gene_core import (
    compute_peak_to_gene_links,
    parse_gtf_tss_detailed,
    parse_peak_coord,
    resolve_gtf_path,
    resolve_indexed_fasta_path,
)
from aria.utils.wnn_diagnostics import sanitize_for_hdf5_uns

_CSV_COLUMNS = (
    "gene",
    "gene_id",
    "peak",
    "correlation",
    "zscore",
    "p_value",
    "padj",
    "p_value_signac",
    "distance_bp",
    "distance_kb",
    "tss_offset_bp",
    "direction",
    "region_class",
    "background_scope",
)

_ASSOCIATIVE_CAVEAT = (
    "Peak-to-gene links measure paired cis accessibility-expression co-variation "
    "controlled against a matched trans-chromosomal background null; correlation "
    "sign and genomic distance are strictly associative and do not establish "
    "regulatory mechanism without independent orthogonal validation."
)


class MissingGTFError(Exception):
    """Raised when a local GTF gene annotation is required but unavailable."""


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
                "before peak-to-gene linking."
            )
    return None, "Peak-to-gene linking requires a valid ATAC peak-by-cell matrix."


def _maybe_map_peaks_from_bed(peak_names: list[str], peaks_path: str) -> list[str]:
    """Attach `chr:start-end` names from a BED file only if `peak_names` lack coordinates."""
    if not peaks_path or not Path(peaks_path).is_file():
        return peak_names
    n_parseable = sum(1 for p in peak_names[: min(len(peak_names), 100)] if parse_peak_coord(p) is not None)
    if n_parseable > 0:
        return peak_names

    bed_coords: list[str] = []
    by_name: dict[str, str] = {}
    with open(peaks_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                s, e = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            coord_str = f"{parts[0].strip()}:{s}-{e}"
            bed_coords.append(coord_str)
            if len(parts) >= 4 and parts[3].strip():
                by_name[parts[3].strip()] = coord_str

    if by_name and sum(1 for p in peak_names if p in by_name) >= max(1, len(peak_names) // 2):
        return [by_name.get(p, p) for p in peak_names]
    if len(bed_coords) == len(peak_names):
        return bed_coords
    return peak_names


def _write_links_csv_atomic(path: Path, links: list[dict[str, Any]]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    os.close(fd)
    tmp_p = Path(tmp_name)
    try:
        with open(tmp_p, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=list(_CSV_COLUMNS), extrasaction="ignore"
            )
            writer.writeheader()
            for row in links:
                writer.writerow(row)
        os.replace(str(tmp_p), str(path))
    finally:
        if tmp_p.exists():
            try:
                tmp_p.unlink()
            except OSError:
                pass
    return str(path)


def integration_peak2gene(params: dict) -> dict[str, Any]:
    h5mu_raw = params.get("h5mu_path") or params.get("data_path")
    rna_files = list(params.get("rna_files") or [])
    atac_files = list(params.get("atac_files") or [])
    peaks_path = str(params.get("peaks_path") or "")
    genome = str(params.get("genome") or "hg38")
    organism = str(params.get("organism") or "Homo sapiens")
    output_dir = Path(str(params.get("output_dir") or "/tmp/aria_p2g"))
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
    min_corr = float(params.get("min_corr", 0.3))
    min_shared_cells = int(params.get("min_shared_cells", 50) or 50)
    min_cells_detected = int(params.get("min_cells_detected", 10) or 10)
    atac_transform = str(params.get("atac_transform") or "raw_counts").strip().lower()
    standard_chroms_only = bool(params.get("standard_chroms_only", False))
    background_null = bool(params.get("background_null", True))
    use_fdr = bool(params.get("use_fdr", True))

    if "max_padj" in params:
        max_padj = (
            float(params["max_padj"]) if params["max_padj"] is not None else None
        )
    else:
        max_padj = 0.05 if (background_null and use_fdr) else None

    max_pvalue_signac = (
        float(params["max_pvalue_signac"])
        if params.get("max_pvalue_signac") is not None
        else None
    )

    if atac_transform not in ("raw_counts", "tfidf"):
        return {
            "status": "error",
            "error_type": "InvalidParameter",
            "details": (
                f"Unsupported atac_transform='{atac_transform}'. "
                "Valid options: 'raw_counts', 'tfidf'."
            ),
        }

    warnings: list[str] = []
    output_dir.mkdir(parents=True, exist_ok=True)

    # Auto-detect paired .h5mu passed through rna_files/atac_files (IntegrationAgent pattern).
    resolved_h5mu: str | None = str(h5mu_raw).strip() if h5mu_raw else None
    if not resolved_h5mu:
        for cand in [*rna_files, *atac_files]:
            if str(cand).lower().endswith(".h5mu") and Path(str(cand)).exists():
                resolved_h5mu = str(cand)
                break

    if resolved_h5mu and not Path(resolved_h5mu).is_file():
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
                "details": "No valid RNA files or paired .h5mu provided for peak-to-gene linking.",
            }
        if not valid_atac:
            return {
                "status": "error",
                "error_type": "ATACMatrixNotLoadable",
                "details": "No valid ATAC files provided for peak-to-gene linking.",
                "warnings": warnings,
            }
    else:
        valid_rna = []
        valid_atac = []

    # Resolve local GTF and optional indexed FASTA before loading heavy matrices.
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
                f"GTF gene annotation not found for genome '{genome}'. "
                "Provide gtf_path explicitly, set ARIA_GTF_PATH, or stage "
                f"~/.aria/genomes/{genome}/annotation.gtf."
            ),
        }

    explicit_fasta = params.get("genome_fasta")
    try:
        fasta_resolved = resolve_indexed_fasta_path(
            str(explicit_fasta) if explicit_fasta else None,
            genome=genome,
        )
    except FileNotFoundError as exc:
        return {
            "status": "error",
            "error_type": "MissingGenomeFasta",
            "details": str(exc),
        }

    try:
        import scipy.sparse as sp

        src_mdata = None
        if resolved_h5mu:
            paired = read_h5mu_paired(
                resolved_h5mu,
                require_raw_counts=False,
            )
            if paired.get("status") != "success":
                return {
                    "status": "error",
                    "error_type": paired.get("error_type", "H5muReadFailed"),
                    "details": paired.get("details", "Failed to read paired .h5mu."),
                    "warnings": warnings,
                }
            common_cells = list(paired["shared_barcodes"])
            if len(common_cells) < min_shared_cells:
                return {
                    "status": "error",
                    "error_type": "InsufficientSharedCells",
                    "details": (
                        f"Only {len(common_cells)} shared cells in paired .h5mu; "
                        f"peak-to-gene linking requires >= {min_shared_cells}."
                    ),
                }

            src_mdata = paired["mdata"]
            src_rna = paired["rna_adata"]
            src_atac = paired["atac_adata"]
            r_map = {str(bc): i for i, bc in enumerate(src_rna.obs_names)}
            a_map = {str(bc): i for i, bc in enumerate(src_atac.obs_names)}
            r_take = [r_map[bc] for bc in common_cells]
            a_take = [a_map[bc] for bc in common_cells]

            rna_is_raw = bool(paired.get("rna_is_raw_counts", False))
            atac_is_raw = bool(paired.get("atac_is_raw_counts", False))
            rna_base = (
                paired["rna_counts_matrix"]
                if paired.get("rna_counts_matrix") is not None
                else src_rna.X
            )
            atac_base = (
                paired["atac_counts_matrix"]
                if paired.get("atac_counts_matrix") is not None
                else src_atac.X
            )
            rna_mat = rna_base[r_take, :]
            atac_mat = atac_base[a_take, :]
            rna_var_names = [str(g) for g in src_rna.var_names]
            atac_var_names = _maybe_map_peaks_from_bed(
                [str(p) for p in src_atac.var_names], peaks_path
            )
            rna_counts_source = str(paired.get("rna_counts_source") or "X")
            atac_counts_source = str(paired.get("atac_counts_source") or "X")
            if not rna_is_raw:
                warnings.append(
                    "RNA modality did not pass integer raw-count verification; "
                    "using values without re-applying log1p(CP10k) normalization."
                )
            if not atac_is_raw:
                warnings.append(
                    "ATAC modality did not pass integer raw-count verification; "
                    "using values as provided."
                )
        else:
            rna_ad = _load_rna_adata(valid_rna[0])
            if rna_ad is None:
                return {
                    "status": "error",
                    "error_type": "UnsupportedRNAFormat",
                    "details": f"Cannot load RNA matrix from {valid_rna[0]}",
                }
            atac_ad, atac_err = _load_atac_adata(valid_atac)
            if atac_ad is None:
                return {
                    "status": "error",
                    "error_type": "ATACMatrixNotLoadable",
                    "details": atac_err or "Peak-to-gene requires a valid ATAC matrix.",
                    "warnings": warnings,
                }

            # Canonical deterministic barcode order preserving RNA obs_names order.
            atac_set = {str(bc) for bc in atac_ad.obs_names}
            seen_bc: set[str] = set()
            common_cells = []
            for bc in (str(x) for x in rna_ad.obs_names):
                if bc in atac_set and bc not in seen_bc:
                    common_cells.append(bc)
                    seen_bc.add(bc)

            if len(common_cells) < min_shared_cells:
                return {
                    "status": "error",
                    "error_type": "InsufficientSharedCells",
                    "details": (
                        f"Only {len(common_cells)} cells shared between RNA and ATAC. "
                        f"Peak-to-gene linking requires at least {min_shared_cells} shared cells."
                    ),
                }

            r_map = {str(bc): i for i, bc in enumerate(rna_ad.obs_names)}
            a_map = {str(bc): i for i, bc in enumerate(atac_ad.obs_names)}
            rna_sub = rna_ad[[r_map[bc] for bc in common_cells], :].copy()
            atac_sub = atac_ad[[a_map[bc] for bc in common_cells], :].copy()

            rna_chk = _select_raw_counts_matrix(
                rna_sub, modality_label="rna", require_raw_counts=False
            )
            atac_chk = _select_raw_counts_matrix(
                atac_sub, modality_label="atac", require_raw_counts=False
            )
            rna_is_raw = bool(rna_chk.get("is_raw_counts", False))
            atac_is_raw = bool(atac_chk.get("is_raw_counts", False))
            rna_mat = (
                rna_chk["matrix"] if rna_chk.get("matrix") is not None else rna_sub.X
            )
            atac_mat = (
                atac_chk["matrix"]
                if atac_chk.get("matrix") is not None
                else atac_sub.X
            )
            rna_var_names = [str(g) for g in rna_sub.var_names]
            atac_var_names = _maybe_map_peaks_from_bed(
                [str(p) for p in atac_sub.var_names], peaks_path
            )
            rna_counts_source = str(rna_chk.get("source") or "X")
            atac_counts_source = str(atac_chk.get("source") or "X")
            if not rna_is_raw:
                warnings.append(
                    "RNA input matrix is already normalized/non-integer; "
                    "skipping NormalizeData/log1p to prevent double normalization."
                )
            if not atac_is_raw:
                warnings.append(
                    "ATAC input matrix is non-integer; correlating values as provided."
                )

        if sp.issparse(rna_mat) and not sp.isspmatrix_csc(rna_mat):
            rna_mat = rna_mat.tocsc()
        if sp.issparse(atac_mat) and not sp.isspmatrix_csc(atac_mat):
            atac_mat = atac_mat.tocsc()

        gene_records, gtf_stats = parse_gtf_tss_detailed(
            gtf_resolved, wanted=set(rna_var_names)
        )
        if not gene_records:
            return {
                "status": "error",
                "error_type": "GTFGeneOverlapEmpty",
                "details": (
                    "Local GTF annotation contained zero genes overlapping RNA "
                    "feature names."
                ),
                "gtf_stats": gtf_stats,
            }

        core_res = compute_peak_to_gene_links(
            rna_matrix=rna_mat,
            rna_var_names=rna_var_names,
            atac_matrix=atac_mat,
            atac_var_names=atac_var_names,
            gene_records=gene_records,
            distance_bp=distance_bp,
            min_distance_bp=min_distance_bp,
            min_corr=min_corr,
            max_padj=max_padj,
            max_pvalue_signac=max_pvalue_signac,
            rna_is_raw_counts=rna_is_raw,
            atac_is_raw_counts=atac_is_raw,
            preprocess="paired",
            atac_transform=atac_transform,
            min_cells_detected=min_cells_detected,
            standard_chroms_only=standard_chroms_only,
            background_null=background_null,
            genome_fasta=fasta_resolved,
            n_bg_bins_per_axis=int(params.get("n_bg_bins_per_axis", 5) or 5),
            n_bg_per_bin=int(params.get("n_bg_per_bin", 200) or 200),
            min_bg_peaks=int(params.get("min_bg_peaks", 30) or 30),
            promoter_bp=promoter_bp,
            distance_mode="interval_overlap_zero",
            seed=seed,
        )

        if core_res.get("status") != "success":
            core_res.setdefault("warnings", warnings)
            return core_res

        warnings.extend(core_res.get("warnings") or [])
        links = list(core_res["links"])
        n_pos = int(core_res["n_positive"])
        n_neg = int(core_res["n_negative"])

        if not links:
            warnings.append(
                f"No peak-gene links passed thresholds (|r| >= {min_corr}"
                + (f", padj <= {max_padj}" if max_padj is not None else "")
                + ")."
            )
        elif n_neg > 0:
            warnings.append(
                f"{n_neg} peaks show NEGATIVE correlation with nearby genes "
                "(discordant accessibility-expression evidence). Flag for "
                "experimental validation before interpreting mechanism."
            )

        csv_path = _write_links_csv_atomic(output_dir / "peak_gene_links.csv", links)

        h5mu_output_path: str | None = None
        want_h5mu_write = bool(
            params.get("write_h5mu", False)
            or params.get("update_h5mu", False)
            or params.get("output_h5mu")
        )
        if want_h5mu_write and src_mdata is not None:
            target_h5mu = Path(
                str(
                    params.get("output_h5mu")
                    or (
                        resolved_h5mu
                        if params.get("update_h5mu", False)
                        else output_dir / "peak2gene_integrated.h5mu"
                    )
                )
            )
            target_h5mu.parent.mkdir(parents=True, exist_ok=True)
            uns_payload = sanitize_for_hdf5_uns(
                {
                    "validation_level": "beta",
                    "seed": int(seed),
                    "n_cells": int(core_res["n_cells"]),
                    "n_genes_tested": int(core_res["n_genes_tested"]),
                    "n_peaks_tested": int(core_res["n_peaks_tested"]),
                    "n_tested": int(core_res["n_tested"]),
                    "n_fdr_significant": int(core_res["n_fdr_significant_all_tested"]),
                    "n_links": len(links),
                    "n_positive_correlations": n_pos,
                    "n_negative_correlations": n_neg,
                    "n_promoter_proximal": int(core_res["n_promoter_proximal"]),
                    "n_distal_cis": int(core_res["n_distal_cis"]),
                    "distance_bp": int(distance_bp),
                    "distance_kb": distance_kb,
                    "min_corr": float(min_corr),
                    "max_padj": max_padj,
                    "atac_transform": atac_transform,
                    "background_method": core_res["background_method"],
                    "background_match": core_res["background_match"],
                    "links_csv_filename": "peak_gene_links.csv",
                    "top_links": links[:50],
                    "ctcf_validated": False,
                    "hic_corroborated": False,
                    "caveat": _ASSOCIATIVE_CAVEAT,
                }
            )
            src_mdata.uns["aria_peak2gene"] = uns_payload
            tmp_h5mu = target_h5mu.parent / (
                f".{target_h5mu.name}.tmp.{os.getpid()}.{uuid.uuid4().hex[:6]}.h5mu"
            )
            try:
                src_mdata.write_h5mu(str(tmp_h5mu))
                os.replace(str(tmp_h5mu), str(target_h5mu))
                h5mu_output_path = str(target_h5mu)
            finally:
                if tmp_h5mu.exists():
                    try:
                        tmp_h5mu.unlink()
                    except OSError:
                        pass

        return {
            "status": "success",
            "validation_level": "beta",
            "method": "paired_cis_peak_gene_background_z_bh_fdr",
            "organism": organism,
            "genome": genome,
            "seed": int(seed),
            "n_cells": int(core_res["n_cells"]),
            "n_genes_tested": int(core_res["n_genes_tested"]),
            "n_peaks_tested": int(core_res["n_peaks_tested"]),
            "n_tested": int(core_res["n_tested"]),
            "n_fdr_significant": int(core_res["n_fdr_significant_all_tested"]),
            "n_links": len(links),
            "n_positive_correlations": n_pos,
            "n_negative_correlations": n_neg,
            "n_promoter_proximal": int(core_res["n_promoter_proximal"]),
            "n_distal_cis": int(core_res["n_distal_cis"]),
            "distance_bp": int(distance_bp),
            "distance_kb": distance_kb,
            "min_corr": float(min_corr),
            "max_padj": max_padj,
            "max_pvalue_signac": max_pvalue_signac,
            "rna_counts_source": rna_counts_source,
            "atac_counts_source": atac_counts_source,
            "rna_normalization": core_res["rna_normalization"],
            "atac_normalization": core_res["atac_normalization"],
            "min_cells_detected": int(core_res["min_cells_detected_effective"]),
            "background_method": core_res["background_method"],
            "background_match": core_res["background_match"],
            "gc_matching": core_res["gc_matching"],
            "gtf_stats": gtf_stats,
            "top_links": links[:20],
            "output_path": csv_path,
            "h5mu_output_path": h5mu_output_path,
            "ctcf_validated": False,
            "hic_corroborated": False,
            "caveat": _ASSOCIATIVE_CAVEAT,
            "warnings": warnings,
        }

    except ImportError as exc:
        return {
            "status": "error",
            "error_type": "MissingDependency",
            "details": (
                f"Peak-to-gene dependencies are required: {exc}. "
                "ARIA never fabricates regulatory links."
            ),
            "warnings": warnings,
        }
    except Exception as exc:
        return {
            "status": "error",
            "error_type": "Peak2GeneFailed",
            "details": str(exc)[:500],
            "warnings": warnings,
        }


if __name__ == "__main__":
    run_script(integration_peak2gene)
