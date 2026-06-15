"""
ARIA Peak-to-Gene Linking Script
----------------------------------
Correlates chromatin accessibility at ATAC peaks with nearby gene
expression to identify regulatory relationships.
Executed inside aria-integration-env by EnvironmentManager.

Method:
  For each peak, find all genes within distance_kb window.
  Compute Pearson correlation between peak accessibility
  and gene expression across all cells.
  Report links where |r| >= min_corr.

Interpretation boundary:
  Positive and negative correlations are descriptive evidence only. ARIA
  reports the sign and magnitude but does not assign regulatory mechanism from
  correlation alone.

High-confidence criteria (stored as Tunnels in ARIAMemory):
  - |r| >= 0.4 (strong correlation)
  - Distance < 250kb
  - Ideally: corroborated by an independent orthogonal assay

Input params:
    rna_files:    list  — scRNA files
    atac_files:   list  — scATAC files or peaks BED
    peaks_path:   str   — path to called peaks (.narrowPeak or .bed)
    genome:       str
    organism:     str
    distance_kb:  int   — search window around TSS (default: 500)
    min_corr:     float — minimum |correlation| threshold (default: 0.3)
    output_dir:   str

Output:
    {
      "status":                 "success",
      "n_links":                int,
      "n_positive_correlations": int,  — open + expressed
      "n_negative_correlations": int,  — open + silent (biologically important)
      "top_links":              [{"gene", "peak", "correlation", "distance_kb"}],
      "output_path":            str,   — CSV of all links
      "ctcf_validated":         bool,
      "hic_corroborated":       bool,
      "warnings":               [str]
    }
"""

from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from aria.scripts._base import run_script


def integration_peak2gene(params: dict) -> dict:
    from pathlib import Path
    import numpy as np

    rna_files   = params.get("rna_files", [])
    atac_files  = params.get("atac_files", [])
    peaks_path  = params.get("peaks_path", "")
    genome      = params.get("genome", "hg38")
    organism    = params.get("organism", "Homo sapiens")
    distance_kb = int(params.get("distance_kb", 500))
    min_corr    = float(params.get("min_corr", 0.3))
    output_dir  = params.get("output_dir", "/tmp/aria_p2g")
    warnings    = []

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Validate files
    valid_rna  = [f for f in rna_files  if Path(f).exists()]
    valid_atac = [f for f in atac_files if Path(f).exists()]
    peaks_ok   = Path(peaks_path).exists() if peaks_path else False

    if not valid_rna:
        return {
            "status":     "error",
            "error_type": "MissingRNA",
            "details":    "No valid RNA files for peak-to-gene linking.",
        }

    try:
        import scanpy as sc
        import anndata as ad
        import pandas as pd
        import numpy as np
        from scipy.stats import pearsonr

        # ── Load RNA ──────────────────────────────────────────────────────
        rna_path = valid_rna[0]
        p = Path(rna_path)
        if p.suffix == ".h5ad":
            rna = sc.read_h5ad(str(p))
        elif p.suffix == ".h5":
            rna = sc.read_10x_h5(str(p))
        elif p.is_dir():
            rna = sc.read_10x_mtx(str(p), var_names="gene_symbols",
                                    cache=True)
        else:
            return {
                "status":     "error",
                "error_type": "UnsupportedRNAFormat",
                "details":    f"Cannot load RNA from {rna_path}",
            }

        # Normalize RNA
        sc.pp.normalize_total(rna, target_sum=1e4)
        sc.pp.log1p(rna)

        # ── Load or create peak accessibility matrix ───────────────────────
        if valid_atac:
            atac = _load_atac_matrix(valid_atac, peaks_path, rna.obs_names)
        else:
            warnings.append(
                "No ATAC files found. Using promoter accessibility proxy."
            )
            atac = None

        if atac is None:
            return {
                "status":     "error",
                "error_type": "ATACMatrixNotLoadable",
                "details":    "Peak-to-gene linking requires a valid ATAC matrix.",
                "warnings":   warnings,
            }

        # Align cells
        common = list(set(rna.obs_names) & set(atac.obs_names))
        if len(common) < 100:
            return {
                "status":     "error",
                "error_type": "InsufficientSharedCells",
                "details":    (
                    f"Only {len(common)} cells shared between RNA and ATAC. "
                    f"Peak-to-gene linking requires at least 100 shared cells."
                ),
            }

        rna_aligned  = rna[common]
        atac_aligned = atac[common]

        # ── Load gene coordinates ──────────────────────────────────────────
        gene_coords = _get_gene_coordinates(rna.var_names, genome, organism)

        # ── Load peak coordinates ──────────────────────────────────────────
        peak_coords = _get_peak_coordinates(peaks_path, atac.var_names)

        # ── Compute peak-gene correlations ────────────────────────────────
        links = []
        n_tested = 0

        # Get dense matrices (only feasible for selected genes/peaks)
        # For full datasets: use sparse matrix operations
        rna_matrix  = rna_aligned.X
        atac_matrix = atac_aligned.X

        # DeepSeek P0 fix: never densify large sparse matrices
        import scipy.sparse as sp
        if sp.issparse(rna_matrix) and not sp.isspmatrix_csc(rna_matrix):
            rna_matrix  = rna_matrix.tocsc()
        if sp.issparse(atac_matrix) and not sp.isspmatrix_csc(atac_matrix):
            atac_matrix = atac_matrix.tocsc()

        def _get_col_sparse(mat, idx: int) -> "np.ndarray":
            """Extract one column without densifying the whole matrix."""
            if sp.issparse(mat):
                return np.asarray(mat.getcol(idx).todense()).ravel()
            return mat[:, idx]

        for gene_idx, gene in enumerate(rna_aligned.var_names[:5000]):  # limit for speed
            if gene not in gene_coords:
                continue

            gene_chrom, gene_start, gene_end = gene_coords[gene]
            gene_expr = rna_matrix[:, gene_idx]

            # Skip if gene not expressed
            if gene_expr.mean() < 0.1:
                continue

            # Find peaks within distance_kb
            for peak_idx, peak in enumerate(atac.var_names):
                if peak not in peak_coords:
                    continue

                peak_chrom, peak_start, peak_end = peak_coords[peak]

                if peak_chrom != gene_chrom:
                    continue

                # Distance from peak center to gene TSS
                peak_center = (peak_start + peak_end) // 2
                dist_bp     = abs(peak_center - gene_start)
                dist_kb     = dist_bp / 1000

                if dist_kb > distance_kb:
                    continue

                n_tested += 1
                peak_acc  = atac_matrix[:, peak_idx]

                if peak_acc.std() < 1e-6:
                    continue

                try:
                    r, p_val = pearsonr(gene_expr, peak_acc)
                except Exception:
                    continue

                if abs(r) >= min_corr:
                    links.append({
                        "gene":        gene,
                        "peak":        peak,
                        "correlation": round(float(r), 4),
                        "p_value":     round(float(p_val), 6),
                        "distance_kb": round(float(dist_kb), 1),
                        "direction":   "positive" if r > 0 else "negative",
                    })

        if not links:
            warnings.append(
                f"No peak-gene links found with |r| >= {min_corr}. "
                f"Consider lowering min_corr or checking data quality."
            )
            return {
                "status":                  "success",
                "n_links":                 0,
                "n_positive_correlations": 0,
                "n_negative_correlations": 0,
                "top_links":               [],
                "ctcf_validated":          False,
                "hic_corroborated":        False,
                "warnings":                warnings,
            }

        # Sort by |correlation|
        links.sort(key=lambda x: abs(x["correlation"]), reverse=True)

        n_pos = sum(1 for l in links if l["direction"] == "positive")
        n_neg = sum(1 for l in links if l["direction"] == "negative")

        if n_neg > 0:
            warnings.append(
                f"{n_neg} peaks show NEGATIVE correlation with nearby genes "
                f"(discordant accessibility-expression evidence). Flag for "
                f"experimental validation before interpreting mechanism."
            )

        # Save to CSV
        output_path = str(Path(output_dir) / "peak_gene_links.csv")
        pd.DataFrame(links).to_csv(output_path, index=False)

        return {
            "status":                  "success",
            "n_links":                 len(links),
            "n_tested":                int(n_tested),
            "n_positive_correlations": int(n_pos),
            "n_negative_correlations": int(n_neg),
            "top_links":               links[:20],
            "output_path":             output_path,
            "ctcf_validated":          False,   # requires ChIP data
            "hic_corroborated":        False,   # requires HiC data
            "warnings":                warnings,
        }

    except ImportError as e:
        return {
            "status":     "error",
            "error_type": "MissingDependency",
            "details":    (
                f"Peak-to-gene dependencies are required: {e}. "
                "Mock peak-to-gene links were removed; ARIA never fabricates "
                "regulatory links."
            ),
            "warnings":   warnings,
        }
    except Exception as e:
        return {
            "status":     "error",
            "error_type": "Peak2GeneFailed",
            "details":    str(e)[:500],
        }


# ── Helper functions ──────────────────────────────────────────────────────────

def _load_atac_matrix(atac_files: list, peaks_path: str,
                       cell_barcodes) -> object:
    """Load ATAC peak-by-cell matrix."""
    try:
        import anndata as ad
        import numpy as np

        # Try loading pre-built matrix (e.g. from cellranger-arc)
        for f in atac_files:
            if f.endswith(".h5ad"):
                return ad.read_h5ad(f)

        # If only fragment files: count fragments in peaks
        # Requires episcanpy or snapatac2
        try:
            import episcanpy.api as epi
            # Full implementation: epi.count_fragments_in_peaks()
            # For now: return None to trigger mock
            return None
        except ImportError:
            return None

    except Exception:
        return None


def _get_gene_coordinates(gene_names, genome: str, organism: str) -> dict:
    """
    Get TSS coordinates for genes from GTF annotation file.

    DeepSeek P0: This function MUST NOT return random coordinates.
    If a GTF is unavailable, it raises MissingGTFError so the caller
    can return a structured error rather than produce fake biology.

    GTF lookup order:
      1. ~/.aria/genomes/{genome}/{genome}.gtf.gz  (auto-downloaded)
      2. ARIA_GTF_PATH environment variable
      3. Explicit gtf_path parameter (future)

    Returns dict: gene -> (chrom, tss_start, tss_end)
    Raises MissingGTFError if no GTF is found.
    """
    import os
    from pathlib import Path

    # Search for GTF in standard locations
    candidates = [
        Path.home() / ".aria" / "genomes" / genome / f"{genome}.gtf.gz",
        Path.home() / ".aria" / "genomes" / genome / f"{genome}.gtf",
        Path(os.environ.get("ARIA_GTF_PATH", "/nonexistent")),
    ]

    gtf_path = next((p for p in candidates if p.exists()), None)

    if gtf_path is None:
        raise MissingGTFError(
            f"GTF annotation not found for genome '{genome}'. "
            f"ARIA cannot compute real gene coordinates without it. "
            f"To fix: run 'aria download-genome {genome}' or set "
            f"ARIA_GTF_PATH=/path/to/{genome}.gtf"
        )

    # Parse GTF for TSS coordinates
    coords = {}
    opener = __import__("gzip").open if str(gtf_path).endswith(".gz") else open
    try:
        with opener(str(gtf_path), "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue
                chrom  = parts[0]
                start  = int(parts[3])
                end    = int(parts[4])
                strand = parts[6]
                attrs  = parts[8]

                # Extract gene name
                gene_name = None
                for attr in attrs.split(";"):
                    attr = attr.strip()
                    if attr.startswith("gene_name"):
                        gene_name = attr.split('"')[1] if '"' in attr else                                     attr.split(" ")[-1]
                        break

                if gene_name and gene_name in gene_names:
                    # TSS = start for +, end for -
                    tss = start if strand == "+" else end
                    coords[gene_name] = (chrom, tss, tss + 1)

                if len(coords) >= len(gene_names):
                    break
    except Exception as e:
        raise MissingGTFError(f"Failed to parse GTF at {gtf_path}: {e}")

    return coords


class MissingGTFError(Exception):
    """Raised when GTF annotation is required but not found."""
    pass


def _get_peak_coordinates(peaks_path: str, peak_names) -> dict:
    """
    Parse peak coordinates from peak names or BED file.
    Peak names typically encode position: chr1:1000-2000
    """
    coords = {}

    # Try parsing from peak names (e.g. "chr1:100000-101000")
    for peak in peak_names:
        try:
            if ":" in peak and "-" in peak:
                chrom, pos = peak.split(":")
                start, end = pos.split("-")
                coords[peak] = (chrom, int(start), int(end))
        except Exception:
            pass

    # If names don't encode positions, try reading BED file
    if not coords and peaks_path and __import__("pathlib").Path(peaks_path).exists():
        try:
            with open(peaks_path) as f:
                for i, line in enumerate(f):
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        peak_id = f"{parts[0]}:{parts[1]}-{parts[2]}"
                        coords[peak_id] = (parts[0], int(parts[1]), int(parts[2]))
                    if i > 500_000:
                        break
        except Exception:
            pass

    return coords


if __name__ == "__main__":
    run_script(integration_peak2gene)
