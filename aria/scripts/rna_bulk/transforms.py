"""Bulk RNA-seq expression transforms: VST, variable-gene selection, TPM, CPM
(A7 split of rna_bulk_de.py; ADR-075 adds CPM for 3'/5'-tag and UMI protocols).

Re-exported from aria.scripts.rna_bulk_de."""
from __future__ import annotations
import warnings

from aria.scripts.rna_bulk.gtf_io import _strip_gene_version


def _run_vst(counts_raw, metadata, warnings: list):
    """
    Variance-Stabilizing Transformation via pydeseq2.

    Returns a DataFrame (genes × samples) of VST-transformed values.
    Falls back to log2(normed+1) if pydeseq2 VST is unavailable.

    Used for: PCA, MDS, heatmaps at the sample level (NOT for DE testing).
    DESeq2 DE still receives raw counts — it has its own internal normalization.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        import pandas as pd
        import numpy as np

        # pydeseq2 expects samples × genes; use intercept-only design since
        # this is a normalization step, not a test for DE.
        dds = DeseqDataSet(
            counts=counts_raw.T.astype(int),
            metadata=metadata,
            design="~1",
            refit_cooks=False,
            quiet=True,
        )

        # Try VST first (preferred — fast, handles large datasets).
        # If not available in the installed version, fall back to rlog,
        # then to log2(normalized_counts+1).
        try:
            dds.fit_size_factors()
            dds.vst_fit(use_design=False)
            vst = dds.vst_transform()
        except (AttributeError, ValueError, RuntimeError):
            try:
                dds.deseq2()
                vst = np.log2(dds.layers["normed_counts"] + 1)
                warnings.append(
                    "pydeseq2 VST unavailable — using log2(size-factor normalized + 1)."
                )
            except Exception as e:
                warnings.append(
                    f"VST/rlog failed ({e}); falling back to raw log2. "
                    f"PCA/MDS may be affected by library size differences."
                )
                vst = np.log2(counts_raw.T.astype(float) + 1)

        # vst is samples × genes; transpose back to genes × samples
        vst_df = pd.DataFrame(
            np.asarray(vst).T,
            index=counts_raw.index,
            columns=counts_raw.columns,
        )
        return vst_df

    except (ImportError, AttributeError, ValueError, RuntimeError) as exc:
        warnings.append(
            f"pydeseq2 VST unavailable ({exc}) — using log2(counts+1) + lib-size scaling. "
            "PCA may be dominated by library size differences."
        )
        import pandas as pd
        import numpy as np
        lib = counts_raw.sum(axis=0)
        positive_lib = lib[lib > 0]
        target = float(positive_lib.median()) if len(positive_lib) else 1.0
        scale_factor = target / lib.replace(0, np.nan)
        normed = (counts_raw * scale_factor).fillna(0.0)
        return pd.DataFrame(
            np.log2(normed.astype(float) + 1),
            index=counts_raw.index,
            columns=counts_raw.columns,
        )


def _select_variable_genes(matrix, n_top: int = 2000,
                             biotype_map: dict | None = None,
                             warnings: list | None = None):
    """
    Select the top-N most variable genes from a VST-transformed matrix.

    If biotype_map is provided, restrict to protein_coding genes first.
    Returns (filtered_matrix, n_protein_coding, n_after_variance_filter).

    Args:
        matrix:      DataFrame (genes × samples) — should be VST-transformed.
        n_top:       number of most-variable genes to keep.
        biotype_map: optional {ensembl_id_no_version: biotype_string}.
        warnings:    list to append advisory messages to.
    """
    if warnings is None:
        warnings = []

    # Strip Ensembl version suffix from matrix index (for biotype lookup)
    if biotype_map:
        def _lookup_biotype(gid):
            return biotype_map.get(_strip_gene_version(str(gid)), "unknown")
        biotypes = matrix.index.map(_lookup_biotype)
        n_pc     = int((biotypes == "protein_coding").sum())
        n_total  = len(matrix)
        pc_frac  = n_pc / max(n_total, 1)

        if n_pc < 500:
            warnings.append(
                f"Only {n_pc} protein_coding genes found in matrix "
                f"({pc_frac:.0%} of {n_total}). Falling back to all biotypes "
                f"for DR — results may be affected by pseudogenes/ncRNAs."
            )
            n_pc = 0
        elif pc_frac < 0.70 and biotype_map:
            warnings.append(
                f"GTF annotation has unusually low protein_coding fraction "
                f"({pc_frac:.0%}). Expected ~70-85% for human/mouse."
            )
            matrix_pc = matrix[biotypes == "protein_coding"]
            matrix    = matrix_pc
        else:
            matrix = matrix[biotypes == "protein_coding"]
    else:
        n_pc = 0

    # Top-N most variable
    variance = matrix.var(axis=1)
    n_keep   = min(n_top, len(matrix))
    top_idx  = variance.nlargest(n_keep).index

    return matrix.loc[top_idx], n_pc, n_keep


def _compute_tpm(counts_raw, gene_lengths: dict, warnings: list):
    """
    Compute TPM (Transcripts Per Million) from raw counts.

    TPM = (reads_per_gene / gene_length_kb) / (sum_of_all / 1e6)

    TPM is NOT used by ARIA for DE or PCA (both use better methods).
    It's computed as a supplementary table for downstream tools that
    require TPM input (e.g., ssGSEA, single-sample deconvolution).

    Args:
        counts_raw:   DataFrame (genes × samples), raw integer counts.
        gene_lengths: {ensembl_id: length_bp} — from GTF exon sum.
                     If missing, returns None with a warning.
    """
    try:
        import pandas as pd
        import numpy as np

        if not gene_lengths:
            warnings.append(
                "Gene lengths not available from GTF — cannot compute TPM. "
                "DE analysis is unaffected (uses raw counts)."
            )
            return None

        # Map matrix rows (possibly versioned IDs) to lengths
        def _get_length(gid):
            clean = _strip_gene_version(str(gid))
            return gene_lengths.get(clean)

        lengths = pd.to_numeric(counts_raw.index.map(_get_length), errors="coerce")
        has_len = np.isfinite(lengths) & (lengths > 0)
        n_with_len = int(has_len.sum())

        if n_with_len < len(counts_raw) * 0.5:
            warnings.append(
                f"TPM: only {n_with_len}/{len(counts_raw)} genes have lengths "
                f"in GTF — TPM table will be incomplete."
            )

        if n_with_len == 0:
            warnings.append("TPM: no genes have valid positive lengths; no table computed.")
            return None

        # Drop genes without lengths (can't TPM-normalize them)
        mat      = counts_raw.loc[has_len]
        lens_kb  = pd.Series(
            lengths[has_len].astype(float) / 1000.0,
            index=mat.index,
        )

        # RPK = reads per kilobase per gene
        rpk = mat.div(lens_kb, axis=0)

        # Scaling factor = sum of RPK per sample / 1e6
        scale = rpk.sum(axis=0) / 1e6

        tpm = rpk.div(scale.replace(0, np.nan), axis=1).fillna(0.0)
        tpm = tpm.round(3)
        return tpm

    except Exception as e:
        warnings.append(f"TPM computation failed: {e}")
        return None


def _compute_cpm(counts_raw, warnings: list):
    """
    Compute CPM (Counts Per Million) from raw integer or UMI counts.

    CPM = reads_or_umis_per_gene / (sample_library_total / 1e6)

    CPM is library-size normalized without transcript-length division, making
    it the valid supplementary expression table for 3'-tag, 5'-tag, inline-UMI,
    and pooled-plate UMI bulk RNA-seq protocols where dividing by exon-union
    length (TPM) introduces severe length bias (ADR-075). NOT used for DE or PCA.
    """
    try:
        import numpy as np

        lib_millions = counts_raw.sum(axis=0).astype(float) / 1e6
        lib_safe = lib_millions.replace(0.0, np.nan)
        cpm = counts_raw.astype(float).div(lib_safe, axis=1).fillna(0.0)
        return cpm.round(3)
    except Exception as e:
        warnings.append(f"CPM computation failed: {e}")
        return None


# ══════════════════════════════════════════════════════════════════════════════


