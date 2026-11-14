# Vignette `04_scatac_multidonor` — Human Hippocampus Lifespan Single-Cell Chromatin Accessibility

- **ARIA Modality Lane:** `scATAC`
- **Public Accession:** `GEO GSE278576`
- **Publication / DOI:** Zemke NR et al., Science (2026) doi:10.1126/science.adt8307 (`10.1126/science.adt8307`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `419.7 s`

## 1. Dataset Characteristics & Public Provenance

6 post-mortem human hippocampal donors (3 young_20_39 [hc29, hc77, hc76] vs 3 aged_80_100 [hc11, hc12, hc19], 500 nuclei/donor = 3,000 single-nucleus ATAC profiles across chr6/7/11/17/19).

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `hc11_scatac.h5mu` | `6,981,721` | `7692a4f42978e6323a09d213...` |
| `hc12_scatac.h5mu` | `11,547,858` | `9319b9b08906f3588bb1f9a9...` |
| `hc19_scatac.h5mu` | `10,239,302` | `0d48c99fc3fdb9552ea21dd1...` |
| `hc29_scatac.h5mu` | `12,513,625` | `54491c5ad56478a53884d82f...` |
| `hc76_scatac.h5mu` | `14,069,361` | `03cef6ee4655dc6caee64108...` |
| `hc77_scatac.h5mu` | `12,258,980` | `b6103fe5b6978c42e513b2ad...` |

**Public Source URLs:**
- `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278576`

## 2. Exact Biological Prompt & Experimental Contrast

```text
Investigate single-nucleus chromatin accessibility across 6 post-mortem human hippocampal donors comparing aged_80_100 (hc11, hc12, hc19) versus young_20_39 (hc29, hc77, hc76) individuals (GEO GSE278576, Homo sapiens hg38). Perform scATAC QC, TF-IDF and LSI dimensionality reduction (discarding depth-correlated component 1), donor-aware chromatin clustering, donor-level pseudobulk differential accessibility comparing aged_80_100 versus young_20_39, JASPAR TF motif enrichment, and cis gene annotation.
```

- **Primary Contrast(s):** `[['aged_80_100', 'young_20_39']]` (`condition_col="condition"`, `replicate_col="donor_id"`) 
- **Execution Flags:** `reproducible_mode=True`, `enable_hypotheses=True`

## 3. Interactive Checkpoint Audit Trail

| Checkpoint | Prompt Summary | Chosen Option |
|---|---|---|
| `CP 1` | 📁 ARIA found the following data: | `Confirm and continue` |
| `CP 2.1` | ARIA proposes the following biological groups: | `Yes — confirm groups` |
| `CP 2.2` | Please confirm the organism and genome assembly: | `Homo sapiens (hg38)` |
| `CP 2.3` | Main experimental factor: proposed 'condition' | `condition` |
| `CP 2.6` | Experimental design ready: | `Yes — proceed` |
| `CP 2` | ARIA Analysis Plan: | `Run recommended + optional supported analyses` |
| `CP 5` | Analysis complete. | `Accept report` |

## 4. Deterministic Pipeline Execution (`RunLedger`)

Total recorded `RunLedger` nodes: **14**.

| Node ID | Modality | Analysis Stage | Label | Status |
|---|---|---|---|---|
| `ledger://chromatin/qc` | `chromatin` | `qc` | Quality control | **`ran`** |
| `ledger://chromatin/scatac_preprocessing` | `chromatin` | `scatac_preprocessing` | scATAC raw preprocessing & matrix construction | **`ran`** |
| `ledger://chromatin/bulk_atac_preprocessing` | `chromatin` | `bulk_atac_preprocessing` | Raw bulk chromatin alignment and filtering | **`not_run`** |
| `ledger://chromatin/dimensionality_reduction` | `chromatin` | `dimensionality_reduction` | Dimensionality reduction (LSI/TF-IDF) | **`ran`** |
| `ledger://chromatin/peak_calling` | `chromatin` | `peak_calling` | Peak calling | **`not_run`** |
| `ledger://chromatin/peak_count_matrix` | `chromatin` | `peak_count_matrix` | Peak-count matrix | **`not_run`** |
| `ledger://chromatin/differential_accessibility` | `chromatin` | `differential_accessibility` | Differential accessibility / binding | **`ran`** |
| `ledger://chromatin/motif_enrichment` | `chromatin` | `motif_enrichment` | TF motif enrichment | **`ran`** |
| `ledger://chromatin/regulatory_layers` | `chromatin` | `regulatory_layers` | Regulatory layers | **`ran`** |
| `ledger://chromatin/matrix_pipeline` | `chromatin` | `matrix_pipeline` | scATAC matrix pipeline | **`not_run`** |
| `ledger://chromatin/peak_annotation` | `chromatin` | `peak_annotation` | Genomic peak annotation | **`skipped`** |
| `ledger://chromatin/peak_ora` | `chromatin` | `peak_ora` | DA-peak functional enrichment | **`skipped`** |
| `ledger://chromatin/differential_tf_footprinting` | `chromatin` | `differential_tf_footprinting` | Differential TF footprinting | **`ran`** |
| `ledger://report/executive_summary` | `report` | `executive_summary` | Executive summary | **`ran`** |

## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine

- **Total Audited Public Claims:** `8` (`{"descriptive": 2, "associative": 6}`)
- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `pass`
- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `0`

### Representative Compiled Claims

- **`[descriptive]` (`executive_summary`):** To address age-related chromatin accessibility changes in the human hippocampus, scATAC-seq analysis of 3,000 cells identified 10,022 differentially accessible peaks across 11 clusters with medium confidence. Motif enrichment analysis of these peaks, tested against 879 JASPAR2024 motifs, revealed significant enrichment in 9 of the 11 cluster groups. However, the absence of FRiP and TSS enrichment metrics limits the assessment of data quality and regulatory specificity. Validation of specific peak-to-gene links and functional characterization of the enriched transcription factors are required to establish causal regulatory mechanisms.
- **`[associative]` (`scatac_preprocessing`):** scATAC preprocessing completed (cohort matrix union): the resulting matrix has 3000 cell(s) x 40057 peak(s).
- **`[descriptive]` (`qc`):** scATAC QC measured 3000 cells and 40057 peaks.
- **`[associative]` (`dimensionality_reduction`):** TF-IDF/LSI over 40057 peaks in 3000 cells produced 11 clusters after removing 1 depth-associated component(s).
- **`[associative]` (`differential_accessibility`):** Per-cluster accessibility testing identified 10022 differentially accessible peaks across 11 cluster(s) at padj<0.05, |log2FC|>1.0.
- **`[associative]` (`differential_accessibility`):** Replicate-level pseudobulk DA ran 1 comparison(s) and found 0 differentially accessible peaks via the shared DESeq2 core.

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `3` total (`3` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `2,981` prompt + `587` completion = `3,568` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`146,703` bytes, Capsule ZIP=`1,765,062` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 04_scatac_multidonor --force
```
