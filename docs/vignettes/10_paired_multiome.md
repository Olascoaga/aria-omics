# Vignette `10_paired_multiome` — Human Hippocampus Same-Cell 10x Multiome RNA+ATAC Cohort (3v3 Donors)

- **ARIA Modality Lane:** `Paired_Multiome`
- **Public Accession:** `GEO GSE278576`
- **Publication / DOI:** Zemke NR et al., Science (2026) doi:10.1126/science.adt8307 (`10.1126/science.adt8307`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `949.4 s`

## 1. Dataset Characteristics & Public Provenance

6 post-mortem human hippocampal donors (3 young_20_39 [hc29, hc77, hc76] vs 3 aged_80_100 [hc11, hc12, hc19], 400 paired same-cell nuclei/donor = 2,400 paired nuclei across chr1/6/7/11/17/19).

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `hc11_paired_multiome.h5mu` | `25,252,046` | `5af99305a2013205471501b9...` |
| `hc12_paired_multiome.h5mu` | `32,123,896` | `e3c1fd0a26087b07548969c4...` |
| `hc19_paired_multiome.h5mu` | `31,063,020` | `d2f34d86000f69f8c47b110a...` |
| `hc29_paired_multiome.h5mu` | `35,558,770` | `0a99dc100734ba0263ff4a61...` |
| `hc76_paired_multiome.h5mu` | `35,757,826` | `b5b57aba7a8b0a95fe823bfe...` |
| `hc77_paired_multiome.h5mu` | `32,506,000` | `5609bbade21bdbc4e5abdd1d...` |

**Public Source URLs:**
- `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278576`

## 2. Exact Biological Prompt & Experimental Contrast

```text
Perform joint same-cell 10x Multiome (scRNA-seq + scATAC-seq) integration across 6 human hippocampal donors comparing aged_80_100 (hc11, hc12, hc19) versus young_20_39 (hc29, hc77, hc76) individuals (GEO GSE278576, Homo sapiens hg38). Execute multi-donor cohort assembly, paired joint RNA+ATAC QC, Weighted Nearest Neighbor (WNN) multimodal graph construction, cis peak-to-gene regulatory linkage, candidate TF->peak->gene eRegulon assembly, and MOFA+ multi-view latent factor decomposition.
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
| `CP 3` | WNN Integration — k Parameter | `Use recommended k=20` |
| `CP 3` | MOFA+ factor count | `Use recommended n_factors=10` |
| `CP 5` | Analysis complete. | `Accept report` |

## 4. Deterministic Pipeline Execution (`RunLedger`)

Total recorded `RunLedger` nodes: **21**.

| Node ID | Modality | Analysis Stage | Label | Status |
|---|---|---|---|---|
| `ledger://chromatin/qc` | `chromatin` | `qc` | Quality control | **`error`** |
| `ledger://chromatin/scatac_preprocessing` | `chromatin` | `scatac_preprocessing` | scATAC raw preprocessing & matrix construction | **`ran`** |
| `ledger://chromatin/bulk_atac_preprocessing` | `chromatin` | `bulk_atac_preprocessing` | Raw bulk chromatin alignment and filtering | **`not_run`** |
| `ledger://chromatin/dimensionality_reduction` | `chromatin` | `dimensionality_reduction` | Dimensionality reduction (LSI/TF-IDF) | **`ran`** |
| `ledger://chromatin/peak_calling` | `chromatin` | `peak_calling` | Peak calling | **`not_run`** |
| `ledger://chromatin/peak_count_matrix` | `chromatin` | `peak_count_matrix` | Peak-count matrix | **`not_run`** |
| `ledger://chromatin/differential_accessibility` | `chromatin` | `differential_accessibility` | Differential accessibility / binding | **`error`** |
| `ledger://chromatin/motif_enrichment` | `chromatin` | `motif_enrichment` | TF motif enrichment | **`ran`** |
| `ledger://chromatin/regulatory_layers` | `chromatin` | `regulatory_layers` | Regulatory layers | **`ran`** |
| `ledger://chromatin/matrix_pipeline` | `chromatin` | `matrix_pipeline` | scATAC matrix pipeline | **`not_run`** |
| `ledger://chromatin/peak_annotation` | `chromatin` | `peak_annotation` | Genomic peak annotation | **`not_run`** |
| `ledger://chromatin/peak_ora` | `chromatin` | `peak_ora` | DA-peak functional enrichment | **`not_run`** |
| `ledger://chromatin/differential_tf_footprinting` | `chromatin` | `differential_tf_footprinting` | Differential TF footprinting | **`ran`** |
| `ledger://integration/preprocessing` | `integration` | `preprocessing` | Multiome raw assembly & cohort union | **`executed`** |
| `ledger://integration/paired_qc` | `integration` | `paired_qc` | Paired multiome QC | **`ran`** |
| `ledger://integration/wnn` | `integration` | `wnn` | Weighted Nearest Neighbors (WNN) | **`ran`** |
| `ledger://integration/mofa` | `integration` | `mofa` | MOFA+ factor decomposition | **`ran`** |
| `ledger://integration/peak_to_gene` | `integration` | `peak_to_gene` | Cis peak-to-gene linkage | **`ran`** |
| `ledger://integration/eregulons` | `integration` | `eregulons` | Candidate eRegulon inference | **`ran`** |
| `ledger://integration/cross_modal_triangulation` | `integration` | `cross_modal_triangulation` | Cross-modal evidence triangulation | **`ran`** |
| `ledger://report/executive_summary` | `report` | `executive_summary` | Executive summary | **`ran`** |

## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine

- **Total Audited Public Claims:** `11` (`{"descriptive": 4, "associative": 7}`)
- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `pass`
- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `0`

### Representative Compiled Claims

- **`[descriptive]` (`executive_summary`):** Paired multiome QC retained 1828 of 2400 shared barcodes across 36601 RNA genes and 55385 ATAC peaks. WNN manifold integration identified 7 joint clusters (WNN vs RNA ARI 0.6702, WNN vs ATAC ARI 0.5765, guardrail status pass). MOFA+ retained 6 of 10 fitted factors (1 confounder factor(s) quarantined); cis peak-to-gene linkage identified 21 BH-FDR associations; candidate eRegulon inference retained 103 associative modules across 3971 triplets. All cross-modal links and candidate eRegulons represent observational co-variation rather than proof of transcriptional causation, and single-modality analyses remain the primary reference.
- **`[associative]` (`scatac_preprocessing`):** scATAC preprocessing completed (cohort matrix union): the resulting matrix has 2400 cell(s) x 55385 peak(s).
- **`[associative]` (`dimensionality_reduction`):** TF-IDF/LSI over 55385 peaks in 2400 cells produced 11 clusters after removing 1 depth-associated component(s).
- **`[associative]` (`regulatory_layers`):** scATAC regulatory-layer analysis ran 1 of 5 optional layer(s); skipped layers are reported with their prerequisites.
- **`[descriptive]` (`preprocessing`):** Multiome intake (paired_cohort_merge, platform paired_multiome) produced a paired RNA+ATAC container with 6 sample(s), 2400 cells, 36601 RNA genes, 55385 consensus ATAC peaks under barcode translation shared_mudata_obs.
- **`[descriptive]` (`paired_qc`):** Paired multiome QC retained 1828 of 2400 shared barcodes (23.83% removed) across 36601 RNA genes and 55385 ATAC peaks.

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `19` total (`19` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `56,585` prompt + `9,192` completion = `65,777` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`1,664,189` bytes, Capsule ZIP=`3,067,000` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 10_paired_multiome --force
```
