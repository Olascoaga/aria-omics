# Vignette `11_cross_modal_concordance` — ENCODE Matched Human Hematopoietic polyA+ RNA-seq x Bulk ATAC-seq (GM12878 vs K562)

- **ARIA Modality Lane:** `Cross_Modal_RNA_x_ATAC`
- **Public Accession:** `ENCODE ENCSR843RJV, ENCSR889TRN, ENCSR594NJP, ENCSR596ACL, ENCSR637XSC, ENCSR483RKN, ENCSR859USB`
- **Publication / DOI:** ENCODE Project Consortium, Nature 583, 699-710 (2020) (`10.1038/s41586-020-2493-4`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `3082.4 s`

## 1. Dataset Characteristics & Public Provenance

N=6 matched biological units (GM12878_rep1..3 vs K562_rep1..3) profiled across both bulk polyA+ RNA-seq and bulk ATAC-seq on GRCh38 chr6, chr7, chr11, chr17, chr19.

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `GM12878_rep1.bam` | `2,110,984,885` | `786eb09af7c40c84be508d23...` |
| `GM12878_rep2.bam` | `1,759,054,627` | `7f7eb74403ff195f8cd96d02...` |
| `GM12878_rep3.bam` | `2,134,120,685` | `32d1685aa2ba66808343b850...` |
| `K562_rep1.bam` | `840,054,452` | `117a1af6cdb1e92177b0ab0b...` |
| `K562_rep2.bam` | `792,289,628` | `a2363caa9b1a92d85a724e7d...` |
| `K562_rep3.bam` | `173,245,367` | `8fd78d5bb93313d2afae8d73...` |
| `encode_gm12878_vs_k562_counts.tsv` | `187,869` | `ed54ed9f85f93a23eb73890f...` |

**Public Source URLs:**
- `https://www.encodeproject.org/files/ENCFF700JWD/@@download/ENCFF700JWD.tsv`
- `https://www.encodeproject.org/files/ENCFF935TRK/@@download/ENCFF935TRK.tsv`
- `https://www.encodeproject.org/files/ENCFF480HNH/@@download/ENCFF480HNH.tsv`
- `https://www.encodeproject.org/files/ENCFF643FJJ/@@download/ENCFF643FJJ.tsv`
- `https://www.encodeproject.org/files/ENCFF939LDL/@@download/ENCFF939LDL.tsv`
- `https://www.encodeproject.org/files/ENCFF241QPH/@@download/ENCFF241QPH.tsv`

## 2. Exact Biological Prompt & Experimental Contrast

```text
Evaluate associative cross-modal regulatory concordance across matched ENCODE K562 erythroleukemia and GM12878 B-lymphoblastoid biological replicates (N=6 shared sample IDs, 3 per lineage) profiled independently by bulk polyA+ RNA-seq and bulk ATAC-seq on GRCh38. Run pyDESeq2 differential expression and differential accessibility (K562 versus GM12878), then quantify cis DE x DA directional concordance, TF expression-motif triangulation, shared pathway convergence, and sample-matched cis Spearman correlation without claiming single-cell pairing.
```

- **Primary Contrast(s):** `[['K562', 'GM12878']]` (`condition_col="condition"`, `replicate_col="replicate"`) 
- **Execution Flags:** `reproducible_mode=True`, `enable_hypotheses=True`

## 3. Interactive Checkpoint Audit Trail

| Checkpoint | Prompt Summary | Chosen Option |
|---|---|---|
| `CP 1` | 📁 ARIA found the following data: | `Confirm and continue` |
| `CP 2.1` | ARIA proposes the following biological groups: | `Yes — confirm groups` |
| `CP 2.2` | Please confirm the organism and genome assembly: | `Homo sapiens (hg38)` |
| `CP 2.3` | Main experimental factor: proposed 'condition' | `condition` |
| `CP 2.5` | Group(s) GM12878, K562 appear to have samples derived from the same biological source (e.g., technical replicates). | `Yes — they are independent biological replicates` |
| `CP 2.6` | Experimental design ready: | `Yes — proceed` |
| `CP 2` | ARIA Analysis Plan: | `Run recommended + optional supported analyses` |
| `CP 5` | Analysis complete. | `Accept report` |

## 4. Deterministic Pipeline Execution (`RunLedger`)

Total recorded `RunLedger` nodes: **27**.

| Node ID | Modality | Analysis Stage | Label | Status |
|---|---|---|---|---|
| `ledger://chromatin/qc` | `chromatin` | `qc` | Quality control | **`ran`** |
| `ledger://chromatin/scatac_preprocessing` | `chromatin` | `scatac_preprocessing` | scATAC raw preprocessing & matrix construction | **`not_run`** |
| `ledger://chromatin/bulk_atac_preprocessing` | `chromatin` | `bulk_atac_preprocessing` | Raw bulk chromatin alignment and filtering | **`not_run`** |
| `ledger://chromatin/dimensionality_reduction` | `chromatin` | `dimensionality_reduction` | Dimensionality reduction (LSI/TF-IDF) | **`not_run`** |
| `ledger://chromatin/peak_calling` | `chromatin` | `peak_calling` | Peak calling | **`ran`** |
| `ledger://chromatin/peak_count_matrix` | `chromatin` | `peak_count_matrix` | Peak-count matrix | **`ran`** |
| `ledger://chromatin/differential_accessibility` | `chromatin` | `differential_accessibility` | Differential accessibility / binding | **`ran`** |
| `ledger://chromatin/motif_enrichment` | `chromatin` | `motif_enrichment` | TF motif enrichment | **`ran`** |
| `ledger://chromatin/regulatory_layers` | `chromatin` | `regulatory_layers` | Regulatory layers | **`not_run`** |
| `ledger://chromatin/matrix_pipeline` | `chromatin` | `matrix_pipeline` | scATAC matrix pipeline | **`not_run`** |
| `ledger://chromatin/peak_annotation` | `chromatin` | `peak_annotation` | Genomic peak annotation | **`ran`** |
| `ledger://chromatin/peak_ora` | `chromatin` | `peak_ora` | DA-peak functional enrichment | **`ran`** |
| `ledger://chromatin/differential_tf_footprinting` | `chromatin` | `differential_tf_footprinting` | Differential TF footprinting | **`not_run`** |
| `ledger://bulk/preprocessing` | `bulk` | `preprocessing` | Raw FASTQ preprocessing & alignment | **`not_run`** |
| `ledger://bulk/qc` | `bulk` | `qc` | Sample QC | **`ran`** |
| `ledger://bulk/differential_expression` | `bulk` | `differential_expression` | Differential expression | **`ran`** |
| `ledger://bulk/pathway_enrichment` | `bulk` | `pathway_enrichment` | Pathway enrichment (ORA/GSEA) | **`ran`** |
| `ledger://integration/preprocessing` | `integration` | `preprocessing` | Multiome raw assembly & cohort union | **`not_run`** |
| `ledger://integration/paired_qc` | `integration` | `paired_qc` | Paired multiome QC | **`not_run`** |
| `ledger://integration/wnn` | `integration` | `wnn` | Weighted Nearest Neighbors (WNN) | **`not_run`** |
| `ledger://integration/mofa` | `integration` | `mofa` | MOFA+ factor decomposition | **`not_run`** |
| `ledger://integration/peak_to_gene` | `integration` | `peak_to_gene` | Cis peak-to-gene linkage | **`not_run`** |
| `ledger://integration/eregulons` | `integration` | `eregulons` | Candidate eRegulon inference | **`not_run`** |
| `ledger://integration/cross_modal_triangulation` | `integration` | `cross_modal_triangulation` | Cross-modal evidence triangulation | **`not_run`** |
| `ledger://cross_modal/concordance` | `cross_modal` | `concordance` | Cross-modal contrast & sample concordance | **`ran`** |
| `ledger://report/executive_summary` | `report` | `executive_summary` | Executive summary | **`ran`** |
| `ledger://report/biological_synthesis` | `report` | `biological_synthesis` | Integrated biological synthesis | **`ran`** |

## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine

- **Total Audited Public Claims:** `20` (`{"descriptive": 4, "associative": 16}`)
- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `pass`
- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `0`

### Representative Compiled Claims

- **`[descriptive]` (`executive_summary`):** ARIA completed the submitted analysis. 12 findings were recorded across all modalities, of which 0 reached HIGH confidence. See the governed claim ledger below for published results.
- **`[descriptive]` (`sample_qc`):** Bulk RNA sample QC evaluated 6 samples with 0 outlier(s) flagged; primary DE retained all samples and sensitivity removed 0.
- **`[descriptive]` (`differential_expression`):** Bulk contrast K562 vs GM12878 had 1833 DE genes.
- **`[associative]` (`pathway_enrichment`):** Pathway enrichment found 80 term(s) for K562 vs GM12878.
- **`[associative]` (`gsea_preranked`):** Preranked GSEA generated ranked pathway support for K562 vs GM12878.
- **`[descriptive]` (`power`):** Approximate bulk RNA power ranged from 29% to 29% across analyzable contrasts.

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `3` total (`3` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `12,158` prompt + `4,422` completion = `16,580` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`4,603,837` bytes, Capsule ZIP=`15,905,759` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 11_cross_modal_concordance --force
```
