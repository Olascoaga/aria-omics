# Vignette `07_cut_and_run_spikein` — Human U2OS SCAF1-AID Acute Auxin Degron CUT&RUN with Heterologous E. coli Spike-In

- **ARIA Modality Lane:** `CUT_AND_RUN`
- **Public Accession:** `GEO GSE273126 / ENA PRJNA1140214`
- **Publication / DOI:** GSE273126 — SCAF1-AID Auxin Degron HA CUT&RUN in Human U2OS Cells (hg38 + E. coli K-12 MG1655 spike-in) (`10.1016/j.molcel.2024.09.018`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `194.8 s`

## 1. Dataset Characteristics & Public Provenance

3 Control vs 3 Auxin (1 uM, 6h) biological replicates + 1 Untagged negative control (first 350,000 paired-end reads/library preserving human + heterologous E. coli MG1655 spike-in reads).

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `Auxin_rep1_cutnrun_R1.fastq.gz` | `7,405,159` | `016cec586f37bb6297646357...` |
| `Auxin_rep1_cutnrun_R2.fastq.gz` | `15,526,897` | `b32678583cade997eb098e67...` |
| `Auxin_rep2_cutnrun_R1.fastq.gz` | `7,357,830` | `0f8746adf353a1ef59718284...` |
| `Auxin_rep2_cutnrun_R2.fastq.gz` | `15,741,691` | `cf55e8d33d3696b39cc2dc6d...` |
| `Auxin_rep3_cutnrun_R1.fastq.gz` | `7,203,249` | `6148d6674f0ae7fa0cfabc1d...` |
| `Auxin_rep3_cutnrun_R2.fastq.gz` | `15,284,380` | `4d01e3f58d0a9784ccbcd0e0...` |
| `Control_rep1_cutnrun_R1.fastq.gz` | `9,009,389` | `beaa11c1d51d852ef436a1fd...` |
| `Control_rep1_cutnrun_R2.fastq.gz` | `17,526,287` | `13b132d700beecf0d4e9d66b...` |
| `Control_rep2_cutnrun_R1.fastq.gz` | `7,577,615` | `c7bdd7d248be09a0bc9e5280...` |
| `Control_rep2_cutnrun_R2.fastq.gz` | `16,158,933` | `5b97f84b415a91d3cd823042...` |
| `Control_rep3_cutnrun_R1.fastq.gz` | `7,521,908` | `0ec1f748bfe2aa3d6e28c037...` |
| `Control_rep3_cutnrun_R2.fastq.gz` | `16,179,460` | `fc35eba3a7374478ef76b55c...` |
| `IgG_Untagged_ctrl_cutnrun_R1.fastq.gz` | `7,510,437` | `81d223dc8179a45686ff3e4b...` |
| `IgG_Untagged_ctrl_cutnrun_R2.fastq.gz` | `15,878,285` | `d0c197dac0afe2418faaea4c...` |

**Public Source URLs:**
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/064/SRR29997464/SRR29997464_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/064/SRR29997464/SRR29997464_2.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/062/SRR29997462/SRR29997462_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/062/SRR29997462/SRR29997462_2.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/060/SRR29997460/SRR29997460_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/060/SRR29997460/SRR29997460_2.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/063/SRR29997463/SRR29997463_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/063/SRR29997463/SRR29997463_2.fastq.gz`
- *(+6 additional mate/replicate URLs in receipt JSON)*

## 2. Exact Biological Prompt & Experimental Contrast

```text
Evaluate acute auxin-induced SCAF1-AID degradation effects on chromatin occupancy in human U2OS cells comparing Auxin versus Control across 3 biological CUT&RUN replicates per group plus an Untagged negative control, using heterologous E. coli K-12 MG1655 spike-in normalization on GRCh38 (GEO GSE273126). Run paired-end FASTQ QC, bwa-mem2 human + E. coli spike-in alignment, nucleosomal fragment and IgG calibration QC, SEACR peak calling, spike-in scaled differential binding (Auxin versus Control), motif enrichment, and cis gene annotation.
```

- **Primary Contrast(s):** `[['Auxin', 'Control']]` (`condition_col="condition"`, `replicate_col="replicate"`) 
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
| `CP 3.5` | Quality Audit found dispatch readiness issue(s) that may compromise your results: | `Proceed anyway` |
| `CP 5` | Analysis complete. | `Accept report` |

## 4. Deterministic Pipeline Execution (`RunLedger`)

Total recorded `RunLedger` nodes: **14**.

| Node ID | Modality | Analysis Stage | Label | Status |
|---|---|---|---|---|
| `ledger://chromatin/qc` | `chromatin` | `qc` | Quality control | **`ran`** |
| `ledger://chromatin/scatac_preprocessing` | `chromatin` | `scatac_preprocessing` | scATAC raw preprocessing & matrix construction | **`not_run`** |
| `ledger://chromatin/bulk_atac_preprocessing` | `chromatin` | `bulk_atac_preprocessing` | Raw bulk chromatin alignment and filtering | **`executed`** |
| `ledger://chromatin/dimensionality_reduction` | `chromatin` | `dimensionality_reduction` | Dimensionality reduction (LSI/TF-IDF) | **`not_run`** |
| `ledger://chromatin/peak_calling` | `chromatin` | `peak_calling` | Peak calling | **`ran`** |
| `ledger://chromatin/peak_count_matrix` | `chromatin` | `peak_count_matrix` | Peak-count matrix | **`ran`** |
| `ledger://chromatin/differential_accessibility` | `chromatin` | `differential_accessibility` | Differential accessibility / binding | **`ran`** |
| `ledger://chromatin/motif_enrichment` | `chromatin` | `motif_enrichment` | TF motif enrichment | **`not_run`** |
| `ledger://chromatin/regulatory_layers` | `chromatin` | `regulatory_layers` | Regulatory layers | **`not_run`** |
| `ledger://chromatin/matrix_pipeline` | `chromatin` | `matrix_pipeline` | scATAC matrix pipeline | **`not_run`** |
| `ledger://chromatin/peak_annotation` | `chromatin` | `peak_annotation` | Genomic peak annotation | **`skipped`** |
| `ledger://chromatin/peak_ora` | `chromatin` | `peak_ora` | DA-peak functional enrichment | **`skipped`** |
| `ledger://chromatin/differential_tf_footprinting` | `chromatin` | `differential_tf_footprinting` | Differential TF footprinting | **`not_run`** |
| `ledger://report/executive_summary` | `report` | `executive_summary` | Executive summary | **`ran`** |

## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine

- **Total Audited Public Claims:** `6` (`{"descriptive": 1, "associative": 1, "weak_mechanistic": 4}`)
- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `pass`
- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `0`

### Representative Compiled Claims

- **`[descriptive]` (`executive_summary`):** ARIA completed the submitted analysis. 3 findings were recorded across all modalities, of which 0 reached HIGH confidence. See the governed claim ledger below for published results.
- **`[associative]` (`bulk_atac_preprocessing`):** Raw bulk chromatin read preprocessing trimmed and aligned 7 sample(s) (paired layout), retaining filtered non-mitochondrial alignments for peak calling.
- **`[weak_mechanistic]` (`qc`):** CUT&RUN QC measured 6 sample(s) with mean mitochondrial fraction 0.0.
- **`[weak_mechanistic]` (`peak_calling`):** SEACR peak calling for CUT&RUN identified 29456 peaks.
- **`[weak_mechanistic]` (`peak_count_matrix`):** CUT&RUN peak counting produced a matrix with 36 peaks across 6 samples.
- **`[weak_mechanistic]` (`differential_accessibility`):** CUT&RUN DESeq2 differential binding ran 1 comparison(s) over 6 replicate samples and found 0 differential peaks.

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `3` total (`3` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `2,910` prompt + `459` completion = `3,369` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`147,875` bytes, Capsule ZIP=`421,082` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 07_cut_and_run_spikein --force
```
