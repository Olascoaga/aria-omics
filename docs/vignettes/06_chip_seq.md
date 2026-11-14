# Vignette `06_chip_seq` — ENCODE Broad Histone H3K27ac ChIP-seq with Matched Input Controls (GM12878 vs K562)

- **ARIA Modality Lane:** `ChIP`
- **Public Accession:** `ENCODE ENCSR000AKC, ENCSR000AKP, ENCSR000AKJ, ENCSR000AKY`
- **Publication / DOI:** ENCODE Project Consortium / Bernstein Lab, Nature 489, 57-74 (2012) (`10.1038/nature11247`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `2025.4 s`

## 1. Dataset Characteristics & Public Provenance

2 GM12878 H3K27ac IP + Input vs 2 K562 H3K27ac IP + Input biological replicates (GRCh38 chr6, chr7, chr11, chr17, chr19).

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `GM12878_H3K27ac_rep1_chip.bam` | `82,445,247` | `2bc7f6174e14ad7a10cb5d27...` |
| `GM12878_H3K27ac_rep2_chip.bam` | `72,657,655` | `92fbb48b4c32969a98a3648a...` |
| `GM12878_Input_ctrl_chip.bam` | `42,121,710` | `53d32b3ec92bbda562ea523b...` |
| `K562_H3K27ac_rep1_chip.bam` | `111,876,247` | `c73fbe46282e8a8d809ff3ae...` |
| `K562_H3K27ac_rep2_chip.bam` | `50,611,631` | `b6773de53f455babafa1ca12...` |
| `K562_Input_ctrl_chip.bam` | `117,699,949` | `03bf17fddfe6bcecacc7bb97...` |

**Public Source URLs:**
- `https://www.encodeproject.org/files/ENCFF804NCH/@@download/ENCFF804NCH.bam`
- `https://www.encodeproject.org/files/ENCFF948GTC/@@download/ENCFF948GTC.bam`
- `https://www.encodeproject.org/files/ENCFF301TVL/@@download/ENCFF301TVL.bam`
- `https://www.encodeproject.org/files/ENCFF879BWC/@@download/ENCFF879BWC.bam`
- `https://www.encodeproject.org/files/ENCFF873ZWP/@@download/ENCFF873ZWP.bam`
- `https://www.encodeproject.org/files/ENCFF829HZY/@@download/ENCFF829HZY.bam`

## 2. Exact Biological Prompt & Experimental Contrast

```text
Identify lineage-specific active enhancer and promoter H3K27ac histone modification regions comparing ENCODE K562 erythroleukemia versus GM12878 B-lymphoblastoid cells using biological replicate ChIP-seq IP BAMs calibrated against matched Input chromatin controls on GRCh38 (ENCSR000AKC, ENCSR000AKP). Execute ChIP-seq complexity and strand cross-correlation QC, Input-controlled broad MACS3 peak calling, consensus peak quantification, replicate-aware pyDESeq2 differential binding (K562 versus GM12878, disclosing n=2 replicate power), TF motif enrichment, and cis gene/pathway ORA.
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
| `CP 3.5` | Quality Audit found dispatch readiness issue(s) that may compromise your results: | `Proceed anyway` |
| `CP 5` | Analysis complete. | `Accept report` |

## 4. Deterministic Pipeline Execution (`RunLedger`)

Total recorded `RunLedger` nodes: **14**.

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
| `ledger://chromatin/differential_tf_footprinting` | `chromatin` | `differential_tf_footprinting` | Differential TF footprinting | **`ran`** |
| `ledger://report/executive_summary` | `report` | `executive_summary` | Executive summary | **`ran`** |

## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine

- **Total Audited Public Claims:** `9` (`{"descriptive": 1, "associative": 8}`)
- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `pass`
- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `0`

### Representative Compiled Claims

- **`[descriptive]` (`executive_summary`):** ARIA completed the submitted analysis. 3 findings were recorded across all modalities, of which 0 reached HIGH confidence. See the governed claim ledger below for published results.
- **`[associative]` (`qc`):** Bulk ATAC QC measured 6 sample(s) with mean mitochondrial fraction 0.0.
- **`[associative]` (`peak_calling`):** MACS3 peak calling for bulk ATAC identified 18814 peaks.
- **`[associative]` (`peak_count_matrix`):** Bulk ATAC peak counting produced a matrix with 10316 peaks across 4 samples.
- **`[associative]` (`differential_accessibility`):** Bulk ATAC DESeq2 differential accessibility ran 1 comparison(s) over 4 replicate samples and found 4643 differentially accessible peaks.
- **`[associative]` (`motif_enrichment`):** TF motif enrichment (JASPAR2024_CORE_vertebrates, hypergeometric) tested 879 motifs and found enriched motifs in 2 of 2 peak group(s).

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `4` total (`4` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `9,791` prompt + `3,437` completion = `13,228` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`226,170` bytes, Capsule ZIP=`1,053,553` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 06_chip_seq --force
```
