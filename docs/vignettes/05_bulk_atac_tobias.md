# Vignette `05_bulk_atac_tobias` — ENCODE Human Hematopoietic Chromatin Accessibility & TF Footprinting (GM12878 vs K562)

- **ARIA Modality Lane:** `bulk_ATAC`
- **Public Accession:** `ENCODE ENCSR637XSC, ENCSR483RKN, ENCSR859USB`
- **Publication / DOI:** ENCODE Project Consortium, Nature 583, 699-710 (2020) (`10.1038/s41586-020-2493-4`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `11868.3 s`

## 1. Dataset Characteristics & Public Provenance

3 GM12878 B-lymphoblastoid replicates vs 3 K562 erythroleukemia replicates (GRCh38 chr6, chr7, chr11, chr17, chr19) enabling inferential pyDESeq2 DA + TOBIAS differential footprinting.

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `GM12878_rep1_atac.bam` | `2,110,984,885` | `786eb09af7c40c84be508d23...` |
| `GM12878_rep2_atac.bam` | `1,759,054,627` | `7f7eb74403ff195f8cd96d02...` |
| `GM12878_rep3_atac.bam` | `2,134,120,685` | `32d1685aa2ba66808343b850...` |
| `K562_rep1_atac.bam` | `840,054,452` | `117a1af6cdb1e92177b0ab0b...` |
| `K562_rep2_atac.bam` | `792,289,628` | `a2363caa9b1a92d85a724e7d...` |
| `K562_rep3_atac.bam` | `173,245,367` | `8fd78d5bb93313d2afae8d73...` |

**Public Source URLs:**
- `https://www.encodeproject.org/files/ENCFF981FXV/@@download/ENCFF981FXV.bam`
- `https://www.encodeproject.org/files/ENCFF962FMH/@@download/ENCFF962FMH.bam`
- `https://www.encodeproject.org/files/ENCFF440GRZ/@@download/ENCFF440GRZ.bam`
- `https://www.encodeproject.org/files/ENCFF512VEZ/@@download/ENCFF512VEZ.bam`
- `https://www.encodeproject.org/files/ENCFF987XOV/@@download/ENCFF987XOV.bam`
- `https://www.encodeproject.org/files/ENCFF557RSA/@@download/ENCFF557RSA.bam`

## 2. Exact Biological Prompt & Experimental Contrast

```text
Compare chromatin accessibility and transcription factor occupancy footprints between ENCODE K562 erythroleukemia cells and GM12878 B-lymphoblastoid cells across 3 biological replicates per lineage on GRCh38 (ENCSR637XSC, ENCSR483RKN, ENCSR859USB). Run bulk ATAC QC, MACS3 narrowPeak calling, consensus peak quantification, replicate-aware pyDESeq2 differential accessibility (K562 versus GM12878), JASPAR TF motif enrichment, cis gene annotation, pathway ORA, and TOBIAS condition-level differential TF footprinting.
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
- **`[associative]` (`peak_calling`):** MACS3 peak calling for bulk ATAC identified 88238 peaks.
- **`[associative]` (`peak_count_matrix`):** Bulk ATAC peak counting produced a matrix with 47022 peaks across 6 samples.
- **`[associative]` (`differential_accessibility`):** Bulk ATAC DESeq2 differential accessibility ran 1 comparison(s) over 6 replicate samples and found 13706 differentially accessible peaks.
- **`[associative]` (`motif_enrichment`):** TF motif enrichment (JASPAR2024_CORE_vertebrates, hypergeometric) tested 879 motifs and found enriched motifs in 2 of 2 peak group(s).

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `4` total (`4` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `10,321` prompt + `5,009` completion = `15,330` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`5,604,032` bytes, Capsule ZIP=`5,756,020` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 05_bulk_atac_tobias --force
```
