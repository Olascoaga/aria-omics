# Vignette `08_cut_and_tag` — Founding Human CUT&Tag Polycomb H3K27me3 Profiling (H1 hESC vs K562 + IgG)

- **ARIA Modality Lane:** `CUT_AND_TAG`
- **Public Accession:** `GEO GSE124557 / ENA PRJNA512492`
- **Publication / DOI:** Kaya-Okur HS, Wu SJ, Codomo CA, Pledger ES, Bryson TD, Henikoff JG, Ahmad K, Henikoff S, Nature Communications 10, 1930 (2019) (`10.1038/s41467-019-09982-5`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `910.0 s`

## 1. Dataset Characteristics & Public Provenance

Full unmodified ENA paired-end FASTQ libraries: 2 H1 hESC vs 2 K562 H3K27me3 biological replicates + matched IgG negative controls.

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `H1_H3K27me3_rep1_cuttag_R1.fastq.gz` | `77,044,488` | `6b92a3104d9e7debfb03d229...` |
| `H1_H3K27me3_rep1_cuttag_R2.fastq.gz` | `80,018,762` | `88f4fcabd955d84665c2d467...` |
| `H1_H3K27me3_rep2_cuttag_R1.fastq.gz` | `78,487,420` | `1200cbcb3f9a4169c68f49f0...` |
| `H1_H3K27me3_rep2_cuttag_R2.fastq.gz` | `82,655,196` | `d5835afaa42bce14e436a15a...` |
| `H1_IgG_ctrl_cuttag_R1.fastq.gz` | `25,672,651` | `85593ab5bd4884b7915b578b...` |
| `H1_IgG_ctrl_cuttag_R2.fastq.gz` | `26,277,923` | `c918ebb697cdbb087f36dd4e...` |
| `K562_H3K27me3_rep1_cuttag_R1.fastq.gz` | `87,346,123` | `2ca99c6bf5ce6325a717b001...` |
| `K562_H3K27me3_rep1_cuttag_R2.fastq.gz` | `89,687,298` | `e38c3e90fae018448583eee5...` |
| `K562_H3K27me3_rep2_cuttag_R1.fastq.gz` | `86,044,983` | `de8332312bfdf69a81587c3a...` |
| `K562_H3K27me3_rep2_cuttag_R2.fastq.gz` | `88,341,165` | `522087d8d2645120200b530f...` |
| `K562_IgG_ctrl_cuttag_R1.fastq.gz` | `38,155,677` | `9c73062d9aa0418a958990c4...` |
| `K562_IgG_ctrl_cuttag_R2.fastq.gz` | `39,282,586` | `d866520697e3bdbeeba523ee...` |

**Public Source URLs:**
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/002/SRR8383482/SRR8383482_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/002/SRR8383482/SRR8383482_2.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/003/SRR8383483/SRR8383483_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/003/SRR8383483/SRR8383483_2.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/000/SRR8383510/SRR8383510_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/000/SRR8383510/SRR8383510_2.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/001/SRR8383511/SRR8383511_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/001/SRR8383511/SRR8383511_2.fastq.gz`
- *(+4 additional mate/replicate URLs in receipt JSON)*

## 2. Exact Biological Prompt & Experimental Contrast

```text
Compare Polycomb H3K27me3 repressive chromatin domains between human K562 erythroleukemia cells and H1 embryonic stem cells using the founding Kaya-Okur et al. 2019 paired-end CUT&Tag FASTQ replicates calibrated against matched IgG negative controls on GRCh38 (GEO GSE124557). Execute fastp QC, Tn5-aware bwa-mem2 alignment, fragment-size and FRiP QC, IgG-calibrated SEACR peak calling, consensus quantification, pyDESeq2 differential binding (K562 versus H1, disclosing n=2 replicate power), and cis gene/pathway annotation.
```

- **Primary Contrast(s):** `[['K562', 'H1']]` (`condition_col="condition"`, `replicate_col="replicate"`) 
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
| `ledger://chromatin/motif_enrichment` | `chromatin` | `motif_enrichment` | TF motif enrichment | **`ran`** |
| `ledger://chromatin/regulatory_layers` | `chromatin` | `regulatory_layers` | Regulatory layers | **`not_run`** |
| `ledger://chromatin/matrix_pipeline` | `chromatin` | `matrix_pipeline` | scATAC matrix pipeline | **`not_run`** |
| `ledger://chromatin/peak_annotation` | `chromatin` | `peak_annotation` | Genomic peak annotation | **`ran`** |
| `ledger://chromatin/peak_ora` | `chromatin` | `peak_ora` | DA-peak functional enrichment | **`ran`** |
| `ledger://chromatin/differential_tf_footprinting` | `chromatin` | `differential_tf_footprinting` | Differential TF footprinting | **`not_run`** |
| `ledger://report/executive_summary` | `report` | `executive_summary` | Executive summary | **`ran`** |

## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine

- **Total Audited Public Claims:** `9` (`{"descriptive": 1, "associative": 4, "weak_mechanistic": 4}`)
- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `pass`
- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `0`

### Representative Compiled Claims

- **`[descriptive]` (`executive_summary`):** ARIA completed the submitted analysis. 3 findings were recorded across all modalities, of which 0 reached HIGH confidence. See the governed claim ledger below for published results.
- **`[associative]` (`bulk_atac_preprocessing`):** Raw bulk chromatin read preprocessing trimmed and aligned 6 sample(s) (paired layout), retaining filtered non-mitochondrial alignments for peak calling.
- **`[weak_mechanistic]` (`qc`):** CUT&TAG QC measured 4 sample(s) with mean mitochondrial fraction 0.0.
- **`[weak_mechanistic]` (`peak_calling`):** SEACR peak calling for CUT&TAG identified 30560 peaks.
- **`[weak_mechanistic]` (`peak_count_matrix`):** CUT&TAG peak counting produced a matrix with 19229 peaks across 4 samples.
- **`[weak_mechanistic]` (`differential_accessibility`):** CUT&TAG DESeq2 differential binding ran 1 comparison(s) over 4 replicate samples and found 1792 differential peaks.

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `4` total (`4` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `9,217` prompt + `5,297` completion = `14,514` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`158,268` bytes, Capsule ZIP=`1,106,616` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 08_cut_and_tag --force
```
