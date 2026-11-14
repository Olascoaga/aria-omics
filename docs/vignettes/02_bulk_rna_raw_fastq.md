# Vignette `02_bulk_rna_raw_fastq` — Human Airway Smooth Muscle Glucocorticoid (Dexamethasone) Response

- **ARIA Modality Lane:** `bulk_RNA_raw`
- **Public Accession:** `GEO GSE52778 / ENA SRP033351`
- **Publication / DOI:** Himes BE et al., PLOS ONE 9(6): e99625 (2014) (`10.1371/journal.pone.0099625`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `1034.7 s`

## 1. Dataset Characteristics & Public Provenance

4 untreated vs 4 dexamethasone-treated primary human airway smooth muscle donor lines (first 250,000 paired-end reads per ENA run).

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `dexamethasone_rep1_R1.fastq.gz` | `13,149,559` | `3611e84670f2cd23c8c78651...` |
| `dexamethasone_rep1_R2.fastq.gz` | `13,229,210` | `bd183ca98113a05440f7e30d...` |
| `dexamethasone_rep2_R1.fastq.gz` | `13,183,802` | `6880ae8f31d05160638e7948...` |
| `dexamethasone_rep2_R2.fastq.gz` | `13,275,526` | `7ca5b4d86afcf325816f4152...` |
| `dexamethasone_rep3_R1.fastq.gz` | `14,341,373` | `6d617ee92a4ed85991991594...` |
| `dexamethasone_rep3_R2.fastq.gz` | `14,461,397` | `4f36ef81542ff6ee907a46b4...` |
| `dexamethasone_rep4_R1.fastq.gz` | `13,237,122` | `bb851b0c1a14a6a687cf324a...` |
| `dexamethasone_rep4_R2.fastq.gz` | `13,326,759` | `c2a41b371eec3489532d2596...` |
| `untreated_rep1_R1.fastq.gz` | `13,153,585` | `5a3de7a4efadee493cb4ef8d...` |
| `untreated_rep1_R2.fastq.gz` | `13,032,053` | `f5a4e6acc497c12b78292800...` |
| `untreated_rep2_R1.fastq.gz` | `13,735,329` | `bc7b1821771218564bf3dc6c...` |
| `untreated_rep2_R2.fastq.gz` | `13,879,452` | `1871199ece7a2e6349d1cbfe...` |
| `untreated_rep3_R1.fastq.gz` | `13,060,742` | `5fbee1a2c1c433505e12970b...` |
| `untreated_rep3_R2.fastq.gz` | `13,225,146` | `8504958d380667e40cec5024...` |
| `untreated_rep4_R1.fastq.gz` | `13,212,798` | `b21425d623094a5dbe80df66...` |
| `untreated_rep4_R2.fastq.gz` | `13,301,174` | `e7cca50406548895acf2a5de...` |

**Public Source URLs:**
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_2.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039512/SRR1039512_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039512/SRR1039512_2.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/006/SRR1039516/SRR1039516_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/006/SRR1039516/SRR1039516_2.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039520/SRR1039520_1.fastq.gz`
- `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039520/SRR1039520_2.fastq.gz`
- *(+8 additional mate/replicate URLs in receipt JSON)*

## 2. Exact Biological Prompt & Experimental Contrast

```text
Analyze the paired-end raw RNA-seq FASTQ libraries from primary human airway smooth muscle cells treated with dexamethasone versus untreated controls across four biological donor lines (GEO GSE52778 / ENA SRP033351, Homo sapiens hg38). Execute fastp read trimming and QC, STAR splice-aware alignment, featureCounts gene quantification, pyDESeq2 differential expression comparing dexamethasone versus untreated, and functional pathway enrichment.
```

- **Primary Contrast(s):** `[['dexamethasone', 'untreated']]` (`condition_col="condition"`, `replicate_col="replicate"`) 
- **Execution Flags:** `reproducible_mode=True`, `enable_hypotheses=True`

## 3. Interactive Checkpoint Audit Trail

| Checkpoint | Prompt Summary | Chosen Option |
|---|---|---|
| `CP 1` | 📁 ARIA found the following data: | `Confirm and continue` |
| `CP 2.1` | ARIA proposes the following biological groups: | `Yes — confirm groups` |
| `CP 2.2` | Please confirm the organism and genome assembly: | `Homo sapiens (hg38)` |
| `CP 2.3` | Main experimental factor: proposed 'condition' | `condition` |
| `CP 2.5` | Group(s) untreated, dexamethasone appear to have samples derived from the same biological source (e.g., technical replicates). | `Yes — they are independent biological replicates` |
| `CP 2.6` | Experimental design ready: | `Yes — proceed` |
| `CP 2` | ARIA Analysis Plan: | `Run recommended + optional supported analyses` |
| `CP 5` | Analysis complete. | `Accept report` |

## 4. Deterministic Pipeline Execution (`RunLedger`)

Total recorded `RunLedger` nodes: **13**.

| Node ID | Modality | Analysis Stage | Label | Status |
|---|---|---|---|---|
| `ledger://bulk/preprocessing` | `bulk` | `preprocessing` | Raw FASTQ preprocessing & alignment | **`ran`** |
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
| `ledger://report/executive_summary` | `report` | `executive_summary` | Executive summary | **`ran`** |
| `ledger://report/biological_synthesis` | `report` | `biological_synthesis` | Integrated biological synthesis | **`ran`** |

## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine

- **Total Audited Public Claims:** `9` (`{"descriptive": 4, "associative": 5}`)
- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `pass`
- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `0`

### Representative Compiled Claims

- **`[descriptive]` (`executive_summary`):** ARIA completed the submitted analysis. 9 findings were recorded across all modalities, of which 0 reached HIGH confidence. See the governed claim ledger below for published results.
- **`[descriptive]` (`preprocessing`):** Raw FASTQ preprocessing (full_length) quantified 16551 genes across 8 samples.
- **`[descriptive]` (`sample_qc`):** Bulk RNA sample QC evaluated 8 samples with 0 outlier(s) flagged; primary DE retained all samples and sensitivity removed 0.
- **`[associative]` (`differential_expression`):** Bulk contrast dexamethasone vs untreated had 18 DE genes.
- **`[associative]` (`pathway_enrichment`):** Pathway enrichment found 12 term(s) for dexamethasone vs untreated.
- **`[associative]` (`gsea_preranked`):** Preranked GSEA generated ranked pathway support for dexamethasone vs untreated.

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `3` total (`3` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `4,017` prompt + `2,795` completion = `6,812` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`2,381,092` bytes, Capsule ZIP=`2,345,324` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 02_bulk_rna_raw_fastq --force
```
