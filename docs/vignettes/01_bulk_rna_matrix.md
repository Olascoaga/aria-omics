# Vignette `01_bulk_rna_matrix` — SEQC/MAQC-III Human Reference Transcriptomes (UHRR vs HBRR, BGI Center)

- **ARIA Modality Lane:** `bulk_RNA`
- **Public Accession:** `GEO GSE47792`
- **Publication / DOI:** SEQC/MAQC-III Consortium, Nature Biotechnology 32, 903-914 (2014) (`10.1038/nbt.2957`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `137.8 s`

## 1. Dataset Characteristics & Public Provenance

4 Universal Human Reference RNA (UHRR) vs 4 Human Brain Reference RNA (HBRR) replicates (ERCC synthetic spike-ins excluded).

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `seqc_uhrr_vs_hbrr_counts.tsv` | `900,008` | `729735ace9d40424e4dd62ce...` |

**Public Source URLs:**
- `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47792`

## 2. Exact Biological Prompt & Experimental Contrast

```text
Compare Human Brain Reference RNA (HBRR) versus Universal Human Reference RNA (UHRR) across the four SEQC/MAQC-III BGI biological replicates per group (GEO GSE47792) on Homo sapiens (hg38). Perform library QC, pyDESeq2 differential expression (HBRR versus UHRR), ORA and prerank GSEA pathway enrichment, and report transcriptome-wide concordance without causal extrapolation.
```

- **Primary Contrast(s):** `[['HBRR', 'UHRR']]` (`condition_col="condition"`, `replicate_col="replicate"`) 
- **Execution Flags:** `reproducible_mode=True`, `enable_hypotheses=True`

## 3. Interactive Checkpoint Audit Trail

| Checkpoint | Prompt Summary | Chosen Option |
|---|---|---|
| `CP 1` | 📁 ARIA found the following data: | `Confirm and continue` |
| `CP 2.1` | ARIA proposes the following biological groups: | `Yes — confirm groups` |
| `CP 2.2` | Please confirm the organism and genome assembly: | `Homo sapiens (hg38)` |
| `CP 2.3` | Main experimental factor: proposed 'condition' | `condition` |
| `CP 2.5` | Group(s) UHRR, HBRR appear to have samples derived from the same biological source (e.g., technical replicates). | `Yes — they are independent biological replicates` |
| `CP 2.6` | Experimental design ready: | `Yes — proceed` |
| `CP 2` | ARIA Analysis Plan: | `Run recommended + optional supported analyses` |
| `CP 5` | Analysis complete. | `Accept report` |

## 4. Deterministic Pipeline Execution (`RunLedger`)

Total recorded `RunLedger` nodes: **6**.

| Node ID | Modality | Analysis Stage | Label | Status |
|---|---|---|---|---|
| `ledger://bulk/preprocessing` | `bulk` | `preprocessing` | Raw FASTQ preprocessing & alignment | **`not_run`** |
| `ledger://bulk/qc` | `bulk` | `qc` | Sample QC | **`ran`** |
| `ledger://bulk/differential_expression` | `bulk` | `differential_expression` | Differential expression | **`ran`** |
| `ledger://bulk/pathway_enrichment` | `bulk` | `pathway_enrichment` | Pathway enrichment (ORA/GSEA) | **`ran`** |
| `ledger://report/executive_summary` | `report` | `executive_summary` | Executive summary | **`ran`** |
| `ledger://report/biological_synthesis` | `report` | `biological_synthesis` | Integrated biological synthesis | **`ran`** |

## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine

- **Total Audited Public Claims:** `8` (`{"descriptive": 3, "associative": 5}`)
- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `pass`
- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `0`

### Representative Compiled Claims

- **`[descriptive]` (`executive_summary`):** ARIA completed the submitted analysis. 8 findings were recorded across all modalities, of which 0 reached HIGH confidence. See the governed claim ledger below for published results.
- **`[descriptive]` (`sample_qc`):** Bulk RNA sample QC evaluated 8 samples with 0 outlier(s) flagged; primary DE retained all samples and sensitivity removed 0.
- **`[associative]` (`differential_expression`):** Bulk contrast HBRR vs UHRR had 9863 DE genes.
- **`[associative]` (`pathway_enrichment`):** Pathway enrichment found 80 term(s) for HBRR vs UHRR.
- **`[associative]` (`gsea_preranked`):** Preranked GSEA generated ranked pathway support for HBRR vs UHRR.
- **`[descriptive]` (`power`):** Approximate bulk RNA power ranged from 57% to 57% across analyzable contrasts.

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `3` total (`3` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `4,720` prompt + `3,078` completion = `7,798` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`5,546,035` bytes, Capsule ZIP=`5,471,822` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 01_bulk_rna_matrix --force
```
