# Vignette `03_scrna_multidonor` — Human Lupus PBMC IFN-beta Stimulation Across 8 Donors

- **ARIA Modality Lane:** `scRNA`
- **Public Accession:** `GEO GSE96583 / Bioconductor EH2259`
- **Publication / DOI:** Kang HM et al., Nature Biotechnology 36, 89-94 (2018) (`10.1038/nbt.4042`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `334.5 s`

## 1. Dataset Characteristics & Public Provenance

8 matched human donors x 2 conditions (control vs ifnb_stimulated, 16 donor-condition libraries, 6,400 singlet PBMCs).

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `kang2018_pbmc_8donors_ifnb.h5ad` | `31,752,496` | `0b6dd0170be9907694c770f6...` |

**Public Source URLs:**
- `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583`
- `https://experimenthub.bioconductor.org/fetch/2292`

## 2. Exact Biological Prompt & Experimental Contrast

```text
Characterize single-cell transcriptional states and interferon-beta (ifnb_stimulated versus control) responses across 8 matched human lupus PBMC donors (GEO GSE96583, Homo sapiens hg38). Perform single-cell QC, doublet detection, donor-aware Harmony batch integration, Leiden clustering, cell-type annotation, donor-level pseudobulk pyDESeq2 differential expression comparing ifnb_stimulated versus control, and pathway enrichment.
```

- **Primary Contrast(s):** `[['ifnb_stimulated', 'control']]` (`condition_col="condition"`, `replicate_col="donor_id"`) 
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
| `CP 3` | 🔬 Parameter: [resolution] for leiden_clustering | `Use recommended (resolution=0.3)` |
| `CP 5` | Analysis complete. | `Accept report` |

## 4. Deterministic Pipeline Execution (`RunLedger`)

Total recorded `RunLedger` nodes: **12**.

| Node ID | Modality | Analysis Stage | Label | Status |
|---|---|---|---|---|
| `ledger://scRNA/qc` | `scRNA` | `qc` | Quality control | **`ran`** |
| `ledger://scRNA/clustering` | `scRNA` | `clustering` | Clustering | **`ran`** |
| `ledger://scRNA/annotation` | `scRNA` | `annotation` | Cell-type annotation | **`ran`** |
| `ledger://scRNA/marker_discovery` | `scRNA` | `marker_discovery` | Per-cluster marker discovery | **`ran`** |
| `ledger://scRNA/pseudobulk_de` | `scRNA` | `pseudobulk_de` | Pseudobulk differential expression | **`ran`** |
| `ledger://scRNA/differential_abundance` | `scRNA` | `differential_abundance` | Differential abundance | **`ran`** |
| `ledger://scRNA/pathway_enrichment` | `scRNA` | `pathway_enrichment` | Pathway enrichment (ORA/GSEA) | **`ran`** |
| `ledger://scRNA/cell_communication` | `scRNA` | `cell_communication` | Cell-cell communication | **`not_run`** |
| `ledger://scRNA/trajectory` | `scRNA` | `trajectory` | Trajectory (PAGA/DPT) | **`not_run`** |
| `ledger://scRNA/rna_velocity` | `scRNA` | `rna_velocity` | RNA velocity | **`not_run`** |
| `ledger://report/executive_summary` | `report` | `executive_summary` | Executive summary | **`ran`** |
| `ledger://report/biological_synthesis` | `report` | `biological_synthesis` | Integrated biological synthesis | **`ran`** |

## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine

- **Total Audited Public Claims:** `15` (`{"descriptive": 7, "associative": 8}`)
- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `pass`
- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `0`

### Representative Compiled Claims

- **`[descriptive]` (`executive_summary`):** ARIA completed the submitted analysis. 7 findings were recorded across all modalities, of which 0 reached HIGH confidence. See the governed claim ledger below for published results.
- **`[descriptive]` (`qc`):** scRNA QC retained 6,015 of 6,317 cells.
- **`[descriptive]` (`marker_discovery`):** Per-cluster marker discovery ranked candidate marker genes across 8 cluster clusters (Wilcoxon cluster-vs-rest, descriptive ranking only).
- **`[associative]` (`differential_abundance`):** 4 cell type(s) shifted in abundance (FDR < 0.1) for ifnb_stimulated_vs_control.
- **`[descriptive]` (`pseudobulk_de`):** B cells ifnb_stimulated_vs_control had 234 per-cluster FDR DE genes.
- **`[descriptive]` (`pseudobulk_de`):** CD16+ NK cells ifnb_stimulated_vs_control had 170 per-cluster FDR DE genes.

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `6` total (`6` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `16,114` prompt + `6,917` completion = `23,031` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`4,551,104` bytes, Capsule ZIP=`6,318,110` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 03_scrna_multidonor --force
```
