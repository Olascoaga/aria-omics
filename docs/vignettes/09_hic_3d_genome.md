# Vignette `09_hic_3d_genome` — Human RPE-1 Whole-Genome Duplication 3D Chromatin Organization (Diploid vs Tetraploid G1)

- **ARIA Modality Lane:** `HiC`
- **Public Accession:** `GEO GSE345433`
- **Publication / DOI:** GSE345433 — Whole genome duplication through mitotic slippage 3D Hi-C in human RPE-1 cells (GRCh38 .mcool) (`NCBI GEO GSE345433`)
- **Organism / Reference:** *Homo sapiens* (`hg38`)
- **Local LLM Reagent:** `qwen3.6-35b-a3b-ud-q4_k_m` via `tantalus` (`$0.00` cloud API cost)
- **End-to-End Status:** `completed` (`verify_reproducible_capsule == "verified_portable"`)
- **Wall-Clock Runtime:** `2003.7 s`

## 1. Dataset Characteristics & Public Provenance

3 Diploid_G1 vs 3 Tetraploid_G1 biological replicates (multi-resolution GRCh38 .mcool matrices).

| Payload File | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `Diploid_rep1.mcool` | `445,936,960` | `e6c1711a604e308f5d01ba4e...` |
| `Diploid_rep2.mcool` | `490,297,183` | `53b2bd3ab1ed2dab44ab53f0...` |
| `Diploid_rep3.mcool` | `438,680,166` | `020488272b6809b653766d2d...` |
| `Tetraploid_rep1.mcool` | `456,144,325` | `0125d5711c9605ec7b37cc4f...` |
| `Tetraploid_rep2.mcool` | `120,133,024` | `0289c4c0bbd3b3326f227191...` |
| `Tetraploid_rep3.mcool` | `468,928,569` | `e0328166b50e445ccdcb1b6f...` |

**Public Source URLs:**
- `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM10003nnn/GSM10003420/suppl/GSM10003420_RPE_2N-MS02_hic_1.mcool`
- `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM10003nnn/GSM10003421/suppl/GSM10003421_RPE_2N-MS02_hic_2.mcool`
- `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM10003nnn/GSM10003422/suppl/GSM10003422_RPE_2N-MS1_hic_1.mcool`
- `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM10003nnn/GSM10003424/suppl/GSM10003424_RPE_4N-MS02_hic_1.mcool`
- `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM10003nnn/GSM10003425/suppl/GSM10003425_RPE_4N-MS02_hic_2.mcool`
- `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM10003nnn/GSM10003426/suppl/GSM10003426_RPE_4N-MS1_hic_1.mcool`

## 2. Exact Biological Prompt & Experimental Contrast

```text
Assess 3D genome reorganization following whole-genome duplication in human RPE-1 G1 cells comparing Tetraploid versus Diploid across 3 biological Hi-C replicates per ploidy state on GRCh38 (GEO GSE345433). Perform ICE contact balancing, P(s) distance-decay polymer scaling and HiCRep stratum-adjusted correlation QC, GC-phased A/B compartment decomposition, Crane multi-scale insulation TAD boundary calling, chromatin loop calling, and replicate-aware differential compartment (Delta E1) and insulation shifts.
```

- **Primary Contrast(s):** `[['Tetraploid', 'Diploid']]` (`condition_col="condition"`, `replicate_col="replicate"`) 
- **Execution Flags:** `reproducible_mode=True`, `enable_hypotheses=True`

## 3. Interactive Checkpoint Audit Trail

| Checkpoint | Prompt Summary | Chosen Option |
|---|---|---|
| `CP 1` | 📁 ARIA found the following data: | `Confirm and continue` |
| `CP 2.1` | ARIA proposes the following biological groups: | `Yes — confirm groups` |
| `CP 2.2` | Please confirm the organism and genome assembly: | `Homo sapiens (hg38)` |
| `CP 2.3` | Main experimental factor: proposed 'condition' | `condition` |
| `CP 2.5` | Group(s) Diploid, Tetraploid appear to have samples derived from the same biological source (e.g., technical replicates). | `Yes — they are independent biological replicates` |
| `CP 2.6` | Experimental design ready: | `Yes — proceed` |
| `CP 2` | ARIA Analysis Plan: | `Run recommended + optional supported analyses` |
| `CP 2` | Hi-C Resolution Selection | `Use recommended resolution (50,000bp)` |
| `CP 5` | Analysis complete. | `Accept report` |
| `CP 3` | Insulation Score Window Size (resolution: 50,000bp) | `Use recommended window (500,000bp)` |
| `CP 5` | Analysis complete. | `Accept report` |
| `CP 5` | Analysis complete. | `Accept report` |

## 4. Deterministic Pipeline Execution (`RunLedger`)

Total recorded `RunLedger` nodes: **7**.

| Node ID | Modality | Analysis Stage | Label | Status |
|---|---|---|---|---|
| `ledger://hic/preprocessing` | `hic` | `preprocessing` | 3D contact preprocessing & cooler build | **`not_run`** |
| `ledger://hic/qc` | `hic` | `qc` | Hi-C contact QC, ICE balancing & P(s)/SCC | **`ran`** |
| `ledger://hic/compartments` | `hic` | `compartments` | A/B compartment decomposition (E1) | **`ran`** |
| `ledger://hic/tads` | `hic` | `tads` | Topologically associating domains (insulation) | **`ran`** |
| `ledger://hic/loops` | `hic` | `loops` | Chromatin loop calling (donut enrichment) | **`ran`** |
| `ledger://hic/differential_3d` | `hic` | `differential_3d` | Differential 3D compartment & insulation | **`ran`** |
| `ledger://report/executive_summary` | `report` | `executive_summary` | Executive summary | **`ran`** |

## 5. Compiled Public Claims (`C1–C4`) & Epistemic Quarantine

- **Total Audited Public Claims:** `6` (`{"descriptive": 4, "associative": 2}`)
- **Speculative Promotion Wall (`assert_no_speculative_promotion`):** `pass`
- **Quarantined `SPECULATIVE` Hypotheses (`S0–S10`):** `0`

### Representative Compiled Claims

- **`[descriptive]` (`executive_summary`):** ARIA completed the submitted analysis. 6 findings were recorded across all modalities, of which 0 reached HIGH confidence. See the governed claim ledger below for published results.
- **`[descriptive]` (`qc`):** Hi-C contact QC and ICE balancing evaluated 6 sample(s) with 463783660 valid contact pairs and cis/trans ratio 0.782.
- **`[descriptive]` (`compartments`):** Hi-C cis eigenvector compartment decomposition at 100000 bp resolution classified 48.21% active and 51.79% inactive compartment bins.
- **`[descriptive]` (`tads`):** Hi-C diamond insulation score analysis identified 568 TAD intervals with median size 1150.0 kb.
- **`[associative]` (`loops`):** Hi-C donut local background enrichment at 10000 bp resolution identified 0 chromatin loop(s).
- **`[associative]` (`differential_3d`):** Hi-C replicate-level differential topology testing across 1 contrast(s) identified 0 differential compartment bin(s), including 0 compartment switch(es), and 2 differential insulation bin(s).

## 6. Local LLM Reagent Accounting & Reproducibility Capsule

- **Served Model:** `qwen3.6-35b-a3b-ud-q4_k_m` (`Qwen3.6-35B-A3B MoE (35B total / 3B active)`)
- **LLM Calls:** `5` total (`5` live, `0` cache hits, `0` fallbacks)
- **Tokens Processed:** `4,966` prompt + `750` completion = `5,716` total tokens
- **Cloud API Cost:** **`$0.00`** (`zero_cost_audit_passed=True`)
- **Capsule Verification:** `verified_portable` (`s14.reproducibility_capsule.v2`, RO-Crate=`True`, SQLite Snapshot=`True`, HTML=`1,434,983` bytes, Capsule ZIP=`2,416,042` bytes)

## 7. One-Command Reproduction

```bash
conda run -n aria-env python scripts/run_public_modality_atlas.py --modality 09_hic_3d_genome --force
```
