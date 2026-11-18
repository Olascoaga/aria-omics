# ARIA Reproducible Public Omics Vignettes (11-Modality Atlas)

Every vignette below documents a complete, unattended end-to-end execution of ARIA v4.8+ on a **100% public, peer-reviewed dataset** using the local **Tantalus** open-weight reagent (`Qwen3.6-35B-A3B-UD-Q4_K_M`, `$0.00` cloud API expenditure), verified via `RO-Crate 1.1` and `s14.reproducibility_capsule.v2`.

| Modality ID | ARIA Lane | Public Accession | Status | Runtime | Claims | Capsule | Vignette |
|---|---|---|---|---:|---:|---|---|
| `01_bulk_rna_matrix` | `bulk_RNA` | `GEO GSE47792` | `completed` | `137.84` | `8` | `verified_portable` | [`01_bulk_rna_matrix.md`](./01_bulk_rna_matrix.md) |
| `02_bulk_rna_raw_fastq` | `bulk_RNA_raw` | `GEO GSE52778 / ENA SRP033351` | `completed` | `1034.71` | `9` | `verified_portable` | [`02_bulk_rna_raw_fastq.md`](./02_bulk_rna_raw_fastq.md) |
| `03_scrna_multidonor` | `scRNA` | `GEO GSE96583 / Bioconductor EH2259` | `completed` | `334.48` | `15` | `verified_portable` | [`03_scrna_multidonor.md`](./03_scrna_multidonor.md) |
| `04_scatac_multidonor` | `scATAC` | `GEO GSE278576` | `completed` | `419.67` | `8` | `verified_portable` | [`04_scatac_multidonor.md`](./04_scatac_multidonor.md) |
| `05_bulk_atac_tobias` | `bulk_ATAC` | `ENCODE ENCSR637XSC, ENCSR483RKN, ENCSR859USB` | `completed` | `11868.3` | `9` | `verified_portable` | [`05_bulk_atac_tobias.md`](./05_bulk_atac_tobias.md) |
| `06_chip_seq` | `ChIP` | `ENCODE ENCSR000AKC, ENCSR000AKP, ENCSR000AKJ, ENCSR000AKY` | `completed` | `2025.37` | `9` | `verified_portable` | [`06_chip_seq.md`](./06_chip_seq.md) |
| `07_cut_and_run_spikein` | `CUT_AND_RUN` | `GEO GSE273126 / ENA PRJNA1140214` | `completed` | `194.79` | `6` | `verified_portable` | [`07_cut_and_run_spikein.md`](./07_cut_and_run_spikein.md) |
| `08_cut_and_tag` | `CUT_AND_TAG` | `GEO GSE124557 / ENA PRJNA512492` | `completed` | `910.01` | `9` | `verified_portable` | [`08_cut_and_tag.md`](./08_cut_and_tag.md) |
| `09_hic_3d_genome` | `HiC` | `GEO GSE345433` | `completed` | `2003.67` | `6` | `verified_portable` | [`09_hic_3d_genome.md`](./09_hic_3d_genome.md) |
| `10_paired_multiome` | `Paired_Multiome` | `GEO GSE278576` | `completed` | `949.35` | `11` | `verified_portable` | [`10_paired_multiome.md`](./10_paired_multiome.md) |
| `11_cross_modal_concordance` | `Cross_Modal_RNA_x_ATAC` | `ENCODE ENCSR843RJV, ENCSR889TRN, ENCSR594NJP, ENCSR596ACL, ENCSR637XSC, ENCSR483RKN, ENCSR859USB` | `completed` | `3082.38` | `20` | `verified_portable` | [`11_cross_modal_concordance.md`](./11_cross_modal_concordance.md) |

## Cross-Model Equivalence & Commercial API Cost Benchmark (`30 / 30` Cells)

Vignette [`12_cross_model_equivalence.md`](./12_cross_model_equivalence.md) evaluates ARIA's model-invariance hypothesis (`Claim 6` / `ADR-011`) across `30` end-to-end `(model × modality)` executions comparing local `Tantalus` (`$0.00`), local WSL system-tool harnesses (`Claude Sonnet 5.5`, `GPT-5.6-Sol`, `Gemini 4 Carbon`, `$0.00` cash spend), and direct commercial REST APIs (`anthropic/claude-sonnet-5-5` `$0.6718` total across `7` workflows; `openai/gpt-5.6-sol` `$0.7890` total across `7` workflows). Full per-cell receipts and live dashboard reside in [`../benchmark_results/cross_model_equivalence/STATUS.md`](../benchmark_results/cross_model_equivalence/STATUS.md).

## 40-Capsule Governed vs. Uncontracted `CodeAct` Benchmark (`ADR-083` Tier 1)

Vignette [`13_governed_vs_uncontracted_benchmark.md`](./13_governed_vs_uncontracted_benchmark.md) evaluates `aria_governed`, `aria_locked_manifest`, `aria_zero_llm_deterministic`, and `uncontracted_codeact_baseline` across `40` deterministic multi-omics benchmark capsules (`5` archetypes × `8` capsules: design/confounding governance, donor pseudobulk vs cell pseudoreplication, epistemic claim & causal licensing, sequencer/modality boundary abstention, and scrubbed-ID `DerivedQueryPlan` relational reasoning), paired with a 20-seed `800`-run sensitivity sweep (`../benchmark_results/v49_comparative_benchmark/comparative_benchmark_v49.json`).
