# Vignette `12_cross_model_equivalence` — Cross-Model Invariance & Commercial API Cost Benchmark

- **Updated (UTC):** `2026-11-16T05:30:11Z`
- **Hypothesis Tested (`Claim 6` / `ADR-011`):** Because ARIA confines LLMs to experimental framing, database-backed label interpretation, and narrative synthesis while audited code executes every statistical test (`pyDESeq2`, `CellTypist`, `LSI`, `WNN`, `Peak2Gene`, `eRegulons`, `MOFA+`, `SEACR`, `cooltools`) and `C1–C4` compiles claims deterministically, **scientific conclusions and epistemic boundaries remain invariant across local open-weight models (`Tantalus`, `$0.00`) and commercial frontier models (`Claude Sonnet 5.5`, `GPT-5.6-Sol`, `Gemini 4 Carbon`).**
- **Cumulative Commercial API Spend:** Anthropic `claude-sonnet-5-5` = **`$0.6718 USD`** · OpenAI `gpt-5.6-sol` = **`$0.7890 USD`** · Local `Tantalus` = **`$0.00 USD`**

## 1. Stage 2A — Local WSL System-Tool Equivalence (`03_scrna`, `04_scatac`, `10_paired_multiome`, `$0.00` Cash Spend)

| Model Family | Transport | Modality | Public Claims | Net ARIA Tokens | Projected List Cost | CLI Raw Input Tokens | Stage Jaccard vs Tantalus | Stat Claim Jaccard vs Tantalus | Table Row Match |
|---|---|---|---:|---:|---:|---:|---:|---:|---|
| `Claude Sonnet 5.5 (WSL System Tool)` | `wsl_cli_shim` | `03_scrna_multidonor` | `15` | `17,659` | `$0.0699` | `22,077` | `1.00` | `0.17` | `1/5` |
| `Claude Sonnet 5.5 (WSL System Tool)` | `wsl_cli_shim` | `04_scatac_multidonor` | `8` | `10,338` | `$0.0391` | `12,149` | `1.00` | `1.00` | `no_bundled_tables` |
| `Claude Sonnet 5.5 (WSL System Tool)` | `wsl_cli_shim` | `10_paired_multiome` | `15` | `43,324` | `$0.1722` | `51,145` | `1.00` | `1.00` | `4/4` |
| `GPT-5.6-Sol (Codex WSL System Tool)` | `wsl_cli_shim` | `03_scrna_multidonor` | `15` | `17,916` | `$0.1297` | `99,622` | `1.00` | `1.00` | `5/5` |
| `GPT-5.6-Sol (Codex WSL System Tool)` | `wsl_cli_shim` | `04_scatac_multidonor` | `8` | `10,096` | `$0.0736` | `65,289` | `1.00` | `1.00` | `no_bundled_tables` |
| `GPT-5.6-Sol (Codex WSL System Tool)` | `wsl_cli_shim` | `10_paired_multiome` | `15` | `35,924` | `$0.2277` | `216,383` | `1.00` | `1.00` | `4/4` |
| `Gemini 4 Carbon Low (Antigravity WSL Tool)` | `wsl_cli_shim` | `03_scrna_multidonor` | `15` | `17,753` | `$0.0614` | `70,313` | `1.00` | `0.17` | `1/5` |
| `Gemini 4 Carbon Low (Antigravity WSL Tool)` | `wsl_cli_shim` | `04_scatac_multidonor` | `8` | `11,306` | `$0.0428` | `45,986` | `1.00` | `1.00` | `no_bundled_tables` |
| `Gemini 4 Carbon Low (Antigravity WSL Tool)` | `wsl_cli_shim` | `10_paired_multiome` | `15` | `49,497` | `$0.1699` | `183,034` | `1.00` | `1.00` | `4/4` |
| `Tantalus Local (Qwen3.6-35B-A3B Q4_K_M)` | `local_tantalus_gguf` | `03_scrna_multidonor` | `15` | `23,031` | `$0.0000` | `—` | `1.00` | `1.00` | `5/5` |
| `Tantalus Local (Qwen3.6-35B-A3B Q4_K_M)` | `local_tantalus_gguf` | `04_scatac_multidonor` | `8` | `3,568` | `$0.0000` | `—` | `1.00` | `1.00` | `no_bundled_tables` |
| `Tantalus Local (Qwen3.6-35B-A3B Q4_K_M)` | `local_tantalus_gguf` | `10_paired_multiome` | `15` | `78,028` | `$0.0000` | `—` | `1.00` | `1.00` | `4/4` |

## 2. Stage 2B — Direct Commercial REST API Cost & Multi-Modality Scaling

| Commercial Provider / Model | Modality | Runtime | Calls | Prompt Tokens | Completion Tokens | Billed Cost (USD) | Stage Jaccard vs Tantalus | Stat Claim Jaccard vs Tantalus | Speculative Wall |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `Claude Sonnet 5.5 (Direct Anthropic REST API)` | `01_bulk_rna_matrix` | `61.7s` | `3` | `6,529` | `5,563` | **`$0.0687`** | `1.00` | `1.00` | `pass` |
| `Claude Sonnet 5.5 (Direct Anthropic REST API)` | `03_scrna_multidonor` | `135.1s` | `6` | `20,064` | `8,133` | **`$0.1215`** | `1.00` | `0.17` | `pass` |
| `Claude Sonnet 5.5 (Direct Anthropic REST API)` | `04_scatac_multidonor` | `267.6s` | `4` | `12,136` | `5,793` | **`$0.0822`** | `1.00` | `1.00` | `pass` |
| `Claude Sonnet 5.5 (Direct Anthropic REST API)` | `08_cut_and_tag` | `801.9s` | `4` | `12,244` | `6,594` | **`$0.0904`** | `1.00` | `1.00` | `pass` |
| `Claude Sonnet 5.5 (Direct Anthropic REST API)` | `09_hic_3d_genome` | `1868.5s` | `3` | `4,242` | `1,811` | **`$0.0266`** | `1.00` | `1.00` | `pass` |
| `Claude Sonnet 5.5 (Direct Anthropic REST API)` | `10_paired_multiome` | `672.5s` | `10` | `43,231` | `10,545` | **`$0.1919`** | `1.00` | `1.00` | `pass` |
| `Claude Sonnet 5.5 (Direct Anthropic REST API)` | `11_cross_modal_concordance` | `2863.7s` | `4` | `11,220` | `6,804` | **`$0.0905`** | `0.62` | `0.37` | `pass` |
| `GPT-5.6-Sol (Direct OpenAI REST API)` | `01_bulk_rna_matrix` | `44.1s` | `3` | `4,381` | `2,988` | **`$0.0773`** | `1.00` | `1.00` | `pass` |
| `GPT-5.6-Sol (Direct OpenAI REST API)` | `03_scrna_multidonor` | `121.8s` | `6` | `14,213` | `4,965` | **`$0.1562`** | `1.00` | `1.00` | `pass` |
| `GPT-5.6-Sol (Direct OpenAI REST API)` | `04_scatac_multidonor` | `245.6s` | `4` | `7,925` | `3,149` | **`$0.0947`** | `1.00` | `1.00` | `pass` |
| `GPT-5.6-Sol (Direct OpenAI REST API)` | `08_cut_and_tag` | `775.6s` | `4` | `8,111` | `2,764` | **`$0.0877`** | `1.00` | `1.00` | `pass` |
| `GPT-5.6-Sol (Direct OpenAI REST API)` | `09_hic_3d_genome` | `2062.3s` | `3` | `2,641` | `1,309` | **`$0.0367`** | `1.00` | `1.00` | `pass` |
| `GPT-5.6-Sol (Direct OpenAI REST API)` | `10_paired_multiome` | `651.0s` | `13` | `30,579` | `5,477` | **`$0.2319`** | `1.00` | `1.00` | `pass` |
| `GPT-5.6-Sol (Direct OpenAI REST API)` | `11_cross_modal_concordance` | `3045.1s` | `3` | `10,210` | `3,185` | **`$0.1045`** | `1.00` | `1.00` | `pass` |

## 3. Empirical Findings & Diagnostic Planning Branches

- **Bit-Identical Statistical & Claim Invariance Across 5/7 Modalities (`1.00` Stage & Claim Jaccard):** Across `01_bulk_rna_matrix` (`8` claims), `04_scatac_multidonor` (`8` claims), `08_cut_and_tag` (`9` claims), `09_hic_3d_genome` (`6` claims, `5/5` SHA-256 identical tables), and `10_paired_multiome` (`15` claims across all `6` models, `4/4` SHA-256 identical tables), local `Tantalus` (`$0.00`), local WSL CLI harnesses (`Claude Sonnet 5.5`, `GPT-5.6-Sol`, `Gemini 4 Carbon`), and direct commercial REST APIs produced **100% identical stage topologies and verbatim compiled statistical claims**, with zero `SPECULATIVE` leakage (`speculative_promotion_guard == "pass"` across `30/30` cells).
- **Diagnostic Free-Planning Branch #1 (`03_scrna_multidonor` Leiden `0.2` vs `0.3`):** All `6/6` models achieved `1.00` Stage Jaccard, `15` compiled claims, and identical `4` CellTypist immune lineages (`CD4 T cells`, `Monocytes`, `B cells`, `NK cells`). `Tantalus` and both `GPT-5.6-Sol` transports selected `resolution=0.3` (`1,252` pseudobulk DE genes, `1.00` Claim Jaccard), whereas both `Claude Sonnet 5.5` transports and `Gemini 4 Carbon` selected `resolution=0.2` (`1,251` pseudobulk DE genes, Pearson $r = 0.9998$ across shared gene $\log_2\text{FC}$s) via `ParameterAdvisor._intent_to_leiden_range`.
- **Diagnostic Free-Planning Branch #2 (`11_cross_modal_concordance` Fail-Closed Contrast Schema):** Local `Tantalus` (`$0.00`) and `openai_api_gpt56` (`$0.1045`) achieved **`1.00` Stage & Statistical Claim Jaccard (`20` claims)** across `BulkRNAAgent` + `ChromatinAgent` + `cross_modal_concordance.py`. `claude_api_sonnet55` emitted a string contrast (`"K562 vs GM12878"`) instead of a 2-element list/dict in `_design_analysis_plan`; `BulkRNAAgent._normalise_explicit_contrasts` refused to guess contrast polarity and halted cleanly at `CP bulk.contrast` (`0` `RunLedger` errors) while `ChromatinAgent` completed all `7` bulk ATAC stages (`8` claims, `0.62` Stage Jaccard).
- **Commercial Cost & Dual Token Accounting:** Direct REST API execution cost **`$0.0266–$0.2319 USD` per complete omics workflow** (`$0.6718` total across `7` Anthropic workflows, avg `$0.096`/run; `$0.7890` total across `7` OpenAI workflows, avg `$0.113`/run) versus **`$0.00 USD`** on local `Tantalus`. Stage 2A dual accounting separates net ARIA prompt/completion tokens from external CLI system-prompt wrapper overhead.

## 4. Archived Self-Contained HTML Reports & Reproducibility Capsules

Every `(model × modality)` execution preserved a complete `report.html`, `methodology.json`, and `reproducible_capsule.zip` (`verified_portable` when inputs reside inside the run bundle, `verified_portable_ext_inputs` when raw sequencing inputs are referenced from the shared immutable public-atlas cache) under `~/.aria/idea2_runs/reports/<model_key>/<modality_id>/`.

