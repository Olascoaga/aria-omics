# ARIA Idea 2 — Cross-Model Equivalence & Commercial API Cost Benchmark

- **Last Updated (UTC):** `2026-11-16T05:30:11Z`
- **Completed Cells:** `30`
- **Actual Commercial API Spend:** Anthropic = **`$0.6718 USD`** · OpenAI = **`$0.7890 USD`** (Hard cap: `$4.50` / provider)

## Live Execution Matrix

| Stage | Model | Modality | Status | Runtime | Claims | Net ARIA Tokens | Billed USD | List-Tariff USD | Stage Jaccard vs Tantalus | Stat Claim Jaccard | Capsule |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `stage2b_commercial` | `claude_api_sonnet55` | `01_bulk_rna_matrix` | **`completed`** | `61.7s` | `8` | `12,092` | `$0.0687` | `$0.0687` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `claude_api_sonnet55` | `03_scrna_multidonor` | **`completed`** | `135.1s` | `15` | `28,197` | `$0.1215` | `$0.1215` | `1.00` | `0.17` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `claude_api_sonnet55` | `04_scatac_multidonor` | **`completed`** | `267.6s` | `8` | `17,929` | `$0.0822` | `$0.0822` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `claude_api_sonnet55` | `08_cut_and_tag` | **`completed`** | `801.9s` | `9` | `18,838` | `$0.0904` | `$0.0904` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `claude_api_sonnet55` | `09_hic_3d_genome` | **`completed`** | `1868.5s` | `6` | `6,053` | `$0.0266` | `$0.0266` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `claude_api_sonnet55` | `10_paired_multiome` | **`completed`** | `672.5s` | `15` | `53,776` | `$0.1919` | `$0.1919` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `claude_api_sonnet55` | `11_cross_modal_concordance` | **`completed`** | `2863.7s` | `8` | `18,024` | `$0.0905` | `$0.0905` | `0.62` | `0.37` | `verified_portable_ext_inputs` |
| `stage2a_wsl` | `claude_wsl_sonnet55` | `03_scrna_multidonor` | **`completed`** | `129.1s` | `15` | `17,659` | `$0.0000` | `$0.0699` | `1.00` | `0.17` | `verified_portable_ext_inputs` |
| `stage2a_wsl` | `claude_wsl_sonnet55` | `04_scatac_multidonor` | **`completed`** | `256.7s` | `8` | `10,338` | `$0.0000` | `$0.0391` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2a_wsl` | `claude_wsl_sonnet55` | `10_paired_multiome` | **`completed`** | `735.4s` | `15` | `43,324` | `$0.0000` | `$0.1722` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2a_wsl` | `codex_wsl_gpt56` | `03_scrna_multidonor` | **`completed`** | `170.6s` | `15` | `17,916` | `$0.0000` | `$0.1297` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2a_wsl` | `codex_wsl_gpt56` | `04_scatac_multidonor` | **`completed`** | `269.6s` | `8` | `10,096` | `$0.0000` | `$0.0736` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2a_wsl` | `codex_wsl_gpt56` | `10_paired_multiome` | **`completed`** | `733.8s` | `15` | `35,924` | `$0.0000` | `$0.2277` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2a_wsl` | `gemini_agy_carbon` | `03_scrna_multidonor` | **`completed`** | `160.6s` | `15` | `17,753` | `$0.0000` | `$0.0614` | `1.00` | `0.17` | `verified_portable_ext_inputs` |
| `stage2a_wsl` | `gemini_agy_carbon` | `04_scatac_multidonor` | **`completed`** | `276.4s` | `8` | `11,306` | `$0.0000` | `$0.0428` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2a_wsl` | `gemini_agy_carbon` | `10_paired_multiome` | **`completed`** | `855.1s` | `15` | `49,497` | `$0.0000` | `$0.1699` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `openai_api_gpt56` | `01_bulk_rna_matrix` | **`completed`** | `44.1s` | `8` | `7,369` | `$0.0773` | `$0.0773` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `openai_api_gpt56` | `03_scrna_multidonor` | **`completed`** | `121.8s` | `15` | `19,178` | `$0.1562` | `$0.1562` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `openai_api_gpt56` | `04_scatac_multidonor` | **`completed`** | `245.6s` | `8` | `11,074` | `$0.0947` | `$0.0947` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `openai_api_gpt56` | `08_cut_and_tag` | **`completed`** | `775.6s` | `9` | `10,875` | `$0.0877` | `$0.0877` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `openai_api_gpt56` | `09_hic_3d_genome` | **`completed`** | `2062.3s` | `6` | `3,950` | `$0.0367` | `$0.0367` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `openai_api_gpt56` | `10_paired_multiome` | **`completed`** | `651.0s` | `15` | `36,056` | `$0.2319` | `$0.2319` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `stage2b_commercial` | `openai_api_gpt56` | `11_cross_modal_concordance` | **`completed`** | `3045.1s` | `20` | `13,395` | `$0.1045` | `$0.1045` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `baseline` | `tantalus_local` | `01_bulk_rna_matrix` | **`completed`** | `137.8s` | `8` | `7,798` | `$0.0000` | `$0.0000` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `baseline` | `tantalus_local` | `03_scrna_multidonor` | **`completed`** | `334.5s` | `15` | `23,031` | `$0.0000` | `$0.0000` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `baseline` | `tantalus_local` | `04_scatac_multidonor` | **`completed`** | `419.7s` | `8` | `3,568` | `$0.0000` | `$0.0000` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `baseline` | `tantalus_local` | `08_cut_and_tag` | **`completed`** | `910.0s` | `9` | `14,514` | `$0.0000` | `$0.0000` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `baseline` | `tantalus_local` | `09_hic_3d_genome` | **`completed`** | `2003.7s` | `6` | `5,716` | `$0.0000` | `$0.0000` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `baseline` | `tantalus_local` | `10_paired_multiome` | **`completed`** | `1041.5s` | `15` | `78,028` | `$0.0000` | `$0.0000` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
| `baseline` | `tantalus_local` | `11_cross_modal_concordance` | **`completed`** | `3082.4s` | `20` | `16,580` | `$0.0000` | `$0.0000` | `1.00` | `1.00` | `verified_portable_ext_inputs` |
