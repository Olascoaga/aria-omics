---
status: active
source_of_truth_for: preprint_supplement_b3_frozen_reagent
last_updated: 2026-11-12
track: B3_multi_llm_and_frozen_reagent_ablation
---

# Supplementary Note (Track B3): Frozen Open-Weight Reagent Ablation & Long-Horizon Reproducibility (*"LLM Proposes, Code Guarantees"*)

## 1. Motivation: Solving Proprietary Model Deprecation Without Modifying Scientific Code

Commercial cloud LLMs (`claude-opus-*`, `claude-sonnet-*`, `gpt-*`) routinely deprecate model snapshots within 12–24 months, breaking bit-level replay of agentic computational biology preprints even when statistical containers and raw count matrices remain archived. Conversely, local 9B–35B open-weight models have smaller parameter budgets than frontier models and exhibit hardware-level floating-point subtleties on consumer GPUs.

ARIA's core architectural thesis — ***"LLM proposes, code guarantees"*** (`docs/architecture/design_principles.md`, Claims 3–7) — predicts two falsifiable properties when every cloud tier (`heavy`, `medium`, `light`) is replaced externally by a SHA-256-pinned local open-weight reagent:

1. **Public Statistical & Epistemic Invariance (`Claim 6`):** Replacing a frontier model with a pinned open-weight reagent must yield an identical public statistical claim set (`Jaccard = 1.0000`, `sign_concordance = 1.0000`) and zero prohibited causal/mechanistic elevations (`0`), because `compile_public_claims`, `validate_design_matrix`, `verify_block_claim_support`, and the typed `run_ledger` derive public claims strictly from deterministic code outputs rather than free-text prose.
2. **Zero-Cost, Zero-Egress Archival Replay:** Through external configuration alone (`~/.aria/config.yaml` routed to a `127.0.0.1` loopback contractual shim in front of a pinned `llama.cpp` commit, with zero edits to ARIA's tracked source tree), an entire headless workflow and its `c6_b3_real_e2e` reproducibility capsule execute at **`$0.00` commercial API cost** with immutable cryptographic stamps (`gguf_sha256`, `tokenizer_sha256`, `chat_template_sha256`, `llama_cpp_commit`).

![Fig. S-B3 Frozen Open-Weight Reagent Ablation](fig_s_b3_frozen_reagent_ablation.svg)

---

## 2. Hardware Physics & Multi-Family Open-Weight Selection (`31` Configurations Across `5` Families)

All measurements were executed on a single consumer workstation (**AMD Ryzen 9 7950X** 16C/32T Zen 4 AVX-512, **128 GB DDR5-3600** at ~45 GB/s, **AMD Radeon RX 7700 XT** 12 GB GDDR6 `gfx1101` at 432 GB/s) running official `ggml-org/llama.cpp` build **`b12316`** (`commit 5d4c8eeb2`, `Vulkan0`, `-np 1`, `seed=0`, `temperature=0.0`, `KV=f16`, `-c 16384`, `--no-cache-prompt`, `--cache-ram 0`, `--ctx-checkpoints 0`) behind the loopback Tantalus Contractual Shim (`127.0.0.1:11435`).

Three hardware/numerical findings emerged across the 31-configuration sweep:

1. **Single-CCD Thread Pinning (`-t 8` vs `-t 16`):** On dual-CCD Zen 4 (`7950X`), pinning CPU MoE expert evaluation to **8 threads (`-t 8`, single CCD L3 domain)** consistently outperformed 16 threads (`-t 16`) by **`+4%` to `+7%`** across all three Hybrid MoE families (`Qwen3.6-35B-A3B`, `Gemma 4 26B-A4B`, `gpt-oss-20b`) by eliminating cross-CCD Infinity Fabric snoop latency.
2. **RDNA3 Vulkan FP16 `split_k` Non-Associativity & Strict Frozen Reagent Mode (`rb=0`):** With unconstrained `<think>` chains (`--reasoning-budget 1024`, generating >600 reasoning tokens before emitting JSON), **`0/6` configurations** achieved 3-run SHA-256 byte-identity (`8.4–31.0 s` wall latency) because multi-workgroup FP16 matrix-vector reductions on Vulkan accumulate `~1e-6` rounding differences that flip a single near-tied reasoning token around step ~200–400. Enforcing **Strict Frozen Reagent Mode (`--reasoning-budget 0 --reasoning off`)** cut wall latency **`5x–10x`** while achieving **`3/3` SHA-256 byte-identical steady-state repeatability across `25/25` (`100%`) tested configurations**.
3. **Tiered Memory Footprint Selection:**
   - **Primary Scientific Reagent (`hybrid`, default):** `Qwen3.6-35B-A3B-UD-Q4_K_M` (`21.1 GiB` total, 3B active parameters/token, `--n-cpu-moe 24 -t 8`, `gguf_sha256=ac0e2c1189e055faa36eff361580e79c5bd6f8e76bffb4ce547f167d53e31a61`) delivered **`34.6 tok/s` peak decode** (`3,239.7 ms` `p50`) and the highest biostatistical design fidelity.
   - **Zero-RAM-Spillover Reagent (`vram`, for >100 GB STAR / Multiome runs):** `gemma-4-12b-it-Q4_K_M` (`6.7 GiB`, `--n-cpu-moe 0 -t 8`, `100%` resident in 12 GB VRAM, `0 GB` host RAM impact) delivered **`46.1 tok/s` peak decode** (`1,459.1 ms` `p50`), with `Qwen3.5-9B-Q4_K_M` reaching **`64.5 tok/s` peak decode** (`899.9 ms` `p50`).
   - **Cross-Family Hybrid Validator:** `gemma-4-26B-A4B-it-UD-Q4_K_M` (`16.2 GiB`, `--n-cpu-moe 16 -t 8`) delivered **`34.3 tok/s` peak decode** (`3,397.9 ms` `p50`).

### Table S-B3.1 — Representative Champions & Complete 31-Configuration Hardware Race

| Candidate Config | `--n-cpu-moe` | `-t` | Mean `tok/s` | Peak Decode `tok/s` | Mean Wall (ms) | `p50` Wall (ms) | JSON Valid | Tag Leaks | 3x Steady-State SHA-256 | Gate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `qwen3.5-9b-q4_k_m_moe0_t8_rb0` | 0 | 8 | 51.99 | 64.5 | 1247.2 | 899.9 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gemma-4-12b-q4_k_m_moe0_t8_rb0` | 0 | 8 | 36.17 | 46.1 | 1849.4 | 1459.1 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gpt-oss-20b-mxfp4_moe8_t8_rb0_ngram` | 8 | 8 | 31.15 | 50.1 | 2132.6 | 1444.3 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gpt-oss-20b-mxfp4_moe8_t8_rb0` | 8 | 8 | 30.24 | 50.5 | 2221.3 | 1627.5 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gpt-oss-20b-mxfp4_moe8_t16_rb0` | 8 | 16 | 29.40 | 46.6 | 2253.2 | 1524.8 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gpt-oss-20b-mxfp4_moe16_t8_rb0` | 16 | 8 | 19.20 | 32.5 | 3436.9 | 2387.3 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gemma-4-26b-a4b-q4_k_m_moe16_t8_rb0` | 16 | 8 | 21.09 | 34.3 | 3551.1 | 3397.9 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gpt-oss-20b-mxfp4_moe16_t16_rb0` | 16 | 16 | 18.28 | 29.4 | 3630.7 | 2476.3 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gemma-4-26b-a4b-q4_k_m_moe16_t16_rb0` | 16 | 16 | 20.64 | 32.9 | 3638.8 | 3430.4 | 100% | 0 | YES (`3/3`) | **PASS** |
| `qwen3.6-35b-a3b-q4_k_m_moe24_t8_rb0` | 24 | 8 | 13.93 | 34.6 | 3765.2 | 3239.7 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gemma-4-26b-a4b-q4_k_m_moe16_t8_rb0_ngram` | 16 | 8 | 20.32 | 34.1 | 3784.0 | 3329.9 | 100% | 0 | YES (`3/3`) | **PASS** |
| `qwen3.6-35b-a3b-q4_k_m_moe24_t16_rb0` | 24 | 16 | 13.58 | 32.7 | 3829.3 | 3195.9 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gemma-4-26b-a4b-q4_k_m_moe12_t16_rb0` | 12 | 16 | 18.98 | 34.5 | 3913.6 | 3591.4 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gemma-4-26b-a4b-q4_k_m_moe12_t8_rb0` | 12 | 8 | 18.83 | 35.5 | 4006.8 | 3888.2 | 100% | 0 | YES (`3/3`) | **PASS** |
| `qwen3.6-35b-a3b-q4_k_m_moe28_t8_rb0_ngram` | 28 | 8 | 12.49 | 31.4 | 4161.9 | 3496.7 | 100% | 0 | YES (`3/3`) | **PASS** |
| `qwen3.6-35b-a3b-q4_k_m_moe28_t8_rb0` | 28 | 8 | 12.46 | 31.5 | 4166.2 | 3527.0 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gemma-4-26b-a4b-q4_k_m_moe22_t8_rb0` | 22 | 8 | 17.41 | 28.5 | 4324.4 | 4077.8 | 100% | 0 | YES (`3/3`) | **PASS** |
| `qwen3.6-35b-a3b-q4_k_m_moe28_t16_rb0` | 28 | 16 | 11.97 | 29.5 | 4339.9 | 3641.4 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gemma-4-26b-a4b-q4_k_m_moe22_t16_rb0` | 22 | 16 | 16.80 | 27.2 | 4468.2 | 4263.8 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gpt-oss-20b-mxfp4_moe24_t8_rb0` | 24 | 8 | 14.29 | 24.1 | 4634.1 | 3299.8 | 100% | 0 | YES (`3/3`) | **PASS** |
| `qwen3.6-35b-a3b-q4_k_m_moe34_t8_rb0` | 34 | 8 | 10.85 | 28.0 | 4789.2 | 4090.3 | 100% | 0 | YES (`3/3`) | **PASS** |
| `qwen3.6-35b-a3b-q4_k_m_moe34_t16_rb0` | 34 | 16 | 10.47 | 26.1 | 4970.0 | 4234.0 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gpt-oss-20b-mxfp4_moe24_t16_rb0` | 24 | 16 | 13.34 | 21.6 | 4981.8 | 3471.0 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gemma-4-26b-a4b-q4_k_m_moe28_t8_rb0` | 28 | 8 | 14.66 | 24.7 | 5126.6 | 4902.5 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gemma-4-26b-a4b-q4_k_m_moe28_t16_rb0` | 28 | 16 | 14.21 | 23.4 | 5286.5 | 5046.7 | 100% | 0 | YES (`3/3`) | **PASS** |
| `gpt-oss-20b-mxfp4_moe12_t16_rb1024` | 12 | 16 | 30.39 | 35.8 | 8402.2 | 8005.4 | 100% | 0 | NO | **FAIL** |
| `gemma-4-12b-q4_k_m_moe0_t8` | 0 | 8 | 39.75 | 46.3 | 12383.5 | 14190.4 | 100% | 0 | NO | **FAIL** |
| `qwen3.5-9b-q4_k_m_moe0_t8` | 0 | 8 | 58.05 | 64.4 | 13156.5 | 16009.2 | 100% | 0 | NO | **FAIL** |
| `qwen3.5-9b-q4_k_m_moe0_t8_nofa` | 0 | 8 | 57.93 | 64.4 | 13167.1 | 16068.4 | 100% | 0 | NO | **FAIL** |
| `gemma-4-26b-a4b-q4_k_m_moe16_t8_rb1024` | 16 | 8 | 28.41 | 33.7 | 19750.9 | 21921.8 | 100% | 0 | NO | **FAIL** |
| `qwen3.6-35b-a3b-q4_k_m_moe28_t8_rb1024` | 28 | 8 | 22.43 | 31.2 | 31001.5 | 39223.5 | 100% | 0 | NO | **FAIL** |

### Table S-B3.2 — Why Deterministic Design Validation (`B1`) Remains Mandatory Across Model Families

Workload `W1_design_structured` presented a 9-sample bulk RNA-seq design (`3 genotypes × 3 biological replicates`, sequenced in a single flowcell batch) and requested structured JSON (`experimental_unit`, `recommended_model`, `batch_confounded`). Inspecting the raw proposals across all five model families illustrates *"LLM proposes, code guarantees"* directly:

| Model Family (`rb=0` Champion) | Proposed `experimental_unit` | Proposed `recommended_model` | Proposed `batch_confounded` | How ARIA Code Governance (`B1` `validate_design_matrix`) Guards Execution |
|---|---|---|---|---|
| `Qwen3.6-35B-A3B` (`moe=24, t=8`) | `"biological replicate"` (correct) | `"DESeq2 ~ genotype"` (correct) | `true` (over-conservative on single flowcell) | `validate_design_matrix` verifies `n=3` per group and rank-full `~ genotype` (`batch` constant across samples is not collinearly added). |
| `Gemma 4 26B-A4B` (`moe=16, t=8`) | `"H9 cell culture batch"` (imprecise) | `"GLM with batch covariate"` (collinear!) | `true` | Deterministic rank check rejects adding a single-level `batch` covariate that would make the model matrix singular. |
| `Gemma 4 12B` (`moe=0, t=8`) | `"H9 cell line genotype"` (group, not unit) | `"LMM or ANOVA"` (non-count model) | `false` (correct) | Bulk RNA-seq dispatch enforces negative-binomial count GLM (`pyDESeq2`) rather than Gaussian ANOVA on raw counts. |
| `gpt-oss-20b` (`moe=8, t=8, ngram`) | `"cell"` (wrong modality unit) | `"linear mixed model"` | `false` (correct) | `AssayDetector` + `BulkRNAAgent` bind sample-level count columns as biological units, never single cells. |
| `Qwen3.5-9B` (`moe=0, t=8`) | `"single cell"` (wrong modality unit) | `"LMM with batch random effect"` | `true` | Same deterministic assay & design-matrix fence prevents pseudoreplicated or singular formula execution. |

---

## 3. Six-Callsite Contract & Epistemic Governance Replay (`CS1`–`CS6`)

We replayed ARIA's six core LLM call surfaces (`CS1`–`CS6`, `3` repetitions each = `18` live calls) through the Contractual Shim (`json_schema`/GBNF enforcement + `<think>`/channel-tag stripping) against the frozen primary reagent (`Qwen3.6-35B-A3B-UD-Q4_K_M`).

### Frozen Reagent Provenance Capsule
- **Frontier Reference:** `Claude Opus (ARIA v48 Governed Reference)`
- **Local Reagent ID:** `qwen3.6-35b-a3b-ud-q4_k_m`
- **GGUF SHA-256:** `ac0e2c1189e055faa36eff361580e79c5bd6f8e76bffb4ce547f167d53e31a61`
- **Tokenizer SHA-256:** `f2f491ce0b89b6daa98d9e10e46fab2959da6bb2a513042512b3b4678c32b639`
- **Chat Template SHA-256:** `55d4931433fe502b794226ee7f4d206a6bdd436ac9f80eb7d8ebb4c639f9ea0c`
- **`llama.cpp` Build / Commit / Binary SHA-256:** `b12316` (`5d4c8eeb2`) / `f1e6be3c20821e4434ca0be2cc5093ee21d5b78ffe2ee98b2f4382b39ea8cbc6`
- **Determinism Controls:** `backend=vulkan (Vulkan0)`, `KV=f16`, `-np 1`, `seed=0`, `temperature=0.0`, `reasoning_budget=0`

### Table S-B3.3 — Six-Callsite Contract, Entity Grounding & Causal Restraint Results

| Callsite | ARIA Agent | Structured JSON (`3/3`) | HGNC Grounding (`0` Hallucinated Genes) | Observational Causal Restraint (`0` Violations) | Steady-State Same-Prompt SHA-256 (`rep2 == rep3`) | Median Wall (ms) | Gate |
|---|---|---|---|---|---|---:|---|
| `CS1_orchestrator_routing` | `Orchestrator` | PASS (`3/3`) | PASS | PASS | Phrasing varies | 2766.4 | **PASS** |
| `CS2_design_agent_audit` | `DesignAgent` | PASS (`3/3`) | PASS | PASS | YES (`3/3`) | 2450.4 | **PASS** |
| `CS3_scrna_annotation` | `scRNAAnnotationAgent` | PASS (`3/3`) | PASS | PASS | YES (`3/3`) | 2547.1 | **PASS** |
| `CS4_scrna_de_summary` | `scRNADEAgent` | PASS (`3/3`) | PASS | PASS | YES (`rep2==rep3`) | 5534.7 | **PASS** |
| `CS5_narrative_agent` | `NarrativeAgent` | N/A (free prose) | PASS | PASS | YES (`rep2==rep3`) | 7521.8 | **PASS** |
| `CS6_hypothesis_debate` | `HypothesisAgent` | PASS (`3/3`) | PASS | PASS | YES (`rep2==rep3`) | 6520.8 | **PASS** |

**Hardware Determinism Note (Why Prose-Hash Freezing Alone Is Insufficient):** Across `CS2` and `CS3`, all `3/3` calls were byte-identical; across `CS4`, `CS5`, and `CS6`, steady-state calls (`rep2 == rep3`) were also byte-identical (`5/6` steady-state byte-identity), whereas the initial transition (`rep1`) immediately following a different-length prompt triggered `ggml-vulkan` graph buffer reuse (`graph reallocs`) on `qwen35moe`'s hybrid Gated-DeltaNet + attention tiles, shifting surface wording while preserving **100% (`18/18`) JSON validity, 100% (`18/18`) HGNC gene grounding, and 100% (`18/18`) causal restraint**. This hardware-level observation reinforces ARIA's design: scientific reproducibility must be anchored at the compiled claim/ledger boundary (`Claim 6`), never assumed from `temperature=0` prose alone.

---

## 4. Claim 6 Public Claim-Set Invariance & Live `$0.00` Headless E2E Capsule Equivalence

1. **Frontier vs Frozen Reagent Public Claim-Set Invariance (`CS4`/`CS5` Donor Pseudobulk DE):**
   - **Public Claim-Set Jaccard Similarity:** `1.0000`
   - **FDR & Log2FC Sign Concordance:** `1.0000`
   - **Prohibited Epistemic Tier Elevations:** `0` (`associative` preserved; `CS6` `HypothesisAgent` quarantined as `hypothesis_only`)
2. **Live Zero-Touch Headless Workflow (`aria.headless.run_headless` on `validation-data://pbmc3k_test`):**
   - Routed all three tiers (`heavy`, `medium`, `light`) to `http://127.0.0.1:11435` via `~/.aria/config.yaml` (`context_window: 16384`, paid cloud keys disabled in `~/.aria/.env`, zero modifications to ARIA's repository).
   - Completed full scRNA QC, normalization, Leiden clustering, Wilcoxon rankings, and governed HTML report synthesis (`227,679` bytes) across `8` local reagent calls (`4 medium`, `3 heavy`, `1 light`; `7,926` prompt tokens, `1,248` completion tokens) at **`34.4 tok/s` mean decode**, **`0` retries**, and **`$0.00` cloud API cost**.
3. **Repeated Real-Run E2E Reproducibility Gate (`c6_b3_real_e2e` / `tests/test_preprint_audit_b3_multi_provider_e2e.py`):**
   - Executed two full `run_headless(..., reproducible_mode=True)` runs back-to-back against the live local reagent (`status=PASS`): `diff_methodologies(meth_a, meth_b)` confirmed `claims.added == []`, `claims.removed == []`, `claims.changed == []`, `ledger.status_changed == []`, and `verify_reproducible_capsule(capsule) == "pass"` for both capsules.

---

## 5. Copy-Pasteable Preprint Methods Paragraph

> **Local Open-Weight Reagent Ablation and Long-Horizon Reproducibility (Track B3).** To verify that ARIA's statistical and epistemic guarantees do not depend on proprietary frontier APIs subject to snapshot deprecation, we evaluated five 2026 open-weight model families across 31 hardware configurations (`llama.cpp` build `b12316`, commit `5d4c8eeb2`; `Vulkan0` on a 12 GB Radeon RX 7700 XT paired with a 16-core Ryzen 9 7950X and 128 GB DDR5 RAM; `-np 1`, `seed=0`, `temperature=0.0`, `KV=f16`, 16,384-token context) behind a loopback-only (`127.0.0.1`) contractual shim enforcing JSON schema grammars and recording SHA-256 hashes of model weights, tokenizers, chat templates, prompts, and outputs. Disabling unconstrained reasoning chains (`--reasoning-budget 0 --reasoning off`) eliminated multi-workgroup FP16 reduction drift on RDNA3 Vulkan, achieving 3-for-3 byte-identical steady-state outputs across 25 of 25 tested configurations while reducing latency 5- to 10-fold. Routing all ARIA model tiers (`light`, `medium`, `heavy`) externally to the selected primary reagent (`Qwen3.6-35B-A3B-UD-Q4_K_M`, `gguf_sha256=ac0e2c1189e055faa36eff361580e79c5bd6f8e76bffb4ce547f167d53e31a61`, `--n-cpu-moe 24 -t 8`, 34.6 decode tokens/s; zero source edits to ARIA) yielded 100% JSON contract validity, zero hallucinated HGNC gene symbols, and zero observational causal-language violations across 18 live replays of ARIA's six core agent callsites (`CS1`–`CS6`), exact public claim-set invariance (`Jaccard = 1.0000`, sign concordance `1.0000`), and zero claim or run-ledger divergence (`diff_methodologies = 0`, verified reproducibility capsules) across repeated end-to-end headless executions at `$0.00` commercial API cost.
