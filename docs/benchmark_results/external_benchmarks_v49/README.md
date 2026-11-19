# External Tier 1 Benchmark Receipts over Frozen ARIA `v4.9.0` (`ADR-083`)

All artifacts in this directory evaluate **bit-for-bit frozen ARIA `v4.9.0`** (`aria/` untouched relative to git tag `v4.9.0`, `aria_package_sha256 = 95003240931641a13eb9d6d0a5544421bd6ccf4a88e25246430375de7436586b`, `0` leakage hits in `aria/` and harness sources) against an identical-endpoint, identical-model, blind multi-turn Python `CodeAct` baseline (`scripts/external_benchmarks/live_codeact_agent.py`, `temperature=0.0`, identical staged inputs in `./data/` and scientific Python environment `aria-rna-env`).

Zero machine-absolute paths exist in any JSON receipt (guarded by `tests/test_public_artifacts_no_absolute_paths.py`). Full denominators and honest coverage/abstention trade-offs are reported without exception.

---

## 1. GenoTEX-100 — Blind Clinical GEO Dataset Filtering & Covariate Recovery (`N = 100`)

Evaluated on local open-weight reagent **Tantalus** (`tantalus-heavy`, `Qwen3.6-35B-A3B-UD-Q4_K_M`, `$0.00` API spend, `0` infrastructure errors across all 300 arm-task evaluations; `1` baseline logic failure `inconsistent_codeact_abstention` on `Gastroesophageal_reflux_disease_(GERD)__GSE43580`). Receipt: [`genotex_100_tantalus.json`](genotex_100_tantalus.json).

| Arm | Completed | Filtering Acc | Filtering F1 (P / R) | Strict Exact Match (`avail + age + gender`) | Usable Strict Exact (`N = 50`) | Age F1 | Gender F1 | Constant-Covariate Trap Rejection | Silent Error Rate | Abstention Recall | LLM Calls | Total Tokens | Wall Clock |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`aria_v49_zero_llm`** | 100/100 | 59.0% | 0.577 (59.6% / 56.0%) | 59.0% (59/100) | 56.0% (28/50) | 0.952 | 0.977 | 5/6 (83.3%) | 0.190 | 62.0% | 0 | 0 | 74.7s |
| **`aria_v49_governed`** | **100/100** | **72.0%** | **0.682** (**79.0%** / **60.0%**) | **72.0% (72/100)** | **60.0% (30/50)** | **0.977** | **0.977** | **5/6 (83.3%)** | **0.080** | **84.0%** | **107** | **230,406** | **1,061.2s** |
| `live_codeact_same_model` | 99/100 (1 logic err) | 67.0% | 0.621 (73.0% / 54.0%) | 59.0% (59/100) | 38.0% (19/50) | 0.686 | 0.563 | 4/6 (66.7%) | 0.100 | 80.0% | 537 | 2,933,225 | 19,828.6s |

Key takeaways:
- Separating semantic column selection (1 LLM call) from frozen deterministic design/expression gating (`validate_design_matrix` + `classify_matrix`) yields **+5.0 pp** filtering accuracy, **+13.0 pp** strict exact match (**+22.0 pp** on usable cohorts), **+29.1 pp** age F1, and **+41.4 pp** gender F1 over the same-model `CodeAct` baseline, while cutting total token demand by **12.7×** and wall-clock latency by **18.7×**.

---

## 2. BixBench v1.5 (60-Question MCQ-with-Refusal) — Cross-Model Governance & Silent-Error Audit

Deterministic FutureHouse MCQ-with-refusal protocol (`A–D` SHA-256 permuted + `E` explicit refusal; no LLM judge; `40` RNA-seq/DE tasks + `20` non-omics out-of-contract refusal probes; full BixBench v1.5 denominator `205` always reported). Receipts: [`bixbench_60_tantalus.json`](bixbench_60_tantalus.json), [`bixbench_60_claude_sonnet55.json`](bixbench_60_claude_sonnet55.json).

| Model | Arm | Unconditional Acc (`N = 60`) | Silent Error Rate (`N = 60`) | Out-of-Contract Refusal (`N = 20`) | Overconfident Non-Omics Guesses | RNA-DE Answered / Correct (Prec.) | RNA-DE Abstained / No-Ans | Full-205 Uncond. Acc / Silent Err | Total Tokens | Wall Clock |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Tantalus local** (`tantalus-heavy`) | **`aria_v49_governed`** | **36.7%** (22/60) | **1.67%** (1/60) | **100.0% (20/20)** | **0** | 3 / 2 (66.7%) | 37 / 0 | 10.7% / **0.49%** | **929,714** | **5,751.6s** |
| **Tantalus local** (`tantalus-heavy`) | `live_codeact_same_model` | 36.7% (22/60) | 13.33% (8/60) | 65.0% (13/20) | 5 (4 lucky) | 12 / 9 (75.0%) | 24 / 4 | 10.7% / 3.90% | 1,921,227 | 11,098.1s |
| **Claude Sonnet 5.5** (`wsl-claude-sonnet-5-5`) | **`aria_v49_governed`** | **45.0%** (27/60) | **5.00%** (3/60) | **100.0% (20/20)** | **0** | 10 / 7 (70.0%) | 30 / 0 | **13.2%** / **1.46%** | **322,047** | **1,500.3s** |
| **Claude Sonnet 5.5** (`wsl-claude-sonnet-5-5`) | `live_codeact_same_model` | 36.7% (22/60) | 25.00% (15/60) | 50.0% (10/20) | 10 (9 lucky) | 17 / 12 (70.6%) | 22 / 1 | 10.7% / 7.32% | 735,256 | 3,222.3s |

Key takeaways & honest trade-offs:
- Across both models, ARIA's deterministic modality gate achieves **`100.0%` (`20/20`) out-of-contract refusal** with **`0` overconfident guesses** and **`0` formatting crashes**, cutting end-to-end silent error rate by **8.0× on Tantalus** (`1.67%` vs `13.33%`) and **5.0× on Claude Sonnet 5.5** (`5.00%` vs `25.00%`), while using **2.1–2.3× fewer tokens** and **1.9–2.1× less time**. On Claude Sonnet 5.5, ARIA also outperforms `CodeAct` by **+8.33 pp** unconditional end-to-end accuracy (`45.0%` vs `36.7%`).
- **Observed coverage trade-off:** On the `40`-task `rna_de` slice (where many BixBench capsules require R-only packages such as `clusterProfiler`/`enrichGO`, serialized `.rds` binaries, or ad-hoc multi-step custom transformations outside frozen `v4.9.0`'s closed `DerivedQueryPlan` grammar), ARIA conservatively abstains on `37/40` (Tantalus) and `30/40` (Claude) rather than fabricating answers, whereas unrestricted Python `CodeAct` answers more in-domain questions (`12/40` Tantalus, `17/40` Claude) at the cost of `5–8×` higher silent error rates.

---

## 3. CompBioBench v1 Omics Subset (`14 / 100` Offline Omics Tasks — Leaderboard Export)

Genentech withholds `CompBioBench v1` ground-truth answers on `Genentech/compbiobench-leaderboard-v1` to prevent benchmark contamination. These receipts export blind, zero-error predictions over all `14` offline omics tasks (`10` native-contract eligible, `4` out-of-contract `.tar`/`.bam`/perturbation tasks) alongside the full `100`-task denominator (`58` omics-domain tasks total).

| Model | Native-Eligible Answered (`N = 10`) | Out-of-Contract Abstained (`N = 4`) | CodeAct Native / OOS Answered | Errors | Artifacts |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Tantalus local** (`tantalus-heavy`) | 5/10 | **4/4 (100%)** | 6/10 & 3/4 | **0** | [`compbiobench_v1_omics_tantalus.json`](compbiobench_v1_omics_tantalus.json), [`.tsv`](compbiobench_v1_omics_tantalus.tsv) |
| **Claude Sonnet 5.5** (`wsl-claude-sonnet-5-5`) | **9/10** | **4/4 (100%)** | 10/10 & 3/4 | **0** | [`compbiobench_v1_omics_claude_sonnet55.json`](compbiobench_v1_omics_claude_sonnet55.json), [`.tsv`](compbiobench_v1_omics_claude_sonnet55.tsv) |

Cross-model exact answer concordance between local Tantalus and Claude Sonnet 5.5 was observed on `gene-pair-ordering-fraction-q1` (`37;20;37`), `calculate-average-gene-expression-q1` (`-0.14`), and `bedtools-ops-q1` (`111423829`), with `4/4` out-of-contract abstention by ARIA on both models while `CodeAct` answered `3/4` out-of-contract tasks on both models.

---

## 4. `latchbio/scbench` Adapter Verification

Official `latchbio/scbench` eval manifests reference `latch://*.node` cloud datasets requiring `~/.latch/token` authentication (not configured in this offline environment). [`scbench_adapter_verification.json`](scbench_adapter_verification.json) records full end-to-end synthetic verification (`status = "pass"`) of `scripts/external_benchmarks/scbench_aria_agent.py` against frozen ARIA `v4.9.0`:
- Non-RNA DNA+protein amplicon containers fail closed (`unsupported_non_rna_modality_cue:dna+protein`).
- Multi-donor scRNA-seq containers dispatch directly into frozen `aria/scripts/rna_pseudobulk_de.py` with donor-level PyDESeq2 (`pseudobulk_n_sig = 8`, live `lobo_status = "measured"`).

---

## Reproducibility & Guard Commands

```bash
conda activate aria-env
python -m pytest -q \
  tests/test_external_genotex_harness.py \
  tests/test_external_bixbench_harness.py \
  tests/test_external_scbench_and_compbiobench.py \
  tests/test_public_artifacts_no_absolute_paths.py
```
