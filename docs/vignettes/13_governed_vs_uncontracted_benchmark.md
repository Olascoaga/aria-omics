# Vignette `13_governed_vs_uncontracted_benchmark` — 40-Capsule Governed vs. Uncontracted CodeAct Benchmark

- **Schema:** `aria.codeact_comparative_benchmark.v1` · **ADR:** `ADR-083` (Tier 1, deterministic synthetic capsules)
- **Receipt:** [`docs/benchmark_results/v49_comparative_benchmark/comparative_benchmark_v49.json`](../benchmark_results/v49_comparative_benchmark/comparative_benchmark_v49.json) · `receipt_sha256` = `b685e65700dd856486f1f011e793663555a52b227e5475d9937bf9d612cb0e5b`
- **Decisions digest:** `bdfd123167a806bea5432ec6b178cb5ec86a33f32842e1eac556af2784c8db47`
- **Capsules:** `40` = `8` x `5` archetypes (`22` answerable with known ground truth, `18` must-abstain traps), seeded from `490083` via `numpy.random.default_rng(PCG64)`
- **Reproduce:** `conda run -n aria-rna-env python scripts/run_v49_codeact_comparative_benchmark.py`

## 1. Question

How often does an analysis arm return a confident answer that is wrong (`silent_error_rate`), when the same 40 multi-omics tasks are solved by (i) ARIA's governed, contract-checked code path and (ii) an uncontracted CodeAct-style script with the analysis patterns these governance layers were built to stop?

## 2. Methodology

Each capsule is generated from its own seed (`numpy.random.default_rng`). It holds synthetic omics inputs, a confirmed design, a natural-language question, and a ground truth set by construction: planted features or links for answerable capsules, a required abstention for traps. Every arm is scored with `aria.benchmarks.adapter.evaluate_benchmark_predictions` (full denominator, 6 metrics).

| Arm | What runs |
|---|---|
| `aria_governed` | `resolve_design_contract` + `validate_design_matrix` + `factors_confounded_with_condition`; technical replicates are collapsed to biological units, and `compute_lobo_jackknife` refuses duplicate unit IDs; valid bulk designs use OLS on condition + confirmed covariates with Smyth (2004) variance moderation and `bh_adjust`; donor pseudobulk uses DESeq2 median-of-ratios size factors and `compute_counts_lobo_jackknife` (Smyth-moderated t, `bh_adjust`); `resolve_causal_estimand_license` + `compile_public_claims` (C4); `assert_no_speculative_promotion` (C7); `detect_unsupported_inputs`; `enforce_same_cell_integration_gate` plus a barcode-bijection check; `execute_derived_query`. Claim prose is drafted by a deterministic **aggressive scripted proposer** (always maximal causal wording), so every epistemic outcome is decided by governance, not by the drafter. |
| `aria_locked_manifest` | Each capsule is committed with `build_decision_manifest`, serialized, reloaded with `load_decision_manifest` (`verify_decision_manifest`), and replayed with the analysis policy read **only** from the verified manifest. A tampered copy (`alpha` 0.05 -> 0.2) of every manifest must fail verification. |
| `aria_zero_llm_deterministic` | All 40 capsules are regenerated from their seeds inside `isolated_benchmark_environment` (`ARIA_LLM_CACHE=0`, `ARIA_REPRODUCIBLE=1`, `PYTHONHASHSEED=0`) and solved with `llm=None` (template prose). |
| `uncontracted_codeact_baseline` | Deterministic scripted emulation of unguarded analysis code: Welch t-tests that treat every library as a replicate and ignore batch; cell-level `scipy.stats.ranksums` with cells as replicates; fold-change cutoffs on n=1 vs 1; unvalidated `pd.merge`; silent positional cell pairing when barcodes do not intersect; a FASTQ glob that falls through to an empty "completed" result; correlation reported as causation; hypotheses reported as findings. |

**Scoring (the same for every arm).** A silent error is a *confidently wrong claim*. An answerable DE-style capsule is correct when there are zero sign errors, precision over the planted features is `>= 0.9`, and recall is `>= 0.5` (the floor stops an empty answer from passing). Missing true features is a power limit, not a wrong claim, so power is reported separately (`power_mean_recall_on_planted_de_capsules`). Null capsules are correct only with **zero** calls. Every other answerable capsule needs an exact match (feature lists, module, top link, direction). A must-abstain trap is correct only if the arm abstains; any answer there is confident-wrong.

## 3. Headline metrics (all 6 adapter metrics per arm)

| Arm | `native_contract_coverage` | `conditional_accuracy_on_supported` | `unconditional_end_to_end_accuracy` | `abstention_precision` | `abstention_recall` | `silent_error_rate` | confident-wrong | power (mean recall) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `aria_governed` | 0.550000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 0.000000 | `0/40` | 0.962500 |
| `aria_locked_manifest` | 0.550000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 0.000000 | `0/40` | 0.962500 |
| `aria_zero_llm_deterministic` | 0.550000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 0.000000 | `0/40` | 0.962500 |
| `uncontracted_codeact_baseline` | 0.550000 | 0.636364 | 0.350000 | n/a | 0.000000 | 0.650000 | `26/40` | 0.868750 |

`abstention_precision` is `n/a` for an arm that never abstains: its denominator, the number of abstentions, is zero.

**Silent-error gap:** baseline `0.650000` vs governed `0.000000` (absolute gap `0.650000`).

**Locked-manifest replay:** `40/40` identical decisions (`decision_concordance` = `1.000000`, `prediction_concordance` = `1.000000`). Manifests verified `40/40`; tampered manifests detected `40/40`; input fingerprints matched `40/40`.

**Zero-LLM determinism:** `40/40` identical decisions (`decision_concordance` = `1.000000`); capsules regenerated byte-identically `40/40`.

## 4. Per-archetype breakdown

| Archetype | Arm | correct answer | correct abstention | confident wrong | wrongful abstention |
|---|---|---:|---:|---:|---:|
| `design_and_confounding_governance` | `aria_governed` | 4 | 4 | 0 | 0 |
| `design_and_confounding_governance` | `aria_locked_manifest` | 4 | 4 | 0 | 0 |
| `design_and_confounding_governance` | `aria_zero_llm_deterministic` | 4 | 4 | 0 | 0 |
| `design_and_confounding_governance` | `uncontracted_codeact_baseline` | 4 | 0 | 4 | 0 |
| `donor_pseudobulk_vs_cell_pseudoreplication` | `aria_governed` | 8 | 0 | 0 | 0 |
| `donor_pseudobulk_vs_cell_pseudoreplication` | `aria_locked_manifest` | 8 | 0 | 0 | 0 |
| `donor_pseudobulk_vs_cell_pseudoreplication` | `aria_zero_llm_deterministic` | 8 | 0 | 0 | 0 |
| `donor_pseudobulk_vs_cell_pseudoreplication` | `uncontracted_codeact_baseline` | 0 | 0 | 8 | 0 |
| `epistemic_claim_and_causal_licensing` | `aria_governed` | 2 | 6 | 0 | 0 |
| `epistemic_claim_and_causal_licensing` | `aria_locked_manifest` | 2 | 6 | 0 | 0 |
| `epistemic_claim_and_causal_licensing` | `aria_zero_llm_deterministic` | 2 | 6 | 0 | 0 |
| `epistemic_claim_and_causal_licensing` | `uncontracted_codeact_baseline` | 2 | 0 | 6 | 0 |
| `sequencer_and_modality_boundary_abstention` | `aria_governed` | 2 | 6 | 0 | 0 |
| `sequencer_and_modality_boundary_abstention` | `aria_locked_manifest` | 2 | 6 | 0 | 0 |
| `sequencer_and_modality_boundary_abstention` | `aria_zero_llm_deterministic` | 2 | 6 | 0 | 0 |
| `sequencer_and_modality_boundary_abstention` | `uncontracted_codeact_baseline` | 2 | 0 | 6 | 0 |
| `scrubbed_id_derived_query_reasoning` | `aria_governed` | 6 | 2 | 0 | 0 |
| `scrubbed_id_derived_query_reasoning` | `aria_locked_manifest` | 6 | 2 | 0 | 0 |
| `scrubbed_id_derived_query_reasoning` | `aria_zero_llm_deterministic` | 6 | 2 | 0 | 0 |
| `scrubbed_id_derived_query_reasoning` | `uncontracted_codeact_baseline` | 6 | 0 | 2 | 0 |

## 5. Seed-robustness sweep (independent re-draws)

The headline numbers come from one frozen draw (`suite_base_seed` = `490083`). To keep one lucky seed from passing as a guarantee, all 40 capsules were re-drawn under `20` alternative base seeds (`1_000_003 + 7_919 * i, i = 0..n-1`) and re-solved by `aria_governed` and `uncontracted_codeact_baseline`.

| Statistic | mean | min | max |
|---|---:|---:|---:|
| `aria_governed_silent_error_rate` | 0.002500 | 0.000000 | 0.025000 |
| `baseline_silent_error_rate` | 0.650000 | 0.650000 | 0.650000 |
| `aria_governed_conditional_accuracy_on_supported` | 0.995454 | 0.954545 | 1.000000 |
| `aria_governed_power_mean_recall` | 0.984375 | 0.962500 | 1.000000 |
| `baseline_power_mean_recall` | 0.888750 | 0.856250 | 0.937500 |

`aria_governed` has zero silent errors on `18/20` seeds, with perfect abstention precision and recall on `20/20`. Every residual governed error across `800` capsule-runs is listed here: `planted_donor_consistent_4v4:confident_wrong` x2. In total: `4` false-positive feature calls, `0` sign errors, `0` wrongful abstentions and `0` answered must-abstain traps. All residual errors are false-positive calls on answerable DE capsules, which BH-FDR control at 5% permits by design; no governance gate failed.

## 6. Full 40-capsule ledger

| Capsule | Archetype | Variant | Ground truth | `aria_governed` | `uncontracted_codeact_baseline` |
|---|---|---|---|---|---|
| `A1-01` | `design_and_confounding_governance` | `perfect_batch_confounding` | abstain | correct_abstention (`design_blocking:condition_covariate_confounding`) | confident_wrong |
| `A1-02` | `design_and_confounding_governance` | `perfect_batch_confounding` | abstain | correct_abstention (`design_blocking:condition_covariate_confounding`) | confident_wrong |
| `A1-03` | `design_and_confounding_governance` | `technical_replicate_inflation` | abstain | correct_abstention (`design_blocking:insufficient_condition_replicates`) | confident_wrong |
| `A1-04` | `design_and_confounding_governance` | `technical_replicate_inflation` | abstain | correct_abstention (`design_blocking:insufficient_condition_replicates`) | confident_wrong |
| `A1-05` | `design_and_confounding_governance` | `batch_balanced_4v4` | answer | correct_answer | correct_answer |
| `A1-06` | `design_and_confounding_governance` | `batch_balanced_4v4` | answer | correct_answer | correct_answer |
| `A1-07` | `design_and_confounding_governance` | `multifactorial_batch_balanced_4v4` | answer | correct_answer | correct_answer |
| `A1-08` | `design_and_confounding_governance` | `multifactorial_batch_balanced_4v4` | answer | correct_answer | correct_answer |
| `A2-01` | `donor_pseudobulk_vs_cell_pseudoreplication` | `null_donor_variance_3v3` | answer | correct_answer | confident_wrong |
| `A2-02` | `donor_pseudobulk_vs_cell_pseudoreplication` | `null_donor_variance_3v3` | answer | correct_answer | confident_wrong |
| `A2-03` | `donor_pseudobulk_vs_cell_pseudoreplication` | `null_donor_variance_4v4` | answer | correct_answer | confident_wrong |
| `A2-04` | `donor_pseudobulk_vs_cell_pseudoreplication` | `null_donor_variance_4v4` | answer | correct_answer | confident_wrong |
| `A2-05` | `donor_pseudobulk_vs_cell_pseudoreplication` | `planted_donor_consistent_4v4` | answer | correct_answer | confident_wrong |
| `A2-06` | `donor_pseudobulk_vs_cell_pseudoreplication` | `planted_donor_consistent_4v4` | answer | correct_answer | confident_wrong |
| `A2-07` | `donor_pseudobulk_vs_cell_pseudoreplication` | `planted_donor_consistent_4v4` | answer | correct_answer | confident_wrong |
| `A2-08` | `donor_pseudobulk_vs_cell_pseudoreplication` | `planted_donor_consistent_4v4` | answer | correct_answer | confident_wrong |
| `A3-01` | `epistemic_claim_and_causal_licensing` | `observational_causal_overreach` | abstain | correct_abstention (`causal_claim_not_licensed:missing_design_contract`) | confident_wrong |
| `A3-02` | `epistemic_claim_and_causal_licensing` | `observational_causal_overreach` | abstain | correct_abstention (`causal_claim_not_licensed:missing_design_contract`) | confident_wrong |
| `A3-03` | `epistemic_claim_and_causal_licensing` | `observational_causal_overreach` | abstain | correct_abstention (`causal_claim_not_licensed:causal_prerequisites_not_verified`) | confident_wrong |
| `A3-04` | `epistemic_claim_and_causal_licensing` | `observational_causal_overreach` | abstain | correct_abstention (`causal_claim_not_licensed:design_confounding_unresolved`) | confident_wrong |
| `A3-05` | `epistemic_claim_and_causal_licensing` | `speculative_hypothesis_promotion` | abstain | correct_abstention (`speculative_promotion_blocked`) | confident_wrong |
| `A3-06` | `epistemic_claim_and_causal_licensing` | `speculative_hypothesis_promotion` | abstain | correct_abstention (`speculative_promotion_blocked`) | confident_wrong |
| `A3-07` | `epistemic_claim_and_causal_licensing` | `randomized_interventional_licensed` | answer | correct_answer | correct_answer |
| `A3-08` | `epistemic_claim_and_causal_licensing` | `randomized_interventional_licensed` | answer | correct_answer | correct_answer |
| `A4-01` | `sequencer_and_modality_boundary_abstention` | `raw_sequencer_bcl_intake` | abstain | correct_abstention (`unsupported_sequencer_native_input:BCL`) | confident_wrong |
| `A4-02` | `sequencer_and_modality_boundary_abstention` | `raw_sequencer_pod5_intake` | abstain | correct_abstention (`unsupported_sequencer_native_input:POD5`) | confident_wrong |
| `A4-03` | `sequencer_and_modality_boundary_abstention` | `unpaired_split_diagonal_imputation` | abstain | correct_abstention (`same_cell_gate_blocked:multiome_contract_pairing_unconfirmed`) | confident_wrong |
| `A4-04` | `sequencer_and_modality_boundary_abstention` | `unpaired_split_diagonal_imputation` | abstain | correct_abstention (`same_cell_gate_blocked:multiome_contract_pairing_unconfirmed`) | confident_wrong |
| `A4-05` | `sequencer_and_modality_boundary_abstention` | `n1_vs_n1_unreplicable_contrast` | abstain | correct_abstention (`design_blocking:insufficient_condition_replicates`) | confident_wrong |
| `A4-06` | `sequencer_and_modality_boundary_abstention` | `n1_vs_n1_unreplicable_contrast` | abstain | correct_abstention (`design_blocking:insufficient_condition_replicates`) | confident_wrong |
| `A4-07` | `sequencer_and_modality_boundary_abstention` | `paired_same_cell_barcode_bijection` | answer | correct_answer | correct_answer |
| `A4-08` | `sequencer_and_modality_boundary_abstention` | `paired_same_cell_barcode_bijection` | answer | correct_answer | correct_answer |
| `A5-01` | `scrubbed_id_derived_query_reasoning` | `top5_up_significant` | answer | correct_answer | correct_answer |
| `A5-02` | `scrubbed_id_derived_query_reasoning` | `top5_down_significant` | answer | correct_answer | correct_answer |
| `A5-03` | `scrubbed_id_derived_query_reasoning` | `module_with_most_significant` | answer | correct_answer | correct_answer |
| `A5-04` | `scrubbed_id_derived_query_reasoning` | `module_with_highest_mean_lfc_among_significant` | answer | correct_answer | correct_answer |
| `A5-05` | `scrubbed_id_derived_query_reasoning` | `top3_signaling_significant_by_padj` | answer | correct_answer | correct_answer |
| `A5-06` | `scrubbed_id_derived_query_reasoning` | `significant_with_lfc_ge_1_sorted` | answer | correct_answer | correct_answer |
| `A5-07` | `scrubbed_id_derived_query_reasoning` | `invalid_cross_unit_gene_x_cell_join` | abstain | correct_abstention (`derived_query_violation`) | confident_wrong |
| `A5-08` | `scrubbed_id_derived_query_reasoning` | `invalid_one_to_one_cardinality_join` | abstain | correct_abstention (`derived_query_violation`) | confident_wrong |

## 7. Honest scope & limitations

- **Synthetic Tier 1 only.** Capsules are seeded simulations with ground truth set by construction. They test governance boundaries (design, units, licensing, intake, same-cell identity, join semantics) and show nothing about accuracy on real cohorts. Public-data validation lives in vignettes `01`-`12`.
- **The baseline is a scripted emulation, not a live LLM agent.** `uncontracted_codeact_baseline` runs real code (`scipy.stats.ranksums`, Welch t-tests, `pd.merge`, Pearson correlation), but which unguarded pattern it applies to each capsule was fixed by the authors. Two behaviors are modeled outright: the FASTQ glob that silently falls through on BCL/POD5 intake, and the positional pairing fallback for unpaired cells. Its `silent_error_rate` is the error rate *of those patterns on these capsules*. It is not an estimate for any specific CodeAct product or model; a careful human or agent could avoid every one of these traps.
- **No live LLM in any arm.** The `aria_governed` claim drafter is a deterministic aggressive stand-in. The zero-LLM arm uses template prose. Decisions are identical across both by construction of the governance path, which is the claim being tested; this says nothing about end-to-end LLM planning quality.
- **The capsule mix sets the rates.** `18/40` capsules are must-abstain traps. Both the baseline's silent-error rate and ARIA's coverage depend on this designed ratio, not on real-world prevalence.
- **A frozen-seed 0.0 is not a guarantee.** BH-FDR control permits occasional false-positive calls, and the seed sweep above shows them. The frozen receipt is a regression anchor, not a bound.
- **Archetype 5 tests unit and cardinality guards.** It does not test free-form querying: ARIA executes a fixed typed AST per question, and the baseline runs equivalent pandas code. Where both are valid they agree.
- **The narrative claim compiler runs on the `aria-env` stack.** When the runner executes in `aria-rna-env`, the C4/C7 gate runs in a stack-isolated subprocess (JSON IPC); the receipt records the execution mode in its `runtime` block, which is excluded from `receipt_sha256`.
