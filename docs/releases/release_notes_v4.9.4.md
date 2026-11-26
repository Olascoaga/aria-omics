# ARIA v4.9.4 Release Notes

`v4.9.4` is a **pre-submission governance, benchmark-rigor, public-freeze hygiene,
and installation-UX hardening release** over `v4.9.3`, closing the engineering
action items identified during adversarial preprint review while keeping all
`9 / 9` production modalities and `3 / 3` cross-cutting multi-modal engines
verified and `100%` CI green.

## What Changed

- **Fail-Closed `DecisionManifest` Commitment & `realized_decisions_sha256` Integrity:**
  - `OrchestratorAgent._dispatch_agents` now halts fail-closed (`reason="decision_manifest_commitment_failed"`) on any `DecisionManifest` verification or commitment failure instead of logging a warning and continuing dispatch.
  - `aria/utils/decision_manifest.py` adds `seal_realized_parameter_decisions()`, cryptographically binding post-lock `realized_parameter_decisions` to the pre-execution `manifest_sha256` via `realized_decisions_sha256` while preserving `manifest_sha256` invariance; `verify_decision_manifest()` and replay dispatch reject unsealed or tampered realized decisions.
- **Strictly Label-Blind `B2` 4-Arm Component Ablation:**
  - `aria/benchmarks/governance_b2.py` now strips every gold annotation metadata key before any arm inspects a candidate `NarrativeBlock`, running a strictly label-blind component ablation across `naive_llm` (emit all), `guards_off` (`verify_block_claim_support` verifier only), `governed` (full `annotate_claim_tiers` + `validate_block` + `verify_block_claim_support` stack), and `template_only` (deterministic quantitative readout policy), reporting exact numerators/denominators and 95% Wilson CIs.
- **Paired Cross-Arm `CP2.5` Detection, Preprocessing Tier Classification & Covariate-Aware Cross-Modal CPM Guard:**
  - `DesignAgent` + `BulkRNAAgent` + `aria/headless.py::default_answer_policy` now detect cross-arm matched biological units at `CP2.5` and model matched donor blocks (`~ donor + condition`), fixing unattended `GSE52778` paired airway runs while failing closed on ambiguous within-group technical collisions.
  - `claim_compiler.py` and `evidence_verifier.py` classify chromatin preprocessing and dimensionality-reduction summaries (`scatac_preprocessing`, `bulk_atac_preprocessing`, `dimensionality_reduction`, `peak_calling`, `peak_count_matrix`, `peak_annotation`, `regulatory_layers`) strictly as `descriptive`.
  - `aria/utils/cross_modal_concordance.py` Layer D adds OLS `sample_covariates` residualization alongside mandatory `cpm_sample_matched_layer_unadjusted_for_covariates` and `unadjusted_sex_chromosome_links_detected` warnings when covariates are absent (empirically verified on 10-donor `GSE278576` human hippocampus Arm B).
- **100% Public & Path-Portable Preprint Freeze (`19 / 19`), Real-Data Pseudoreplication Null & Scope-Gated `CodeAct` Control:**
  - Replaced the restricted chr21 H9 smoke lane (`c1_h9_fastq_e2e`) with `c1_public_fastq_e2e` running public paired `GSE52778` raw FASTQs (`fastp`/`STAR`/`featureCounts`/paired `pyDESeq2` `~ donor + condition`).
  - Removed all machine-local paths from `aria/benchmarks/preprint_freeze.py` and `scripts/run_scatac_multisample_da_concordance.py`, and added `ARIA_PROVENANCE_EXCLUDE_PATH` in `aria/version.py` so sequential freeze writes preserve clean `git_dirty=false` provenance.
  - Added `scripts/run_kang_real_pseudoreplication_null.py` (all 35 unique balanced 4v4 donor splits $\times$ 6 eligible PBMC cell types = 210 real-data null tests on unstimulated `Kang GSE96583`: `112` pseudobulk FPs vs `4,198`/`4,248` naive cell-level Welch/Wilcoxon FPs, a `37.5×–37.9×` reduction) and `scripts/external_benchmarks/run_bixbench_scope_gated_control.py` (zero-API deterministic `CodeAct + ARIA scope gate` ablation decoupling OOS refusal from closed-grammar grounding over frozen `v4.9.0` `BixBench-60` receipts).
- **One-Command Lock-First Environment Provisioning (`aria setup`) & `count_classifier` Path Fix:**
  - Added `aria setup [profile]` (`aria/setup_cli.py`, wired in `aria/tui.py` and `install.sh`), defaulting to `aria setup rna` (1 environment covering Bulk RNA matrices & scRNA-seq) and supporting `raw-rna`, `chromatin`, `tobias`, `multiome`, `hic`, `all`, `--list`, and `--dry-run` via committed `@EXPLICIT` `.linux-64.lock` + `.pip.lock` files.
  - Updated `SetupAgent._ensure_env()` to prefer exact lockfiles over `.yml` solving and stopped requesting `aria-tobias-env` for `ChIP` / `CUT_AND_RUN` / `CUT_AND_TAG`.
  - Fixed `aria/utils/count_classifier.py::_tool_signature_score` to match normalization hints on the filename basename only so random temporary parent directories never falsely reject raw integer count matrices.

## Compatibility & Invariants

- **Zero diff against `v4.9.0`** on `aria/utils/stage_manifest.py` and `docs/benchmark_results/preprint_v1/human/*.csv`.
- External Tier 1 blind benchmarks (`docs/benchmark_results/external_benchmarks_v49/genotex_100_*.json`, `bixbench_60_{tantalus,claude_sonnet55}.json`, `compbiobench_v1_*.json`) remain untouched and explicitly attributed to immutable tag `v4.9.0`.
- All `19 / 19` `automated_v1_gate` machine freeze lanes are 100% public/synthetic and freshly verified on clean `v4.9.4` (`git_dirty=false`, `0` machine blockers); the `3 / 22` independent `human_gold` reader-study lanes remain explicitly preserved as `deferred_to_v2_reader_study` (`ADR-083`).
