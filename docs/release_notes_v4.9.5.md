# v4.9.5 — systematic tech-debt audit (D1–D10): silent failures, data alignment, contracts, convergence, test oracles, CI coverage and a debt ratchet

`v4.9.5` is the batched remediation release closing ARIA's **systematic tech-debt
audit (slices D1–D10)** on top of `v4.9.4` (`8a0cb68..4e70292`, 123 non-merge
commits plus the `9b6e41e` D2/D3 merge). All `9 / 9` single-assay modalities and
`3 / 3` cross-cutting multi-modal engines remain at their existing readiness
tiers, with **`0` diff against `v4.9.0`** on `aria/utils/stage_manifest.py` and
`docs/benchmark_results/preprint_v1/human/*.csv`.

---

## Scope and Method

- **Auditor / verifier / fixer separation:** Every audit slice ran read-only to
  emit structured findings; an independent verifier re-executed each finding's
  `repro_test` on clean HEAD in its target Conda environment before triage;
  approved fixes then landed separately on `fix/d1-silent-failures` and
  `test/d2-d3-oracles` (merged cleanly at `9b6e41e`).
- **Executable evidence bar:** No finding reached `CONFIRMED` without a
  self-contained failing `repro_test` on HEAD. Unconfirmed findings remained
  `PLAUSIBLE` (or `REJECTED`) and were left untouched unless guarded alongside
  confirmed defects.
- **FAIL→PASS regression proof:** Every fix round verified its permanent
  regression tests against clean `git archive` snapshots of its starting base vs.
  HEAD across the relevant Conda environments (`aria-env`, `aria-rna-env`,
  `aria-chromatin-env`, `aria-hic-env`, `aria-integration-env`):
  - **D1:** `52 / 52` FAIL→PASS (`19` Round 1 + `14` Round 2 + `19` Round 3).
  - **D4:** `25 / 25` FAIL→PASS across all 5 Conda environments.
  - **D7:** `26 / 26` FAIL→PASS across 3 Conda environments.
  - **D8:** `27 / 27` FAIL→PASS across all 5 Conda environments.
  - **D5 + D9:** `36 / 36` FAIL→PASS in `aria-env` (`27` D5 + `9` D9) plus `1 / 1`
    in `aria-rna-env` (`D9-001`).
  - **D2:** all `10 / 10` real bugs exposed as strict xfails (`D2-BUG-1..10`)
    flipped `XFAIL → PASS` (`0` xfails remaining in `tests/oracle/`).

### Slice Totals (from Audit Charter §8)

| Slice | Scope | Findings (`C / H / M / L`) | Disposition | Metric Before → After | Commits |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **D1** | Silent failures in scientific kernels | `7 / 12 / 5 / 1` (`20 CONFIRMED`, `5 PLAUSIBLE`) | `20 / 20 CONFIRMED` fixed (Track A R1–R2 + Track B `W-WARN` R3); `5 PLAUSIBLE` untouched | Silent broad handlers `360 → 346` (`11` annotated auxiliary broad catches) | `aa2de54`..`3b22849` |
| **D2** | Statistical oracles + mutation testing | `10` real bugs exposed (`D2-BUG-1..10`: `5` new, `5` corroborated with D4/D8) | `10 / 10` fixed (`0` strict xfails remain in `tests/oracle/`); `182` oracle items in `10` modules | Core mutation score `70.2% → 82.75%` (`11/12` modules `≥ 80%`; `peak2gene_core` `79.2%` raw / `>87%` non-equiv) | `0c8ca78`..`5c96c60`, `9b6e41e`, `0098611`..`297ee90` |
| **D3** | Test quality & CI coverage | `4 / 13 / 5 / 0` (`22`) | `22 / 22` addressed (`27 / 27` manual mutations killed: `22` D3 + `5` heavy-module) | CI native file collection `46.8% → 100%` (`126/269 → 281/281` at D3 close; dynamic `tests/` on HEAD) | `3ea58d5`, `9946bdd`, `48b5bad`, `5c96c60`, `af1aad5`, `bed54c0`, `1884076`, `4e70292` |
| **D4** | Silent data misalignment | `2 / 8 / 7 / 0` `CONFIRMED` (`17`) + `2 PLAUSIBLE` | `19 / 19` fixed + `67` kernel `B905` + `13` in-scope `F841` | Kernel `B905` (`zip()` without `strict=`) `67 → 0`; `F841` `13 → 0` | `73e9f7b`..`d284dcb` |
| **D5** | Silent failures in claims & orchestration | `5 / 17 / 1 / 0` `CONFIRMED` (`23`) + `2 PLAUSIBLE` | All `25` findings + all remaining classification sites handled (`0` unannotated in D5 scope) | Silent broad handlers `334 → 294` (`# aria: allow-broad-except` `11 → 162`) | `5b7af89`..`d57ad7b` |
| **D6** | Monolithic function decomposition plan | `19` targets (`17` `≥ 400L`, `15` `CC ≥ 100`) | Plan delivered & reviewed (`118` pure + `29` I/O extractions, `60` pre-refactor oracles); **execution deferred** | Planned: `≥ 400L` `17 → 0`, `CC ≥ 100` `15 → 0` | Plan only (`D6_2026-11-26_plan.{md,json}`) |
| **D7** | Contracts between layers | `2 / 15 / 4 / 0` `CONFIRMED` (`21`) | `21 / 21` fixed across `665` audited contract rows | `22` `orphan_write`, `13` `masked_read`, `2` `type_mismatch`, `4` positional resolved | `4856aa3`..`cf88bc8`, `876f3ea` |
| **D8** | Diverged duplication across modalities | `1 / 11 / 10 / 1` `CONFIRMED` (`23`) + `2 DUPLICATE` of D1 | `23 / 23` fixed (`2` duplicates skipped) + absorbed cluster hygiene (`CL-028/030/031/034/040`) | `48` clusters: `29` divergent → converged onto `aria/utils/`; **`-131` net LOC** in `aria/` | `b79b9bb`..`ea49754`, `50d99cd` |
| **D9** | Hygiene, dead code & determinism | `0 / 3 / 6 / 0` `CONFIRMED` (`9`) + `1 PLAUSIBLE` | All `10` fixed + safe deletions + `PYTHONHASHSEED`-invariant verifier/B2 sorting | `F401` `101 → 18` live re-exports, `F811` `16 → 0`, `F821` `3 → 0`, `vulture` `6 → 2` | `ad4ae7d`..`8a70a4d` |
| **D10** | Structural debt ratchet | — | Stdlib-AST ratchet guard (`tests/test_debt_ratchet.py`, `tests/debt_ratchet.json`, slack `0`, ADR to raise) | Unannotated broad `507 → 313`; silent `148 → 99`; kernel `zip` w/o `strict` `67 → 0`; fns `> 300L` `33 → 35` | `642f3e0` |

---

## Behavior Changes Users Can Notice

- **Fail-closed I/O and cohort integrity in scientific kernels (`D1`, `D4`, `D5`):**
  - Corrupt or unreadable inputs now halt the stage fail-closed instead of silently
    returning empty/partial results or dropping samples: corrupt `.hic` blocks in
    `hic_preprocess.py` (`D1-001`), truncated GTFs in `rna_bulk/gtf_io.py` and
    `chromatin_peak_annotation.py` (`D1-003`, `D1-021`), corrupt `fragments.tsv.gz`
    and unreadable peak files in `chromatin_qc.py` / `chromatin_peaks.py` (`D1-008`,
    `D1-015`, `D1-022`), corrupt per-sample bulk count files (`D1-006`), corrupt
    BAM/CRAM in bulk chromatin QC (`D1-018`), corrupt Hi-C samples (`D1-007`),
    corrupt 3rd modalities in `integration_mofa.py` (`D1-013`), corrupt background
    `.h5ad` in scRNA (`D5-014`), corrupt SQLite reads in `NarrativeAgent` (`D5-019`),
    and unreadable bulk/chromatin metadata CSVs (with `utf-8` → `utf-8-sig` →
    `cp1252` → `latin-1` encoding detection in `ChromatinAgent`, `D5-008`, `D5-015`).
  - Failed bulk sample QC (`D1-009`, `642867f`, `1bfdff8`) now reports `status="error"`
    and `outliers=None` uniformly across all 6 consumers via
    `aria/utils/bulk_qc_contract.py` (`sample_qc_failed`, `sample_qc_size_ratio`)
    while preserving real library `size_ratio` when computable rather than emitting
    `"0 outlier(s)"` / `"1.0x"` at `HIGH` confidence.
  - Auxiliary diagnostic paths (LOBO ×4, PCA/MDS rendering, bulk sample QC, VST
    library-size fallback `D1-023`, apeGLM shrinkage, and Enrichr network calls
    `D1-025`) intentionally retain annotated broad catches (`# aria: allow-broad-except`)
    so diagnostic failures record errors/warnings without aborting primary DE/DA
    (`b231fef`). In `hic_topology.py`, plain-gzip FASTAs and read-only FASTA
    directories without `.fai` fall back to GTF TSS-density phasing with an explicit
    `phasing_fallback_reason` warning, while post-open truncated-sequence `fetch()`
    corruption fails closed (`79a0d37`).
  - Additional fail-closed alignment guards: duplicate sample IDs in bulk metadata
    and counts (`D4-004`), missing `obs` metadata / non-invariant within-donor
    pseudobulk covariates (`D4-006`), mismatched `scale_factors` in SEACR/broad
    domains (`D4-007`), zero control contig overlap in broad domains (`D4-008`),
    missing Hi-C files or `sample_ids` length mismatch (`D4-009`), duplicate
    barcodes in scATAC label transfer and `peak2gene` (`D4-010`), duplicate
    `replicate_ids` within condition in TOBIAS (`D4-015`), verified/reindexed
    `adata.raw.X` `obs`/`var` alignment (`D4-002`), index-reindexed `pd.Series`
    `size_factors` in LOBO (`D4-003`), ambiguous multi-locus RNA symbol exclusion
    instead of splicing in `cross_modal_concordance.load_rna_table` (`D4-011`),
    `fin_both_mask` applied in paired WNN modality weight audits (`D4-013`), and
    `verify_reproducible_capsule` failing closed on unverified environment metadata
    (`D5-021`).
- **Bulk multi-file `count_union_policy` defaults to `"error"` (`D4-001`, `D4-005`, `73e9f7b`, `1624980`):**
  - Multi-file bulk count loading (`aria/scripts/rna_bulk/counts.py`, `rna_bulk_de.py`,
    `BulkRNAAgent`, `script_contracts.py`) no longer silently outer-joins mismatched
    gene sets with `fillna(0)`. Ensembl gene version suffixes are normalized before
    comparison; mismatched gene universes now **fail closed by default**
    (`count_union_policy="error"`, `CountGeneSetMismatch`). Callers that explicitly
    want non-identical gene sets merged must opt in to `count_union_policy="intersection"`
    or `"zero_fill"`, which `BulkRnaNarrator` surfaces as a warning caveat. Headered
    per-sample count files and colliding subdirectory stems (`dirA/counts.tsv`,
    `dirB/counts.tsv`) are now parsed and disambiguated cleanly (`D4-005`).
- **Single 0-based half-open BED / 1-based closed GTF conversion — TSS distances shift by 1 bp (`D4-014`, `D4-016`, `D4-017`, `D8-004`–`D8-007`, `D8-011`):**
  - All GTF parsing, peak-coordinate parsing, and UCSC/Ensembl contig normalization
    converged onto pure-stdlib `aria/utils/genomic_coords.py` (`gtf_to_bed_interval`,
    `gtf_tss_0based`, `parse_gtf_tss_records`, `parse_gtf_gene_models`,
    `parse_gtf_symbol_map`, `parse_peak_coord`, `canonical_contig`, `ensembl_contig`).
  - 1-based closed GTF coordinates are converted **once at parse time** to 0-based
    half-open BED `[start - 1, end)`: `+`-strand TSS shifts from `S` to `S - 1` and
    `-`-strand TSS shifts from `E` to `E - 1`. Reported `distance_to_tss` and
    `tss_offset_bp` values therefore shift by **1 bp** compared to `v4.9.4`, and
    exon-boundary peaks are classified accurately (`Exonic` instead of `Intronic`).
  - `peak2gene_core.compute_peak_gc_from_indexed_fasta` no longer subtracts `1` a
    second time from already 0-based BED starts (`D4-017`, `62ee32a`). Contig
    normalization (`chr1` ↔ `1`, preserving unplaced scaffolds such as `chrUn_*`)
    now applies consistently across CRAM counting, broad-domain controls, consensus
    peak merging, and IP-vs-control JSD (`220fff2`, `62ee32a`, `7cb6a7b`).
- **scRNA mitochondrial detection no longer fabricates `0%` (`D8-002`, `b79b9bb`, `ea49754`):**
  - `aria/utils/ambient_qc.py::detect_mt_mask` is now shared by `rna_qc.py` and
    `integration_paired_qc.py`, resolving mitochondrial genes across both `var_names`
    and `var` symbol columns (`VAR_SYMBOL_COLS`) for Ensembl-indexed datasets while
    preserving the case-insensitive `^mt-` fallback warning.
  - When zero mitochondrial genes match in `rna_qc.py`, ARIA no longer fabricates
    `0%` mitochondrial content or runs a no-op filter: it omits `qc_vars=["mt"]`,
    skips mitochondrial filtering, and explicitly records `mt_measured=False`,
    `mito_fraction=None`, `mt_threshold_used=None`, `mito_qc_degraded=True`, and
    `metrics_not_computed=["mito_fraction (...)"]`, surfaced truthfully in Methods
    and via `warning_router.py`.
- **No positional design fallback (`D7-001`, `4856aa3`):**
  - `BulkRNAAgent._apply_design` removed its silent positional fallback when matrix
    column names do not match design sample IDs (which could invert DE `log2FC`
    signs when column order differed from design order). Explicit `sample_id_map`,
    `sample_aliases`, and GEO SOFT `GSM↔SRR` (`sra_accessions`) maps are resolved;
    unmatched IDs without an explicit mapping now **fail closed** with a `ValueError`
    naming unmatched IDs on both sides.
- **`sample_covariates` now forwarded to cross-modal concordance (`D7-003`, `cc40b8a`) and layer contracts aligned (`D7`):**
  - `OrchestratorAgent._run_cross_modal_concordance` now forwards every optional
    input declared in `SCRIPT_CONTRACTS["aria/scripts/cross_modal_concordance.py"]`
    — including `sample_covariates` — to Layer D SVD OLS residualization instead of
    silently dropping them.
  - Shared `resolve_design_contract` (`aria/utils/design_matrix.py`) now drives
    `devils_advocate`, `audit_agent`, `design_intelligence`, `decision_manifest`
    (retaining `batch_covariate` in `design["covariates"]`, `D7-013`),
    `falsification_certificates` (`D7-008`), `hypothesis/adapters/bulk_rna` (`D7-012`),
    `ChromatinAgent` TOBIAS `paired_replicates` (`D7-009`), and `IntegrationAgent`
    MOFA+ `condition_col`/`batch_key` (`D7-021`). Additional contract fixes: actual
    fitted formula in `rna_bulk_de.py` `design_used` (`D7-007`), `devils_advocate`
    reading `qc["scrublet"]["ran"]` directly (`D7-002`), `cross_modal_triangulation`
    reading `da_peaks_by_cluster` (`D7-006`), `hypothesis/adapters/hic` falling
    through to `qc["n_samples"]` (`D7-011`), `hypothesis/adapters/cross_modal` using
    `LOW_POWER_MIN_MATCHED_SAMPLES = 10` and `pearson_r` polarity matching `padj`
    (`D7-014`, `D7-015`), and `chromatin_qc` standardizing on `subnucleosomal_fraction`
    / `mononucleosomal_fraction` (`D7-019`).
- **Per-narrator isolation and governance-pass failures surfaced in the report (`D5-001`–`D5-025`, `5b7af89`..`d57ad7b`):**
  - `NarrativeRegistry.collect_blocks` isolates per-narrator exceptions into a
    modality-scoped `<modality>.narrator_error` block while healthy narrators
    continue, and `NarrativeAgent` no longer swallows outer block-collection errors
    (`D5-001`). Crashed agents emit structured error dicts routed through `W-WARN`,
    `NarrativeRegistry`, and `RunLedger` (`D5-005`), and `BiologicalSynthesisAgent`
    slice failures emit `integration.synthesis_error.<slice>` error blocks (`D5-006`).
  - Failures in governance passes — `RunLedger` (`D5-002`), `devils_advocate`
    (`D5-003`), robustness multiverse / falsification certificates (now decoupled,
    `D5-004`), `compile_public_claims` (`D5-009`, `D5-010`, tagged
    `compiler_internal_error` with withheld counts and explicit unclassified tier
    badges `D5-025`), triangulation and claim-linkage (`D5-018`, `D5-023`),
    speculative hypotheses (`D5-020`), `AuditAgent` checks (`D5-013`),
    `experiment_view` read models (`D5-022`), and figure/artifact staging (`D5-024`,
    recorded in `figures._errors`) — now surface explicitly in `methodology.json`
    and HTML error banners instead of failing silently.
- **LOBO errors reported as `error` in falsification certificates (`D1-014`, `1377472`, `0c3f224`, `3838a77`):**
  - `build_falsification_certificates` no longer masks failed LOBO jackknife runs
    (`status="error"`) as `"eligible_not_executed"`. Failed contrasts now report
    `status="error"` with `n_contrasts_errored` and preserved `error_details`.
- **apeGLM shrinkage status disclosed truthfully in Methods (`D1-005` part 2, `6eee920`; `D8-001` / `D2-BUG-6`, `8f0713b`):**
  - `rna_pseudobulk_de.py` now records realized `lfc_shrinkage` counts (`applied`,
    `n_applied`, `n_failed`, `failed_reasons`), and `scrna/_common.py::_lfc_shrinkage_clause`
    discloses full, partial, or unshrunken MLE status truthfully in Methods instead
    of claiming `"apeGLM-shrunken"` whenever shrinkage was merely requested.
  - Spike-in DESeq2 (`aria/utils/epigenomic_spikein_de.py`) now delegates to
    canonical `rna_bulk/deseq2.py::_run_deseq2(..., external_size_factors=...)`
    with `C(condition, Treatment(reference=...))`, fixing apeGLM coefficient
    lookup on alphabetical contrasts (`KO` vs `WT`).
- **Centralized warning router (`W-WARN`) caveats (`e1d4792`, `0a885b9`, `0670449`, `247c0e6`, `e53d4a4`):**
  - Pure-stdlib `aria/agents/narrative/warning_router.py` (`296` lines, zero I/O,
    zero `except Exception`) runs at the single `NarrativeRegistry.collect_blocks`
    choke point (plus `cross_modal_synthesis.py` sample-matched blocks), retiring
    `9` manual warning loops in `HiCNarrator` and `ChromatinNarrator` and providing
    uniform routing across all 5 narrators without introducing any new `status`
    value.
  - Structured signals (`error` on failed stages with `status in {"error", "failed"}`,
    `*_error`, `figures._errors`, `failed_databases` e.g. partial/total Enrichr
    outage `D1-025`, `failed_samples`/`failed_modalities`, `*_degraded`,
    `phasing_fallback_reason`, `lfc_shrinkage` discrepancy,
    `covariate_adjustment.unadjusted_confounding_risk`,
    `sex_chromosome_significant_fraction >= 0.20`, `pass_qc=False` with
    `qc_failures`, `pairing_status='partial'` with dropped donors,
    `batch_condition_confounded`, `n_nonconverged_excluded > 0`) emit structured
    caveats, may cap `confidence` (`high → medium`), or synthesize orphan
    `limitation`/`error` blocks.
  - Free-text stage `warnings` route strictly by structural position (zero keyword
    matching) with `severity="warning"` and never alter `confidence` or `status`;
    root warnings attach only to the QC/summary block or fall back to
    `blocks[0].metadata["degradation_signals"]` (`severity="info"`). Clean E2E
    runs now surface previously dropped chromatin `qc`/`peak_ora` warnings and
    multiome `metrics_not_computed` (`info`) while collapsing duplicate warnings
    (`0` status or confidence changes across `62` audited E2E blocks).
- **Canonical BH-FDR and statistical bug fixes (`D8-016`, `D2-BUG-1..4`, `D2-BUG-8..10`, `D8-014`, `D8-015`):**
  - All `7` local BH copies converged onto `aria/utils/stats.py::bh_correct`
    (`2d7c995`), masking `NaN`/`Inf` inputs before adjustment so finite p-values
    are adjusted with the true finite denominator $m$ instead of poisoned to
    all-`NaN` (`rna_pseudobulk_de._global_bh`, `synthetic_de._bh_adjust`).
  - `check_calibration` fails closed on `NaN`/`Inf`/`bool` metrics (`0098611`);
    Layer D `sample_matched_layer` drops strictly constant nonzero features before
    `log2` (`e110515`); `fit_copula_idr_pair` equalizes exact tie groups for
    permutation invariance and replicate-swap symmetry (`2771948`);
    `peak2gene_core` quantile binning uses tie-safe average ranks invariant to
    peak order (`58b21bc`); `-log10(p)` underflow plot capping converged onto
    `stats.adaptive_neglog10` (`min_cap=5.0`, `3f54024`, `ea49754`); CPM/log2-CPM
    converged onto `stats.compute_cpm` / `compute_log2_cpm` (`c91f692`).
- **Narrowed pseudobulk warning filter & global log2 filter removed (`D9-001`, `D9-002`, `ad4ae7d`, plus `D9-004`–`D9-009`):**
  - Removed blanket `warnings.filterwarnings("ignore")` in `rna_pseudobulk_de.py`
    (narrowed to benign categories only; real numerical `RuntimeWarning`s are now
    captured in per-comparison `fit_warnings` and top-level `warnings`) and
    removed the global `.*log2.*` `RuntimeWarning` filter in `aria/scripts/_base.py`.
  - Also fixed `LLMProvider.from_config` to honor `ARIA_CONFIG_FILE` and forward
    `extra_params` (`20afb49`), unified GTF env resolution across `ARIA_GTF` /
    `ARIA_GTF_PATH` + `genome_dir()`, wired `propose_illumina_scatac_manifest` into
    `DataAuditAgent` CP1 previews, removed the undocumented
    `ARIA_ALLOW_BULK_METADATA_INFERENCE` bypass, and surfaced `n_ambiguous_direction`
    in cross-modal TF synthesis (`2aebfb6`).

---

## Test, CI and Guard Infrastructure

- **Statistical oracle suite (`tests/oracle/`, `0c8ca78`, `5ee6d22`, `3b9a19c`, `78f5558`, `ca8532f`):**
  Added `10` analytical oracle modules (`182` items) comparing ARIA's core
  statistics against independent references (`statsmodels`, `scipy.stats`,
  `cooler`/`cooltools`, closed-form Smyth OLS, Signac TF-IDF/`log1p_cp10k`,
  bivariate Gaussian copula E-step) and verifying sample/feature permutation and
  level-renaming invariance. All `10 / 10` exposed `D2-BUG-1..10` strict xfails
  now pass (`0` xfails remaining in `tests/oracle/`).
- **Mutation score `70.2% → 82.75%` (`D2`):**
  `mutmut 3.8.1` across `12` core statistical modules (`6,945` scored mutants)
  improved from `70.22%` (`4,877 / 6,945`) to **`82.75%` (`5,747 / 6,945`,
  `+12.53 pp`)** — `11 / 12` modules `≥ 80.1%` and `peak2gene_core.py` at `79.21%`
  raw (`87.3%` excluding equivalent/cosmetic mutants) — complemented by **`27 / 27`
  (100%)** empirical manual mutations killed (`22 / 22` D3 + `5 / 5` heavy-module
  `pydeseq2`/`cooltools` mutations).
- **CI collection `46.8% → 100%` with new Conda lanes (`48b5bad`, `af1aad5`, `bed54c0`, `1884076`, `4e70292`):**
  Replaced hand-maintained file allowlists in `.github/workflows/ci.yml` with
  dynamic `pytest tests/` across `pr-checks`, `rna-env-tests`, new
  `hic-env-tests` (`envs/aria-hic-env.yml`), new `integration-env-tests`
  (`envs/aria-integration-env.yml`), and real in-container test runs in all 4
  Docker images under `release-repro`. Native test-file collection jumped from
  `126 / 269` (`46.8%`, `143` uncollected) to **`100%`** (`281 / 281` at D3 close;
  all `307 / 307` non-ignored native files on HEAD; `6` legacy standalone scripts
  in `tests/conftest.py:collect_ignore`). Resolved the 4 Conda-lane blockers
  surfaced by full-suite execution (`mudata_io._stage_cache` cold-miss lockfile,
  pseudobulk composition collinearity fixture, `episcanpy` `tests/` package
  shadowing, and `[tui]` extra gating on cockpit smoke).
- **`tests/test_ci_collects_all_tests.py` (`48b5bad`):**
  AST + YAML guard failing CI whenever any `tests/**/test_*.py` file is neither
  collected by a CI job nor listed with a concrete reason in
  `EXPLICIT_CI_EXCLUSIONS` matching `tests/conftest.py:collect_ignore`.
- **`tests/test_debt_ratchet.py` + `tests/debt_ratchet.json` (`D10`, `642f3e0`):**
  Zero-dependency stdlib-AST guard enforcing monotonic debt reduction across
  `aria/` (`slack = 0`; raising any ceiling requires an explicit ADR ID):
  - `unannotated_broad_except`: **`507 → 313`** (limit `313`; `162` justified
    handlers explicitly annotated `# aria: allow-broad-except <reason>`).
  - `silent_broad_except`: **`148 → 99`** (limit `99`).
  - `zip_without_strict` (`aria/scripts/` + `aria/utils/`): **`67 → 0`** (limit `0`).
  - `functions_over_300_lines`: **`33 → 35`** (limit `35`; targeted by deferred D6).
- **`tests/test_no_agent_imports_scripts.py` (`1bfdff8`):**
  AST layering guard enforcing `0` top-level `aria.scripts` imports under
  `aria/agents/` and freezing the `4` pre-existing lazy edges so shared
  agent/script contracts remain in pure utility modules (`aria/utils/`).

---

## Receipts and Preprint Impact

- **v49 comparative receipt regenerated (`876f3ea`):**
  Regenerated `docs/benchmark_results/v49_comparative_benchmark/comparative_benchmark_v49.json`
  (and `docs/vignettes/13_governed_vs_uncontracted_benchmark.md`) via the
  canonical command after `D7-013` fixed `decision_manifest._normalize_design` so
  `batch` is recorded in `design.covariates` on the 6 batch capsules. All `40`
  capsule outcomes, arm metrics, 20-seed sweep scores, and tamper detection
  (**`40 / 40`**) are byte-identical; only the 6 capsule manifest digests and
  receipt SHA changed (`b29d5911… → b685e657…`).
- **B2 byte-reproducible across `PYTHONHASHSEED` (`690b961`):**
  Sorting token sets in `evidence_verifier.py` and `*_rejection_reasons` in
  `governance_b2.py` eliminated hash-seed-dependent list ordering. `score_b2()`
  now hashes identically (`9eefa7a7…`) under `PYTHONHASHSEED=0` and
  `PYTHONHASHSEED=1`. Relative to `v4.9.4` (`fed6e138…`), the diff is strictly
  alphabetical ordering inside 6 `per_claim` rejection-reason lists; every score,
  CI, tier, emission, and verdict is **identical**.
- **Historical `adr079` and `adr080` receipts NOT recomputed:**
  `docs/benchmark_results/multiome_concordance/adr079_multidonor_multiome_validation.json`
  and `adr080_cross_modal_concordance_validation.json` remain historical receipts
  and were **not** recomputed for `v4.9.5`; code paths exercised by those
  harnesses changed in `v4.9.5` via `D4-011` (ambiguous duplicate RNA symbol
  exclusion in `load_rna_table`), `D4-017` (0-based BED start fix in
  `compute_peak_gc_from_indexed_fasta`), and `D7-003` (`sample_covariates`
  forwarding in `_run_cross_modal_concordance`). External Tier 1 benchmarks
  (`docs/benchmark_results/external_benchmarks_v49/*`) remain untouched and
  attributed to immutable tag `v4.9.0`.
- **19-lane `preprint_v1` freeze:**
  All `19 / 19` `automated_v1_gate` machine freeze lanes under
  `docs/benchmark_results/preprint_v1/` will be refreshed from a clean commit for
  `v4.9.5` (`3 / 22` human-gold lanes remain `deferred_to_v2_reader_study`).

---

## Deferred

- **D6 monolithic function decomposition execution:** The complete 6-wave plan
  for all `19` targets (`17` functions `≥ 400` lines, `15` blocks `CC ≥ 100`;
  `118` pure zero-I/O extractions + `29` I/O adapters) is delivered and reviewed,
  with execution deferred post-`v4.9.5` so the `60` missing pre-refactor
  analytical/golden-SHA-256 oracle tests are written first and line ranges are
  regenerated on top of the F1–F4 audit changes.
- **B2 / C4 manuscript decision:** Pending editorial decision on updating the
  preprint manuscript (and `verify_paper_claims.py`) to adopt the honest `v4.9.4`
  label-blind B2 component ablation numbers (`false_narrative_rate = 0.1875`,
  `9/48` vs. historical `0.1042`, `5/48`; causal overreach `6/6` / `5/6 → 0/6` at
  `16/18` clean coverage) and document the three named deterministic verifier
  boundaries and `missing_caveat_rate` metric definition.

---

## Known Limitations

- **Remaining broad catches and long functions:** `313` unannotated broad
  `except` handlers (`99` silent) outside the D1/D5 audit scopes and `35`
  functions `> 300` lines remain in `aria/`; all are now capped with `slack = 0`
  by `tests/test_debt_ratchet.py`.
- **Free-text warnings on zero-block stages (`N-R3-1`):** A free-text `warnings`
  entry on a non-root stage when an agent emitted zero matching blocks cannot be
  attached to a non-existent block; structured failure/degradation fields still
  synthesize explicit orphan `limitation`/`error` blocks.
- **Open-time FASTA `OSError` ambiguity (`D1-020` / `R3`):** Because
  `pysam.FastaFile(...)` raises identical `OSError` at open/index time for valid
  plain-gzip FASTAs, read-only directories without `.fai`, and unparseable
  headers, open-time `OSError` falls back to GTF phasing with an explicit
  `phasing_fallback_reason` warning, whereas post-open truncated-sequence
  `fetch()` errors fail closed.
- **Pre-existing non-audit boundaries unchanged:** `3 / 22` human-gold reader
  study lanes remain `deferred_to_v2_reader_study`; air-gapped mode remains
  `Partial` for the pre-CP1 biological-question LLM parse (`A1b`); the 3 named
  deterministic `evidence_verifier.py` boundaries (English word numbers,
  superlative `top_term` checks, qualitative hype lexicon) remain documented
  boundaries; `6` legacy standalone scripts remain in
  `tests/conftest.py:collect_ignore`; and `4` pre-existing lazy `aria.scripts`
  imports in `aria/agents/` remain frozen by `tests/test_no_agent_imports_scripts.py`.
