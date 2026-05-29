---
status: active
source_of_truth_for: pre_atac_remediation
last_updated: 2026-05-28
---

# ARIA Senior Audit + Pre-ATAC Remediation Plan (2026-05-28)

Audit performed at HEAD `c611e45` (post-v4.5.2 hardening) from a dual lens:
senior software engineer + senior computational-biology/bioinformatics
researcher. Every finding below was grounded in reading the implementation;
one hypothesis was falsified by an executed probe (see F-SCI-COMPOSITION).

**Mandate from Samael:** solve these before moving to v4.6 scATAC.

This file is the authoritative remediation tracker. Each item carries:
severity, evidence (`file:line`), and the exact execution steps. Status is
updated in place as items land.

---

## Verification note (why this audit trusts itself)

The composition-covariate concern was initially flagged as a HIGH bug
(continuous log-ratio passed to `pydeseq2` in `design_factors` without
`continuous_factors` → suspected one-hot rank-deficiency). A synthetic probe
in `aria-rna-env` disproved it: this pydeseq2 version auto-detects numeric
dtype as continuous (design matrix = `[Intercept, condition[T.STIM],
_aria_composition_log_ratio]`, one column). Reclassified to a MINOR
robustness item (F-SCI-COMPOSITION). Lesson recorded: verify before asserting.

---

## Findings

### Software / architecture

- **F-ENG-E2E [HIGH] — the checkpoint state machine has no integration test.**
  The orchestrator → `DesignAgent` (CP2.1–2.6) → CP2-plan → dispatch control
  path is the most complex, most user-facing code and is untested.
  `tests/test_pytest_smoke.py:999` only tests `_format_plan_summary` via
  `__new__`. Fragile, real subtleties are uncovered: `checkpoint=2` means
  three different things (design-confirm store, plan escalation, post-audit),
  and `options[0]` is NOT the safe answer at CP2.3 (factor: `options[0]`
  can be `"genotype"` while the correct factor is `"stim"`/`"stimulation"`).
  Evidence: `aria/agents/orchestrator_agent.py:168-289`,
  `aria/agents/design_agent.py:165-391,467-503`.

- **F-ENG-REPRO [MEDIUM-HIGH] — reproducibility is documented, not enforced.**
  No `.github/workflows/`, no `Dockerfile`/Apptainer. Lockfiles are
  snapshot-from-installed (no real solve) and only `aria-rna-env` has a lock;
  chromatin/hic/integration have none. `--reproducible` + provenance hashes
  are good but depend on one machine. Evidence: `envs/` contents,
  absence of CI config.

- **F-ENG-BUS [MEDIUM] — the message bus is in-memory only.**
  `MessageBus._log` is a `deque` with no persistence
  (`aria/bus/message_bus.py:91`). A crash mid-run loses findings/decisions of
  completed stages; resume only covers file-valid h5ad stages, not the
  narrative/finding layer. Post-mortem audit and replay are incomplete.

- **F-ENG-PERF [MINOR-MEDIUM] — avoidable per-call overhead.**
  `EnvironmentManager.run_in_stack` calls `check_environments()` (which spawns
  `conda env list --json`) on every script dispatch
  (`aria/utils/environment_manager.py:400`). `conda run` re-activates the env
  per call. `bus.get_pending_checkpoints()`/`get_log()` scan the full deque
  (up to 100k) on every 0.5s poll tick.

- **F-ENG-GODFILES [MEDIUM] — maintainability hotspots.**
  `narrative_agent.py` 2723, `rna_bulk_de.py` 2431, `_narrative_scrna.py`
  2200, `scrna_agent.py` 1983 LOC. High merge/regression risk.

- **F-ENG-TESTGAPS [MEDIUM] — untested modules.**
  `setup_agent` (719 LOC, env detection), `rna_quantify`/`rna_align` (raw
  FASTQ path), `base_agent`, `rna_inject_condition`, `rna_figure_paga` have no
  test reference.

### Scientific / bioinformatics

- **F-SCI-POWER [MEDIUM] — reported power is optimistic vs the decision
  threshold.** `power_estimation.py` computes Wald-NB power at nominal
  α=0.05, but significance is declared via **global BH** over ~127k tests
  (`rna_pseudobulk_de.py:504`). The effective per-test α is far below 0.05, so
  reported power (8–13% in the d72 PBMC run) overstates the true power against
  the applied bar. Must be disclosed in Methods or reconciled.

- **F-SCI-LOGNORM [MEDIUM] — count "recovery" is not surfaced in provenance.**
  `rna_pseudobulk_de.py:165-235` reverses `round(expm1(x)·lib/1e4)` when
  `raw.X` is log-normalized. Clever, but: (a) the return dict does NOT record
  that counts were reverse-engineered → reports may not disclose it;
  (b) invalid if the input was scaled/regressed (ScaleData), only weakly
  guarded by an 85%-integer probe on 200 cells; (c) assumes scale_factor=1e4.
  Borders on the "no fabricated scientific outputs" principle.

- **F-SCI-CAUSAL [MEDIUM] — the causal-language guard is a 10-phrase
  blocklist.** `aria/agents/narrative/validators.py:15-26` scans only
  `claim` + evidence labels for a fixed list ("drives", "binds to", …). It
  misses most causal phrasing ("induces", "triggers", "causes", "controls",
  "leads to", "master regulator", "upstream of", "promotes", "responsible
  for") and does NOT scan the composed prose / LLM integrated interpretation,
  which is exactly where causal overreach appears. Weakest possible
  implementation of a core project principle.

- **F-SCI-FDR [LOW-MEDIUM] — global FDR pooling is opinionated.**
  Pooling every gene×celltype×contrast p-value into one BH family
  (`rna_pseudobulk_de.py:495-515`) is defensible but conservative for small
  blocks and couples unrelated cell types. Rationale should be documented; an
  IHW / per-contrast / hierarchical option would strengthen it.

- **F-SCI-LIANA [MINOR] — `n_perms=100` default is low** for stable LIANA
  permutation specificity p-values (`aria/scripts/rna_cellcomm.py:68`);
  1000 is typical.

- **F-SCI-COMPOSITION [MINOR] — composition covariate relies on pydeseq2
  dtype auto-detection.** Works today (verified) but is fragile to version /
  dtype changes; passing `continuous_factors=[COMPOSITION_COL]` explicitly is
  defensive. `rna_pseudobulk_de.py:431-436`.

### Missing capability

- No ambient-RNA correction (SoupX/decontX); only Scrublet doublets.
- No integration-QA metrics (LISI/kBET/silhouette) to validate Harmony.
- No effect-size shrinkage (apeglm/ashr) or s-values; raw DESeq2 LFCs only.
- No real RNA velocity (scVelo); PAGA/DPT only (honestly caveated).
- scATAC/bulk-ATAC/Hi-C/integration scaffolded (roadmap v4.6–v4.8).

---

## Remediation plan — execution steps per proposal

Status legend: ☐ planned · ◐ in progress · ✅ done.

### P1 — Headless runner + E2E test of the control path  (closes F-ENG-E2E)

Highest value/risk ratio. Reuses the validated `/tmp/aria_pbmc_headless.py`
driver written for the PBMC rerun.

Steps:
1. ✅ Add `aria/headless.py`: a maintained non-interactive runner. It mirrors
   `tui.run_analysis` but resolves checkpoints with an auditable answer policy
   that follows ARIA's own inference (parse proposed factor, accept inferred
   groups/organism, recommended+optional plan) — no hardcoded biology
   (ADR-011). Exposes `run_headless(data_dir, question, answer_policy=...)`.
2. ✅ Add `tests/test_headless_design_e2e.py`: drive the full `DesignAgent`
   state machine (start_design → groups → organism → factor → batch/pseudorep
   auto-skip → confirm → `_build_design`) from a synthetic high-confidence
   `inferred_design`, asserting the final design (condition/replicate/groupby,
   formula, n_total). Pure Python, no LLM, no subprocess. Neutral synthetic
   labels only.
3. ✅ Validate: `pytest -q tests/test_headless_design_e2e.py` + full smoke.

### P3 — Honesty layer hardening  (closes F-SCI-CAUSAL)

Steps:
1. ✅ Broaden `CAUSAL_PATTERNS` to cover the common causal verbs/nouns above.
2. ✅ Scan composed prose: have `compose_prose`/`render_blocks` route the final
   per-block prose through the causal guard, not just `claim`+evidence.
3. ✅ Add `tests/` cases: a block whose prose contains "induces"/"master
   regulator" gets a caveat + confidence downgrade.

### P4a — Lognorm-recovery provenance  (closes F-SCI-LOGNORM)

Steps:
1. ✅ `rna_pseudobulk_de.py`: add `count_source`
   (`raw_counts`|`X_counts`|`recovered_from_lognorm`), `lognorm_recovered`
   bool, and `norm_scale_factor_used` to the return dict.
2. ✅ scRNA narrator emits a visible caveat when `lognorm_recovered` is true.
3. ✅ Test: synthetic log-norm h5ad → return dict flags recovery.

### P4b — Power honesty  (closes F-SCI-POWER)

Steps:
1. ✅ Surface `power_alpha` and an explicit note that power is computed at the
   nominal α, while significance uses global BH (so power is an upper bound).
   Add to the pseudobulk/bulk power payload and Methods prose.
2. ✅ Test: power payload includes the disclosure field.

### P-ENG-PERF — cache env detection  (closes F-ENG-PERF, env part)

Steps:
1. ✅ Cache `check_environments()` result on the `EnvironmentManager` instance
   (invalidate on alias-file mtime change). Keep correctness for SetupAgent.
2. ✅ Covered by existing env tests + a new cache test.

### P2 — CI + container artifacts  (closes F-ENG-REPRO, partially)

Cannot run CI/Docker build in this environment; delivered as committed
artifacts that run on the GitHub remote.

Steps:
1. ✅ `.github/workflows/ci.yml`: matrix that sets up `aria-env`, runs
   `compileall` + `tests/test_pytest_smoke.py` and the narrative suite; a
   second job (allowed-to-fail until envs are containerized) for the
   pydeseq2-gated tests.
2. ✅ `Dockerfile` (or `Apptainer.def`) skeleton per the rna stack, documented
   in `docs/`. Build is deferred to Samael's machine / CI.

### Deferred to a dedicated v4.5.x patch (documented, not executed now)

These need package installs / heavy validation and must not be done blind:

- ☐ **F-SCI-FDR**: add IHW/hierarchical FDR option behind a flag; document the
  global-pooling rationale in Methods.
- ☐ **Effect-size shrinkage** (apeglm/ashr) + s-values in `rna_pseudobulk_de`
  and `rna_bulk_de`.
- ☐ **F-ENG-BUS**: persist bus events to the SQLite memory as an event log so
  resume reconstructs narrative state.
- ☐ **Integration QA** (LISI/kBET/silhouette) + a report "analysis confidence
  score" closing the DesignIntelligence loop.
- ☐ **Ambient RNA** correction (decontX/SoupX) in `rna_qc`.
- ☐ **F-ENG-GODFILES**: split `narrative_agent.py` / `_narrative_scrna.py`.
- ☐ **F-ENG-TESTGAPS**: native pytest for `setup_agent` + FASTQ path
  (overlaps the existing P2 test-debt item).
- ☐ **F-SCI-LIANA**: raise default `n_perms` to 1000 (perf-validate first).

---

## Execution log (2026-05-28)

Executed and validated this session (all green, zero regressions —
`tests/test_pytest_smoke.py` 86 passed/4 skipped plus the new files):

- **P1 ✅** `aria/headless.py` (maintained non-interactive runner) +
  `tests/test_headless_design_e2e.py` (2 tests: full CP2.1–2.6 walk asserting
  the confirmed design uses the inferred factor not `options[0]`; reject-at-
  confirm cancels). Closes F-ENG-E2E.
- **P3 ✅** Broadened `CAUSAL_PATTERNS` (~25 phrases), added
  `find_causal_language`, and made `render_blocks` scan the FINAL composed
  prose (not just claim/evidence). `tests/test_causal_guard.py` (5 tests).
  Closes F-SCI-CAUSAL.
- **P4a ✅** `rna_pseudobulk_de` now returns `count_source`,
  `lognorm_recovered`, `norm_scale_factor_used`, `lognorm_lib_size_col`; the
  scRNA narrator emits a recovery caveat. New narrator test. Closes
  F-SCI-LOGNORM.
- **P4b ✅** `rna_pseudobulk_de` returns a `power` disclosure block
  (nominal-alpha vs applied global-BH threshold = upper bound). Closes
  F-SCI-POWER (data contract; methods prose surfacing tracked).
- **P-ENG-PERF ✅** `EnvironmentManager.check_environments()` is memoized on
  env-aliases mtime. `tests/test_env_manager_cache.py`. Closes F-ENG-PERF
  (env part).
- **P2 ◐** `.github/workflows/ci.yml` (light lane runs compileall + audit-gate
  tests + smoke; heavy lane = micromamba rna env, allowed-to-fail) and a
  `Dockerfile` (micromamba + `aria-rna-env.yml`). Build/CI run on the GitHub
  remote, not verifiable locally. Partially closes F-ENG-REPRO.

Headless-runner bug found and fixed during the PBMC rerun: the live loop
filtered `bus.get_log(experiment_id)`, but agent STATUS messages do not carry
`experiment_id`, so the narrative "Report saved" completion signal was hidden
and the runner timed out despite a finished report. Fixed to poll the full bus
log (like `tui._live_analysis_loop`) plus a filesystem fallback.

## PBMC v4.5.2 rerun outcome (BLOCKER FIXED — full report rerun still useful)

The headless PBMC rerun produced
`~/.aria/reports/aria_20260528_165233_interferonbeta_myeloidcells_lymphoidcell_-1db/`.

Provenance is clean: `git_sha=0569689` (clean tree), `aria_version=4.5.2`,
input SHA-256 `af0696e9…` (matches canonical pbmc.h5ad), conda/pip locks
embedded, LLM usage 3 haiku calls / $0.0145 (real calls, honestly accounted).

**But the report is THIN and does NOT reproduce the v4.4 (`cbcde8e`/d72)
depth.** Only two success narrative blocks: QC (13,807/13,836 cells, 13
`obs['cluster']` groups) and LIANA (50 L-R interactions, low confidence).
**Pseudobulk DE, differential abundance/composition, pathway ORA, and
trajectory did NOT run** — the decisions log has no pseudobulk/composition
decision (d72 had `composition_covariate=ON`, `groupby=cluster; condition=stim;
replicate=Donor`).

### DIAGNOSED 2026-05-28 — root cause + fix

Isolation verdict (evidence-based, not a rerun guess):

- **NOT cell-focus steering.** `methodology.design_intelligence.profiles[0].focus
  == []` — no focus subset was applied.
- **NOT a v4.4→v4.5.2 code regression.** The gate logic is unchanged; the
  behavior difference is driven entirely by question phrasing.
- **Root cause = silent plan/dispatch disagreement (latent in v4.4 too).**
  `scRNAAgent._needs_pseudobulk` re-gated pseudobulk DE on
  `PSEUDOBULK_KEYWORDS` matched against the free-text question, and that gate
  **overrode** DesignIntelligence's recommendation + the user's CP2 approval.
  The blocker question ("interferon-beta transcriptional response programs /
  myeloid-lymphoid signaling networks") contains NONE of the keywords
  (`compare`, `versus`, `differential`, `control`, `treatment`, …), so
  `_needs_pseudobulk` returned False and DE was dropped **without a logged
  decision or a reported skip** — while DesignIntelligence had recommended
  `Donor/sample-level pseudobulk DE: condition=stim, replicate=Donor,
  groupby=cluster` and the user selected "recommended + optional". A
  d72-style question hits 4 keywords, which is why the v4.4 report ran DE.

**Fix (committed):** `_needs_pseudobulk` now honors the approved plan — when
there is an explicit obs design AND DesignIntelligence recommended pseudobulk,
it runs regardless of question keywords. The keyword gate remains only as a
fallback for designs with no DI recommendation (e.g. filename-inferred groups
on a purely descriptive question). Regression guard:
`tests/test_pseudobulk_gate.py` (3 tests).

Follow-up worth considering: a run that silently drops a recommended/approved
analysis is a report-integrity violation — the dispatch layer should at least
emit a visible "planned but not run" finding whenever it diverges from the
approved plan (candidate hardening, tracked).

## External AI audit (auditoria.txt) — adjudication 2026-05-28

Two external AIs audited ARIA **statically from GitHub (not executed, stale/
partial snapshot)**. Each claim below was checked against the real code. The
methodology caveat matters: several headline P0s are factually wrong because
the snapshot missed files that exist on `main`.

### REJECTED — false against the current tree (with evidence)

- **"orchestrator references missing agents design_agent/audit_agent/
  raw_ingestion_agent; central flow broken" (IA2 P0 2.1)** — FALSE. All exist
  (`design_agent.py` 761 LOC, `audit_agent.py` 441, `raw_ingestion_agent.py`,
  `design_intelligence.py` 371). The design state machine is now covered by
  `tests/test_headless_design_e2e.py`.
- **"DesignIntelligence absent" (5.3) / "no section narrators" (5.7)** —
  FALSE. `design_intelligence.py` and `narrative/narrators/{scrna,bulk_rna}.py`
  exist (the v4.5.2 Narrative Kernel already does evidence-card narrators).
- **"SQLite free concurrent writes will corrupt state" (IA1 1.1)** —
  OVERSTATED. `memory.py` already uses `check_same_thread=False` + `RLock`
  serialization + `WAL`, and dispatch runs in a single daemon thread. The real
  residual is durability, tracked as F-ENG-BUS.
- **"Quantum superposition agent architecture" (IA1 novel #4)** — REJECT the
  framing (buzzword); a modest reframe is accepted as X19.

Lesson: the external audits would not have produced these false P0s if a
registry-integrity test existed — which is itself a good accepted item (X3).

### CROSS-REF — already done or already in this plan

- lognorm recovery should be LOW/INSUFFICIENT, not a log line (IA1 sci #2) →
  **DONE this session (P4a)**; can additionally force a confidence downgrade.
- bus not durable (IA2 2.6) → F-ENG-BUS (deferred).
- ambient RNA, LISI/kBET, effect-size shrinkage, per-sample QC, FASTQ-path
  tests (IA2 4.4) → already in this file's deferred list.
- reproducibility must be transitive/cryptographic (IA1 systemic #3) → the
  Docker image digest in provenance is that answer; P2 (CI+Docker) started.
- pseudobulk design/dispatch must honor the plan → the PBMC-blocker fix.

### ACCEPTED — new items folded into the remediation plan

P0/cheap integrity (a "v4.5.3 Integrity & Trust" mini-milestone before ATAC):

- ✅ **X1** Centralized version in `aria/version.py`; `aria.__version__`,
  `aria.llm.__version__`, `setup.py`, TUI, and `install.sh` now read the same
  source. Version metadata bumped to `4.5.3` locally.
- ✅ **X2** Installer no longer writes API keys to `~/.bashrc`; keys are stored
  only in `~/.aria/.env` with `chmod 600` and exported to the current installer
  process for verification.
- ✅ **X3** Added registry-integrity checks/tests: import every
  `AGENT_REGISTRY` entry, validate modality validation metadata, and check
  required script contracts.
- ✅ **X4** Chromatin modalities are marked `scaffold` and dispatch-gated until
  v4.6+ scripts land; direct calls to planned missing chromatin scripts return
  structured `script_not_implemented` blockers.
- ✅ **X17** Added tiered `aria doctor` (`--smoke`/`--synthetic`/`--benchmark`)
  plus installer smoke verification. Installer no longer claims full readiness
  after mock-only integration checks.

Scientific contracts & validation:

- ✅ **X5** Typed IPC contracts (pydantic `ScriptContract`: input/output schema,
  required files, `validation_level`) at the `EnvironmentManager` boundary.
  Implemented in `aria.utils.script_contracts`; input/output contract failures
  now return structured `InvalidScriptParams` or `IncompatibleScriptContract`
  errors instead of obscure downstream failures.
- ☐ **X6** Synthetic ground-truth benchmark (splatter/scDesign3) with known
  true DE genes — numerical-accuracy regression, not just flow (IA1 missing #1,
  IA2 5.6). Overlaps the P2 test-debt.
- ✅ **X7** Design-matrix sanity validator before DESeq2: rank deficiency,
  batch↔condition confounding, n=1 blocks, continuous-vs-miscategorized factor
  (IA1 sci #1, IA2 4.4/5.3). Implemented in
  `aria.utils.design_matrix`, surfaced in DesignIntelligence/AuditAgent, and
  enforced before bulk/scRNA DESeq2 calls.
- ☐ **X8** Integration QA as red flags: negative silhouette / poor LISI/kBET
  must raise a flagged AuditAgent finding (overcorrection), not a passive
  report line (IA1 sci #3). Folds into the deferred integration-QA item.
- ☐ **X9** Annotation-reuse marker-coherence check: when reusing `obs` labels,
  verify canonical markers are consistent (avoid mislabeled atlas/species
  annotations) (IA1 sci #4). Mind ADR-011: use a generic, documented marker
  reference, flag-only.
- ☐ **X10** Privacy firewall: `ARIA_AIR_GAPPED` app-level network block + LLM
  prompt redactor (sample-name → token map, kept locally) + human-data cache
  policy (IA1 missing #2, IA2 2.5/5.5).
- ☐ **X11** DebateCouncil cost governance: cap rounds, gate debates to
  high-impact findings, hard budget ceiling (IA1 systemic #1).
- ☐ **X12** ParameterAdvisor bias guardrail: historical suggestions are
  provenance, not authority; never silently override scientific defaults
  (IA1 systemic #2).
- ☐ **X13** DataAudit scan limits: max files/depth/total-size, symlink policy,
  sampling for inference on huge directories (IA1 systemic #4, IA2 2.7).
- ☐ **X18** Spatial-transcriptomics scaffold (Visium/MERFISH); the `spatial`
  env alias already exists in `EnvironmentManager.STACKS` (IA1 missing #4).

Flagship / deferred research (post-integrity, high value):

- ☐ **X14** **Claim Compiler** — classify each biological claim as
  descriptive / associative / weak-mechanistic / strong-mechanistic /
  causal-experimental, with a per-claim manifest (claim_id → evidence paths →
  code commit → limitations → confidence). Evolves the hardened causal guard +
  the NarrativeBlock kernel (already has evidence/claim/caveats/confidence) into
  the real thing. Both AIs converge on this; it is the single highest-value
  differentiator (IA1 novel #3, IA2 4.5/5.1).
- ☐ **X15** Shadow / multiverse analysis: rerun key decisions (QC strictness,
  Harmony vs scVI, DESeq2 vs edgeR, thresholds) on a subset; report robustness
  ("conclusions stable across 5/6 choices") (IA2 5.4).
- ☐ **X16** Evidence/knowledge graph per run (dataset→samples→design→contrasts
  →scripts→artifacts→results→claims→report) for navigation, partial reanalysis,
  and auto-methods (IA2 5.2). Overlaps X14 + bus persistence.
- ☐ **X19** Calibrated-uncertainty hypotheses (reframed from "quantum"): when
  competing cell-state hypotheses exist, report a probability distribution +
  the discriminating experiment, instead of forcing a binary DebateCouncil
  verdict (IA1 novel #4, reframed).
- ☐ **X20** Knowledge-synthesis agents: Null-Result-of-Interest + LitCross
  (local PubMed embeddings) (IA1 novel #1/#2). Most speculative; research track.

Net: the external audits add real value in **X1–X14** (integrity, contracts,
privacy, claim compiler). Their P0 "broken core" narrative is rejected as a
stale-snapshot artifact; the registry-integrity test (X3) is the durable
defense against that class of false alarm.

## Validation gates (run after each executed item)

```bash
/home/medusa/anaconda3/envs/aria-env/bin/python -m compileall -q aria
/home/medusa/anaconda3/envs/aria-env/bin/python -m pytest -q tests/test_pytest_smoke.py
# narrative + new targeted tests
```

## v4.5.3 Integrity & Trust execution log (2026-05-28)

Implemented locally after commit `27d4b15`:

- X1-X4 and X17 closed as described above.
- New files: `aria/version.py`, `aria/doctor.py`,
  `aria/utils/registry_integrity.py`,
  `tests/test_registry_integrity.py`, `tests/test_doctor.py`,
  `tests/test_chromatin_dispatch_gate.py`.
- Validation:
  - `python -m compileall -q aria` -> pass
  - `python -m pytest -q tests/test_registry_integrity.py tests/test_doctor.py tests/test_chromatin_dispatch_gate.py` -> 8 passed
  - `python -m pytest -q tests/test_chromatin_agent.py tests/test_env_manager_cache.py` -> 1 passed
  - `python -m pytest -q tests/test_pytest_smoke.py` -> 86 passed / 4 skipped
  - Combined targeted + smoke -> 94 passed / 4 skipped
  - `python -c "import sys; sys.argv=['aria','doctor','--smoke']; from aria.tui import main; main()"` -> passed with one warning: HiC is still scaffold but dispatch-enabled from earlier behavior.
  - `git diff --check` -> pass

## X7 design-matrix validator execution log (2026-05-28)

Implemented on top of `22c93ac`:

- New `aria/utils/design_matrix.py` validates sample-level DE designs for:
  rank-deficient model matrices, complete condition-covariate confounding,
  insufficient condition replicates, n=1 condition x covariate cells, numeric
  condition factors, binary numeric covariates, and text-stored continuous
  covariates.
- `AuditAgent` now runs the validator before dispatch when confirmed design
  metadata can be mapped to count-matrix columns; blocking design issues enter
  CP3.5 rather than failing inside DESeq2.
- `DesignIntelligence` surfaces early design-matrix warnings from confirmed
  groups plus batch/covariate maps.
- `rna_bulk_de.py` validates each contrast before DESeq2 and returns
  `InvalidDesignMatrix` instead of opaque model errors.
- `rna_pseudobulk_de.py` validates each cell-group contrast before DESeq2,
  skips invalid blocks with structured `design_matrix_invalid`, records
  `design_check`, and passes explicit `continuous_factors` to pydeseq2 when
  the installed API supports it.
- `ScrnaNarrator` surfaces design-matrix warnings as visible pseudobulk caveats.

Validation:

- `python -m compileall -q aria` -> pass
- `python -m pytest -q tests/test_design_matrix_validator.py tests/test_narrator_scrna.py` -> 8 passed
- Narrative + X7 suite (`test_narrative_types`, `test_narrative_validators`,
  `test_narrator_scrna`, `test_narrator_bulk`,
  `test_narrative_render_blocks`, `test_design_matrix_validator`) -> 23 passed
- Integrity suite (`test_registry_integrity`, `test_doctor`,
  `test_chromatin_dispatch_gate`) -> 8 passed
- `python -m pytest -q tests/test_pytest_smoke.py` -> 86 passed / 4 skipped
- `git diff --check` -> pass

## X5 typed IPC contract execution log (2026-05-28)

Implemented after `e616d30`:

- Added pydantic-backed `aria/utils/script_contracts.py` with `ScriptContract`,
  `ContractField`, `ContractIssue`, contract version `1.0`, required input
  fields, required file checks, success-output checks, and validation levels.
- `EnvironmentManager.run_in_stack` now validates registered script inputs
  before writing IPC JSON or invoking conda, validates output JSON after script
  completion, attaches `ipc_contract` metadata to valid outputs, and returns
  structured contract errors.
- Registered contracts for critical RNA and chromatin scripts:
  `rna_qc`, `rna_clustering`, `rna_bulk_de`, `rna_diff_abundance`,
  `rna_pseudobulk_de`, `rna_pathway_per_cluster`, `rna_cellcomm`,
  `rna_trajectory`, `chromatin_qc`, and `chromatin_peaks`.
- `registry_integrity` now checks IPC contract script existence and version
  coherence.
- `setup.py` and `install.sh` now include `pydantic>=2.0.0`.

Validation:

- `python -m compileall -q aria` -> pass
- `python -m pytest -q tests/test_script_contracts.py tests/test_registry_integrity.py tests/test_doctor.py` -> 10 passed
- `python -m pytest -q tests/test_script_contracts.py tests/test_registry_integrity.py tests/test_doctor.py tests/test_chromatin_dispatch_gate.py tests/test_design_matrix_validator.py tests/test_narrator_scrna.py` -> 20 passed
- `python -m pytest -q tests/test_pytest_smoke.py` -> 86 passed / 4 skipped
- `aria doctor --smoke` via TUI -> pass with the expected tracked HiC scaffold warning
- `git diff --check` -> pass
- Legacy `python tests/test_environment_manager.py` was also attempted and
  failed for pre-existing environmental reasons: scanpy/numba cache locator in
  `aria-env` and leftover temp JSON files in `~/.aria/workspace`; no X5
  contract failure was indicated by that legacy script.

### HiC scaffold-dispatch — DECIDED 2026-05-28 (accept as tracked warning)

`MODALITY_VALIDATION["HiC"]` is `level="scaffold"` with `dispatch_enabled=True`.
Reviewed and **accepted as a tracked warning**, not silenced: unlike the
chromatin family (whose `chromatin_motifs.py` / `chromatin_differential.py` are
absent, hence hard-blocked), all HiC scripts exist and run
(`hic_inspect`, `hic_qc_and_balance`, `hic_topology`), so prior QC behavior is
preserved. `registry_integrity.check_registry_integrity` emits a
`scaffold_dispatch_enabled` **warning** so the state stays visible and auditable.
Caveat carried forward: `hic_topology` emits TADs/loops, which are unvalidated
biological conclusions under a scaffold level — the proper resolution
(inspect+QC only, gating topology until HiC validation) is deferred to the HiC
validation milestone. Recorded as ADR-012 in `DECISIONS.md`.

Status: X1, X2, X3, X4, X17 are CLOSED in this v4.5.3 work (committed, not
tagged — no `v4.5.3` tag was created). Remaining external-audit items
(X5–X16, X18–X20) stay open in the plan above.
