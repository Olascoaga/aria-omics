---
status: active
source_of_truth_for: senior_audit_2026_05_29
last_updated: 2026-05-29
supersedes: []
related:
  - audit/2026-05-28_senior_audit.md
---

# ARIA Senior Audit (2026-05-29) — full code/architecture/docs review

Audit performed at HEAD `482ad79` (`Document v4.5.4 baseline`; `v4.5.4` tag at
`b7cd67f`, one commit behind HEAD) from a dual lens: senior software engineer +
senior computational-biology/bioinformatics researcher. Every finding below is
grounded in reading the implementation, with `file:line` evidence.

This file is a remediation tracker in the same style as
`audit/2026-05-28_senior_audit.md`. It does **not** supersede that file; it
records findings that survived the pre-ATAC freeze plus newly discovered ones.
Items marked **[NEW]** are not in the 2026-05-28 tracker.

Status legend: ☐ planned · ◐ in progress · ✅ done.

---

## Verdict

In its **validated path (scRNA / bulk RNA → report)** ARIA is solid and honest:
provenance, dual FDR, power reconciled to the decision rule, claim compiler,
narrative validators. The risk is **uneven maturity**: the "brain"
infrastructure (ParameterAdvisor, DebateCouncil, IntegrationAgent) is written in
production prose but internally contains **checkpoints that do not block,
fabricated metrics, and uncontrolled non-determinism**. Several of those cracks
directly contradict the five inviolable principles ("LLM proposes / code
guarantees / user decides", "no fabricated scientific outputs", "auditable
decisions").

---

## Findings — what is WRONG (bugs / principle violations)

### B1 [HIGH] [NEW] — ✅ fixed: in-agent parameter checkpoints do not block; they are decorative

`scrna_agent.py:812-840` calls `advise_leiden_resolution`, then
`publish_escalation(checkpoint=3, options=["Use recommended","Enter custom
resolution","Skip clustering"])`, then **immediately** runs
`run_in_stack(... resolution=decision.chosen_value ...)`. There is no
`Event`/wait. The user's response (custom / skip) can never take effect. Same
fire-and-forget pattern in `integration_agent.py:274-301` (WNN k) and the MOFA
path. CP1/CP2/CP3-thresholds/CP3.5 **are** awaited because they happen in the
orchestrator before dispatch; the CP3s published from inside the dispatch
thread are not. This violates "LLM proposes / **code guarantees** / **user
decides**": the user believes they choose the Leiden resolution, but the
advisor value always runs.

Evidence: `aria/agents/scrna_agent.py:812-845`,
`aria/agents/integration_agent.py:259-301`,
`aria/agents/orchestrator_agent.py:218-304` (async resolution model).

Status: ✅ fixed in the post-`v4.5.4` B1 remediation. `MessageBus` now exposes a
condition-backed `wait_for_checkpoint_resolution`; `BaseAgent` exposes
`publish_blocking_escalation` for in-dispatch parameter checkpoints; `scRNAAgent`
blocks Leiden resolution until user approval/custom override/skip; and
`IntegrationAgent` blocks WNN k and MOFA+ factor count before script execution.
Agent-parameter CP3 messages carry `agent_parameter_checkpoint=True`, so
`OrchestratorAgent.on_checkpoint_resolved` resolves them without invoking the
threshold CP3 dispatch path. The TUI live loop no longer marks ESCALATION
messages as seen before the checkpoint handler can display them. Validation:
`aria-env` `compileall` pass; B1 regressions for blocking/custom Leiden,
WNN skip-before-script, and internal CP3 no-dispatch passed; IntegrationAgent
legacy suite 24/24 passed; full `tests/test_pytest_smoke.py` 90 passed /
4 skipped.

### B2 [HIGH] [NEW] — Fabricated metrics shown to the user as if measured

`parameter_advisor.py:158-192` `MetricEvaluator.wnn_k` returns RNA/ATAC weights
from a linear heuristic (`0.6 - k*0.003`), flagged `is_estimated: True`, with a
comment "*real WNN evaluation planned for v0.3*". But
`integration_agent._format_wnn_checkpoint` (`:718-731`) prints them as
`RNA_weight=0.54 ATAC_weight=0.46` without disclosing they are invented. This is
fabricated scientific output presented to the user — a direct breach of "no
fabricated scientific outputs". Latent today (integration is outside the
validated path) but a landmine for v4.8.

### B3 [HIGH] [NEW] — `modularity` is always 0.0 → 30% of the clustering score is dead weight

`parameter_advisor.py:115-119` has `pass  # full implementation in production`
inside the igraph try block, so modularity is never computed. The score is
`sil*0.6 + mod*0.3` (`:568`) with `mod` always 0. The resolution recommendation
is de facto silhouette-only but is presented as modularity-informed. Functional
bug + honesty issue.

### B4 [MEDIUM] [NEW] — `IntegrationAgent` is not gated by `MODALITY_VALIDATION`

X4 gated chromatin/HiC scaffolds at the orchestrator
(`orchestrator_agent.py:117-159`), but `integration_agent` is dispatched on
`if n_mods >= 2 or plan.get("integration_needed")`
(`orchestrator_agent.py:650-664`) with no validation-level gate. The scripts
`integration_{wnn,mofa,peak2gene}.py` **exist** (not stubs) and are unvalidated
v4.8 roadmap. If the user has `aria-integration-env`, this runs unvalidated
science and **publishes findings with `Confidence.HIGH/MEDIUM`**
(`integration_agent.py:354,499`). The pre-ATAC integrity gate does not cover
integration.

### B5 [MEDIUM] [NEW] — Causal guard over-fires and degrades legitimate findings

`validators.py:15-54` `CAUSAL_PATTERNS` includes `"regulates"`, `"controls"`,
`"promotes"`, `"inhibits"`, `"activates"`, `"mediates"`. The guard also scans
`evidence.label` (`:129-135`), and GO/ORA term names literally contain
"regulation of…", "positive regulation of…". Result: purely descriptive ORA
blocks get a causal caveat and a **one-level confidence downgrade** (`:150`).
This is a systematic false positive that pollutes otherwise-honest reports. The
guard must distinguish "term name in evidence" from "author's asserted claim".

### B6 [MEDIUM] [NEW] — The canonical operational memory is in `.gitignore`

`.gitignore` excludes `memory/`, `CLAUDE.md`, `AGENTS.md`. That means
`PROJECT_STATE.md`, the roadmap plans, and `memory/audit/*` — the declared
source of truth — are **not version-controlled**: no history, no diff, no
recovery if the Medusa machine is lost. For a project whose entire discipline is
"auditable decisions", the decision log is the one thing without version audit.
Likely intentional for privacy, but then it needs a private mirror repo.

### B7 [MINOR] [NEW] — Housekeeping

- `codigo_aria.txt` (1.6 MB code dump) is **untracked and not ignored** → it is
  the cause of `git_dirty=True` in provenance reports.
- `_bh_correct` is duplicated in `rna_pseudobulk_de.py:72` and
  `rna_diff_abundance.py:347`.
- `env_manager` (`environment_manager.py:509`) and `bus`
  (`message_bus.py:165`) are instantiated at module import → import-time side
  effect (creates `~/.aria/workspace`).

### B8 [HIGH] [NEW-2] — ✅ fixed: v4.5.4 power disclosure is stale for the per-cluster default

(Surfaced by the second AI audit; verified.) `rna_pseudobulk_de.py:692-704` emits
a `power` block whose `note` states unconditionally: *"Significance is declared
with global BH-FDR across all blocks; the empirical per-test cutoff is
effective_alpha_global."* But the v4.5.4 default is `fdr_strategy="per_cluster"`
(`:146`), under which significance is `padj_local` per block — **not** global BH.
`power_estimate_at_effective_alpha` is also computed at the global effective
alpha (`:616-624`) regardless of strategy. So under the shipped default the power
disclosure and effective-alpha narrative describe a decision rule that is not the
one applied. A reviewer reading Methods would be misled. The `note`/effective-α
must branch on `fdr_strategy`. This refines R1 with a concrete, current bug.
Evidence: `aria/scripts/rna_pseudobulk_de.py:146,616-624,692-704`,
`memory/DECISIONS.md:319` (ADR-015).

Status: ✅ fixed in the post-`v4.5.4` B8 remediation. The script now records
per-block `effective_alpha_local`, `effective_alpha_global`,
`effective_alpha_primary`, and `effective_alpha_strategy`; computes
`power_estimate_at_effective_alpha` at the primary strategy-specific cutoff; and
branches the top-level `power.applied_threshold` / `power.note` via
`_power_disclosure_for_strategy`. For the default `per_cluster` strategy,
`effective_alpha_global` is explicitly labeled as a secondary whole-experiment
diagnostic rather than the applied decision rule. Validation: `compileall` pass;
`tests/test_pytest_smoke.py::test_global_fdr_is_more_conservative_than_local_family`,
`tests/test_pytest_smoke.py::test_pseudobulk_power_disclosure_tracks_primary_fdr_strategy`,
and `tests/test_pytest_smoke.py::test_power_estimate_monotone_in_n` all passed;
`aria-env` full `tests/test_pytest_smoke.py` passed 87 / skipped 4. Base Python
smoke failed from environment problems unrelated to B8 (`litellm` absent and
NumPy 2 ABI conflicts in compiled packages).

### B9 [HIGH] [NEW-2] — chromatin_qc.py ships placeholder metrics (fabricated if ungated)

(Second AI audit; verified.) `chromatin_qc.py` returns invented QC numbers:
`_compute_tss_enrichment` returns `base + hash(frag_file) % 30 / 10` — a function
of **file size and a string hash**, not TSS signal (`:320-335`); `_estimate_frip`
returns the constant `0.35` (`:400-404`); `_compute_fragment_sizes` sets
`n_barcodes = len(set())` ≡ 0 (`:368`). Real fragment-size/mito parsing exists
alongside, so it is a mix of real and fabricated. Today this is mitigated because
scATAC is `dispatch_enabled=False` (`orchestrator_agent.py:130`), but the code is
**not** ready for v4.6: unblocking it without rewriting these helpers would emit
fabricated science (breach of ADR-002). The v4.6 first step ("validate
chromatin_qc.py standalone") must replace these, not wrap them.

### B10 [HIGH] [NEW-2] — Bulk RNA silently rounds any matrix to integers (no raw-count guard)

(Second AI audit; verified.) `rna_bulk_de.py:_load_counts` ends with
`counts = counts.round().astype(int)` (`:694`) with **no** raw-vs-normalized
detection. TPM/CPM/FPKM/log-normalized inputs become pseudo-counts handed to
DESeq2 with no provenance and no caveat. This is the bulk analog of the pseudobulk
log-norm issue (R7), but **worse**: pseudobulk has recovery provenance +
`confidence=low` cap (v4.5.4), bulk has nothing. Needs a raw-count classifier and
a hard-refuse path mirroring `rna_pseudobulk_de`'s `allow_lognorm_recovery=False`.
Evidence: `aria/scripts/rna_bulk_de.py:632,694`.

### B11 [MEDIUM] [NEW-2] — Partial ADR-011 violations (hardcoded biological content in runtime)

(Second AI audit; verified, with nuance.) Confirmed hardcoded biology in runtime
paths:
- `design_intelligence.py:317-321` — cell-type alias map `{"microglia":
  {"Microglia"}, "opc": {"OPC"}, "astro"...}` (brain/hippocampus types — exactly
  the over-fitting class ADR-011 was written to prevent).
- `data_audit_agent.py:651-656` — keyword list contains `"ifn "`, `"ifn-"` even
  though the adjacent comment claims *"no 'interferon'"*. Comment contradicts code.
- `rna_bulk_de.py:2143-2150` — `_mock_pathways` fabricates "immune response" /
  "inflammatory response" GO terms when gseapy is unavailable (verify it is gated
  by `mocks_allowed()`; if not, it is silent fake science, ADR-002).

Nuance / partial REJECT: `data_audit_agent.py:443` `human_markers` (a fixed gene
set) is defensible — it is used for **species inference**, a technical detection,
not a biological claim. Recommend keeping it but documenting it as an explicit
ADR-011 exception (like the PBMC golden fixture), rather than removing it.

### B12 [MINOR-MEDIUM] [NEW-2] — SetupAgent ↔ EnvironmentManager contradiction + orphaned alias path

(Second AI audit; verified + extended.) `setup_agent.py:7` states the philosophy
*"No detection of what the user already has. No aliases. No clever matching."*
But `environment_manager._resolve_env` reads `~/.aria/env_aliases.json` first and
resolves stacks through it (`:450-476`). Worse: a grep finds **no writer** for
`env_aliases.json` anywhere in `aria/` — the alias-resolution branch is effectively
dead unless a user hand-creates the file. Either wire SetupAgent to write it
(and fix the docstring) or remove the branch.

---

## Findings — what can FAIL (runtime reliability)

### R1 [HIGH] [NEW] — Narrative reproducibility is not guaranteed

No `temperature` or `seed` on any LLM call (`provider.py:249-265`, confirmed by
grep across `aria/llm/*.py`). `--reproducible` covers code+data+locks, but the
prose, parameter justifications, and **DebateCouncil verdicts/confidence** are
stochastic. The prompt cache masks this *within* a machine; two clean runs (or
two reviewers) can get different confidence levels for the same evidence. For a
publication-grade tool the reported confidence level should be deterministic.
Cheap fix: `temperature=0` + fixed `seed`, both recorded in provenance.

### R2 [MEDIUM] [NEW] — DebateCouncil is dormant where it matters

`DebateCouncil.resolve` is invoked **only** in `chromatin_agent`,
`genome_arch_agent`, and `integration_agent` (grep) — the three
scaffold/non-dispatchable modalities. The validated path (scRNA / bulk RNA) has
no debate. The flagship anti-sycophancy mechanism does not protect the only
science ARIA ships today. (claim_compiler + causal guard cover part of it, but
the "alternative hypothesis first" reasoning never runs on validated results.)

### R3 [MEDIUM] [NEW] — LLM calls have no timeout

`completion(**kwargs)` (`provider.py:265`) passes no `timeout`. A hung provider
blocks the dispatch thread indefinitely; there is no watchdog. Tier fallback
only triggers on exception, not on a hang.

### R4 [MEDIUM] [NEW] — Silent model degradation + stale default model IDs

A transient 429 on the HEAVY model silently drops to sonnet→haiku
(`provider.py:160-170`) without the report recording "this section ran on a
fallback model". Reasoning quality drops with no trace. Defaults are also stale
(`claude-opus-4-7`, `gemini/gemini-1.5-pro` at `provider.py:73,77`) → the first
HEAVY call may always fail and run on fallback by default.

### R5 [MEDIUM] [NEW] — Orphaned grandchild processes on timeout

`environment_manager.py:219-281`: `subprocess.run(timeout=...)` kills the direct
child (`conda run`) but grandchildren (the real `python` + BLAS/numba threads)
can survive on large datasets. No `start_new_session` / process-group kill.

### R6 [MEDIUM] — Crash loses bus state (F-ENG-BUS, still open)

`message_bus.py:91` is an in-memory `deque(maxlen=100k)`; `memory.py` does not
persist the bus log. A mid-run crash loses findings/decisions of completed
stages (resume only covers file-valid h5ad stages). Also
`get_findings`/`get_pending_checkpoints` do an O(n) scan of the full deque on
every 0.5s TUI poll tick — the bus part of F-ENG-PERF is still unfixed (only the
env-detection part was memoized).

**[NEW-2] cross-run contamination via the global bus.** `headless.py:184-194`
deliberately polls the **full** bus log unfiltered by `experiment_id` (because
some STATUS messages carry no `experiment_id`). The global singleton `bus`
(`message_bus.py:165`) means two concurrent headless runs in one process would
read each other's events. Single-run headless is fine and the choice is
documented, but concurrency is unsafe. A per-run bus instance (or mandatory
`experiment_id` on every message) closes this.

### R7 [MINOR-MEDIUM] — Fragile log-norm recovery probe

`rna_pseudobulk_de.py:162-189`: the integer-likeness probe uses `mat[:200]`
(first 200 cells, **not random**); if the h5ad is ordered by cell type or
condition the probe is biased. Assumes `log1p` and `scale=1e4`; another
normalization yields plausible-but-wrong "recovered" counts. Well mitigated by
the v4.5.4 cap to `confidence=low`, but the sampling should be randomized.

---

## Findings — what is MISSING

### Scientific (comp-bio lens)

- **C1 [MEDIUM] [NEW] — Differential abundance uses Poisson without
  overdispersion.** `rna_diff_abundance.py:240-246`: Poisson GLM with offset
  assumes mean=variance; cell-type counts across replicates are overdispersed →
  **anti-conservative (over-calls significance)**. Because this gates the
  composition covariate in pseudobulk DE (`scrna_agent` turns on
  `composition_covariate` when `any_significant`), a false positive here turns on
  spurious correction there. Field standard: NB / quasi-Poisson, or `propeller`
  (logit + empirical Bayes, Phipson 2022).
- **C2 [MEDIUM] [NEW] — ORA background is global, not per-cluster.**
  `rna_pseudobulk_de.py:200-205`: `background_genes` = genes expressed in *any*
  cell of the whole dataset. Per-cluster ORA should use the genes
  tested/expressed in that cell type's pseudobulk as the universe. A global
  background inflates per-cluster enrichment.
- **C3 [LOW-MEDIUM] [NEW] — Composition-covariate statistical tension.**
  `rna_pseudobulk_de.py:340-348`: the covariate is `log(fraction of the cell
  type itself)`. When abundance shifts with condition (exactly when it is turned
  on), the covariate is **collinear with the condition factor** → variance
  inflation / partial absorption of the real signal. The 2026-05-28 audit closed
  this as "minor" after a probe; from the science lens, regressing a cell type's
  own proportion in its own DE is methodologically questionable (muscat/Crowell
  handle composition via normalization / separate abundance modeling, not a
  self-proportion covariate).
- C4 — No effect-size shrinkage (apeglm/ashr), no s-values; raw LFCs. No ambient
  RNA correction (SoupX/decontX), only Scrublet. No integration QA (LISI/kBET) —
  `integration_qc.py` exists but consumes only silhouettes. LIANA `n_perms=100`
  is still low (`rna_cellcomm.py`). No donor-level QC, no mixed-models /
  beta-binomial for composition (aligns with C1). (Carried from prior audit +
  second AI audit.)
- **C8 [HIGH] [NEW-2] — The v4.6 scATAC entry path does not match its planned
  input.** `V46_SCATAC_PLAN.md:11` names a `.h5mu` (MuData) validation input, but
  `DataAuditAgent._scan_directory` (`data_audit_agent.py:289-296`) does not list
  `.h5mu` in its extension set, and `chromatin_qc.py` requires a fragment file
  (`"fragments" in f.lower() or f.endswith(".tsv.gz")`, `:94-100`) — not MuData.
  The first step of v4.6 will not detect the input or will fall into an
  unvalidated path. v4.6 must add `.h5mu` detection + a real MuData reader before
  any chromatin QC runs.

### Engineering

- C5 — **God-files** unsplit: `narrative_agent.py` 2742, `rna_bulk_de.py` 2455,
  `_narrative_scrna.py` 2211, `scrna_agent.py` 2043. `tests/test_pytest_smoke.py`
  is 2919 LOC and hard to maintain.
- C6 — Privacy (X10, open): params with `data_path` to human data are written in
  cleartext to `~/.aria/workspace/input_*.json`, and prompts (with sample/gene
  names) are sent to cloud APIs. No air-gapped firewall + redactor. The LLM cache
  has no TTL / no version salt (`provider.py:294-300`).
- C7 — Test gaps: `setup_agent` (719 LOC), `rna_quantify`/`rna_align` (FASTQ
  path), `integration_*`, `genome_arch`. Much of the suite relies on mocks.

---

## Documentation / process findings

- **D1 [MEDIUM] [NEW-2] — Duplicate ADR number.** `memory/DECISIONS.md` has TWO
  `ADR-013`: line 233 "DESeq2 Designs Must Pass Preflight Matrix Validation" and
  line 286 "Claims Are Evidence-Tiered". Renumber one (e.g. the design-matrix ADR
  to a free number) and fix back-references in `PROJECT_STATE.md`/`START_HERE.md`.
- **D2 [LOW] — HEAD vs tag drift in memory.** `START_HERE.md`/`PROJECT_STATE.md`
  say "current HEAD is the `v4.5.4` tagged commit", but `git describe` is
  `v4.5.4-1-g482ad79` (HEAD is one docs commit ahead of the tag `b7cd67f`).
  Correct the memory wording without moving the tag. (Both AI audits agree.)
- **D3 [MEDIUM] — Canonical memory is untracked.** See B6: `memory/` is in
  `.gitignore`, so DECISIONS/PROJECT_STATE/audits have no version history.

## Novel proposals (real added value)

1. **[P-CHK] Checkpoint contract that closes B1 at the root.** Make
   `publish_escalation` inside the dispatch thread return a `Future`/`Event` the
   agent **must** await (`bus.await_resolution(msg_id, timeout)`), or move every
   parameter decision to pre-dispatch (as thresholds already are). Add an E2E
   test that fails if an agent runs a script after an unresolved checkpoint. Turns
   decorative checkpoints into real guarantees — the literal embodiment of "code
   guarantees".

2. **[P-LEDGER] Planned-vs-Run ledger.** A deterministic per-run manifest
   `{planned, dispatched, completed, skipped+reason}` that the NarrativeAgent
   **must** render. Any plan↔execution divergence becomes visible (would have
   caught the PBMC thin-report bug before the rerun). Closes the dispatch
   integrity gap noted in the 2026-05-28 audit.

3. **[P-DET] Deterministic LLM + model provenance in the report.**
   `temperature=0`, fixed `seed`, and record per section which model/tier actually
   answered (and whether it was a cache hit / fallback) in `methodology.json`.
   Makes the *narrative* reproducible, not just the compute. Closes R1 + R4.

4. **[P-CLAIM2] Quantitative stats-evidence gate in the Claim Compiler.** Today
   it classifies by analysis type; add thresholds: an "associative" claim
   requires `n_sig>0 ∧ power_at_effective_alpha>X ∧ ¬low_power ∧
   ¬lognorm_recovered`. Auto-downgrade to descriptive/insufficient when the
   numbers don't hold. Closes the gap between "the design type licenses it" and
   "the numbers support it".

5. **[P-DEVIL] Cheap DebateCouncil on the validated path via the Claim
   Compiler.** Instead of 3 expensive LLM rounds, a deterministic single-shot
   devil's-advocate over each associative/strong-mechanistic claim: "what is the
   alternative technical explanation (batch, ambient, doublets, composition) and
   what evidence rules it out?". Reuses the kernel without round cost (also
   attacks X11). Closes R2.

6. **[P-MULTIVERSE] Lightweight robustness multiverse (X15) on two axes you
   already compute.** Re-run pseudobulk DE under `fdr_strategy ∈ {per_cluster,
   global}` and `composition_covariate ∈ {on, off}`, report "N genes stable
   across all 4 combinations". Nearly free (both FDR families are already
   computed) and yields a publishable robustness metric.

---

## Reconciliation with the second AI static audit (2026-05-29b)

A second AI ran an independent static audit (no suite run, no file changes). Each
of its claims was checked against the real code here. Adjudication:

**ACCEPTED + verified (folded in above as `[NEW-2]`):**
- FDR/power inconsistency in v4.5.4 → **B8** (its #1; sharp, current, verified).
- v4.6 `.h5mu` path mismatch → **C8** (its #2; verified, real readiness gap).
- chromatin placeholder metrics → **B9** (its #3; verified — TSS via hash/size,
  FRiP constant, n_barcodes≡0).
- bulk `_load_counts` rounds non-raw matrices → **B10** (its #4; verified).
- partial ADR-011 violations → **B11** (its #5; verified, with the
  `human_markers`→species-inference nuance partially rejected).
- headless reads the global bus unfiltered → folded into **R6**.
- SetupAgent ↔ EnvironmentManager contradiction → **B12** (verified + extended:
  no writer for `env_aliases.json` exists in `aria/`).
- duplicate ADR-013 → **D1**; HEAD-vs-tag drift → **D2** (both AIs agree).

**ACCEPTED but already tracked (here or in the 2026-05-28 file):**
- bus not durable → R6 / F-ENG-BUS. IPC contracts shallow (no FASTQ/HiC/
  integration) → X5 partial. DataAudit no scan limits (`rglob("*")`,
  `data_audit_agent.py:288`) → X13 (its evidence confirms it). Legacy tests can
  false-green → P2 test debt. HiC TAD/loop gating → ADR-012. Missing science
  (ambient RNA, LISI/kBET, shrinkage, s-values, beta-binomial composition) →
  C1/C4. Planned-but-not-run manifesto → **P-LEDGER** (we converged
  independently). Multiverse → **P-MULTIVERSE** (converged). ChromatinNarrator
  from day one of v4.6 → already mandated by `V46_SCATAC_PLAN.md`.

**Where the two audits converge (high confidence, fix first):** non-blocking
checkpoints (B1), LLM non-determinism (R1), the v4.5.4 FDR/power text (B8),
planned-vs-run ledger, multiverse, privacy/air-gapped mode.

**Net new from the second audit:** B8, B9, B10, B11, C8, D1 — all verified, all
concrete. Its strongest unique catch is **B8** (the v4.5.4 power text is wrong for
the shipped default — a live scientific-honesty regression introduced by the same
commit that was meant to improve honesty).

**Its proposals mapped to this file:** Run Evidence Graph → X16; Planned-vs-
Executed Ledger → P-LEDGER; Raw Matrix Classifier → unifies B10 + R7 (new
**P-RAWCLASS** below); Scientific Calibration Suite → extends X6; Privacy/
Air-Gapped → X10 / C6; Multiverse → P-MULTIVERSE; ChromatinNarrator → v4.6 plan.

7. **[P-RAWCLASS] Shared raw-count classifier (closes B10 + R7).** One detector
   (`is_raw_counts` / `classify_matrix → {raw, cpm, tpm, lognorm, scaled}`) used by
   both `rna_bulk_de._load_counts` and `rna_pseudobulk_de`, with a single
   hard-refuse + provenance contract (`count_source`, confidence cap). Removes the
   duplicated, divergent count handling across the two DE scripts.

---

## Remediation plan — suggested order

| Priority | Items | Why |
|---|---|---|
| **Before v4.6** | B2/B3 (or scaffold-gate IntegrationAgent B4), R1 | Inviolable principles; B1 and B8 are fixed. |
| **v4.5.5 honesty patch** | B5, B10/P-RAWCLASS, B11, C1, C2, R7, D1/D2, B7 | Fix science/precision/docs without a new modality. |
| **v4.6 readiness gate** | C8, B9, B12 | Must precede chromatin work: `.h5mu` detection, real QC (no placeholders), env-alias path. |
| **During/after v4.6** | R3-R6, C4-C7, P-LEDGER/P-DET/P-CLAIM2/P-DEVIL/P-MULTIVERSE, D3 | More effort or depend on installing envs. |

Cheapest highest-value remaining: **R1** (`temperature=0` + seed) plus
**B2/B3** or **B4**. B1 (non-blocking checkpoints) and B8 (v4.5.4 power/FDR text)
are fixed. These items are direct
contradictions of the project's identity and fix in a few lines.
**Do not start v4.6** until C8 + B9 are addressed — the planned `.h5mu` input is
undetectable today and chromatin QC would emit fabricated metrics.

---

## Notes on prior-audit items still open at this HEAD

- F-ENG-BUS (persistence) → R6 here, still open.
- F-ENG-PERF → env part memoized; **bus-scan part still open** (R6).
- X10 privacy, X11 debate cost, X12 ParameterAdvisor bias, X13 DataAudit scan
  limits, X15/X16/X18-X20 → still open (X12 confirmed real: `_choose_best`
  historical bonus at `parameter_advisor.py:660-673` can override objective
  score; `_recall_similar_decisions:610-637` leaks parameter choices across
  unrelated datasets/experiments).
- God-file split, apeglm shrinkage, IHW/hierarchical FDR, FASTQ-path tests →
  still open (C4/C5).
