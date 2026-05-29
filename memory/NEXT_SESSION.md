---
status: active
source_of_truth_for: next_session
last_updated: 2026-05-29
supersedes:
  - archive/sessions/2026-05-14_close.md
---

# NEXT SESSION

## ACTIVE MANDATE (2026-05-29) — clear the senior audit BEFORE v4.6

Samael's directive: **resolve the entire `memory/audit/2026-05-29_senior_audit.md`
tracker before starting v4.6 scATAC.** That file is the authoritative source of
truth for these items (severity, `file:line` evidence, fixes). Do not start v4.6
implementation until every item below is closed or explicitly deferred with
Samael's sign-off.

This audit reconciles a second independent AI static audit; `[NEW-2]` items came
from it and were code-verified. Suggested execution order (full detail in the
audit file):

**Stage 1 — Before anything else (cheap, principle-level + live regression):**
- ✅ **B1** — in-agent parameter checkpoints do not block (Leiden res, WNN k,
  MOFA). Make them await resolution or move to pre-dispatch. E2E test.
- ✅ **B8** — v4.5.4 `power` block text says "global BH-FDR" but default is
  per-cluster → Methods is wrong for the shipped default. Branch on `fdr_strategy`.
- ☐ **R1** — set `temperature=0` + fixed `seed` on all LLM calls; record in
  provenance (narrative/confidence are non-deterministic today).
- ☐ **B2/B3** — ParameterAdvisor fabricated WNN metrics shown as measured;
  `modularity` always 0.0. Fix or stop displaying as real.
- ☐ **B4** — gate `IntegrationAgent` by `MODALITY_VALIDATION` (scaffold) like
  chromatin/HiC.

**Stage 2 — v4.5.5 honesty/precision patch (no new modality):**
- ☐ **B5** — causal guard over-fires on GO/ORA term names; downgrades honest
  blocks. Distinguish "term name in evidence" from asserted claim.
- ☐ **B10 + R7 → P-RAWCLASS** — shared raw-count classifier + hard-refuse for
  bulk (`_load_counts` rounds any matrix to int) and pseudobulk; random-sample the
  log-norm probe.
- ☐ **B11** — ADR-011: remove/justify hardcoded biology (microglia/OPC aliases,
  `ifn` token vs comment, `_mock_pathways`); document `human_markers` as an
  explicit species-inference exception.
- ☐ **C1** — differential abundance Poisson → NB/quasi-Poisson or propeller
  (overdispersion; it gates the composition covariate).
- ☐ **C2** — ORA background should be per-cluster, not global.
- ☐ **D1** — duplicate `ADR-013` in DECISIONS.md (renumber one).
- ☐ **D2** — memory says HEAD = tagged v4.5.4 but git is `v4.5.4-1-g482ad79`;
  correct wording WITHOUT moving the tag.
- ☐ **B7** — housekeeping: gitignore/remove `codigo_aria.txt`, dedupe
  `_bh_correct`, import-time side effects.

**Stage 3 — v4.6 readiness gate (MUST precede chromatin work):**
- ☐ **C8** — DataAudit does not detect `.h5mu`; `chromatin_qc` needs fragments.
  Add `.h5mu` detection + real MuData reader.
- ☐ **B9** — `chromatin_qc.py` placeholder metrics (TSS via hash/size, FRiP
  constant 0.35, `n_barcodes`≡0). Replace with real QC, don't wrap.
- ☐ **B12** — SetupAgent ("No aliases") vs EnvironmentManager (reads
  `env_aliases.json`, no writer exists). Wire it or remove the branch.

**Stage 4 — during/after, larger or env-dependent:**
- ☐ **R3** LLM call timeout · **R4** model-degradation provenance + refresh stale
  default model IDs · **R5** process-group kill on subprocess timeout · **R6**
  bus durability + per-run bus (headless reads global bus unfiltered) ·
  **C4** ambient RNA / LISI-kBET / shrinkage / s-values · **C5** god-file splits ·
  **C6/X10** privacy / air-gapped mode + LLM cache TTL · **C7** test gaps ·
  **D3/B6** take `memory/` out of `.gitignore` (or private mirror) so DECISIONS/
  audits get version history.
- Proposals to land alongside: **P-CHK**, **P-LEDGER** (planned-vs-run manifest),
  **P-DET**, **P-CLAIM2**, **P-DEVIL**, **P-MULTIVERSE**, **P-RAWCLASS**.

**Start by:** open `memory/audit/2026-05-29_senior_audit.md`, continue Stage 1
with **R1** or **B2/B3** (highest remaining value), run the validation gates
(`compileall` + `tests/test_pytest_smoke.py` + targeted tests), and tick the box
here as each item closes. Mark `[NEW-2]` items as verified-fixed in the audit
file's status legend (☐ → ✅).

## Last Completed (most recent first)

- **B1 audit remediation (2026-05-29)** — internal parameter checkpoints now
  block dispatch-thread execution until the user/headless runner resolves them.
  Leiden resolution, WNN k, and MOFA+ factor count use blocking checkpoints with
  recommended/custom/skip handling; internal CP3 resolution no longer triggers
  threshold CP3 redispatch; the TUI live loop leaves ESCALATION messages visible
  for the checkpoint handler. Validation: `aria-env` `compileall` pass; B1
  regressions 3 passed; IntegrationAgent legacy 24/24 passed; full
  `tests/test_pytest_smoke.py` 90 passed / 4 skipped.
- **B8 audit remediation (2026-05-29)** — pseudobulk power disclosure now follows
  the actual `fdr_strategy`. The default per-cluster path uses per-block
  `effective_alpha_primary` for `power_estimate_at_effective_alpha` and labels
  `effective_alpha_global` as a secondary whole-experiment diagnostic; the global
  strategy still uses `effective_alpha_global` as primary. Validation:
  `aria-env` `compileall` pass; `tests/test_pytest_smoke.py` 87 passed /
  4 skipped. The base Python smoke remains unsuitable here because it lacks
  `litellm` and has NumPy 2 ABI conflicts in compiled scientific packages.
- **`v4.5.4` Scientific-Honesty Hardening (2026-05-29)** — tagged + pushed on
  top of `v4.5.3`. Per-cluster FDR is now the pseudobulk default
  (`fdr_strategy`, ADR-015); power reconciled with the global-BH decision rule
  (`effective_alpha_global` + `power_estimate_at_effective_alpha`, ADR-016);
  log-norm-recovered DE capped at `confidence="low"` (ADR-016). Validation:
  aria-env 128 passed / 4 skipped, aria-rna-env pydeseq2-gated 3 passed
  (recall=1.0/FDR=0 preserved), diff-check clean. **Behavior change: pseudobulk
  significant-gene counts differ from pre-v4.5.4 runs.** Details in
  `PROJECT_STATE.md` and `memory/audit/2026-05-28_senior_audit.md`.
- **`v4.5.3` tagged + pushed** at `bab6fbd` — fixes F-NEW-3 (the pre-ATAC
  integrity freeze X1-X9/X14/X17 had no tag).
- Senior dual-lens audit ran 2026-05-29; remaining open items (F-SCI-FDR
  hierarchical/IHW upgrade, F-ENG-BUS persistence, F-NEW-1 single-env test
  story + real Docker, F-NEW-2 composition-covariate justification, HiC
  topology gating) are documented in the audit file. v4.6 scATAC remains next.

## Last Completed

- 4.3 line closed.
- `v4.4` Publication Readiness is closed and tagged at `cbcde8e`
  (`v4.4 release: align version and tool provenance`).
- `v4.5` Raw Ingestion is closed at `1d54cc0`
  (`v4.5 raw ingestion bridge`), with patch tag `v4.5.1` at `a0b33dd`
  (`v4.5.1 add gated kb ingestion execution`).
- Pre-memory-reorganization docs commit: `2e0eb39`.
- Final code commit: `d3de169`.
- `v4.3.12` tag remains at `3a0c40e`.
- `v4.3.12.post1` tag remains at `805e0b2`.
- Existing tags must not be moved.
- Current HEAD is post-`v4.5.4` audit remediation on `origin/main`;
  `aria.__version__` remains `4.5.4`. Tag `v4.5.4` remains at `b7cd67f` and
  must not be moved. Verify the exact hash with
  `git log --oneline --decorate -5`.
- **PBMC Stage C blocker rerun reviewed.** Four blockers
  documented in `memory/roadmap/V44_PBMC_BLOCKERS.md`:
  1. ✅ CLOSED in `ba4e21e`: pseudobulk uses `replicate × condition` for
     paired designs and adds donor as a covariate for balanced paired DESeq2
     blocks.
  2. ✅ CLOSED in `ba4e21e`: scRNA intermediates are written to ARIA
     workspace paths and DataAudit ignores ARIA-generated contamination.
  3. ✅ CLOSED in `7a73c22`: embedding figures fall back to t-SNE or any
     other 2D `X_*` embedding and label the figure honestly.
  4. ✅ CLOSED in `7a73c22`: `methodology.json` persists `llm_usage`, and
     usage collection includes a small timestamp grace window.
- The clean PBMC rerun at
  `/home/medusa/.aria/reports/aria_20260516_115010_monocytes_tcells_bcells_-87a/report.html`
  passes the seven-question audit for commit `7a73c22`: exact commit,
  dependency locks, composition correction, local/global BH, replicate/power,
  input SHA-256, and LLM usage are visible; paired pseudobulk now produces 11
  analyzable significant blocks.
- `cbcde8e` finalizes release metadata: ARIA/setup/TUI now report `4.4.0`,
  and `methodology.json["tools"]` falls back to committed lockfiles for
  packages missing from the report writer's active Python env.
- Final PBMC Stage C report reviewed:
  `/home/medusa/.aria/reports/aria_20260520_124352_monocytes_tcells_bcells_-d72/report.html`.
  It was generated from `cbcde8e`, shows `aria_version=4.4.0`, embeds
  conda/pip lockfiles, reports lockfile-derived `pydeseq2=0.5.4` and
  `gseapy=1.1.13`, and passes the seven-question peer-reviewer audit:
  exact commit, exact dependencies, composition correction, local/global BH,
  min replicates/power, input SHA-256, and LLM token/cost accounting.
  `git_dirty=True` is attributable to unrelated untracked `codigo_aria.txt`;
  tracked files had no diff.
- Last stable tag: `v4.5.4`
  (`v4.5.4 scientific-honesty hardening`; per-cluster pseudobulk FDR default).
- Previous stable tag: `v4.5.3` (`bab6fbd`,
  pre-ATAC integrity freeze).
- Earlier stable tag: `v4.5.2`
  (`v4.5.2 narrative kernel`).
- Previous stable tag before that: `v4.5.1` (`a0b33dd`,
  `v4.5.1 add gated kb ingestion execution`).
- Previous base tag: `v4.5` (`1d54cc0`,
  `v4.5 raw ingestion bridge`).
- P0 audit fixes landed in `05f6f4e`
  (`Close P0 audit findings for v4.3.19`).
- P1 audit fixes landed in `3054318` on top of `v4.5.1` (untagged) —
  see "Audit Follow-ups After v4.5.1" below for the full list.
- **`v4.5.4` is the current pre-v4.6 reporting baseline.** It includes the
  `v4.5.2` Narrative Kernel plus `v4.5.3` integrity freeze and `v4.5.4`
  scientific-honesty defaults. scRNA and bulk RNA report sections are composed
  from validated `NarrativeBlock` objects via modality narrators.
- Last previous tagged release: `v4.3.18` (`9dc48aa`,
  `Make design intelligence choices explicit`).
- Maintenance tags after `v4.3.12.post1`:
  - `v4.3.13`: initial TUI long-prompt / covariate display patch.
  - `v4.3.14`: explicit multiline TUI prompt with `END` sentinel and h5ad
    human-species inference.
  - `v4.3.15`: focused scRNA h5ad pre-QC subsetting for requested cell groups.
  - `v4.3.16`: tightened scRNA focus, LIANA p-value display, integrated
    interpretation prompt sanitization, and richer scRNA decision logging.
  - `v4.3.17`: cross-modality Design Intelligence layer before compute.
  - `v4.3.18`: CP2 explicitly chooses recommended-only vs
    recommended+optional supported analyses, and cell-focus detection only
    uses explicit focus/restriction clauses.
  - `v4.3.19` (`eece228`; P0 code in `05f6f4e`): four P0 audit fixes.
    Pseudobulk DE narrative header derives label from `condition_col` instead
    of hardcoded "Age-associated"; chromatin QC fallback gated by
    `mocks_allowed()`; `rna_cellcomm` no longer swallows non-ImportError
    LIANA failures; `hic_inspect` reads `RAM_ESTIMATES_GB` from
    `aria.utils.hic_constants` (removes script→agent import). Release docs and
    version metadata are aligned to `4.3.19`.

## Audit Follow-ups After v4.5.1

The P1 items from the 2026-05-15 audit landed in commit `3054318`
(`Post-v4.5.1 audit hardening: close P1 findings`) on `main`, untagged on
top of `v4.5.1`:

- ✅ **P1 — bus dispatch robustness**: `MessageBus.publish` logs receiver
  exceptions and continues fan-out.
- ✅ **P1 — resume validity**: RNA `.h5ad` resume paths now validate cached
  h5ads with `aria.utils.safe_h5ad.h5ad_is_readable`; RNA script readers use
  `read_h5ad` for explicit corrupt-file errors.
- ✅ **P1 — failed-run TOCTOU**: failed-run eviction uses `failed/.lock` with
  `fcntl.flock` and ignores the lockfile during FIFO deletion.
- ✅ **New audit finding — decisions checkpoint type**:
  `memory.decisions.checkpoint` is now `TEXT`; existing INTEGER tables migrate
  in place and `store_decision` normalizes checkpoint IDs to strings.
- ✅ **New audit finding — kb hash TOCTOU**: `execute_kb_count` converts
  index/t2g hashing failures into structured blockers instead of crashing.
- ⚠️ **P2 — tests partly hardened**: the smoke wrapper for the legacy bulk
  script now fails if the script prints nonzero internal failures. Full
  conversion of `tests/test_bulk_rna.py` / `tests/test_scrna.py` to native
  pytest plus golden DE fixtures remains open.

Validation at `3054318`:

- `python -m compileall -q aria` -> pass
- `python -m pytest -q tests/test_pytest_smoke.py` -> 85 passed, 4 skipped
- `python -c "import aria; print(aria.__version__)"` -> `4.5.1`

## Narrative Depth After v4.5.1

Implemented locally after the P1 hardening commit, not yet committed:

- `aria/agents/_narrative_scrna.py` adds deterministic result-local helpers
  for:
  - top pseudobulk DE blocks with global/local FDR, top up/down genes,
    composition correction, power estimates, ORA support, and caveats;
  - differential-abundance overlap with DE blocks;
  - pathway support per block;
  - LIANA top non-autocrine ligand-receptor candidates with rank metric;
  - PAGA/DPT trajectory depth and interpretation limits.
- `build_scrna_integrated_interpretation` now reuses those helpers, so the
  final synthesis is deeper but still grounded in structured outputs.
- New regression test:
  `tests/test_pytest_smoke.py::test_scrna_narrative_adds_per_result_depth`.

Validation:

- `python -m compileall -q aria` -> pass
- Targeted narrative tests -> 4 passed
- `python -m pytest -q tests/test_pytest_smoke.py` -> 86 passed, 4 skipped
- `git diff --check` -> pass

Follow-up recommendation: before starting chromatin narrative in v4.6, mirror
this contract in `_narrative_chromatin.py`: each peak/accessibility/motif
result should receive local evidence, caveats, and synthesis rather than only
aggregate counts.

## Narrative Kernel v4.5.2

Implemented for the v4.5.2 closeout:

- New `aria.agents.narrative` package with:
  - `types.py` (`EvidenceItem`, `Caveat`, `NarrativeBlock`);
  - `protocols.py` (`ModalityNarrator`);
  - `registry.py` (`NarrativeRegistry`);
  - `validators.py` (claim/evidence, failed-analysis, low-confidence,
    causal-language, PAGA/DPT, and file-reference rules);
  - `render_blocks.py` (HTML composer).
- `ScrnaNarrator` emits blocks for QC, marker-discovery errors, composition,
  pseudobulk DE, ORA, LIANA, and trajectory.
- `BulkRnaNarrator` emits blocks for QC, contrasts, pathways, and power.
- `NarrativeAgent` uses block rendering for modalities with blocks and keeps
  legacy fallback for chromatin, Hi-C, and integration until their narrators
  exist.
- `methodology.json` persists serialized `narrative_blocks`.
- Offline scRNA harness reports now persist input SHA-256 records via
  `rna_narrative_adapter.py`.
- Version metadata is aligned to `4.5.2`.

Validation recorded so far:

- `python -m compileall -q aria` -> pass
- Narrative kernel tests:
  `python -m pytest -q tests/test_narrative_types.py tests/test_narrative_validators.py tests/test_narrator_scrna.py tests/test_narrator_bulk.py tests/test_narrative_render_blocks.py`
  -> 16 passed
- `python -m pytest -q tests/test_pytest_smoke.py` -> 86 passed, 4 skipped

Before starting v4.6 implementation, verify whether the final PBMC report
rerun was completed for v4.5.2 and recorded in `PROJECT_STATE.md`.

## Anti-Hardcode Follow-up After v4.5.2

After review, the new narrative-kernel synthetic tests were neutralized so
they no longer use PBMC/interferon/cell-type/gene labels. This is now governed
by `ADR-011`: runtime and generic tests must not hardcode biological content;
real names are reserved for explicit golden-dataset regressions and historical
validation notes.

## Start By

```bash
git status --short
git log --oneline --decorate -5
git tag --list 'v4.5*'
python -c "import aria; print(aria.__version__)"
```

Then read:

1. `memory/PROJECT_STATE.md`
2. `memory/AGENT_PROTOCOL.md`
3. A task-specific file only if needed.

## Known Caveats

- 4.3 is closed; do not reopen it for new features.
- `Memoria/`, if present, is a temporary local export folder and should not be
  treated as canonical.
- `memory/archive/` is history only.
- Reports generated before `805e0b2` / `d3de169` may still show older wording;
  new reports should use the fixed templates.
- TUI prompts can now be pasted as multi-line text; verify in an actual terminal
  by ending the question with a line containing only `END`.
- Design Intelligence now appears in the analysis-plan checkpoint and logs a
  `2.pre` decision. It recommends, offers optional analyses, and flags
  unsupported/not-recommended analyses before compute.
- In CP2, option 1 runs only recommended analyses; option 2 adds optional
  supported analyses. Unsupported/not-recommended items remain blocked unless
  the plan is modified explicitly.
- scRNA cell focus should not be inferred from incidental or negative mentions;
  it should require explicit focus/restriction language.
- Network may be restricted; LiteLLM cost-map warnings are non-fatal.
- Matplotlib may use a temporary cache if the user config dir is not writable.

## v4.4 sprint closed

The 2026-05-15 expert audit changed the order. `v4.4` is no longer scATAC —
it became a **Publication Readiness** sprint that closed scientific and
reproducibility blockers before any new modality. On 2026-05-20, Raw
Ingestion was promoted to `v4.5` and scATAC moved to `v4.6`.

- `v4.4` — Publication Readiness, closed at `cbcde8e`. Full plan in
  `memory/roadmap/V44_PUBLICATION_READINESS.md`. Tier 1 (5 scientific) +
  Tier 2 (5 reproducibility) items. Done criterion: the seven-question
  peer-reviewer test passes on a real scRNA rerun.
- `v4.5` — Raw Ingestion (`memory/roadmap/V45_RAW_INGESTION_PLAN.md`),
  closed at `1d54cc0`.
- `v4.6` — scATAC E2E (`memory/roadmap/V46_SCATAC_PLAN.md`), next.
- `v4.7` — bulk ATAC E2E (`memory/roadmap/V47_BULK_ATAC_PLAN.md`).
- `v4.8` — multimodal integration (`memory/roadmap/V48_INTEGRATION_PLAN.md`).

The P1/P2 audit items listed above remain valid follow-up work for `v4.4.x`
or later patches; `v4.4` itself is scoped tightly to Tier 1 + Tier 2.

### Sprint progress snapshot

- ✅ **T1.3** (commit `e8fc846`) — default `min_replicates_per_condition` 3;
  low_power_warning at n=2 surfaced in pseudobulk + bulk DE; narrative
  caveats wired; DesignIntelligence downgrades DE to `optional` at n=2.
  4 new pytest cases. See V44 plan for details.
- ✅ **T1.1** (commit `52bf703`) — differential abundance + composition
  correction. New `rna_diff_abundance.py` (Poisson-offset GLM with
  Fisher fallback), `composition_covariate` param in pseudobulk DE,
  scrna_agent runs DA before DE and toggles the covariate on
  `any_significant`, new "Cell-type abundance" narrative section, DE
  header reports composition-corrected count. 4 new pytest cases.
- ✅ **T1.2** (`ab246ca`, wording tightened in `efa7136`) — global FDR across all
  cell-type × comparison blocks. Pseudobulk rows now carry `padj_local`,
  `padj_global`, and global-FDR counts; narrative reports local vs global
  and ORA uses global-significant genes.
- ✅ **T1.4** (`ab246ca`) — explicit ORA background.
  scRNA and bulk ORA pass dataset-expressed genes; reports surface
  `background_size` and `background_source`.
- ✅ **T1.5** (`ab246ca`) — power statements. NB-Wald
  approximate power is surfaced for pseudobulk and bulk and summarized in
  Methods.
- ✅ **T2.1-T2.5** (`ab246ca` + `efa7136`) — provenance block,
  input SHA-256, nested params hashes, lockfile embed/warning,
  `methodology.json`, `--reproducible`, memory snapshot, and LLM token/cost
  usage table implemented.
- ⚠️ **Stage C first rerun reviewed** —
  `/home/medusa/.aria/reports/aria_20260515_193650_oligodendrocytes_opc_maturemyelinatingol_-824/report.html`.
  Do not tag from this report: it lacks lockfiles, LLM token/cost usage, and
  visible per-stage params hashes; global-FDR wording was also confusing.
  Fixes for all but lockfile generation landed in `efa7136`.
- ✅ **PBMC rerun reviewed** —
  `/home/medusa/.aria/reports/aria_20260516_115010_monocytes_tcells_bcells_-87a/report.html`.
  The paired-design, t-SNE, lockfile, params-hash, global-FDR, power, input
  SHA, and LLM usage checks passed, but it predated the final `4.4.0`
  metadata rerun.
- ✅ **Final PBMC v4.4 report reviewed** —
  `/home/medusa/.aria/reports/aria_20260520_124352_monocytes_tcells_bcells_-d72/report.html`.
  The report was generated from `cbcde8e`, shows `aria_version=4.4.0`, and
  passes all seven reviewer questions. It is the report that justifies the
  public `v4.4` tag.

## If Resuming Mid-Sprint

The plan in `memory/roadmap/V44_PUBLICATION_READINESS.md` is the single
source of truth for what is done and what is next. Each item has a Touch
list (files), a Test list (pytest cases), and a Reviewer check. Closed
items have a `✅ COMPLETED` header with the exact files touched and the
test counts; the original spec is preserved below each closure header.

Execution order (do NOT improvise):

1. **Stage A — scientific blockers** (Tier 1): code done.
2. **Stage B — reproducibility blockers** (Tier 2): code done.
3. **Stage C — validation**: PBMC rerun from `7a73c22` is clean and
   pseudobulk DE is non-zero/analyzable. The final release report from
   `cbcde8e` passed the same seven-question answers and shows
   `aria_version=4.4.0`; `v4.4` is tagged at `cbcde8e`.

Validation gates between every item:

```bash
python -m compileall -q aria
python -m pytest -q tests/test_pytest_smoke.py
```

Latest validation at `cbcde8e`: `compileall` clean,
`python -m pytest -q tests/test_pytest_smoke.py` -> 72 passed, 4 skipped,
and the two pydeseq2-gated paired/composition tests pass in `aria-rna-env`
with `NUMBA_CACHE_DIR=/tmp/numba_cache`.
Each item adds new tests; never remove. Note
that some pseudobulk tests require pydeseq2 (skip in `aria-env`, pass in
`aria-rna-env`); that is expected.

Lockfile note: `ac48599` pivoted lockfile generation to
snapshot-from-installed; the scRNA/rna lock pair exists and is sufficient
for the PBMC Stage C rerun.

## v4.5 Raw Ingestion Closeout

Implemented in `1d54cc0` and completed in `a0b33dd`:

- deterministic 10X matrix-triplet discovery and validation;
- `scanpy.read_10x_mtx` conversion to canonical workspace `.h5ad`;
- source/output SHA-256 and ingestion parameter hashes;
- `RawIngestionAgent` integration after SetupAgent and before scRNA dispatch;
- report/methodology Raw Ingestion provenance section;
- FASTQ grouping/planning with explicit blockers for missing chemistry,
  reference/index hash, and `kb-python` tooling;
- gated `execute_kb_count` for the fully explicit FASTQ path: FASTQs,
  index path/hash, transcript-to-gene path/hash, chemistry, output dir, and
  installed `kb` are required before execution;
- version metadata aligned to `4.5.1`.

Validation:

- `python -m compileall -q aria` -> pass
- `python -m pytest -q tests/test_pytest_smoke.py` -> 79 passed, 4 skipped
- `python -c "import aria; print(aria.__version__)"` -> `4.5.1`
- `git diff --check` -> pass

## Next Milestone — v4.6 scATAC

Before starting v4.6, note the post-v4.5.2 hardening was committed and pushed
in `c611e45` (untagged, `aria.__version__` = `4.5.2` at that point) on top of
`4bf9741`:

- GEO/SRA metadata design inference was tightened to prefer experimental
  characteristic keys and preserve GSM/SRR IDs plus title aliases for count
  matrix mapping.
- Bulk preranked GSEA was fixed for count matrices whose gene IDs are already
  symbols and now receives explicit `NarrativeBlock` coverage.
- Narrative rendering now composes deterministic prose from blocks via
  `aria.agents.narrative.compose_prose`, with evidence tables as audit support.
- The versionable dependency/impact graph lives at
  `docs/architecture/code_graph.md`.
- `memory/architecture/PROJECT_ARCHITECTURE.md` was refreshed locally; the
  tracked durable architecture docs are under `docs/architecture/`.

Validation from the hardening pass (re-confirmed at `c611e45` on 2026-05-28
with the `aria-env` interpreter; narrative/smoke tests need `litellm`, present
in `aria-env`, absent in `aria-rna-env`):

- `python -m pytest -q tests/test_pathway_viz.py tests/test_narrator_bulk.py tests/test_narrative_render_blocks.py tests/test_geo_design.py` -> 10 passed
- narrative kernel suite -> 17 passed (was 16 at v4.5.2)
- `python -m pytest -q tests/test_pytest_smoke.py` -> 86 passed, 4 skipped
- `python -m compileall -q aria` -> pass
- `git diff --check` -> pass

Read `memory/roadmap/V46_SCATAC_PLAN.md` first. scATAC must reuse v4.4
provenance/methodology guarantees and v4.5 raw-ingestion patterns rather than
inventing separate input handling.

## Pre-ATAC remediation (2026-05-28) — DO NOT skip before scATAC

A senior audit (`memory/audit/2026-05-28_senior_audit.md`) is the authoritative
tracker. Samael's mandate: solve before v4.6. Executed this session (108
passed / 4 skipped, zero regressions): P1 headless runner + design E2E test
(`aria/headless.py`, `tests/test_headless_design_e2e.py`); P3 causal-guard
hardening (`tests/test_causal_guard.py`); P4a lognorm-recovery provenance; P4b
power-disclosure; P-ENG-PERF env-detection cache (`tests/test_env_manager_cache.py`);
P2 partial CI/Dockerfile.

**PBMC v4.5.2 thin-report blocker — diagnosed, fixed, confirmed.** The thin
report was caused by `scRNAAgent._needs_pseudobulk` silently re-gating an
approved DesignIntelligence pseudobulk plan on free-text keywords. Commit
`27d4b15` fixes the gate so explicit obs design + DI-recommended pseudobulk
runs regardless of keyword phrasing. Confirmation rerun produced
`pseudobulk_de.csv` with 3,376 DE genes across 11 STIM_vs_CTRL blocks plus
differential abundance. A longer-timeout full HTML rerun is still useful for
release evidence, but the scientific-depth regression is closed.

Deferred remediation items (steps in the audit file): IHW/hierarchical FDR,
apeglm shrinkage, bus persistence, integration QA (LISI/kBET), ambient-RNA
correction, god-file splits, FASTQ-path tests.

**External AI audit adjudicated** (`auditoria.txt`, two static-GitHub audits):
verdicts + items X1–X20 are in the audit file. Their P0 "missing
agents / broken core" claims are REJECTED as stale-snapshot artifacts
(design/audit/raw_ingestion/design_intelligence agents and narrative narrators
all exist). VERIFIED-real accepts form a cheap **"v4.5.3 Integrity & Trust"**
mini-milestone before ATAC; it was tagged at `bab6fbd` as `v4.5.3`: X1 central
version source, X2 no API-key writes to `.bashrc`, X3 registry-integrity
checks, X4 chromatin scaffold dispatch gate, and X17 tiered `aria doctor`
(with X5-X9/X14 closed later in the same freeze). Validation: compileall pass,
targeted integrity tests 8 passed, smoke 86 passed / 4 skipped, combined
targeted + smoke 94 passed / 4 skipped, diff-check pass.

**X7 design-matrix validator is closed.** `aria.utils.design_matrix` now
validates rank deficiency, complete condition-covariate confounding,
insufficient condition replicates, n=1 design cells, numeric condition factors,
binary numeric covariates, and text-stored continuous covariates. The same
contract is surfaced in DesignIntelligence/AuditAgent and enforced in bulk RNA
and scRNA pseudobulk scripts before DESeq2. Validation: compileall pass,
X7+narrator 8 passed, narrative+X7 23 passed, integrity 8 passed, smoke 86
passed / 4 skipped, diff-check pass.

**X5 typed IPC contracts are closed.** `aria.utils.script_contracts` defines
pydantic `ScriptContract`s; `EnvironmentManager` validates registered script
inputs before conda, validates output JSON after execution, returns
`InvalidScriptParams` / `IncompatibleScriptContract`, and attaches
`ipc_contract` metadata to valid outputs. `registry_integrity` checks contract
script existence/version coherence. Validation: compileall pass, X5+
registry/doctor 10 passed, combined X5+integrity+X7+narrator 20 passed, smoke
86 passed / 4 skipped, doctor smoke pass with expected HiC warning, diff-check
pass.

**X6 synthetic-DE benchmark is closed.** `aria/benchmarks/synthetic_de.py`
simulates an NB dataset with known true-DE genes and runs ARIA's real
pseudobulk DE, scoring recall / empirical FDR vs tolerances; pydeseq2-gated
test + `aria doctor --benchmark` wiring. Measured recall=1.000 /
empirical_fdr=0.000 (aria-rna-env); aria-env suite 98 passed / 5 skipped.

**X14 Claim Compiler (flagship) is closed.**
`aria/agents/narrative/claim_compiler.py` tiers every claim
(descriptive→associative→weak/strong mechanistic→causal_experimental) from
structured evidence and caps the licensed language; observational omics caps at
associative, causal needs an interventional design; methodology.json gains a
`claims` manifest, HTML shows an evidence-tier badge, policy in ADR-013.
Validation: claim-compiler+guard+narrative 20 passed; full smoke 109 passed /
4 skipped.

**X8/X9 scientific QC are closed.** `aria/utils/integration_qc.py` flags
integration overcorrection / residual batch from silhouettes;
`aria/utils/annotation_qc.py` flags reused obs labels as unverified or lacking
distinct markers (data-driven, no hardcoded marker map, ADR-011). Surfaced via
`scrna_agent` findings + a `scrna.data_quality` narrative block. Validation:
qc 8 passed, narrator+qc 13 passed, full smoke 109 passed / 4 skipped.

Pre-ATAC integrity-freeze status: **X1–X9, X14, X17 are CLOSED.** Remaining
audit items (lower priority, can run alongside or after v4.6): X10 privacy
firewall, X11/X12 DebateCouncil cost + ParameterAdvisor bias guardrails, X13
DataAudit scan limits, X18 spatial scaffold, and the research-track
X15/X16/X19/X20. The original deferred engineering items (bus persistence,
apeglm shrinkage, IHW FDR, god-file splits, FASTQ-path tests) also remain in
the audit file. The pre-ATAC integrity freeze is essentially complete; v4.6
scATAC can begin, inheriting these contracts.
