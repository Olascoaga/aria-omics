---
status: active
source_of_truth_for: next_session
last_updated: 2026-05-22
supersedes:
  - archive/sessions/2026-05-14_close.md
---

# NEXT SESSION

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
- Current HEAD: v4.5.2 narrative-kernel release commit after closeout.
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
- Last stable tag: `v4.5.2`
  (`v4.5.2 narrative kernel`).
- Previous stable tag: `v4.5.1` (`a0b33dd`,
  `v4.5.1 add gated kb ingestion execution`).
- Previous base tag: `v4.5` (`1d54cc0`,
  `v4.5 raw ingestion bridge`).
- P0 audit fixes landed in `05f6f4e`
  (`Close P0 audit findings for v4.3.19`).
- P1 audit fixes landed in `3054318` on top of `v4.5.1` (untagged) —
  see "Audit Follow-ups After v4.5.1" below for the full list.
- **`v4.5.2` Narrative Kernel is the current pre-v4.6 reporting baseline.**
  scRNA and bulk RNA report sections are composed from validated
  `NarrativeBlock` objects via modality narrators. See "Narrative Kernel
  v4.5.2" below.
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

Read `memory/roadmap/V46_SCATAC_PLAN.md` first. scATAC must reuse v4.4
provenance/methodology guarantees and v4.5 raw-ingestion patterns rather than
inventing separate input handling.
