---
status: active
source_of_truth_for: project_state
last_updated: 2026-05-22
supersedes:
  - archive/sessions/2026-05-14_close.md
  - archive/handoffs/CODEX_ARIA_HANDOFF_2026-05-12.md
---

# ARIA Project State

## Repository

- Local repo: `/home/medusa/Samael/ARIA`
- Branch: `main`
- Remote: `origin/main`
- Verify current HEAD with `git log --oneline --decorate -5`.
- Last pre-memory-reorganization docs commit: `2e0eb39`
  (`Document final v4.3.12 closeout`)
- Last pre-maintenance code commit: `d3de169`
  (`Remove dataset-specific narrative guardrails`)
- Current HEAD: verify with `git log --oneline --decorate -5`; the latest
  pushed line is v4.5.2 plus post-release anti-hardcode hardening.
- v4.5.2 release commit: `ca8b169`
  (`v4.5.2 narrative kernel`).
- Last stable tag: `v4.5.2`
  (`v4.5.2 narrative kernel`).
- Previous stable tag: `v4.5.1` (`a0b33dd`,
  `v4.5.1 add gated kb ingestion execution`).
- Previous base tag: `v4.5` (`1d54cc0`,
  `v4.5 raw ingestion bridge`).
- P0 audit fixes landed in `05f6f4e`
  (`Close P0 audit findings for v4.3.19`).
- P1 audit fixes landed in `3054318` (untagged) — safe h5ad helper, bus
  fan-out try/except, `failed/` eviction flock, `decisions.checkpoint` TEXT
  migration, kb hash blocker, RNA script resume validity guards.
- Post-v4.5.2 anti-hardcode hardening neutralized new narrative-kernel
  synthetic tests and added ADR-011 / design-principle language: runtime and
  generic tests must not hardcode biological content unless explicitly marked
  as named golden-dataset regressions.
- Tags:
  - `v4.3.12` -> `3a0c40e`
  - `v4.3.12.post1` -> `805e0b2`
  - `v4.3.13` -> `9094cb1`
  - `v4.3.14` -> `f4679d7`
  - `v4.3.15` -> `bf1b65b`
  - `v4.3.16` -> `5812741`
  - `v4.3.17` -> `e93707f`
  - `v4.3.18` -> `9dc48aa`
  - `v4.3.19` -> `eece228`
  - `v4.4` -> `cbcde8e`
  - `v4.5` -> `1d54cc0`
  - `v4.5.1` -> `a0b33dd`
  - `v4.5.2` -> narrative kernel release commit

Do not move existing tags. Use a new patch tag if one is ever needed.

## Current Release Boundary

The 4.3 line remains closed to broad new modality features. `v4.4` is the
latest stable scientific/reproducibility baseline for RNA workflows and report
provenance. `v4.3.12` remains the historical stable baseline for bulk RNA and
scRNA workflows; `v4.3.12.post1` marks report-fidelity fixes; `d3de169`
removed dataset-specific runtime narrative guardrails after the post1 tag.
`v4.3.18` made Design Intelligence choices explicit. `v4.3.19` closes four
P0 audit findings against documented principles and aligns release docs/version
metadata.

## v4.3.13-v4.3.19 Maintenance Patches

- `v4.3.13`: initial TUI long-prompt / covariate display patch.
- `v4.3.14`: explicit multiline TUI prompt with `END` sentinel and h5ad
  human-species inference.
- `v4.3.15`: focused scRNA h5ad pre-QC subsetting for requested cell groups.
- `v4.3.16`: tightened scRNA focus, LIANA p-value display, integrated
  interpretation prompt sanitization, and richer scRNA decision logging.
- `v4.3.17`: cross-modality Design Intelligence layer before compute.
- `v4.3.18`: CP2 explicitly chooses recommended-only vs
  recommended+optional supported analyses, and cell-focus detection only uses
  explicit focus/restriction clauses.
- `v4.3.19` (`eece228`; P0 code in `05f6f4e`): four P0 audit fixes — dynamic
  pseudobulk DE narrative header (was hardcoded "Age-associated"); chromatin
  QC fallback now gated by `mocks_allowed()`; `rna_cellcomm` no longer
  swallows non-ImportError LIANA failures; `hic_inspect` reads
  `RAM_ESTIMATES_GB` from `aria.utils.hic_constants` instead of importing an
  agent (script→agent layer violation removed). Release docs, TUI, `setup.py`,
  and `aria/llm.__version__` are aligned to `4.3.19`.

Design Intelligence:

- Runs after experimental design confirmation and before the final plan
  checkpoint.
- Produces recommended, optional/supported, unsupported/not-recommended, and
  warning items.
- Is rules-first and currently covers scRNA, bulk RNA, chromatin/scATAC,
  Hi-C, and multimodal integration feasibility.
- Logs a `2.pre` decision so reports expose the design assessment.
- scRNA uses the assessment to avoid unsupported LIANA or PAGA/DPT runs when
  the selected design does not justify them.
- Option 1 in CP2 runs only recommended analyses; option 2 adds optional
  supported analyses. Unsupported/not-recommended items stay blocked unless the
  user modifies the plan explicitly.

Validation:

- `python -m compileall -q aria` -> pass
- `python -m pytest -q tests/test_pytest_smoke.py` -> 40 passed

## Production-Like Validated

- Bulk RNA-seq count-matrix workflow.
- scRNA single-sample workflow.
- scRNA multi-sample workflow.
- Processed `.h5ad` pseudobulk workflow using `obs` design metadata.
- Narrative reports with figures and supplementary TSV exports.

## Validated / Beta

- Bulk RNA FASTQ preprocessing.
- PAGA + DPT trajectory context.
- LIANA cell-cell communication.
- GEO/SRA connector path.

## Scaffolded / Roadmap

- scATAC-seq.
- Bulk ATAC-seq.
- ChIP-seq / CUT&RUN / CUT&TAG.
- Hi-C / Micro-C.
- WNN / MOFA+ / peak-to-gene integration.

## Final v4.3.12 Closeout

Core fixes:

- DataAudit infers condition, replicate, groupby, and covariates from `.h5ad`
  `obs`.
- `rna_qc.py` handles processed h5ads using existing `obs` QC metrics.
- ARIA filters stale intermediate h5ads when real source h5ads exist.
- Existing `obs` labels, such as `obs['subclass']`, can be reused; Leiden and
  CellTypist are skipped when a trusted groupby column is provided.
- Pseudobulk donor-level DE is the primary inferential layer for scRNA
  condition contrasts.

Report fixes:

- Reports describe reused `obs` labels truthfully instead of claiming Leiden.
- `subclasss` typo fixed.
- PAGA/DPT language avoids active-differentiation claims without velocity or
  time-course data.
- Per-cluster Wilcoxon timeout is surfaced explicitly.
- HTML report version is read from `aria.__version__`.
- Methods and decisions are HTML-escaped.
- Raw Python error dicts are converted to human-readable report text.
- LIANA table shows rank, metric value, and metric name.
- Bulk interpretation anti-causality guardrails are generic.
- Runtime code and prompts contain no dataset-specific guardrails.

## Latest Validation

After the latest v4.4 release commit (`cbcde8e` — version + methodology
tool provenance):

- `python -m compileall -q aria` -> pass
- `python -m pytest -q tests/test_pytest_smoke.py` -> 72 passed, 4 skipped
- `NUMBA_CACHE_DIR=/tmp/numba_cache /home/medusa/anaconda3/envs/aria-rna-env/bin/python -m pytest -q tests/test_pytest_smoke.py::test_pseudobulk_paired_design_uses_donor_condition_pseudosamples tests/test_pytest_smoke.py::test_pseudobulk_composition_covariate_flag_propagates_into_blocks` -> 2 passed
- `python tests/test_bulk_rna.py` -> 30 passed, 0 failed
- `python -c "import aria; print(aria.__version__)"` -> `4.4.0`
- `git diff --check` -> pass
- `envs/aria-rna-env.linux-64.lock` produced (136 conda packages,
  @EXPLICIT format, snapshot-from-installed)
- `envs/aria-rna-env.pip.lock` produced (86 pip packages)
- Live verification on `/home/medusa/Samael/datasets/pbmc.h5ad`
  (Kang et al. IFN-β, 8 donors): condition_col=stim, replicate_col=Donor,
  groupby_col=cluster, comparisons=[[STIM,CTRL]], confidence=high.

## Latest Reviewed scRNA Report

`/home/medusa/.aria/reports/aria_20260514_143352_oligodendrocytes_opcs_microglia_-009/report.html`

Reviewed outputs:

- 295,033 starting cells.
- 242,405 cells retained after QC.
- 18 subclass groups reused from `obs['subclass']`.
- Pseudobulk DE: 79 analyzable blocks, 57 significant.
- Pathway ORA: 38/57 DE blocks with enrichment.
- LIANA: 50 non-autocrine interactions.
- PAGA/DPT generated.
- RNA velocity skipped because no spliced/unspliced layers were available.

The reviewed HTML predates the final report-fidelity hotfixes, so newly
generated reports should be cleaner than that specific file.

## v4.4 Publication Readiness Closeout

Decided 2026-05-15 after the expert audit: ARIA paused new modalities and
spent `v4.4` closing the scientific and reproducibility gaps that block peer
review. Full step-by-step plan lives in
`memory/roadmap/V44_PUBLICATION_READINESS.md`. The sprint is closed and
tagged after the final real PBMC Stage C report from `cbcde8e` passed the
seven-question peer-reviewer test.

### v4.4 sprint progress

| Item | Status     | Notes |
|------|------------|-------|
| T1.3 | ✅ done    | Commit `e8fc846`. Default n=3; low_power_warning at n=2; DI downgrade. 4 new tests; 42 pytest + 30 legacy bulk + 2 rna-env-gated all green. |
| T1.1 | ✅ done    | Commit `52bf703`. New `rna_diff_abundance.py` (Poisson-offset GLM + Fisher fallback); `composition_covariate` param in pseudobulk DE; scrna_agent wires DA before DE and gates the covariate on `any_significant`; new "Cell-type abundance" narrative section; DE header reports composition-corrected count. 4 new tests; 45 pytest + 30 legacy + 3 rna-env-gated all green. |
| T1.2 | ✅ done | Commit `ab246ca`. `rna_pseudobulk_de.py` now carries `padj_local`, `padj_global`, `multiple_testing.n_tests_global`, and global-FDR counts; narrative and ORA use global-FDR significant genes by default. Wording tightened in `efa7136`. |
| T1.4 | ✅ done | Commit `ab246ca`. scRNA and bulk ORA pass dataset-expressed backgrounds; reports surface `background_size` / `background_source`; gseapy background incompatibility falls back visibly. |
| T1.5 | ✅ done | Commit `ab246ca`. `aria.utils.power_estimation` adds NB-Wald approximate power; pseudobulk and bulk blocks surface `power_estimate_at_lfc_min`; methods report power ranges. |
| T2.1-T2.5 | ✅ done | `ab246ca` implemented provenance block, input SHA-256, params hashes, lockfile embed/warning, `methodology.json`, `--reproducible`, and memory snapshot. `efa7136` tightened Stage C report blockers: nested params hashes, LLM token/cost report table, clearer global-FDR wording, scRNA params-hash propagation, and `idr` moved to conda for the chromatin env. `ac48599` closed the residual: `inputs` are now persisted in `methodology.json`, and `generate_locks.sh` switched to snapshot-from-installed (`conda list --explicit` + `pip freeze`) so the solver hang no longer blocks the v4.4 closeout. Final PBMC Stage C real-data validation passed on 2026-05-20. |

### Stage C first rerun review

Reviewed:
`/home/medusa/.aria/reports/aria_20260515_193650_oligodendrocytes_opc_maturemyelinatingol_-824/report.html`.

Result: **do not tag `v4.4` from this report.** It answers git commit,
composition correction, local/global multiple testing, replicate/power, and
input SHA-256. It does not answer exact dependency lockfiles, LLM token/cost
usage, or per-stage parameter hashes in the HTML. The report also used
confusing global-FDR wording ("remained significant") that was corrected in
`efa7136`.

### Stage C closeout fixes (commit `ac48599`)

A second-pass audit caught two more issues that landed in `ac48599`:

1. `methodology.json` did NOT persist `inputs` even though the HTML
   rendered the SHA-256 table. Fixed: `_build_methodology_json` now
   returns the `input_files` list under the `inputs` key, with two new
   pytest cases (`_persists_input_hashes` and an extended
   `_emitted_with_required_keys`).
2. The conda-lock-based approach hung not only on `aria-chromatin-env`
   but on any env mixing bioconda + pip on this machine. Pivoted to
   `conda list --explicit` snapshot-from-installed in
   `scripts/generate_locks.sh`, which (a) does not invoke the solver,
   (b) writes the same `@EXPLICIT` lockfile format conda-lock would,
   (c) pairs every conda lock with a `pip freeze` sibling
   (`envs/<env>.pip.lock`). Only `aria-rna-env` is installed locally,
   so only its lock pair was produced; `hic`, `integration`, and
   `chromatin` are skipped with explicit "deferred" messages because
   v4.4 Stage C is scRNA-only and does not need them.
   `_build_lockfile_section` was extended to embed pip locks next to
   conda locks and to note "no pip packages" rather than silence when
   a conda lock has no pip sibling.

After `ac48599` the seven-question test now has answers wired for all 7
points; the next rerun should pass review.

### PBMC IFN-β Stage C rerun review (2026-05-15 22:22)

A second TUI rerun against `/home/medusa/Samael/datasets/pbmc.h5ad`
landed at
`/home/medusa/.aria/reports/aria_20260515_222216_monocytes_tcells_bcells_-474/`.
Provenance is correct (git_sha = `0957d32`, input SHA preserved,
`inputs` key populated in `methodology.json`, lockfile pair embedded).
**Four blockers prevent the v4.4 tag** — see
`memory/roadmap/V44_PBMC_BLOCKERS.md` for full forensics and the
proposed execution order. Headline: pseudobulk DE silently fails on
paired designs (Donor as the aggregation key collapses STIM and CTRL
aliquots into one pseudosample), ARIA writes intermediates inside the
user's data_dir, `rna_figure_umap.py` crashes on t-SNE-only h5ads, and
`methodology.json["llm_usage"]` is null despite real LLM calls.

Bug 1 and Bug 2 are closed in `ba4e21e`. Pseudobulk and differential
abundance now use a `replicate × condition` pseudosample key for paired
designs, and balanced paired DESeq2 blocks add donor as a covariate.
scRNAAgent now passes explicit per-experiment workspace output directories
to QC, integration, clustering, annotation, DE, ORA, pseudobulk, trajectory,
and LIANA; `rna_clustering.py` honors `output_dir`; DataAudit filters
ARIA-generated intermediates across modalities.

Bug 3 and Bug 4 are closed in `7a73c22`. Embedding figures now fall back
from UMAP to t-SNE or any other 2D `X_*` embedding and label figures/captions
with the actual embedding type. `methodology.json` now persists the same
`llm_usage` summary rendered in HTML, and usage collection applies a small
timestamp grace window to include near-start LLM calls.

PBMC rerun reviewed 2026-05-16:
`/home/medusa/.aria/reports/aria_20260516_115010_monocytes_tcells_bcells_-87a/report.html`.
The rerun was generated from clean `7a73c22` provenance and answers the
seven-question audit in the HTML: exact commit, conda/pip lockfiles,
composition correction (`~ stim + Donor` for balanced paired blocks), local
and global BH (`n_tests_global=127643`), replicate/power summaries, input
SHA-256, and LLM usage (3 cache hits, 0 tokens/cost). Pseudobulk now produced
11 analyzable significant cell-type blocks and ORA ran on all 11.

Release metadata closeout landed in `cbcde8e`: ARIA/package/TUI version now
reports `4.4.0`, and `methodology.json["tools"]` falls back to committed
lockfiles for packages such as `pydeseq2` and `gseapy` when the report writer
runs outside `aria-rna-env`.

Final PBMC Stage C report reviewed 2026-05-20:
`/home/medusa/.aria/reports/aria_20260520_124352_monocytes_tcells_bcells_-d72/report.html`.
The report was generated from `cbcde8e` and passes the seven-question audit in
the HTML: `aria_version=4.4.0`, exact git SHA
`cbcde8ea4b2e114ef0829bbab3105793ddbc3395`, embedded conda/pip lockfiles,
lockfile-derived tool versions (`pydeseq2=0.5.4`, `gseapy=1.1.13`), input
SHA-256
`af0696e90f9abd4904e370b4abf14a0f88cd9c3df7a8455ef7beb102369e0d9a`,
per-stage parameter hashes, composition correction enabled for 11/11
analyzable pseudobulk blocks, local and global BH
(`n_tests_global=127643`), min replicates per condition 3, approximate power
8%-13%, and LLM usage (3 cache hits, 0 tokens, $0 estimated cost). Pseudobulk
DE produced 11 analyzable significant blocks, 3,376 global-FDR genes, ORA ran
on all 11 blocks with 14,044 dataset-expressed background genes, and t-SNE
figures rendered with honest labels. `git_dirty=True` in provenance is
explained by an unrelated untracked `codigo_aria.txt`; tracked repository
files had no diff at report generation review.

### Generalization regression: PBMC IFN-β (commit `0957d32`)

A separate audit revealed that ARIA had been over-fit to the hippocampus
dataset and failed on the public Kang et al. PBMC IFN-β stimulation
dataset (`/home/medusa/Samael/datasets/pbmc.h5ad`). Root causes were all
inside the obs-inference layer:

1. Vocabulary was too narrow — `stim`, `stimulation`, `perturbation`,
   `state`, `intervention` were not in the condition priority list.
2. Column matching was case-sensitive — `Donor` (capitalized) missed.
3. Groupby priority did not include `cluster`.
4. `_usable_replicate_col` allowed n=2 levels, which let
   Seurat's `orig.ident` (a duplicate of the condition) get picked as
   the replicate column.
5. CP 2.1 manual-group example was hardcoded to hippocampus donor
   barcodes.
6. CP 2.3 factor menu would surface LLM hallucinated sentences as
   factor options (e.g. "need file names... paste your rna-seq file").
7. `tui.py` "Other" example mentioned a specific bird species.

All seven were fixed in `0957d32` with three new pytest cases (PBMC
golden fixture, sanitizer guard, case-insensitive replicate). Live
re-verification confirms ARIA now infers stim/Donor/cluster with high
confidence and no warnings on the real Kang dataset.

The v4.4 publication-readiness bar is now: any new scRNA dataset whose
obs follows standard naming conventions (single-cell vocabulary covered
by the widened priority lists) should pass CP1/CP2.1/CP2.3 without
manual override. The PBMC fixture acts as the regression guard.

Summary of scope:

**Tier 1 — scientific corrections (blocking):**

- T1.1 Differential abundance + composition correction for scRNA pseudobulk.
- T1.2 Global FDR across all cell-type × comparison blocks (BH).
- T1.3 Default `min_replicates_per_condition = 3`; escalated warning at n=2.
- T1.4 Pathway ORA with explicit dataset-expressed background.
- T1.5 Automated power and effect-size statement in methods section.

**Tier 2 — reproducibility (blocking):**

- T2.1 Provenance block embedded in every HTML report (git SHA, versions).
- T2.2 SHA-256 hashes for inputs and per-stage params.
- T2.3 conda-lock files generated and embedded.
- T2.4 `methodology.json` exported alongside the HTML.
- T2.5 `--reproducible` mode (deterministic byte-identity + memory snapshot).

Done criteria passed on the final real scRNA rerun generated from `cbcde8e`;
`v4.4` is tagged at that commit. Raw Ingestion becomes `v4.5`
(see `V45_RAW_INGESTION_PLAN.md`) so users can start from Cell Ranger-style
matrices or FASTQs without leaving ARIA; scATAC moves to `v4.6`
(see `V46_SCATAC_PLAN.md`) so it inherits both v4.4 and v4.5 guarantees.

Tier 3 (operability) and Tier 4 (design depth) are explicitly out of scope
for `v4.4` and documented at the bottom of the plan file.

## Roadmap After v4.4

- `v4.5` — Raw Ingestion (`V45_RAW_INGESTION_PLAN.md`) closed in `1d54cc0`.
  Cell Ranger-style 10X matrix triplets are ingested to canonical workspace
  `.h5ad` files with hashes/provenance; reports expose raw ingestion records;
  FASTQ inputs are detected/grouped and blocked until chemistry, reference,
  index hash, and `kb-python` tooling are explicit.
- `v4.6` — scATAC E2E (`V46_SCATAC_PLAN.md`, moved from v4.5), next.
- `v4.7` — bulk ATAC E2E (`V47_BULK_ATAC_PLAN.md`, moved from v4.6).
- `v4.8` — multimodal integration WNN + MOFA+ + peak2gene
  (`V48_INTEGRATION_PLAN.md`, moved from v4.7).

## v4.5 Raw Ingestion Closeout

Commit/tag: `v4.5.1` -> `a0b33dd`
(`v4.5.1 add gated kb ingestion execution`). Base `v4.5` remains at
`1d54cc0` and was not moved.

Implemented:

- `aria.utils.raw_ingestion`: deterministic 10X triplet discovery,
  gzip/Matrix Market/dimension validation, `scanpy.read_10x_mtx` conversion,
  `.h5ad` output hashing, FASTQ grouping/blocker plan, explicit `kb count`
  command builder, and gated `execute_kb_count` that runs only after FASTQs,
  index path/hash, transcript-to-gene path/hash, chemistry, output dir, and
  `kb` tooling are explicit.
- `aria.agents.raw_ingestion_agent`: scans the data directory after setup,
  writes canonical `.h5ad` files to `~/.aria/workspace/<run_id>/ingested/`,
  updates `exp_context["modalities"]["scRNA"]`, and records ingestion
  provenance.
- `OrchestratorAgent`: runs RawIngestionAgent after SetupAgent and before
  modality dispatch.
- `NarrativeAgent`: writes `raw_ingestion` into `methodology.json` and renders
  a Raw Ingestion provenance table in HTML.
- TUI labels for Raw Ingest.
- Version metadata aligned to `4.5.1`.

Validation:

- `python -m compileall -q aria` -> pass
- `python -m pytest -q tests/test_pytest_smoke.py` -> 79 passed, 4 skipped
- `python -c "import aria; print(aria.__version__)"` -> `4.5.1`
- `git diff --check` -> pass

## Post-v4.5.1 Audit Hardening

Closed the live P1 audit findings from 2026-05-15 before starting `v4.6`.
Code landed in commit `3054318` (`Post-v4.5.1 audit hardening: close P1
findings`, untagged) on top of `v4.5.1`:

- `MessageBus.publish` now logs receiver exceptions and continues dispatching
  to later subscribers.
- RNA `.h5ad` resume logic now validates cached h5ads before returning cached
  summaries; RNA script readers use `aria.utils.safe_h5ad.read_h5ad` for
  explicit corrupt-file errors.
- Failed-run eviction under `workspace/failed/` now uses `failed/.lock` with
  `fcntl.flock`.
- `memory.decisions.checkpoint` is normalized to `TEXT`, with an in-place
  migration for existing SQLite tables.
- `execute_kb_count` now reports index/t2g hashing failures as structured
  blockers instead of raising through the orchestrator.
- The legacy bulk RNA smoke wrapper now rejects printed nonzero internal
  failure counts. Native pytest conversion and golden DE fixtures remain open
  P2 follow-up work.

Validation at `3054318`:

- `python -m compileall -q aria` -> pass
- `python -m pytest -q tests/test_pytest_smoke.py` -> 85 passed, 4 skipped
- `python -c "import aria; print(aria.__version__)"` -> `4.5.1`

## Post-v4.5.1 Narrative Depth Work

Status: implemented locally on top of `3054318`, not yet committed or tagged.

Rationale: the scRNA NarrativeAgent was already peer-reviewable for
reproducibility and honesty, but its main prose often stopped at aggregate
counts. The narrative now adds deterministic, result-local interpretation
before broad synthesis.

Implemented:

- `aria/agents/_narrative_scrna.py` now describes top pseudobulk DE blocks
  with global/local FDR counts, top up/down genes, composition-correction
  status, approximate power, matched ORA support, and caveats.
- Differential abundance is connected back to matched DE blocks so reports
  explain when composition shifts were modeled by a log-proportion covariate.
- Pathway, LIANA, and trajectory sections now add local interpretation:
  ORA terms per block, top non-autocrine ligand-receptor candidates with rank
  metric, PAGA strength, and DPT ordering limits.
- `build_scrna_integrated_interpretation` uses the same deterministic helpers
  so the final synthesis is deeper without giving the LLM authority over
  which analyses ran.
- `tests/test_pytest_smoke.py::test_scrna_narrative_adds_per_result_depth`
  locks the new behavior with a synthetic DE/ORA/LIANA/trajectory fixture.

Validation after the narrative-depth work:

- `python -m compileall -q aria` -> pass
- Targeted narrative tests:
  `python -m pytest -q tests/test_pytest_smoke.py::test_scrna_narrative_adds_per_result_depth tests/test_pytest_smoke.py::test_integrated_interpretation_sanitizes_long_prompt tests/test_pytest_smoke.py::test_scrna_narrative_reports_predefined_groupby_not_leiden tests/test_pytest_smoke.py::test_cellcomm_table_labels_rank_metric_not_generic_score`
  -> 4 passed
- `python -m pytest -q tests/test_pytest_smoke.py` -> 86 passed, 4 skipped
- `git diff --check` -> pass

## v4.5.2 Narrative Kernel Closeout

`v4.5.2` promotes the post-v4.5.1 narrative-depth work into a structured
Narrative Kernel before v4.6 scATAC begins.

Implemented:

- New `aria.agents.narrative` package:
  - `types.py`: `EvidenceItem`, `Caveat`, and `NarrativeBlock` dataclasses;
  - `protocols.py`: `ModalityNarrator` protocol;
  - `registry.py`: narrator registration and block collection;
  - `validators.py`: integrity validators for claims, evidence, failed
    analyses, low/insufficient visibility, causal-language downgrade,
    PAGA/DPT caveats, and file existence;
  - `render_blocks.py`: HTML composer for narrative blocks.
- New narrators:
  - `ScrnaNarrator` wraps `_narrative_scrna.py` and emits blocks for QC,
    marker-discovery errors, composition, pseudobulk DE, ORA, LIANA, and
    trajectory.
  - `BulkRnaNarrator` emits QC, contrast, pathway, and power blocks.
- `NarrativeAgent` now uses the registry for scRNA and bulk RNA. If blocks
  exist for a modality, the block composer renders that modality; modalities
  without blocks keep the legacy fallback path.
- `methodology.json` now includes serialized `narrative_blocks`.
- `rna_narrative_adapter.py` now persists input file SHA-256 records for
  offline harness-rendered reports.
- Version metadata aligned to `4.5.2`.
- Documentation updated in `docs/architecture/reporting_and_outputs.md` and
  `docs/release_notes_v4.5.2.md`.

Validation:

- `python -m compileall -q aria` -> pass
- `python -m pytest -q tests/test_narrative_types.py tests/test_narrative_validators.py tests/test_narrator_scrna.py tests/test_narrator_bulk.py tests/test_narrative_render_blocks.py` -> 16 passed
- `python -m pytest -q tests/test_pytest_smoke.py` -> 86 passed, 4 skipped

Release caveat: the final PBMC report rerun must be recorded here before the
`v4.5.2` tag is considered fully closed.

## Post-v4.5.2 Narrative/GEO/GSEA Hardening

Implemented on top of `4bf9741`, pending commit at the time of this memory
update:

- `render_blocks.py` no longer presents block findings primarily as raw
  claim/evidence tables. New `aria.agents.narrative.compose_prose` composes
  deterministic prose per `NarrativeBlock` analysis type, while structured
  evidence remains collapsible audit support.
- `BulkRnaNarrator` now emits explicit preranked GSEA narrative blocks when
  `rna_bulk_de.py` produced GSEA tables or running-sum/top-table figures.
- `rna_pathway_viz.py` now keeps already-symbolic gene identifiers when no
  Ensembl-to-symbol `symbol_map` is available, fixing the preranked-GSEA path
  that previously dropped all genes for symbol-indexed count matrices.
- GEO/SRA design inference now prefers experimental characteristic keys
  (`condition`, `treatment`, `group`, `genotype`, etc.), uses GSM/SRR IDs as
  canonical samples, and stores sample-title aliases so bulk count matrix
  columns can map by accession or readable title.
- `DesignAgent` preserves external `sample_aliases`, and `BulkRNAAgent`
  uses them while applying confirmed design to count matrices.
- New versionable impact map: `docs/architecture/code_graph.md`, linked from
  `docs/README.md` and `docs/architecture/overview.md`.
- `memory/architecture/PROJECT_ARCHITECTURE.md` was refreshed locally with
  current approximate line counts, the narrative kernel package, and visible
  debt notes for `narrative_agent.py` and `_narrative_scrna.py`.

Validation:

- `python -m pytest -q tests/test_pathway_viz.py tests/test_narrator_bulk.py tests/test_narrative_render_blocks.py tests/test_geo_design.py` -> 10 passed
- `python -m pytest -q tests/test_pytest_smoke.py` -> 86 passed, 4 skipped
- `python -m compileall -q aria` -> pass
- `git diff --check` -> pass
