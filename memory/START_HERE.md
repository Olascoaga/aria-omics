---
status: active
source_of_truth_for: entrypoint
last_updated: 2026-05-22
---

# ARIA START HERE

Read this file first in every ARIA session.

## Current Status

- The 4.3 line is closed.
- Current branch: `main`.
- Verify current HEAD with `git log --oneline --decorate -5`.
- Last pre-memory-reorganization docs commit: `2e0eb39`
  (`Document final v4.3.12 closeout`).
- Last pre-v4.3 maintenance baseline code change: `d3de169`
  (`Remove dataset-specific narrative guardrails`).
- Current HEAD: v4.5.2 narrative-kernel release commit after closeout.
- Last stable tag: `v4.5.2`
  (`v4.5.2 narrative kernel`).
- Previous stable tag: `v4.5.1` (`a0b33dd`,
  `v4.5.1 add gated kb ingestion execution`).
- Previous base tag: `v4.5` (`1d54cc0`,
  `v4.5 raw ingestion bridge`).
- P0 audit fixes landed in `05f6f4e`
  (`Close P0 audit findings for v4.3.19`).
- P1 audit fixes landed in `3054318`
  (`Post-v4.5.1 audit hardening: close P1 findings`): safe h5ad readability
  helper, bus fan-out resilience, failed/ eviction flock,
  `decisions.checkpoint` TEXT migration, kb hash blocker, and RNA-script
  resume guards. Untagged.
- `v4.5.2` Narrative Kernel: `NarrativeAgent` now composes scRNA and bulk RNA
  report sections from validated `NarrativeBlock` objects via modality
  narrators. `methodology.json` persists serialized narrative blocks.
  Validation: `compileall` pass, narrative-kernel tests 16 passed, full
  `tests/test_pytest_smoke.py` 86 passed / 4 skipped. PBMC rerun closeout
  should be checked before starting v4.6 implementation work.
- Last previous tagged release: `v4.3.18` (`9dc48aa`).
- `v4.3.12` tag: `3a0c40e`.
- `v4.3.12.post1` tag: `805e0b2`.
- `v4.4` Publication Readiness is closed and tagged after a real PBMC Stage C
  report generated from `cbcde8e` passed the seven-question peer-reviewer
  audit. Full plan:
  `memory/roadmap/V44_PUBLICATION_READINESS.md`. Tier 1 (composition
  correction, global FDR, n=3 default, ORA background, power statement) +
  Tier 2 (provenance, hashes, conda lock, methodology.json, reproducible
  mode).
- **v4.5 Raw Ingestion is closed and tagged (`v4.5.1`).** It adds deterministic
  ingestion for Cell Ranger-style 10X matrix triplets, canonical workspace
  `.h5ad` generation with hashes/provenance, report/methodology surfacing,
  a FASTQ detection/planning bridge that blocks on missing chemistry,
  references, index hashes, and `kb-python` tooling, plus gated `kb count`
  execution when all deterministic inputs are explicit. Plan:
  `memory/roadmap/V45_RAW_INGESTION_PLAN.md`.
- scATAC is next as `v4.6` (`memory/roadmap/V46_SCATAC_PLAN.md`) and must
  inherit both v4.4 publication readiness and v4.5 raw-ingestion guarantees.
- **Sprint progress:**
  - T1.3 closed in commit `e8fc846` (n=3 default + low_power_warning + DI).
  - T1.1 closed in commit `52bf703` (`rna_diff_abundance.py` + composition
    covariate in pseudobulk DE + scrna_agent wiring + "Cell-type abundance"
    narrative section).
  - T1.2, T1.4, T1.5, and T2.1-T2.5 implemented in `ab246ca`
    (global FDR, explicit ORA background, power estimates,
    provenance/hashes/lockfile embed, methodology.json, and reproducible
    mode).
  - Stage C first rerun reviewed:
    `/home/medusa/.aria/reports/aria_20260515_193650_oligodendrocytes_opc_maturemyelinatingol_-824/report.html`.
    It does **not** close `v4.4`: lockfiles are missing, LLM usage is not in
    that report, per-stage parameter hashes are not shown, and global-FDR
    wording needed tightening.
  - Follow-up fixes landed in `efa7136`: nested parameter hashes in
    provenance, LLM token/cost accounting and report table, clearer global
    FDR wording, scRNA params hash propagation, and `idr` moved from pip to
    conda in the chromatin env. Validation is green, but conda-lock solving
    for `aria-chromatin-env.yml` hung and no lockfiles were generated.
  - Closeout in `ac48599`: `inputs` persisted in `methodology.json`,
    `scripts/generate_locks.sh` pivoted to snapshot-from-installed
    (conda list --explicit + pip freeze) — solver no longer in the loop.
    `envs/aria-rna-env.linux-64.lock` (136 pkgs) and `.pip.lock` (86 pkgs)
    committed. hic / integration / chromatin envs are not installed
    locally and are skipped with explicit deferred messages; the v4.4
    Stage C rerun is scRNA-only and does not need them.
  - **PBMC Stage C rerun reviewed 2026-05-15 22:22 — v4.4 NOT ready.**
    Four blockers documented in `memory/roadmap/V44_PBMC_BLOCKERS.md`.
  - Bug 1 + Bug 2 closed in `ba4e21e`: paired designs now aggregate
    pseudobulk by `replicate × condition` and add donor as a DESeq2
    covariate for balanced paired blocks; scRNA intermediates are written
    under ARIA workspace paths instead of the user's data directory; the
    audit layer ignores ARIA-generated contamination.
  - Bug 3 + Bug 4 closed in `7a73c22`: embedding figures fall back from
    UMAP to t-SNE or any other 2D `X_*` embedding with honest labels, and
    `methodology.json` now persists `llm_usage` with a timestamp grace
    window for near-start LLM calls.
  - PBMC rerun reviewed:
    `/home/medusa/.aria/reports/aria_20260516_115010_monocytes_tcells_bcells_-87a/report.html`.
    Data dir was clean, pseudobulk paired design produced 11 analyzable
    significant blocks, t-SNE fallback rendered honestly, lockfiles/input
    SHA/params hashes/global FDR/power/LLM usage were visible. The report
    was produced from `7a73c22` and still showed `aria_version=4.3.19`;
    release metadata was finalized in `cbcde8e` by bumping ARIA/TUI/setup to
    `4.4.0` and making `methodology.json["tools"]` fall back to lockfiles
    when the report process is not running inside the RNA env.
  - Final PBMC Stage C report reviewed 2026-05-20:
    `/home/medusa/.aria/reports/aria_20260520_124352_monocytes_tcells_bcells_-d72/report.html`.
    It was generated from `cbcde8e`, shows `aria_version=4.4.0`, embeds
    conda/pip lockfiles and lockfile-derived `pydeseq2=0.5.4` /
    `gseapy=1.1.13`, reports composition correction, local/global BH
    (`n_tests_global=127643`), replicate/power, input SHA-256, per-stage
    parameter hashes, and LLM usage (3 cache hits, 0 tokens/cost).
    `git_dirty=True` is attributable to an unrelated untracked
    `codigo_aria.txt`; tracked files had no diff.

## Do Not

- Do not move existing tags.
- Before starting v4.6 scATAC, read `memory/roadmap/V46_SCATAC_PLAN.md`;
  it inherits v4.4 and v4.5 guarantees.
- Do not add dataset-specific runtime guardrails.
- Do not hardcode genes, perturbations, datasets, or report-rescue phrases.
- Do not allow silent mocks or fake scientific outputs.
- Do not treat files in `memory/archive/` as source of truth.

## Read Order

1. `memory/PROJECT_STATE.md`
2. `memory/NEXT_SESSION.md`
3. `memory/AGENT_PROTOCOL.md`
4. Task-specific file only if needed.

## Source-Of-Truth Hierarchy

1. `PROJECT_STATE.md` for current version, tags, validation, and status.
2. `NEXT_SESSION.md` for the next live handoff.
3. `DECISIONS.md` for accepted architecture and scientific-policy decisions.
4. `AGENT_PROTOCOL.md` for collaboration rules.
5. `REFERENCE_PATHS.md` for local paths.
6. `datasets/`, `architecture/`, and `roadmap/` for task-specific context.
7. `archive/` for history only.

If a fact appears in two active files and conflicts, update the wrong file
immediately. Each type of truth should have one owner.
