# CLAUDE.md — ARIA Operating Manual

Last revised: 2026-05-12 evening handoff, TUI rerun in progress.

This file is the handoff manual for Claude or any other coding agent working on
ARIA. Keep it current after every meaningful architecture, testing, or workflow
change.

## Human Context

Samael Olascoaga is a bioinformatics researcher in Mexico.

Working style:

- Use Spanish for discussion, planning, and status updates unless he switches.
- Keep code, comments, docstrings, and technical docs in English.
- Be direct and concrete. Diagnose before coding.
- Avoid marketing language and vague optimism.
- When something is broken, identify the root cause and fix it.
- If there are options, present them briefly and recommend one.

## Product Definition

ARIA means Agentic Research Intelligence for -omics Analysis.

ARIA is a local-first multi-agent system that takes a biological question and
omics data, confirms the experiment design, runs modality-specific analyses,
records methodology decisions, and generates a publication-style HTML report.

The differentiator is the semantic layer:

- biological intent parsing,
- experimental design confirmation,
- parameter decisions with rationale,
- guarded interpretation,
- structured reports grounded in real output files.

The pipeline stack alone is not the product. The product is the supervised,
auditable, biologically aware orchestration around the pipeline.

## Core Engineering Principles

### 1. LLM Proposes, Code Guarantees

Never trust model prose for scientific invariants.

Examples of things code must guarantee:

- contrast generation,
- sample-to-group mapping,
- threshold application,
- missing-result handling,
- pathway input gene lists,
- report input schema,
- production-vs-dev mock behavior.

### 2. No Silent Fake Science

Mock outputs are development tools, not production behavior.

Mocks are allowed only when one of these is explicit:

- params include `allow_mock=true`,
- params include `allow_mocks=true`,
- environment has `ARIA_ALLOW_MOCKS=1`,
- environment has `ARIA_DEV_MODE=true`.

Otherwise missing tools/dependencies must return structured errors.

This rule currently applies to the hardened core paths:

- `pydeseq2` in bulk DE,
- `gseapy` pathway fallback,
- `fastp`,
- `STAR`,
- `featureCounts`,
- Hi-C cooler/hic-straw paths,
- WNN,
- MOFA+,
- peak-to-gene.

### 3. Missing Results Stay Missing

If enrichment, plots, annotation, or integration fails, pass the exact warning
forward. Do not let NarrativeAgent or any LLM invent plausible biology.

### 4. Resume Logic Is File-Validity Based

Heavy steps must check real output validity:

- BAM validity,
- gzip integrity,
- mtime comparisons,
- expected columns,
- parseable JSON,
- non-empty output files.

Do not rely on in-memory flags.

### 5. Methodology Must Be Auditable

Every important user or code decision should be visible later:

- thresholds,
- design formula,
- normalization,
- filtering,
- covariates,
- batch decisions,
- quality warnings,
- skipped analyses and why.

## Current Version

Current repo state: v4.3.12-dev.

Packaged/released version strings are still v4.3.11 in:

- `aria/__init__.py`
- `aria/llm/__init__.py`
- `setup.py`
- `aria/tui.py`
- `install.sh`
- `README.md`
- `docs/INSTALLATION.md`

Do not bump/tag until the hippocampus TUI rerun validates the new h5ad obs
design path end-to-end.

## What Was Done In v4.3.11

### Bulk RNA Regression Fixed

`tests/test_bulk_rna.py` was failing 3 checks:

- full synthetic bulk pipeline,
- multi-contrast bulk pipeline,
- DESeq2 design factor propagation.

Root causes:

- sample QC could mark all synthetic samples as outliers and remove them before
  DE, leaving no valid groups;
- pyDESeq2 API changed from `design_factors=...` to `design="~ factor"`;
- mocks were being used implicitly when dependencies were absent.

Fixes:

- Added `_prune_outliers_for_design()` in `aria/scripts/rna_bulk_de.py`.
- Preserved minimum two-replicate group structure before dropping QC outliers.
- Added pyDESeq2 compatibility path for both new and old APIs.
- Added explicit mock gating via `mocks_allowed()`.

Validation:

- `python tests/test_bulk_rna.py` now passes 30/30.

### Pytest Collection Fixed

Most tests were legacy executable scripts with top-level code and `sys.exit()`.
Running `pytest` imported them, executed them during collection, and aborted.

Fixes:

- Added `tests/conftest.py` to ignore legacy scripts during pytest collection.
- Added `tests/test_pytest_smoke.py` with subprocess wrappers.

Validation:

- `python -m pytest -q` passes.

### Production Mock Hardening

Added `mocks_allowed(params)` to `aria/scripts/_base.py`.

Updated scripts so production runs fail loudly instead of silently generating
fake output when dependencies are missing:

- `aria/scripts/rna_bulk_de.py`
- `aria/scripts/rna_fastq_qc.py`
- `aria/scripts/rna_align.py`
- `aria/scripts/rna_quantify.py`
- `aria/scripts/hic_qc_and_balance.py`
- `aria/scripts/hic_topology.py`
- `aria/scripts/integration_wnn.py`
- `aria/scripts/integration_mofa.py`
- `aria/scripts/integration_peak2gene.py`

Tests that intentionally need mocks pass `allow_mock=true`.

### Memory Updated

The previous `ARIA_CONTEXT.md` described v3.9 and v4.0 roadmap. It was no
longer true. It has been replaced with current v4.3.11 memory.

## What Was Done After v4.3.11

### h5ad obs Design Repair

The TUI run at
`/home/medusa/.aria/reports/aria_20260512_131402_oligodendrocytes_opcs_microglia_-e5f/report.html`
completed the scRNA workflow but did not run young-vs-old pseudobulk DE. The
report explicitly said age-stratified differential analyses were absent.

Root cause:

- `DataAuditAgent` did not classify `.h5ad` as scRNA or inspect `adata.obs`.
- `DesignAgent` accepted inferred design only when it came from GEO metadata.
- `scRNAAgent` pseudobulk always injected a condition column from
  `design.groups`, even when the h5ad already had condition/replicate/cell-type
  columns.

Fixes:

- `aria/agents/data_audit_agent.py` now recognizes `.h5ad` and infers design
  hints from `obs`: condition, biological replicate, groupby/cell-type,
  covariates, groups, and comparisons.
- `aria/agents/design_agent.py` now preserves inferred design from either GEO or
  h5ad obs, seeds the main factor checkpoint, and carries
  `design["pseudobulk"]` downstream.
- `aria/agents/scrna_agent.py` now runs pseudobulk directly on native obs
  columns when `design["pseudobulk"]["from_obs"]` is true.
- `tests/test_pytest_smoke.py` covers h5ad obs design inference and direct
  obs-based pseudobulk parameter wiring.

Validation:

- `python -m compileall -q aria`
- `python -m pytest -q` passes 5 tests.

### Processed h5ad QC Repair

The follow-up TUI report at
`/home/medusa/.aria/reports/aria_20260512_165405_oligodendrocytes_opcs_microglia_-0d2/report.html`
proved h5ad obs design inference reached CP1/CP2, but scRNA failed immediately:
the report said zero cells passed QC.

Root cause:

- The input h5ad was a processed Seurat-style object with 295,033 cells, 2,000
  features, scaled/log-normalized `X`, and real QC metrics in `obs`
  (`nFeature_RNA`, `nCount_RNA`, `percent.mt`).
- `rna_qc.py` treated `X` as raw counts, recalculated Scanpy QC metrics, derived
  a destructive `min_genes=1998`, and filtered all cells.
- The script then returned success with `n_cells_after=0`, letting the report
  become the first visible failure point.

Fixes:

- `aria/scripts/rna_qc.py` now detects h5ads with existing obs QC metrics and
  filters cells using those metrics instead of recalculating from processed `X`.
- It skips Scrublet in that mode because Scrublet requires raw counts.
- It returns structured `NoCellsAfterQC` errors if any QC path leaves zero
  cells, rather than saving an empty h5ad as success.
- `tests/test_pytest_smoke.py` includes a processed-h5ad regression test with
  negative/scaled `X` and Seurat-style QC columns.

Latest validation:

- The follow-up h5ad TUI report completed analytically: QC retained
  242,405/295,033 cells, pseudobulk DE ran across 18 subclasses, pathway ORA,
  LIANA communication, and PAGA/DPT trajectory outputs were present.
- The remaining defects were in NarrativeAgent: the executive summary
  contradicted the body by saying completed analyses were absent, scRNA prose
  was too shallow, and the report `tables/` directory stayed empty for
  scRNA-only runs.
- NarrativeAgent now uses deterministic scRNA summaries when structured
  pseudobulk/ORA/LIANA/trajectory outputs exist, includes those outputs in the
  LLM concrete-results block, adds an Integrated Interpretation section, and
  exports scRNA supplementary TSVs for QC, cell types, markers, pseudobulk DE,
  pathway enrichment, LIANA interactions, PAGA connections, and DPT
  pseudotime.
- A hardcode audit removed dataset-specific narrative wording. Remaining
  age/young/old strings in production code are either generic design heuristics
  or UI examples; dataset-specific values are confined to tests and handoff
  documentation.
- Leave pre-existing untracked `audit.txt` and `pathways_per_cluster.csv`
  untouched unless Samael explicitly decides what to do with them.

## Architecture Map

### Agents

- `orchestrator_agent.py`
  Central flow: parse question, run DataAudit, hand off to DesignAgent,
  confirm plan, run AuditAgent, dispatch modality agents, call NarrativeAgent.

- `data_audit_agent.py`
  Scans input directory, classifies modalities, infers organism/genome, builds
  experiment context, triggers checkpoint 1.

- `design_agent.py`
  Interactive design state machine: groups, organism, factor, batch,
  pseudoreplication, confirmation.

- `audit_agent.py`
  Pre-dispatch quality linter. Currently strongest for bulk RNA; checks
  replicate correlation, PCA batch dominance, STAR alignment rates.

- `bulk_rna_agent.py`
  Bulk RNA-seq orchestration, including raw FASTQ preprocessing when needed.

- `scrna_agent.py`
  scRNA path: QC, concat, Harmony, clustering, CellTypist, DE, pseudobulk,
  trajectory, cell communication.

- `chromatin_agent.py`
  ATAC/ChIP/CUT&RUN/CUT&TAG scaffold. Next major target.

- `genome_arch_agent.py`
  Hi-C / Micro-C orchestration: QC/balance, compartments, TADs, loops.

- `integration_agent.py`
  WNN/MOFA+/peak-to-gene integration scaffold. Defer production hardening until
  standalone modalities are stronger.

- `narrative_agent.py`
  HTML report and methods. Must stay grounded in structured outputs.

### Infrastructure

- `environment_manager.py`
  Runs scripts in Conda stacks via JSON files. Archives failed runs.

- `memory.py`
  SQLite persistent memory: wings, halls, rooms, findings, tunnels, decisions.
  Thread-hardened with lock + WAL.

- `message_bus.py`
  In-process pub/sub with checkpoint log. Thread-hardened: dispatches outside
  lock.

- `provider.py`
  LiteLLM abstraction with tier routing and fallback.

## Current Validation Commands

Run these after changes:

```bash
python -m compileall -q aria
python tests/test_bulk_rna.py
python -m pytest -q
```

Current expected result:

- compileall: no output, exit 0
- bulk RNA: 30 passed / 0 failed
- pytest: 14 passed

Known non-fatal warnings in this environment:

- LiteLLM may warn that remote model cost map cannot be fetched because network
  is restricted.
- Matplotlib may create a temporary cache under `/tmp` because
  `/home/medusa/.config/matplotlib` is not writable.
- `blitzgsea` and optional pathway plotting dependencies may be absent.

## Known Limitations

### Tests

The test suite is not fully modern.

Current state:

- legacy scripts still exist,
- pytest wrappers keep collection stable,
- broad CI-quality coverage is not done.

Next test work:

- convert legacy scripts into pytest modules,
- remove top-level execution from tests,
- replace manual counters with asserts,
- add markers: `slow`, `external_data`, `requires_conda_stack`,
- add small deterministic fixtures.

### Bulk RNA

Now passes synthetic validation. Still needs a real production dependency check
for RNA stack installation and possibly a small fixture that runs with real
`pydeseq2` when available.

### scRNA

The code now has a specific repair for `.h5ad.obs` experimental design metadata
feeding pseudobulk DE. Still needs the real hippocampus TUI rerun and cleaner
pytest coverage for report shape normalization.

Current v4.3.12 hardening also includes large-dataset scRNA guards:

- QC, concat, Harmony integration, Leiden resolution advice, and clustering can
  resume from file-backed summaries.
- Resume for QC, concat, clustering, and Harmony is valid only when the cached
  summary has matching parameter/manifest signatures; legacy or mismatched
  caches are rerun.
- Harmony is skipped above a configurable cell limit to avoid OOM and the skip
  reason is propagated to findings/reports.
- Leiden resolution advice and clustering use a deterministic sketch when cell
  counts exceed the configured limit, and reports flag sketch-based results.
- If CellTypist cannot produce usable labels, ARIA applies conservative
  marker-panel labels where possible and marks them as `marker_fallback` so the
  report treats them as curation targets, not definitive identities.

### Chromatin/scATAC

Next real milestone. Current roadmap target: v4.4.

Need:

- scATAC QC,
- LSI,
- clustering,
- differential accessibility,
- motif enrichment,
- report module,
- tests with small fixtures.

### Integration

Do not over-invest yet. WNN/MOFA+/peak-to-gene should wait until scATAC/bulk
ATAC/Hi-C standalone workflows are reliable.

## Priority Roadmap

### Immediate

1. Rerun the hippocampus TUI case and inspect whether pseudobulk DE executes.
2. If the rerun passes, package/tag the h5ad design repair as v4.3.12.
3. Decide whether `audit.txt` and `pathways_per_cluster.csv` should be ignored,
   committed as artifacts, or removed by user instruction.

### v4.4

scATAC end-to-end:

- input detection hardening,
- fragment/peak matrix QC,
- LSI + clustering,
- marker peaks,
- differential accessibility,
- motif enrichment,
- NarrativeAgent chromatin/scATAC section,
- focused tests.

### v4.5

Bulk ATAC end-to-end:

- peak count matrix,
- DA via DESeq2-style count model,
- FRiP/TSS QC,
- pathway/motif summary,
- report integration.

### v4.6

IntegrationAgent:

- WNN only after scRNA + scATAC are standalone stable,
- MOFA+ after at least two robust modalities,
- peak-to-gene after ATAC + RNA inputs are reliable,
- DebateCouncil only for interpretation, not data fabrication.

## How To Work In This Repo

Before coding:

```bash
git status --short
rg "target_symbol_or_error"
```

When editing:

- Keep changes scoped.
- Prefer existing patterns.
- Use `apply_patch`.
- Do not revert user changes.
- Do not add broad abstractions unless they remove real duplication or risk.

When working on scripts:

- Scripts in `aria/scripts/` are subprocess-only.
- They should take params dict and return JSON-serializable dicts.
- They must not access the message bus.
- They must return structured errors.
- If a mock is needed, gate it through `mocks_allowed(params)`.

When working on reports:

- Preserve missing-result warnings.
- Do not hide failures.
- Do not ask LLMs to infer from absent data.
- Keep figures both embedded and written as files when possible.

## Files To Read First By Task

Bulk RNA:

- `aria/agents/bulk_rna_agent.py`
- `aria/scripts/rna_bulk_de.py`
- `tests/test_bulk_rna.py`
- `aria/agents/narrative_agent.py`

scRNA:

- `aria/agents/data_audit_agent.py`
- `aria/agents/design_agent.py`
- `aria/agents/scrna_agent.py`
- `aria/agents/_narrative_scrna.py`
- `tests/test_scrna_e2e.py`
- relevant `aria/scripts/rna_*.py`

Chromatin/scATAC:

- `aria/agents/chromatin_agent.py`
- `aria/scripts/chromatin_qc.py`
- `aria/scripts/chromatin_peaks.py`
- `tests/test_chromatin_agent.py`

Hi-C:

- `aria/agents/genome_arch_agent.py`
- `aria/scripts/hic_qc_and_balance.py`
- `aria/scripts/hic_topology.py`
- `tests/test_genome_arch_agent.py`

Infrastructure:

- `aria/utils/environment_manager.py`
- `aria/memory/memory.py`
- `aria/bus/message_bus.py`
- `aria/llm/provider.py`

## Final Reminder

ARIA should be scientifically conservative. It is better to stop with a precise
dependency/data/design error than to produce a polished but invalid report.
