# CLAUDE.md — ARIA Operating Manual

Last revised: v4.3.11 hardening pass.

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

Current repo state: v4.3.11.

Version strings were aligned in:

- `aria/__init__.py`
- `aria/llm/__init__.py`
- `setup.py`
- `aria/tui.py`
- `install.sh`
- `README.md`
- `docs/INSTALLATION.md`

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
- pytest: 2 passed

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

The README claims substantial scRNA E2E validation. Code is organized and
appears actively hardened. Still needs cleaner pytest coverage and explicit
fixtures for report shape normalization.

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

1. Commit or package v4.3.11 hardening changes.
2. Run full local verification once more after any doc/manual edits.
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

