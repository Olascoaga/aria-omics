# ARIA Project Memory

Last revised: 2026-05-12 evening handoff, TUI rerun in progress.

This is the compact memory to paste into a new AI/code-agent session before
working on ARIA. For the full operating manual, read `CLAUDE.md`.

## Collaborator

Samael Olascoaga is a bioinformatics researcher in Mexico. Spanish is preferred
for planning and project discussion; code, docstrings, comments, and technical
artifacts are usually English. He values direct diagnosis, concrete fixes, and
low-fluff communication.

## What ARIA Is

ARIA is a local-first multi-agent bioinformatics system for omics analysis. The
goal is not merely to run pipelines, but to add a semantic layer:

- understand the biological question,
- detect data modalities,
- confirm experimental design interactively,
- run the correct modality agents,
- record methodology decisions,
- produce a publication-style HTML report grounded in real outputs.

Core rule: LLM proposes, code guarantees.

## Architecture Snapshot

- `aria/agents/orchestrator_agent.py`: checkpoint flow and agent dispatch.
- `aria/agents/data_audit_agent.py`: scans data, detects modalities, builds
  experiment context.
- `aria/agents/design_agent.py`: interactive design confirmation.
- `aria/agents/audit_agent.py`: pre-dispatch quality linter.
- `aria/agents/bulk_rna_agent.py`: bulk RNA-seq orchestration.
- `aria/agents/scrna_agent.py`: scRNA orchestration.
- `aria/agents/chromatin_agent.py`: ATAC/ChIP/CUT&RUN/CUT&TAG scaffold.
- `aria/agents/genome_arch_agent.py`: Hi-C / 3D genome orchestration.
- `aria/agents/integration_agent.py`: WNN/MOFA+/peak-to-gene scaffold.
- `aria/agents/narrative_agent.py`: HTML report generation.
- `aria/scripts/`: subprocess-only analytical scripts; no bus access.
- `aria/utils/environment_manager.py`: Conda stack isolation + JSON IPC.
- `aria/memory/memory.py`: SQLite memory for decisions/findings.
- `aria/bus/message_bus.py`: in-process pub/sub and checkpoint log.

## Non-Negotiable Design Rules

1. LLM proposes; code guarantees.
   Contrasts, thresholds, sample mappings, missing-result warnings, and report
   inputs must be enforced deterministically.

2. No silent fake science.
   Mock outputs are allowed only with explicit development opt-in:
   `allow_mock=true`, `allow_mocks=true`, `ARIA_ALLOW_MOCKS=1`, or
   `ARIA_DEV_MODE=true`. Production analysis must fail loudly when required
   tools or dependencies are missing.

3. Resume logic validates files.
   Heavy steps skip only when outputs are present and valid, not because a flag
   says a step ran.

4. Missing results remain missing.
   Feed NarrativeAgent explicit warnings and empty-result markers. Do not let an
   LLM invent biology or causes for missing enrichment/plots.

5. Methodology is auditable.
   User decisions, thresholds, design formulas, normalization choices, filters,
   and warnings belong in memory and report methods.

## Current State: v4.3.12-dev

v4.3.11 fixed the bulk RNA and production-run hardening layer:

- Bulk RNA synthetic validation repaired: `tests/test_bulk_rna.py` passes 30/30.
- Pytest collection repaired with wrappers: `python -m pytest -q` passes.
- Bulk RNA outlier pruning now preserves minimum replicate structure before
  removing samples from DE.
- pyDESeq2 compatibility supports both newer `design="~ factor"` and older
  `design_factors=...` APIs.
- Core production mocks now require explicit dev/test opt-in.
- Version strings aligned to v4.3.11 across package, LLM module, setup, TUI,
  installer, README, and installation docs.
- The old Claude memory was replaced because it described v3.9/v4.0 and was no
  longer a source of truth.

The 2026-05-12 TUI run at
`/home/medusa/.aria/reports/aria_20260512_131402_oligodendrocytes_opcs_microglia_-e5f/report.html`
completed QC/Harmony/clustering/CellTypist/trajectory/LIANA, but still did not
run age-stratified pseudobulk differential expression. Root cause: ARIA treated
`.h5ad` as data only, not as a carrier of experimental design metadata. CP2 used
manual filename/sample mapping while scRNA pseudobulk required injected
condition columns, so existing `obs` fields such as age group, donor, sex, and
cell subclass were ignored.

New repair in progress:

- `DataAuditAgent` now classifies `.h5ad` as scRNA and inspects `adata.obs`
  before CP1.
- It infers condition, biological replicate, cell-type/groupby, covariates, and
  comparisons from common obs columns such as `age_group`, `donor_id`,
  `subclass`, and `Gender`.
- `DesignAgent` now accepts inferred designs from either GEO metadata or
  `.h5ad.obs`, not only GEO. It preserves the inferred main factor and passes a
  `design["pseudobulk"]` block downstream.
- `scRNAAgent` can run pseudobulk directly from native obs columns when
  `design["pseudobulk"]["from_obs"]` is true, skipping condition injection.
- Added pytest coverage for `.h5ad.obs` design inference and for direct
  obs-based pseudobulk parameterization.
- Follow-up from the 2026-05-12 16:54 TUI report:
  `/home/medusa/.aria/reports/aria_20260512_165405_oligodendrocytes_opcs_microglia_-0d2/report.html`
  showed CP1/CP2 h5ad obs design was working, but QC filtered all cells. Root
  cause: `rna_qc.py` reapplied X-based Scanpy filters to a preprocessed Seurat
  h5ad whose `X` is scaled/log-normalized and whose real QC metrics live in
  `obs` (`nFeature_RNA`, `nCount_RNA`, `percent.mt`). `n_genes_by_counts`
  collapsed to a flat 1998-feature value, leading to a destructive second
  filter. `rna_qc.py` now uses existing h5ad obs QC metrics when present,
  skips invalid Scrublet-on-processed-X, and returns `NoCellsAfterQC` instead
  of success if filtering ever leaves zero cells.

Validated locally after this repair:

- `python -m compileall -q aria`
- `python -m pytest -q` (5 passed, includes the legacy bulk wrapper)

Evening handoff, 2026-05-12:

- Samael started a new TUI rerun after the processed-h5ad QC repair, then
  reported the full h5ad run at
  `/home/medusa/.aria/reports/aria_20260512_191451_oligodendrocytes_opcs_microglia_-bb7/report.html`.
- That report succeeded analytically: QC retained 242,405/295,033 cells,
  pseudobulk DE ran across 18 subclasses, pathway ORA, LIANA communication,
  and PAGA/DPT trajectory outputs were present.
- Remaining defect from that report was narrative quality. The executive
  summary falsely said DE, trajectory, and L-R outputs were absent even though
  the body contained them, and the scRNA findings section compressed the rich
  result set into a shallow paragraph.
- NarrativeAgent was repaired so scRNA executive summaries are deterministic
  when structured pseudobulk/ORA/LIANA/trajectory outputs exist, the LLM
  concrete-results block includes those analyses, the scRNA findings prose is
  more detailed, and a final Integrated Interpretation section combines DE,
  pathways, communication, and trajectory.
- NarrativeAgent also now exports scRNA supplementary TSV tables before HTML
  rendering. The old staging path only copied bulk RNA contrast tables, so
  scRNA-only reports produced an empty `tables/` directory. New scRNA reports
  can include TSVs for QC per sample, cell-type annotations, cluster markers,
  pseudobulk DE summaries, pseudobulk DE genes, pathway enrichment, LIANA
  interactions, PAGA connections, and DPT pseudotime by group, with links in
  the scRNA report section.
- A pre-commit hardcode audit removed dataset-specific narrative wording.
  Remaining age/young/old strings in production code are generic design
  heuristics or UI examples; hippocampus-specific values are confined to tests
  and handoff documentation.
- ARIA now auto-loads API keys from `~/.aria/.env` when `LLMProvider` is
  constructed. The loader accepts both `KEY=value` and `export KEY=value`,
  preserves variables already exported in the terminal, and can be redirected
  with `ARIA_ENV_FILE`.
- Current working tree intentionally has uncommitted v4.3.12-dev changes in:
  `ARIA_CONTEXT.md`, `CLAUDE.md`, `aria/agents/data_audit_agent.py`,
  `aria/agents/design_agent.py`, `aria/agents/narrative_agent.py`,
  `aria/agents/_narrative_scrna.py`, `aria/agents/scrna_agent.py`,
  `aria/scripts/rna_qc.py`, `aria/utils/env_loader.py`,
  `tests/test_narrative_agent.py`, and `tests/test_pytest_smoke.py`.
- Pre-existing untracked files remain untouched: `audit.txt` and
  `pathways_per_cluster.csv`.

## Known Limitations

- Most tests are still legacy executable scripts rather than idiomatic pytest
  modules. They are preserved and wrapped for now.
- scRNA appears actively hardened and documented, but broad pytest coverage is
  still shallow.
- The `.h5ad.obs` design inference is heuristic. It covers common columns but
  still needs validation on the hippocampus TUI rerun and more real public h5ads.
- Chromatin/scATAC is the next real standalone modality to harden.
- IntegrationAgent should remain deferred until standalone scATAC/bulk ATAC/Hi-C
  paths are reliable.
- Optional plotting/enrichment paths can still degrade when external packages
  such as `gseapy` or `blitzgsea` are missing; warnings must stay explicit.

## Priority Roadmap

1. Rerun or regenerate the hippocampus TUI report after the NarrativeAgent
   repair. Expected new behavior: the executive summary no longer claims that
   completed DE/trajectory/L-R analyses are missing, scRNA findings include
   detailed pseudobulk/pathway/communication/trajectory prose, and the report
   ends the findings section with an Integrated Interpretation. The report
   directory should also contain non-empty scRNA TSVs under `tables/`.

2. If the rerun succeeds, tag this as v4.3.12 or fold it into v4.4 prework.

3. v4.4 scATAC end-to-end:
   QC, LSI, clustering, differential accessibility, motif enrichment, and
   NarrativeAgent report integration.

4. Test modernization:
   Convert legacy test scripts into pytest modules with fixtures, asserts,
   markers for slow/external-data tests, and no `sys.exit()` during import.

5. Dependency hardening:
   Make env YAMLs and installer match every advertised production workflow.

6. Narrative risk reduction:
   Keep reports driven by structured outputs and explicit warnings only.

7. IntegrationAgent productionization:
   WNN/MOFA+/peak-to-gene after standalone modalities are validated.
