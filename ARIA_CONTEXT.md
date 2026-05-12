# ARIA Project Memory

Last revised: v4.3.11 hardening pass.

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

## Current State: v4.3.11

What was just fixed:

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

Validated locally after the hardening pass:

- `python -m compileall -q aria`
- `python tests/test_bulk_rna.py`
- `python -m pytest -q`

## Known Limitations

- Most tests are still legacy executable scripts rather than idiomatic pytest
  modules. They are preserved and wrapped for now.
- scRNA appears actively hardened and documented, but broad pytest coverage is
  still shallow.
- Chromatin/scATAC is the next real standalone modality to harden.
- IntegrationAgent should remain deferred until standalone scATAC/bulk ATAC/Hi-C
  paths are reliable.
- Optional plotting/enrichment paths can still degrade when external packages
  such as `gseapy` or `blitzgsea` are missing; warnings must stay explicit.

## Priority Roadmap

1. v4.4 scATAC end-to-end:
   QC, LSI, clustering, differential accessibility, motif enrichment, and
   NarrativeAgent report integration.

2. Test modernization:
   Convert legacy test scripts into pytest modules with fixtures, asserts,
   markers for slow/external-data tests, and no `sys.exit()` during import.

3. Dependency hardening:
   Make env YAMLs and installer match every advertised production workflow.

4. Narrative risk reduction:
   Keep reports driven by structured outputs and explicit warnings only.

5. IntegrationAgent productionization:
   WNN/MOFA+/peak-to-gene after standalone modalities are validated.

