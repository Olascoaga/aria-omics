---
status: active
source_of_truth_for: project_state
last_updated: 2026-05-14
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
- Last code commit: `d3de169` (`Remove dataset-specific narrative guardrails`)
- Tags:
  - `v4.3.12` -> `3a0c40e`
  - `v4.3.12.post1` -> `805e0b2`

Do not move existing tags. Use a new patch tag if one is ever needed.

## Current Release Boundary

The 4.3 line is closed. `v4.3.12` is the stable baseline for bulk RNA and
scRNA workflows. `v4.3.12.post1` marks report-fidelity fixes; `d3de169` removed
dataset-specific runtime narrative guardrails after the post1 tag.

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

After the final code commit `d3de169`:

- `python -m compileall -q aria` -> pass
- `python -m pytest -q` -> 29 passed
- `python tests/test_narrative_agent.py` -> 23 passed
- `python tests/test_bulk_rna.py` -> 30 passed
- Runtime search over `aria/` and the new smoke tests found no hardcoded
  validation-gene guardrails.

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

## Next Milestone

Possible next milestone: `v4.4 scATAC`, but only after Samael explicitly asks.

Starting point if approved:

- Input: `/home/medusa/Samael/Erosion/data_inputs/muon_processed/hc11_paired.h5mu`
- Create `aria-chromatin-env` with a lockfile first.
- Validate `chromatin_qc.py` standalone before adding new scripts.
- Missing pieces: `chromatin_lsi_clustering.py`, `chromatin_diffacc.py`,
  `chromatin_motifs.py`, `_narrative_chromatin.py`, NarrativeAgent wiring.
