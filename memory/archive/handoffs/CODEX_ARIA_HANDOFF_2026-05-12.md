---
status: archived
source_of_truth_for: none
superseded_by: ../../PROJECT_STATE.md
last_updated: 2026-05-14
---

Historical record only. Not source of truth. See `../../PROJECT_STATE.md` and
`../../NEXT_SESSION.md`.

# ARIA Handoff - 2026-05-14 Final 4.3 Closeout

Repo: `/home/medusa/Samael/ARIA`

Status: **4.3 line closed**. User wrote "cerramos 4.13.2"; treat that as
`v4.3.12`. Current `main` is pushed at `2e0eb39` (documentation-only final
closeout). Last code change is:

- `d3de169 Remove dataset-specific narrative guardrails`
- `2e0eb39 Document final v4.3.12 closeout` is the current HEAD
- `805e0b2` is tagged `v4.3.12.post1`
- `3a0c40e` is tagged `v4.3.12`

Do not move existing tags. If a future formal tag is needed, use a new patch
tag rather than rewriting `v4.3.12` or `v4.3.12.post1`.

## Human Context

Samael prefers Spanish for planning/status and direct, concrete engineering
diagnosis. Keep code/comments/docs in English unless there is a reason not to.
He wants closure before new features. Next session should not start scATAC
unless he explicitly asks; first verify the repo is clean and the 4.3 closeout
state is understood.

## Non-Negotiables

- LLM proposes; code guarantees scientific invariants.
- No silent fake science: mocks require explicit dev/test opt-in.
- Resume logic must validate real files and matching parameters/manifest.
- Missing results stay missing. Reports must propagate warnings and
  empty-result markers instead of inventing biology.
- Methodology and warnings must be auditable.
- Runtime code and prompts must remain dataset-agnostic. No hardcoded genes,
  perturbations, validation datasets, or analysis-specific rescue phrases in
  narrative guardrails.

## What Closed In v4.3.12

Core closeout:

- Processed `.h5ad` path stabilized for the 40-donor GSE278576 hippocampus
  dataset.
- DataAudit infers condition, replicate, groupby, and covariates from `obs`.
- `rna_qc.py` handles processed h5ads using existing obs QC metrics.
- ARIA filters stale intermediate h5ads when real source h5ads exist.
- Existing `obs['subclass']` labels can be reused; Leiden and CellTypist are
  skipped when a trusted groupby column is provided.
- Pseudobulk donor-level DE is the main inferential layer for condition
  contrasts.
- scRNA reports export figures and TSV supplements.

Narrative/report fixes:

- Reports now say "reused obs['subclass']" rather than falsely reporting
  Leiden clustering when pre-existing annotations are reused.
- Typo `subclasss` fixed.
- PAGA/DPT language says exploratory manifold connectivity, not proof of
  active differentiation.
- Per-cluster Wilcoxon marker discovery timeout is surfaced in executive
  summary/body instead of being silently omitted.
- Report HTML uses `aria.__version__`, not `v0.2`.
- Methods and decisions blocks are HTML-escaped so `<10 counts` does not
  disappear.
- Raw Python error dicts are translated into human-readable report text.
- LIANA table shows rank and metric value explicitly; generic "Score" column
  removed to avoid misleading all-zero displays.
- Bulk RNA interpretation gets a generic anti-causality guardrail.
- Final fix removed dataset-specific guardrails and examples from runtime code
  (`BMAL1`, `REV-ERB`, circadian/pluripotency references no longer appear in
  `aria/` runtime code or the new smoke tests).

## Latest Validation

Final closeout validation after `d3de169`:

- `python -m compileall -q aria` -> pass
- `python -m pytest -q` -> 29 passed
- `python tests/test_narrative_agent.py` -> 23 passed
- `python tests/test_bulk_rna.py` -> 30 passed
- `rg` over `aria` and `tests/test_pytest_smoke.py` found no runtime mentions
  of `BMAL1`, `REV-ERB`, `REVERB`, `circadian`, `pluripotency`, or
  `transcriptional repressor`.

Known non-fatal warnings:

- LiteLLM cannot fetch the remote cost map because network is restricted.
- Matplotlib may use a temp cache because `/home/medusa/.config/matplotlib`
  is not writable.
- Optional pathway plotting dependencies such as `blitzgsea` may be absent.

## Latest Reviewed Report

`/home/medusa/.aria/reports/aria_20260514_143352_oligodendrocytes_opcs_microglia_-009/report.html`

Observed outputs:

- 295,033 starting cells -> 242,405 retained after QC.
- 18 subclass groups reused from `obs['subclass']`.
- Pseudobulk DE: 79 analyzable group x comparison blocks, 57 significant.
- Top DE blocks: Oligo/Astro/OPC `80-100_vs_40-59`.
- Pathway ORA: 38/57 DE blocks with enrichment.
- LIANA: 50 non-autocrine interactions.
- PAGA/DPT generated; RNA velocity skipped because no spliced/unspliced layers.

The report was reviewed before final report-fidelity hotfixes, so newly
generated reports should be cleaner than this specific HTML.

## Repo Hygiene

- `ARIA_CONTEXT.md` and `CLAUDE.md` remain private/ignored.
- `audit.txt` and `pathways_per_cluster.csv` are ignored local artifacts.
- Worktree should be clean at handoff unless the user asks for further docs.

## Next Work If/When User Restarts

Recommended next milestone: `v4.4 scATAC`, but only after a fresh explicit
request.

Likely starting point:

- Input: `/home/medusa/Samael/Erosion/data_inputs/muon_processed/hc11_paired.h5mu`
- Create `aria-chromatin-env` with a lockfile first.
- Validate `chromatin_qc.py` standalone before adding new scripts.
- Missing pieces: `chromatin_lsi_clustering.py`, `chromatin_diffacc.py`,
  `chromatin_motifs.py`, `_narrative_chromatin.py`, and NarrativeAgent wiring.
