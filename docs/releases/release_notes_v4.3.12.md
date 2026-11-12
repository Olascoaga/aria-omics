# ARIA v4.3.12 Release Notes

Date: 2026-05-14

## Scope

`v4.3.12` closes the current bulk RNA + scRNA hardening cycle. The release is
focused on stability, report fidelity, processed `.h5ad` handling, and
large-dataset behavior. It intentionally does not add new modality features.

## Stabilized Paths

- Bulk RNA-seq count-matrix workflow with regression tests.
- scRNA single-sample and multi-sample workflows.
- Processed `.h5ad` pseudobulk workflows using design metadata from `obs`.
- Narrative reports with embedded figures and supplementary TSV exports.
- PAGA/DPT and LIANA reporting as beta analysis layers.

## scRNA Hippocampus Rerun

The reviewed GSE278576 hippocampus rerun used a processed 40-donor `.h5ad`:

- 295,033 starting cells.
- 242,405 cells retained after QC.
- Existing `obs['subclass']` annotations reused as the grouping layer.
- Donor-level pseudobulk DE used as the primary inferential evidence.
- 79 group x comparison pseudobulk blocks analyzable.
- 57 blocks with significant DE.
- 38/57 DE blocks with pathway enrichment.
- 50 non-autocrine LIANA interactions.
- PAGA/DPT trajectory summaries generated for exploratory manifold context.

## Scientific Policy

- Pseudobulk donor-level contrasts carry the main inferential weight for
  between-condition scRNA claims.
- Cell-level per-cluster DE is optional on atlas-scale inputs. A timeout in
  `rna_de_per_cluster.py` does not invalidate a run when pseudobulk, pathway,
  communication, and trajectory outputs complete.
- PAGA/DPT is reported as manifold ordering/connectivity, not proof of active
  differentiation. Velocity or time-course data are required for stronger
  lineage claims.
- Reused input annotations are reported as reused `obs` groupings, not as newly
  inferred Leiden clusters.
- Partial analysis failures are surfaced explicitly in the report. On
  atlas-scale scRNA inputs, per-cluster Wilcoxon marker discovery may time out;
  this is reported as unavailable rather than silently omitted.
- Bulk RNA integrated interpretations must avoid unsupported causal language.
  Differential expression and pathway enrichment can support hypotheses, but
  not direct transcriptional causality without orthogonal validation.
- Narrative guardrails must be dataset-agnostic. No gene, perturbation, or
  validation dataset should be hardcoded into runtime prompts or report
  post-processing.

## Final 4.3 Closeout Hotfixes

After the initial `v4.3.12` tag, two report-fidelity hotfixes closed the 4.3
branch:

- `805e0b2` / `v4.3.12.post1`: fixed report version display, escaped Methods
  and decisions HTML, rendered raw error dictionaries as human-readable
  messages, clarified LIANA rank metrics, and added a generic anti-causality
  guardrail for bulk RNA interpretation.
- `d3de169`: removed dataset-specific narrative guardrails and prompt examples
  from runtime code. The report fidelity rules are now generic and no longer
  reference any particular gene, perturbation, or validation experiment.

## Repository Hygiene

- Private agent memory files remain ignored.
- Local analysis scratch outputs such as `audit.txt` and
  `pathways_per_cluster.csv` are ignored.
- The public documentation and roadmap now mark `v4.3.12` as the stable
  baseline for the current cycle.

## Verification

Release closeout validation:

- `python -m compileall -q aria`
- `python -m pytest -q` -> 29 passed after final closeout
- `python tests/test_narrative_agent.py` -> 23 passed
- `python tests/test_bulk_rna.py` -> 30 passed
