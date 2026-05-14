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

## Repository Hygiene

- Private agent memory files remain ignored.
- Local analysis scratch outputs such as `audit.txt` and
  `pathways_per_cluster.csv` are ignored.
- The public documentation and roadmap now mark `v4.3.12` as the stable
  baseline for the current cycle.

## Verification

Release closeout validation:

- `python -m compileall -q aria`
- `python -m pytest -q`
- `python tests/test_narrative_agent.py`

