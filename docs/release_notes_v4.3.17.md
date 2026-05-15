# ARIA v4.3.17 Release Notes

`v4.3.17` adds the first cross-modality Design Intelligence layer before v4.4.

## Added

- `DesignIntelligence` evaluates detected modalities, confirmed design, user
  intent, and available h5ad metadata before computation.
- The analysis-plan checkpoint now shows recommended, optional, unsupported, and
  warning items before the user launches compute.
- scRNA design intelligence detects condition, replicate, cell grouping,
  covariates, focused cell groups, velocity feasibility, LIANA feasibility,
  trajectory feasibility, and pseudobulk suitability.
- Initial bulk RNA, chromatin, Hi-C, and integration feasibility profiles provide
  modality-wide recommendations and limitations for current and roadmap agents.
- scRNAAgent honors Design Intelligence blocks for analyses such as LIANA or
  PAGA/DPT when they are not supported by the selected design.

## Validation

- `python -m compileall -q aria`
- `python -m pytest -q tests/test_pytest_smoke.py`
