# ARIA v4.3.19 Release Notes

`v4.3.19` closes four P0 audit findings from the v4.3 maintenance cycle.

## Fixes

- scRNA report text now derives the pseudobulk DE narrative header from the
  active condition column instead of hardcoding `Age-associated`.
- Chromatin QC basic-stat fallback is gated by explicit mock permission.
  Production runs with missing chromatin dependencies now fail loudly.
- LIANA cell-cell communication failures are no longer swallowed as fallback
  results. Only missing LIANA imports fall back to mean-expression scoring.
- Hi-C inspection scripts no longer import constants from an agent module.
  `RAM_ESTIMATES_GB` now lives in `aria.utils.hic_constants`.

## Validation

- `python -m compileall -q aria`
- `python -m pytest -q tests/test_pytest_smoke.py`
