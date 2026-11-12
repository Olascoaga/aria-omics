# ARIA v4.3.15 Release Notes

`v4.3.15` is a pre-v4.4 scRNA efficiency patch.

## Fixes

- scRNA runs now materialize a focused h5ad before QC when the user explicitly
  asks to focus on one or a few trusted `obs` cell groups, such as `OPC` and
  `Oligo` under `obs['subclass']`.
- Downstream QC, clustering/annotation reuse, pseudobulk DE, trajectory, LIANA,
  figures, and reporting operate on the focused h5ad instead of all unrelated
  cell types.
- The focus step skips automatically when no specific focus is detected or when
  the detected focus covers every available group.

## Validation

- `python -m compileall -q aria`
- `python -m pytest -q tests/test_pytest_smoke.py`
