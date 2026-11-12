# ARIA v4.3.16 Release Notes

`v4.3.16` is a pre-v4.4 report-integrity and scRNA focus patch.

## Fixes

- scRNA cell focusing no longer expands generic `neuron(s)` mentions into
  hippocampal neuronal subclasses. Focus is now limited to explicit cell-type
  matches and conservative aliases such as `oligodendrocyte`, `OPC`, and
  `microglía`.
- Integrated scRNA interpretation now uses a concise biological question instead
  of copying long prompt instructions into the report.
- LIANA/CellPhone p-values reported as zero are treated as unavailable and
  rendered as `—` instead of a misleading `0`.
- scRNAAgent now logs additional durable decisions for cell focusing, trusted
  obs grouping reuse, pseudobulk design, trajectory grouping, and LIANA grouping.

## Validation

- `python -m compileall -q aria`
- `python -m pytest -q tests/test_pytest_smoke.py`
