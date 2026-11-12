# ARIA v4.3.18 Release Notes

`v4.3.18` tightens Design Intelligence control and cell-focus specificity.

## Fixes

- The analysis-plan checkpoint now asks whether to run only recommended analyses
  or to add optional supported analyses, making Design Intelligence suggestions
  actionable instead of ambiguous.
- Optional supported analyses are propagated through `exp_context`; scRNA can
  add optional LIANA or PAGA/DPT only when the user selects that mode and Design
  Intelligence does not mark them unsupported.
- scRNA cell focusing now uses only explicit focus/restriction clauses such as
  `focus only on astrocytes`, `obs['subclass'] == 'Microglia'`, `solo`, or
  `únicamente`.
- Mentions outside explicit focus clauses, including negative instructions such
  as “do not run oligodendrocyte trajectory,” no longer expand the focused h5ad.

## Validation

- `python -m compileall -q aria`
- `python -m pytest -q tests/test_pytest_smoke.py`
