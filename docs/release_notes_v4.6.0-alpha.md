# ARIA v4.6.0-alpha Release Notes

`v4.6.0-alpha` opens the scATAC dispatch lane with explicit acknowledgement
rather than silent autonomous execution.

## scATAC Alpha Gate

- `ChromatinAgent` is marked alpha and dispatch-enabled for scATAC.
- The orchestrator validation gate now treats `scATAC` as `alpha` with
  `dispatch_enabled=True`.
- The readiness matrix keeps scATAC yellow (`requires_ack`), so CP3.5 explicit
  acknowledgement is required before dispatch.
- bulk ATAC, ChIP, CUT&RUN, and CUT&TAG remain scaffolded and blocked.

## Honesty Boundary

- scATAC alpha is not yet declared fully autonomous or publication-grade.
- A live orchestrator/TUI validation run remains required before promoting the
  lane beyond alpha.
- The version badge and release-note guard are synchronized from
  `aria.version.__version__`.

## Validation

- Focused chromatin dispatch/readiness/registry tests passed.
- Bulk RNA ADR-011 refactor tests passed.
- Version metadata smoke passed after adding these release notes and regenerating
  the README badge.
- Full smoke passed with writable Numba/Matplotlib caches: 100 passed / 4
  skipped.
- Report provenance now prefers RNA lockfiles for RNA/report tools, preventing
  the chromatin lock's `gseapy` version from shadowing the RNA lock.
