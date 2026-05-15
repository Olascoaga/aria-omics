---
status: active
source_of_truth_for: next_session
last_updated: 2026-05-14
supersedes:
  - archive/sessions/2026-05-14_close.md
---

# NEXT SESSION

## Last Completed

- 4.3 line closed.
- Pre-memory-reorganization docs commit: `2e0eb39`.
- Final code commit: `d3de169`.
- `v4.3.12` tag remains at `3a0c40e`.
- `v4.3.12.post1` tag remains at `805e0b2`.
- Existing tags must not be moved.

## Start By

```bash
git status --short
git log --oneline --decorate -5
git tag --list 'v4.3.12*'
python -c "import aria; print(aria.__version__)"
```

Then read:

1. `memory/PROJECT_STATE.md`
2. `memory/AGENT_PROTOCOL.md`
3. A task-specific file only if needed.

## Known Caveats

- 4.3 is closed; do not reopen it for new features.
- `Memoria/`, if present, is a temporary local export folder and should not be
  treated as canonical.
- `memory/archive/` is history only.
- Reports generated before `805e0b2` / `d3de169` may still show older wording;
  new reports should use the fixed templates.
- Network may be restricted; LiteLLM cost-map warnings are non-fatal.
- Matplotlib may use a temporary cache if the user config dir is not writable.

## If Samael Asks For v4.4

Do not start from implementation. First verify:

- `hc11_paired.h5mu` exists.
- `aria-chromatin-env` status.
- Existing `chromatin_qc.py`, `chromatin_agent.py`, and `chromatin_peaks.py`.
- `memory/roadmap/V44_SCATAC_PLAN.md` if it exists, otherwise create it before
  coding.
