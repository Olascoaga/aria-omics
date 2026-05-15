---
status: active
source_of_truth_for: entrypoint
last_updated: 2026-05-14
---

# ARIA START HERE

Read this file first in every ARIA session.

## Current Status

- The 4.3 line is closed.
- Current branch: `main`.
- Verify current HEAD with `git log --oneline --decorate -5`.
- Last pre-memory-reorganization docs commit: `2e0eb39`
  (`Document final v4.3.12 closeout`).
- Last code change: `d3de169` (`Remove dataset-specific narrative guardrails`).
- `v4.3.12` tag: `3a0c40e`.
- `v4.3.12.post1` tag: `805e0b2`.
- Working tree should be clean at handoff.

## Do Not

- Do not move existing tags.
- Do not start v4.4 scATAC unless Samael explicitly asks.
- Do not add dataset-specific runtime guardrails.
- Do not hardcode genes, perturbations, datasets, or report-rescue phrases.
- Do not allow silent mocks or fake scientific outputs.
- Do not treat files in `memory/archive/` as source of truth.

## Read Order

1. `memory/PROJECT_STATE.md`
2. `memory/NEXT_SESSION.md`
3. `memory/AGENT_PROTOCOL.md`
4. Task-specific file only if needed.

## Source-Of-Truth Hierarchy

1. `PROJECT_STATE.md` for current version, tags, validation, and status.
2. `NEXT_SESSION.md` for the next live handoff.
3. `DECISIONS.md` for accepted architecture and scientific-policy decisions.
4. `AGENT_PROTOCOL.md` for collaboration rules.
5. `REFERENCE_PATHS.md` for local paths.
6. `datasets/`, `architecture/`, and `roadmap/` for task-specific context.
7. `archive/` for history only.

If a fact appears in two active files and conflicts, update the wrong file
immediately. Each type of truth should have one owner.
