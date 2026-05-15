---
status: active
source_of_truth_for: agent_protocol
last_updated: 2026-05-14
supersedes:
  - archive/sessions/2026-05-14_close.md
---

# ARIA Agent Protocol

## Collaboration Style

- Use Spanish for planning, status, and conversational explanations.
- Use English in code comments, docstrings, commit messages, and technical docs.
- Diagnose before coding.
- Prefer concrete file/line evidence over speculation.
- Keep changes minimal and scoped.
- Do not rewrite unrelated modules.
- Do not move existing tags.
- Do not force-push.
- Own mistakes briefly and fix them; no long apologies.

## When Samael Reports A Bug

1. Inspect relevant code and outputs first.
2. Identify the root cause.
3. State the smallest safe fix.
4. Implement only the needed change.
5. Run targeted tests, then broader tests when risk warrants it.
6. Summarize clearly.

## When A Design Choice Exists

Offer options when useful:

```text
Option A (simple): ...
Option B (robust): ...
Option C (overbuilt): ...

Recommendation: B because ...
```

If Samael gives explicit approval, execute without re-litigating the plan.

## Session Close Protocol

At the end of a sprint:

1. Update `memory/PROJECT_STATE.md` if project state changed.
2. Update `memory/DECISIONS.md` only for new durable decisions.
3. Update `memory/NEXT_SESSION.md` with the live handoff.
4. Archive detailed historical notes under `memory/archive/sessions/`.
5. Update `memory/START_HERE.md` only if the entrypoint or major status changed.

Do not touch every memory file reflexively.

## Session Start Protocol

Run:

```bash
git status --short
git log --oneline --decorate -5
```

Then read:

1. `memory/START_HERE.md`
2. `memory/PROJECT_STATE.md`
3. `memory/NEXT_SESSION.md`
4. `memory/AGENT_PROTOCOL.md`

Read task-specific context only when needed.

## Report Integrity Rules

- Missing results stay missing.
- Low and insufficient confidence findings remain visible.
- Methods must be copy-pasteable and HTML-safe.
- Integrated interpretations must separate observed DE/pathway evidence from
  causal hypotheses.
- No dataset-specific runtime guardrails.
