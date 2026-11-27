# v4.9.7: release-repro patch, the first tag with the full CI green

`v4.9.7` contains **no analytical code changes** relative to `v4.9.5`/`v4.9.6`. It fixes the
tag-only `release-repro` job, which builds the five Docker images and runs in-container
smoke tests. Every scientific behavior, receipt and preprint number is unchanged. For the
audit itself see [`releases/release_notes_v4.9.5.md`](releases/release_notes_v4.9.5.md).

## Why
On the `v4.9.6` tag, `main`-equivalent jobs were green, but `release-repro` failed in the
in-container steps that the `D3` CI rework had added. That job runs only on tags, so it had
never been executed before:

| Problem | Fix |
|---|---|
| The chromatin, ATAC-seq and TOBIAS images ship no `pytest`; only the RNA env does. | Install `pytest` in the throwaway container for those steps only. Published images are unchanged. |
| `bash -lc` (a login shell) re-sources `/etc/profile` and drops the PATH that the micromamba entrypoint activates, so `python` was not found. | Use `bash -c`. |
| Worker images (by design) do not ship the orchestration layer (`litellm`/`pydantic`: `aria.agents`, script contracts, env manager, run ledger). Module-level agent imports broke test collection there. | Agent, contract and ledger imports in those test files are lazy and gated with `importorskip`. 0 assertions changed. The same tests still run in `pr-checks` and `aria-env`. |
| The CI collection guard (`test_ci_collects_all_tests`) only unwrapped `bash -lc`. | It now unwraps `bash -c` too. Path validation is unchanged. |

`release-repro` can now also run through `workflow_dispatch`, so a branch can validate the
Docker images before a release is tagged.

## Verification
- `fix/release-repro` was validated on GitHub through `workflow_dispatch` (run `45930304429`):
  - `release-repro` succeeded: 5 image builds, the 4 in-container smoke steps and the
    CycloneDX SBOM;
  - `pr-checks`, the RNA, Hi-C and integration lanes all succeeded.
- `secret-scan` under `workflow_dispatch` scans the full history and reports 62 pre-existing
  false positives in benchmark receipts: `model_key` (model identifiers), `module_key` and
  `tokenizer_sha256` (content hashes). Push and tag scans cover only new commits and are
  green. Whether to add a gitleaks allowlist is left to the maintainer.
- The in-container steps were reproduced locally in the matching conda worker envs:
  chromatin 65 passed, ATAC 44 passed, TOBIAS 38 passed.
- The preprint `v1` automated gate was refreshed (19/19), and `verify_paper_claims.py` was
  unchanged (46/47, B2 decision pending).
