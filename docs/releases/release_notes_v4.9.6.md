# v4.9.6: CI patch for the v4.9.5 tech-debt audit release

`v4.9.6` contains **no analytical code changes** relative to `v4.9.5`. It makes the CI that
`v4.9.5` introduced (full `pytest tests/` collection in every lane) pass on GitHub. Every
scientific behavior, receipt and preprint number of `v4.9.5` is unchanged. See
[`release_notes_v4.9.5.md`](release_notes_v4.9.5.md) for the audit itself.

## Why
The `v4.9.5` tag was the first run of the reworked CI (`D3`) on GitHub. Four jobs failed, and
none of the failures was scientific:

| Job | Cause | Fix |
|---|---|---|
| `pr-checks` (43 failed) | The pure-pip lane has no scanpy, seaborn, numba or conda stacks, and some tests did not skip without them. | `pytest.importorskip(...)` gates and a shared `require_aria_stack()` helper (`tests/conftest.py`). Every gated test still runs and passes in a lane that has its dependencies. |
| `rna-env-tests` (2 failed) | Two contract tests dispatch to the `rna` conda stack, which is not named `aria-rna-env` on CI. | Gated by `require_aria_stack("rna")`; they run where the stack exists. |
| `integration-env-tests` (13 failed) | `envs/aria-integration-env.yml` did not pin snapatac2's `pyarrow`/`polars`. pip resolved `pyarrow 26` (requires NumPy ≥ 2) into the `numpy<2` env, so fragment import failed and the stage honest-skipped (`fragment_import_failed`). | Pinned `polars==1.38.1`, `polars-runtime-32==1.38.1` and `pyarrow==23.0.1` to the validated pip lock. |
| `release-repro` (1 failed) | The reversed-contrast oracle required `padj` parity to `rtol=1e-4`. The RNA image (numpy 2.2.6, scipy 1.14.1) differs by up to ~2%, because pydeseq2 IRLS starts from another parametrization. | The oracle now requires identical calls away from the threshold, plus `log10(padj)` agreement. The exact LFC sign and magnitude checks are unchanged and still kill the LFC-sign mutant. |

**Zero assertions were changed** by the gating commits
(`git diff v4.9.5 -- tests/ | grep -E "^[-+]\s*assert"` is empty apart from the documented oracle
rewrite).

## Verification
The five CI lanes were reproduced locally before pushing:
- `pr-checks`: a clean Python 3.11 venv with the exact pip set and no conda on PATH;
- the RNA and Hi-C lanes: cloned conda envs with `pip install -e .`, both with conda visible
  and with conda hidden;
- the integration lane, run the same way;
- the full `aria-env` suite: 3295 passed, 0 new skips.

The preprint `v1` automated gate was refreshed for `v4.9.6` (19/19 lanes), and
`verify_paper_claims.py` was unchanged from `v4.9.5` (46/47, B2 decision pending).
