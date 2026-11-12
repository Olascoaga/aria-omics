# CLI, Control Center & Headless Execution Guide

ARIA exposes three execution surfaces — the **Textual Control Center**, the **classic Rich terminal flow**, and the **programmatic headless runner (`run_headless`)** — plus diagnostic and capsule-auditing subcommands. All surfaces drive the exact same `OrchestratorAgent`, checkpoint contracts, `EnvironmentManager`, and `C1–C4` claim compiler.

---

## 1. Interactive Launch (`aria`)

```bash
conda activate aria-env

aria                  # Launch Textual Control Center (when TTY + textual installed)
aria --classic-tui    # Force classic Rich terminal prompts
aria --reproducible   # Interactive reproducible mode (disables Textual cockpit, exports capsule inputs)
aria --hypotheses     # Opt in to the quarantined HypothesisAgent (SPECULATIVE appendix)
```

See [Control Center Guide](CONTROL_CENTER.md) for keyboard shortcuts, intake validation, session-aware clipboard handling (`Ctrl+V` / SSH safety), and live `RunLedger` / `Decisions` / `Artifacts` tabs.

---

## 2. Ledger, Export & Capsule Subcommands

Every completed report directory (`~/.aria/reports/aria_<timestamp>_<slug>_<id>/`) contains structured provenance in `methodology.json` that can be diffed, packaged, and independently verified:

```bash
# Compare structured claims, tiers, ledger linkage, and provenance between two reports
aria diff <REPORT_DIR_A> <REPORT_DIR_B>

# Export RO-Crate metadata (ro-crate-metadata.json) and a portable SHA-256 capsule ZIP
aria export <REPORT_DIR>

# Verify a reproducibility capsule's hashes, single-experiment SQLite scope, and claim set
aria reproduce <CAPSULE.zip>
```

> **Note:** `aria reproduce` cryptographically verifies the capsule manifest, single-experiment `memory_snapshot.sqlite`, lockfile hashes, and public claim set, and emits a reproduction command template. It never silently re-executes compute from prose.

---

## 3. System & Environment Diagnostics (`aria doctor`)

Run `aria doctor` (or `aria-doctor`) before sensitive or archival analyses:

```bash
aria doctor --smoke       # Verify core imports, registry integrity, and mock pipeline health
aria doctor --llm         # Verify configured LLM provider connectivity and model routing
aria doctor --secrets     # Scan environment and ~/.aria/.env permissions for secret hygiene
aria doctor --benchmark   # Run synthetic DE + label-permutation W-CALIB negative controls
```

To verify that installed scientific Conda environments match the committed exact Linux lockfiles (`envs/*.linux-64.lock` and `envs/*.pip.lock`) and required binaries:

```bash
python -m aria.utils.environment_audit
python -m aria.utils.environment_audit --env aria-rna-env
python -m aria.utils.environment_audit --env aria-integration-env
```

---

## 4. Environment Variables

| Variable | Effect |
| :--- | :--- |
| `ARIA_AIR_GAPPED=1` | **Set before process launch.** Blocks all cloud LLM calls, Enrichr/network enrichment, GEO/SRA downloads, and subprocess network egress; requires a local LLM endpoint (e.g., Ollama). |
| `ARIA_NO_TUI=1` | Disables the Textual Control Center and uses the classic Rich CLI. |
| `ARIA_LLM_CACHE=0` | Disables on-disk prompt caching under `~/.aria/llm_cache/`. |
| `ARIA_CLIPBOARD_SOURCE` | Controls Ctrl+V clipboard backend (`windows`, `x11`, `wayland`, `mac`, or `off`). Remote SSH sessions default to terminal paste safety (`off`). |
| `ARIA_ALLOW_FILENAME_FALLBACK=1` | Development-only opt-in allowing sample-name group guessing when explicit metadata is absent (stamps `LOW` confidence / not publication-grade). |
| `ARIA_ALLOW_EXPERIMENTAL_HIC=1` | Experimental-only opt-in for the scaffolded Hi-C / Micro-C lane (stamps report as non-publication-grade). |

---

## 5. Programmatic Headless Execution (`run_headless`)

For batch automation, reproducible benchmarks, and CI end-to-end tests, call `aria.headless.run_headless`:

```python
from pathlib import Path
from aria.headless import run_headless

def study_policy(checkpoint_id: str, payload: dict) -> str:
    # Provide explicit study-specific checkpoint confirmations
    if checkpoint_id == "cp3_5_audit":
        return "Proceed anyway"  # Explicit acknowledgement for beta ATAC / paired Multiome
    return "1"

result = run_headless(
    data_path=Path("/path/to/study_directory"),
    question="Compare chromatin accessibility and paired gene expression across hippocampal subclasses.",
    answer_policy=study_policy,
    reproducible_mode=True,
    enable_hypotheses=False,
)

print(result["status"], result["report_dir"])
```

Every policy-selected checkpoint answer is persisted to `ARIAMemory` (`decisions` table with `made_by="headless_policy"`) and audited in `methodology.json` and `memory_snapshot.sqlite`.
