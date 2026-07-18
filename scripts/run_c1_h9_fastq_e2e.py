#!/usr/bin/env python3
"""c1_h9_fastq_e2e preprint-freeze harness — FASTQ-to-report bulk-RNA E2E.

Clean-checkout, end-to-end evidence for Claim 1. Drives the REAL raw-RNA
pipeline (STAR alignment -> featureCounts -> pyDESeq2 -> F3-governed report)
over the local H9 paired-end FASTQ inputs through ``aria.headless.run_headless``
under a frozen, deterministic checkpoint policy, then publishes a portable
reproducibility capsule under the freeze output root.

Honesty contract (see memory/DECISIONS.md and the preprint audit invariants):

  * No canned success. STAR, featureCounts and pyDESeq2 run for real in their
    own scientific environments (dispatched by ``EnvironmentManager`` into
    ``aria-rnaseq-env`` / ``aria-rna-env``). This harness only orchestrates the
    real run and packages the resulting evidence.
  * The genome/index is auto-resolved by ARIA's own setup agent to the managed
    ``~/.aria/genomes/hg38``. STAR builds the managed index when it is absent or
    incomplete and reuses it thereafter. Nothing about the reference is
    hard-coded here.
  * The LLM is a hermetic, deterministic in-process double (no network, no
    cost). ARIA's scientific results and the F3-governed report are
    evidence-derived, not free-text; the double only satisfies incidental
    proposal calls so the run is fully offline and reproducible.
  * The frozen H9 checkpoint policy supplies the independently known sample
    grouping at CP2.1 and authorizes the two explicit contrasts at CP2. This is
    dataset-specific validation metadata in the harness, not runtime inference:
    DesignAgent still parses, validates and confirms the groups, replicate
    structure and factor through the real checkpoint state machine.
  * If the completed report lacks a differential-expression table or figure, the
    harness fails loudly rather than emitting a partial capsule.

The declared, receipt-covered artifacts are exactly the five files this harness
emits under ``--output-dir``:

    capsule.json  report.html  methodology.json  de_results.tsv  fig1_h9_bulk_de.svg

``capsule.json`` additionally records, by name + sha256, every figure and table
the report produced, so nothing the run emitted is hidden even though only the
canonical DE table/figure are copied out.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

CAPSULE_SCHEMA = "aria.preprint_freeze.c1_h9_capsule.v1"
LANE_ID = "c1_h9_fastq_e2e"

DEFAULT_DATA_DIR = "~/Samael/H9-RNA/raw_fastq"
# The question is the user's own input; naming the knockout conditions here is
# legitimate intake, not a hidden policy injection. The design is inferred by
# ARIA and confirmed by the deterministic checkpoint policy.
DEFAULT_QUESTION = (
    "In H9 human bulk RNA-seq, identify differentially expressed genes for the "
    "BMAL1 knockout versus wild type and for the REV-ERBa knockout versus wild "
    "type."
)

# Environments whose committed lockfiles pin the tool stack that produced this
# evidence: the orchestrator env plus the two scientific envs the pipeline
# dispatches into.
CAPSULE_ENVIRONMENTS = ("aria-env", "aria-rnaseq-env", "aria-rna-env")
FREEZE_N_CPUS = 30
H9_GROUPS = {
    "B": ("B1", "B2", "B3"),
    "R": ("R1", "R2", "R3"),
    "WT": ("WT1", "WT2", "WT3"),
}
H9_CONTRASTS = (
    {"numerator": "B", "denominator": "WT"},
    {"numerator": "R", "denominator": "WT"},
)
H9_MANUAL_GROUP_ASSIGNMENT = "; ".join(
    f"{group}={','.join(samples)}" for group, samples in H9_GROUPS.items()
)

# DE table / figure discovery in the completed report directory. Patterns are
# DE-biased and fall back conservatively; a completed bulk-DE report always
# emits both, so a miss is a real failure, not a silent skip.
_DE_TABLE_PATTERNS = (
    "*deseq2*.tsv", "*deseq*.tsv", "*differential*.tsv", "*_de.tsv",
    "*de_results*.tsv", "*deseq2*.csv", "*differential*.csv", "*_de.csv",
)
_DE_FIGURE_PATTERNS = (
    "*volcano*.svg", "*deseq2*.svg", "*differential*.svg", "*_de*.svg",
    "*ma_plot*.svg", "*volcano*.png", "*_de*.png",
)


class FrozenLLMDouble:
    """Hermetic, deterministic stand-in for :class:`LLMProvider`.

    Returns stable neutral text for any incidental proposal call so a headless
    run is fully offline and reproducible. It never contacts a provider, never
    records usage and never injects dataset-specific biology. It satisfies the
    per-execution provider seam the orchestrator uses (``for_execution``).
    """

    marker = "[frozen-llm-double: proposal suppressed for reproducible freeze]"
    intent_response = json.dumps({
        "analysis_type": "differential",
        "biological_entities": [],
        "comparison": "",
        "key_modalities_needed": ["RNA"],
        "complexity": "moderate",
        "summary": "Differential expression analysis.",
    }, sort_keys=True)
    plan_response = json.dumps({
        "steps": [{
            "order": 1,
            "agent": "bulk_rna_agent",
            "analysis": "Bulk RNA differential expression",
            "depends_on": [],
            "can_parallel": False,
        }],
        "contrasts": list(H9_CONTRASTS),
        "integration_needed": False,
        "integration_type": "none",
        "estimated_complexity": "medium",
        "rationale": (
            "Run the two pre-specified contrasts against the shared reference."
        ),
    }, sort_keys=True)

    def for_execution(self, experiment_id, usage_log, egress_policy):
        return self

    def complete(self, prompt, system="", tier=None, max_tokens=1024,
                 messages=None):
        if "analysis_type" in prompt and "key_modalities_needed" in prompt:
            return self.intent_response
        if '"steps"' in prompt and '"contrasts"' in prompt:
            return self.plan_response
        return self.marker

    def complete_heavy(self, prompt, system="", max_tokens=2048, messages=None):
        return self.marker

    def complete_medium(self, prompt, system="", max_tokens=1024, messages=None):
        return self.marker

    def complete_light(self, prompt, system="", max_tokens=512, messages=None):
        return self.marker

    def get_active_model(self, tier=None):
        return None


def frozen_h9_answer_policy(cp_num, question: str, options: list) -> str:
    """Confirm the frozen H9 design through ARIA's real checkpoint contract."""
    from aria.headless import default_answer_policy

    if cp_num == 2.1:
        return H9_MANUAL_GROUP_ASSIGNMENT
    if cp_num == 2.3:
        condition = next(
            (option for option in options
             if option.strip().lower() == "condition"),
            None,
        )
        if condition:
            return condition
    if cp_num == 2:
        recommended = next(
            (option for option in options
             if option.strip().lower() == "run recommended plan only"),
            None,
        )
        if recommended:
            return recommended
    return default_answer_policy(cp_num, question, options)


def _sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _find_first(report_dir: Path, subdir: str, patterns) -> Path | None:
    directory = report_dir / subdir
    if not directory.is_dir():
        directory = report_dir
    for pattern in patterns:
        hits = sorted(directory.glob(pattern))
        if hits:
            return hits[0]
    return None


def _inventory(report_dir: Path, subdir: str) -> list[dict]:
    """Name + sha256 + size for every artifact the report emitted in ``subdir``."""
    directory = report_dir / subdir
    entries: list[dict] = []
    if directory.is_dir():
        for path in sorted(directory.rglob("*")):
            if path.is_file():
                entries.append({
                    "name": path.relative_to(directory).as_posix(),
                    "size_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                })
    return entries


def _env_locks(repo_root: Path) -> list[dict]:
    """Committed lockfile identity for each environment used by this lane."""
    locks: list[dict] = []
    for env_name in CAPSULE_ENVIRONMENTS:
        entry: dict = {"env_name": env_name, "lock_file": None, "sha256": None}
        for fname in (f"{env_name}.linux-64.lock", f"{env_name}.pip.lock"):
            lock_path = repo_root / "envs" / fname
            if lock_path.is_file():
                entry["lock_file"] = lock_path.relative_to(repo_root).as_posix()
                entry["sha256"] = _sha256(lock_path)
                break
        locks.append(entry)
    return locks


def _scan_de_summary(value, out: list[dict]) -> None:
    """Best-effort, schema-agnostic extraction of DE contrast summaries."""
    if isinstance(value, dict):
        if any(k in value for k in ("n_significant", "n_upregulated",
                                    "n_downregulated")):
            summary = {
                k: value[k]
                for k in ("name", "contrast", "n_significant",
                          "n_upregulated", "n_downregulated")
                if k in value
            }
            if summary:
                out.append(summary)
        for nested in value.values():
            _scan_de_summary(nested, out)
    elif isinstance(value, list):
        for nested in value:
            _scan_de_summary(nested, out)


def publish(repo_root: Path, report_dir: Path, output_dir: Path,
            *, experiment_id: str, decisions: list) -> dict:
    """Copy the canonical artifacts out and assemble the portable capsule.

    Returns the capsule payload. Raises loudly if the completed report is
    missing its DE table or DE figure.
    """
    from aria.version import collect_version_metadata

    report_dir = Path(report_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report_html = report_dir / "report.html"
    methodology = report_dir / "methodology.json"
    if not report_html.is_file():
        raise RuntimeError(f"completed run has no report.html at {report_dir}")
    if not methodology.is_file():
        raise RuntimeError(f"completed run has no methodology.json at {report_dir}")

    de_table = _find_first(report_dir, "tables", _DE_TABLE_PATTERNS)
    de_figure = _find_first(report_dir, "figures", _DE_FIGURE_PATTERNS)
    if de_table is None:
        raise RuntimeError(
            "completed report has no differential-expression table; refusing to "
            "emit a partial capsule (searched tables/ for DESeq2/DE outputs)"
        )
    if de_figure is None:
        raise RuntimeError(
            "completed report has no differential-expression figure; refusing to "
            "emit a partial capsule (searched figures/ for volcano/MA/DE plots)"
        )

    copies = {
        "report.html": report_html,
        "methodology.json": methodology,
        "de_results.tsv": de_table,
        "fig1_h9_bulk_de.svg": de_figure,
    }
    artifacts: list[dict] = []
    for out_name, src in copies.items():
        dst = output_dir / out_name
        shutil.copyfile(src, dst)
        artifacts.append({
            "name": out_name,
            "source_name": src.name,
            "size_bytes": dst.stat().st_size,
            "sha256": _sha256(dst),
        })

    de_summary: list[dict] = []
    try:
        _scan_de_summary(json.loads(methodology.read_text(encoding="utf-8")),
                         de_summary)
    except Exception:
        de_summary = []

    version = collect_version_metadata(repo_root)
    capsule = {
        "schema_version": CAPSULE_SCHEMA,
        "lane_id": LANE_ID,
        "claim": "claim_1",
        "aria_version": version.get("aria_version"),
        "git_commit": version.get("git_commit"),
        "git_tree_sha": version.get("git_tree_sha"),
        "experiment_id": experiment_id,
        "data_source": "validation-data://h9_bulk_fastq",
        "pipeline": "raw FASTQ -> fastp/STAR -> featureCounts -> pyDESeq2 -> "
                    "F3-governed report",
        "llm": "frozen-deterministic-double (hermetic, no network)",
        "checkpoint_policy": "c1_h9_fastq_e2e.frozen_h9_answer_policy",
        "design_contract": {
            "groups": {
                group: list(samples) for group, samples in H9_GROUPS.items()
            },
            "contrasts": list(H9_CONTRASTS),
            "factor": "condition",
        },
        "requested_cpus": FREEZE_N_CPUS,
        "checkpoint_decisions": [
            {"checkpoint": cp, "choice": choice} for cp, choice in decisions
        ],
        "environments": _env_locks(repo_root),
        "de_summary": de_summary,
        "canonical_artifacts": artifacts,
        "report_figures": _inventory(report_dir, "figures"),
        "report_tables": _inventory(report_dir, "tables"),
    }
    capsule_path = output_dir / "capsule.json"
    capsule_path.write_text(
        json.dumps(capsule, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return capsule


def run(data_dir: str, output_dir: str, question: str, timeout: float) -> dict:
    from aria.headless import run_headless

    data_path = str(Path(data_dir).expanduser())
    result = run_headless(
        data_path,
        question,
        policy=frozen_h9_answer_policy,
        reproducible_mode=True,
        enable_hypotheses=False,
        timeout=timeout,
        context_overrides={"n_cpus": FREEZE_N_CPUS},
        llm_provider=FrozenLLMDouble(),
    )
    if result.status != "completed" or not result.report_path:
        raise RuntimeError(
            f"headless E2E did not complete: status={result.status!r} "
            f"report_path={result.report_path!r}"
        )
    report_dir = Path(result.report_path).parent
    return publish(
        ROOT, report_dir, Path(output_dir),
        experiment_id=result.experiment_id,
        decisions=result.decisions,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    parser.add_argument(
        "--output-dir",
        default="docs/benchmark_results/preprint_v1/claim_1/h9_e2e",
    )
    parser.add_argument("--question", default=DEFAULT_QUESTION)
    parser.add_argument("--timeout", type=float, default=28_800.0)
    args = parser.parse_args(argv)

    capsule = run(args.data_dir, args.output_dir, args.question, args.timeout)
    print(json.dumps({
        "lane_id": capsule["lane_id"],
        "experiment_id": capsule["experiment_id"],
        "de_summary": capsule["de_summary"],
        "canonical_artifacts": [a["name"] for a in capsule["canonical_artifacts"]],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
