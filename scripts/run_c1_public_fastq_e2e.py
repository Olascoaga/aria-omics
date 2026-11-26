#!/usr/bin/env python3
"""c1_public_fastq_e2e preprint-freeze harness — public paired GSE52778 FASTQ-to-report E2E.

Clean-checkout, 100% public end-to-end evidence for Claim 1 (v4.9.4+).
Drives ARIA's REAL raw bulk RNA-seq pipeline (fastp/STAR alignment ->
featureCounts -> paired pyDESeq2 ``~ donor + condition`` -> F3-governed report)
over the public human airway smooth muscle glucocorticoid dataset
(GEO GSE52778 / ENA SRP033351, Himes et al., PLOS ONE 2014; 4 matched human
donors x 2 conditions: untreated vs dexamethasone, 8 paired-end FASTQ pairs)
through ``aria.headless.run_headless`` under a frozen deterministic checkpoint
policy, then publishes a portable reproducibility capsule under the freeze root.

Honesty contract:
  * 100% public input data (``~/.aria/public_atlas_datasets/02_bulk_rna_raw_fastq/inputs``
    or ``$ARIA_GSE52778_FASTQ_DIR``). Zero unpublished or restricted lab cohorts.
  * No canned success: STAR, featureCounts, and paired pyDESeq2 run for real in
    their pinned scientific environments (``aria-rnaseq-env`` / ``aria-rna-env``).
  * CP 2.5 explicitly exercises ARIA's cross-arm paired biological replicate
    detection and selects ``Paired biological replicates — model matched
    donor/unit block``, yielding ``~ donor + condition`` in DESeq2.
  * The LLM is a hermetic deterministic in-process double (no network, $0.00).
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

CAPSULE_SCHEMA = "aria.preprint_freeze.c1_public_fastq_capsule.v1"
LANE_ID = "c1_public_fastq_e2e"

DEFAULT_DATA_DIR = (
    os.environ.get("ARIA_GSE52778_FASTQ_DIR", "").strip()
    or "~/.aria/public_atlas_datasets/02_bulk_rna_raw_fastq/inputs"
)
DEFAULT_QUESTION = (
    "Analyze paired-end raw bulk RNA-seq FASTQ libraries from primary human "
    "airway smooth muscle cells treated with dexamethasone versus untreated "
    "controls across four matched human donor lines (GEO GSE52778 / ENA "
    "SRP033351, Homo sapiens hg38). Identify differentially expressed genes "
    "for dexamethasone versus untreated using a paired donor-blocked design."
)

CAPSULE_ENVIRONMENTS = ("aria-env", "aria-rnaseq-env", "aria-rna-env")
FREEZE_N_CPUS = 30
GSE52778_GROUPS = {
    "dexamethasone": (
        "dexamethasone_rep1",
        "dexamethasone_rep2",
        "dexamethasone_rep3",
        "dexamethasone_rep4",
    ),
    "untreated": (
        "untreated_rep1",
        "untreated_rep2",
        "untreated_rep3",
        "untreated_rep4",
    ),
}
GSE52778_CONTRASTS = (
    {"numerator": "dexamethasone", "denominator": "untreated"},
)
GSE52778_MANUAL_GROUP_ASSIGNMENT = "; ".join(
    f"{group}={','.join(samples)}" for group, samples in GSE52778_GROUPS.items()
)

_DE_TABLE_PATTERNS = (
    "*deseq2*.tsv", "*deseq*.tsv", "*differential*.tsv", "*_de.tsv",
    "*_de_genes.tsv", "*de_results*.tsv", "*deseq2*.csv",
    "*differential*.csv", "*_de.csv", "*_de_genes.csv",
)
_DE_FIGURE_PATTERNS = (
    "*volcano*.svg", "*deseq2*.svg", "*differential*.svg", "*_de*.svg",
    "*ma_plot*.svg", "*volcano*.png", "*_de*.png",
)


class FrozenLLMDouble:
    """Hermetic, deterministic in-process stand-in for LLMProvider."""

    marker = "[frozen-llm-double: proposal suppressed for reproducible freeze]"
    intent_response = json.dumps({
        "analysis_type": "differential",
        "biological_entities": [],
        "comparison": "dexamethasone vs untreated",
        "key_modalities_needed": ["RNA"],
        "complexity": "moderate",
        "summary": "Paired bulk RNA-seq differential expression (dexamethasone vs untreated).",
    }, sort_keys=True)
    plan_response = json.dumps({
        "steps": [{
            "order": 1,
            "agent": "bulk_rna_agent",
            "analysis": "Bulk RNA differential expression",
            "depends_on": [],
            "can_parallel": False,
        }],
        "contrasts": list(GSE52778_CONTRASTS),
        "integration_needed": False,
        "integration_type": "none",
        "estimated_complexity": "medium",
        "rationale": (
            "Run paired donor-blocked differential expression for dexamethasone versus untreated."
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


def frozen_gse52778_answer_policy(cp_num, question: str, options: list) -> str:
    """Confirm the public paired GSE52778 design through ARIA's real checkpoint contract."""
    from aria.headless import default_answer_policy

    if cp_num == 2.1:
        return GSE52778_MANUAL_GROUP_ASSIGNMENT
    if cp_num == 2.2:
        hg38 = next(
            (o for o in options if "homo sapiens" in o.lower() and "hg38" in o.lower()),
            None,
        )
        if hg38:
            return hg38
    if cp_num == 2.3:
        condition = next(
            (option for option in options if option.strip().lower() == "condition"),
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
    # CP 2.5 intentionally delegates to default_answer_policy so the real
    # paired-donor block policy is exercised end-to-end!
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
        hits = sorted(directory.rglob(pattern))
        if hits:
            return hits[0]
    return None


def _portable_json_paths(value, replacements: dict[str, str]):
    if isinstance(value, dict):
        return {
            key: _portable_json_paths(nested, replacements)
            for key, nested in value.items()
        }
    if isinstance(value, list):
        return [_portable_json_paths(nested, replacements) for nested in value]
    if isinstance(value, str) and Path(value).is_absolute():
        portable = Path(value).name or "artifact"
        replacements[value] = portable
        return portable
    return value


def _portable_report_html(text: str, replacements: dict[str, str]) -> str:
    for source, portable in sorted(
        replacements.items(), key=lambda item: len(item[0]), reverse=True
    ):
        text = text.replace(source, portable)
        text = text.replace(html.escape(source, quote=True), portable)
    if re.search(r"(?<!:)\/(?:home|Users|mnt|tmp|var|opt|private)\/", text):
        raise RuntimeError(
            "completed report still contains a machine-absolute path after "
            "portable publication"
        )
    return text


def _inventory(report_dir: Path, subdir: str) -> list[dict]:
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
    if isinstance(value, dict):
        if any(k in value for k in ("n_significant", "n_upregulated", "n_downregulated")):
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
    """Copy the canonical artifacts out and assemble the portable capsule."""
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
            "emit a partial capsule"
        )
    if de_figure is None:
        raise RuntimeError(
            "completed report has no differential-expression figure; refusing to "
            "emit a partial capsule"
        )

    methodology_payload = json.loads(methodology.read_text(encoding="utf-8"))
    path_replacements: dict[str, str] = {}
    portable_methodology = _portable_json_paths(
        methodology_payload, path_replacements
    )
    portable_report = _portable_report_html(
        report_html.read_text(encoding="utf-8"), path_replacements
    )

    copies = {
        "de_results.tsv": de_table,
        "fig1_gse52778_bulk_de.svg": de_figure,
    }
    artifacts: list[dict] = []
    portable_copies = {
        "report.html": (portable_report, report_html.name),
        "methodology.json": (
            json.dumps(portable_methodology, indent=2, sort_keys=True) + "\n",
            methodology.name,
        ),
    }
    for out_name, (content, source_name) in portable_copies.items():
        dst = output_dir / out_name
        dst.write_text(content, encoding="utf-8")
        artifacts.append({
            "name": out_name,
            "source_name": source_name,
            "size_bytes": dst.stat().st_size,
            "sha256": _sha256(dst),
        })
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
        _scan_de_summary(portable_methodology, de_summary)
    except Exception:
        de_summary = []

    design_used = (portable_methodology.get("design") or {}).get("design_formula")
    version = collect_version_metadata(repo_root)
    capsule = {
        "schema_version": CAPSULE_SCHEMA,
        "lane_id": LANE_ID,
        "claim": "claim_1",
        "aria_version": version.get("aria_version"),
        "git_commit": version.get("git_commit"),
        "git_tree_sha": version.get("git_tree_sha"),
        "git_dirty": version.get("git_dirty"),
        "experiment_id": experiment_id,
        "data_source": "public-atlas://02_bulk_rna_raw_fastq/GSE52778",
        "accession": "GEO GSE52778 / ENA SRP033351",
        "doi": "10.1371/journal.pone.0099625",
        "pipeline": "raw FASTQ -> fastp/STAR -> featureCounts -> paired pyDESeq2 (~ donor + condition) -> F3-governed report",
        "llm": "frozen-deterministic-double (hermetic, no network)",
        "checkpoint_policy": "c1_public_fastq_e2e.frozen_gse52778_answer_policy",
        "design_formula_recorded": design_used,
        "design_contract": {
            "groups": {
                group: list(samples) for group, samples in GSE52778_GROUPS.items()
            },
            "contrasts": list(GSE52778_CONTRASTS),
            "factor": "condition",
            "expected_paired_formula": "~ donor + condition",
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
        policy=frozen_gse52778_answer_policy,
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
        default="docs/benchmark_results/preprint_v1/claim_1/public_fastq_e2e",
    )
    parser.add_argument("--question", default=DEFAULT_QUESTION)
    parser.add_argument("--timeout", type=float, default=14_400.0)
    args = parser.parse_args(argv)

    capsule = run(args.data_dir, args.output_dir, args.question, args.timeout)
    print(json.dumps({
        "lane_id": capsule["lane_id"],
        "experiment_id": capsule["experiment_id"],
        "design_formula_recorded": capsule.get("design_formula_recorded"),
        "de_summary": capsule["de_summary"],
        "canonical_artifacts": [a["name"] for a in capsule["canonical_artifacts"]],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
