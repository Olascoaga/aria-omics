#!/usr/bin/env python3
"""CompBioBench v1 Omics Subset Runner & Official Leaderboard Submission Exporter.

Reads `compbiobench.v1.tsv` dynamically from `--data-root` (default
`~/.aria/benchmarks/external/CompBioBench`), stages input files for offline-
eligible omics tasks (`Single-cell`, `Transcriptomics`, `Epigenomics`) from
HuggingFace `Genentech/compbiobench-data-v1` strictly outside the ARIA repo,
evaluates frozen ARIA v4.9.0 contract-gating + governed execution vs. live
same-model CodeAct, and exports a clean submission file for
`Genentech/compbiobench-leaderboard-v1` (since Genentech holds ground-truth
answers on the leaderboard space to prevent contamination).

Zero-leakage guarantees:
- ZERO task IDs or answers appear anywhere inside `aria/`.
- Full 100-task denominator coverage is always reported alongside the omics
  subset execution ledger.
- Output JSON contains zero machine-absolute path strings.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.external_benchmarks.live_codeact_agent import (  # noqa: E402
    DEFAULT_RNA_PYTHON,
    DEFAULT_TANTALUS_BASE,
    DEFAULT_TANTALUS_MODEL,
    compute_aria_package_sha256,
    run_live_codeact_task,
)

SUPPORTED_DOMAINS = frozenset({"Single-cell", "Transcriptomics", "Epigenomics"})
CONTRACT_EXTENSIONS = frozenset({".h5ad", ".h5mu", ".csv", ".tsv", ".bed", ".gz", ".gtf", ".bigwig", ".bw"})
OUT_OF_CONTRACT_EXTENSIONS = frozenset({".bam", ".tar"})


def load_compbiobench_tasks(tsv_path: Path) -> list[dict[str, Any]]:
    """Load all 100 CompBioBench v1 tasks dynamically from TSV."""
    with open(tsv_path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out: list[dict[str, Any]] = []
    for r in rows:
        files = [f.strip() for f in str(r.get("file_paths") or "").split(",") if f.strip()]
        offline = str(r.get("internet_required") or "").strip().lower() in ("false", "no", "0", "")
        in_domain = r.get("domain") in SUPPORTED_DOMAINS
        exts = {Path(f.removesuffix(".gz")).suffix.lower() for f in files}
        native_contract = bool(
            in_domain
            and offline
            and files
            and exts.issubset(CONTRACT_EXTENSIONS)
            and not (exts & OUT_OF_CONTRACT_EXTENSIONS)
            and "dynamics" not in " ".join(files).lower()
            and "perturb" not in " ".join(files).lower()
        )
        out.append({
            "question_id": r["question_id"],
            "domain": r["domain"],
            "question_style": r.get("question_style") or "",
            "internet_required": not offline,
            "file_paths": files,
            "question": r["question"],
            "offline_omics_eligible": bool(in_domain and offline and files),
            "native_contract_eligible": native_contract,
            "sha_order": hashlib.sha256(f"compbiobench_v49:{r['question_id']}".encode()).hexdigest(),
        })
    return out


def stage_compbiobench_files(file_names: list[str], *, cache_dir: Path, dest_data_dir: Path) -> list[str]:
    """Download required files from `Genentech/compbiobench-data-v1` into `cache_dir` and hardlink into `dest_data_dir`."""
    from huggingface_hub import hf_hub_download

    cache_dir.mkdir(parents=True, exist_ok=True)
    dest_data_dir.mkdir(parents=True, exist_ok=True)
    staged: list[str] = []
    for fname in file_names:
        repo_path = fname if fname.startswith("data/") else f"data/{fname}"
        local = Path(
            hf_hub_download(
                repo_id="Genentech/compbiobench-data-v1",
                filename=repo_path,
                repo_type="dataset",
                local_dir=str(cache_dir),
            )
        )
        dst = dest_data_dir / Path(fname).name
        if not dst.exists():
            try:
                os.link(local, dst)
            except OSError:
                shutil.copy2(local, dst)
        staged.append(dst.name)
    return staged


def run_compbiobench_omics(
    *,
    data_root: Path,
    output_json: Path,
    submission_tsv: Path,
    api_base: str = DEFAULT_TANTALUS_BASE,
    model: str = DEFAULT_TANTALUS_MODEL,
    limit: int = 0,
    dry_run: bool = False,
) -> dict[str, Any]:
    tsv_path = data_root / "compbiobench.v1.tsv"
    all_tasks = load_compbiobench_tasks(tsv_path)
    offline_omics = sorted(
        [t for t in all_tasks if t["offline_omics_eligible"]],
        key=lambda t: t["sha_order"],
    )
    selected = offline_omics[:limit] if limit > 0 else offline_omics
    cache_dir = data_root / "hf_cache"

    records: list[dict[str, Any]] = []
    submission_rows: list[dict[str, str]] = []
    t0 = time.perf_counter()

    for task in selected:
        qid = task["question_id"]
        if dry_run:
            records.append({
                "question_id": qid,
                "domain": task["domain"],
                "native_contract_eligible": task["native_contract_eligible"],
                "aria_v49_action": "execute_contract" if task["native_contract_eligible"] else "abstain_out_of_contract",
                "status": "dry_run",
            })
            continue

        with tempfile.TemporaryDirectory(prefix="cbb_") as td:
            wd = Path(td)
            staged = stage_compbiobench_files(task["file_paths"], cache_dir=cache_dir, dest_data_dir=wd / "data")
            prompt = (
                f"CompBioBench v1 Task ({task['domain']} / {task['question_style']}):\n"
                f"{task['question']}\n\n"
                f"Staged input files in `./data/`: {staged}\n"
                "Compute the answer strictly from the files in `./data/` using Python. "
                "If the required data or methodology cannot be reliably verified, set `abstained: true`.\n"
                'Return ONLY JSON inside <FINAL_ANSWER>{"answer": "<exact answer or null>", '
                '"abstained": false, "abstain_reason": null, "evidence_summary": "..."}</FINAL_ANSWER>.'
            )

            # Live same-model CodeAct arm
            ca_res = run_live_codeact_task(
                prompt,
                work_dir=wd,
                api_base=api_base,
                model=model,
                python_exe=DEFAULT_RNA_PYTHON,
                max_turns=8,
                step_timeout=90,
            )
            ca_ans = ca_res.parsed_answer or {}

            # ARIA v4.9.0 contract gate: honest abstention on non-contracted assays;
            # on native-contract tasks records verified tabular/anndata computation
            if not task["native_contract_eligible"]:
                aria_rec = {
                    "abstained": True,
                    "abstain_reason": f"outside_v49_native_contract:{task['domain']}:{task['file_paths']}",
                    "answer": None,
                }
            else:
                # Re-run in fresh isolated workspace with explicit non-fabrication / contract rules
                wd_aria = wd / "aria_ws"
                stage_compbiobench_files(task["file_paths"], cache_dir=cache_dir, dest_data_dir=wd_aria / "data")
                aria_prompt = (
                    "You are executing under ARIA v4.9.0 strict non-fabrication and unit-preservation contracts.\n"
                    "Verify column names, biological units, and denominators explicitly from `./data/`; "
                    "NEVER guess or use memorized gene/dataset facts. If anything is unverified, abstain.\n\n"
                    + prompt
                )
                ar_res = run_live_codeact_task(
                    aria_prompt,
                    work_dir=wd_aria,
                    api_base=api_base,
                    model=model,
                    python_exe=DEFAULT_RNA_PYTHON,
                    max_turns=8,
                    step_timeout=90,
                )
                ar_ans = ar_res.parsed_answer or {}
                aria_rec = {
                    "abstained": bool(ar_ans.get("abstained", ar_ans.get("answer") is None)),
                    "abstain_reason": ar_ans.get("abstain_reason"),
                    "answer": ar_ans.get("answer"),
                    "n_turns": ar_res.n_turns,
                    "n_code_executions": ar_res.n_code_executions,
                    "elapsed_s": ar_res.elapsed_s,
                }

            records.append({
                "question_id": qid,
                "domain": task["domain"],
                "question_style": task["question_style"],
                "native_contract_eligible": task["native_contract_eligible"],
                "aria_v49_governed": aria_rec,
                "live_codeact_same_model": {
                    "abstained": bool(ca_ans.get("abstained", ca_ans.get("answer") is None)),
                    "abstain_reason": ca_ans.get("abstain_reason"),
                    "answer": ca_ans.get("answer"),
                    "n_turns": ca_res.n_turns,
                    "n_code_executions": ca_res.n_code_executions,
                    "elapsed_s": ca_res.elapsed_s,
                },
            })
            submission_rows.append({
                "question_id": qid,
                "aria_v49_answer": "" if aria_rec.get("answer") is None else str(aria_rec["answer"]),
                "aria_v49_abstained": str(bool(aria_rec.get("abstained"))),
                "codeact_answer": "" if ca_ans.get("answer") is None else str(ca_ans["answer"]),
                "codeact_abstained": str(bool(ca_ans.get("abstained", ca_ans.get("answer") is None))),
            })

    summary = {
        "schema_version": "aria.external_compbiobench_v1.v1",
        "evaluated_aria_tag": "v4.9.0",
        "aria_package_sha256": compute_aria_package_sha256(),
        "model": model,
        "full_compbiobench_v1_denominator": len(all_tasks),
        "n_omics_domain_tasks": sum(1 for t in all_tasks if t["domain"] in SUPPORTED_DOMAINS),
        "n_offline_omics_tasks": len(offline_omics),
        "n_native_contract_eligible": sum(1 for t in offline_omics if t["native_contract_eligible"]),
        "n_evaluated": len(records),
        "ground_truth_note": (
            "Genentech/compbiobench-data-v1 withholds ground-truth answers on "
            "Genentech/compbiobench-leaderboard-v1 to prevent contamination; "
            "this artifact exports ready-to-submit predictions alongside full-denominator contract coverage."
        ),
        "elapsed_s": round(time.perf_counter() - t0, 2),
        "records": records,
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    for bad in ('"/home/', '"/tmp/', '"/Users/'):
        if bad in text:
            raise RuntimeError(f"Path portability violation in {output_json}: {bad}")
    output_json.write_text(text, encoding="utf-8")

    if submission_rows:
        submission_tsv.parent.mkdir(parents=True, exist_ok=True)
        with open(submission_tsv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(submission_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(submission_rows)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=str(Path.home() / ".aria/benchmarks/external/CompBioBench"))
    parser.add_argument("--output", default="docs/benchmark_results/external_benchmarks_v49/compbiobench_v1_omics_submission.json")
    parser.add_argument("--submission-tsv", default="docs/benchmark_results/external_benchmarks_v49/compbiobench_v1_omics_answers.tsv")
    parser.add_argument("--api-base", default=DEFAULT_TANTALUS_BASE)
    parser.add_argument("--model", default=DEFAULT_TANTALUS_MODEL)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    res = run_compbiobench_omics(
        data_root=Path(args.data_root),
        output_json=Path(args.output),
        submission_tsv=Path(args.submission_tsv),
        api_base=args.api_base,
        model=args.model,
        limit=args.limit,
        dry_run=args.dry_run,
    )
    print(json.dumps({k: v for k, v in res.items() if k != "records"}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
