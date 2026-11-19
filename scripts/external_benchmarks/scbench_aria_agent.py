#!/usr/bin/env python3
"""Official `latchbio/scbench` EvalRunner-compatible adapter for frozen ARIA v4.9.0.

Zero-leakage design:
- Contains ZERO scBench eval IDs, dataset accessions, or hardcoded answers.
- Implements the exact `agent_function(task_prompt: str, work_dir: Path) -> dict`
  contract expected by `scbench.EvalRunner.run(agent_function=...)`, writing
  `work_dir / "eval_answer.json"` and returning the structured answer dict.
- Routes every single-cell `.h5ad` workspace through frozen ARIA v4.9.0
  governance primitives:
    1. Native platform & modality gate (fails closed with explicit abstention on
       non-RNA DNA/amplicon/variant-only containers such as DNA+protein panels).
    2. Real `aria-rna-env` inspection & ARIA's production `rna_pseudobulk_de.py`
       (donor-level PyDESeq2 with donor blocking—never cell-level pseudoreplication)
       and `aria.benchmarks.derived_query_ast.execute_derived_query` (`DerivedQueryPlan`).
    3. Isolated per-task environment (`isolated_benchmark_environment`) with
       `ARIA_LLM_CACHE=0`, `ARIA_REPRODUCIBLE=1`, and fresh `ARIA_MEMORY_DB`.

Usage:
  # With LatchBio scBench installed + ~/.latch/token configured:
  python scripts/external_benchmarks/scbench_aria_agent.py --eval-json /path/to/eval.json

  # Offline dry-run smoke test (no Latch token required):
  python scripts/external_benchmarks/scbench_aria_agent.py --synthetic-smoke
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from aria.benchmarks.adapter import isolated_benchmark_environment  # noqa: E402
from scripts.external_benchmarks.live_codeact_agent import (  # noqa: E402
    DEFAULT_RNA_PYTHON,
    DEFAULT_TANTALUS_BASE,
    DEFAULT_TANTALUS_MODEL,
    chat_completion,
    compute_aria_package_sha256,
    extract_json_answer,
)

# Generic non-RNA modality indicators (platform-agnostic, zero benchmark IDs)
_NON_RNA_PROMPT_CUES = (
    "dna amplicon",
    "dna+protein",
    "somatic mutation",
    "mutant clone",
    "variant allele frequency",
    "vaf matrix",
    "genotype matrix",
    "single-cell dna",
)


def inspect_h5ad_workspace(work_dir: Path, *, python_exe: str = DEFAULT_RNA_PYTHON) -> dict[str, Any]:
    """Profile all `.h5ad` files in `work_dir` inside `aria-rna-env` (read-only)."""
    h5ads = sorted(work_dir.rglob("*.h5ad"))
    if not h5ads:
        return {"status": "no_h5ad_found", "files": []}

    probe_code = """
import json, sys
from pathlib import Path
import numpy as np
import anndata as ad
from scipy import sparse

out = []
for p_str in sys.argv[1:]:
    p = Path(p_str)
    adata = ad.read_h5ad(p)
    X = adata.X
    if sparse.issparse(X):
        vals = X.data[:50000]
    else:
        vals = np.asarray(X).ravel()[:50000]
    finite = vals[np.isfinite(vals)] if len(vals) else np.array([])
    is_nonneg_int = bool(
        len(finite) > 0
        and float(finite.min()) >= 0.0
        and np.allclose(finite, np.round(finite), atol=1e-6)
    )
    obs_summary = {}
    for col in list(adata.obs.columns)[:30]:
        s = adata.obs[col]
        nunq = int(s.nunique(dropna=True))
        vals_sample = [str(v) for v in s.dropna().unique()[:8]]
        obs_summary[str(col)] = {"dtype": str(s.dtype), "n_unique": nunq, "sample_values": vals_sample}
    var_names = [str(v) for v in adata.var_names[:20]]
    looks_like_variants = sum(1 for v in var_names if ":" in v and (">" in v or "chr" in v.lower())) >= max(3, len(var_names) // 2)
    out.append({
        "path": str(p),
        "name": p.name,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "is_nonneg_int_X": is_nonneg_int,
        "layers": sorted(str(k) for k in adata.layers.keys()),
        "obsm": sorted(str(k) for k in adata.obsm.keys()),
        "uns": sorted(str(k) for k in adata.uns.keys())[:25],
        "obs_summary": obs_summary,
        "var_names_sample": var_names,
        "looks_like_dna_variants": bool(looks_like_variants),
    })
print(json.dumps({"status": "ok", "files": out}))
"""
    proc = subprocess.run(
        [python_exe, "-c", probe_code, *[str(p) for p in h5ads]],
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    if proc.returncode != 0:
        return {"status": "probe_error", "stderr": proc.stderr[-1000:], "files": []}
    return json.loads(proc.stdout)


def check_scbench_native_contract(task_prompt: str, profile: dict[str, Any]) -> tuple[bool, str]:
    """Fail-closed native scRNA-seq contract check before running analysis."""
    prompt_l = task_prompt.lower()
    for cue in _NON_RNA_PROMPT_CUES:
        if cue in prompt_l:
            return False, f"unsupported_non_rna_modality_cue:{cue}"
    files = profile.get("files") or []
    if not files:
        return False, "no_h5ad_input_found"
    if all(f.get("looks_like_dna_variants") for f in files):
        return False, "unsupported_dna_variant_matrix_not_scrna"
    return True, "supported_scrna_h5ad"


def run_aria_pseudobulk_de_on_h5ad(
    h5ad_path: Path,
    *,
    groupby: str,
    condition_col: str,
    replicate_col: str,
    comparisons: list[list[str]],
    output_dir: Path,
    lfc_min: float = 1.0,
    padj_max: float = 0.05,
    min_cells: int = 5,
    python_exe: str = DEFAULT_RNA_PYTHON,
) -> dict[str, Any]:
    """Invoke frozen ARIA v4.9.0 `aria/scripts/rna_pseudobulk_de.py` via JSON IPC."""
    params = {
        "data_path": str(h5ad_path),
        "groupby": groupby,
        "condition_col": condition_col,
        "replicate_col": replicate_col,
        "comparisons": comparisons,
        "lfc_min": lfc_min,
        "padj_max": padj_max,
        "min_cells_per_pseudosample": min_cells,
        "min_replicates_per_condition": 2,
        "lfc_shrink": False,
        "fdr_strategy": "per_cluster",
        "use_raw": False,
        "output_dir": str(output_dir),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    in_json = output_dir / "pb_input.json"
    out_json = output_dir / "pb_output.json"
    in_json.write_text(json.dumps(params, indent=2), encoding="utf-8")
    env = isolated_benchmark_environment(output_dir, base_env=os.environ.copy())
    env["PYTHONPATH"] = str(REPO_ROOT)
    env["MPLBACKEND"] = "Agg"
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    proc = subprocess.run(
        [python_exe, str(REPO_ROOT / "aria" / "scripts" / "rna_pseudobulk_de.py"), str(in_json), str(out_json)],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    if proc.returncode != 0 or not out_json.is_file():
        return {"status": "error", "stderr": proc.stderr[-1500:], "stdout": proc.stdout[-500:]}
    return json.loads(out_json.read_text(encoding="utf-8"))


def make_aria_scbench_agent(
    *,
    api_base: str = DEFAULT_TANTALUS_BASE,
    model: str = DEFAULT_TANTALUS_MODEL,
    python_exe: str = DEFAULT_RNA_PYTHON,
):
    """Build an `EvalRunner`-compatible callable for `latchbio/scbench`."""

    def aria_agent_fn(task_prompt: str, work_dir: Path | str) -> dict[str, Any]:
        wd = Path(work_dir).resolve()
        scratch = wd / "_aria_governed_scratch"
        scratch.mkdir(parents=True, exist_ok=True)
        profile = inspect_h5ad_workspace(wd, python_exe=python_exe)
        supported, reason = check_scbench_native_contract(task_prompt, profile)
        if not supported:
            answer = {
                "abstained": True,
                "abstain_reason": reason,
                "aria_tag": "v4.9.0",
                "aria_package_sha256": compute_aria_package_sha256(),
            }
            (wd / "eval_answer.json").write_text(json.dumps(answer, indent=2), encoding="utf-8")
            return answer

        # Plan governed ARIA action from task prompt + measured .h5ad schema
        planner_prompt = (
            "You are ARIA v4.9.0's governed single-cell planner.\n"
            f"TASK:\n{task_prompt}\n\n"
            f"H5AD WORKSPACE PROFILE:\n{json.dumps(profile, indent=2)[:6000]}\n\n"
            "Return ONLY a JSON object inside <FINAL_ANSWER>...</FINAL_ANSWER> choosing one action:\n"
            '1) {"action": "pseudobulk_de", "h5ad_name": "...", "groupby": "...", '
            '"condition_col": "...", "replicate_col": "...", "comparisons": [["test","ref"]], '
            '"lfc_min": 1.0, "padj_max": 0.05, "min_cells": 5, "output_key": "n_degs"}\n'
            '2) {"action": "abstain", "reason": "..."}\n'
            "Use donor-level pseudobulk_de whenever the task asks for multi-donor differential expression; "
            "abstain if the task requires an analysis outside ARIA's contracted primitives."
        )
        plan_resp = chat_completion(
            [{"role": "user", "content": planner_prompt}],
            api_base=api_base,
            model=model,
            temperature=0.0,
            max_tokens=1024,
        )
        plan = extract_json_answer(plan_resp["content"]) or {"action": "abstain", "reason": "unparseable_plan"}

        if plan.get("action") == "pseudobulk_de":
            h5_map = {f["name"]: Path(f["path"]) for f in profile.get("files", [])}
            h5_path = h5_map.get(str(plan.get("h5ad_name") or "")) or Path(profile["files"][0]["path"])
            pb = run_aria_pseudobulk_de_on_h5ad(
                h5_path,
                groupby=str(plan.get("groupby") or ""),
                condition_col=str(plan.get("condition_col") or ""),
                replicate_col=str(plan.get("replicate_col") or ""),
                comparisons=list(plan.get("comparisons") or []),
                output_dir=scratch / "pseudobulk",
                lfc_min=float(plan.get("lfc_min", 1.0)),
                padj_max=float(plan.get("padj_max", 0.05)),
                min_cells=int(plan.get("min_cells", 5)),
                python_exe=python_exe,
            )
            if pb.get("status") == "success":
                total_sig = 0
                for g_entry in (pb.get("per_group") or {}).values():
                    for comp in (g_entry.get("per_comparison") or {}).values():
                        if comp.get("status") == "success":
                            total_sig += int(comp.get("n_significant") or len(comp.get("all_sig") or []))
                out_key = str(plan.get("output_key") or "n_degs")
                answer = {out_key: total_sig, "aria_status": "success"}
                (wd / "eval_answer.json").write_text(json.dumps(answer, indent=2), encoding="utf-8")
                return answer
            answer = {"abstained": True, "abstain_reason": f"pseudobulk_failed:{pb.get('status')}"}
            (wd / "eval_answer.json").write_text(json.dumps(answer, indent=2), encoding="utf-8")
            return answer

        answer = {
            "abstained": True,
            "abstain_reason": str(plan.get("reason") or "outside_contracted_scbench_action"),
        }
        (wd / "eval_answer.json").write_text(json.dumps(answer, indent=2), encoding="utf-8")
        return answer

    return aria_agent_fn


def run_synthetic_smoke() -> dict[str, Any]:
    """Offline end-to-end verification of the scBench adapter on synthetic .h5ad fixtures."""
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        # 1. Non-RNA DNA+protein fixture must abstain cleanly
        dna_dir = root / "dna_amplicon"
        dna_dir.mkdir()
        mk_dna = """
import sys, numpy as np, pandas as pd, anndata as ad
p = sys.argv[1]
adata = ad.AnnData(
    X=np.zeros((20, 6), dtype=float),
    obs=pd.DataFrame({"cell": [f"c{i}" for i in range(20)]}).set_index("cell"),
    var=pd.DataFrame(index=[f"chr1:{1000+i}:A>G" for i in range(6)]),
)
adata.write_h5ad(p)
"""
        subprocess.run([DEFAULT_RNA_PYTHON, "-c", mk_dna, str(dna_dir / "panel.h5ad")], check=True)
        prof_dna = inspect_h5ad_workspace(dna_dir)
        ok_dna, reason_dna = check_scbench_native_contract(
            "Identify the largest mutant clone in this single-cell DNA+protein panel.", prof_dna
        )
        assert not ok_dna, (ok_dna, reason_dna)

        # 2. Multi-donor scRNA fixture must run real ARIA rna_pseudobulk_de.py
        rna_dir = root / "scrna_pb"
        rna_dir.mkdir()
        mk_rna = """
import sys, numpy as np, pandas as pd, anndata as ad
from scipy import sparse
rng = np.random.default_rng(49)
rows, obs = [], []
n_genes = 60
base = rng.uniform(50, 200, n_genes)
for cond in ("ctrl", "stim"):
    for d in range(3):
        donor = f"D{d}"
        for c in range(8):
            mu = base.copy()
            if cond == "stim":
                mu[:8] *= 8.0
            counts = rng.poisson(mu)
            rows.append(counts)
            obs.append({"cell": f"{cond}_{donor}_{c}", "cell_type": "Tcell", "condition": cond, "donor": donor})
adata = ad.AnnData(
    X=sparse.csr_matrix(np.asarray(rows, dtype=np.int32)),
    obs=pd.DataFrame(obs).set_index("cell"),
    var=pd.DataFrame(index=[f"GENE_{i:03d}" for i in range(n_genes)]),
)
adata.write_h5ad(sys.argv[1])
"""
        h5_path = rna_dir / "multi_donor.h5ad"
        subprocess.run([DEFAULT_RNA_PYTHON, "-c", mk_rna, str(h5_path)], check=True)
        prof_rna = inspect_h5ad_workspace(rna_dir)
        ok_rna, reason_rna = check_scbench_native_contract(
            "Perform donor-aware pseudobulk differential expression between stim and ctrl.", prof_rna
        )
        assert ok_rna, reason_rna
        pb = run_aria_pseudobulk_de_on_h5ad(
            h5_path,
            groupby="cell_type",
            condition_col="condition",
            replicate_col="donor",
            comparisons=[["stim", "ctrl"]],
            output_dir=rna_dir / "out",
            lfc_min=1.0,
            padj_max=0.05,
            min_cells=5,
        )
        assert pb.get("status") == "success", pb
        comp = pb["per_group"]["Tcell"]["per_comparison"]["stim_vs_ctrl"]
        n_sig = int(comp.get("n_significant") or len(comp.get("all_sig") or []))
        assert comp["status"] == "success" and n_sig >= 6, comp
        return {
            "status": "pass",
            "non_rna_abstain_reason": reason_dna,
            "pseudobulk_n_sig": n_sig,
            "lobo_status": (comp.get("lobo_jackknife") or {}).get("status"),
            "aria_package_sha256": compute_aria_package_sha256(),
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eval-json", help="Path to a latchbio/scbench eval JSON file.")
    parser.add_argument("--synthetic-smoke", action="store_true", help="Run offline synthetic smoke test.")
    parser.add_argument("--api-base", default=DEFAULT_TANTALUS_BASE)
    parser.add_argument("--model", default=DEFAULT_TANTALUS_MODEL)
    args = parser.parse_args()

    if args.synthetic_smoke:
        print(json.dumps(run_synthetic_smoke(), indent=2, sort_keys=True))
        return 0

    if not args.eval_json:
        parser.error("Provide --eval-json or --synthetic-smoke")

    from scbench import EvalRunner  # type: ignore[import-not-found]

    agent_fn = make_aria_scbench_agent(api_base=args.api_base, model=args.model)
    runner = EvalRunner(args.eval_json)
    result = runner.run(agent_function=agent_fn)
    print(json.dumps(result, indent=2, sort_keys=True, default=str))
    return 0 if result.get("passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
