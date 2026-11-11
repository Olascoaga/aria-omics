#!/usr/bin/env python3
"""V48 Slice I3 — End-to-end WNN anti-degradation guardrail validation runner.

Executes inside ``aria-chromatin-env``:
  1. ``aria/scripts/integration_paired_qc.py`` (V48 I1 joint QC + optional
     bijective reference obs label join via ``h5py``).
  2. ``aria/scripts/integration_wnn.py`` (V48 I2 deterministic WNN + inline I3
     guardrail persisted into ``mdata.uns['aria_wnn_diagnostics']``).
  3. Two seeded runs of ``aria/scripts/integration_wnn_diagnostics.py`` (V48 I3
     standalone IPC script) verifying bit-exact determinism and full 5-axis
     anti-degradation diagnostics across both unsupervised per-modality rare
     clusters and optional external reference rare groups.

Emits a sanitized JSON artifact with zero machine-absolute paths (enforced by
``tests/test_public_artifacts_no_absolute_paths.py``).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _resolve_conda_bin() -> str:
    found = shutil.which("conda")
    if found:
        return found
    candidate = Path(sys.executable).resolve().parents[2] / "bin" / "conda"
    if candidate.is_file():
        return str(candidate)
    return "conda"


def _run_ipc_script(
    env: str,
    script_rel: str,
    params: dict[str, Any],
    work_dir: Path,
    tag: str,
) -> dict[str, Any]:
    in_json = work_dir / f"{tag}_in.json"
    out_json = work_dir / f"{tag}_out.json"
    in_json.write_text(json.dumps(params, indent=2), encoding="utf-8")
    proc = subprocess.run(
        [
            _resolve_conda_bin(),
            "run",
            "--no-capture-output",
            "-n",
            env,
            "python",
            str(ROOT / script_rel),
            str(in_json),
            str(out_json),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not out_json.is_file():
        raise RuntimeError(
            f"{script_rel} exited with {proc.returncode}:\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(out_json.read_text(encoding="utf-8"))


def _sanitize_public_value(val: Any) -> Any:
    """Strip machine-absolute filesystem prefixes from public manifest values."""
    if isinstance(val, dict):
        return {str(k): _sanitize_public_value(v) for k, v in val.items()}
    if isinstance(val, list):
        return [_sanitize_public_value(v) for v in val]
    if isinstance(val, str):
        p = Path(val)
        if val.startswith(
            ("/home/", "/Users/", "/root/", "/tmp/", "/mnt/", "/private/", "/data/")
        ):
            try:
                return str(p.resolve().relative_to(ROOT.resolve()))
            except Exception:
                return p.name
        return val
    return val


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--h5mu",
        required=True,
        help="Input paired .h5mu (raw counts or pre-QC paired MuData).",
    )
    ap.add_argument(
        "--reference-h5ad",
        default="",
        help="Optional reference .h5ad for bijective /obs label join in I1.",
    )
    ap.add_argument(
        "--audit-group-col",
        default="",
        help="Optional obs column to audit rare reference groups in I3.",
    )
    ap.add_argument(
        "--sample-id",
        default="",
        help="Optional sample prefix for bijective reference barcode join.",
    )
    ap.add_argument(
        "--reference-sample-col",
        default="",
        help="Optional reference sample column for bijective join.",
    )
    ap.add_argument(
        "--keep-raw-group-names",
        action="store_true",
        help="Keep raw audit group names instead of neutral rank IDs in JSON.",
    )
    ap.add_argument("--chromatin-env", default="aria-chromatin-env")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--k", type=int, default=20)
    ap.add_argument("--resolution", type=float, default=0.5)
    ap.add_argument("--skip-paired-qc", action="store_true")
    ap.add_argument("--work-dir", default="")
    ap.add_argument(
        "--output-json",
        default=str(
            ROOT
            / "docs"
            / "benchmark_results"
            / "multiome_concordance"
            / "v48_i3_wnn_guardrail_hc11.json"
        ),
    )
    args = ap.parse_args()

    cleanup_tmp: tempfile.TemporaryDirectory[str] | None = None
    if args.work_dir:
        work_dir = Path(args.work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
    else:
        cleanup_tmp = tempfile.TemporaryDirectory(prefix="v48_i3_guardrail_")
        work_dir = Path(cleanup_tmp.name)

    try:
        input_h5mu = Path(args.h5mu)
        qc_summary: dict[str, Any] | None = None
        if args.skip_paired_qc:
            wnn_input_h5mu = input_h5mu
        else:
            qc_dir = work_dir / "i1_qc"
            qc_params: dict[str, Any] = {
                "h5mu_path": str(input_h5mu),
                "output_dir": str(qc_dir),
            }
            if args.sample_id:
                qc_params["sample_id"] = str(args.sample_id)
            if args.reference_h5ad and args.audit_group_col:
                qc_params["reference_h5ad"] = str(args.reference_h5ad)
                qc_params["reference_label_col"] = str(args.audit_group_col)
                if args.reference_sample_col:
                    qc_params["reference_sample_col"] = str(args.reference_sample_col)
            qc_res = _run_ipc_script(
                args.chromatin_env,
                "aria/scripts/integration_paired_qc.py",
                qc_params,
                work_dir=work_dir,
                tag="i1_qc",
            )
            if qc_res.get("status") != "success":
                raise RuntimeError(f"I1 paired QC failed: {qc_res}")
            wnn_input_h5mu = Path(str(qc_res["output_path"]))
            lj = dict(qc_res.get("label_join") or {})
            if not args.keep_raw_group_names and isinstance(lj.get("label_counts"), dict):
                raw_counts = lj.pop("label_counts")
                lj["n_unique_labels"] = len(raw_counts)
                lj["label_size_distribution"] = sorted(
                    int(v) for v in raw_counts.values()
                )
            qc_summary = {
                "n_cells_before": qc_res["n_cells_before"],
                "n_cells_after": qc_res["n_cells_after"],
                "n_genes": qc_res["n_genes"],
                "n_peaks": qc_res["n_peaks"],
                "label_join": lj,
            }

        wnn_dir = work_dir / "i2_wnn"
        wnn_params: dict[str, Any] = {
            "h5mu_path": str(wnn_input_h5mu),
            "output_dir": str(wnn_dir),
            "k": int(args.k),
            "resolution": float(args.resolution),
            "seed": int(args.seed),
        }
        if args.audit_group_col:
            wnn_params["audit_group_col"] = str(args.audit_group_col)

        wnn_res = _run_ipc_script(
            args.chromatin_env,
            "aria/scripts/integration_wnn.py",
            wnn_params,
            work_dir=work_dir,
            tag="i2_wnn",
        )
        if wnn_res.get("status") != "success":
            raise RuntimeError(f"I2 WNN failed: {wnn_res}")

        wnn_h5mu = Path(str(wnn_res["output_path"]))
        diag_params: dict[str, Any] = {
            "h5mu_path": str(wnn_h5mu),
            "seed": int(args.seed),
            "k": int(args.k),
        }
        if args.audit_group_col:
            diag_params["audit_group_col"] = str(args.audit_group_col)

        diag_run1 = _run_ipc_script(
            args.chromatin_env,
            "aria/scripts/integration_wnn_diagnostics.py",
            diag_params,
            work_dir=work_dir,
            tag="i3_diag_run1",
        )
        diag_run2 = _run_ipc_script(
            args.chromatin_env,
            "aria/scripts/integration_wnn_diagnostics.py",
            diag_params,
            work_dir=work_dir,
            tag="i3_diag_run2",
        )
        if diag_run1.get("status") != "success" or diag_run2.get("status") != "success":
            raise RuntimeError(f"I3 diagnostics failed: {diag_run1}")

        bit_exact_deterministic = (
            json.dumps(diag_run1["diagnostics"], sort_keys=True)
            == json.dumps(diag_run2["diagnostics"], sort_keys=True)
        )
        inline_matches_standalone = (
            json.dumps(wnn_res["anti_degradation"], sort_keys=True)
            == json.dumps(diag_run1["diagnostics"], sort_keys=True)
        )

        diag_out = json.loads(json.dumps(diag_run1["diagnostics"]))
        if not args.keep_raw_group_names:
            ref_block = diag_out.get("rare_population_retention", {}).get(
                "reference_rare_groups"
            )
            if isinstance(ref_block, dict) and isinstance(ref_block.get("groups"), list):
                groups_sorted = sorted(
                    ref_block["groups"],
                    key=lambda g: (int(g.get("n_cells", 0)), str(g.get("group_id", ""))),
                )
                for rank_idx, grp_entry in enumerate(groups_sorted, start=1):
                    grp_entry["group_id"] = (
                        f"ref_rare_{rank_idx:02d}_n{int(grp_entry['n_cells'])}"
                    )
                ref_block["groups"] = groups_sorted

        payload = _sanitize_public_value(
            {
                "slice": "V48_I3_wnn_anti_degradation_guardrail",
                "principle": "lens_not_replacement",
                "dataset_artifact": input_h5mu.name,
                "chromatin_env": args.chromatin_env,
                "seed": int(args.seed),
                "k": int(args.k),
                "resolution": float(args.resolution),
                "paired_qc_summary": qc_summary,
                "determinism": {
                    "standalone_repeat_bit_exact": bit_exact_deterministic,
                    "inline_vs_standalone_bit_exact": inline_matches_standalone,
                },
                "guardrail_status": diag_run1["guardrail_status"],
                "issue_codes": diag_run1["issue_codes"],
                "issues": diag_run1["issues"],
                "diagnostics": diag_out,
            }
        )

        out_path = Path(args.output_json)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(
            json.dumps(
                {
                    "status": "ok",
                    "guardrail_status": payload["guardrail_status"],
                    "issue_codes": payload["issue_codes"],
                    "determinism": payload["determinism"],
                    "output_json": _sanitize_public_value(str(out_path)),
                }
            )
        )
        return 0
    finally:
        if cleanup_tmp is not None:
            cleanup_tmp.cleanup()


if __name__ == "__main__":
    raise SystemExit(main())
