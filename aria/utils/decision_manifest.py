"""Committed Pre-Execution DecisionManifest & Analytical Forking-Path Ledger (P0).

Prevents silent agentic HARKing and cross-provider plan drift by:
1. Serializing and SHA-256 committing the full pre-outcome analytical
   specification (`DecisionManifest`) BEFORE any modality compute stage runs.
2. Recording an append-only per-study `AnalyticalForkLedger` (`fork_ledger.jsonl`)
   keyed by input content fingerprint so every parameter/contrast/covariate
   branch increments an explicit, auditable fork counter while pure deterministic
   replays keep `fork_index` invariant.
3. Supporting zero-LLM replay (`replay_decision_manifest`) where `OrchestratorAgent`
   and `ParameterAdvisor` consume the committed manifest without issuing LLM
   planning or justification calls.
"""

from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path, PureWindowsPath
from typing import Any

from aria.utils.provenance import hash_params

DECISION_MANIFEST_SCHEMA = "aria.decision_manifest.v1"
FORK_LEDGER_SCHEMA = "aria.fork_ledger.v1"


def _is_private_absolute_path(value: str) -> bool:
    text = str(value or "")
    if not text:
        return False
    try:
        return (
            Path(text).is_absolute()
            or PureWindowsPath(text).is_absolute()
            or text.startswith("~/")
        )
    except (OSError, ValueError):
        return False


def _portable_path_ref(value: str, content_sha256: Any = None) -> str:
    digest = str(content_sha256 or "").strip()
    if digest and digest != "unavailable":
        return f"input://sha256/{digest}"
    path_digest = hashlib.sha256(str(value).encode("utf-8")).hexdigest()
    return f"private-path://sha256/{path_digest}"


def _sanitize_json_value(value: Any) -> Any:
    """Recursively replace non-finite floats with None for strict RFC-8259 JSON."""
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {k: _sanitize_json_value(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_sanitize_json_value(v) for v in value]
    if isinstance(value, tuple):
        return [_sanitize_json_value(v) for v in value]
    if isinstance(value, set):
        return [_sanitize_json_value(v) for v in sorted(value, key=str)]
    item = getattr(value, "item", None)
    if callable(item) and type(value).__module__.startswith("numpy"):
        try:
            return _sanitize_json_value(item())
        except Exception:
            pass
    return value


def _sanitize_portable_value(value: Any, *, content_sha256: Any = None) -> Any:
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    item = getattr(value, "item", None)
    if callable(item) and type(value).__module__.startswith("numpy"):
        try:
            return _sanitize_portable_value(item(), content_sha256=content_sha256)
        except Exception:
            pass
    if isinstance(value, dict):
        item_sha = value.get("sha256")
        return {
            key: _sanitize_portable_value(
                item,
                content_sha256=item_sha if key == "path" else None,
            )
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_sanitize_portable_value(item) for item in value]
    if isinstance(value, tuple):
        return [_sanitize_portable_value(item) for item in value]
    if isinstance(value, Path):
        value = str(value)
    if isinstance(value, str) and _is_private_absolute_path(value):
        return _portable_path_ref(value, content_sha256)
    return value


def _sanitize_free_text_paths(text: str, ctx: dict[str, Any]) -> str:
    out = str(text or "")
    if not out:
        return out
    replacements: dict[str, str] = {}
    for rec in ctx.get("input_files") or []:
        if isinstance(rec, dict) and rec.get("path"):
            raw = str(rec["path"])
            if _is_private_absolute_path(raw):
                replacements[raw] = _portable_path_ref(raw, rec.get("sha256"))
    for paths in (ctx.get("modalities") or {}).values():
        if isinstance(paths, list):
            for raw_p in paths:
                raw = str(raw_p)
                if _is_private_absolute_path(raw) and raw not in replacements:
                    replacements[raw] = _portable_path_ref(raw)
    data_dir = ctx.get("data_dir")
    if isinstance(data_dir, (str, Path)):
        raw_dir = str(data_dir)
        if _is_private_absolute_path(raw_dir) and raw_dir not in replacements:
            replacements[raw_dir] = _portable_path_ref(raw_dir)
    for private, public in sorted(
        replacements.items(), key=lambda kv: len(kv[0]), reverse=True
    ):
        out = out.replace(private, public)
    if _is_private_absolute_path(out):
        return _portable_path_ref(out)
    return out


def _threshold_float(primary: Any, fallback: Any, default: float) -> float:
    for candidate in (primary, fallback, default):
        if candidate is None or isinstance(candidate, bool):
            continue
        try:
            val = float(candidate)
        except (TypeError, ValueError):
            continue
        if math.isfinite(val):
            return val
    return float(default)


@contextmanager
def _advisory_lock(lock_path: Path):
    fh = None
    try:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        fh = open(lock_path, "a+", encoding="utf-8")
        try:
            import fcntl
            fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
        except (ImportError, OSError):
            pass
        yield
    finally:
        if fh is not None:
            try:
                import fcntl
                fcntl.flock(fh.fileno(), fcntl.LOCK_UN)
            except (ImportError, OSError):
                pass
            fh.close()


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(
            _sanitize_json_value(payload),
            fh,
            indent=2,
            sort_keys=True,
            default=str,
            allow_nan=False,
        )
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def compute_study_input_fingerprint(exp_context: dict[str, Any]) -> str:
    """Compute a path-portable fingerprint of the study's biological input files."""
    input_records = exp_context.get("input_files") or []
    canonical_inputs: list[dict[str, Any]] = []
    if isinstance(input_records, list) and input_records:
        for rec in input_records:
            if not isinstance(rec, dict):
                continue
            sha = str(rec.get("sha256") or "")
            raw_path = str(rec.get("path") or "")
            basename = Path(raw_path).name if raw_path else ""
            canonical_inputs.append({
                "modality": str(rec.get("modality") or "unknown"),
                "sha256": sha if sha and sha != "unavailable" else "",
                "basename": basename if not sha or sha == "unavailable" else "",
                "size_bytes": int(rec.get("size_bytes") or 0) if not sha or sha == "unavailable" else 0,
            })
    if not canonical_inputs:
        modalities = exp_context.get("modalities") or {}
        if isinstance(modalities, dict):
            for mod in sorted(modalities):
                for raw_path in sorted(str(p) for p in (modalities.get(mod) or [])):
                    p = Path(raw_path)
                    size = int(p.stat().st_size) if p.exists() and p.is_file() else 0
                    canonical_inputs.append({
                        "modality": str(mod),
                        "sha256": "",
                        "basename": p.name,
                        "size_bytes": size,
                    })
    canonical_inputs.sort(
        key=lambda item: (
            item["modality"],
            item["sha256"],
            item["basename"],
            item["size_bytes"],
        )
    )
    return hash_params({"inputs": canonical_inputs})


def _normalize_design(design: dict[str, Any] | None) -> dict[str, Any]:
    d = dict(design or {})
    contrasts = d.get("plan_contrasts") or d.get("contrasts") or d.get("comparisons") or []
    norm_contrasts = []
    for c in contrasts:
        if isinstance(c, dict):
            norm_contrasts.append({
                "numerator": str(
                    c.get("numerator")
                    or c.get("test")
                    or c.get("case")
                    or c.get("treatment")
                    or c.get("group1")
                    or ""
                ),
                "denominator": str(
                    c.get("denominator")
                    or c.get("reference")
                    or c.get("baseline")
                    or c.get("control")
                    or c.get("group2")
                    or ""
                ),
                "factor": str(c.get("factor") or d.get("condition_col") or d.get("main_factor") or d.get("factor") or ""),
            })
        elif isinstance(c, (list, tuple)) and len(c) >= 2:
            norm_contrasts.append({
                "numerator": str(c[0]),
                "denominator": str(c[1]),
                "factor": str(d.get("condition_col") or d.get("main_factor") or d.get("factor") or ""),
            })
        else:
            norm_contrasts.append({"contrast": str(c)})
    groups = d.get("groups")
    norm_groups = (
        {str(k): sorted(str(s) for s in (v or [])) for k, v in sorted(groups.items())}
        if isinstance(groups, dict)
        else {}
    )
    analysis_groups = d.get("analysis_groups")
    norm_analysis_groups = (
        {str(k): sorted(str(s) for s in (v or [])) for k, v in sorted(analysis_groups.items())}
        if isinstance(analysis_groups, dict)
        else {}
    )
    replicates = d.get("replicates")
    norm_replicates = (
        {str(k): int(v) for k, v in sorted(replicates.items())}
        if isinstance(replicates, dict)
        else {}
    )
    return {
        "organism": str(d.get("organism") or ""),
        "genome": str(d.get("genome") or ""),
        "main_factor": str(d.get("main_factor") or d.get("factor") or d.get("condition_col") or ""),
        "factor": str(d.get("factor") or d.get("main_factor") or d.get("condition_col") or ""),
        "formula": str(d.get("formula") or ""),
        "batch_col": str(d.get("batch_col") or d.get("batch_covariate") or d.get("batch") or ""),
        "donor_col": str(d.get("donor_col") or d.get("biological_unit_col") or ""),
        "covariates": sorted(str(x) for x in (d.get("covariates") or []) if x),
        "paired": bool(d.get("paired", False)),
        "groups": norm_groups,
        "analysis_groups": norm_analysis_groups,
        "replicates": norm_replicates,
        "replicate_handling": dict(d.get("replicate_handling") or {}),
        "contrasts": norm_contrasts,
        "causal_estimands": list(d.get("causal_estimands") or []),
    }


def _normalize_plan(plan: dict[str, Any] | None, exp_context: dict[str, Any]) -> dict[str, Any]:
    p = dict(plan or {})
    steps = []
    for step in p.get("steps") or []:
        if isinstance(step, dict):
            steps.append({
                "order": step.get("order"),
                "agent": str(step.get("agent") or ""),
                "analysis": str(step.get("analysis") or ""),
                "depends_on": list(step.get("depends_on") or []),
            })
    return {
        "steps": steps,
        "contrasts": list(p.get("contrasts") or []),
        "integration_needed": bool(p.get("integration_needed", False)),
        "integration_type": str(p.get("integration_type") or "none"),
        "design_intelligence_selection": str(
            exp_context.get("design_intelligence_selection") or "recommended_only"
        ),
        "run_optional_supported": bool(exp_context.get("run_optional_supported", False)),
        "enable_hypotheses": bool(exp_context.get("enable_hypotheses", False)),
        "reproducible_mode": bool(exp_context.get("reproducible_mode", False)),
    }


def _diff_specification_fields(prev_spec: dict[str, Any], cur_spec: dict[str, Any]) -> list[str]:
    changed: list[str] = []
    for top_key in sorted(set(prev_spec) | set(cur_spec)):
        pv = prev_spec.get(top_key)
        cv = cur_spec.get(top_key)
        if pv == cv:
            continue
        if isinstance(pv, dict) and isinstance(cv, dict):
            for sub_key in sorted(set(pv) | set(cv)):
                if pv.get(sub_key) != cv.get(sub_key):
                    changed.append(f"{top_key}.{sub_key}")
        else:
            changed.append(top_key)
    return changed


def build_decision_manifest(
    *,
    experiment_id: str,
    exp_context: dict[str, Any],
    plan: dict[str, Any] | None,
    biological_intent: dict[str, Any] | None,
    commitment_status: str = "locked_before_execution",
    locked_parameter_decisions: dict[str, Any] | None = None,
    replayed_from_manifest_sha256: str | None = None,
) -> dict[str, Any]:
    """Build a canonical pre-execution `DecisionManifest` dictionary."""
    ctx = exp_context or {}
    intent = dict(biological_intent or {})
    modalities_map = ctx.get("modalities") or {}
    modalities_summary = {
        str(k): len(v) if isinstance(v, list) else 1
        for k, v in sorted(modalities_map.items())
    }
    reproducible = bool(ctx.get("reproducible_mode", False))
    raw_q = _sanitize_free_text_paths(
        str(
            intent.get("raw_user_question")
            or intent.get("user_question")
            or ctx.get("user_question")
            or ""
        ),
        ctx,
    )
    summary_text = _sanitize_free_text_paths(
        str(intent.get("summary") or raw_q[:200]),
        ctx,
    )

    core_spec = {
        "study_input_fingerprint": compute_study_input_fingerprint(ctx),
        "organism": str(ctx.get("organism") or ""),
        "genome": str(ctx.get("genome") or ""),
        "modalities": modalities_summary,
        "design": _sanitize_portable_value(_normalize_design(ctx.get("design"))),
        "thresholds": {
            "global_padj": _threshold_float(
                ctx.get("global_padj"), (plan or {}).get("padj_threshold"), 0.05
            ),
            "global_lfc": _threshold_float(
                ctx.get("global_lfc"), (plan or {}).get("lfc_threshold"), 1.0
            ),
            "lfc_threshold_provenance": str(
                ctx.get("lfc_threshold_provenance")
                or (plan or {}).get("lfc_threshold_provenance")
                or "default"
            ),
        },
        "analysis_plan": _sanitize_portable_value(_normalize_plan(plan, ctx)),
        "parameter_grid_policy": {
            "leiden_grid": "canonical_union_12pt" if reproducible else "intent_window_4pt",
            "historical_bonus_disabled": reproducible,
        },
        "locked_parameter_decisions": _sanitize_portable_value(
            dict(sorted((locked_parameter_decisions or {}).items()))
        ),
    }
    specification_sha256 = hash_params(core_spec)

    committed_ts = (
        "<timestamp redacted for byte-identity>"
        if reproducible
        else datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    )
    manifest_body = {
        "schema_version": DECISION_MANIFEST_SCHEMA,
        "experiment_id": str(experiment_id),
        "commitment_status": str(commitment_status),
        "committed_before_outcome_compute": True,
        "committed_at_utc": committed_ts,
        "specification_sha256": specification_sha256,
        "replayed_from_manifest_sha256": replayed_from_manifest_sha256,
        "biological_intent": {
            "analysis_type": str(intent.get("analysis_type") or "differential"),
            "comparison": _sanitize_free_text_paths(str(intent.get("comparison") or ""), ctx),
            "raw_user_question": raw_q,
            "summary": summary_text,
        },
        **core_spec,
    }
    hashable_body = {
        k: v for k, v in manifest_body.items()
        if k not in {"experiment_id", "committed_at_utc"}
    }
    manifest_body["manifest_sha256"] = hash_params(hashable_body)
    return manifest_body


def verify_decision_manifest(manifest: dict[str, Any]) -> tuple[bool, str]:
    """Verify cryptographic integrity of a serialized `DecisionManifest`."""
    if not isinstance(manifest, dict):
        return False, "manifest_not_a_dict"
    if manifest.get("schema_version") != DECISION_MANIFEST_SCHEMA:
        return False, f"unsupported_schema:{manifest.get('schema_version')}"
    recorded_manifest_sha = str(manifest.get("manifest_sha256") or "")
    if not recorded_manifest_sha:
        return False, "missing_manifest_sha256"
    body = {
        k: v for k, v in manifest.items()
        if k not in {
            "manifest_sha256",
            "experiment_id",
            "committed_at_utc",
            "realized_parameter_decisions",
            "fork_ledger_summary",
        }
    }
    expected_manifest_sha = hash_params(body)
    if recorded_manifest_sha != expected_manifest_sha:
        return False, "manifest_sha256_mismatch"

    core_keys = (
        "study_input_fingerprint",
        "organism",
        "genome",
        "modalities",
        "design",
        "thresholds",
        "analysis_plan",
        "parameter_grid_policy",
        "locked_parameter_decisions",
    )
    core_spec = {k: manifest.get(k) for k in core_keys}
    if hash_params(core_spec) != manifest.get("specification_sha256"):
        return False, "specification_sha256_mismatch"
    return True, "verified"


def load_decision_manifest(source: str | Path | dict[str, Any]) -> dict[str, Any]:
    """Load and verify a `DecisionManifest` from dict, `decision_manifest.json`, or `methodology.json`."""
    if isinstance(source, dict):
        payload = dict(source)
    else:
        p = Path(source).expanduser().resolve()
        if p.is_dir():
            dm_file = p / "decision_manifest.json"
            meth_file = p / "methodology.json"
            target = dm_file if dm_file.is_file() else meth_file
        else:
            target = p
        raw = json.loads(target.read_text(encoding="utf-8"))
        if isinstance(raw, dict) and raw.get("schema_version") == DECISION_MANIFEST_SCHEMA:
            payload = raw
        elif isinstance(raw, dict) and isinstance(raw.get("decision_manifest"), dict):
            payload = dict(raw["decision_manifest"])
        else:
            raise ValueError(f"No valid DecisionManifest found in {target}")

    ok, reason = verify_decision_manifest(payload)
    if not ok:
        raise ValueError(f"DecisionManifest integrity check failed: {reason}")
    return payload


def resolve_fork_ledger_path(memory_db_path: str | None = None) -> Path:
    """Resolve the append-only fork ledger JSONL path, colocated with scoped `memory.db`."""
    explicit = os.environ.get("ARIA_FORK_LEDGER_PATH", "").strip()
    if explicit:
        return Path(explicit).expanduser()
    if memory_db_path and memory_db_path != ":memory:":
        db_p = Path(memory_db_path).expanduser()
        return db_p.parent / f"{db_p.stem}_fork_ledger.jsonl"
    env_db = os.environ.get("ARIA_MEMORY_DB", "").strip()
    if env_db and env_db != ":memory:":
        db_p = Path(env_db).expanduser()
        return db_p.parent / f"{db_p.stem}_fork_ledger.jsonl"
    return Path.home() / ".aria" / "fork_ledger.jsonl"


def commit_decision_manifest(
    manifest: dict[str, Any],
    *,
    workspace_dir: str | Path | None = None,
    ledger_path: str | Path | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Atomically persist `decision_manifest.json` and register it in `fork_ledger.jsonl`."""
    ok, reason = verify_decision_manifest(manifest)
    if not ok:
        raise ValueError(f"Refusing to commit invalid DecisionManifest: {reason}")

    lpath = Path(ledger_path) if ledger_path is not None else resolve_fork_ledger_path()
    lock_path = lpath.with_name(f".{lpath.name}.lock")
    study_fp = manifest["study_input_fingerprint"]
    spec_sha = manifest["specification_sha256"]

    core_keys = (
        "study_input_fingerprint",
        "organism",
        "genome",
        "modalities",
        "design",
        "thresholds",
        "analysis_plan",
        "parameter_grid_policy",
        "locked_parameter_decisions",
    )
    cur_core = {k: manifest.get(k) for k in core_keys}

    with _advisory_lock(lock_path):
        prior_study_entries: list[dict[str, Any]] = []
        if lpath.is_file():
            for line in lpath.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                if isinstance(rec, dict) and rec.get("study_input_fingerprint") == study_fp:
                    prior_study_entries.append(rec)

        distinct_specs: dict[str, dict[str, Any]] = {}
        for entry in prior_study_entries:
            s_sha = str(entry.get("specification_sha256") or "")
            if s_sha and s_sha not in distinct_specs:
                distinct_specs[s_sha] = entry

        if spec_sha in distinct_specs:
            matched = distinct_specs[spec_sha]
            fork_index = int(matched.get("fork_index", 0))
            parent_manifest_sha = matched.get("parent_manifest_sha256")
            fork_changes: list[str] = list(matched.get("fork_changes") or [])
            is_exact_spec_replay = True
            total_distinct_specs = len(distinct_specs)
        elif not distinct_specs:
            fork_index = 0
            parent_manifest_sha = None
            fork_changes = []
            is_exact_spec_replay = False
            total_distinct_specs = 1
        else:
            fork_index = len(distinct_specs)
            latest_prev = prior_study_entries[-1]
            parent_manifest_sha = latest_prev.get("manifest_sha256")
            prev_core = latest_prev.get("core_specification") or {}
            fork_changes = (
                _diff_specification_fields(prev_core, cur_core)
                if isinstance(prev_core, dict) and prev_core
                else ["analytical_specification"]
            )
            is_exact_spec_replay = False
            total_distinct_specs = len(distinct_specs) + 1

        effective_status = str(manifest.get("commitment_status") or "locked_before_execution")
        if is_exact_spec_replay:
            matched_status = str(matched.get("commitment_status") or "")
            if fork_index > 0 and effective_status == "locked_before_execution" and matched_status:
                effective_status = matched_status
        elif fork_index > 0 and effective_status == "locked_before_execution":
            effective_status = "post_hoc_fork"

        runs_observed = len(prior_study_entries) + 1
        ledger_record = {
            "schema_version": FORK_LEDGER_SCHEMA,
            "recorded_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "experiment_id": manifest["experiment_id"],
            "study_input_fingerprint": study_fp,
            "specification_sha256": spec_sha,
            "manifest_sha256": manifest["manifest_sha256"],
            "commitment_status": effective_status,
            "fork_index": fork_index,
            "total_analytical_forks": total_distinct_specs,
            "is_exact_spec_replay": is_exact_spec_replay,
            "parent_manifest_sha256": parent_manifest_sha,
            "fork_changes": fork_changes,
            "ledger_runs_observed": runs_observed,
            "core_specification": cur_core,
        }
        lpath.parent.mkdir(parents=True, exist_ok=True)
        with open(lpath, "a", encoding="utf-8") as fh:
            fh.write(
                json.dumps(
                    _sanitize_json_value(ledger_record),
                    sort_keys=True,
                    default=str,
                    allow_nan=False,
                )
                + "\n"
            )
            fh.flush()
            os.fsync(fh.fileno())

    reproducible = bool(
        (manifest.get("analysis_plan") or {}).get("reproducible_mode", False)
    )
    fork_summary = {
        "schema_version": FORK_LEDGER_SCHEMA,
        "study_input_fingerprint": study_fp,
        "specification_sha256": spec_sha,
        "manifest_sha256": manifest["manifest_sha256"],
        "commitment_status": effective_status,
        "fork_index": fork_index,
        "total_analytical_forks": total_distinct_specs,
        "is_exact_spec_replay": is_exact_spec_replay,
        "parent_manifest_sha256": parent_manifest_sha,
        "fork_changes": fork_changes,
        "ledger_runs_observed": runs_observed,
    }
    manifest_fork_summary = dict(fork_summary)
    if reproducible:
        # Keep sequential replay telemetry in the external append-only fork
        # ledger and returned fork_summary, but omit run-count / replay-order
        # fields inside decision_manifest.json so Run 1 and exact Replay Run 2
        # remain byte-identical under reproducible_mode=True.
        manifest_fork_summary.pop("ledger_runs_observed", None)
        manifest_fork_summary.pop("is_exact_spec_replay", None)
    enriched_manifest = {
        **manifest,
        "fork_ledger_summary": manifest_fork_summary,
    }

    if workspace_dir is not None:
        ws_path = Path(workspace_dir) / "decision_manifest.json"
        _atomic_write_json(ws_path, enriched_manifest)

    return enriched_manifest, fork_summary
