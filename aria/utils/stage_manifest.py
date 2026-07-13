"""Per-stage content-addressed resume manifests (A4).

Pipeline resume must not reuse a stage's output when its inputs, parameters,
reference, or tool version changed. The FASTQ QC / alignment resume gates used to
check only output existence and integrity, so a changed input FASTQ, GTF /
reference, tool version, or parameter silently reused a stale output.

Each stage now writes a small manifest recording a content signature of its
inputs plus an order-invariant params hash and the tool version. The resume gate
reuses a stage only when :func:`stage_is_current` confirms the manifest still
matches; any changed input / param / reference / version selectively invalidates
that stage (and, through the DAG, its downstream consumers whose input — the
upstream output — then changes too).

Signature strategy (Samael, 2026-07-12): hybrid content-addressing. Every input
carries a content sha256, but on resume an unchanged file is recognised by its
``(size, mtime_ns)`` fast path and is NOT re-hashed. Only a size change, or a
same-size ``mtime`` change whose re-hash actually differs, invalidates — so a
multi-GB sequencing file is not re-hashed on every resume check while a real edit
is still caught.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from aria.utils.provenance import hash_file, hash_params

SCHEMA = "aria.stage_manifest.v1"


def _file_signature(path: str | Path) -> dict[str, Any]:
    st = os.stat(path)
    return {
        "path": str(path),
        "size": int(st.st_size),
        "mtime_ns": int(st.st_mtime_ns),
        "sha256": hash_file(path),
    }


def build_stage_manifest(
    *,
    stage: str,
    inputs: list[tuple[str, str | None]],
    params: dict | None,
    tool_version: str | None,
    extra: dict | None = None,
) -> dict[str, Any]:
    """Assemble a stage manifest.

    ``inputs`` is a list of ``(role, path)`` pairs; a ``None`` path (e.g. the R2
    of a single-end sample) is skipped. Each present input is hashed once here.
    """
    sig_inputs = []
    for role, path in inputs:
        if path is None:
            continue
        signature = _file_signature(path)
        signature["role"] = role
        sig_inputs.append(signature)
    manifest = {
        "schema": SCHEMA,
        "stage": str(stage),
        "tool_version": str(tool_version or ""),
        "params_sha256": hash_params(params or {}),
        "inputs": sig_inputs,
    }
    if extra:
        manifest["extra"] = extra
    return manifest


def write_stage_manifest(manifest_path: str | Path, manifest: dict) -> None:
    """Atomically write a stage manifest (serialize → tmp+fsync → os.replace)."""
    path = Path(manifest_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, sort_keys=True)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _input_matches(recorded: dict, current_path: str) -> tuple[bool, str]:
    if not current_path or not os.path.exists(current_path):
        return False, f"input_missing:{current_path}"
    st = os.stat(current_path)
    if int(st.st_size) != int(recorded.get("size", -1)):
        return False, f"input_changed:{current_path}"
    if int(st.st_mtime_ns) == int(recorded.get("mtime_ns", -1)):
        return True, ""  # fast path: size + mtime unchanged, trust the hash
    # Same size, different mtime → re-hash to distinguish a touch from an edit.
    if hash_file(current_path) != recorded.get("sha256"):
        return False, f"input_changed:{current_path}"
    return True, ""


def stage_is_current(
    manifest_path: str | Path,
    *,
    inputs: list[tuple[str, str | None]],
    params: dict | None,
    tool_version: str | None,
) -> tuple[bool, str]:
    """Return ``(current, reason)`` for a stage's resume decision.

    ``current`` is True only when the stored manifest exists and its tool
    version, params hash, and every input's content signature still match the
    supplied stage inputs. A missing/corrupt manifest, changed tool version,
    changed param, changed/missing input, or a changed input set → ``False`` with
    a concrete reason.
    """
    path = Path(manifest_path)
    if not path.exists():
        return False, "no_manifest"
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False, "corrupt_manifest"

    if str(manifest.get("tool_version") or "") != str(tool_version or ""):
        return False, "tool_version_changed"
    if manifest.get("params_sha256") != hash_params(params or {}):
        return False, "params_changed"

    recorded_by_role = {s.get("role"): s for s in manifest.get("inputs", [])}
    provided = [(role, p) for role, p in inputs if p is not None]
    if set(recorded_by_role) != {role for role, _ in provided}:
        return False, "input_set_changed"

    for role, current_path in provided:
        ok, why = _input_matches(recorded_by_role[role], current_path)
        if not ok:
            return False, why
    return True, "current"
