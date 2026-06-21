"""W-LEDGER export: a machine-readable reproducible capsule per report + `aria diff`.

The report's ``methodology.json`` already carries the full run evidence graph:
provenance (version/commit/dirty/workflow_hash/seeds), input hashes, the claims
manifest (each claim linked to an evidence card AND a run-ledger node), the run
ledger, devil's-advocate, robustness, calibration, tools, and decisions. This
module serializes that into:

- an **RO-Crate 1.1** metadata file (`ro-crate-metadata.json`, W3C-PROV-flavored
  JSON-LD) describing the run as a provenance graph — DOI/repository ready;
- a **reproducible capsule** ZIP bundling the report directory + the crate;
- `aria diff reportA reportB` — a structured comparison of two run ledgers
  (provenance, executed analyses, claims, calibration).

Pure bookkeeping over structured JSON; no LLM, no biology, no network.
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path
from typing import Any

from aria.utils.reference_integrity import sha256_file

RO_CRATE_CONTEXT = "https://w3id.org/ro/crate/1.1/context"
_METHODOLOGY_NAME = "methodology.json"
_CAPSULE_MANIFEST = "capsule_manifest.json"
_CAPSULE_SCHEMA_VERSION = "s14.reproducibility_capsule.v1"


# ── loading ──────────────────────────────────────────────────────────────────

def load_methodology(source: str | Path) -> dict:
    """Load a methodology dict from a methodology.json path or a report dir."""
    p = Path(source)
    if p.is_dir():
        p = p / _METHODOLOGY_NAME
    return json.loads(p.read_text())


def _provenance(methodology: dict) -> dict:
    return (methodology or {}).get("provenance", {}) or {}


def _version(methodology: dict) -> str:
    prov = _provenance(methodology)
    return str(prov.get("version") or prov.get("aria_version") or "unknown")


def _commit(methodology: dict) -> str:
    prov = _provenance(methodology)
    return str(prov.get("git_commit") or prov.get("git_sha") or "unknown")


def _calibration(methodology: dict) -> dict:
    # Calibration may live at the top level or inside provenance.
    cal = (methodology or {}).get("calibration")
    if not isinstance(cal, dict):
        cal = _provenance(methodology).get("calibration")
    return cal if isinstance(cal, dict) else {}


# ── RO-Crate / W3C-PROV export ───────────────────────────────────────────────

def build_ro_crate(methodology: dict, report_dir: Path | None = None) -> dict:
    """Build an RO-Crate 1.1 metadata graph (JSON-LD) for one report.

    The crate describes the ARIA run as a provenance graph: a CreateAction
    (the workflow run, stamped with version/commit/workflow_hash/seeds), the input
    File entities with their SHA-256, the report.html + methodology.json outputs,
    and one contextual entity per claim that links to its evidence card and its
    run-ledger node. Deterministic; asserts nothing not already in methodology.
    """
    prov = _provenance(methodology)
    version = _version(methodology)
    commit = _commit(methodology)

    graph: list[dict[str, Any]] = [
        {
            "@type": "CreativeWork",
            "@id": "ro-crate-metadata.json",
            "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
            "about": {"@id": "./"},
        },
    ]

    has_part: list[dict] = []
    output_files = ["report.html", _METHODOLOGY_NAME, "ro-crate-metadata.json"]
    for name in output_files:
        has_part.append({"@id": name})

    # Input File entities (with content hashes).
    inputs = methodology.get("inputs", []) or []
    for i, inp in enumerate(inputs):
        if not isinstance(inp, dict):
            continue
        fid = str(inp.get("path") or f"#input-{i}")
        has_part.append({"@id": fid})
        graph.append({
            "@type": "File",
            "@id": fid,
            "name": Path(str(inp.get("path", fid))).name,
            "sha256": inp.get("sha256"),
            "contentSize": inp.get("bytes"),
            "encodingFormat": inp.get("modality"),
            "role": "input",
        })

    # Root dataset.
    graph.append({
        "@type": "Dataset",
        "@id": "./",
        "name": f"ARIA report (v{version}, commit {commit[:12]})",
        "description": (
            "ARIA analysis run evidence graph: provenance, input hashes, "
            "claims linked to evidence cards and run-ledger nodes, and "
            "numerical calibration."
        ),
        "datePublished": prov.get("timestamp_utc"),
        "version": version,
        "hasPart": has_part,
        "mentions": [{"@id": "#aria-run"}],
    })

    # Output File entities.
    for name in output_files:
        graph.append({
            "@type": "File",
            "@id": name,
            "name": name,
            "role": "output",
        })

    # The workflow run as a PROV-flavored CreateAction.
    cal = _calibration(methodology)
    graph.append({
        "@type": "CreateAction",
        "@id": "#aria-run",
        "name": "ARIA analysis run",
        "agent": {"@id": "#aria"},
        "instrument": {"@id": "#aria"},
        "endTime": prov.get("timestamp_utc"),
        "ariaVersion": version,
        "gitCommit": commit,
        "gitDirty": prov.get("git_dirty"),
        "workflowHash": prov.get("workflow_hash"),
        "imageDigest": prov.get("image_digest"),
        "seeds": methodology.get("seeds", {}),
        "calibrationStatus": cal.get("status"),
        "object": [{"@id": str(i.get("path"))}
                   for i in inputs if isinstance(i, dict) and i.get("path")],
        "result": [{"@id": n} for n in output_files],
    })
    graph.append({
        "@type": "SoftwareApplication",
        "@id": "#aria",
        "name": "ARIA",
        "version": version,
        "softwareVersion": version,
    })

    # One contextual entity per claim, linked to its evidence card + ledger node.
    for claim in methodology.get("claims", []) or []:
        if not isinstance(claim, dict):
            continue
        cid = str(claim.get("claim_id") or "")
        if not cid:
            continue
        graph.append({
            "@type": "Claim",
            "@id": f"#claim-{cid}",
            "text": claim.get("text"),
            "claimTier": claim.get("tier"),
            "modality": claim.get("modality"),
            "analysis": claim.get("analysis"),
            "evidenceCard": claim.get("evidence_card_id"),
            "ledgerNode": claim.get("ledger_node_id"),
            "ledgerStatus": claim.get("ledger_status"),
            "isPartOf": {"@id": "./"},
        })

    return {"@context": RO_CRATE_CONTEXT, "@graph": graph}


def write_ro_crate(report_dir: str | Path) -> Path:
    """Write ``ro-crate-metadata.json`` into a report dir from its methodology."""
    rd = Path(report_dir)
    methodology = load_methodology(rd)
    crate = build_ro_crate(methodology, rd)
    out = rd / "ro-crate-metadata.json"
    out.write_text(json.dumps(crate, indent=2, default=str))
    return out


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _file_entry(path: Path, arcname: str, role: str) -> dict[str, Any]:
    return {
        "path": arcname,
        "role": role,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _repo_snapshot_files(repo_root: Path) -> list[tuple[Path, str, str]]:
    files: list[tuple[Path, str, str]] = []
    for pattern in ("*.lock", "*.pip.lock"):
        for path in sorted((repo_root / "envs").glob(pattern)):
            if path.is_file():
                rel = path.relative_to(repo_root).as_posix()
                files.append((path, f"repository/{rel}", "lockfile"))
    for rel in (
        "requirements.lock",
        "docs/architecture/code_graph.md",
        "docs/architecture/graphify/README.md",
        "docs/architecture/graphify/GRAPH_REPORT.md",
        "docs/architecture/graphify/manifest.json",
        "docs/architecture/graphify/graph.json",
    ):
        path = repo_root / rel
        if path.is_file():
            role = "graph_snapshot" if rel.startswith("docs/architecture") else "lockfile"
            files.append((path, f"repository/{rel}", role))
    return files


def _capsule_reproduction_plan(methodology: dict) -> dict[str, Any]:
    prov = _provenance(methodology)
    env = prov.get("environment") if isinstance(prov.get("environment"), dict) else {}
    inputs = methodology.get("inputs", []) or []
    return {
        "status": "verification_ready",
        "automatic_rerun": False,
        "reason": (
            "The capsule can verify report, input, code, graph, and environment "
            "identity. Automatic re-execution requires the original data paths "
            "and biological question/run policy; use the command template after "
            "confirming those inputs."
        ),
        "required_git_commit": _commit(methodology),
        "required_workflow_hash": prov.get("workflow_hash"),
        "required_env_lock_file": env.get("env_lock_file"),
        "required_env_lock_sha256": env.get("env_lock_sha256"),
        "input_hashes": [
            {
                "path": item.get("path"),
                "sha256": item.get("sha256"),
                "size_bytes": item.get("size_bytes") or item.get("bytes"),
                "modality": item.get("modality"),
            }
            for item in inputs if isinstance(item, dict)
        ],
        "command_template": (
            "aria --reproducible  # or use aria.headless.run_headless(...) with "
            "the same input files, biological question, and checkpoint policy"
        ),
    }


def build_capsule_manifest(
    report_dir: str | Path,
    *,
    repo_root: str | Path | None = None,
) -> tuple[dict[str, Any], list[tuple[Path, str, str]]]:
    """Return the S14 capsule manifest and source files to zip."""
    rd = Path(report_dir)
    root = Path(repo_root).resolve() if repo_root else _repo_root()
    methodology = load_methodology(rd)
    write_ro_crate(rd)

    sources: list[tuple[Path, str, str]] = []
    for path in sorted(rd.rglob("*")):
        if path.is_file():
            sources.append((path, path.relative_to(rd.parent).as_posix(), "report_artifact"))
    sources.extend(_repo_snapshot_files(root))

    files = [_file_entry(path, arcname, role) for path, arcname, role in sources]
    manifest = {
        "schema_version": _CAPSULE_SCHEMA_VERSION,
        "report_dir": rd.name,
        "aria_version": _version(methodology),
        "git_commit": _commit(methodology),
        "git_dirty": _provenance(methodology).get("git_dirty"),
        "workflow_hash": _provenance(methodology).get("workflow_hash"),
        "provenance": _provenance(methodology),
        "environment": _provenance(methodology).get("environment", {}),
        "ro_crate": "ro-crate-metadata.json",
        "files": files,
        "reproduction": _capsule_reproduction_plan(methodology),
    }
    return manifest, sources


def write_reproducible_capsule(report_dir: str | Path,
                               out_zip: str | Path | None = None,
                               *,
                               repo_root: str | Path | None = None) -> Path:
    """Bundle a report dir, RO-Crate, lockfiles, graph snapshot, and checksums."""
    rd = Path(report_dir)
    manifest, sources = build_capsule_manifest(rd, repo_root=repo_root)
    out = Path(out_zip) if out_zip else rd.parent / f"{rd.name}_capsule.zip"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for path, arcname, _role in sources:
            if path.resolve() != out.resolve():
                z.write(path, arcname)
        z.writestr(
            _CAPSULE_MANIFEST,
            json.dumps(manifest, indent=2, sort_keys=True, default=str),
        )
    return out


def _zip_sha256(z: zipfile.ZipFile, name: str) -> str:
    import hashlib

    digest = hashlib.sha256()
    with z.open(name, "r") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _methodology_from_zip(z: zipfile.ZipFile, manifest: dict) -> dict:
    report_dir = str(manifest.get("report_dir") or "")
    candidates = [
        f"{report_dir}/{_METHODOLOGY_NAME}" if report_dir else "",
        _METHODOLOGY_NAME,
    ]
    for name in candidates:
        if name and name in z.namelist():
            return json.loads(z.read(name).decode("utf-8"))
    for name in z.namelist():
        if name.endswith(f"/{_METHODOLOGY_NAME}"):
            return json.loads(z.read(name).decode("utf-8"))
    return {}


def verify_reproducible_capsule(
    capsule_zip: str | Path,
    *,
    compare_report: str | Path | None = None,
    repo_root: str | Path | None = None,
) -> dict[str, Any]:
    """Verify a capsule ZIP and optionally diff it against a reproduced report."""
    root = Path(repo_root).resolve() if repo_root else _repo_root()
    result: dict[str, Any] = {
        "capsule": str(capsule_zip),
        "status": "pass",
        "files": {"checked": 0, "missing": [], "mismatched": []},
        "inputs": {"checked": 0, "missing": [], "mismatched": []},
        "environment": {},
        "git": {},
        "diff": None,
    }
    with zipfile.ZipFile(capsule_zip, "r") as z:
        if _CAPSULE_MANIFEST not in z.namelist():
            result["status"] = "fail"
            result["files"]["missing"].append(_CAPSULE_MANIFEST)
            return result
        manifest = json.loads(z.read(_CAPSULE_MANIFEST).decode("utf-8"))
        result["manifest"] = {
            "schema_version": manifest.get("schema_version"),
            "git_commit": manifest.get("git_commit"),
            "workflow_hash": manifest.get("workflow_hash"),
        }
        for entry in manifest.get("files", []) or []:
            if not isinstance(entry, dict):
                continue
            name = str(entry.get("path") or "")
            expected = entry.get("sha256")
            if not name:
                continue
            if name not in z.namelist():
                result["files"]["missing"].append(name)
                continue
            observed = _zip_sha256(z, name)
            result["files"]["checked"] += 1
            if expected and observed != expected:
                result["files"]["mismatched"].append({
                    "path": name, "expected": expected, "observed": observed,
                })

        methodology = _methodology_from_zip(z, manifest)

    for item in (methodology.get("inputs", []) or []):
        if not isinstance(item, dict):
            continue
        path = item.get("path")
        expected = item.get("sha256")
        if not path or not expected or expected == "unavailable":
            continue
        p = Path(str(path)).expanduser()
        if not p.is_file():
            result["inputs"]["missing"].append(str(path))
            continue
        observed = sha256_file(p)
        result["inputs"]["checked"] += 1
        if observed != expected:
            result["inputs"]["mismatched"].append({
                "path": str(path), "expected": expected, "observed": observed,
            })

    try:
        from aria.version import collect_environment_metadata, collect_version_metadata
        current_env = collect_environment_metadata(root)
        expected_env = (_provenance(methodology).get("environment") or {})
        result["environment"] = {
            "expected_env_lock_file": expected_env.get("env_lock_file"),
            "expected_env_lock_sha256": expected_env.get("env_lock_sha256"),
            "current_env_lock_file": current_env.get("env_lock_file"),
            "current_env_lock_sha256": current_env.get("env_lock_sha256"),
            "matches": (
                bool(expected_env.get("env_lock_sha256"))
                and current_env.get("env_lock_sha256") == expected_env.get("env_lock_sha256")
            ),
        }
        current_version = collect_version_metadata(root)
        expected_commit = _commit(methodology)
        result["git"] = {
            "expected_commit": expected_commit,
            "current_commit": current_version.get("git_commit"),
            "matches": current_version.get("git_commit") == expected_commit,
            "current_dirty": current_version.get("git_dirty"),
        }
    except Exception as exc:
        result["environment"] = {"status": "unverified", "error": str(exc)}

    if compare_report:
        result["diff"] = diff_methodologies(methodology, load_methodology(compare_report))

    if (result["files"]["missing"] or result["files"]["mismatched"]
            or result["inputs"]["mismatched"]):
        result["status"] = "fail"
    elif result["inputs"]["missing"]:
        result["status"] = "warning"
    return result


# ── aria diff ────────────────────────────────────────────────────────────────

def _ledger_index(methodology: dict) -> dict[str, dict]:
    entries = (methodology.get("run_ledger", {}) or {}).get("entries", []) or []
    return {e["node_id"]: e for e in entries
            if isinstance(e, dict) and e.get("node_id")}


def _claim_index(methodology: dict) -> dict[str, dict]:
    return {str(c.get("claim_id")): c
            for c in (methodology.get("claims", []) or [])
            if isinstance(c, dict) and c.get("claim_id")}


def diff_methodologies(a: dict, b: dict) -> dict:
    """Structured diff of two report methodologies over their run ledgers.

    Compares provenance (version/commit/dirty/workflow_hash/seeds/input hashes),
    executed analyses (ledger node status), the claim set (added/removed/tier or
    ledger-status changes), and calibration status. ``identical`` is true when no
    tracked field differs.
    """
    prov_a, prov_b = _provenance(a), _provenance(b)
    provenance: dict[str, Any] = {}
    for key in ("version", "git_commit", "git_sha", "git_dirty",
                "workflow_hash", "image_digest"):
        va, vb = prov_a.get(key), prov_b.get(key)
        if va != vb:
            provenance[key] = {"a": va, "b": vb}
    if a.get("seeds") != b.get("seeds"):
        provenance["seeds"] = {"a": a.get("seeds"), "b": b.get("seeds")}

    # Input hashes.
    def _hashes(m):
        return {str(i.get("path")): i.get("sha256")
                for i in (m.get("inputs", []) or []) if isinstance(i, dict)}
    ha, hb = _hashes(a), _hashes(b)
    input_changes = {p: {"a": ha.get(p), "b": hb.get(p)}
                     for p in (set(ha) | set(hb)) if ha.get(p) != hb.get(p)}
    if input_changes:
        provenance["inputs"] = input_changes

    # Ledger nodes.
    la, lb = _ledger_index(a), _ledger_index(b)
    ledger = {
        "only_in_a": sorted(set(la) - set(lb)),
        "only_in_b": sorted(set(lb) - set(la)),
        "status_changed": [
            {"node_id": n, "a": la[n].get("status"), "b": lb[n].get("status")}
            for n in sorted(set(la) & set(lb))
            if la[n].get("status") != lb[n].get("status")
        ],
    }

    # Claims.
    ca, cb = _claim_index(a), _claim_index(b)
    changed = []
    for cid in sorted(set(ca) & set(cb)):
        before, after = ca[cid], cb[cid]
        delta = {}
        for key in ("tier", "ledger_status", "ledger_node_id"):
            if before.get(key) != after.get(key):
                delta[key] = {"a": before.get(key), "b": after.get(key)}
        if delta:
            changed.append({"claim_id": cid, **delta})
    claims = {
        "added": sorted(set(cb) - set(ca)),
        "removed": sorted(set(ca) - set(cb)),
        "changed": changed,
    }

    cal_a, cal_b = _calibration(a), _calibration(b)
    calibration = {}
    if cal_a.get("status") != cal_b.get("status"):
        calibration = {"status": {"a": cal_a.get("status"), "b": cal_b.get("status")}}

    identical = not (provenance or ledger["only_in_a"] or ledger["only_in_b"]
                     or ledger["status_changed"] or claims["added"]
                     or claims["removed"] or claims["changed"] or calibration)
    return {
        "identical": identical,
        "provenance": provenance,
        "ledger": ledger,
        "claims": claims,
        "calibration": calibration,
    }


def format_diff(diff: dict) -> str:
    """Render a structured diff as human-readable text."""
    if diff.get("identical"):
        return "Reports are identical across tracked provenance, ledger, and claims."
    lines = ["Report diff (A vs B):"]
    prov = diff.get("provenance", {})
    if prov:
        lines.append("  Provenance:")
        for key, val in prov.items():
            lines.append(f"    {key}: {val.get('a')} -> {val.get('b')}"
                         if isinstance(val, dict) and "a" in val
                         else f"    {key}: {val}")
    ledger = diff.get("ledger", {})
    if ledger.get("only_in_a"):
        lines.append(f"  Analyses only in A: {', '.join(ledger['only_in_a'])}")
    if ledger.get("only_in_b"):
        lines.append(f"  Analyses only in B: {', '.join(ledger['only_in_b'])}")
    for ch in ledger.get("status_changed", []):
        lines.append(f"  Analysis {ch['node_id']}: {ch['a']} -> {ch['b']}")
    claims = diff.get("claims", {})
    if claims.get("added"):
        lines.append(f"  Claims added in B: {', '.join(claims['added'])}")
    if claims.get("removed"):
        lines.append(f"  Claims removed in B: {', '.join(claims['removed'])}")
    for ch in claims.get("changed", []):
        lines.append(f"  Claim {ch['claim_id']} changed: "
                     f"{ {k: v for k, v in ch.items() if k != 'claim_id'} }")
    cal = diff.get("calibration", {})
    if cal.get("status"):
        lines.append(f"  Calibration status: {cal['status']['a']} -> {cal['status']['b']}")
    return "\n".join(lines)


# ── CLI (routed from aria.tui:main) ──────────────────────────────────────────

def cli_main(argv: list[str]) -> int:
    """`aria diff`, `aria export`, and `aria reproduce` over run ledgers."""
    import argparse

    parser = argparse.ArgumentParser(prog="aria")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_diff = sub.add_parser("diff", help="Diff two report run ledgers.")
    p_diff.add_argument("report_a", help="Report dir or methodology.json")
    p_diff.add_argument("report_b", help="Report dir or methodology.json")
    p_diff.add_argument("--json", action="store_true", help="Emit the structured diff as JSON.")

    p_exp = sub.add_parser("export", help="Write the RO-Crate + reproducible capsule.")
    p_exp.add_argument("report_dir", help="Report directory containing methodology.json")
    p_exp.add_argument("--zip", dest="zip_out", default=None, help="Capsule ZIP output path")
    p_exp.add_argument("--no-zip", action="store_true", help="Only write ro-crate-metadata.json")

    p_rep = sub.add_parser("reproduce", help="Verify a reproducible capsule.")
    p_rep.add_argument("capsule_zip", help="Capsule ZIP created by `aria export`")
    p_rep.add_argument(
        "--compare-report",
        default=None,
        help="Optional reproduced report dir/methodology.json to diff against the capsule.",
    )
    p_rep.add_argument("--json", action="store_true", help="Emit verification as JSON.")

    args = parser.parse_args(argv)

    if args.cmd == "diff":
        d = diff_methodologies(load_methodology(args.report_a),
                               load_methodology(args.report_b))
        print(json.dumps(d, indent=2, default=str) if args.json else format_diff(d))
        return 0

    if args.cmd == "export":
        rd = Path(args.report_dir)
        crate = write_ro_crate(rd)
        print(f"wrote {crate}")
        if not args.no_zip:
            cap = write_reproducible_capsule(rd, args.zip_out)
            print(f"wrote {cap}")
        return 0

    if args.cmd == "reproduce":
        result = verify_reproducible_capsule(
            args.capsule_zip,
            compare_report=args.compare_report,
        )
        if args.json:
            print(json.dumps(result, indent=2, default=str))
        else:
            print(format_reproduce_result(result))
        return 0 if result.get("status") in {"pass", "warning"} else 1

    return 2


def format_reproduce_result(result: dict[str, Any]) -> str:
    """Human-readable capsule verification report."""
    lines = [f"Capsule verification: {result.get('status', 'unknown')}"]
    files = result.get("files", {})
    lines.append(f"  Files checked: {files.get('checked', 0)}")
    if files.get("missing"):
        lines.append(f"  Missing files: {', '.join(files['missing'])}")
    if files.get("mismatched"):
        lines.append(f"  Mismatched files: {len(files['mismatched'])}")

    inputs = result.get("inputs", {})
    lines.append(f"  Inputs checked: {inputs.get('checked', 0)}")
    if inputs.get("missing"):
        lines.append(f"  Inputs missing locally: {', '.join(inputs['missing'])}")
    if inputs.get("mismatched"):
        lines.append(f"  Input hash mismatches: {len(inputs['mismatched'])}")

    env = result.get("environment", {})
    if env:
        lines.append(
            "  Env lock: "
            f"{env.get('current_env_lock_sha256')} "
            f"({'match' if env.get('matches') else 'no match or unrecorded'})"
        )
    git = result.get("git", {})
    if git:
        lines.append(
            "  Git commit: "
            f"{git.get('current_commit')} "
            f"({'match' if git.get('matches') else 'differs'})"
        )
    if result.get("diff"):
        lines.append(format_diff(result["diff"]))
    return "\n".join(lines)
