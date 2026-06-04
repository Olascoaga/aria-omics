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

RO_CRATE_CONTEXT = "https://w3id.org/ro/crate/1.1/context"
_METHODOLOGY_NAME = "methodology.json"


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


def write_reproducible_capsule(report_dir: str | Path,
                               out_zip: str | Path | None = None) -> Path:
    """Bundle a report dir (+ its RO-Crate metadata) into a reproducible ZIP."""
    rd = Path(report_dir)
    # Ensure the crate exists/refreshes before bundling.
    if (rd / _METHODOLOGY_NAME).exists():
        write_ro_crate(rd)
    out = Path(out_zip) if out_zip else rd.parent / f"{rd.name}_capsule.zip"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for path in sorted(rd.rglob("*")):
            if path.is_file() and path.resolve() != out.resolve():
                z.write(path, path.relative_to(rd.parent))
    return out


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
    """`aria diff A B` and `aria export <reportDir> [--zip out.zip]`."""
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

    return 2
