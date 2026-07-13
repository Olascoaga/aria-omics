#!/usr/bin/env python
"""Reproduce the E6 accession-to-DataAudit public-study validation.

The four accessions are deliberately small processed-data studies spanning the
four ARIA modalities.  GEOConnector still exercises network download, archive
extraction/decompression, content checks, manifest publication, and cache reuse.
Raw SRA fallback is covered deterministically in the E6 guard suite because it
requires an external SRA Toolkit installation and public raw runs can be large.

Usage:
    conda run -n aria-env python scripts/run_e6_geo_sra_validation.py

By default the cache and report live outside the repository under
``~/.aria/workspace/e6_geo_sra_validation``.  Pass ``--output`` when refreshing
the versioned benchmark artifact after reviewing a live run.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

STUDIES = {
    "GSE183948": "bulk_RNA",
    "GSE202494": "scRNA",
    "GSE96769": "scATAC",
    "GSE142660": "bulk_ATAC",
}


class _ValidationMemory:
    """Minimal memory surface required by DataAuditAgent's classification run."""

    @staticmethod
    def create_wing(*_args, **_kwargs):
        return None

    @staticmethod
    def create_hall(*_args, **_kwargs):
        return None


def _git(*args: str) -> str:
    result = subprocess.run(
        ["git", *args], cwd=REPO, text=True, capture_output=True, check=True
    )
    return result.stdout.strip()


def _audit_result(result: dict) -> dict[str, int]:
    from aria.agents.data_audit_agent import DataAuditAgent

    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent.memory = _ValidationMemory()
    agent.publish_status = lambda *_args, **_kwargs: None
    agent.publish_escalation = lambda *_args, **_kwargs: None
    audited = agent.run(
        f"e6_{result['accession'].lower()}",
        {
            "data_dir": result["local_dir"],
            "user_question": "validate public accession ingestion",
            "geo_metadata": result,
        },
    )
    if audited.get("status") != "awaiting_checkpoint":
        raise RuntimeError(
            f"DataAudit rejected {result['accession']}: {audited.get('error')}"
        )
    modalities = audited["exp_context"]["modalities"]
    return {
        modality: len(paths)
        for modality, paths in sorted(modalities.items())
        if paths
    }


def _validate_study(connector, accession: str, expected: str) -> dict:
    from aria.utils.atomic_retrieval import validate_retrieval_manifest

    result = connector.fetch(accession, status_callback=print)
    root = Path(result["local_dir"]).resolve()
    validation = validate_retrieval_manifest(root, expected_accession=accession)
    if validation.get("status") != "valid":
        raise RuntimeError(
            f"manifest rejected {accession}: {validation.get('errors')}"
        )

    manifest = validation["manifest"]
    records = {row["path"]: row for row in manifest["files"]}
    analyzable = []
    for bucket, paths in sorted((result.get("files") or {}).items()):
        for raw_path in paths:
            relative = Path(raw_path).resolve().relative_to(root).as_posix()
            record = records.get(relative)
            if record is None:
                raise RuntimeError(
                    f"declared analyzable payload is not manifested: {relative}"
                )
            analyzable.append({
                "bucket": bucket,
                "path": relative,
                "size": record["size"],
                "sha256": record["sha256"],
            })
    if not analyzable:
        raise RuntimeError(f"{accession} produced no analyzable payload")

    modalities = _audit_result(result)
    if expected not in modalities:
        raise RuntimeError(
            f"{accession}: expected DataAudit modality {expected}, got {modalities}"
        )
    if set(modalities) != {expected}:
        raise RuntimeError(
            f"{accession}: unexpected extra DataAudit modalities: {modalities}"
        )
    return {
        "accession": accession,
        "expected_modality": expected,
        "connector_data_type": result.get("data_type"),
        "retrieval_status": result.get("retrieval_status"),
        "source_accessions": result.get("source_accessions") or [accession],
        "manifest_status": validation["status"],
        "manifest_file_count": len(records),
        "analyzable_payload_count": len(analyzable),
        "analyzable_payloads": analyzable,
        "data_audit_modalities": modalities,
    }


def main() -> int:
    from aria import __version__
    from aria.connectors.geo_connector import GEOConnector

    workspace = Path.home() / ".aria" / "workspace" / "e6_geo_sra_validation"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=workspace / "cache")
    parser.add_argument("--output", type=Path, default=workspace / "report.json")
    args = parser.parse_args()

    connector = GEOConnector(cache_dir=str(args.cache_dir.expanduser()))
    studies = [
        _validate_study(connector, accession, expected)
        for accession, expected in STUDIES.items()
    ]
    report = {
        "artifact": "E6 atomic GEO/SRA accession ingestion validation",
        "aria_version": __version__,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "provenance": {
            "git_commit": _git("rev-parse", "HEAD"),
            "git_dirty": bool(_git("status", "--short")),
            "network": "NCBI GEO FTP over HTTPS",
        },
        "scope": {
            "validates": [
                "four real public GEO studies spanning bulk_RNA, scRNA, scATAC, and bulk_ATAC",
                "atomic retrieval manifests with byte size and SHA-256 verification",
                "archive extraction or gzip decompression where supplied upstream",
                "published connector result accepted as analyzable input by DataAudit",
            ],
            "does_not_validate": [
                "downstream biological analysis of the retrieved payloads",
                "live raw-SRA conversion; external SRA Toolkit behavior is guarded deterministically in tests",
            ],
        },
        "studies": studies,
        "summary": {
            "study_count": len(studies),
            "modalities": sorted({row["expected_modality"] for row in studies}),
            "all_manifests_valid": all(
                row["manifest_status"] == "valid" for row in studies
            ),
            "all_data_audit_transitions_valid": all(
                row["expected_modality"] in row["data_audit_modalities"]
                for row in studies
            ),
        },
    }
    args.output = args.output.expanduser()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
