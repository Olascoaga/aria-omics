"""
ARIA scRNA FASTQ quantification via kallisto-bustools.

Executed inside the RNA stack by EnvironmentManager. This script is a thin IPC
wrapper around the deterministic raw-ingestion helper so FASTQ quantification
runs under ARIA's normal subprocess, timeout, failed-run archive, and contract
controls.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from aria.scripts._base import run_script
from aria.utils.raw_ingestion import execute_kb_count


def rna_kb_count(params: dict) -> dict:
    result = execute_kb_count(params)
    if result.get("status") == "blocked":
        blockers = [str(b) for b in result.get("blockers", [])]
        return {
            **result,
            "status": "error",
            "error_type": "KbInputBlocked",
            "details": "; ".join(blockers) or "kb count inputs are incomplete.",
        }
    if result.get("status") == "failed":
        return {
            **result,
            "status": "error",
            "error_type": "KbCountFailed",
            "details": result.get("details") or (
                f"kb count exited with return code {result.get('returncode')}."
            ),
        }
    return result


if __name__ == "__main__":
    run_script(rna_kb_count)
