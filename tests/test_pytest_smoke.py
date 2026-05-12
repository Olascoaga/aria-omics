from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_package_compiles():
    result = subprocess.run(
        [sys.executable, "-m", "compileall", "-q", "aria"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr or result.stdout


def test_bulk_rna_legacy_script_passes():
    env = os.environ.copy()
    env["ARIA_ALLOW_MOCKS"] = "1"
    result = subprocess.run(
        [sys.executable, "tests/test_bulk_rna.py"],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr
