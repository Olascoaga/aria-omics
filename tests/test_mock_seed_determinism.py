"""A2 (audit 2026-06-11): the `_mock_*` helpers seeded with `hash(name)`, which is
randomized per process (PYTHONHASHSEED), so mock outputs were not reproducible
across runs. They live only behind `allow_mock` (ADR-002), so this is cosmetic —
but "seed pinned everywhere it matters" should hold. The fix uses a deterministic
hash (zlib.crc32).

Failing-first: under the old `hash(name)` seed, two subprocesses with different
PYTHONHASHSEED produce different mock outputs.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

_SNIPPETS = {
    "align": (
        "from aria.scripts.rna_align import _mock_alignment;"
        "import json;print(json.dumps(_mock_alignment('sampleX','/tmp/x_'),"
        "sort_keys=True))"
    ),
    "fastp": (
        "from aria.scripts.rna_fastq_qc import _mock_fastp_result;"
        "import json;print(json.dumps(_mock_fastp_result('sampleX','r1','r2',"
        "'o1','o2',True),sort_keys=True))"
    ),
}


def _run(snippet: str, hashseed: str) -> str:
    env = {**os.environ, "PYTHONHASHSEED": hashseed, "PYTHONPATH": str(ROOT)}
    proc = subprocess.run(
        [sys.executable, "-c", snippet],
        capture_output=True, text=True, env=env, cwd=str(ROOT),
    )
    assert proc.returncode == 0, proc.stderr
    return proc.stdout.strip()


def test_mock_alignment_deterministic_across_hashseeds():
    a = _run(_SNIPPETS["align"], "0")
    b = _run(_SNIPPETS["align"], "12345")
    assert a and json.loads(a) == json.loads(b)


def test_mock_fastp_deterministic_across_hashseeds():
    a = _run(_SNIPPETS["fastp"], "0")
    b = _run(_SNIPPETS["fastp"], "12345")
    assert a and json.loads(a) == json.loads(b)
