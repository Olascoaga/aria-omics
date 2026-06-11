"""A1 (audit 2026-06-11): public, preprint-citable benchmark artifacts must not
embed machine-absolute filesystem paths.

`docs/benchmark_results/a1_external/*.json` shipped
`/home/medusa/Samael/ARIA/...` `results_tsv` values — a privacy leak (author home
path) and a reproducibility break (not relocatable). This guard fails on any
absolute system path inside a committed benchmark JSON, so the regression cannot
return. Failing-first: it is RED until the artifacts are relativized.
"""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCH_DIR = ROOT / "docs" / "benchmark_results"

# Prefixes that indicate a machine-absolute filesystem path (not a URL route or
# a repo-relative path). Conservative on purpose to avoid false positives.
_ABSOLUTE_PREFIXES = (
    "/home/", "/Users/", "/root/", "/tmp/", "/mnt/", "/private/", "/data/",
)


def _iter_strings(obj):
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for value in obj.values():
            yield from _iter_strings(value)
    elif isinstance(obj, list):
        for value in obj:
            yield from _iter_strings(value)


def _looks_absolute(s: str) -> bool:
    return s.startswith(_ABSOLUTE_PREFIXES)


def test_public_benchmark_json_has_no_absolute_paths():
    if not BENCH_DIR.exists():
        return
    offenders: list[tuple[str, str]] = []
    for json_file in sorted(BENCH_DIR.rglob("*.json")):
        try:
            data = json.loads(json_file.read_text(encoding="utf-8"))
        except (ValueError, UnicodeDecodeError):
            continue
        for value in _iter_strings(data):
            if _looks_absolute(value):
                offenders.append((str(json_file.relative_to(ROOT)), value))
    assert not offenders, (
        "absolute machine paths in public benchmark artifacts "
        "(use repo-relative paths): " + "; ".join(f"{f}: {v}" for f, v in offenders)
    )


# ── root-cause guard: the A1 generator relativizes paths before serializing ──

def test_repo_relative_normalizes_paths():
    from aria.scripts.benchmark_a1_external_comparators import (
        _relativize_paths, _repo_relative, _REPO_ROOT,
    )
    under_repo = str(_REPO_ROOT / "docs/benchmark_results/a1_external/x.tsv")
    assert _repo_relative(under_repo) == "docs/benchmark_results/a1_external/x.tsv"
    # already relative / non-path / outside-repo
    assert _repo_relative("docs/benchmark_results/x.tsv") == "docs/benchmark_results/x.tsv"
    assert _repo_relative("DESeq2") == "DESeq2"
    assert _repo_relative("/etc/hosts") == "hosts"   # outside repo -> basename
    assert _repo_relative(106) == 106                # non-string passes through
    # recursive over a nested manifest-like object
    obj = {"methods": {"deseq2": {"results_tsv": under_repo, "n": 106}}}
    out = _relativize_paths(obj)
    assert out["methods"]["deseq2"]["results_tsv"] == (
        "docs/benchmark_results/a1_external/x.tsv"
    )
    assert out["methods"]["deseq2"]["n"] == 106
