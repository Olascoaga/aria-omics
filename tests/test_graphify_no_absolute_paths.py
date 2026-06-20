"""F11 (preprint audit 2026-06-19): committed Graphify artifacts must not embed a
machine-absolute filesystem path, and the README snapshot pointer must reference a
real commit reachable from HEAD (not a stale/garbage hash).

The Graphify regen used to rewrite the temp corpus/run_out roots to the ABSOLUTE repo
root, leaking `/home/<user>/Samael/ARIA/...` into manifest/report/html on every regen
(a privacy leak + a reproducibility break — not relocatable). This guard fails on any
absolute path under `docs/architecture/graphify/`, so the regression cannot return.
"""
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRAPHIFY_DIR = ROOT / "docs" / "architecture" / "graphify"

# Machine-absolute path prefixes (not repo-relative paths or URL routes).
_ABSOLUTE_PREFIXES = (
    "/home/", "/Users/", "/root/", "/tmp/", "/mnt/", "/private/", "/data/",
)


def test_graphify_artifacts_have_no_absolute_paths():
    if not GRAPHIFY_DIR.exists():
        return
    offenders: list[str] = []
    for f in sorted(GRAPHIFY_DIR.rglob("*")):
        if not f.is_file() or f.suffix not in {".json", ".md", ".html"}:
            continue
        text = f.read_text(encoding="utf-8", errors="ignore")
        for pref in _ABSOLUTE_PREFIXES:
            idx = text.find(pref)
            if idx != -1:
                offenders.append(
                    f"{f.relative_to(ROOT)}: ...{text[idx:idx + 48]!r}...")
                break
    assert not offenders, (
        "machine-absolute paths in committed Graphify artifacts "
        "(must be repo-relative): " + "; ".join(offenders))


def test_graphify_readme_commit_pointer_is_reachable_from_head():
    readme = GRAPHIFY_DIR / "README.md"
    if not readme.exists():
        return
    m = re.search(r"(?m)^- Commit: `([0-9a-f]{7,40})`",
                  readme.read_text(encoding="utf-8"))
    assert m, "README has no parseable '- Commit: `<sha>`' snapshot pointer"
    sha = m.group(1)
    # The pointer must resolve to a real commit object...
    kind = subprocess.run(
        ["git", "-C", str(ROOT), "cat-file", "-t", sha],
        capture_output=True, text=True)
    assert kind.returncode == 0 and kind.stdout.strip() == "commit", (
        f"README commit pointer {sha} does not resolve to a real commit")
    # ...and be reachable from HEAD (in this branch's history, not a stale/foreign
    # hash). The regen script stamps the build HEAD, so the committed pointer is HEAD
    # or an ancestor (the graph is committed one step after it is built).
    anc = subprocess.run(
        ["git", "-C", str(ROOT), "merge-base", "--is-ancestor", sha, "HEAD"],
        capture_output=True, text=True)
    assert anc.returncode == 0, (
        f"README commit pointer {sha} is not reachable from HEAD (stale/foreign)")
