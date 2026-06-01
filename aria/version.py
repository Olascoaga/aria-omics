"""Single source of truth for ARIA package version and build stamp."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

__version__ = "4.5.4"

VERSION_SOURCE = "aria.version.__version__"


def _repo_root(repo_root: str | Path | None = None) -> Path:
    if repo_root is not None:
        return Path(repo_root).resolve()
    return Path(__file__).resolve().parents[1]


def _git_text(args: list[str], repo_root: str | Path | None = None) -> str:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=_repo_root(repo_root),
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return "unknown"


def _git_bytes(args: list[str], repo_root: str | Path | None = None) -> bytes:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=_repo_root(repo_root),
            capture_output=True,
            timeout=5,
            check=False,
        )
        if result.returncode == 0:
            return result.stdout
    except Exception:
        pass
    return b""


def _fallback_source_hash(repo_root: Path) -> str:
    digest = hashlib.sha256()
    roots = [
        repo_root / "aria",
        repo_root / "scripts",
        repo_root / "setup.py",
        repo_root / "README.md",
    ]
    for root in roots:
        if not root.exists():
            continue
        paths = [root] if root.is_file() else sorted(root.rglob("*"))
        for path in paths:
            if not path.is_file():
                continue
            rel = path.relative_to(repo_root).as_posix()
            digest.update(rel.encode("utf-8"))
            digest.update(b"\0")
            try:
                digest.update(path.read_bytes())
            except Exception:
                digest.update(b"unreadable")
            digest.update(b"\0")
    return digest.hexdigest()


def workflow_hash(repo_root: str | Path | None = None) -> str:
    """Return a deterministic hash for the code workflow behind this checkout."""
    root = _repo_root(repo_root)
    git_sha = _git_text(["rev-parse", "HEAD"], root)
    git_tree_sha = _git_text(["rev-parse", "HEAD^{tree}"], root)
    if git_sha == "unknown" and git_tree_sha == "unknown":
        return _fallback_source_hash(root)

    status = _git_text(["status", "--porcelain"], root)
    diff_tracked = _git_bytes(["diff", "--binary", "HEAD", "--", "."], root)
    diff_cached = _git_bytes(["diff", "--cached", "--binary", "--", "."], root)
    payload = {
        "version": __version__,
        "git_sha": git_sha,
        "git_tree_sha": git_tree_sha,
        "git_dirty": bool(status and status != "unknown"),
        "git_status_sha256": hashlib.sha256(status.encode("utf-8")).hexdigest(),
        "git_diff_sha256": hashlib.sha256(diff_tracked + diff_cached).hexdigest(),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def version_badge_url(color: str = "blue") -> str:
    """Return the README badge URL derived from the package version."""
    return f"https://img.shields.io/badge/version-{__version__}-{color}"


def collect_version_metadata(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Collect version and source-control metadata for reports and releases."""
    root = _repo_root(repo_root)
    status = _git_text(["status", "--porcelain"], root)
    git_sha = _git_text(["rev-parse", "HEAD"], root)
    return {
        "aria_version": __version__,
        "version": __version__,
        "version_source": VERSION_SOURCE,
        "git_sha": git_sha,
        "git_commit": git_sha,
        "git_dirty": bool(status and status != "unknown"),
        "git_tree_sha": _git_text(["rev-parse", "HEAD^{tree}"], root),
        "git_describe": _git_text(["describe", "--tags", "--dirty", "--always"], root),
        "workflow_hash": workflow_hash(root),
        "workflow_hash_algorithm": (
            "sha256(version+git_sha+git_tree_sha+git_dirty+git_status+git_diff)"
        ),
    }
