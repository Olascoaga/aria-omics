"""F6 (preprint audit 2026-06-19): re-stamp provenance onto the V47 bulk ATAC
benchmark artifacts.

The two committed V47 artifacts
(``docs/benchmark_results/bulk_atac/v47_k562_gm12878_da_validation.json`` and
``..._motif_enrichment.json``) were authored by hand during development and shipped
with NO provenance — no ``git_commit``, ``aria_version``, tool versions, seed, or
run date — even though the standard helper ``collect_version_metadata`` exists and
``synthetic_atac_da`` already uses it. That is a reproducibility/credibility gap for
a preprint-citable artifact.

This module is a re-emission *stamper*, not a re-computation runner (Samael's call):
the scientific body of each artifact is the **measured record** from the validated
runs and is preserved verbatim; we only strip any prior provenance block and append a
fresh one (current repo state via ``collect_version_metadata`` + the source commit
where the numbers were actually computed + the tool versions used + determinism note
+ run date). Re-running the full MACS3->bedtools->DESeq2->motif pipeline from the real
ENCODE BAMs is separate preprint-freeze work; it is heavy (multi-env, prior OOM) and
would change the numbers, so it is intentionally out of scope here.

No DE/DA core touched. No fabrication: tool versions/seed are the values documented in
the source runs; the snapatac2 numeric version was not recorded in the source run and
is labelled as such rather than invented.
"""
from __future__ import annotations

import datetime as _dt
import json
from pathlib import Path
from typing import Any

from aria.version import __version__, collect_version_metadata

# ── what we stamp ───────────────────────────────────────────────────────────

# The benchmark directory whose JSONs must all carry provenance (F6 test surface).
BENCH_SUBDIR = Path("docs/benchmark_results/bulk_atac")

# Artifact filename -> the commit where its scientific values were actually computed
# (the "measured-at" commit; distinct from the current HEAD where we re-stamp).
SOURCE_COMMITS: dict[str, str] = {
    "v47_k562_gm12878_da_validation.json": "5c11fc0",
    "v47_k562_gm12878_motif_enrichment.json": "5a407be",
}

# Tool versions actually used in the source runs (documented in the artifact bodies /
# session memory). Not queried from the live env: provenance must record what produced
# the numbers, not whatever happens to be installed now.
TOOL_VERSIONS: dict[str, dict[str, str]] = {
    "v47_k562_gm12878_da_validation.json": {
        "macs3": "3.0.4",
        "bedtools": "2.31.1",
        "pydeseq2": "0.5.4",
    },
    "v47_k562_gm12878_motif_enrichment.json": {
        "engine": "snapatac2.tl.motif_enrichment (shared scATAC chromatin_motifs lane)",
        "snapatac2": "unrecorded_in_source_run",
        "motif_collection": "JASPAR2024_CORE_vertebrates (release 2024, 879 motifs)",
    },
}

# Keys we own and (re)write, so stamping is idempotent regardless of prior state.
_PROVENANCE_KEYS = ("aria_version", "provenance")

_DETERMINISM_NOTE = (
    "Replicate-gated DESeq2 DA and hypergeometric motif enrichment are deterministic; "
    "no RNG seed is used."
)
_STAMP_NOTE = (
    "Scientific values are the measured record from `source_commit`; this provenance "
    "block was re-stamped (not recomputed) at `git_commit`. Full re-computation from "
    "the real ENCODE BAMs is separate preprint-freeze work."
)


def _strip_provenance(body: dict[str, Any]) -> dict[str, Any]:
    """Return ``body`` without any previously stamped provenance keys (idempotency)."""
    return {k: v for k, v in body.items() if k not in _PROVENANCE_KEYS}


def stamp_body(
    body: dict[str, Any],
    artifact_name: str,
    *,
    repo_root: str | Path | None = None,
    now: _dt.datetime | None = None,
) -> dict[str, Any]:
    """Return a copy of ``body`` with a fresh provenance block appended.

    The scientific content is preserved verbatim; only the owned provenance keys are
    replaced. ``artifact_name`` selects the recorded source commit + tool versions.
    """
    stamped = _strip_provenance(body)
    when = (now or _dt.datetime.now(_dt.timezone.utc)).replace(microsecond=0)
    stamped["aria_version"] = __version__
    stamped["provenance"] = {
        **collect_version_metadata(repo_root),
        "stamped_by": "aria/benchmarks/bulk_atac_v47_provenance.py",
        "source_commit": SOURCE_COMMITS.get(artifact_name),
        "tool_versions": TOOL_VERSIONS.get(artifact_name, {}),
        "seed": None,
        "determinism": _DETERMINISM_NOTE,
        "generated_utc": when.isoformat().replace("+00:00", "Z"),
        "note": _STAMP_NOTE,
    }
    return stamped


def stamp_artifacts(
    bench_dir: str | Path | None = None,
    *,
    repo_root: str | Path | None = None,
    now: _dt.datetime | None = None,
) -> dict[str, Any]:
    """Re-stamp provenance on every known V47 artifact found in ``bench_dir``.

    Reads each artifact, preserves its measured body, appends fresh provenance, and
    writes it back (insertion order preserved + trailing newline). Returns a manifest
    of what was stamped. Missing files are reported, never fabricated.
    """
    root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[2]
    directory = Path(bench_dir) if bench_dir is not None else root / BENCH_SUBDIR
    stamped: list[dict[str, Any]] = []
    missing: list[str] = []
    for name in sorted(SOURCE_COMMITS):
        path = directory / name
        if not path.exists():
            missing.append(name)
            continue
        body = json.loads(path.read_text(encoding="utf-8"))
        out = stamp_body(body, name, repo_root=root, now=now)
        path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
        stamped.append({
            "artifact": name,
            "source_commit": out["provenance"]["source_commit"],
            "git_commit": out["provenance"]["git_commit"],
            "aria_version": out["aria_version"],
        })
    return {
        "status": "ok" if stamped and not missing else ("partial" if stamped else "empty"),
        "bench_dir": str(directory),
        "stamped": stamped,
        "missing": missing,
    }
