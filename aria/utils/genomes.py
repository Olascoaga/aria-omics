"""
ARIA reference-genome resolver (v4.6 scATAC motif scanning).
-----------------------------------------------------------
Motif scanning needs a reference genome FASTA. ARIA must resolve it AUTOMATICALLY
from the assembly it already infers (DesignAgent → organism/genome), never by
asking a non-expert user to set an env var. This module is the resolution policy
(light, no snapatac2/heavy deps, importable anywhere):

    explicit path (caller)                       — power user / agent
    → ARIA_GENOME_FASTA                          — power-user env fallback
    → managed store  ARIA_GENOME_DIR/<assembly>/ — curated local genomes + manifest
    → snapatac2 auto-fetch by assembly           — resolved in the script, GOVERNED
      (snap.genome.<attr>, cached; egress only when ARIA_AIR_GAPPED is off)
    → guided checkpoint / honest skip            — only when all of the above fail

The default store is ``~/.aria/genomes`` (mirrors ``ARIA_MOTIF_DIR`` /
``ARIA_GMT_DIR``). The snapatac2 auto-fetch (the part that needs the chromatin
env) is performed by ``aria/scripts/chromatin_motifs.py`` using ``snapatac2_attr``;
this module only names the attribute, it never imports snapatac2.
"""

from __future__ import annotations

import os
from pathlib import Path

from aria.utils.reference_integrity import reference_is_usable, verify_reference_file

GENOME_DIR_ENV = "ARIA_GENOME_DIR"
GENOME_FASTA_ENV = "ARIA_GENOME_FASTA"

_FASTA_GLOBS = ("*.fa", "*.fasta", "*.fa.gz", "*.fasta.gz", "*.fna", "*.fna.gz")

# Assembly name → the set of directory names to look under in the managed store
# (UCSC vs GENCODE/Ensembl naming for the same assembly are interchangeable).
_DIR_ALIASES = {
    "hg38": {"hg38", "grch38"}, "grch38": {"hg38", "grch38"},
    "hg19": {"hg19", "grch37"}, "grch37": {"hg19", "grch37"},
    "mm10": {"mm10", "grcm38"}, "grcm38": {"mm10", "grcm38"},
    "mm39": {"mm39", "grcm39"}, "grcm39": {"mm39", "grcm39"},
}

# Assembly name → the snapatac2 `snap.genome.<attr>` attribute for auto-fetch.
_SNAP_ATTR = {
    "hg38": "hg38", "grch38": "GRCh38",
    "hg19": "hg19", "grch37": "GRCh37",
    "mm10": "mm10", "grcm38": "GRCm38",
    "mm39": "GRCm39", "grcm39": "GRCm39",
}


def _norm(assembly: str | None) -> str:
    return str(assembly or "").strip().lower()


def genome_dir() -> Path:
    """Resolve the managed genome directory (``ARIA_GENOME_DIR`` →
    ``$ARIA_HOME/genomes`` → ``~/.aria/genomes``)."""
    override = os.environ.get(GENOME_DIR_ENV)
    if override:
        return Path(override).expanduser()
    home = os.environ.get("ARIA_HOME")
    base = Path(home).expanduser() if home else (Path.home() / ".aria")
    return base / "genomes"


def genome_fasta_from_env() -> str | None:
    """Return a local genome FASTA path from ``ARIA_GENOME_FASTA`` if set and the
    file exists, else None. Power-user fallback; no network, no fabrication."""
    path = os.environ.get(GENOME_FASTA_ENV)
    if not path:
        return None
    p = Path(path).expanduser()
    return str(p) if p.is_file() else None


def _candidate_dirs(assembly: str) -> set[str]:
    n = _norm(assembly)
    return _DIR_ALIASES.get(n, {n}) if n else set()


def _candidate_fasta_paths(assembly: str | None) -> list[Path]:
    base = genome_dir()
    candidates: list[Path] = []
    for d in _candidate_dirs(assembly):
        sub = base / d
        if sub.is_dir():
            for pat in _FASTA_GLOBS:
                candidates.extend(sorted(sub.glob(pat)))
        for pat in _FASTA_GLOBS:
            candidates.append(base / pat.replace("*", d))
    return candidates


def local_genome_fasta_with_integrity(
    assembly: str | None,
) -> tuple[str | None, dict | None]:
    """Find a usable curated FASTA and return ``(path, integrity)``.

    A declared SHA-256 mismatch makes the candidate unusable. Missing manifests
    remain accepted for backwards compatibility, but are reported as unchecked.
    """
    for cand in _candidate_fasta_paths(assembly):
        if not cand.is_file():
            continue
        integrity = verify_reference_file(cand)
        if reference_is_usable(integrity):
            return str(cand), integrity
        return None, integrity
    return None, None


def local_genome_fasta(assembly: str | None) -> str | None:
    """Find a curated, checksum-usable FASTA for ``assembly`` under the managed store.

    Looks for ``ARIA_GENOME_DIR/<dir>/<*.fa[.gz]>`` and ``ARIA_GENOME_DIR/<dir>.fa``
    across the assembly's alias dir names. No network, no fabrication."""
    path, _integrity = local_genome_fasta_with_integrity(assembly)
    return path


def resolve_local_genome_fasta_with_integrity(
    assembly: str | None,
) -> tuple[str | None, str | None, dict | None]:
    """Resolve a local genome FASTA and include checksum integrity metadata."""
    env = genome_fasta_from_env()
    if env:
        integrity = verify_reference_file(env)
        if reference_is_usable(integrity):
            return env, "env:ARIA_GENOME_FASTA", integrity
        return None, "env:ARIA_GENOME_FASTA", integrity
    local, integrity = local_genome_fasta_with_integrity(assembly)
    if local:
        return local, "managed:ARIA_GENOME_DIR", integrity
    return None, None, integrity


def resolve_local_genome_fasta(assembly: str | None) -> tuple[str | None, str | None]:
    """Resolve a LOCAL genome FASTA (no network): env fallback, then the managed
    store keyed by assembly. Returns ``(path|None, source|None)``."""
    path, source, _integrity = resolve_local_genome_fasta_with_integrity(assembly)
    return path, source


def snapatac2_attr(assembly: str | None) -> str | None:
    """Return the ``snap.genome`` attribute name for ``assembly`` (for governed
    auto-fetch in the chromatin env), or None when the assembly is unknown."""
    return _SNAP_ATTR.get(_norm(assembly))
