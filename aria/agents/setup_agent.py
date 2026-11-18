"""
ARIA SetupAgent
---------------
Guarantees the computational environment is ready before any analysis.
The user never touches conda, pip, STAR, or GTF files.

Philosophy: ARIA owns its environments. No detection of what the user
already has. No aliases. No clever matching.

ARIA ships with YAMLs for every stack in envs/.
SetupAgent installs them if missing. That's it.

Lifecycle:
  First run  → installs envs + downloads genome + builds index (~1h)
  Every run  → checks if ready, skips if yes (~2 seconds)
  User sees  → a progress bar, nothing else

What SetupAgent does NOT do:
  - Search for existing tools on the user's system
  - Try to reuse user's conda environments
  - Ask the user anything about their setup
"""

from __future__ import annotations

import json
import logging
import os
import platform
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence
from aria.llm.provider import LLMProvider
from aria.memory.memory import ARIAMemory
from aria.utils.environment_specs import (
    environment_names_by_stack,
    setup_environment_definitions,
)

log = logging.getLogger("aria.setup")


# ── Reference genome registry ─────────────────────────────────────────────────

GENOME_REGISTRY = {
    "hg38": {
        "organism":       "Homo sapiens",
        "aliases":        ["human", "homo sapiens", "hg38", "grch38"],
        "fasta_url":      (
            "https://ftp.ensembl.org/pub/release-112/fasta/"
            "homo_sapiens/dna/"
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        ),
        "gtf_url": (
            "https://ftp.ensembl.org/pub/release-112/gtf/"
            "homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
        ),
        "star_sa_bases":  14,
        "size_gb":        3.1,
    },
    "mm39": {
        "organism":       "Mus musculus",
        "aliases":        ["mouse", "mus musculus", "mm39", "grcm39", "murine"],
        "fasta_url":      (
            "https://ftp.ensembl.org/pub/release-112/fasta/"
            "mus_musculus/dna/"
            "Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
        ),
        "gtf_url": (
            "https://ftp.ensembl.org/pub/release-112/gtf/"
            "mus_musculus/Mus_musculus.GRCm39.112.gtf.gz"
        ),
        "star_sa_bases":  14,
        "size_gb":        2.7,
    },
    "mm10": {
        "organism":       "Mus musculus",
        "aliases":        ["mm10", "grcm38"],
        "fasta_url":      (
            "https://ftp.ensembl.org/pub/release-102/fasta/"
            "mus_musculus/dna/"
            "Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
        ),
        "gtf_url": (
            "https://ftp.ensembl.org/pub/release-102/gtf/"
            "mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
        ),
        "star_sa_bases":  14,
        "size_gb":        2.7,
    },
    "rn7": {
        "organism":       "Rattus norvegicus",
        "aliases":        ["rat", "rattus norvegicus", "rn7"],
        "fasta_url":      (
            "https://ftp.ensembl.org/pub/release-112/fasta/"
            "rattus_norvegicus/dna/"
            "Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.fa.gz"
        ),
        "gtf_url": (
            "https://ftp.ensembl.org/pub/release-112/gtf/"
            "rattus_norvegicus/"
            "Rattus_norvegicus.mRatBN7.2.112.gtf.gz"
        ),
        "star_sa_bases":  14,
        "size_gb":        2.9,
    },
    "danRer11": {
        "organism":       "Danio rerio",
        "aliases":        ["zebrafish", "danio rerio", "danrer11"],
        "fasta_url":      (
            "https://ftp.ensembl.org/pub/release-112/fasta/"
            "danio_rerio/dna/"
            "Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
        ),
        "gtf_url": (
            "https://ftp.ensembl.org/pub/release-112/gtf/"
            "danio_rerio/Danio_rerio.GRCz11.112.gtf.gz"
        ),
        "star_sa_bases":  14,
        "size_gb":        1.4,
    },
}

# ── Conda environment definitions ─────────────────────────────────────────────
# Derived from the single registry; dependency content lives only in envs/*.yml.
ARIA_ENVS = setup_environment_definitions()


class SetupAgent(BaseAgent):
    """
    Provisions everything ARIA needs before analysis starts.
    Runs as Checkpoint 0, before DataAuditAgent.
    Silent on subsequent runs when everything is already in place.
    """

    name        = "setup_agent"
    description = "Installs environments and downloads references on first run."

    ARIA_HOME = Path.home() / ".aria"

    def __init__(self, memory: ARIAMemory,
                 llm: LLMProvider,
                 api_key: str = None):
        super().__init__(memory, llm, api_key)
        self.genomes_dir = self.ARIA_HOME / "genomes"
        self.ARIA_HOME.mkdir(exist_ok=True)
        self.genomes_dir.mkdir(exist_ok=True)

    # ── Main entry point ──────────────────────────────────────────────────

    def run(self, experiment_id: str, context: dict) -> dict:
        exp_ctx    = context.get("exp_context", {})
        intent     = context.get("biological_intent", {})
        organism   = exp_ctx.get("organism", "")
        has_fastq  = self._has_fastq(exp_ctx)
        modalities = exp_ctx.get("modalities", {})

        self.publish_status(experiment_id,
                            "Checking environment...", 0.0)

        actions  = []
        warnings = []

        # ── 1. Conda ──────────────────────────────────────────────────────
        conda_cmd = self._ensure_conda(experiment_id, actions, warnings)
        if conda_cmd is None:
            return self._fatal(experiment_id, warnings,
                "Could not find or install conda. "
                "Please install Miniforge from: "
                "https://github.com/conda-forge/miniforge/releases/latest"
            )

        # ── 2. Conda environments ─────────────────────────────────────────
        needed_envs = self._needed_envs(modalities, has_fastq, exp_ctx=exp_ctx)
        for env_name in needed_envs:
            self._ensure_env(
                env_name, conda_cmd,
                experiment_id, actions, warnings
            )

        # ── 3. Reference genome + GTF ─────────────────────────────────────
        audited_genome = exp_ctx.get("genome")
        genome_key = self._detect_genome(
            organism, intent, audited_genome=audited_genome
        )
        genome_cfg = {}

        if genome_key:
            self.publish_status(
                experiment_id,
                f"Checking reference genome ({genome_key})...",
                0.6,
            )
            genome_cfg = self._ensure_genome(
                genome_key, experiment_id, actions, warnings
            )

            # ── 4. STAR index (bulk RNA FASTQ lanes only) ────────────────
            needs_star = has_fastq and any(
                m in modalities for m in ("bulk_RNA", "bulk_RNA_raw")
            )
            if needs_star and genome_cfg.get("fasta_ready"):
                self.publish_status(
                    experiment_id,
                    f"Checking STAR index ({genome_key})...",
                    0.8,
                )
                self._ensure_star_index(
                    genome_key, genome_cfg,
                    conda_cmd, experiment_id, actions, warnings
                )

            # ── 5. kb-python index + t2g (scRNA / Multiome FASTQ lanes) ──
            needs_kb = has_fastq and (
                "scRNA" in modalities
                or bool(exp_ctx.get("multiome_fastq_manifest"))
                or (
                    isinstance(exp_ctx.get("paired_multiome_contract"), dict)
                    and exp_ctx["paired_multiome_contract"].get("source") == "paired_fastq_manifest"
                )
            )
            if needs_kb and genome_cfg.get("fasta_ready") and genome_cfg.get("gtf_ready"):
                self.publish_status(
                    experiment_id,
                    f"Checking kb-python transcriptome index ({genome_key})...",
                    0.88,
                )
                ensure_kb = getattr(self, "_ensure_kb_index", None)
                if callable(ensure_kb):
                    ensure_kb(
                        genome_key, genome_cfg,
                        conda_cmd, experiment_id, actions, warnings
                    )
        elif str(audited_genome or "").strip().lower() not in {"", "unknown", "none"}:
            genome_cfg = self._local_staged_genome_cfg(
                str(audited_genome).strip(), warnings
            )

        self.publish_status(experiment_id, "Environment ready.", 1.0)

        summary = (
            f"Setup complete. {len(actions)} action(s): "
            f"{', '.join(actions[:3])}{'...' if len(actions) > 3 else ''}."
            if actions else "Environment ready — all tools present."
        )

        self.publish_finding(
            experiment_id,
            {"summary": summary, "actions": actions,
             "genome_config": genome_cfg},
            Confidence.HIGH,
        )

        return {
            "status":        "done",
            "genome_config": genome_cfg,
            "actions":       actions,
            "warnings":      warnings,
        }

    # ── Conda ─────────────────────────────────────────────────────────────

    def _ensure_conda(self, experiment_id: str,
                       actions: list, warnings: list) -> Optional[str]:
        """Return path to conda/mamba, installing Miniforge if needed."""
        # Check PATH and common locations
        for cmd in ["mamba", "conda"]:
            if shutil.which(cmd):
                return cmd

        for p in [
            Path.home() / "miniforge3"  / "bin" / "conda",
            Path.home() / "anaconda3"   / "bin" / "conda",
            Path.home() / "miniconda3"  / "bin" / "conda",
            Path("/opt/conda/bin/conda"),
        ]:
            if p.exists():
                return str(p)

        # Not found — install Miniforge
        self.publish_status(
            experiment_id,
            "Installing Miniforge (first-time setup, ~2 min)...",
            0.05,
        )
        result = self._install_miniforge()
        if result == "ok":
            actions.append("installed_miniforge")
            conda = str(Path.home() / "miniforge3" / "bin" / "conda")
            return conda if Path(conda).exists() else None

        warnings.append(f"Miniforge install failed: {result}")
        return None

    def _install_miniforge(self) -> str:
        """Download and install Miniforge silently."""
        # Preprint audit A1: refuse the Miniforge download/install under air-gap;
        # an air-gapped run expects a pre-provisioned conda, not a network install.
        from aria.utils.privacy import egress_allowed
        if not egress_allowed():
            return ("ARIA_AIR_GAPPED is enabled; refusing to download/install "
                    "Miniforge. Provide a pre-installed conda in air-gapped mode.")
        system  = platform.system().lower()
        machine = platform.machine().lower()

        url_map = {
            ("linux",  "x86_64"):  "Miniforge3-Linux-x86_64.sh",
            ("linux",  "aarch64"): "Miniforge3-Linux-aarch64.sh",
            ("darwin", "arm64"):   "Miniforge3-MacOSX-arm64.sh",
            ("darwin", "x86_64"):  "Miniforge3-MacOSX-x86_64.sh",
        }
        arch    = "aarch64" if "aarch" in machine else \
                  "arm64"   if "arm"   in machine else "x86_64"
        fname   = url_map.get((system, arch))
        if not fname:
            return f"Unsupported platform: {system}/{machine}"

        url     = f"https://github.com/conda-forge/miniforge/releases/latest/download/{fname}"
        script  = Path("/tmp/miniforge_install.sh")
        install = Path.home() / "miniforge3"

        try:
            subprocess.run(
                ["curl", "-fsSL", "-o", str(script), url],
                check=True, timeout=300,
            )
            subprocess.run(
                ["bash", str(script), "-b", "-p", str(install)],
                check=True, timeout=600,
            )
            # Add to PATH for this process
            new_path = str(install / "bin")
            os.environ["PATH"] = new_path + ":" + os.environ.get("PATH", "")
            return "ok"
        except subprocess.CalledProcessError as e:
            return str(e)
        except subprocess.TimeoutExpired:
            return "timeout"
        except FileNotFoundError:
            return "curl not found"

    def _list_envs(self, conda_cmd: str) -> set[str]:
        """List installed conda environment names."""
        try:
            r = subprocess.run(
                [conda_cmd, "env", "list", "--json"],
                capture_output=True, text=True, timeout=30,
            )
            data = json.loads(r.stdout)
            return {Path(e).name for e in data.get("envs", [])}
        except Exception:
            return set()

    def _needed_envs(
        self,
        modalities: dict,
        has_fastq: bool,
        exp_ctx: Optional[dict] = None,
    ) -> list[str]:
        """Which ARIA environments does this experiment need?"""
        by_stack = environment_names_by_stack()
        envs = []
        modality_set = set(modalities)
        ctx = exp_ctx or {}
        has_multiome_manifest = bool(ctx.get("multiome_fastq_manifest"))
        if has_fastq and ("scRNA" in modality_set or has_multiome_manifest):
            envs.append(by_stack["ingestion"])
        if has_fastq and ({"bulk_RNA", "bulk_RNA_raw"} & modality_set):
            envs.append(by_stack["rnaseq"])
        if {"scRNA", "bulk_RNA"} & modality_set or has_multiome_manifest:
            envs.append(by_stack["rna"])
        chromatin = {
            "scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"
        } & modality_set
        if chromatin or has_multiome_manifest:
            envs.extend([
                by_stack["chromatin"],
                by_stack["rna"],
                by_stack["tobias"],
            ])
            if has_fastq or has_multiome_manifest:
                envs.extend([by_stack["rnaseq"], by_stack["atacseq"]])
        if {"HiC", "Micro-C"} & modality_set:
            envs.append(by_stack["hic"])
            if has_fastq:
                envs.extend([by_stack["rnaseq"], by_stack["atacseq"]])
        has_paired_h5mu = any(
            str(f).lower().endswith(".h5mu")
            for files in modalities.values()
            for f in (files or [])
        )
        has_paired_detection = any(
            isinstance(rec, dict)
            and str((rec.get("evidence") or {}).get("format") or "")
            in {"h5mu", "10x_multiome_h5", "10x_multiome_mtx"}
            for rec in (ctx.get("assay_detections") or [])
        )
        if (
            has_paired_h5mu
            or has_paired_detection
            or has_multiome_manifest
            or ({"scRNA", "scATAC"} <= modality_set)
        ):
            envs.append(by_stack["chromatin"])
            envs.append(by_stack["integration"])
        return list(dict.fromkeys(envs))

    def _ensure_env(self, env_name: str, conda_cmd: str,
                     experiment_id: str,
                     actions: list, warnings: list):
        """Create aria conda env from bundled YAML if not present."""
        installed = self._list_envs(conda_cmd)
        if env_name in installed:
            log.debug(f"Env '{env_name}' already installed.")
            return

        env_def = ARIA_ENVS.get(env_name)
        if not env_def:
            warnings.append(f"No definition found for env '{env_name}'")
            return

        self.publish_status(
            experiment_id,
            f"Installing {env_name} ({env_def['description']})...",
            0.2,
        )
        log.info(f"Creating conda env: {env_name}")

        yml_path = Path(env_def["yml_path"])
        if not yml_path.exists():
            warnings.append(
                f"Versioned environment spec missing for {env_name}: {yml_path}"
            )
            return

        try:
            r = subprocess.run(
                [conda_cmd, "env", "create",
                 "-f", str(yml_path), "--name", env_name, "-y"],
                capture_output=True, text=True,
                timeout=1800,   # 30 min
            )
            if r.returncode == 0:
                actions.append(f"created:{env_name}")
                log.info(f"Created env: {env_name}")
            else:
                warnings.append(
                    f"Failed to create {env_name}: {r.stderr[-200:]}"
                )
        except subprocess.TimeoutExpired:
            warnings.append(f"Timeout creating {env_name} (>30 min)")
        except Exception as e:
            warnings.append(f"Error creating {env_name}: {e}")

    # ── Reference genome ─────────────────────────────────────────────────

    _CANONICAL_ASSEMBLY_ALIASES = {
        "hg38": "hg38",
        "grch38": "hg38",
        "hg19": "hg19",
        "grch37": "hg19",
        "mm39": "mm39",
        "grcm39": "mm39",
        "mm10": "mm10",
        "grcm38": "mm10",
        "rn7": "rn7",
        "mratbn7.2": "rn7",
        "danrer11": "danRer11",
        "grcz11": "danRer11",
        "dm6": "dm6",
        "ce11": "ce11",
        "saccer3": "sacCer3",
    }

    def _detect_genome(self, organism: str,
                        intent: dict,
                        audited_genome: Optional[str] = None) -> Optional[str]:
        """Resolve genome key honoring DataAuditAgent's audited assembly first.

        Never substitutes a different assembly of the same species when an
        explicit ``audited_genome`` (e.g. ``mm10`` or ``hg19``) was declared.
        """
        raw_audited = str(audited_genome or "").strip()
        if raw_audited.lower() not in {"", "unknown", "none"}:
            canonical = self._CANONICAL_ASSEMBLY_ALIASES.get(
                raw_audited.lower(), raw_audited
            )
            return canonical if canonical in GENOME_REGISTRY else None

        text = " ".join([
            str(organism or ""),
            str((intent or {}).get("summary", "")),
            str((intent or {}).get("user_question", "")),
        ]).lower()

        # Check explicit assembly tokens BEFORE generic species names so
        # e.g. "mouse mm10" resolves to mm10 rather than generic mouse -> mm39.
        for token, canonical in self._CANONICAL_ASSEMBLY_ALIASES.items():
            if token in text:
                return canonical if canonical in GENOME_REGISTRY else None

        for key, info in GENOME_REGISTRY.items():
            if any(alias in text for alias in info["aliases"]):
                return key

        # Fallback heuristics
        if any(k in text for k in ["human", "homo", "patient", "clinical"]):
            return "hg38"
        if any(k in text for k in ["mouse", "murine", "mus musculus"]):
            return "mm39"
        if any(k in text for k in ["rat", "rattus"]):
            return "rn7"

        return None

    def _local_staged_genome_cfg(self, audited_genome: str,
                                  warnings: list) -> dict:
        """Inspect ~/.aria/genomes/<key> for locally staged non-registry assemblies."""
        canonical = self._CANONICAL_ASSEMBLY_ALIASES.get(
            audited_genome.lower(), audited_genome
        )
        gdir = self.genomes_dir / canonical
        fasta = gdir / "genome.fa"
        fasta_gz = gdir / "genome.fa.gz"
        gtf = gdir / "annotation.gtf"
        gtf_gz = gdir / "annotation.gtf.gz"
        star_idx = gdir / "star_index"
        fasta_ready = fasta.exists() and fasta.stat().st_size > 0
        gtf_ready = gtf.exists() and gtf.stat().st_size > 0
        index_ready = (star_idx / "SA").exists()
        if not fasta_ready:
            warnings.append(
                f"Audited assembly '{canonical}' is not auto-provisioned by "
                f"SetupAgent and was not found under {gdir}; refusing to "
                f"substitute another assembly."
            )
        return {
            "genome_key":  canonical,
            "genome_dir":  str(gdir),
            "fasta":       str(fasta),
            "fasta_gz":    str(fasta_gz),
            "gtf":         str(gtf),
            "gtf_gz":      str(gtf_gz),
            "star_index":  str(star_idx),
            "fasta_ready": fasta_ready,
            "gtf_ready":   gtf_ready,
            "index_ready": index_ready,
            "strand":      "auto",
        }

    def _ensure_genome(self, genome_key: str,
                        experiment_id: str,
                        actions: list,
                        warnings: list) -> dict:
        """
        Download FASTA + GTF from Ensembl if not in ~/.aria/genomes/.
        Decompresses FASTA (STAR cannot read .gz FASTAs).
        Returns genome_cfg dict for downstream agents.
        """
        info      = GENOME_REGISTRY[genome_key]
        gdir      = self.genomes_dir / genome_key
        gdir.mkdir(exist_ok=True)

        fasta_gz  = gdir / "genome.fa.gz"
        fasta     = gdir / "genome.fa"        # decompressed (STAR needs this)
        gtf_gz    = gdir / "annotation.gtf.gz"

        cfg = {
            "genome_key":  genome_key,
            "genome_dir":  str(gdir),
            "fasta":       str(fasta),           # decompressed
            "fasta_gz":    str(fasta_gz),
            "gtf":         str(gdir / "annotation.gtf"),   # decompressed
            "gtf_gz":      str(gtf_gz),
            "star_index":  str(gdir / "star_index"),
            "kb_dir":      str(gdir / "kb_index"),
            "kb_index":    str(gdir / "kb_index" / "index.idx"),
            "kb_t2g":      str(gdir / "kb_index" / "t2g.txt"),
            "fasta_ready": False,
            "gtf_ready":   False,
            "index_ready": False,
            "kb_ready":    False,
            "strand":      "auto",               # trigger strand auto-detection
        }

        # ── FASTA: download if missing, then decompress ──────────────────
        if not (fasta_gz.exists() and fasta_gz.stat().st_size > 1_000_000):
            self.publish_status(
                experiment_id,
                f"Downloading {genome_key} genome "
                f"(~{info['size_gb']:.0f} GB, one-time)...",
                0.62,
            )
            err = self._download(info["fasta_url"], fasta_gz)
            if err:
                warnings.append(f"FASTA download failed: {err}")
                return cfg
            actions.append(f"downloaded_fasta:{genome_key}")

        # Decompress if not already done
        if not (fasta.exists() and fasta.stat().st_size > 1_000_000):
            self.publish_status(
                experiment_id,
                f"Decompressing {genome_key} FASTA (required by STAR)...",
                0.66,
            )
            err = self._decompress(fasta_gz, fasta)
            if err:
                warnings.append(f"FASTA decompression failed: {err}")
                return cfg
            actions.append(f"decompressed_fasta:{genome_key}")

        cfg["fasta_ready"] = True

        # ── GTF: download, decompress (STAR + featureCounts both prefer plain) ──
        gtf     = gdir / "annotation.gtf"

        if not (gtf_gz.exists() and gtf_gz.stat().st_size > 100_000):
            self.publish_status(
                experiment_id,
                f"Downloading {genome_key} annotation (GTF)...",
                0.72,
            )
            err = self._download(info["gtf_url"], gtf_gz)
            if err:
                warnings.append(f"GTF download failed: {err}")
                return cfg
            actions.append(f"downloaded_gtf:{genome_key}")

        if not (gtf.exists() and gtf.stat().st_size > 100_000):
            self.publish_status(
                experiment_id,
                f"Decompressing {genome_key} GTF...",
                0.74,
            )
            err = self._decompress(gtf_gz, gtf)
            if err:
                warnings.append(f"GTF decompression failed: {err}")
            else:
                actions.append(f"decompressed_gtf:{genome_key}")

        if gtf.exists() and gtf.stat().st_size > 100_000:
            cfg["gtf"] = str(gtf)          # use decompressed
            cfg["gtf_gz"] = str(gtf_gz)
            cfg["gtf_ready"] = True

        return cfg

    def _decompress(self, src_gz: Path, dest: Path) -> Optional[str]:
        """Decompress a .gz file. Uses gunzip via subprocess (streams, low memory)."""
        try:
            # gunzip keeps the source; we use -c and redirect for safety.
            # For a ~900 MB hg38 FASTA this takes ~30 s on a typical machine.
            with open(dest, "wb") as out:
                r = subprocess.run(
                    ["gunzip", "-c", str(src_gz)],
                    stdout=out, stderr=subprocess.PIPE, timeout=600,
                )
            if r.returncode != 0:
                dest.unlink(missing_ok=True)
                return r.stderr.decode("utf-8", errors="ignore")[:200]
            if not dest.exists() or dest.stat().st_size == 0:
                return "decompression produced empty file"
            return None
        except subprocess.TimeoutExpired:
            dest.unlink(missing_ok=True)
            return "gunzip timed out (>10 min)"
        except FileNotFoundError:
            # Fallback: pure Python gzip
            return self._decompress_python(src_gz, dest)
        except Exception as e:
            dest.unlink(missing_ok=True)
            return str(e)

    def _decompress_python(self, src_gz: Path, dest: Path) -> Optional[str]:
        """Pure-Python decompression fallback (no gunzip binary)."""
        import gzip as _gzip
        try:
            with _gzip.open(src_gz, "rb") as f_in, \
                 open(dest, "wb") as f_out:
                # Stream in 16 MB chunks
                while True:
                    chunk = f_in.read(16 * 1024 * 1024)
                    if not chunk:
                        break
                    f_out.write(chunk)
            return None
        except Exception as e:
            dest.unlink(missing_ok=True)
            return str(e)

    def _download(self, url: str, dest: Path) -> Optional[str]:
        """Download a file. Returns None on success, error string on failure."""
        # Preprint audit A1: setup egress must honor the air-gap gate too, not just
        # the LLM layer. Fail closed — refuse the download and surface an honest
        # reason instead of silently shelling out to curl/wget.
        from aria.utils.privacy import egress_allowed
        if not egress_allowed():
            return (f"ARIA_AIR_GAPPED is enabled; refusing to download {url}. "
                    f"Stage the reference locally or disable air-gapped mode.")
        for downloader in [
            ["curl", "-fsSL", "--progress-bar", "-o", str(dest), url],
            ["wget", "-q", "--show-progress", "-O", str(dest), url],
        ]:
            if not shutil.which(downloader[0]):
                continue
            try:
                r = subprocess.run(downloader, timeout=7200)
                if r.returncode == 0 and dest.exists() \
                        and dest.stat().st_size > 0:
                    return None
            except subprocess.TimeoutExpired:
                return "timeout (>2h)"
            except Exception as e:
                continue

        return "curl and wget both failed or not found"

    # ── Reference Index Build Lock + STAR / kb-python Indexes ─────────────

    @staticmethod
    def _reference_build_lock(lock_path: Path):
        """Serialize concurrent reference-index builds sharing a genome directory."""
        from contextlib import contextmanager

        @contextmanager
        def _ctx():
            fh = None
            try:
                lock_path.parent.mkdir(parents=True, exist_ok=True)
                fh = open(lock_path, "a+", encoding="utf-8")
                try:
                    import fcntl
                    fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
                except (ImportError, OSError):
                    pass
                yield
            finally:
                if fh is not None:
                    try:
                        import fcntl
                        fcntl.flock(fh.fileno(), fcntl.LOCK_UN)
                    except (ImportError, OSError):
                        pass
                    fh.close()

        return _ctx()

    def _ensure_star_index(self, genome_key: str,
                            genome_cfg: dict,
                            conda_cmd: str,
                            experiment_id: str,
                            actions: list,
                            warnings: list):
        """Build STAR index if not present. ~30 min, one time only."""
        index_dir = Path(genome_cfg["star_index"])
        ready_files = ["SA", "SAindex", "Genome"]
        lock_file = Path(genome_cfg.get("genome_dir", index_dir.parent)) / ".star_build.lock"

        with self._reference_build_lock(lock_file):
            if (index_dir.exists()
                    and all((index_dir / f).exists() for f in ready_files)):
                genome_cfg["index_ready"] = True
                return

            if not genome_cfg.get("fasta_ready") or not genome_cfg.get("gtf_ready"):
                warnings.append(
                    "STAR index cannot be built: genome files not ready"
                )
                return

            self.publish_status(
                experiment_id,
                f"Building STAR index for {genome_key} "
                f"(~30 min, one-time setup)...",
                0.83,
            )
            index_dir.mkdir(parents=True, exist_ok=True)

            n_threads = str(min(os.cpu_count() or 4, 16))
            sa_bases  = str(GENOME_REGISTRY[genome_key]["star_sa_bases"])

            cmd = [
                conda_cmd, "run", "-n", "aria-rnaseq-env",
                "STAR",
                "--runMode",             "genomeGenerate",
                "--genomeDir",           str(index_dir),
                "--genomeFastaFiles",    genome_cfg["fasta"],
                "--sjdbGTFfile",         genome_cfg["gtf"],
                "--runThreadN",          n_threads,
                "--genomeSAindexNbases", sa_bases,
            ]

            try:
                r = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=7200
                )
                if r.returncode == 0:
                    genome_cfg["index_ready"] = True
                    actions.append(f"built_star_index:{genome_key}")
                    log.info(f"STAR index built: {index_dir}")
                else:
                    warnings.append(
                        f"STAR index build failed: {r.stderr[-200:]}"
                    )
            except subprocess.TimeoutExpired:
                warnings.append("STAR index build timed out (>2h)")
            except Exception as e:
                warnings.append(f"STAR index error: {e}")

    def _ensure_kb_index(self, genome_key: str,
                         genome_cfg: dict,
                         conda_cmd: str,
                         experiment_id: str,
                         actions: list,
                         warnings: list):
        """Auto-build `kb ref` transcriptome index + `t2g.txt` (+ SHA-256 digests) once per assembly."""
        from aria.utils.provenance import hash_file

        kb_dir = Path(
            genome_cfg.get("kb_dir")
            or (Path(genome_cfg.get("genome_dir", ".")) / "kb_index")
        )
        idx_path = Path(genome_cfg.get("kb_index") or (kb_dir / "index.idx"))
        t2g_path = Path(genome_cfg.get("kb_t2g") or (kb_dir / "t2g.txt"))
        cdna_path = kb_dir / "cdna.fa"
        lock_file = Path(genome_cfg.get("genome_dir", kb_dir.parent)) / ".kb_build.lock"

        with self._reference_build_lock(lock_file):
            if (
                idx_path.is_file()
                and idx_path.stat().st_size > 0
                and t2g_path.is_file()
                and t2g_path.stat().st_size > 0
            ):
                genome_cfg["kb_index"] = str(idx_path)
                genome_cfg["kb_t2g"] = str(t2g_path)
                genome_cfg["kb_index_sha256"] = hash_file(idx_path)
                genome_cfg["kb_t2g_sha256"] = hash_file(t2g_path)
                genome_cfg["kb_ready"] = True
                return

            if not genome_cfg.get("fasta_ready") or not genome_cfg.get("gtf_ready"):
                warnings.append(
                    "kb-python index cannot be built: genome FASTA/GTF not ready"
                )
                return

            self.publish_status(
                experiment_id,
                f"Building kb-python transcriptome index for {genome_key} (one-time setup)...",
                0.90,
            )
            kb_dir.mkdir(parents=True, exist_ok=True)
            n_threads = str(min(os.cpu_count() or 4, 16))
            cmd = [
                conda_cmd, "run", "-n", "aria-ingestion-env",
                "kb", "ref",
                "-i", str(idx_path),
                "-g", str(t2g_path),
                "-f1", str(cdna_path),
                "-t", n_threads,
                str(genome_cfg["fasta"]),
                str(genome_cfg["gtf"]),
            ]
            try:
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
                if (
                    r.returncode == 0
                    and idx_path.is_file()
                    and idx_path.stat().st_size > 0
                    and t2g_path.is_file()
                    and t2g_path.stat().st_size > 0
                ):
                    genome_cfg["kb_index"] = str(idx_path)
                    genome_cfg["kb_t2g"] = str(t2g_path)
                    genome_cfg["kb_index_sha256"] = hash_file(idx_path)
                    genome_cfg["kb_t2g_sha256"] = hash_file(t2g_path)
                    genome_cfg["kb_ready"] = True
                    actions.append(f"built_kb_index:{genome_key}")
                    log.info(f"kb-python index built: {idx_path}")
                else:
                    warnings.append(
                        f"kb-python index build failed: {(r.stderr or r.stdout)[-200:]}"
                    )
            except subprocess.TimeoutExpired:
                warnings.append("kb-python index build timed out (>2h)")
            except Exception as exc:
                warnings.append(f"kb-python index error: {exc}")

    # ── Helpers ───────────────────────────────────────────────────────────

    def _has_fastq(self, exp_ctx: dict) -> bool:
        if exp_ctx.get("multiome_fastq_manifest"):
            return True
        all_files = [
            f
            for files in exp_ctx.get("modalities", {}).values()
            if isinstance(files, list)
            for f in files
        ]
        return any(
            str(f).lower().endswith(
                (".fastq.gz", ".fq.gz", ".fastq", ".fq")
            )
            for f in all_files
        )

    def _fatal(self, experiment_id: str,
                warnings: list, message: str) -> dict:
        self.publish_finding(
            experiment_id,
            {"summary": f"Setup failed: {message}"},
            Confidence.INSUFFICIENT,
        )
        return {
            "status":     "error",
            "error_type": "SetupFailed",
            "details":    message,
            "warnings":   warnings,
        }

    def receive(self, message):
        pass
