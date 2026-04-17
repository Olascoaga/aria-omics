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
# These are the authoritative env specs ARIA uses.
# Users never edit these. ARIA installs them on first run.

ARIA_ENVS = {
    "aria-rnaseq-env": {
        "description": "Raw RNA-seq: fastp, STAR, featureCounts, samtools",
        "stack":       "rnaseq",
        "yml": """name: aria-rnaseq-env
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - fastp>=1.0
  - star>=2.7.11
  - subread>=2.1
  - samtools>=1.21
  - multiqc>=1.25
  - fastqc>=0.12
  - pandas>=1.5
  - numpy<2.0.0
""",
    },
    "aria-rna-env": {
        "description": "scRNA-seq and bulk DE: scanpy, pydeseq2, gseapy",
        "stack":       "rna",
        "yml": """name: aria-rna-env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - scanpy>=1.10
  - anndata>=0.10
  - leidenalg
  - igraph
  - umap-learn
  - numba
  - llvmlite
  - pip:
    - pydeseq2
    - gseapy
    - scrublet
""",
    },
    "aria-chromatin-env": {
        "description": "Chromatin analysis: pysam, MACS3, episcanpy",
        "stack":       "chromatin",
        "yml": """name: aria-chromatin-env
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - pysam>=0.21
  - samtools>=1.21
  - numpy<2.0.0
  - pandas>=1.5
  - pip:
    - macs3
    - episcanpy
""",
    },
    "aria-hic-env": {
        "description": "3D genome: cooler, cooltools",
        "stack":       "hic",
        "yml": """name: aria-hic-env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - cooler>=0.9
  - h5py>=3.8
  - hdf5=1.14.*
  - numpy<2.0.0
  - pandas>=1.5
  - pip:
    - cooltools
""",
    },
}


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
        needed_envs = self._needed_envs(modalities, has_fastq)
        for env_name in needed_envs:
            self._ensure_env(
                env_name, conda_cmd,
                experiment_id, actions, warnings
            )

        # ── 3. Reference genome + GTF ─────────────────────────────────────
        genome_key = self._detect_genome(organism, intent)
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

            # ── 4. STAR index ─────────────────────────────────────────────
            if has_fastq and genome_cfg.get("fasta_ready"):
                self.publish_status(
                    experiment_id,
                    f"Checking STAR index ({genome_key})...",
                    0.8,
                )
                self._ensure_star_index(
                    genome_key, genome_cfg,
                    conda_cmd, experiment_id, actions, warnings
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

    def _needed_envs(self, modalities: dict, has_fastq: bool) -> list[str]:
        """Which ARIA environments does this experiment need?"""
        envs = []
        if has_fastq:
            envs.append("aria-rnaseq-env")
        if {"scRNA", "bulk_RNA"} & set(modalities):
            envs.append("aria-rna-env")
        if {"scATAC", "bulk_ATAC", "ChIP",
            "CUT_AND_RUN", "CUT_AND_TAG"} & set(modalities):
            envs.append("aria-chromatin-env")
        if {"HiC", "Micro-C"} & set(modalities):
            envs.append("aria-hic-env")
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

        import tempfile
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yml", delete=False
        ) as f:
            f.write(env_def["yml"])
            yml_path = f.name

        try:
            r = subprocess.run(
                [conda_cmd, "env", "create",
                 "-f", yml_path, "--name", env_name, "-y"],
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
        finally:
            Path(yml_path).unlink(missing_ok=True)

    # ── Reference genome ─────────────────────────────────────────────────

    def _detect_genome(self, organism: str,
                        intent: dict) -> Optional[str]:
        """Infer genome key from organism name or biological question."""
        text = (organism + " " + intent.get("summary", "")).lower()

        for key, info in GENOME_REGISTRY.items():
            if any(alias in text for alias in info["aliases"]):
                return key

        # Fallback heuristics
        if any(k in text for k in ["human", "homo", "patient", "clinical"]):
            return "hg38"
        if any(k in text for k in ["mouse", "murine", "mus", "bmal", "circadian"]):
            return "mm39"
        if any(k in text for k in ["rat", "rattus"]):
            return "rn7"

        return None

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
            "fasta_ready": False,
            "gtf_ready":   False,
            "index_ready": False,
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

    # ── STAR index ────────────────────────────────────────────────────────

    def _ensure_star_index(self, genome_key: str,
                            genome_cfg: dict,
                            conda_cmd: str,
                            experiment_id: str,
                            actions: list,
                            warnings: list):
        """Build STAR index if not present. ~30 min, one time only."""
        index_dir = Path(genome_cfg["star_index"])
        ready_files = ["SA", "SAindex", "Genome"]

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

    # ── Helpers ───────────────────────────────────────────────────────────

    def _has_fastq(self, exp_ctx: dict) -> bool:
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
