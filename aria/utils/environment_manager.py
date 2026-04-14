"""
ARIA EnvironmentManager
-----------------------
Isolates bioinformatics tool execution in dedicated Conda environments.

Prevents agents from importing libraries directly, eliminating the
"dependency hell" problem (e.g. C-level conflicts between scanpy,
cooler, MACS3, pysam across different stacks).

IPC strategy: JSON file-based communication.
  - Agent writes params  -> input_{id}.json
  - Subprocess runs      -> conda run -n <env> python <script> in out
  - Agent reads results  <- output_{id}.json

Stdout/stderr are NOT used for data transfer — bioinformatics tools
pollute stdout with C library warnings that break pipe-based parsing.
"""

from __future__ import annotations

import json
import logging
import subprocess
import uuid
from pathlib import Path
from typing import Any

log = logging.getLogger("aria.env")


class EnvironmentManager:
    """
    Manages isolated subprocess execution per analytical stack.

    Each stack maps to a dedicated Conda environment with pinned
    dependency versions, preventing cross-contamination between
    tools like scanpy, cooler, MACS3, and pysam.

    Usage:
        result = env_manager.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_qc.py",
            params={"data_path": "/data/sample.h5ad", "organism": "Homo sapiens"}
        )
        if result["status"] == "success":
            print(result["n_cells_after_qc"])
    """

    # Analytical stack -> Conda environment name
    STACKS: dict[str, str] = {
        "rna":         "aria-rna-env",
        "chromatin":   "aria-chromatin-env",
        "hic":         "aria-hic-env",
        "integration": "aria-integration-env",
        "spatial":     "aria-rna-env",   # spatial uses same env as RNA
    }

    # Per-stack execution timeouts (seconds)
    # HiC and integration can take hours on large datasets
    TIMEOUTS: dict[str, int] = {
        "rna":         3600,    # 1 hour
        "chromatin":   7200,    # 2 hours
        "hic":         14400,   # 4 hours
        "integration": 7200,    # 2 hours
        "spatial":     3600,    # 1 hour
    }

    FALLBACK_ENV = "aria-env"   # Base environment used when stack env is missing

    def __init__(self, workspace_dir: str = "~/.aria/workspace"):
        self.workspace = Path(workspace_dir).expanduser()
        self.workspace.mkdir(parents=True, exist_ok=True)
        self._conda_ok = self._verify_conda()

    # ── Public interface ──────────────────────────────────────────────────

    def run_in_stack(
        self,
        stack:       str,
        script_path: str,
        params:      dict[str, Any],
        timeout:     int = None,
    ) -> dict[str, Any]:
        """
        Execute a Python script inside the appropriate isolated Conda env.

        Args:
            stack:       Analytical stack key ("rna", "chromatin", "hic", etc.)
            script_path: Path to the analysis script (must follow _base.py contract)
            params:      Parameters dict written to input JSON
            timeout:     Override default timeout in seconds (None = use stack default)

        Returns:
            dict with at minimum {"status": "success" | "error"}
            On success: additional keys depend on the script
            On error:   {"status": "error", "error_type": str, "details": str}
        """
        if stack not in self.STACKS:
            return {
                "status":     "error",
                "error_type": "UnknownStack",
                "details":    f"Unknown stack: '{stack}'. Valid: {list(self.STACKS)}",
            }

        env_name    = self._resolve_env(stack)
        run_id      = str(uuid.uuid4())[:8]
        input_file  = self.workspace / f"input_{run_id}.json"
        output_file = self.workspace / f"output_{run_id}.json"
        max_time    = timeout or self.TIMEOUTS.get(stack, 3600)

        try:
            # 1. Write parameters to input file
            with open(input_file, "w") as f:
                json.dump(params, f)

            # 2. Build command
            cmd = [
                "conda", "run",
                "--no-capture-output",   # let heavy C logs go to system stderr
                "-n", env_name,
                "python",
                str(Path(script_path).resolve()),
                str(input_file),
                str(output_file),
            ]

            log.debug(f"Running {Path(script_path).name} in {env_name} "
                      f"(timeout={max_time}s)")

            # 3. Execute
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=max_time,
            )

            # 4. Check for subprocess failure
            if process.returncode != 0:
                log.error(
                    f"Subprocess failed in {env_name} "
                    f"(exit {process.returncode}): {process.stderr[-500:]}"
                )
                return {
                    "status":     "error",
                    "error_type": "SubprocessFailed",
                    "exit_code":  process.returncode,
                    "details":    process.stderr[-1000:],
                }

            # 5. Read structured output
            if not output_file.exists():
                return {
                    "status":     "error",
                    "error_type": "MissingOutput",
                    "details":    (
                        f"Script {Path(script_path).name} exited cleanly "
                        f"but produced no output JSON."
                    ),
                }

            with open(output_file, "r") as f:
                result = json.load(f)

            return result

        except subprocess.TimeoutExpired:
            log.error(f"Stack '{stack}' timed out after {max_time}s")
            return {
                "status":     "error",
                "error_type": "Timeout",
                "details":    f"Execution exceeded {max_time}s limit.",
            }

        except Exception as e:
            log.error(f"EnvironmentManager exception: {e}")
            return {
                "status":     "error",
                "error_type": type(e).__name__,
                "details":    str(e),
            }

        finally:
            # Always clean up temp files
            for f in (input_file, output_file):
                if f.exists():
                    try:
                        f.unlink()
                    except OSError:
                        pass

    def check_environments(self) -> dict[str, bool]:
        """
        Check which ARIA Conda environments are installed on this system.

        Returns:
            dict mapping stack name -> True if env exists, False otherwise
        """
        if not self._conda_ok:
            return {stack: False for stack in self.STACKS}

        try:
            result = subprocess.run(
                ["conda", "env", "list", "--json"],
                capture_output=True, text=True, check=True,
            )
            data         = json.loads(result.stdout)
            installed    = {Path(e).name for e in data.get("envs", [])}
            return {
                stack: self.STACKS[stack] in installed
                for stack in self.STACKS
            }
        except Exception as e:
            log.warning(f"Could not enumerate Conda environments: {e}")
            return {stack: False for stack in self.STACKS}

    def get_status_report(self) -> dict:
        """
        Return a human-readable status report for the TUI.
        Shows which environments are ready and which need installation.
        """
        envs    = self.check_environments()
        ready   = [s for s, ok in envs.items() if ok]
        missing = [s for s, ok in envs.items() if not ok]
        return {
            "conda_available": self._conda_ok,
            "environments":    envs,
            "ready_stacks":    ready,
            "missing_stacks":  missing,
            "fallback_active": bool(missing),
        }

    # ── Private methods ───────────────────────────────────────────────────

    def _resolve_env(self, stack: str) -> str:
        """
        Return the Conda env name for a stack.
        Falls back to FALLBACK_ENV if the dedicated env is not installed.
        """
        env_name = self.STACKS[stack]
        available = self.check_environments()

        if not available.get(stack, False):
            log.warning(
                f"Environment '{env_name}' not found. "
                f"Falling back to '{self.FALLBACK_ENV}'. "
                f"Install with: conda env create -f envs/aria-{stack}-env.yml"
            )
            return self.FALLBACK_ENV

        return env_name

    def _verify_conda(self) -> bool:
        """Check that conda is installed and accessible in PATH."""
        try:
            subprocess.run(
                ["conda", "--version"],
                capture_output=True, check=True,
            )
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            log.warning(
                "Conda not found in PATH. "
                "Isolated environment execution will not work. "
                "Install Miniforge: https://github.com/conda-forge/miniforge"
            )
            return False


# Global instance — imported by all agents
env_manager = EnvironmentManager()
