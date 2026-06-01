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
import os
import shutil
import signal
import subprocess
import uuid
import fcntl
from pathlib import Path
from typing import Any
from aria.utils.provenance import hash_params
from aria.utils.privacy import redact_sensitive_params
from aria.utils.script_contracts import ContractIssue, contract_for_script

shutil_rmtree = shutil.rmtree

log = logging.getLogger("aria.env")


# ARIA package root — this file is at /aria/utils/environment_manager.py
# so the package root is two levels up. Used to resolve script_path values
# like "aria/scripts/rna_qc.py" reliably regardless of the caller's CWD.
# (Without this, Path("aria/scripts/...").resolve() depends on os.getcwd(),
# which produces wrong paths like /Samael/ARIA/aria/scripts/aria/scripts/...
# when ARIA is launched from inside aria/ or aria/scripts/.)
_ARIA_PACKAGE_ROOT = Path(__file__).resolve().parent.parent.parent


def _resolve_script_path(script_path: str) -> Path:
    """
    Resolve a script_path string to an absolute Path.

    Acceptable inputs:
      - "aria/scripts/rna_qc.py"      (relative to ARIA package root)
      - "scripts/rna_qc.py"           (relative — also tried under aria/)
      - "/abs/path/to/script.py"      (absolute — used as-is)

    Returns the first form that resolves to an existing file. If none
    exist, returns the package-root-relative form (so the error message
    points at a sensible location).
    """
    p = Path(script_path)
    if p.is_absolute():
        return p

    candidates = [
        _ARIA_PACKAGE_ROOT / script_path,           # aria/scripts/rna_qc.py
        _ARIA_PACKAGE_ROOT / "aria" / script_path,  # in case "scripts/..." form
        Path(script_path).resolve(),                # CWD-relative as last resort
    ]
    for c in candidates:
        if c.exists():
            return c
    return candidates[0]   # for the error message


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
        "rnaseq":      "aria-rnaseq-env",   # raw FASTQ processing: fastp/STAR/featureCounts
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
        "hic":         14400,
        "rnaseq":      10800,   # 3h — STAR alignment can be slow   # 4 hours
        "integration": 7200,    # 2 hours
        "spatial":     3600,    # 1 hour
    }

    FALLBACK_ENV = "aria-env"   # Base environment used when stack env is missing

    # Number of failed runs to keep on disk for postmortem inspection.
    # Files live in <workspace>/failed/ — older runs are evicted FIFO.
    MAX_FAILED_RUNS = 20

    def __init__(self, workspace_dir: str = "~/.aria/workspace"):
        self.workspace = Path(workspace_dir).expanduser()
        self.workspace.mkdir(parents=True, exist_ok=True)
        (self.workspace / "failed").mkdir(exist_ok=True)
        self._conda_ok = self._verify_conda()
        # Cache for check_environments() — enumerating conda envs spawns a
        # subprocess, and run_in_stack used to do it on every dispatch. The env
        # set is fixed for the life of the process, so this is memoized once.
        self._env_cache: dict[str, bool] | None = None
        self._env_cache_key: str | None = None

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

        run_id      = str(uuid.uuid4())[:8]
        input_file  = self.workspace / f"input_{run_id}.json"
        output_file = self.workspace / f"output_{run_id}.json"
        params_sha256 = hash_params(params)
        max_time    = max(timeout or self.TIMEOUTS.get(stack, 3600), 120)  # min 2 min

        result: dict = {"status": "error"}
        succeeded = False
        try:
            # 1. Build command — resolve script path against package root
            # (CWD-based resolution caused duplicate path bugs like
            # /aria/scripts/aria/scripts/foo.py when launched from aria/scripts/)
            resolved_script = _resolve_script_path(script_path)
            if not resolved_script.exists():
                result = {
                    "status":     "error",
                    "error_type": "ScriptNotFound",
                    "details": (
                        f"Script not found: {script_path}\n"
                        f"Resolved to: {resolved_script}\n"
                        f"ARIA package root: {_ARIA_PACKAGE_ROOT}"
                    ),
                }
                return result

            contract = contract_for_script(resolved_script)
            if contract is not None:
                issues = contract.validate_params(params)
                if issues:
                    result = self._contract_error(
                        stage="input",
                        contract_path=contract.script_path,
                        issues=issues,
                    )
                    return result

            # 2. Write parameters to input file after schema validation. The
            # live IPC file needs real paths, but failed-run archives redact it
            # by default below (C6/X10 privacy).
            with open(input_file, "w") as f:
                json.dump(params, f)
            try:
                input_file.chmod(0o600)
            except OSError:
                pass
            (self.workspace / f"params_{run_id}.sha256").write_text(
                params_sha256 + "\n"
            )

            env_name = self._resolve_env(stack)
            cmd = [
                "conda", "run",
                "--no-capture-output",   # let heavy C logs go to system stderr
                "-n", env_name,
                "python",
                str(resolved_script),
                str(input_file),
                str(output_file),
            ]

            log.debug(f"Running {resolved_script.name} in {env_name} "
                      f"(timeout={max_time}s)")

            # 3. Execute in a NEW SESSION so the whole process tree can be killed
            # on timeout (R5, audit 2026-05-29). `conda run` spawns `python`,
            # which spawns BLAS/numba threads; `subprocess.run(timeout=)` only
            # kills the direct child (`conda run`), orphaning the grandchildren
            # on large datasets. start_new_session makes the child a process
            # group leader so a timeout can reap the entire group.
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                start_new_session=True,
            )
            try:
                stdout, stderr = proc.communicate(timeout=max_time)
            except subprocess.TimeoutExpired:
                self._terminate_process_tree(proc)
                # Drain whatever the killed tree already wrote, then surface the
                # timeout through the handler below.
                try:
                    proc.communicate(timeout=10)
                except Exception:
                    pass
                raise
            process = subprocess.CompletedProcess(
                cmd, proc.returncode, stdout, stderr
            )

            # 4. Check for subprocess failure
            if process.returncode != 0:
                log.error(
                    f"Subprocess failed in {env_name} "
                    f"(exit {process.returncode}): {process.stderr[-500:]}"
                )
                result = {
                    "status":     "error",
                    "error_type": "SubprocessFailed",
                    "exit_code":  process.returncode,
                    "details":    process.stderr[-1000:],
                }
                return result

            # 5. Read structured output
            if not output_file.exists():
                result = {
                    "status":     "error",
                    "error_type": "MissingOutput",
                    "details":    (
                        f"Script {Path(script_path).name} exited cleanly "
                        f"but produced no output JSON."
                    ),
                }
                return result

            with open(output_file, "r") as f:
                result = json.load(f)
            result.setdefault("params_sha256", params_sha256)
            if contract is not None:
                issues = contract.validate_result(result)
                if issues:
                    result = self._contract_error(
                        stage="output",
                        contract_path=contract.script_path,
                        issues=issues,
                        params_sha256=params_sha256,
                    )
                    return result
                result.setdefault("ipc_contract", {
                    "script_path": contract.script_path,
                    "contract_version": contract.contract_version,
                    "validation_level": contract.validation_level,
                })

            succeeded = (result.get("status") == "success")
            return result

        except subprocess.TimeoutExpired:
            log.error(f"Stack '{stack}' timed out after {max_time}s")
            result = {
                "status":     "error",
                "error_type": "Timeout",
                "details":    f"Execution exceeded {max_time}s limit.",
            }
            return result

        except Exception as e:
            log.error(f"EnvironmentManager exception: {e}")
            result = {
                "status":     "error",
                "error_type": type(e).__name__,
                "details":    str(e),
            }
            return result

        finally:
            if succeeded:
                # Happy path: drop the temp files immediately.
                for f in (input_file, output_file):
                    if f.exists():
                        try:
                            f.unlink()
                        except OSError:
                            pass
                sha_file = self.workspace / f"params_{run_id}.sha256"
                if sha_file.exists():
                    try:
                        sha_file.unlink()
                    except OSError:
                        pass
            else:
                # Failure: preserve input + output for postmortem.
                result.setdefault("params_sha256", params_sha256)
                self._archive_failed_run(run_id, stack, input_file, output_file, result)

    @staticmethod
    def _contract_error(
        stage: str,
        contract_path: str,
        issues: list[ContractIssue],
        params_sha256: str | None = None,
    ) -> dict[str, Any]:
        issue_dicts = [issue.model_dump() for issue in issues]
        version_issue = any(issue.field in {
            "_aria_contract_version",
            "contract_version",
            "ipc_contract_version",
        } for issue in issues)
        result = {
            "status": "error",
            "error_type": (
                "IncompatibleScriptContract"
                if version_issue or stage == "output"
                else "InvalidScriptParams"
            ),
            "details": (
                f"{stage} IPC contract validation failed for "
                f"{contract_path}: "
                + "; ".join(issue.message for issue in issues)
            ),
            "contract_stage": stage,
            "script_path": contract_path,
            "contract_issues": issue_dicts,
        }
        if params_sha256:
            result["params_sha256"] = params_sha256
        return result

    def _archive_failed_run(self, run_id: str, stack: str,
                             input_file: Path, output_file: Path,
                             result: dict):
        """
        Move input/output JSON of a failed run to <workspace>/failed/ so the
        user can inspect what went wrong. Evicts the oldest archived runs to
        cap disk usage at MAX_FAILED_RUNS.
        """
        failed_dir = self.workspace / "failed"
        run_dir    = failed_dir / f"{stack}_{run_id}"
        # P1-8 (c) / ADR-024: archived failed-run artifacts are redacted by
        # default so a sensitive run never leaves raw paths/params/outputs on
        # disk. The same opt-in that preserves raw inputs preserves raw
        # output/error for debugging.
        preserve_raw = (
            os.environ.get("ARIA_PRESERVE_FAILED_INPUTS", "0")
            .strip().lower() in {"1", "true", "yes", "on"}
        )

        def _archive_json(payload, dest_raw: Path, dest_redacted: Path):
            """Write `payload` redacted (default) or raw (opt-in)."""
            if preserve_raw:
                dest_raw.write_text(json.dumps(payload, indent=2),
                                    encoding="utf-8")
                return
            try:
                redacted = redact_sensitive_params(payload)
                dest_redacted.write_text(json.dumps(redacted, indent=2),
                                         encoding="utf-8")
            except Exception:
                dest_redacted.write_text(json.dumps({"redacted": True}, indent=2),
                                         encoding="utf-8")

        try:
            run_dir.mkdir(parents=True, exist_ok=True)
            if input_file.exists():
                if preserve_raw:
                    input_file.replace(run_dir / "input.json")
                else:
                    try:
                        with input_file.open("r", encoding="utf-8") as fh:
                            raw_params = json.load(fh)
                        _archive_json(raw_params, run_dir / "input.json",
                                      run_dir / "input.redacted.json")
                    except Exception:
                        (run_dir / "input.redacted.json").write_text(
                            json.dumps({"redacted": True}, indent=2),
                            encoding="utf-8",
                        )
                    try:
                        input_file.unlink()
                    except OSError:
                        pass
            if output_file.exists():
                if preserve_raw:
                    output_file.replace(run_dir / "output.json")
                else:
                    try:
                        with output_file.open("r", encoding="utf-8") as fh:
                            raw_output = json.load(fh)
                        _archive_json(raw_output, run_dir / "output.json",
                                      run_dir / "output.redacted.json")
                    except Exception:
                        (run_dir / "output.redacted.json").write_text(
                            json.dumps({"redacted": True}, indent=2),
                            encoding="utf-8",
                        )
                    try:
                        output_file.unlink()
                    except OSError:
                        pass
            # Drop the error summary next to the JSON for easy triage (redacted
            # by default — the result dict can echo paths/params).
            _archive_json(result, run_dir / "error.json",
                          run_dir / "error.redacted.json")
            sha_file = self.workspace / f"params_{run_id}.sha256"
            if sha_file.exists():
                sha_file.replace(run_dir / "params.sha256")
        except Exception as e:
            log.warning(f"Could not archive failed run {run_id}: {e}")
            return

        # FIFO eviction so the failed/ dir stays bounded.
        lock_file = None
        try:
            lock_path = failed_dir / ".lock"
            lock_file = lock_path.open("a+")
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
            archived = sorted(
                (p for p in failed_dir.iterdir()
                 if p.is_dir() and p.name != ".lock"),
                key=lambda p: p.stat().st_mtime,
            )
            while len(archived) > self.MAX_FAILED_RUNS:
                oldest = archived.pop(0)
                shutil_rmtree(oldest)
        except Exception as e:
            log.debug(f"Failed-run eviction skipped: {e}")
        finally:
            if lock_file is not None:
                try:
                    fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
                    lock_file.close()
                except Exception:
                    pass

    def check_environments(self) -> dict[str, bool]:
        """
        Check which ARIA Conda environments are installed on this system.

        Returns:
            dict mapping stack name -> True if env exists, False otherwise
        """
        if not self._conda_ok:
            return {stack: False for stack in self.STACKS}

        # Memoize per process (F-ENG-PERF): enumerating conda envs spawns a
        # subprocess, and run_in_stack used to do it on every dispatch. ARIA's
        # env set does not change mid-run (SetupAgent installs before dispatch,
        # by name, with no aliasing — see B12), so a stable cache key is enough.
        cache_key = "process"
        if self._env_cache is not None and self._env_cache_key == cache_key:
            return dict(self._env_cache)

        try:
            result = subprocess.run(
                ["conda", "env", "list", "--json"],
                capture_output=True, text=True, check=True,
            )
            data         = json.loads(result.stdout)
            installed    = {Path(e).name for e in data.get("envs", [])}
            envs = {
                stack: self.STACKS[stack] in installed
                for stack in self.STACKS
            }
            self._env_cache = dict(envs)
            self._env_cache_key = cache_key
            return envs
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
        Resolve which conda environment to use for a stack.

        Resolution order:
          1. Preferred env name (STACKS[stack]) if installed
          2. FALLBACK_ENV (aria-env / base)

        SetupAgent owns ARIA's environments by name and does NOT do tool-aware
        detection or aliasing of pre-existing user environments (B12, audit
        2026-05-29). There is no `env_aliases.json` writer in the codebase, so
        the former alias-resolution branch was dead; it was removed so this
        method's behavior matches SetupAgent's documented "no aliases" policy.
        """
        env_name  = self.STACKS.get(stack, self.FALLBACK_ENV)
        available = self.check_environments()

        if not available.get(stack, False):
            log.warning(
                f"Environment '{env_name}' not found. "
                f"Falling back to '{self.FALLBACK_ENV}'. "
                f"ARIA will attempt automatic setup on next run."
            )
            return self.FALLBACK_ENV

        return env_name

    @staticmethod
    def _terminate_process_tree(proc: "subprocess.Popen") -> None:
        """
        Kill the entire process group led by `proc` (R5). The child was started
        with `start_new_session=True`, so it is its own group leader and
        `os.killpg` reaches `conda run` plus the real `python` and its
        BLAS/numba threads. SIGTERM first, then SIGKILL if it lingers; falls
        back to killing the direct child if the group is already gone.
        """
        try:
            pgid = os.getpgid(proc.pid)
        except (ProcessLookupError, OSError):
            pgid = None

        if pgid is not None:
            try:
                os.killpg(pgid, signal.SIGTERM)
                try:
                    proc.wait(timeout=5)
                    return
                except subprocess.TimeoutExpired:
                    os.killpg(pgid, signal.SIGKILL)
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        pass
                    return
            except (ProcessLookupError, PermissionError, OSError):
                pass

        # Fallback: at least kill the direct child.
        try:
            proc.kill()
        except Exception:
            pass

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


class _LazyEnvironmentManager:
    """Create the global EnvironmentManager only when it is first used."""

    def __init__(self):
        self._instance: EnvironmentManager | None = None
        self._lock = None

    def _get(self) -> EnvironmentManager:
        if self._instance is None:
            import threading
            if self._lock is None:
                self._lock = threading.Lock()
            with self._lock:
                if self._instance is None:
                    self._instance = EnvironmentManager()
        return self._instance

    def __getattr__(self, name: str):
        return getattr(self._get(), name)


# Global accessor imported by agents; construction is lazy to avoid import-time
# workspace creation and environment probing.
env_manager = _LazyEnvironmentManager()
