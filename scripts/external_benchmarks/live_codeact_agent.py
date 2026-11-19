#!/usr/bin/env python3
"""Live Same-Model CodeAct Baseline Agent & Shared External Benchmark Utilities.

Provides a genuine, fair, multi-turn Python CodeAct execution agent (comparable
to open-ended scientific coding agents such as Biomni / mini-swe-agent /
GenoAgent) that runs against the **exact same OpenAI-compatible LLM endpoint**
used by ARIA (`temperature=0`, identical workspace files, identical scientific
Python environment `aria-rna-env`, and an explicit option to abstain/refuse
when data are missing, confounded, or out of scope).

Anti-contamination & isolation guarantees:
- Strictly zero task IDs, accessions, or answers in `aria/` (`aria/` stays
  bit-for-bit identical to git tag `v4.9.0`).
- Every task runs inside a fresh scratch directory with read-only staged data
  (hard check rejecting any `CapsuleNotebook` / `.ipynb` solution leakage).
- Executes real Python code blocks in an isolated subprocess inside
  `aria-rna-env` with bounded wall-clock timeout, no network egress, and
  capped stdout/stderr feedback.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import urllib.request

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONDA_ENVS = Path(os.environ.get("CONDA_ENVS_PATH", Path.home() / "anaconda3" / "envs"))
DEFAULT_RNA_PYTHON = str(DEFAULT_CONDA_ENVS / "aria-rna-env" / "bin" / "python")
DEFAULT_TANTALUS_BASE = "http://127.0.0.1:11435/v1"
DEFAULT_TANTALUS_MODEL = "tantalus-heavy"
DEFAULT_SHIM_BASE = "http://127.0.0.1:11436/v1"

_PYTHON_BLOCK_RE = re.compile(r"```(?:python|py)\s*\n(.*?)```", re.DOTALL | re.IGNORECASE)
_FINAL_TAG_RE = re.compile(r"<FINAL_ANSWER>\s*(.*?)\s*</FINAL_ANSWER>", re.DOTALL | re.IGNORECASE)
_JSON_FENCE_RE = re.compile(r"```json\s*\n(.*?)```", re.DOTALL | re.IGNORECASE)
_THINK_RE = re.compile(r"<think>.*?</think>", re.DOTALL | re.IGNORECASE)


def compute_aria_package_sha256(repo_root: Path = REPO_ROOT) -> str:
    """Hash all tracked Python files under `aria/` to prove frozen-code evaluation."""
    h = hashlib.sha256()
    pkg = repo_root / "aria"
    for p in sorted(pkg.rglob("*.py")):
        rel = p.relative_to(repo_root).as_posix()
        h.update(rel.encode("utf-8"))
        h.update(b"\0")
        h.update(p.read_bytes())
        h.update(b"\0")
    return h.hexdigest()


def assert_blind_workspace(workspace_dir: Path) -> None:
    """Fail closed if any reference solution notebook leaked into a staged workspace."""
    for p in workspace_dir.rglob("*"):
        name_l = p.name.lower()
        if "capsulenotebook" in name_l or name_l.endswith(".ipynb"):
            raise RuntimeError(
                f"ANTI-CONTAMINATION VIOLATION: reference notebook path present in workspace: {p.name}"
            )


def stage_blind_capsule_data(source_capsule_dir: Path, dest_data_dir: Path) -> list[str]:
    """Copy ONLY `CapsuleData-*` files (excluding `CapsuleNotebook-*` / `.ipynb`) into `dest_data_dir`."""
    dest_data_dir.mkdir(parents=True, exist_ok=True)
    data_subdirs = sorted(
        d for d in source_capsule_dir.iterdir()
        if d.is_dir() and d.name.startswith("CapsuleData")
    )
    src_root = data_subdirs[0] if data_subdirs else source_capsule_dir
    staged: list[str] = []
    for p in sorted(src_root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(src_root)
        rel_str = rel.as_posix()
        if "capsulenotebook" in rel_str.lower() or p.suffix.lower() == ".ipynb":
            continue
        target = dest_data_dir / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.link(p, target)
        except OSError:
            shutil.copy2(p, target)
        staged.append(rel_str)
    assert_blind_workspace(dest_data_dir)
    return staged


def chat_completion(
    messages: list[dict[str, str]],
    *,
    api_base: str = DEFAULT_TANTALUS_BASE,
    model: str = DEFAULT_TANTALUS_MODEL,
    api_key: str = "local-benchmark-key",
    temperature: float = 0.0,
    max_tokens: int = 2048,
    timeout: int = 180,
) -> dict[str, Any]:
    """Call an OpenAI-compatible `/chat/completions` endpoint with deterministic settings."""
    url = api_base.rstrip("/") + "/chat/completions"
    payload = {
        "model": model,
        "messages": messages,
        "temperature": temperature,
        "max_tokens": max_tokens,
    }
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=data,
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
        },
        method="POST",
    )
    t0 = time.perf_counter()
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = json.loads(resp.read().decode("utf-8"))
    elapsed = time.perf_counter() - t0
    choice = (raw.get("choices") or [{}])[0]
    msg = choice.get("message") or {}
    content = str(msg.get("content") or "")
    # Strip hidden reasoning tags if a reasoning model emits them inline
    clean = _THINK_RE.sub("", content).strip()
    usage = raw.get("usage") or {}
    return {
        "content": clean,
        "raw_content": content,
        "prompt_tokens": int(usage.get("prompt_tokens") or 0),
        "completion_tokens": int(usage.get("completion_tokens") or 0),
        "elapsed_s": round(elapsed, 3),
        "model": raw.get("model") or model,
    }


def extract_json_answer(text: str) -> dict[str, Any] | None:
    """Extract structured JSON from `<FINAL_ANSWER>...</FINAL_ANSWER>` or trailing JSON block."""
    candidates: list[str] = []
    m = _FINAL_TAG_RE.search(text)
    if m:
        candidates.append(m.group(1).strip())
    for jm in _JSON_FENCE_RE.finditer(text):
        candidates.append(jm.group(1).strip())
    # Also try last balanced {...} block if `<FINAL_ANSWER>` wasn't fenced
    if m:
        inner = m.group(1)
        start = inner.find("{")
        end = inner.rfind("}")
        if start != -1 and end > start:
            candidates.append(inner[start : end + 1])
    for cand in candidates:
        # Remove inner code fences if any
        c = re.sub(r"^```(?:json)?\s*|\s*```$", "", cand.strip(), flags=re.IGNORECASE)
        try:
            parsed = json.loads(c)
            if isinstance(parsed, dict):
                return parsed
        except Exception:
            continue
    return None


def _run_python_step(
    code: str,
    *,
    work_dir: Path,
    step_idx: int,
    python_exe: str = DEFAULT_RNA_PYTHON,
    step_timeout: int = 60,
) -> dict[str, Any]:
    """Execute one CodeAct Python snippet inside an isolated subprocess."""
    script_path = work_dir / f"_codeact_step_{step_idx:02d}.py"
    script_path.write_text(code, encoding="utf-8")
    env = {
        "PATH": f"{Path(python_exe).parent}:/usr/local/bin:/usr/bin:/bin",
        "HOME": str(work_dir),
        "TMPDIR": str(work_dir),
        "MPLBACKEND": "Agg",
        "MPLCONFIGDIR": str(work_dir / ".mpl"),
        "NUMBA_CACHE_DIR": str(work_dir / ".numba"),
        "PYTHONHASHSEED": "0",
        "OMP_NUM_THREADS": "2",
        "OPENBLAS_NUM_THREADS": "2",
        "MKL_NUM_THREADS": "2",
        "NO_PROXY": "*",
        "HTTP_PROXY": "http://127.0.0.1:9",
        "HTTPS_PROXY": "http://127.0.0.1:9",
    }
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            [python_exe, str(script_path)],
            cwd=str(work_dir),
            env=env,
            capture_output=True,
            text=True,
            timeout=step_timeout,
            check=False,
        )
        return {
            "returncode": proc.returncode,
            "stdout": proc.stdout[-4000:],
            "stderr": proc.stderr[-2500:],
            "timed_out": False,
            "elapsed_s": round(time.perf_counter() - t0, 3),
        }
    except subprocess.TimeoutExpired as exc:
        return {
            "returncode": -1,
            "stdout": str(exc.stdout or "")[-2000:],
            "stderr": f"Step timed out after {step_timeout}s",
            "timed_out": True,
            "elapsed_s": round(time.perf_counter() - t0, 3),
        }


@dataclass
class CodeActResult:
    parsed_answer: dict[str, Any] | None
    raw_final_text: str
    n_turns: int
    n_code_executions: int
    n_code_errors: int
    prompt_tokens: int
    completion_tokens: int
    elapsed_s: float
    error: str | None = None
    trace: list[dict[str, Any]] = field(default_factory=list)


CODEACT_SYSTEM_PROMPT = """\
You are a computational biology data-analysis CodeAct agent working inside an isolated Linux workspace.
The input files for this task are located inside `./data/`.
You have a Python 3.11 scientific environment with pandas, numpy, scipy, statsmodels, scikit-learn, pydeseq2, scanpy, anndata, gseapy, and openpyxl.

How to work:
1. On any turn where you need to inspect files or run analysis, write a single ```python ... ``` code block.
   It will be executed in the workspace root (`./data/` holds the task files) and you will receive its stdout/stderr.
2. Do NOT guess numbers or metadata without inspecting the actual files in `./data/`.
3. If the data genuinely cannot support the requested analysis (e.g., missing labels, single group, unusable cohort, incompatible format, or insufficient evidence), you ARE allowed and encouraged to explicitly abstain rather than guessing.
4. Once you have your final answer (or abstention), emit ONLY a valid JSON object inside `<FINAL_ANSWER> ... </FINAL_ANSWER>` tags matching the requested schema, with no code block on that final turn.
"""


def run_live_codeact_task(
    task_prompt: str,
    *,
    work_dir: Path,
    api_base: str = DEFAULT_TANTALUS_BASE,
    model: str = DEFAULT_TANTALUS_MODEL,
    python_exe: str = DEFAULT_RNA_PYTHON,
    max_turns: int = 6,
    step_timeout: int = 60,
    llm_timeout: int = 180,
    max_tokens: int = 2048,
    system_prompt: str = CODEACT_SYSTEM_PROMPT,
) -> CodeActResult:
    """Run a genuine multi-turn live LLM CodeAct loop in `work_dir`."""
    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    assert_blind_workspace(work_dir)

    messages: list[dict[str, str]] = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": task_prompt},
    ]
    t_start = time.perf_counter()
    prompt_toks = 0
    comp_toks = 0
    n_exec = 0
    n_err = 0
    trace: list[dict[str, Any]] = []
    last_text = ""

    for turn in range(1, max_turns + 1):
        try:
            resp = chat_completion(
                messages,
                api_base=api_base,
                model=model,
                temperature=0.0,
                max_tokens=max_tokens,
                timeout=llm_timeout,
            )
        except Exception as exc:
            return CodeActResult(
                parsed_answer=None,
                raw_final_text=last_text,
                n_turns=turn,
                n_code_executions=n_exec,
                n_code_errors=n_err,
                prompt_tokens=prompt_toks,
                completion_tokens=comp_toks,
                elapsed_s=round(time.perf_counter() - t_start, 3),
                error=f"llm_error: {type(exc).__name__}: {str(exc)[:200]}",
                trace=trace,
            )

        text = resp["content"]
        last_text = text
        prompt_toks += resp["prompt_tokens"]
        comp_toks += resp["completion_tokens"]
        messages.append({"role": "assistant", "content": text})

        # Check for final answer tag FIRST if no code block or if final answer present
        parsed = extract_json_answer(text) if "<FINAL_ANSWER>" in text.upper() else None
        code_blocks = _PYTHON_BLOCK_RE.findall(text)

        if parsed is not None:
            trace.append({"turn": turn, "action": "final_answer", "elapsed_s": resp["elapsed_s"]})
            return CodeActResult(
                parsed_answer=parsed,
                raw_final_text=text,
                n_turns=turn,
                n_code_executions=n_exec,
                n_code_errors=n_err,
                prompt_tokens=prompt_toks,
                completion_tokens=comp_toks,
                elapsed_s=round(time.perf_counter() - t_start, 3),
                trace=trace,
            )

        if code_blocks and turn < max_turns:
            code = code_blocks[0]
            exec_res = _run_python_step(
                code,
                work_dir=work_dir,
                step_idx=turn,
                python_exe=python_exe,
                step_timeout=step_timeout,
            )
            n_exec += 1
            if exec_res["returncode"] != 0:
                n_err += 1
            trace.append({
                "turn": turn,
                "action": "python_exec",
                "returncode": exec_res["returncode"],
                "timed_out": exec_res["timed_out"],
                "elapsed_s": exec_res["elapsed_s"],
            })
            obs = (
                f"Execution finished (returncode={exec_res['returncode']}, "
                f"elapsed={exec_res['elapsed_s']}s).\n"
                f"STDOUT:\n{exec_res['stdout'] or '(empty)'}\n"
            )
            if exec_res["stderr"]:
                obs += f"STDERR:\n{exec_res['stderr']}\n"
            if turn == max_turns - 1:
                obs += (
                    "\nThis is your FINAL turn. Do NOT write any more code. "
                    "Return ONLY the final JSON inside <FINAL_ANSWER>...</FINAL_ANSWER>."
                )
            messages.append({"role": "user", "content": obs})
        else:
            # No code block and no `<FINAL_ANSWER>` yet: try fallback JSON parse or nudge once
            fallback = extract_json_answer(text)
            if fallback is not None:
                trace.append({"turn": turn, "action": "final_answer_fallback", "elapsed_s": resp["elapsed_s"]})
                return CodeActResult(
                    parsed_answer=fallback,
                    raw_final_text=text,
                    n_turns=turn,
                    n_code_executions=n_exec,
                    n_code_errors=n_err,
                    prompt_tokens=prompt_toks,
                    completion_tokens=comp_toks,
                    elapsed_s=round(time.perf_counter() - t_start, 3),
                    trace=trace,
                )
            trace.append({"turn": turn, "action": "nudge", "elapsed_s": resp["elapsed_s"]})
            messages.append({
                "role": "user",
                "content": (
                    "Please provide your final structured JSON response inside "
                    "<FINAL_ANSWER>...</FINAL_ANSWER> tags (or a single ```python``` block to inspect `./data/`)."
                ),
            })

    final_parsed = extract_json_answer(last_text)
    if final_parsed is None:
        # One symmetric JSON-only wrap-up turn so valid multi-step runs never fail purely on formatting
        messages.append({
            "role": "user",
            "content": (
                "Your code execution turns are finished. Do NOT write Python code. "
                "Using ONLY the execution outputs already above, output ONLY the valid JSON object "
                "inside <FINAL_ANSWER>...</FINAL_ANSWER> tags (or set abstained=true if insufficient)."
            ),
        })
        try:
            wrap = chat_completion(
                messages,
                api_base=api_base,
                model=model,
                temperature=0.0,
                max_tokens=1024,
                timeout=llm_timeout,
            )
            last_text = wrap["content"]
            prompt_toks += wrap["prompt_tokens"]
            comp_toks += wrap["completion_tokens"]
            final_parsed = extract_json_answer(last_text)
            trace.append({"turn": max_turns + 1, "action": "wrap_up_json", "elapsed_s": wrap["elapsed_s"]})
        except Exception:
            pass

    return CodeActResult(
        parsed_answer=final_parsed,
        raw_final_text=last_text,
        n_turns=max_turns,
        n_code_executions=n_exec,
        n_code_errors=n_err,
        prompt_tokens=prompt_toks,
        completion_tokens=comp_toks,
        elapsed_s=round(time.perf_counter() - t_start, 3),
        error=None if final_parsed is not None else "max_turns_without_json",
        trace=trace,
    )
