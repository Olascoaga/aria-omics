#!/usr/bin/env python3
"""Loopback OpenAI-compatible HTTP shim wrapping local WSL system CLIs (`claude`, `codex`, `agy`).

Exposes `GET /health` and `POST /v1/chat/completions` strictly on `127.0.0.1` so
`aria.llm.provider.LLMProvider` can benchmark ARIA end-to-end across Claude,
OpenAI/Codex, and Google Gemini (Antigravity `agy`) in Pure Text-Completion Mode
(empty temporary working directory, zero file/tool access, paid env keys stripped)
at `$0.00` out-of-pocket API spend.

Dual Token & Cost Accounting:
- Returns **net ARIA prompt/completion tokens** (`tiktoken` `cl100k_base` over
  ARIA's exact `messages` and `response`) in the OpenAI `usage` response so
  ARIA's `llm_usage.jsonl` measures pure ARIA token demand without CLI-harness
  system-prompt inflation.
- Appends full per-call dual telemetry (both net ARIA tokens + projected list-price
  USD AND raw CLI `usage` + `cli_reported_cost_usd`) to `--ledger-path`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any


# Official published per-1M-token USD list tariffs (matched to LiteLLM cost map)
PUBLISHED_TARIFF_PER_MTOK: dict[str, tuple[float, float]] = {
    "wsl-claude-sonnet-5-5": (2.00, 10.00),
    "wsl-codex-gpt-5-6-sol": (4.00, 20.00),
    "wsl-agy-gemini-4-carbon-low": (1.25, 10.00),
}

SUPPORTED_MODELS = tuple(PUBLISHED_TARIFF_PER_MTOK.keys())

_LEDGER_LOCK = threading.Lock()
_TIKTOKEN_ENC = None


def _count_tokens(text: str) -> int:
    global _TIKTOKEN_ENC
    if not text:
        return 0
    try:
        if _TIKTOKEN_ENC is None:
            import tiktoken

            _TIKTOKEN_ENC = tiktoken.get_encoding("cl100k_base")
        return len(_TIKTOKEN_ENC.encode(text, disallowed_special=()))
    except Exception:
        return max(1, len(text) // 4)


def _clean_subprocess_env() -> dict[str, str]:
    stripped = {
        "ANTHROPIC_API_KEY",
        "OPENAI_API_KEY",
        "GEMINI_API_KEY",
        "GOOGLE_API_KEY",
    }
    return {k: v for k, v in os.environ.items() if k not in stripped}


def _extract_prompts(messages: list[dict[str, Any]]) -> tuple[str, str]:
    sys_parts: list[str] = []
    conv_parts: list[str] = []
    for m in messages or []:
        role = str(m.get("role") or "user").lower()
        content = m.get("content") or ""
        if isinstance(content, list):
            content = "\n".join(
                str(part.get("text") or "") if isinstance(part, dict) else str(part)
                for part in content
            )
        text = str(content)
        if role == "system":
            if text.strip():
                sys_parts.append(text)
        elif role == "assistant":
            conv_parts.append(f"[Previous Assistant Response]\n{text}")
        else:
            conv_parts.append(text)
    return "\n\n".join(sys_parts).strip(), "\n\n".join(conv_parts).strip()


def _format_combined_prompt(system_text: str, user_text: str) -> str:
    header = (
        "[PURE TEXT COMPLETION MODE: Do NOT call any tools, shell commands, or "
        "file readers. Act strictly as a stateless text/JSON language model and "
        "reply directly to the prompt below.]"
    )
    if system_text:
        return (
            f"{header}\n\n"
            f"<SYSTEM_INSTRUCTIONS>\n{system_text}\n</SYSTEM_INSTRUCTIONS>\n\n"
            f"<USER_PROMPT>\n{user_text}\n</USER_PROMPT>"
        )
    return f"{header}\n\n{user_text}"


def _invoke_claude_wsl(system_text: str, user_text: str, timeout_s: float) -> tuple[str, dict[str, Any]]:
    claude_bin = shutil.which("claude") or "/usr/bin/claude"
    safe_user = user_text or "Proceed."
    if safe_user.startswith("-"):
        safe_user = " \n" + safe_user
    cmd = [
        claude_bin,
        "-p",
        safe_user,
        "--model",
        "claude-sonnet-5-5",
        "--effort",
        "low",
        "--tools",
        "",
        "--safe-mode",
        "--no-session-persistence",
        "--output-format",
        "json",
    ]
    if system_text:
        safe_sys = (" \n" + system_text) if system_text.startswith("-") else system_text
        cmd.extend(["--system-prompt", safe_sys])

    with tempfile.TemporaryDirectory(prefix="aria_wsl_claude_") as td:
        proc = subprocess.run(
            cmd,
            cwd=td,
            env=_clean_subprocess_env(),
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
    if proc.returncode != 0:
        raise RuntimeError(
            f"claude -p failed (rc={proc.returncode}): {(proc.stderr or proc.stdout)[:300]}"
        )
    data = json.loads(proc.stdout)
    if data.get("is_error"):
        raise RuntimeError(f"claude -p returned error payload: {str(data)[:300]}")
    text = str(data.get("result") or "")
    u = data.get("usage") or {}
    raw_meta = {
        "cli": "claude",
        "upstream_model": "claude-sonnet-5-5",
        "cli_input_tokens": int(u.get("input_tokens") or 0),
        "cli_cache_creation_tokens": int(u.get("cache_creation_input_tokens") or 0),
        "cli_cache_read_tokens": int(u.get("cache_read_input_tokens") or 0),
        "cli_output_tokens": int(u.get("output_tokens") or 0),
        "cli_thinking_tokens": int(
            ((u.get("output_tokens_details") or {}).get("thinking_tokens")) or 0
        ),
        "cli_reported_cost_usd": float(data.get("total_cost_usd") or 0.0),
        "num_turns": int(data.get("num_turns") or 1),
    }
    return text, raw_meta


def _invoke_codex_wsl(system_text: str, user_text: str, timeout_s: float) -> tuple[str, dict[str, Any]]:
    codex_bin = shutil.which("codex") or "/usr/bin/codex"
    combined = _format_combined_prompt(system_text, user_text)
    with tempfile.TemporaryDirectory(prefix="aria_wsl_codex_") as td:
        cmd = [
            codex_bin,
            "exec",
            combined,
            "-m",
            "gpt-5.6-sol",
            "--ephemeral",
            "--skip-git-repo-check",
            "--ignore-user-config",
            "--ignore-rules",
            "-s",
            "read-only",
            "-C",
            td,
            "-c",
            'model_reasoning_effort="low"',
            "--json",
        ]
        proc = subprocess.run(
            cmd,
            cwd=td,
            env=_clean_subprocess_env(),
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
    if proc.returncode != 0:
        raise RuntimeError(
            f"codex exec failed (rc={proc.returncode}): {(proc.stderr or proc.stdout)[:300]}"
        )
    last_text = ""
    usage: dict[str, Any] = {}
    turns = 0
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            ev = json.loads(line)
        except json.JSONDecodeError:
            continue
        etype = ev.get("type")
        if etype == "item.completed":
            item = ev.get("item") or {}
            if item.get("type") == "agent_message" and item.get("text"):
                last_text = str(item["text"])
        elif etype == "turn.completed":
            turns += 1
            usage = ev.get("usage") or usage
    if not last_text:
        raise RuntimeError(f"codex exec produced no agent_message: {proc.stdout[:300]}")
    raw_meta = {
        "cli": "codex",
        "upstream_model": "gpt-5.6-sol",
        "cli_input_tokens": int(usage.get("input_tokens") or 0),
        "cli_cache_creation_tokens": int(usage.get("cache_write_input_tokens") or 0),
        "cli_cache_read_tokens": int(usage.get("cached_input_tokens") or 0),
        "cli_output_tokens": int(usage.get("output_tokens") or 0),
        "cli_thinking_tokens": int(usage.get("reasoning_output_tokens") or 0),
        "cli_reported_cost_usd": None,
        "num_turns": turns or 1,
    }
    return last_text, raw_meta


def _invoke_agy_wsl(system_text: str, user_text: str, timeout_s: float) -> tuple[str, dict[str, Any]]:
    agy_bin = shutil.which("agy") or str(Path.home() / ".local" / "bin" / "agy")
    combined = _format_combined_prompt(system_text, user_text)
    with tempfile.TemporaryDirectory(prefix="aria_wsl_agy_") as td:
        cmd = [
            agy_bin,
            "-p",
            combined,
            "--model",
            "gemini-4-carbon-low",
            "--disable-slash-commands",
            "--output-format",
            "json",
        ]
        proc = subprocess.run(
            cmd,
            cwd=td,
            env=_clean_subprocess_env(),
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
    if proc.returncode != 0:
        raise RuntimeError(
            f"agy -p failed (rc={proc.returncode}): {(proc.stderr or proc.stdout)[:300]}"
        )
    # Locate the JSON object in stdout (skipping any spinner carriage returns if present)
    raw_out = proc.stdout.strip()
    idx = raw_out.find("{")
    if idx < 0:
        raise RuntimeError(f"agy -p returned non-JSON stdout: {raw_out[:300]}")
    data = json.loads(raw_out[idx:])
    if data.get("status") != "SUCCESS":
        raise RuntimeError(f"agy -p non-SUCCESS status: {str(data)[:300]}")
    text = str(data.get("response") or "")
    u = data.get("usage") or {}
    raw_meta = {
        "cli": "agy",
        "upstream_model": "gemini-4-carbon-low",
        "cli_input_tokens": int(u.get("input_tokens") or 0),
        "cli_cache_creation_tokens": 0,
        "cli_cache_read_tokens": int(u.get("cache_read_tokens") or 0),
        "cli_output_tokens": int(u.get("output_tokens") or 0),
        "cli_thinking_tokens": int(u.get("thinking_tokens") or 0),
        "cli_reported_cost_usd": None,
        "num_turns": int(data.get("num_turns") or 1),
    }
    return text, raw_meta


def _dispatch_cli(model_id: str, system_text: str, user_text: str, timeout_s: float) -> tuple[str, dict[str, Any]]:
    last_err: Exception | None = None
    for attempt in range(2):
        try:
            if model_id == "wsl-claude-sonnet-5-5":
                return _invoke_claude_wsl(system_text, user_text, timeout_s)
            if model_id == "wsl-codex-gpt-5-6-sol":
                return _invoke_codex_wsl(system_text, user_text, timeout_s)
            if model_id == "wsl-agy-gemini-4-carbon-low":
                return _invoke_agy_wsl(system_text, user_text, timeout_s)
            raise ValueError(f"Unsupported shim model: {model_id!r}")
        except Exception as exc:
            last_err = exc
            if attempt == 0:
                time.sleep(2.0)
    assert last_err is not None
    raise last_err


class ShimHandler(BaseHTTPRequestHandler):
    ledger_path: Path = Path.home() / ".aria" / "idea2_runs" / "shim_telemetry.jsonl"
    call_timeout_s: float = 240.0

    def log_message(self, fmt: str, *args: Any) -> None:
        return

    def _send_json(self, code: int, payload: dict[str, Any]) -> None:
        raw = json.dumps(payload).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def do_GET(self) -> None:
        if self.path.startswith("/health") or self.path.startswith("/v1/models"):
            self._send_json(
                200,
                {
                    "status": "ok",
                    "shim": "wsl_cli_openai_shim",
                    "data": [{"id": m, "object": "model"} for m in SUPPORTED_MODELS],
                },
            )
            return
        self._send_json(404, {"error": "not_found"})

    def do_POST(self) -> None:
        if not self.path.startswith("/v1/chat/completions"):
            self._send_json(404, {"error": "not_found"})
            return
        try:
            n_bytes = int(self.headers.get("Content-Length") or 0)
            body = json.loads(self.rfile.read(n_bytes).decode("utf-8"))
            raw_model = str(body.get("model") or "").strip()
            model_id = raw_model.split("/")[-1]
            messages = body.get("messages") or []
            sys_text, usr_text = _extract_prompts(messages)

            t0 = time.monotonic()
            text, raw_meta = _dispatch_cli(
                model_id=model_id,
                system_text=sys_text,
                user_text=usr_text,
                timeout_s=self.call_timeout_s,
            )
            elapsed_s = round(time.monotonic() - t0, 3)

            net_in = _count_tokens(sys_text) + _count_tokens(usr_text) + 6
            net_out = _count_tokens(text)
            in_rate, out_rate = PUBLISHED_TARIFF_PER_MTOK.get(model_id, (0.0, 0.0))
            projected_cost_usd = round(
                (net_in * in_rate + net_out * out_rate) / 1_000_000.0, 6
            )

            call_id = "chatcmpl-" + hashlib.sha256(
                f"{time.time_ns()}:{model_id}:{net_in}:{net_out}".encode("utf-8")
            ).hexdigest()[:16]

            ledger_row = {
                "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                "call_id": call_id,
                "shim_model": model_id,
                "elapsed_seconds": elapsed_s,
                "net_aria_prompt_tokens": net_in,
                "net_aria_completion_tokens": net_out,
                "net_aria_total_tokens": net_in + net_out,
                "projected_api_list_cost_usd": projected_cost_usd,
                "tariff_per_mtok_usd": {"input": in_rate, "output": out_rate},
                **raw_meta,
            }
            self.ledger_path.parent.mkdir(parents=True, exist_ok=True)
            with _LEDGER_LOCK:
                with open(self.ledger_path, "a", encoding="utf-8") as lf:
                    lf.write(json.dumps(ledger_row) + "\n")

            self._send_json(
                200,
                {
                    "id": call_id,
                    "object": "chat.completion",
                    "created": int(time.time()),
                    "model": model_id,
                    "choices": [
                        {
                            "index": 0,
                            "message": {"role": "assistant", "content": text},
                            "finish_reason": "stop",
                        }
                    ],
                    "usage": {
                        "prompt_tokens": net_in,
                        "completion_tokens": net_out,
                        "total_tokens": net_in + net_out,
                    },
                },
            )
        except Exception as exc:
            self._send_json(
                500,
                {
                    "error": {
                        "message": f"{type(exc).__name__}: {str(exc)[:400]}",
                        "type": "shim_execution_error",
                    }
                },
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Loopback OpenAI shim for local WSL CLIs.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=11436)
    parser.add_argument(
        "--ledger-path",
        type=Path,
        default=Path.home() / ".aria" / "idea2_runs" / "shim_telemetry.jsonl",
    )
    parser.add_argument("--timeout", type=float, default=240.0)
    args = parser.parse_args()

    if args.host != "127.0.0.1":
        raise ValueError("wsl_cli_openai_shim must bind strictly to 127.0.0.1")

    ShimHandler.ledger_path = args.ledger_path
    ShimHandler.call_timeout_s = args.timeout
    server = ThreadingHTTPServer((args.host, args.port), ShimHandler)
    print(
        f"[wsl_cli_openai_shim] Listening on http://{args.host}:{args.port} "
        f"(models={list(SUPPORTED_MODELS)})",
        flush=True,
    )
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
