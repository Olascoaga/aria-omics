"""Read the OS clipboard for the Textual front door.

Textual's built-in Ctrl+V (``Input``/``TextArea`` ``action_paste``) pastes from
the app's OWN in-process clipboard (``App.clipboard``), which is empty unless the
user copied something inside the running app. Terminals that pass Ctrl+V through
to the application (instead of performing their own bracketed paste) therefore
paste nothing — the symptom Samael hit in Windows Terminal and MobaXterm, while
the classic cooked-mode Rich TUI pasted fine.

This reads the REAL OS clipboard so the intake can paste external text (data
paths, biological questions). It is:

- local-only — no network egress, so it is unaffected by air-gapped mode;
- best-effort — tries several backends and returns ``None`` if none work;
- non-raising — a missing/failing backend never crashes the UI.
"""

from __future__ import annotations

import os
import shutil
import subprocess

# Clipboard read commands, ordered by likelihood on ARIA's target environments:
# WSL2 (Windows clipboard via interop) first, then Wayland, X11, and macOS.
_BACKENDS: list[list[str]] = [
    ["powershell.exe", "-NoProfile", "-Command", "Get-Clipboard"],
    ["pwsh.exe", "-NoProfile", "-Command", "Get-Clipboard"],
    ["wl-paste", "--no-newline"],
    ["xclip", "-selection", "clipboard", "-o"],
    ["xsel", "-b", "-o"],
    ["pbpaste"],
]

# Absolute fallbacks for tools that WSL interop can execute even when the Windows
# paths are not on a (login-stripped) PATH. Only used when `shutil.which` misses.
_FALLBACK_PATHS: dict[str, list[str]] = {
    "powershell.exe": [
        "/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe",
    ],
    "pwsh.exe": [
        "/mnt/c/Program Files/PowerShell/7/pwsh.exe",
    ],
}


def _resolve(cmd: str) -> str | None:
    """Resolve a backend command to a runnable path (PATH first, then fallbacks)."""
    hit = shutil.which(cmd)
    if hit:
        return hit
    for path in _FALLBACK_PATHS.get(cmd, []):
        if os.path.isfile(path):
            return path
    return None


def clipboard_backend() -> str | None:
    """Return the name of the first available clipboard backend, or ``None``."""
    for argv in _BACKENDS:
        if _resolve(argv[0]) is not None:
            return argv[0]
    return None


def read_clipboard(*, timeout: float = 2.0) -> str | None:
    """Return the OS clipboard text, or ``None`` if it cannot be read.

    Trailing newlines and CRLF line endings (notably from
    ``powershell.exe Get-Clipboard``) are normalized so a pasted path or
    question does not carry stray carriage returns or a terminal newline.
    """
    for argv in _BACKENDS:
        resolved = _resolve(argv[0])
        if resolved is None:
            continue
        try:
            proc = subprocess.run(
                [resolved, *argv[1:]], capture_output=True, text=True,
                timeout=timeout,
            )
        except Exception:
            continue
        if proc.returncode != 0 or not proc.stdout:
            continue
        text = proc.stdout.replace("\r\n", "\n").replace("\r", "\n")
        if text.endswith("\n"):
            text = text[:-1]
        if text:
            return text
    return None
