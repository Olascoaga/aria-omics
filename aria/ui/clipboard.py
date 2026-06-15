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

import shutil
import subprocess

# Clipboard read commands, ordered by likelihood on ARIA's target environments:
# WSL2 (Windows clipboard via interop) first, then Wayland, X11, and macOS.
_BACKENDS: list[list[str]] = [
    ["powershell.exe", "-NoProfile", "-Command", "Get-Clipboard"],
    ["wl-paste", "--no-newline"],
    ["xclip", "-selection", "clipboard", "-o"],
    ["xsel", "-b", "-o"],
    ["pbpaste"],
]


def clipboard_backend() -> str | None:
    """Return the name of the first available clipboard backend, or ``None``."""
    for argv in _BACKENDS:
        if shutil.which(argv[0]) is not None:
            return argv[0]
    return None


def read_clipboard(*, timeout: float = 2.0) -> str | None:
    """Return the OS clipboard text, or ``None`` if it cannot be read.

    Trailing newlines and CRLF line endings (notably from
    ``powershell.exe Get-Clipboard``) are normalized so a pasted path or
    question does not carry stray carriage returns or a terminal newline.
    """
    for argv in _BACKENDS:
        if shutil.which(argv[0]) is None:
            continue
        try:
            proc = subprocess.run(
                argv, capture_output=True, text=True, timeout=timeout
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
