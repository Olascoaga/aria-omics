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

Clipboard *source* matters. In a local WSL2/X11 session, text copied in a Linux
app lands in the X11 clipboard, while ``powershell.exe Get-Clipboard`` reads the
Windows clipboard. Reading the wrong one pastes stale, unrelated text.

In a remote SSH session (e.g. MobaXterm connected to a server), ARIA cannot read
the client's local clipboard. OS clipboard backends run on the remote host and
would paste the remote machine's clipboard instead. In that topology Ctrl+V is
disabled by default; use the terminal's paste operation (bracketed paste) or set
``ARIA_CLIPBOARD_SOURCE`` explicitly if reading the remote OS clipboard is
intended.
"""

from __future__ import annotations

import os
import shutil
import subprocess

# Backend command groups by clipboard source.
_WINDOWS: list[list[str]] = [
    ["powershell.exe", "-NoProfile", "-Command", "Get-Clipboard"],
    ["pwsh.exe", "-NoProfile", "-Command", "Get-Clipboard"],
]
_WAYLAND: list[list[str]] = [
    ["wl-paste", "--no-newline"],
]
# X11 has two selections: CLIPBOARD (explicit copy) and PRIMARY (text selection,
# which is how terminal copy-on-select often behaves). Try CLIPBOARD first.
_X11: list[list[str]] = [
    ["xclip", "-selection", "clipboard", "-o"],
    ["xsel", "-b", "-o"],
    ["xclip", "-selection", "primary", "-o"],
    ["xsel", "-p", "-o"],
]
_MAC: list[list[str]] = [
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


# Named clipboard sources for the ARIA_CLIPBOARD_SOURCE override.
_SOURCE_GROUPS: dict[str, list[list[str]]] = {
    "windows": _WINDOWS,
    "x11": _X11,
    "wayland": _WAYLAND,
    "mac": _MAC,
}
_SOURCE_ORDER = ("windows", "x11", "wayland", "mac")
_DISABLED_SOURCES = {"off", "none", "disabled", "terminal"}


def _remote_ssh_session() -> bool:
    """True when ARIA is running on a host reached through SSH.

    In this topology the process-local clipboard backends read the remote host,
    not the user's local terminal/client clipboard.
    """
    return bool(os.environ.get("SSH_CONNECTION") or os.environ.get("SSH_CLIENT"))


def _ordered_backends() -> list[list[str]]:
    """Backend order matching the active GUI session's clipboard.

    A copy made in the user's environment lands in the clipboard of whatever GUI
    session they are in, so prefer that source and keep the others as fallback.

    ``ARIA_CLIPBOARD_SOURCE`` (windows|x11|wayland|mac) forces a primary source,
    keeping the rest as fallback. ``off``/``none``/``terminal`` disables OS
    clipboard reads. ``auto`` / unset keeps the session-derived default.
    """
    forced = os.environ.get("ARIA_CLIPBOARD_SOURCE", "").strip().lower()
    if forced in _DISABLED_SOURCES:
        return []
    if forced and forced != "auto" and forced in _SOURCE_GROUPS:
        rest = [c for name in _SOURCE_ORDER if name != forced
                for c in _SOURCE_GROUPS[name]]
        return _SOURCE_GROUPS[forced] + rest
    if _remote_ssh_session():
        return []
    if os.environ.get("WAYLAND_DISPLAY"):
        return _WAYLAND + _X11 + _WINDOWS + _MAC
    if os.environ.get("DISPLAY"):
        return _X11 + _WINDOWS + _WAYLAND + _MAC
    return _WINDOWS + _WAYLAND + _X11 + _MAC


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
    for argv in _ordered_backends():
        if _resolve(argv[0]) is not None:
            return argv[0]
    return None


def unavailable_message() -> str:
    """Explain why Ctrl+V cannot read an OS clipboard in this session."""
    forced = os.environ.get("ARIA_CLIPBOARD_SOURCE", "").strip().lower()
    if forced in _DISABLED_SOURCES:
        return (
            "OS clipboard reads are disabled. Use the terminal paste operation "
            "(right-click, Shift+Insert, or MobaXterm Paste)."
        )
    if _remote_ssh_session() and not (forced and forced in _SOURCE_GROUPS):
        return (
            "Remote SSH sessions cannot read your local clipboard. Use the "
            "terminal paste operation (right-click, Shift+Insert, or MobaXterm "
            "Paste), or set ARIA_CLIPBOARD_SOURCE explicitly to read the remote "
            "host clipboard."
        )
    return (
        "No clipboard tool found. Install xclip/wl-clipboard, or use the "
        "terminal paste operation (Shift+Insert or right-click)."
    )


def read_clipboard(*, timeout: float = 2.0) -> str | None:
    """Return the OS clipboard text, or ``None`` if it cannot be read.

    Backends are tried in session-appropriate order; the first that yields
    non-empty text wins. Trailing newlines and CRLF line endings (notably from
    ``powershell.exe Get-Clipboard``) are normalized so a pasted path or question
    does not carry stray carriage returns or a terminal newline.
    """
    for argv in _ordered_backends():
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
