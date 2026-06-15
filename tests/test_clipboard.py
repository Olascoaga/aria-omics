"""OS clipboard reader for the Textual front door (offline, mocked backends).

Textual's built-in Ctrl+V pastes from the empty in-app clipboard; ARIA reads the
real OS clipboard instead. These tests pin backend selection (session-aware
ordering, PATH + absolute WSL fallbacks) and normalization without touching a
real clipboard.
"""

from __future__ import annotations

import pytest

from aria.ui import clipboard


class _Proc:
    def __init__(self, returncode: int, stdout: str):
        self.returncode = returncode
        self.stdout = stdout


@pytest.fixture
def no_display(monkeypatch):
    """Force the Windows-first order (headless WSL): no X11/Wayland session."""
    monkeypatch.delenv("DISPLAY", raising=False)
    monkeypatch.delenv("WAYLAND_DISPLAY", raising=False)


def test_read_clipboard_none_when_no_backend(monkeypatch, no_display):
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: None)
    monkeypatch.setattr(clipboard.os.path, "isfile", lambda path: False)
    assert clipboard.clipboard_backend() is None
    assert clipboard.read_clipboard() is None


def test_read_clipboard_normalizes_crlf_and_trailing_newline(monkeypatch, no_display):
    monkeypatch.setattr(
        clipboard.shutil, "which",
        lambda name: "/fake/ps" if name == "powershell.exe" else None,
    )
    monkeypatch.setattr(clipboard.os.path, "isfile", lambda path: False)
    monkeypatch.setattr(
        clipboard.subprocess, "run",
        lambda *a, **k: _Proc(0, "C:\\data\\exp\r\nsecond\r\n"),
    )
    assert clipboard.clipboard_backend() == "powershell.exe"
    assert clipboard.read_clipboard() == "C:\\data\\exp\nsecond"


def test_read_clipboard_uses_absolute_powershell_fallback(monkeypatch, no_display):
    # powershell.exe is NOT on PATH (login-stripped), but the WSL absolute path
    # exists and interop can run it.
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: None)
    ps_path = "/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe"
    monkeypatch.setattr(clipboard.os.path, "isfile", lambda path: path == ps_path)
    seen = {}

    def fake_run(argv, **kwargs):
        seen["argv0"] = argv[0]
        return _Proc(0, "/data/from/windows\r\n")

    monkeypatch.setattr(clipboard.subprocess, "run", fake_run)
    assert clipboard.clipboard_backend() == "powershell.exe"
    assert clipboard.read_clipboard() == "/data/from/windows"
    assert seen["argv0"] == ps_path  # ran the resolved absolute path


def test_x11_clipboard_preferred_over_windows_when_display_set(monkeypatch):
    # MobaXterm/X11 session: a Linux copy lands in X11, not the Windows clipboard.
    monkeypatch.setenv("DISPLAY", "localhost:11.0")
    monkeypatch.delenv("WAYLAND_DISPLAY", raising=False)
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: "/usr/bin/" + name)
    monkeypatch.setattr(clipboard.os.path, "isfile", lambda path: False)

    def fake_run(argv, **kwargs):
        if "xclip" in argv[0]:
            return _Proc(0, "from-x11\n")
        return _Proc(0, "stale-windows\n")

    monkeypatch.setattr(clipboard.subprocess, "run", fake_run)
    assert clipboard.clipboard_backend() == "xclip"
    assert clipboard.read_clipboard() == "from-x11"


def test_read_clipboard_skips_failing_backend(monkeypatch, no_display):
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: "/fake/" + name)
    monkeypatch.setattr(clipboard.os.path, "isfile", lambda path: False)

    def fake_run(argv, **kwargs):
        # First backend (powershell, Windows-first order) fails; next wins.
        if argv[0].endswith("powershell.exe"):
            return _Proc(1, "")
        return _Proc(0, "pasted\n")

    monkeypatch.setattr(clipboard.subprocess, "run", fake_run)
    assert clipboard.read_clipboard() == "pasted"


def test_read_clipboard_swallows_backend_errors(monkeypatch, no_display):
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: "/fake/" + name)
    monkeypatch.setattr(clipboard.os.path, "isfile", lambda path: False)

    def boom(*a, **k):
        raise OSError("backend exploded")

    monkeypatch.setattr(clipboard.subprocess, "run", boom)
    assert clipboard.read_clipboard() is None
