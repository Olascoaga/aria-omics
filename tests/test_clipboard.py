"""OS clipboard reader for the Textual front door (offline, mocked backends).

Textual's built-in Ctrl+V pastes from the empty in-app clipboard; ARIA reads the
real OS clipboard instead. These tests pin backend selection (PATH + absolute
WSL fallbacks) and normalization without touching a real clipboard.
"""

from __future__ import annotations

from aria.ui import clipboard


class _Proc:
    def __init__(self, returncode: int, stdout: str):
        self.returncode = returncode
        self.stdout = stdout


def test_read_clipboard_none_when_no_backend(monkeypatch):
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: None)
    # No absolute fallback exists either.
    monkeypatch.setattr(clipboard.os.path, "isfile", lambda path: False)
    assert clipboard.clipboard_backend() is None
    assert clipboard.read_clipboard() is None


def test_read_clipboard_normalizes_crlf_and_trailing_newline(monkeypatch):
    # Only powershell.exe is "available"; it returns CRLF + trailing newline.
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


def test_read_clipboard_uses_absolute_powershell_fallback(monkeypatch):
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


def test_read_clipboard_skips_failing_backend(monkeypatch):
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: "/fake/" + name)
    monkeypatch.setattr(clipboard.os.path, "isfile", lambda path: False)

    def fake_run(argv, **kwargs):
        # First backend (powershell) fails; the next that yields output wins.
        if argv[0].endswith("powershell.exe"):
            return _Proc(1, "")
        return _Proc(0, "pasted\n")

    monkeypatch.setattr(clipboard.subprocess, "run", fake_run)
    assert clipboard.read_clipboard() == "pasted"


def test_read_clipboard_swallows_backend_errors(monkeypatch):
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: "/fake/" + name)
    monkeypatch.setattr(clipboard.os.path, "isfile", lambda path: False)

    def boom(*a, **k):
        raise OSError("backend exploded")

    monkeypatch.setattr(clipboard.subprocess, "run", boom)
    assert clipboard.read_clipboard() is None
