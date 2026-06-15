"""OS clipboard reader for the Textual front door (offline, mocked backends).

Textual's built-in Ctrl+V pastes from the empty in-app clipboard; ARIA reads the
real OS clipboard instead. These tests pin the backend selection + normalization
without touching a real clipboard.
"""

from __future__ import annotations

from aria.ui import clipboard


class _Proc:
    def __init__(self, returncode: int, stdout: str):
        self.returncode = returncode
        self.stdout = stdout


def test_read_clipboard_none_when_no_backend(monkeypatch):
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: None)
    assert clipboard.clipboard_backend() is None
    assert clipboard.read_clipboard() is None


def test_read_clipboard_normalizes_crlf_and_trailing_newline(monkeypatch):
    # Only powershell.exe is "available"; it returns CRLF + trailing newline.
    monkeypatch.setattr(
        clipboard.shutil, "which",
        lambda name: "/fake/ps" if name == "powershell.exe" else None,
    )
    monkeypatch.setattr(
        clipboard.subprocess, "run",
        lambda *a, **k: _Proc(0, "C:\\data\\exp\r\nsecond\r\n"),
    )
    assert clipboard.clipboard_backend() == "powershell.exe"
    assert clipboard.read_clipboard() == "C:\\data\\exp\nsecond"


def test_read_clipboard_skips_failing_backend(monkeypatch):
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: "/fake/" + name)

    def fake_run(argv, **kwargs):
        # First backend (powershell) fails; the next that yields output wins.
        if argv[0] == "powershell.exe":
            return _Proc(1, "")
        return _Proc(0, "pasted\n")

    monkeypatch.setattr(clipboard.subprocess, "run", fake_run)
    assert clipboard.read_clipboard() == "pasted"


def test_read_clipboard_swallows_backend_errors(monkeypatch):
    monkeypatch.setattr(clipboard.shutil, "which", lambda name: "/fake/" + name)

    def boom(*a, **k):
        raise OSError("backend exploded")

    monkeypatch.setattr(clipboard.subprocess, "run", boom)
    assert clipboard.read_clipboard() is None
