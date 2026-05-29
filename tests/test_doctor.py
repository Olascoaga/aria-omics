from pathlib import Path
import re


def test_doctor_smoke_passes_without_benchmark_data(monkeypatch, tmp_path):
    from aria.doctor import run_doctor

    monkeypatch.setenv("HOME", str(tmp_path))

    code, messages = run_doctor("smoke")

    assert code == 0
    assert any("Result: passed" in message for message in messages)


def test_aria_versions_share_single_source():
    import aria
    import aria.llm
    from aria.version import __version__

    assert aria.__version__ == __version__
    assert aria.llm.__version__ == __version__


def test_installer_does_not_write_api_keys_to_bashrc():
    install_text = Path("install.sh").read_text(encoding="utf-8")

    assert re.search(r"ANTHROPIC_API_KEY.*\.bashrc", install_text) is None
    assert re.search(r"GEMINI_API_KEY.*\.bashrc", install_text) is None
    assert re.search(r"GOOGLE_API_KEY.*\.bashrc", install_text) is None
    assert ">> \"$HOME/.bashrc\"" in install_text
