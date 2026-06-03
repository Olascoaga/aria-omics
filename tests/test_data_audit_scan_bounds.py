from __future__ import annotations

import sys
import types

litellm_stub = types.ModuleType("litellm")
litellm_stub.completion = lambda *args, **kwargs: None
sys.modules.setdefault("litellm", litellm_stub)

from aria.agents.data_audit_agent import DataAuditAgent, DataAuditScanLimits


def _agent() -> DataAuditAgent:
    agent = DataAuditAgent.__new__(DataAuditAgent)
    agent.publish_status = lambda *args, **kwargs: None
    return agent


def test_scan_directory_respects_max_files(tmp_path):
    for i in range(5):
        (tmp_path / f"sample_{i}_R1.fastq.gz").write_text("reads\n")

    agent = _agent()
    files = agent._scan_directory(
        tmp_path,
        DataAuditScanLimits(max_files=2, max_entries=100, max_depth=2),
    )

    assert len(files) == 2
    assert agent._last_scan_report["truncated"] is True
    assert agent._last_scan_report["truncated_reason"] == "max_files"
    assert "max_files" in agent._scan_report_warnings()[0]


def test_scan_directory_respects_depth_and_skips_symlinks(tmp_path):
    deep = tmp_path / "level1" / "level2"
    deep.mkdir(parents=True)
    (deep / "counts.tsv").write_text("gene\ta\n")
    (tmp_path / "root_counts.tsv").write_text("gene\ta\n")
    (tmp_path / "linked_counts.tsv").symlink_to(tmp_path / "root_counts.tsv")

    agent = _agent()
    files = agent._scan_directory(
        tmp_path,
        DataAuditScanLimits(max_files=20, max_entries=100, max_depth=1),
    )

    assert tmp_path / "root_counts.tsv" in files
    assert deep / "counts.tsv" not in files
    assert tmp_path / "linked_counts.tsv" not in files
    assert agent._last_scan_report["skipped_deep_dirs"] == 1
    assert agent._last_scan_report["skipped_symlinks"] == 1
    warnings = "\n".join(agent._scan_report_warnings())
    assert "symlinked paths" in warnings
    assert "deeper than the configured limit" in warnings


def test_scan_directory_respects_immediate_time_budget(tmp_path):
    (tmp_path / "counts.tsv").write_text("gene\ta\n")

    agent = _agent()
    files = agent._scan_directory(
        tmp_path,
        DataAuditScanLimits(max_files=20, max_entries=100, max_depth=2,
                            max_seconds=0.0),
    )

    assert files == []
    assert agent._last_scan_report["truncated"] is True
    assert agent._last_scan_report["truncated_reason"] == "scan_timeout"


def test_scan_directory_respects_configured_file_size_limit(tmp_path):
    small = tmp_path / "small_counts.tsv"
    large = tmp_path / "large_counts.tsv"
    small.write_text("x")
    large.write_text("abcdef")

    agent = _agent()
    files = agent._scan_directory(
        tmp_path,
        DataAuditScanLimits(max_files=20, max_entries=100, max_depth=2,
                            max_file_size_bytes=1),
    )

    assert small in files
    assert large not in files
    assert agent._last_scan_report["skipped_large_files"] == 1
    assert "size limit" in "\n".join(agent._scan_report_warnings())


def test_scan_limits_can_come_from_context_or_environment(tmp_path, monkeypatch):
    for i in range(3):
        (tmp_path / f"counts_{i}.tsv").write_text("gene\ta\n")

    monkeypatch.setenv("ARIA_DATA_AUDIT_MAX_FILES", "1")
    agent = _agent()
    files = agent._scan_directory(tmp_path)

    assert len(files) == 1
    assert agent._last_scan_report["truncated_reason"] == "max_files"
    assert agent._last_scan_report["limits"]["max_files"] == 1

    files = agent._scan_directory(tmp_path, limits={"max_files": 2})

    assert len(files) == 2
    assert agent._last_scan_report["limits"]["max_files"] == 2
