"""Stage 4 R6: bus durability + per-run isolation + indexed reads.

- A persisted bus survives a process crash and can be replayed.
- get_findings / get_pending_checkpoints serve from indices and can be scoped to
  one experiment so concurrent runs sharing the global bus do not cross-read.
- The eviction-aware indices stay consistent when the deque overflows.
"""

from aria.bus.message_bus import (
    MessageBus, Message, MessageType, Confidence,
)


def _finding(exp, sender="a"):
    return Message(sender=sender, type=MessageType.FINDING,
                   confidence=Confidence.HIGH, experiment_id=exp,
                   payload={"k": "v"})


def _escalation(exp, cp=1):
    return Message(sender="a", type=MessageType.ESCALATION,
                   confidence=Confidence.MEDIUM, experiment_id=exp,
                   checkpoint=cp, payload={"checkpoint": cp})


# ── indexed + scoped reads ───────────────────────────────────────────────────

def test_get_findings_scoped_to_experiment():
    bus = MessageBus()
    bus.publish(_finding("exp-A"))
    bus.publish(_finding("exp-A"))
    bus.publish(_finding("exp-B"))
    assert len(bus.get_findings("exp-A")) == 2
    assert len(bus.get_findings("exp-B")) == 1
    assert bus.get_findings("exp-missing") == []


def test_pending_checkpoints_isolated_between_runs():
    bus = MessageBus()
    a = _escalation("exp-A")
    b = _escalation("exp-B")
    bus.publish(a)
    bus.publish(b)
    # global view sees both; scoped view sees only its own run
    assert len(bus.get_pending_checkpoints()) == 2
    scoped = bus.get_pending_checkpoints(experiment_id="exp-A")
    assert [m.id for m in scoped] == [a.id]
    # resolving one removes it from pending (index reflects payload mutation)
    bus.resolve_checkpoint(a.id, {"choice": "ok"})
    assert bus.get_pending_checkpoints(experiment_id="exp-A") == []
    assert len(bus.get_pending_checkpoints(experiment_id="exp-B")) == 1


def test_indices_stay_consistent_under_eviction():
    bus = MessageBus(max_log_size=3)
    f1 = _finding("exp-A")
    bus.publish(f1)
    bus.publish(_escalation("exp-A"))
    bus.publish(_finding("exp-A"))
    # next publish overflows the deque (maxlen=3) and evicts f1
    bus.publish(_finding("exp-A"))
    findings = bus.get_findings("exp-A")
    assert f1 not in findings              # evicted finding is de-indexed
    assert len(findings) == 2              # the two newest findings remain


# ── durability ───────────────────────────────────────────────────────────────

def test_persisted_bus_survives_and_replays(tmp_path):
    log_path = tmp_path / "run1" / "bus_log.jsonl"
    bus = MessageBus(persist_path=str(log_path))
    bus.publish(_finding("exp-A", sender="qc"))
    bus.publish(_escalation("exp-A", cp=2))
    bus.publish(_finding("exp-B", sender="de"))

    assert log_path.exists()
    # Simulate a crash: a brand-new process has no in-memory log, only the file.
    replayed = MessageBus.replay(str(log_path))
    assert len(replayed) == 3
    assert {r["experiment_id"] for r in replayed} == {"exp-A", "exp-B"}
    # scoped replay for postmortem of a single run
    only_a = MessageBus.replay(str(log_path), experiment_id="exp-A")
    assert len(only_a) == 2
    assert all(r["experiment_id"] == "exp-A" for r in only_a)


def test_replay_missing_file_is_empty():
    assert MessageBus.replay("/nonexistent/path/bus_log.jsonl") == []


def test_enable_persistence_after_construction(tmp_path):
    log_path = tmp_path / "late" / "bus_log.jsonl"
    bus = MessageBus()
    bus.publish(_finding("exp-A"))          # before persistence — not written
    bus.enable_persistence(str(log_path))
    bus.publish(_finding("exp-A"))          # after — written
    replayed = MessageBus.replay(str(log_path))
    assert len(replayed) == 1
