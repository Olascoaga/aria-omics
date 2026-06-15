"""U2 DesignDraft tests (pure, no Textual).

Covers the draft logic and the load-bearing contract: the editor's JSON output
must be parsed back by DesignAgent._parse_manual_groups (the existing manual
assignment path) so the editor needs no new design code path.
"""

from __future__ import annotations

import json

from aria.ui.design_model import DesignDraft


def test_from_proposed_and_to_groups_roundtrip():
    proposed = {"treated": ["s1", "s2"], "control": ["s3", "s4"]}
    d = DesignDraft.from_proposed(proposed)
    assert d.samples == ["s1", "s2", "s3", "s4"]
    assert d.groups == ["treated", "control"]
    assert d.group_of("s1") == "treated"
    assert d.to_groups() == proposed


def test_cycle_round_robin():
    d = DesignDraft.from_proposed({"a": ["s1"], "b": ["s2"]})
    assert d.cycle("s1") == "b"     # a -> b
    assert d.cycle("s1") == "a"     # b -> a (wrap)
    assert d.group_of("s1") == "a"


def test_assign_registers_new_group_and_drops_empty():
    d = DesignDraft.from_proposed({"a": ["s1"], "b": ["s2"]})
    d.assign("s2", "a")             # b now empty -> dropped
    groups = d.to_groups()
    assert groups == {"a": ["s1", "s2"]}
    # new group via assign
    d.assign("s1", "c")
    assert "c" in d.groups
    assert d.to_groups()["c"] == ["s1"]


def test_is_valid_requires_two_nonempty_groups():
    d = DesignDraft.from_proposed({"a": ["s1"], "b": ["s2"]})
    assert d.is_valid()[0] is True
    d.assign("s2", "a")             # collapse to one group
    ok, reason = d.is_valid()
    assert ok is False
    assert "2" in reason


def test_to_json_is_valid_json():
    d = DesignDraft.from_proposed({"treated": ["s1"], "control": ["s2"]})
    obj = json.loads(d.to_json())
    assert obj == {"treated": ["s1"], "control": ["s2"]}


def test_json_roundtrips_through_design_agent_manual_path():
    # The editor submits JSON; DesignAgent._parse_manual_groups must accept it.
    from aria.agents.design_agent import DesignAgent

    d = DesignDraft.from_proposed({"treated": ["s1", "s2"], "control": ["s3"]})
    d.cycle("s2")                   # move s2 treated -> control

    agent = object.__new__(DesignAgent)   # bypass __init__/memory; only need state
    agent._parsed_samples = [
        {"stem": s} for s in ("s1", "s2", "s3")
    ]
    parsed = agent._parse_manual_groups(d.to_json())
    # to_groups preserves sample order (s2 before s3), so control == [s2, s3].
    assert parsed == {"treated": ["s1"], "control": ["s2", "s3"]}
