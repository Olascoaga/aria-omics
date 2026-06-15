"""Pure design-draft model for the U2 tabular design editor.

A small, Textual-free state object over a ``{group: [samples]}`` mapping (the
shape ``DesignAgent`` proposes at checkpoint 2.1 and accepts back as JSON via
``_parse_manual_groups``). The Textual editor (:mod:`aria.ui.design_editor`)
drives this; keeping the logic here means it unit-tests in the standard env and
the editor stays a thin view.

The serialized form is a JSON ``{group: [samples]}`` object, which the design
agent's manual-assignment path parses directly — so the editor submits through
the SAME governed checkpoint resolution, never a new design code path.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field


@dataclass
class DesignDraft:
    samples: list[str] = field(default_factory=list)
    assignment: dict[str, str] = field(default_factory=dict)  # sample -> group
    groups: list[str] = field(default_factory=list)           # ordered group names

    @classmethod
    def from_proposed(cls, proposed: dict[str, list[str]]) -> "DesignDraft":
        groups: list[str] = []
        samples: list[str] = []
        assignment: dict[str, str] = {}
        for group, members in (proposed or {}).items():
            g = str(group)
            if g not in groups:
                groups.append(g)
            for s in members or []:
                s = str(s)
                if s not in assignment:
                    samples.append(s)
                assignment[s] = g
        return cls(samples=samples, assignment=assignment, groups=groups)

    def group_of(self, sample: str) -> str:
        return self.assignment.get(sample, "")

    def add_group(self, name: str) -> None:
        name = str(name).strip()
        if name and name not in self.groups:
            self.groups.append(name)

    def assign(self, sample: str, group: str) -> None:
        """Assign ``sample`` to ``group`` (registering the group if new)."""
        if sample not in self.assignment:
            self.samples.append(sample)
        self.add_group(group)
        self.assignment[sample] = str(group).strip()

    def cycle(self, sample: str) -> str:
        """Move ``sample`` to the next known group (round-robin). Returns it."""
        if not self.groups:
            return self.group_of(sample)
        cur = self.assignment.get(sample)
        if cur in self.groups:
            nxt = self.groups[(self.groups.index(cur) + 1) % len(self.groups)]
        else:
            nxt = self.groups[0]
        self.assignment[sample] = nxt
        return nxt

    def to_groups(self) -> dict[str, list[str]]:
        """Render ``{group: [samples]}`` in group order, sample order preserved,
        dropping any group left with no samples."""
        out: dict[str, list[str]] = {g: [] for g in self.groups}
        for s in self.samples:
            g = self.assignment.get(s)
            if g is None:
                continue
            out.setdefault(g, []).append(s)
        return {g: members for g, members in out.items() if members}

    def to_json(self) -> str:
        return json.dumps(self.to_groups())

    def is_valid(self) -> tuple[bool, str]:
        """A comparison needs at least two non-empty groups."""
        groups = self.to_groups()
        if len(groups) < 2:
            return False, "need at least 2 non-empty groups for a comparison"
        return True, ""
