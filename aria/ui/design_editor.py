"""U2 — tabular design editor (Textual modal).

Opened from the cockpit when the design-groups checkpoint (CP2.1) is pending. It
edits the proposed ``{group: [samples]}`` mapping in a table and, on submit,
emits it as JSON through the caller's ``on_submit`` — which the cockpit wires to
the SAME ``orchestrator.on_checkpoint_resolved`` the manual path already uses
(``DesignAgent._parse_manual_groups`` accepts JSON). The editor is a view; all
state logic lives in the Textual-free :class:`aria.ui.design_model.DesignDraft`.
"""

from __future__ import annotations

from typing import Callable

from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import DataTable, Footer, Input, Static

from aria.ui.design_model import DesignDraft


class DesignEditorScreen(ModalScreen):
    """Modal table editor for experimental groups (CP2.1)."""

    CSS = """
    DesignEditorScreen { align: center middle; }
    #design-editor { width: 80%; height: 80%; border: round $accent; padding: 1 2; }
    #design-table { height: 1fr; }
    #design-input { display: none; }
    """

    BINDINGS = [
        Binding("c", "cycle", "Cycle group"),
        Binding("n", "new_group", "New group"),
        Binding("s", "submit", "Submit"),
        Binding("escape", "cancel", "Cancel"),
    ]

    def __init__(self, proposed_groups: dict,
                 on_submit: Callable[[str], None]):
        super().__init__()
        self.draft = DesignDraft.from_proposed(proposed_groups)
        self._on_submit = on_submit

    def compose(self):
        with Vertical(id="design-editor"):
            yield Static(
                "Tabular design editor — [c] cycle group · [n] new group · "
                "[s] submit · [esc] cancel",
                id="design-help")
            yield DataTable(id="design-table")
            yield Static("", id="design-status")
            yield Input(placeholder="new group name", id="design-input")
            yield Footer()

    def on_mount(self) -> None:
        table = self.query_one("#design-table", DataTable)
        table.cursor_type = "row"
        table.add_column("Sample", key="sample")
        table.add_column("Group", key="group")
        for s in self.draft.samples:
            table.add_row(s, self.draft.group_of(s), key=s)
        table.focus()
        self._refresh_status()

    # ── helpers ──────────────────────────────────────────────────────────────
    def _current_sample(self):
        if not self.draft.samples:
            return None
        row = self.query_one("#design-table", DataTable).cursor_row
        if row is None or not (0 <= row < len(self.draft.samples)):
            return None
        return self.draft.samples[row]

    def _refresh_status(self) -> None:
        groups = self.draft.to_groups()
        parts = ", ".join(f"{g}:{len(m)}" for g, m in groups.items()) \
            or "no groups"
        self.query_one("#design-status", Static).update(f"Groups → {parts}")

    # ── actions ──────────────────────────────────────────────────────────────
    def action_cycle(self) -> None:
        s = self._current_sample()
        if s is None:
            return
        new = self.draft.cycle(s)
        self.query_one("#design-table", DataTable).update_cell(s, "group", new)
        self._refresh_status()

    def action_new_group(self) -> None:
        inp = self.query_one("#design-input", Input)
        inp.display = True
        inp.focus()

    def on_input_submitted(self, event: Input.Submitted) -> None:
        name = (event.value or "").strip()
        inp = self.query_one("#design-input", Input)
        inp.value = ""
        inp.display = False
        if name:
            s = self._current_sample()
            if s is not None:
                self.draft.assign(s, name)
                self.query_one("#design-table", DataTable).update_cell(
                    s, "group", name)
        self.query_one("#design-table", DataTable).focus()
        self._refresh_status()

    def action_submit(self) -> None:
        ok, reason = self.draft.is_valid()
        if not ok:
            self.query_one("#design-status", Static).update(
                f"[red]Cannot submit: {reason}[/]")
            return
        payload = self.draft.to_json()
        try:
            self._on_submit(payload)
        except Exception:
            pass
        self.dismiss(payload)

    def action_cancel(self) -> None:
        self.dismiss(None)
