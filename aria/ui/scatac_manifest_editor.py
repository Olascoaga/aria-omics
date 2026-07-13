"""Typed JSON editor for the CP1 scATAC FASTQ library manifest."""

from __future__ import annotations

import json
from typing import Callable

from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Footer, Static, TextArea

from aria.utils.scatac_fastq_manifest import resolve_scatac_fastq_manifest


EMPTY_SCATAC_MANIFEST = {
    "schema_version": "1",
    "libraries": [{
        "library_id": "",
        "sample_id": "",
        "donor_id": "",
        "fastqs": {"R1": "", "R2": "", "R3": ""},
        "barcode_whitelist": "",
        "metadata": {},
    }],
}


def parse_scatac_manifest_edit(text: str, data_dir: str) -> tuple[dict | None, list[str]]:
    """Parse and fully validate a user-edited manifest payload."""
    try:
        raw = json.loads(text)
    except json.JSONDecodeError as exc:
        return None, [f"invalid JSON: {exc}"]
    result = resolve_scatac_fastq_manifest(
        {"data_dir": data_dir, "scatac_fastq_manifest": raw},
        require_paths=True,
    )
    return result.get("manifest"), list(result.get("errors") or [])


class ScatacManifestEditorScreen(ModalScreen):
    """Modal CP1 editor; scientific state remains in the shared manifest validator."""

    CSS = """
    ScatacManifestEditorScreen { align: center middle; }
    #manifest-editor { width: 90%; height: 90%; border: round $accent; padding: 1 2; }
    #manifest-json { height: 1fr; }
    """

    BINDINGS = [
        Binding("ctrl+s", "submit", "Validate + submit"),
        Binding("escape", "cancel", "Cancel"),
    ]

    def __init__(self, manifest: dict | None, data_dir: str,
                 on_submit: Callable[[dict], None]):
        super().__init__()
        self.manifest = manifest or EMPTY_SCATAC_MANIFEST
        self.data_dir = data_dir
        self._on_submit = on_submit

    def compose(self):
        with Vertical(id="manifest-editor"):
            yield Static(
                "10x/scATAC manifest — one row per library; R1/R2/R3, "
                "whitelist, sample and donor are required. Ctrl+S validates."
            )
            yield TextArea(
                json.dumps(self.manifest, indent=2, sort_keys=True),
                id="manifest-json",
            )
            yield Static("", id="manifest-status")
            yield Footer()

    def action_submit(self) -> None:
        text = self.query_one("#manifest-json", TextArea).text
        manifest, errors = parse_scatac_manifest_edit(text, self.data_dir)
        if errors or manifest is None:
            self.query_one("#manifest-status", Static).update(
                "[red]Cannot submit: " + "; ".join(errors) + "[/]"
            )
            return
        try:
            self._on_submit(manifest)
        except Exception:
            return
        self.dismiss(manifest)

    def action_cancel(self) -> None:
        self.dismiss(None)
