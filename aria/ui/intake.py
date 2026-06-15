"""Textual front door for ARIA analysis intake.

This is intentionally thin: it collects the same data/accession + biological
question that the classic Rich intake asks for, then hands the existing context
back to ``aria.tui``. It owns no analysis logic.
"""

from __future__ import annotations

from dataclasses import dataclass

from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.events import Paste
from textual.widgets import Button, Footer, Input, Static, TextArea

from aria.ui.brand import ARIA_BANNER, TAGLINE


@dataclass(frozen=True)
class IntakeResult:
    """User-provided startup fields for a new ARIA analysis."""

    data_input: str
    question: str


class AriaIntakeApp(App):
    """Minimal Textual intake screen shown before a run starts."""

    CSS = """
    #body { height: 1fr; }
    #left { width: 48; padding: 1 2; border-right: solid $accent; }
    #main { width: 1fr; padding: 1 3; }
    #brand-logo { text-style: bold; color: $accent; margin-bottom: 1; }
    #brand-tagline { color: $text-muted; margin-bottom: 1; }
    #app-version { color: $text-muted; margin-bottom: 1; }
    #form-title { text-style: bold; margin-bottom: 1; }
    #memory { height: 1fr; color: $text-muted; }
    #experiments { height: auto; color: $text-muted; }
    #data-input { margin-top: 1; margin-bottom: 1; }
    #question-input { height: 10; margin-top: 1; margin-bottom: 1; }
    #status { height: 1; color: $warning; margin-top: 1; }
    #actions { height: auto; margin-top: 1; }
    Button { margin-right: 2; }
    """

    BINDINGS = [
        Binding("ctrl+s", "start", "Start"),
        Binding("escape", "cancel", "Exit"),
    ]

    def __init__(
        self,
        *,
        startup_context: str = "",
        experiments: list[dict] | None = None,
        version: str = "",
    ):
        super().__init__()
        self.startup_context = startup_context
        self.experiments = experiments or []
        self.version = version
        self.submitted: IntakeResult | None = None
        self.status_message = ""

    def compose(self) -> ComposeResult:
        yield Horizontal(
            Vertical(
                Static(ARIA_BANNER, id="brand-logo"),
                Static(TAGLINE, id="brand-tagline"),
                Static(self.version, id="app-version"),
                Static(self._memory_text(), id="memory"),
                Static(self._experiments_text(), id="experiments"),
                id="left",
            ),
            Vertical(
                Static("New analysis", id="form-title"),
                Static("Data directory or GEO/SRA accession"),
                Input(
                    placeholder="/data/my_experiment or GSE183948",
                    id="data-input",
                ),
                Static("Biological question"),
                TextArea(id="question-input"),
                Static("", id="status"),
                Horizontal(
                    Button("Start", id="start", variant="primary"),
                    Button("Exit", id="exit"),
                    id="actions",
                ),
                id="main",
            ),
            id="body",
        )
        yield Footer()

    def on_mount(self) -> None:
        self.query_one("#data-input", Input).focus()

    def on_paste(self, event: Paste) -> None:
        """Route terminal paste into the focused intake field.

        Textual widgets have paste handlers, but terminal/emulator paste events
        can arrive at the app level. Handling it here keeps the startup form
        usable for long paths and multi-line biological questions.
        """
        if self._paste_into_focused(event.text):
            event.prevent_default()
            event.stop()

    def _memory_text(self) -> str:
        text = (self.startup_context or "").strip()
        if not text or "No experiments" in text:
            return "No prior experiment context."
        return text

    def _experiments_text(self) -> str:
        if not self.experiments:
            return ""
        lines = ["", "Recent experiments"]
        for exp in self.experiments[-5:]:
            exp_id = str(exp.get("id", ""))[:8]
            name = str(exp.get("name", "untitled"))[:28]
            organism = exp.get("organism") or "?"
            updated = str(exp.get("updated_at", ""))[:10]
            lines.append(f"{exp_id}  {name}  {organism}  {updated}")
        return "\n".join(lines)

    def _set_status(self, message: str) -> None:
        self.status_message = message
        self.query_one("#status", Static).update(message)

    def _paste_into_focused(self, text: str) -> bool:
        focused = self.focused
        if isinstance(focused, Input):
            focused.insert_text_at_cursor(text.strip("\r\n"))
            return True
        if isinstance(focused, TextArea):
            focused.insert(text)
            return True
        return False

    def _values(self) -> tuple[str, str]:
        data_input = self.query_one("#data-input", Input).value.strip()
        question = self.query_one("#question-input", TextArea).text.strip()
        return data_input, question

    def action_start(self) -> None:
        data_input, question = self._values()
        if not data_input:
            self._set_status("Enter a data directory or accession.")
            self.query_one("#data-input", Input).focus()
            return
        if not question:
            self._set_status("Enter a biological question.")
            self.query_one("#question-input", TextArea).focus()
            return
        result = IntakeResult(data_input=data_input, question=question)
        self.submitted = result
        self.exit(result)

    def action_cancel(self) -> None:
        self.submitted = None
        self.exit(None)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "start":
            self.action_start()
        elif event.button.id == "exit":
            self.action_cancel()


def launch_intake(
    *,
    startup_context: str = "",
    experiments: list[dict] | None = None,
    version: str = "",
) -> IntakeResult | None:
    """Run the Textual intake and return submitted fields, or ``None``."""

    app = AriaIntakeApp(
        startup_context=startup_context,
        experiments=experiments,
        version=version,
    )
    return app.run()
