# ARIA Control Center

The Control Center is ARIA's Textual terminal interface. It presents the same
orchestrator state and resolves the same checkpoints as the classic Rich flow;
it does not own analysis logic or scientific state.

## Start it

Install ARIA with the optional UI dependency in the active core environment:

```bash
conda activate aria-env
python -m pip install -e '.[tui]'
aria
```

The Control Center opens only when Textual is installed and both stdin and
stdout are attached to a terminal. Verify the condition directly:

```bash
python -c 'from aria.ui.cockpit import cockpit_available; print(cockpit_available())'
```

ARIA uses the classic interface when any of these conditions applies:

- Textual is not installed;
- input or output is not a TTY;
- `ARIA_NO_TUI` is set;
- `aria --classic-tui` is used;
- `aria --reproducible` is used.

## Intake screen

The first screen collects:

- a local data directory or a GEO/SRA accession;
- the biological question;
- an explicit opt-in for speculative hypotheses, off by default.

Use `Ctrl+S` to start or `Esc` to exit. Local paths and accessions are validated
before the screen closes. The left panel shows recent experiment context and
resolved report bundles when they exist; it does not claim that an ephemeral
in-memory run can be resumed after process loss.

The intake and run views live inside one Textual application. Data audit begins
on a worker thread so the interface remains responsive during the transition.

## Run views and keys

The footer shows the active bindings. The center panel switches between these
read-only projections of the run state:

| Key | View | Content |
|---|---|---|
| `o` | Overview | Current phase, progress, checkpoint, and report state |
| `g` | Agents | Agent execution progress |
| `d` | Decisions | Active checkpoint and available choices |
| `f` | Findings | Evidence-tiered findings, including older-item counts |
| `a` | Artifacts | Report, methods, figures, tables, and compiled claim linkage |
| `u` | Resources | Runtime resource information |
| `l` | Ledger | Planned versus executed analyses |
| `r` | Readiness | Registry-derived modality readiness cards |
| `q` | Quit | Exit the application |

When a checkpoint becomes pending, the Control Center automatically focuses the
Decisions view. Press `1` through `6` to choose the corresponding visible option.
After completion, it focuses Artifacts when a report bundle is available.

Two context-specific editors are available:

- `e` opens the group-design editor at checkpoint 2.1 when proposed groups are
  present;
- `m` opens the typed scATAC FASTQ manifest editor at data-audit checkpoint 1
  when FASTQ input requires it.

Both editors submit into the existing checkpoint resolver. They are not
alternative design or manifest pipelines.

## Checkpoint contract

The Control Center can display and resolve the supervised checkpoints for data
confirmation, group mapping, organism, factor, batch, replicate structure,
assembled design, analysis plan, thresholds, readiness acknowledgement, and
final review. The exact sequence depends on the detected modalities and audit
findings.

For beta modalities, readiness acknowledgement is a real execution gate. The UI
cannot bypass it. Scaffolded modalities remain dispatch-disabled even if their
files appear in the intake directory.

## Artifacts and history

The artifact browser reads the completed report directory and can expose:

- `report.html`;
- `methodology.json`;
- files under `figures/` and `tables/`;
- public claims and their run-ledger linkage.

It is empty until a real report path exists. It does not invent a separate
artifact graph or imply that a missing output ran successfully.

Recent-experiment history is derived from ARIA's local memory and on-disk report
bundles. A displayed report is a review/resume point for the artifact, not a
promise that the original process and MessageBus can be reconstructed.

## Privacy

The Control Center does not change provider or egress policy. If a run must be
air-gapped from its first operation, configure a local provider and set the flag
before launch:

```bash
export ARIA_AIR_GAPPED=1
aria
```

The sensitivity checkpoint occurs after the initial question parse. See
[Installation: Privacy and air-gapped launch](INSTALLATION.md#privacy-and-air-gapped-launch)
for the full boundary.

## Troubleshooting

If `aria` shows the classic interface unexpectedly:

```bash
which aria
python -c 'import aria; print(aria.__version__, aria.__file__)'
python -c 'from aria.ui.cockpit import cockpit_available; print(cockpit_available())'
python -m pip install -e '.[tui]'
```

Confirm that `which aria` and `aria.__file__` belong to the same intended
environment and checkout. A clone at `v4.7.0` and a clone at current `main` both
report package version 4.7.0, so compare `git rev-parse HEAD` as well.
