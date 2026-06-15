"""Presentation-agnostic read-model over the MessageBus (U0).

Both the interactive TUI (:mod:`aria.tui`) and the headless runner
(:mod:`aria.headless`) independently polled the process-global MessageBus and
re-derived the same run state: de-duplicating messages by id, classifying
STATUS/FINDING/ESCALATION, detecting completion (``narrative_agent`` progress
``>= 1.0``), and extracting the report path. A Textual cockpit (U1) would have
been a third copy of that logic.

This module turns the bus state (plus an optional :class:`ExperimentSession`
for the planned-vs-run ledger) into one immutable, serializable snapshot. It is
pure: it reads the bus, never writes, and imports no UI toolkit. Checkpoint
*resolution* stays with the orchestrator; this only *reads* pending checkpoints.

Normalization note (latent-bug fix): agents publish STATUS text under
``payload["message"]`` (``base_agent.publish_status``), while the legacy
consumers read ``payload["status"]``. The snapshot reads either key, so status
text and the ``"Report saved: <path>"`` completion line are surfaced correctly
regardless of which key carried them.
"""

from __future__ import annotations

import re
import os
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

from aria.bus.message_bus import bus, Message, MessageType

# Canonical checkpoint titles. Single source of truth for any presentation
# layer (the legacy `aria.tui` keeps its own copy until U1 converges on this).
CHECKPOINT_TITLES: dict[int | float, str] = {
    1:   "Data Audit Results",
    2:   "Analysis Plan",
    2.1: "Experimental Groups",
    2.2: "Organism",
    2.3: "Experimental Factor",
    2.4: "Batch Effects",
    2.5: "Pseudoreplication Check",
    2.6: "Design Confirmation",
    3:   "Quality Control / Parameter Decision",
    4:   "Preliminary Findings",
    5:   "Final Report Ready",
}

# Pre-dispatch design checkpoints (CP1 audit, CP2.x design). Used for phase.
_AUDIT_CHECKPOINTS = {1}
_DESIGN_CHECKPOINTS = {2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6}

_REPORT_SAVED_RE = re.compile(r"Report saved:\s*(\S+)")

_CONFIDENCE_BUCKETS = ("HIGH", "MEDIUM", "LOW", "INSUFFICIENT")


def status_text(payload: Optional[dict]) -> str:
    """Normalize a STATUS payload's human-readable text.

    Agents write the text under ``message`` (``base_agent.publish_status``); the
    historical consumers read ``status``. Prefer an explicit ``status`` (legacy
    key) when present, otherwise fall back to ``message``.
    """
    if not payload:
        return ""
    return str(payload.get("status") or payload.get("message") or "")


# ── View dataclasses (immutable, serializable) ───────────────────────────────

@dataclass(frozen=True)
class ProgressEvent:
    ts: datetime
    sender: str
    text: str
    progress: float


@dataclass(frozen=True)
class FindingView:
    ts: datetime
    sender: str
    confidence: str          # HIGH | MEDIUM | LOW | INSUFFICIENT
    summary: str
    payload: dict = field(default_factory=dict)


@dataclass(frozen=True)
class CheckpointView:
    message_id: str
    number: int | float | str
    title: str
    question: str
    options: list[str]
    resolved: bool = False
    decision: Any = None
    # Structured escalation context (e.g. CP2.1 carries
    # {"proposed_groups": {group: [samples]}}) for richer editors like U2.
    context: dict = field(default_factory=dict)


@dataclass(frozen=True)
class LedgerNodeView:
    modality: str
    analysis: str
    label: str
    status: str              # ran | skipped | not_run | error
    planned: bool
    divergence: bool
    reason: Optional[str]
    node_id: str


@dataclass(frozen=True)
class ModalityCardView:
    modality: str
    validation_level: str    # validated | alpha | beta | scaffold | unknown
    status: str              # green | yellow | red
    dispatch_policy: str     # allowed | requires_ack | blocked
    reason: str
    findings: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class ResourceView:
    category: str            # env | geneset | motif | genome | celltypist | privacy
    name: str
    status: str              # ready | missing | blocked | pending | info
    detail: str
    path: Optional[str] = None
    action: Optional[str] = None


@dataclass(frozen=True)
class ArtifactView:
    category: str            # report | methodology | figure | table | claim
    name: str
    # present | missing for files; ok | warn | violation for claims.
    status: str
    detail: str
    path: Optional[str] = None


@dataclass(frozen=True)
class ExperimentSnapshot:
    experiment_id: str
    phase: str               # audit | design | dispatch | report | done
    progress: float
    last_status: Optional[ProgressEvent]
    findings_by_confidence: dict[str, list[FindingView]]
    pending_checkpoint: Optional[CheckpointView]
    ledger: list[LedgerNodeView]
    report_path: Optional[str]
    done: bool
    elapsed_s: float
    silent_s: float
    # U3 — per-modality readiness cards (from exp_context); defaulted so older
    # constructions stay valid.
    readiness: list[ModalityCardView] = field(default_factory=list)
    # U4 — local resources and egress policy, presentation-only.
    resources: list[ResourceView] = field(default_factory=list)
    # U6 — report artifacts (report/figs/tables/methodology/claims), read-only.
    artifacts: list[ArtifactView] = field(default_factory=list)


# ── Internal helpers ─────────────────────────────────────────────────────────

def _confidence_label(msg: Message) -> str:
    raw = msg.payload.get("confidence")
    if raw is None:
        raw = getattr(getattr(msg, "confidence", None), "value", None)
    label = str(raw or "medium").upper()
    return label if label in _CONFIDENCE_BUCKETS else "MEDIUM"


def _detect_completion(log: list[Message]) -> tuple[bool, Optional[str]]:
    """Return (done, report_path) from a run's STATUS history.

    Completion is the narrative agent reaching progress ``>= 1.0`` (the report
    is written), or an orchestrator "complete" at ``>= 1.0``. The report path is
    taken from a ``report_path`` payload key or the ``"Report saved: <path>"``
    status line — both now read correctly via :func:`status_text`.
    """
    done = False
    report_path: Optional[str] = None
    for m in log:
        payload = m.payload or {}
        if "report_path" in payload and payload["report_path"]:
            report_path = str(payload["report_path"])
        if m.type != MessageType.STATUS:
            continue
        progress = payload.get("progress") or 0
        try:
            progress = float(progress)
        except (TypeError, ValueError):
            progress = 0.0
        text = status_text(payload)
        if progress >= 1.0 and m.sender == "narrative_agent":
            done = True
            mm = _REPORT_SAVED_RE.search(text)
            if mm:
                report_path = mm.group(1)
        elif (progress >= 1.0 and m.sender == "orchestrator"
              and "complete" in text.lower()):
            done = True
    return done, report_path


def _phase(progress: float, done: bool,
           pending: Optional[CheckpointView]) -> str:
    if done:
        return "done"
    if pending is not None:
        n = pending.number
        if n in _AUDIT_CHECKPOINTS:
            return "audit"
        if n in _DESIGN_CHECKPOINTS:
            return "design"
        if n in (3, 4):
            return "dispatch"
        if n == 5:
            return "report"
    if progress < 0.1:
        return "audit"
    if progress < 0.7:
        return "dispatch"
    return "report"


def _ledger_views(session: Any) -> list[LedgerNodeView]:
    if session is None:
        return []
    exp_ctx = getattr(session, "exp_context", None) or {}
    agent_results = getattr(session, "agent_results", None) or {}
    if not exp_ctx and not agent_results:
        return []
    try:
        from aria.agents.narrative.run_ledger import build_run_ledger
        ledger = build_run_ledger(exp_ctx, agent_results)
    except Exception:
        return []
    out: list[LedgerNodeView] = []
    for e in ledger.get("entries", []) or []:
        if not isinstance(e, dict):
            continue
        out.append(LedgerNodeView(
            modality=str(e.get("modality", "")),
            analysis=str(e.get("analysis", "")),
            label=str(e.get("label", "")),
            status=str(e.get("status", "not_run")),
            planned=bool(e.get("planned")),
            divergence=bool(e.get("divergence")),
            reason=e.get("reason"),
            node_id=str(e.get("node_id", "")),
        ))
    return out


def _readiness_views(session: Any) -> list[ModalityCardView]:
    if session is None:
        return []
    exp_ctx = getattr(session, "exp_context", None) or {}
    cards = exp_ctx.get("readiness_cards")
    if not cards:
        cards = (exp_ctx.get("capability_matrix") or {}).get("cards")
    if not isinstance(cards, dict):
        return []
    out: list[ModalityCardView] = []
    for modality, card in cards.items():
        if not isinstance(card, dict):
            continue
        findings = [
            str(f.get("message", ""))
            for f in card.get("findings", []) or []
            if isinstance(f, dict) and f.get("message")
        ]
        out.append(ModalityCardView(
            modality=str(card.get("modality", modality)),
            validation_level=str(card.get("validation_level", "unknown")),
            status=str(card.get("status", "green")),
            dispatch_policy=str(card.get("dispatch_policy", "allowed")),
            reason=str(card.get("reason", "")),
            findings=findings,
        ))
    out.sort(key=lambda c: c.modality)
    return out


def _all_input_files(exp_ctx: dict) -> list[str]:
    out: list[str] = []
    for files in (exp_ctx.get("modalities") or {}).values():
        if isinstance(files, list):
            out.extend(str(f) for f in files)
    return out


def _has_fastq(exp_ctx: dict) -> bool:
    return any(
        f.lower().endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))
        for f in _all_input_files(exp_ctx)
    )


def _needed_envs(exp_ctx: dict) -> list[str]:
    modalities = set((exp_ctx.get("modalities") or {}).keys())
    has_fastq = _has_fastq(exp_ctx)
    envs: list[str] = []
    if has_fastq and "scRNA" in modalities:
        envs.append("aria-ingestion-env")
    if has_fastq and ({"bulk_RNA", "bulk_RNA_raw"} & modalities):
        envs.append("aria-rnaseq-env")
    if {"scRNA", "bulk_RNA"} & modalities:
        envs.append("aria-rna-env")
    if {"scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"} & modalities:
        envs.append("aria-chromatin-env")
    if {"HiC", "Micro-C"} & modalities:
        envs.append("aria-hic-env")
    return list(dict.fromkeys(envs))


def _path_exists(path: str | os.PathLike | None) -> bool:
    if not path:
        return False
    try:
        return Path(path).expanduser().exists()
    except Exception:
        return False


def _count_children_with_suffix(base: Path, suffix: str) -> int:
    try:
        if not base.exists():
            return 0
        return sum(1 for p in base.glob(f"*/*{suffix}") if p.is_file())
    except Exception:
        return 0


def _resource_setup_views(exp_ctx: dict, agent_results: dict) -> list[ResourceView]:
    setup = agent_results.get("setup_agent") or {}
    needed = _needed_envs(exp_ctx)
    if not needed:
        return [ResourceView(
            category="env", name="Conda stacks", status="pending",
            detail="No modality-specific stack selected yet.",
            action="Confirm input modalities to compute required envs.",
        )]

    status = "ready" if setup.get("status") == "done" else "pending"
    detail = "SetupAgent completed environment checks." if status == "ready" else (
        "SetupAgent has not completed for this run yet."
    )
    warnings = setup.get("warnings") or []
    if warnings:
        detail = f"{detail} {len(warnings)} warning(s) recorded."
    return [ResourceView(
        category="env", name=env_name, status=status, detail=detail,
        action=None if status == "ready" else "SetupAgent will verify/install it.",
    ) for env_name in needed]


def _resource_geneset_view() -> ResourceView:
    try:
        from aria.utils import ora
        base = ora.genesets_dir()
        n_libs = _count_children_with_suffix(base, ".gmt")
        if n_libs:
            return ResourceView(
                category="geneset", name="Local GMT libraries", status="ready",
                detail=f"{n_libs} versioned GMT librar(ies) staged.",
                path=str(base),
            )
        return ResourceView(
            category="geneset", name="Local GMT libraries", status="missing",
            detail="No local GMT libraries staged for offline ORA.",
            path=str(base), action="Run scripts/fetch_genesets.py explicitly.",
        )
    except Exception as exc:
        return ResourceView(
            category="geneset", name="Local GMT libraries", status="missing",
            detail=f"Could not inspect GMT directory: {exc}",
        )


def _resource_motif_view(exp_ctx: dict) -> ResourceView:
    try:
        from aria.utils import motifs
        collection = exp_ctx.get("motif_collection") or motifs.DEFAULT_COLLECTION
        base = motifs.motifs_dir()
        local = motifs.load_local_motif_collection(collection)
        if local is not None:
            meme_path, version = local
            return ResourceView(
                category="motif", name=f"Motifs: {collection}", status="ready",
                detail=f"{version.get('n_motifs', 0)} motifs staged.",
                path=meme_path,
            )
        n_collections = _count_children_with_suffix(base, ".meme")
        detail = "Requested collection is not staged."
        if n_collections:
            detail += f" {n_collections} other collection(s) are present."
        return ResourceView(
            category="motif", name=f"Motifs: {collection}", status="missing",
            detail=detail, path=str(base),
            action="Run scripts/fetch_motifs.py explicitly.",
        )
    except Exception as exc:
        return ResourceView(
            category="motif", name="Motif collections", status="missing",
            detail=f"Could not inspect motif directory: {exc}",
        )


def _resource_genome_view(exp_ctx: dict) -> ResourceView:
    assembly = exp_ctx.get("genome") or exp_ctx.get("assembly")
    if not assembly:
        return ResourceView(
            category="genome", name="Reference genome", status="pending",
            detail="No assembly selected yet.",
            action="Confirm organism/genome during design.",
        )
    try:
        from aria.utils import genomes
        fasta, source = genomes.resolve_local_genome_fasta(str(assembly))
        if fasta:
            return ResourceView(
                category="genome", name=f"Reference genome: {assembly}",
                status="ready", detail=f"Resolved locally via {source}.",
                path=fasta,
            )
        snap_attr = genomes.snapatac2_attr(str(assembly))
        action = (
            f"Stage FASTA under {genomes.GENOME_DIR_ENV}/{assembly}."
        )
        if snap_attr:
            action += " Governed snapatac2 fetch may be offered later."
        return ResourceView(
            category="genome", name=f"Reference genome: {assembly}",
            status="missing", detail="No local FASTA resolved.",
            path=str(genomes.genome_dir()), action=action,
        )
    except Exception as exc:
        return ResourceView(
            category="genome", name=f"Reference genome: {assembly}",
            status="missing", detail=f"Could not inspect genome store: {exc}",
        )


def _celltypist_cache_dirs() -> list[Path]:
    return [
        Path.home() / ".celltypist" / "data" / "models",
        Path.home() / ".celltypist" / "models",
    ]


def _celltypist_model_cached(model_name: str) -> bool:
    explicit = Path(str(model_name)).expanduser()
    if explicit.is_file():
        return True
    return any((d / str(model_name)).is_file() for d in _celltypist_cache_dirs())


def _celltypist_model_from_context(exp_ctx: dict) -> tuple[str, str]:
    model = (exp_ctx.get("celltypist_model")
             or exp_ctx.get("celltypist_model_path"))
    organism = str(exp_ctx.get("organism") or "")
    hint = str(exp_ctx.get("celltypist_tissue_hint")
               or exp_ctx.get("tissue_hint")
               or exp_ctx.get("tissue")
               or "").lower().strip()
    try:
        from aria.scripts.rna_celltypist import _resolve_celltypist_model
        resolved, source, _warning = _resolve_celltypist_model(
            model, organism, hint)
        return str(resolved), str(source)
    except Exception:
        if model:
            return str(model), "explicit"
        return "Immune_All_Low.pkl", "default_immune_fallback"


def _resource_celltypist_view(exp_ctx: dict) -> list[ResourceView]:
    if "scRNA" not in set((exp_ctx.get("modalities") or {}).keys()):
        return []
    model, source = _celltypist_model_from_context(exp_ctx)
    cached = _celltypist_model_cached(model)
    status = "ready" if cached else "missing"
    detail = f"Model source: {source}; cache {'hit' if cached else 'miss'}."
    action = None if cached else "Stage model locally or allow governed model-hub egress."
    return [ResourceView(
        category="celltypist", name=f"CellTypist: {model}",
        status=status, detail=detail,
        path=str(Path(model).expanduser()) if _path_exists(model) else None,
        action=action,
    )]


def _resource_privacy_view(exp_ctx: dict) -> ResourceView:
    enabled = bool(exp_ctx.get("air_gapped"))
    if not enabled:
        try:
            from aria.utils.privacy import air_gapped_enabled
            enabled = air_gapped_enabled()
        except Exception:
            enabled = False
    return ResourceView(
        category="privacy", name="Network egress", status="blocked" if enabled else "info",
        detail=("Air-gapped mode is active; network fetches are blocked."
                if enabled else "Network egress is allowed unless a run opts into air-gap."),
        action=("Use local staged resources/caches."
                if enabled else "Use air-gapped mode for sensitive data."),
    )


def _resource_views(session: Any) -> list[ResourceView]:
    exp_ctx = getattr(session, "exp_context", None) or {}
    agent_results = getattr(session, "agent_results", None) or {}
    resources: list[ResourceView] = []
    resources.extend(_resource_setup_views(exp_ctx, agent_results))
    resources.append(_resource_geneset_view())
    resources.append(_resource_motif_view(exp_ctx))
    resources.append(_resource_genome_view(exp_ctx))
    resources.extend(_resource_celltypist_view(exp_ctx))
    resources.append(_resource_privacy_view(exp_ctx))
    return resources


# ── U6 artifact browser ──────────────────────────────────────────────────────

# A compiled claim's run-ledger status maps to an artifact-browser status. A
# violation (non-descriptive claim citing a not-run/error node) is detected
# upstream by `link_claims_to_ledger` and surfaced here verbatim, so U6 reuses
# the report-builder's evidence graph rather than re-deriving tier logic.
_CLAIM_LEDGER_STATUS = {
    "ran": "ok",
    "skipped": "warn",
    "not_run": "warn",
    "error": "warn",
    "no_ledger_node": "warn",
}


def _read_methodology(report_dir: Path) -> Optional[dict]:
    path = report_dir / "methodology.json"
    if not path.is_file():
        return None
    try:
        import json
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _list_dir_artifacts(report_dir: Path, sub: str,
                        category: str) -> list[ArtifactView]:
    base = report_dir / sub
    out: list[ArtifactView] = []
    try:
        if not base.is_dir():
            return out
        for p in sorted(base.iterdir()):
            if not p.is_file():
                continue
            try:
                size = p.stat().st_size
            except Exception:
                size = 0
            out.append(ArtifactView(
                category=category, name=p.name, status="present",
                detail=f"{size:,} bytes", path=str(p),
            ))
    except Exception:
        return out
    return out


def _claim_artifact_views(methodology: Optional[dict]) -> list[ArtifactView]:
    if not methodology:
        return []
    claims = methodology.get("claims")
    if not isinstance(claims, list):
        return []
    linkage = ((methodology.get("run_ledger") or {}).get("claim_linkage")) or {}
    violation_ids = {
        v.get("claim_id")
        for v in (linkage.get("violations") or [])
        if isinstance(v, dict)
    }
    out: list[ArtifactView] = []
    for c in claims:
        if not isinstance(c, dict):
            continue
        cid = c.get("claim_id")
        ledger_status = str(c.get("ledger_status") or "no_ledger_node")
        if cid in violation_ids:
            status = "violation"
        else:
            status = _CLAIM_LEDGER_STATUS.get(ledger_status, "warn")
        tier = str(c.get("tier") or "—")
        conf = str(c.get("confidence") or "—")
        text = str(c.get("text") or "").strip()
        detail = f"tier={tier} · {conf} · ledger={ledger_status}"
        if text:
            detail = f"{detail} · {text}"
        out.append(ArtifactView(
            category="claim", name=str(cid or "claim"),
            status=status, detail=detail,
        ))
    return out


def _artifact_views(report_path: Optional[str]) -> list[ArtifactView]:
    """Surface the on-disk report bundle: report/methodology/figures/tables/claims.

    Read-only and honest: returns ``[]`` until the report path is known, and
    reflects exactly what is on disk. Claims and their ledger status are read
    from the report-builder's ``methodology.json`` (which already ran
    ``link_claims_to_ledger``), so U6 never invents a separate artifact graph.
    """
    if not report_path:
        return []
    report_file = Path(report_path).expanduser()
    report_dir = report_file.parent
    if not report_dir.is_dir():
        return []

    out: list[ArtifactView] = []
    report_present = report_file.is_file()
    out.append(ArtifactView(
        category="report", name=report_file.name,
        status="present" if report_present else "missing",
        detail=("Narrative HTML report." if report_present
                else "Report file not found on disk."),
        path=str(report_file) if report_present else None,
    ))

    meth_path = report_dir / "methodology.json"
    meth_present = meth_path.is_file()
    out.append(ArtifactView(
        category="methodology", name="methodology.json",
        status="present" if meth_present else "missing",
        detail=("Provenance, thresholds, tools, claims, run-ledger."
                if meth_present else "Methodology manifest not found."),
        path=str(meth_path) if meth_present else None,
    ))

    out.extend(_list_dir_artifacts(report_dir, "figures", "figure"))
    out.extend(_list_dir_artifacts(report_dir, "tables", "table"))
    out.extend(_claim_artifact_views(_read_methodology(report_dir)))
    return out


# ── Public API ───────────────────────────────────────────────────────────────

def build_snapshot(experiment_id: str, *,
                   session: Any = None,
                   start_time: Optional[datetime] = None,
                   now: Optional[datetime] = None) -> ExperimentSnapshot:
    """Build an immutable, UI-agnostic snapshot of one experiment's run state.

    Reads the process-global bus, scoped to ``experiment_id``. ``session`` (an
    :class:`ExperimentSession`) is optional and only used to reconcile the
    planned-vs-run ledger. ``start_time`` anchors ``elapsed_s``; ``now`` is
    injectable for deterministic tests.
    """
    now = now or datetime.now()
    log = bus.get_log(experiment_id)

    # Latest STATUS message -> progress + last_status banner.
    last_status: Optional[ProgressEvent] = None
    progress = 0.0
    last_msg_ts: Optional[datetime] = None
    for m in log:
        if last_msg_ts is None or m.timestamp > last_msg_ts:
            last_msg_ts = m.timestamp
        if m.type == MessageType.STATUS:
            payload = m.payload or {}
            p = payload.get("progress") or 0
            try:
                p = float(p)
            except (TypeError, ValueError):
                p = 0.0
            last_status = ProgressEvent(
                ts=m.timestamp, sender=m.sender,
                text=status_text(payload), progress=p,
            )
            progress = p

    # Findings grouped by confidence (served from the bus index).
    grouped: dict[str, list[FindingView]] = {b: [] for b in _CONFIDENCE_BUCKETS}
    for m in bus.get_findings(experiment_id):
        grouped[_confidence_label(m)].append(FindingView(
            ts=m.timestamp, sender=m.sender,
            confidence=_confidence_label(m),
            summary=str((m.payload or {}).get("summary", "")),
            payload=dict(m.payload or {}),
        ))

    # First unresolved checkpoint (both legacy loops act on pending[0]).
    pending: Optional[CheckpointView] = None
    pendings = [
        m for m in bus.get_pending_checkpoints(experiment_id=experiment_id)
        if not (m.payload or {}).get("resolved")
    ]
    if pendings:
        m = pendings[0]
        payload = m.payload or {}
        num = payload.get("checkpoint", m.checkpoint)
        pending = CheckpointView(
            message_id=m.id,
            number=num if num is not None else "?",
            title=CHECKPOINT_TITLES.get(num, f"Checkpoint {num}"),
            question=str(payload.get("question", "")),
            options=list(payload.get("options", ["Continue", "Cancel"])),
            resolved=bool(payload.get("resolved")),
            decision=payload.get("user_decision"),
            context=dict(payload.get("context") or {}),
        )

    done, report_path = _detect_completion(log)
    phase = _phase(progress, done, pending)

    elapsed_s = (now - start_time).total_seconds() if start_time else 0.0
    silent_s = (now - last_msg_ts).total_seconds() if last_msg_ts else 0.0

    return ExperimentSnapshot(
        experiment_id=experiment_id,
        phase=phase,
        progress=progress,
        last_status=last_status,
        findings_by_confidence=grouped,
        pending_checkpoint=pending,
        ledger=_ledger_views(session),
        report_path=report_path,
        done=done,
        elapsed_s=elapsed_s,
        silent_s=silent_s,
        readiness=_readiness_views(session),
        resources=_resource_views(session),
        artifacts=_artifact_views(report_path),
    )
