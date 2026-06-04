"""ARIA report-section + provenance HTML builders (P2-8 follow-up).

Pure, deterministic helpers extracted from `narrative_agent.py` (all were
`@staticmethod`). They are re-aliased on `NarrativeAgent` so call sites and the
public surface are unchanged. No LLM, no `self` — just structured-data ->
HTML/text rendering and provenance collection."""

from __future__ import annotations

import html as _html
import re
from pathlib import Path


def _avg_pct_passed(qc_info: dict) -> float:
        """Average percentage of reads passing QC across samples."""
        samples = qc_info.get("samples", [])
        if not samples:
            return 0.0
        pct_list = [s.get("pct_passed", 0) for s in samples]
        return sum(pct_list) / max(len(pct_list), 1)


def _collect_tool_versions(packages: tuple[str, ...]) -> dict:
        tools = {}
        try:
            from importlib.metadata import version, PackageNotFoundError
            for pkg in packages:
                try:
                    tools[pkg] = version(pkg)
                except PackageNotFoundError:
                    tools[pkg] = "not installed"
        except Exception:
            tools = {pkg: "not installed" for pkg in packages}

        locked = _tool_versions_from_lockfiles(packages)
        for pkg, locked_version in locked.items():
            if tools.get(pkg) in (None, "not installed"):
                tools[pkg] = locked_version
            elif tools[pkg] != locked_version:
                tools[f"{pkg}_lock"] = locked_version
        return tools


def _tool_versions_from_lockfiles(packages: tuple[str, ...]) -> dict:
        root = Path(__file__).resolve().parents[2]
        package_set = {p.lower(): p for p in packages}
        versions = {}
        pip_locks = sorted((root / "envs").glob("*.pip.lock"))
        for lock in pip_locks:
            try:
                for line in lock.read_text(encoding="utf-8").splitlines():
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    if "==" in stripped:
                        name, ver = stripped.split("==", 1)
                    else:
                        continue
                    canonical = package_set.get(name.strip().lower())
                    if canonical and canonical not in versions:
                        versions[canonical] = ver.strip()
            except Exception:
                continue

        conda_locks = sorted((root / "envs").glob("*.linux-64.lock"))
        for lock in conda_locks:
            try:
                for line in lock.read_text(encoding="utf-8").splitlines():
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#") or stripped == "@EXPLICIT":
                        continue
                    parsed = _package_version_from_conda_url(
                        stripped, package_set
                    )
                    if parsed:
                        canonical, ver = parsed
                        versions.setdefault(canonical, ver)
            except Exception:
                continue
        return versions


def _package_version_from_conda_url(
        url: str, package_set: dict[str, str]
    ) -> tuple[str, str] | None:
        import re
        filename = url.rsplit("/", 1)[-1]
        filename = re.sub(r"\.(conda|tar\.bz2)$", "", filename)
        for lower_name, canonical in package_set.items():
            prefix = f"{lower_name}-"
            if filename.lower().startswith(prefix):
                rest = filename[len(prefix):]
                version = rest.split("-", 1)[0]
                return canonical, version
        return None


def _build_calibration_badge(calibration: dict | None) -> str:
    """W-CALIB: render the numerical-calibration badge for the report.

    ``calibration`` is the manifest produced by
    ``aria.benchmarks.run_calibration_suite`` (recall + empirical FDR + the
    label-permutation negative-control false-positive rate, for the bulk and
    pseudobulk DE paths). It is a property of the ARIA build, not of the user's
    dataset, so it is only present when a real calibration run was attached.

    Honesty contract (no fabrication): when no manifest is attached, the badge
    states that calibration was NOT measured in this report environment and that
    the DE paths are covered by ARIA's calibration test suite (the CI release
    gate) — it shows NO metric. When a manifest IS attached it shows exactly the
    measured numbers and their pass/fail status.
    """
    if not isinstance(calibration, dict) or not calibration.get("measured"):
        return (
            "<h3>Numerical Calibration (W-CALIB)</h3>"
            "<p><span class='badge'>not measured in this run</span> "
            "DE numerical calibration (recall, empirical FDR, and label-permutation "
            "negative controls on the bulk &amp; pseudobulk paths) is exercised by "
            "ARIA's calibration test suite / CI release gate, not re-run for each "
            "report. No calibration metric is asserted for this run.</p>"
        )

    status = str(calibration.get("status", "")).lower()
    css = {"pass": "high", "fail": "low", "error": "insuff"}.get(status, "medium")
    label = {"pass": "PASS", "fail": "FAIL", "error": "ERROR"}.get(status, status.upper() or "—")
    summary = calibration.get("summary", {}) or {}

    def _fmt(key: str) -> str:
        v = summary.get(key)
        try:
            return f"{float(v):.3f}"
        except (TypeError, ValueError):
            return "—"

    rows = [
        ("Bulk recall", _fmt("bulk_recall")),
        ("Bulk empirical FDR", _fmt("bulk_empirical_fdr")),
        ("Bulk null false-positive rate", _fmt("bulk_null_fpr")),
        ("Pseudobulk recall", _fmt("pseudobulk_recall")),
        ("Pseudobulk empirical FDR", _fmt("pseudobulk_empirical_fdr")),
        ("Pseudobulk null false-positive rate", _fmt("pseudobulk_null_fpr")),
    ]
    body = "".join(
        f"<tr><td>{_html.escape(name)}</td><td><code>{_html.escape(val)}</code></td></tr>"
        for name, val in rows
    )
    seed = _html.escape(str(calibration.get("seed", "")))
    return (
        "<h3>Numerical Calibration (W-CALIB)</h3>"
        f"<p><span class='badge {css}'>{label}</span> measured on synthetic "
        "ground-truth data (recovery + label-permutation negative control). "
        f"Seed <code>{seed}</code>. The null false-positive rate is the empirical "
        "type-I rate when condition labels are permuted; it should sit at or below "
        "the nominal alpha.</p>"
        "<table><tr><th>Metric</th><th>Value</th></tr>"
        + body
        + "</table>"
    )


def _build_run_ledger_section(run_ledger: dict | None) -> str:
        """P-LEDGER: render the planned-vs-run manifest. Any analysis the plan
        called for that did not run is flagged as a divergence."""
        entries = (run_ledger or {}).get("entries", []) or []
        if not entries:
            return ""
        n_div = (run_ledger or {}).get("n_divergences", 0)
        rows = []
        for e in entries:
            planned = "yes" if e.get("planned") else "no"
            status = str(e.get("status", ""))
            reason = e.get("reason")
            flag = " ⚠ planned but not run" if e.get("divergence") else ""
            detail = f" ({_html.escape(str(reason))})" if reason else ""
            rows.append(
                "<tr>"
                f"<td>{_html.escape(str(e.get('label', e.get('analysis', ''))))}</td>"
                f"<td>{_html.escape(planned)}</td>"
                f"<td><code>{_html.escape(status)}</code>{detail}"
                f"{_html.escape(flag)}</td>"
                "</tr>"
            )
        header = (
            f"<h3>Run Ledger (planned vs executed)</h3>"
            f"<p><em>{n_div} plan/execution divergence(s).</em></p>"
        )
        return (
            header
            + "<table><tr><th>Analysis</th><th>Planned</th>"
              "<th>Status</th></tr>"
            + "".join(rows)
            + "</table>"
        )


def _build_raw_ingestion_section(agent_results: dict,
                                     exp_ctx_records: list | None = None) -> str:
        records = []
        raw = (agent_results or {}).get("raw_ingestion_agent", {}) or {}
        records.extend(raw.get("records", []) or [])
        records.extend(exp_ctx_records or [])
        if not records:
            return ""
        rows = []
        for rec in records:
            source = rec.get("source_directory") or rec.get("mode", "")
            output = rec.get("output_h5ad") or ""
            output_hash = rec.get("output_sha256") or ""
            blockers = "; ".join(rec.get("blockers", [])[:3]) if rec.get("blockers") else ""
            rows.append(
                "<tr>"
                f"<td>{_html.escape(str(rec.get('mode', '')))}</td>"
                f"<td><code>{_html.escape(str(source))}</code></td>"
                f"<td><code>{_html.escape(str(output))}</code></td>"
                f"<td><code>{_html.escape(str(output_hash))}</code></td>"
                f"<td>{_html.escape(blockers)}</td>"
                "</tr>"
            )
        return (
            "<h3>Raw Ingestion</h3>"
            "<table><tr><th>Mode</th><th>Source</th><th>Generated h5ad</th>"
            "<th>Output SHA-256</th><th>Blockers</th></tr>"
            + "".join(rows)
            + "</table>"
        )


def _collect_param_hashes(obj, prefix: str = "") -> list[tuple[str, str]]:
        rows: list[tuple[str, str]] = []
        if isinstance(obj, dict):
            digest = obj.get("params_sha256")
            if digest:
                rows.append((prefix or "root", str(digest)))
            for key, val in obj.items():
                if key == "params_sha256":
                    continue
                child_prefix = f"{prefix}.{key}" if prefix else str(key)
                rows.extend(_collect_param_hashes(val, child_prefix))
        elif isinstance(obj, list):
            for idx, val in enumerate(obj):
                rows.extend(
                    _collect_param_hashes(val, f"{prefix}[{idx}]")
                )
        seen = set()
        unique = []
        for label, digest in rows:
            pair = (label, digest)
            if pair not in seen:
                seen.add(pair)
                unique.append(pair)
        return unique


def _build_lockfile_section() -> str:
        root = Path(__file__).resolve().parents[2]
        env_dir = root / "envs"
        # Pair every conda lockfile with its sibling pip lock (if present).
        # `conda create --name X --file <env>.linux-64.lock` reproduces the
        # conda side byte-by-byte; `pip install -r <env>.pip.lock` covers
        # the pip side. Both are needed for full peer-reviewable reproduction.
        lockfiles = sorted(env_dir.glob("*.linux-64.lock"))
        pip_locks = sorted(env_dir.glob("*.pip.lock"))
        if not lockfiles and not pip_locks:
            return (
                "<div class='warning'>No conda lockfiles found in "
                "<code>envs/*.linux-64.lock</code>. Run "
                "<code>scripts/generate_locks.sh</code> before tagging the release.</div>"
            )
        blocks = []
        for lock in lockfiles:
            try:
                text = lock.read_text(encoding="utf-8")[:20000]
            except Exception as exc:
                text = f"Could not read lockfile: {exc}"
            blocks.append(
                f"<details><summary>{_html.escape(lock.name)}</summary>"
                f"<pre>{_html.escape(text)}</pre></details>"
            )
        for pip_lock in pip_locks:
            try:
                text = pip_lock.read_text(encoding="utf-8")[:20000]
            except Exception as exc:
                text = f"Could not read pip lockfile: {exc}"
            blocks.append(
                f"<details><summary>{_html.escape(pip_lock.name)} "
                f"(pip side; use after the conda lock)</summary>"
                f"<pre>{_html.escape(text)}</pre></details>"
            )
        # If a conda lock has no sibling pip lock, note that pip was empty
        # (rather than missing) so a reviewer doesn't suspect omission.
        for lock in lockfiles:
            sibling_pip = env_dir / lock.name.replace(".linux-64.lock", ".pip.lock")
            if not sibling_pip.exists():
                blocks.append(
                    f"<div><em>{_html.escape(lock.stem)}: no pip packages "
                    f"in env (pip lock not emitted).</em></div>"
                )
        return "".join(blocks)


def _build_slug(intent: dict, exp_ctx: dict) -> str:
        """
        Build a short URL-safe slug from biological entities or the question.
        Examples:
          {entities: ['GeneA', 'GeneB']} -> "genea_geneb"
          {summary: "Effect of treatment on cells"} -> "effect_treatment_cells"
        """
        import re as _re
        # Prefer biological entities (most specific)
        entities = intent.get("biological_entities", []) or []
        if entities:
            parts = [_re.sub(r"[^a-z0-9]+", "", str(e).lower())
                     for e in entities[:3]]
            parts = [p for p in parts if p and p not in
                     ("cells", "cell", "h9", "h1", "wildtype", "wt")]
            if parts:
                return "_".join(parts)[:40]

        # Fallback: first few words of the summary or user_question
        text = (str(intent.get("summary", "")) or
                str(exp_ctx.get("user_question", ""))).lower()
        if text:
            words = _re.findall(r"[a-z0-9]+", text)
            stop  = {"the", "a", "an", "of", "to", "in", "and", "or", "for",
                     "with", "is", "are", "what", "how", "this", "that",
                     "vs", "versus", "between", "from", "on"}
            keep  = [w for w in words if w not in stop and len(w) > 2][:3]
            if keep:
                return "_".join(keep)[:40]
        return ""


def _plain_text_to_html(text: str) -> str:
        """Escape plain report text and preserve line breaks."""
        escaped = _html.escape(str(text or ""))
        return escaped.replace("\n", "<br>")


def _format_finding_summary(finding: dict) -> str:
        """Human-readable finding text; no raw Python dicts in reports."""
        summary = finding.get("summary")
        if isinstance(summary, dict):
            if summary.get("error"):
                text = f"Skipped: {summary['error']}"
            else:
                text = "; ".join(
                    f"{k}={v}" for k, v in summary.items()
                    if v is not None
                )
        elif summary:
            text = str(summary)
        elif finding.get("error"):
            text = f"Skipped: {finding['error']}"
        else:
            text = str(finding)

        if "Integration requires at least 2 modalities" in text:
            text = (
                "Integration was skipped because fewer than two modalities "
                "were available for this run."
            )
        text = re.sub(r"\s+", " ", text).strip()
        if len(text) > 200:
            text = text[:197].rstrip() + "..."
        return _html.escape(text)


def _guard_bulk_interpretation(text: str) -> str:
        """Downgrade unsupported causal/mechanistic wording in bulk reports."""
        replacements = {
            r"\bacts as an? ([a-z -]+)\b":
                r"is consistent with a \1-associated expression pattern",
            r"\bsuppressing\b":
                "associated with lower expression of",
            r"\bposition ([^.]+?) as hierarchical gatekeepers\b":
                r"support the hypothesis that \1 influence",
            r"\bhierarchical gatekeepers\b":
                "candidate upstream modulators",
            r"\b([A-Za-z0-9_.-]+) enforces ([A-Za-z0-9_ -]+)\b":
                r"\1 perturbation is associated with \2-related changes",
            r"\b([A-Za-z0-9_.-]+) drives ([A-Za-z0-9_ -]+)\b":
                r"\1 perturbation is associated with \2",
        }
        guarded = text
        for pattern, repl in replacements.items():
            guarded = re.sub(pattern, repl, guarded, flags=re.IGNORECASE)

        causal_terms = (
            "These interpretations are hypotheses from differential "
            "expression and pathway enrichment, not direct evidence of "
            "transcriptional causality; orthogonal perturbation, chromatin, "
            "or binding assays are required to distinguish direct regulation "
            "from secondary state changes."
        )
        if "not direct evidence of transcriptional causality" not in guarded:
            guarded = f"{guarded}\n\n{causal_terms}"
        return guarded
