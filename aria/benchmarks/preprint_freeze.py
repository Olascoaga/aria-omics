"""Fail-closed inventory for the ARIA preprint-v1 evidence freeze.

Historical benchmark files are useful context, but they never count as freshly
regenerated freeze evidence. A lane becomes ``verified`` only through a receipt
under the dedicated ``preprint_v1`` root, bound to the current clean commit and
the hashes of that lane's emitted artifacts.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
from typing import Any, Mapping

from aria.version import collect_version_metadata


SCHEMA_VERSION = "aria.preprint_freeze.inventory.v1"
CLAIM_IDS = tuple(f"claim_{i}" for i in range(1, 8))


def _lane(
    lane_id: str,
    claims: tuple[str, ...],
    requirement: str,
    command: str | None,
    environment: str,
    resources: tuple[str, ...] = (),
    artifacts: tuple[str, ...] = (),
    *,
    implementation: str = "available",
    evidence_kind: str = "artifact",
) -> dict[str, Any]:
    return {
        "lane_id": lane_id,
        "claims": list(claims),
        "requirement": requirement,
        "command": command,
        "environment": environment,
        "resources": list(resources),
        "expected_artifacts": list(artifacts),
        "implementation": implementation,
        "evidence_kind": evidence_kind,
        "required_for_freeze": True,
    }


LANES: tuple[dict[str, Any], ...] = (
    _lane(
        "c1_a1_synthetic_bulk_de", ("claim_1",),
        "Fresh A1 synthetic-truth bulk-DE calibration and LFC frontier.",
        "conda run -n aria-rna-env python scripts/run_a1_bulk_de_benchmark.py "
        "--output-dir docs/benchmark_results/preprint_v1/claim_1/a1_synthetic "
        "--manifest-name a1_bulk_de.json --figure-name fig1_a1_bulk_de.svg",
        "aria-rna-env", ("env:aria-rna-env",),
        ("claim_1/a1_synthetic/a1_bulk_de.json",
         "claim_1/a1_synthetic/fig1_a1_bulk_de.svg"),
    ),
    _lane(
        "c1_a1_external_comparators", ("claim_1",),
        "Fresh DESeq2/edgeR/limma comparator execution on the identical A1 matrix.",
        "python scripts/run_a1_external_comparators.py --output-dir "
        "docs/benchmark_results/preprint_v1/claim_1/a1_external "
        "--manifest-name a1_external_comparators.json",
        "aria-env -> aria-bench-env", ("env:aria-env", "env:aria-bench-env"),
        (
            "claim_1/a1_external/a1_external_comparators.json",
            "claim_1/a1_external/inputs/a1_counts.tsv",
            "claim_1/a1_external/inputs/a1_metadata.tsv",
            "claim_1/a1_external/inputs/a1_truth.tsv",
            "claim_1/a1_external/r_outputs/r_comparators.json",
            "claim_1/a1_external/r_outputs/deseq2.tsv",
            "claim_1/a1_external/r_outputs/edgeR_QLF.tsv",
            "claim_1/a1_external/r_outputs/limma_voom.tsv",
        ),
    ),
    _lane(
        "c1_seqc_maqc_taqman", ("claim_1",),
        "Fresh ARIA-vs-TaqMan SEQC/MAQC reference validation.",
        "conda run -n aria-rna-env python scripts/run_a1_seqc_maqc_benchmark.py "
        "--bundle ~/.aria/benchmarks/seqc_maqc_BGI --output-dir "
        "docs/benchmark_results/preprint_v1/claim_1/seqc",
        "aria-rna-env", ("env:aria-rna-env", "data:seqc_bgi"),
        (
            "claim_1/seqc/a1_seqc_maqc_v4.5.5.json",
            "claim_1/seqc/fig1_a1_seqc_maqc_v4.5.5.svg",
        ),
    ),
    _lane(
        "c1_seqc_multisite", ("claim_1",),
        "Fresh five-site SEQC log2FC reproducibility matrix.",
        "conda run -n aria-rna-env python scripts/run_a1_seqc_multisite_benchmark.py "
        "--output-dir docs/benchmark_results/preprint_v1/claim_1/seqc_multisite",
        "aria-rna-env", ("env:aria-rna-env", "data:seqc_multisite"),
        (
            "claim_1/seqc_multisite/a1_seqc_multisite_v4.5.5.json",
            "claim_1/seqc_multisite/fig1_a1_seqc_multisite_v4.5.5.svg",
        ),
    ),
    _lane(
        "c1_ercc_dose_response", ("claim_1",),
        "Fresh ERCC fold-change and dynamic-range recovery.",
        "conda run -n aria-rna-env python scripts/run_a1_ercc_dose_response.py "
        "--bundle ~/.aria/benchmarks/seqc_maqc_BGI --output-dir "
        "docs/benchmark_results/preprint_v1/claim_1/ercc",
        "aria-rna-env", ("env:aria-rna-env", "data:seqc_ercc"),
        (
            "claim_1/ercc/a1_ercc_dose_response_v4.5.5.json",
            "claim_1/ercc/fig1_a1_ercc_dose_response_v4.5.5.svg",
        ),
    ),
    _lane(
        "c1_h9_fastq_e2e", ("claim_1",),
        "Clean-checkout FASTQ-to-report bulk-RNA E2E with frozen checkpoint policy.",
        "conda run -n aria-env python scripts/run_c1_h9_fastq_e2e.py "
        "--output-dir docs/benchmark_results/preprint_v1/claim_1/h9_e2e",
        "aria-env -> aria-rnaseq-env + aria-rna-env",
        ("env:aria-env", "env:aria-rnaseq-env", "env:aria-rna-env", "data:h9_fastq"),
        (
            "claim_1/h9_e2e/capsule.json",
            "claim_1/h9_e2e/report.html",
            "claim_1/h9_e2e/methodology.json",
            "claim_1/h9_e2e/de_results.tsv",
            "claim_1/h9_e2e/fig1_h9_bulk_de.svg",
        ),
        evidence_kind="e2e_capsule",
    ),
    _lane(
        "c2_a2_synthetic_pseudobulk", ("claim_2",),
        "Fresh donor-aware recovery and donor-heterogeneity null.",
        "conda run -n aria-rna-env python scripts/run_a2_pseudobulk_benchmark.py "
        "--output-dir docs/benchmark_results/preprint_v1/claim_2/a2_synthetic "
        "--manifest-name a2_pseudobulk.json --figure-name fig2_a2_pseudobulk.svg "
        "--artifact-version preprint_v1",
        "aria-rna-env", ("env:aria-rna-env",),
        ("claim_2/a2_synthetic/a2_pseudobulk.json",
         "claim_2/a2_synthetic/fig2_a2_pseudobulk.svg"),
    ),
    _lane(
        "c2_kang_muscat", ("claim_2",),
        "Fresh ARIA-vs-muscat/edgeR comparison on identical Kang pseudobulk matrices.",
        "conda run -n aria-rna-env python scripts/run_a2_external_comparators.py "
        "--export-dir ~/.aria/benchmarks/kang_muscat --output-dir "
        "docs/benchmark_results/preprint_v1/claim_2/kang_muscat",
        "aria-rna-env", ("env:aria-rna-env", "data:kang_muscat"),
        (
            "claim_2/kang_muscat/a2_kang_muscat_v4.5.5.json",
            "claim_2/kang_muscat/fig2_a2_kang_muscat_v4.5.5.svg",
        ),
    ),
    _lane(
        "c2_pairing_and_technical_replicates", ("claim_2",),
        "Regression proof for partial pairing and explicit technical replicate units.",
        "conda run -n aria-env python -m pytest "
        "tests/test_preprint_audit_b1_technical_replicates.py "
        "tests/test_preprint_audit_b3_pairing.py -q",
        "aria-env", ("env:aria-env",), evidence_kind="test_receipt",
    ),
    _lane(
        "c2_scatac_donor_aware", ("claim_2",),
        "Fresh donor-aware scATAC DA concordance on the real multi-donor matrix.",
        "conda run -n aria-rna-env python "
        "scripts/run_scatac_multisample_da_concordance.py "
        "--data-dir ~/Samael/Single/prueba/muon_processed "
        "--young hc937,hc5579,hc13344,hc935,hc5614 "
        "--old hc26,hc98,hc35,hc40,hc73 "
        "--rna-env aria-rna-env --bench-env aria-bench-env "
        "--output-dir docs/benchmark_results/preprint_v1/claim_2/scatac_donor "
        "--manifest-name scatac_donor_da_concordance.json "
        "--work-dir /tmp/c2_scatac_donor",
        "aria-rna-env + aria-bench-env",
        ("env:aria-rna-env", "env:aria-bench-env", "data:gse278576"),
        ("claim_2/scatac_donor/scatac_donor_da_concordance.json",),
        evidence_kind="external_concordance",
    ),
    _lane(
        "c3_b1_design_governance", ("claim_3",),
        "Fresh adversarial design inference/block/escalation confusion matrix.",
        "conda run -n aria-env python scripts/run_b1_design_governance.py "
        "--output-dir docs/benchmark_results/preprint_v1/claim_3/b1_design",
        "aria-env", ("env:aria-env",),
        (
            "claim_3/b1_design/b1_design_governance_v4.5.5.json",
            "claim_3/b1_design/fig3_b1_design_governance_v4.5.5.svg",
        ),
    ),
    _lane(
        "c3_blind_multifactorial_corpus", ("claim_3",),
        "Independent blind gold for a multifactorial held-out design corpus.",
        None, "human review", ("human:design_gold",),
        implementation="missing", evidence_kind="human_gold",
    ),
    _lane(
        "c4_b2_governance_ablation", ("claim_4",),
        "Fresh four-arm governed/guards-off/naive/template claim benchmark.",
        "conda run -n aria-env python scripts/run_b2_claim_narrative.py "
        "--output-dir docs/benchmark_results/preprint_v1/claim_4/b2_claim",
        "aria-env", ("env:aria-env",),
        (
            "claim_4/b2_claim/b2_claim_narrative_v4.5.5.json",
            "claim_4/b2_claim/fig4_b2_claim_narrative_v4.5.5.svg",
        ),
    ),
    _lane(
        "c4_b2_multi_annotator_gold", ("claim_4",),
        "At least two independent annotators plus adjudication and agreement metrics.",
        "python scripts/run_b2_annotation_kit.py score "
        "docs/benchmark_results/preprint_v1/human/b2_annotator_1.csv "
        "docs/benchmark_results/preprint_v1/human/b2_annotator_2.csv "
        "--adjudicator docs/benchmark_results/preprint_v1/human/b2_adjudicator.csv "
        "--out docs/benchmark_results/preprint_v1/claim_4/b2_human_gold.json",
        "aria-env", ("env:aria-env", "human:b2_annotators"),
        ("claim_4/b2_human_gold.json",), evidence_kind="human_gold",
    ),
    _lane(
        "c4_report_e2e_false_narrative", ("claim_4",),
        "E2E benchmark through bus, compiler, report and independent human scoring.",
        None, "aria-env + human review", ("env:aria-env", "human:b2_annotators"),
        implementation="missing", evidence_kind="e2e_human_gold",
    ),
    _lane(
        "c5_b4_null_narrative", ("claim_5",),
        "Fresh null-evidence narrative benchmark.",
        "conda run -n aria-env python scripts/run_b4_null_narrative.py "
        "--output-dir docs/benchmark_results/preprint_v1/claim_5/b4_null",
        "aria-env", ("env:aria-env",),
        (
            "claim_5/b4_null/b4_null_narrative_v4.5.5.json",
            "claim_5/b4_null/fig5_b4_null_narrative_v4.5.5.svg",
        ),
    ),
    _lane(
        "c5_semantic_negation", ("claim_5",),
        "Semantic positive/null contradiction and negation regression.",
        "conda run -n aria-env python -m pytest "
        "tests/test_preprint_audit_c2_semantic_facts.py -q",
        "aria-env", ("env:aria-env",), evidence_kind="test_receipt",
    ),
    _lane(
        "c5_multimodal_null_permutations", ("claim_5",),
        "Fresh RNA and ATAC label-permutation null reports through the public compiler.",
        None, "aria-rna-env + aria-chromatin-env",
        ("env:aria-rna-env", "env:aria-chromatin-env"),
        implementation="missing", evidence_kind="multimodal_e2e",
    ),
    _lane(
        "c6_b3_in_process_matrix", ("claim_6",),
        "Fresh 3-provider x 3-repetition fail-closed-monotone compiler matrix.",
        "conda run -n aria-env python -m pytest "
        "tests/test_preprint_audit_b3_multi_provider.py -q",
        "aria-env", ("env:aria-env",), evidence_kind="test_receipt",
    ),
    _lane(
        "c6_b3_real_e2e", ("claim_6",),
        "Repeated real-run ledger/claim/capsule equivalence on a frozen dataset.",
        "ARIA_B3_E2E_DATA=~/aria-data/pbmc3k_test conda run -n aria-env python "
        "-m pytest tests/test_preprint_audit_b3_multi_provider_e2e.py -q",
        "aria-env", ("env:aria-env", "data:pbmc3k", "manual:configured_provider"),
        evidence_kind="e2e_capsule",
    ),
    _lane(
        "c7_speculative_quarantine", ("claim_7",),
        "Proof that speculative text cannot enter claims, ledger export or capsules.",
        "conda run -n aria-env python -m pytest tests/test_hypothesis_quarantine.py "
        "tests/test_hypothesis_report_section.py tests/test_ledger_export.py "
        "tests/test_preprint_audit_c1_claim_boundary.py -q",
        "aria-env", ("env:aria-env",), evidence_kind="test_receipt",
    ),
)


def _resource_defaults(output_root: Path) -> dict[str, tuple[str, Path | None]]:
    home = Path.home()
    env_root = Path(os.environ.get("CONDA_ENVS_PATH", home / "anaconda3" / "envs"))
    return {
        "env:aria-env": ("conda-env://aria-env", env_root / "aria-env"),
        "env:aria-rna-env": ("conda-env://aria-rna-env", env_root / "aria-rna-env"),
        "env:aria-bench-env": ("conda-env://aria-bench-env", env_root / "aria-bench-env"),
        "env:aria-rnaseq-env": ("conda-env://aria-rnaseq-env", env_root / "aria-rnaseq-env"),
        "env:aria-chromatin-env": ("conda-env://aria-chromatin-env", env_root / "aria-chromatin-env"),
        "data:seqc_bgi": ("benchmark-data://seqc_maqc_BGI", home / ".aria/benchmarks/seqc_maqc_BGI/counts.tsv"),
        "data:seqc_multisite": ("benchmark-data://seqc_maqc_{BGI,CNL,MAY,AGR,NVS}", home / ".aria/benchmarks/seqc_maqc_NVS/counts.tsv"),
        "data:seqc_ercc": ("benchmark-data://seqc_maqc_BGI/ercc_truth.tsv", home / ".aria/benchmarks/seqc_maqc_BGI/ercc_truth.tsv"),
        "data:kang_muscat": ("benchmark-data://kang_muscat/clusters.json", home / ".aria/benchmarks/kang_muscat/clusters.json"),
        "data:h9_fastq": ("validation-data://h9_bulk_fastq", home / "Samael/H9-RNA/raw_fastq"),
        "data:pbmc3k": ("validation-data://pbmc3k_test", home / "aria-data/pbmc3k_test"),
        "data:gse278576": ("validation-data://gse278576_scatac_donor_h5mu", home / "Samael/Single/prueba/muon_processed/hc937_paired.h5mu"),
        "human:b2_annotators": ("freeze-input://human/b2_annotator_{1,2}.csv+adjudicator", output_root / "human/b2_adjudicator.csv"),
        "human:design_gold": ("freeze-input://human/design_gold.csv", output_root / "human/design_gold.csv"),
        "manual:configured_provider": ("manual-confirmation://configured-llm-provider", None),
    }


def probe_resources(
    output_root: str | Path,
    *,
    overrides: Mapping[str, bool | None] | None = None,
) -> dict[str, dict[str, Any]]:
    """Probe availability without serializing machine-absolute paths or secrets."""
    output_root = Path(output_root)
    overrides = dict(overrides or {})
    result: dict[str, dict[str, Any]] = {}
    for resource_id, (path_ref, path) in _resource_defaults(output_root).items():
        if resource_id in overrides:
            available = overrides[resource_id]
        elif path is None:
            available = None
        else:
            available = path.exists()
        state = (
            "available" if available is True else
            "missing" if available is False else
            "manual_confirmation_required"
        )
        result[resource_id] = {
            "available": available,
            "state": state,
            "path_ref": path_ref,
        }
    return result


def _receipt_status(
    lane: Mapping[str, Any], output_root: Path, provenance: Mapping[str, Any]
) -> str | None:
    receipt = output_root / "receipts" / f"{lane['lane_id']}.json"
    if not receipt.is_file():
        return None
    try:
        payload = json.loads(receipt.read_text(encoding="utf-8"))
    except Exception:
        return "invalid_receipt"
    if payload.get("status") != "pass":
        return "failed_receipt"
    if any(
        payload.get(key) != provenance.get(key)
        for key in ("source_snapshot_hash", "workflow_hash", "aria_version")
    ):
        return "stale_receipt"
    if (
        payload.get("schema_version") != "aria.preprint_freeze.receipt.v1"
        or payload.get("lane_id") != lane["lane_id"]
        or tuple(payload.get("claims") or ()) != tuple(lane["claims"])
        or payload.get("command") != lane.get("command")
        or payload.get("environment") != lane.get("environment")
        or payload.get("returncode") != 0
    ):
        return "invalid_receipt"

    artifact_entries = payload.get("artifacts")
    if not isinstance(artifact_entries, list):
        return "invalid_receipt"
    expected = tuple(lane.get("expected_artifacts") or ())
    if tuple(item.get("path") for item in artifact_entries) != expected:
        return "artifact_mismatch"
    resolved_root = output_root.resolve()
    for item in artifact_entries:
        relative = item.get("path")
        if not isinstance(relative, str):
            return "invalid_receipt"
        path = output_root / relative
        try:
            path.resolve().relative_to(resolved_root)
        except ValueError:
            return "invalid_receipt"
        if (
            not path.is_file()
            or item.get("size_bytes") != path.stat().st_size
            or item.get("sha256") != _sha256(path)
        ):
            return "artifact_mismatch"
    return "verified"


def _sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_tree_clean(repo_root: Path, output_root: Path) -> bool:
    """Check tracked/source cleanliness while allowing the dedicated output root."""
    try:
        rel = output_root.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        rel = ""
    cmd = ["git", "status", "--porcelain", "--", "."]
    if rel:
        cmd.append(f":(exclude){rel}")
    proc = subprocess.run(
        cmd, cwd=repo_root, capture_output=True, text=True, check=False,
    )
    return proc.returncode == 0 and not proc.stdout.strip()


def _source_snapshot_hash(repo_root: Path, output_root: Path) -> str:
    """Hash the indexed source snapshot while excluding freeze outputs."""
    try:
        rel = output_root.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        rel = ""
    cmd = ["git", "ls-files", "--stage", "-z", "--", "."]
    if rel:
        cmd.append(f":(exclude){rel}")
    proc = subprocess.run(cmd, cwd=repo_root, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError("could not hash the indexed source snapshot")
    return hashlib.sha256(proc.stdout).hexdigest()


def _safe_test_summary(stdout: str) -> str | None:
    for line in reversed(str(stdout or "").splitlines()):
        line = line.strip()
        if re.fullmatch(
            r"[.A-Za-z0-9 ,%()_-]*(?:passed|failed|skipped)"
            r"[.A-Za-z0-9 ,%()_-]*",
            line,
        ):
            return line[-300:]
    return None


def _sanitize_public_json_artifact(path: Path) -> None:
    """Remove machine-local Conda prefixes from a public JSON artifact.

    The environment name and lock hash remain in provenance. Other absolute
    paths are intentionally untouched so the repository-wide public-artifact
    guard can still expose unexpected leaks instead of hiding them.
    """
    payload = json.loads(path.read_text(encoding="utf-8"))
    changed = False

    def visit(value: Any) -> None:
        nonlocal changed
        if isinstance(value, dict):
            if value.get("conda_prefix") is not None:
                value["conda_prefix"] = None
                changed = True
            for nested in value.values():
                visit(nested)
        elif isinstance(value, list):
            for nested in value:
                visit(nested)

    visit(payload)
    if changed:
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )


def execute_lane(
    repo_root: str | Path,
    output_root: str | Path,
    lane_id: str,
    *,
    timeout: int = 14_400,
) -> dict[str, Any]:
    """Execute one inventory lane and atomically publish a current-commit receipt.

    Only statically registered, implemented lanes may run. A zero exit is not
    enough for artifact lanes: every declared output must exist and be hashed.
    The source tree must be clean aside from the dedicated freeze output root.
    """
    repo_root = Path(repo_root).resolve()
    output_root = Path(output_root)
    lane = next((dict(item) for item in LANES if item["lane_id"] == lane_id), None)
    if lane is None:
        raise KeyError(f"unknown preprint-freeze lane: {lane_id}")
    if lane["implementation"] != "available" or not lane.get("command"):
        raise RuntimeError(f"lane {lane_id} has no executable implementation")
    if not _source_tree_clean(repo_root, output_root):
        raise RuntimeError(
            "source tree is dirty outside docs/benchmark_results/preprint_v1"
        )

    version = collect_version_metadata(repo_root)
    inventory_path = output_root / "inventory.json"
    if not inventory_path.is_file():
        raise RuntimeError("generate inventory.json from a clean checkout first")
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    baseline = inventory.get("provenance") or {}
    if baseline.get("git_dirty") is not False:
        raise RuntimeError("inventory was not generated from a clean checkout")
    if baseline.get("git_commit") != version.get("git_commit"):
        raise RuntimeError("inventory commit does not match current HEAD")

    tokens = shlex.split(str(lane["command"]))
    env = os.environ.copy()
    while tokens and re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", tokens[0]):
        key, value = tokens.pop(0).split("=", 1)
        env[key] = os.path.expanduser(value)
    tokens = [os.path.expanduser(token) for token in tokens]
    started = datetime.now(timezone.utc).isoformat()
    proc = subprocess.run(
        tokens,
        cwd=repo_root,
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"lane {lane_id} failed with exit {proc.returncode}: "
            f"{proc.stderr[-2000:]}"
        )

    artifacts = []
    for relative in lane.get("expected_artifacts") or []:
        path = output_root / relative
        if not path.is_file():
            raise RuntimeError(
                f"lane {lane_id} returned zero but did not produce {relative}"
            )
        if path.suffix.lower() == ".json":
            _sanitize_public_json_artifact(path)
        artifacts.append({
            "path": relative,
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        })
    receipt = {
        "schema_version": "aria.preprint_freeze.receipt.v1",
        "lane_id": lane_id,
        "claims": lane["claims"],
        "status": "pass",
        "git_commit": baseline.get("git_commit"),
        "git_tree_sha": baseline.get("git_tree_sha"),
        "source_snapshot_hash": baseline.get("source_snapshot_hash"),
        "workflow_hash": baseline.get("workflow_hash"),
        "aria_version": baseline.get("aria_version"),
        "environment": lane["environment"],
        "command": lane["command"],
        "started_utc": started,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "returncode": proc.returncode,
        "stdout_sha256": hashlib.sha256(proc.stdout.encode("utf-8")).hexdigest(),
        "stderr_sha256": hashlib.sha256(proc.stderr.encode("utf-8")).hexdigest(),
        "test_summary": _safe_test_summary(proc.stdout),
        "artifacts": artifacts,
    }
    receipts = output_root / "receipts"
    receipts.mkdir(parents=True, exist_ok=True)
    write_inventory(receipt, receipts / f"{lane_id}.json")
    return receipt


def build_inventory(
    repo_root: str | Path,
    output_root: str | Path,
    *,
    resource_overrides: Mapping[str, bool | None] | None = None,
) -> dict[str, Any]:
    repo_root = Path(repo_root).resolve()
    output_root = Path(output_root)
    version = collect_version_metadata(repo_root)
    source_clean = _source_tree_clean(repo_root, output_root)
    source_snapshot_hash = _source_snapshot_hash(repo_root, output_root)
    git_describe = str(version.get("git_describe") or "")
    if source_clean and git_describe.endswith("-dirty"):
        git_describe = git_describe[:-6]
    freeze_workflow_payload = "\0".join((
        str(version.get("aria_version") or ""),
        source_snapshot_hash,
    ))
    provenance = {
        "aria_version": version.get("aria_version"),
        "version_source": version.get("version_source"),
        "git_commit": version.get("git_commit"),
        "git_tree_sha": version.get("git_tree_sha"),
        "source_snapshot_hash": source_snapshot_hash,
        "source_snapshot_hash_algorithm": (
            "sha256(git-ls-files-stage-excluding-output)"
        ),
        "git_describe": git_describe,
        "git_dirty": not source_clean,
        "workflow_hash": hashlib.sha256(
            freeze_workflow_payload.encode("utf-8")
        ).hexdigest(),
        "workflow_hash_algorithm": "sha256(version+source_snapshot_hash)",
    }
    resources = probe_resources(output_root, overrides=resource_overrides)
    lanes: list[dict[str, Any]] = []
    for spec in LANES:
        lane = dict(spec)
        receipt_state = _receipt_status(lane, output_root, provenance)
        resource_states = [resources[r]["available"] for r in lane["resources"]]
        if receipt_state is not None:
            status = receipt_state
        elif lane["implementation"] != "available":
            status = "blocked_missing_implementation"
        elif any(state is False for state in resource_states):
            status = "blocked_missing_resources"
        elif any(state is None for state in resource_states):
            status = "blocked_manual_confirmation"
        else:
            status = "ready_to_run"
        lane["status"] = status
        lanes.append(lane)

    claims = []
    for claim_id in CLAIM_IDS:
        claim_lanes = [lane for lane in lanes if claim_id in lane["claims"]]
        claims.append({
            "claim_id": claim_id,
            "lane_ids": [lane["lane_id"] for lane in claim_lanes],
            "status": (
                "verified" if claim_lanes and all(
                    lane["status"] == "verified" for lane in claim_lanes
                ) else "blocked"
            ),
        })

    blockers = [
        {"lane_id": lane["lane_id"], "status": lane["status"]}
        for lane in lanes if lane["required_for_freeze"]
        and lane["status"] != "verified"
    ]
    freeze_ready = (
        not provenance.get("git_dirty")
        and not blockers
        and all(claim["status"] == "verified" for claim in claims)
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "output_root": "docs/benchmark_results/preprint_v1",
        "provenance": provenance,
        "resources": resources,
        "claims": claims,
        "lanes": lanes,
        "freeze_gate": {
            "ready": freeze_ready,
            "n_required_lanes": sum(
                1 for lane in lanes if lane["required_for_freeze"]
            ),
            "n_verified_lanes": sum(
                1 for lane in lanes if lane["status"] == "verified"
            ),
            "blockers": blockers,
            "tag_action_authorized": False,
        },
    }


def write_inventory(payload: Mapping[str, Any], path: str | Path) -> Path:
    """Atomically publish a path-portable inventory JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
    if re.search(r'"/(?:home|Users|tmp)/', text):
        raise ValueError("inventory contains a machine-absolute path")
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)
    return path


__all__ = [
    "CLAIM_IDS", "LANES", "SCHEMA_VERSION", "build_inventory",
    "execute_lane", "probe_resources", "write_inventory",
]
