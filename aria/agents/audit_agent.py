"""
ARIA AuditAgent (v4.1) — Quality Linter
-----------------------------------------
Runs synchronously after experimental design is confirmed, before any
computationally expensive dispatch (STAR, DESeq2, etc.).

Philosophy: "Inspect before investing."
Only flags issues where the finding would change the user's decision.

Three checks in v4.1 (bulk RNA-seq focus):
  1. Replicate correlation  — detect outliers and sample swaps in count matrix
  2. PCA batch dominance    — flag when PC1 is driven by batch, not condition
  3. STAR alignment rate    — flag low unique-mapping % from STAR Log.final.out

Each finding has a severity:
  "blocking" → requires user acknowledgement before proceeding (CP 3.5)
  "warning"  → logged and surfaced in report, analysis proceeds automatically
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Optional

from aria.agents.base_agent import BaseAgent
from aria.bus.message_bus import Confidence

log = logging.getLogger("aria.audit")

# Thresholds
_CORR_WARN_THR   = 0.85   # within-group Pearson r below this → outlier warning
_ALIGN_WARN_THR  = 70.0   # STAR unique-mapping % below this → alignment warning
_PC1_BATCH_THR   = 0.50   # PC1 variance fraction above this → check batch dominance


class AuditAgent(BaseAgent):
    name        = "audit_agent"
    description = "Quality linter: inspects data before dispatch."

    def run(self, experiment_id: str, context: dict) -> dict:
        return self.run_audit(context.get("exp_context", context), experiment_id)

    # ── Public entry point ────────────────────────────────────────────────

    def run_audit(self, exp_context: dict, experiment_id: str) -> dict:
        """
        Run all applicable checks.

        Returns:
            {
              "status":   "clean" | "warnings" | "blocking",
              "findings": [
                  {
                    "severity":       "blocking" | "warning",
                    "check":          str,
                    "message":        str,   # user-facing, specific
                    "recommendation": str,
                  },
                  ...
              ]
            }
        """
        self.publish_status(experiment_id, "Running quality audit...", 0.05)
        findings = []

        files  = exp_context.get("files", [])
        design = exp_context.get("design", {})

        # ── 1. Replicate correlation ──────────────────────────────────────
        try:
            corr_findings = self._check_replicate_correlation(files, design)
            findings.extend(corr_findings)
        except Exception as e:
            log.warning(f"Replicate correlation check failed: {e}")

        # ── 2. PCA batch dominance ────────────────────────────────────────
        if design.get("batch_factor"):
            try:
                batch_findings = self._check_pca_batch_dominance(files, design)
                findings.extend(batch_findings)
            except Exception as e:
                log.warning(f"PCA batch dominance check failed: {e}")

        # ── 3. STAR alignment rate ────────────────────────────────────────
        try:
            star_findings = self._check_star_alignment(files)
            findings.extend(star_findings)
        except Exception as e:
            log.warning(f"STAR alignment check failed: {e}")

        has_blocking = any(f["severity"] == "blocking" for f in findings)
        has_warnings = any(f["severity"] == "warning" for f in findings)

        status = "blocking" if has_blocking else ("warnings" if has_warnings else "clean")

        for f in findings:
            log.warning(f"[audit:{f['severity']}] {f['check']}: {f['message']}")

        if not findings:
            self.publish_status(experiment_id, "Quality audit passed — no issues found.", 0.08)
        else:
            n_block = sum(1 for f in findings if f["severity"] == "blocking")
            n_warn  = sum(1 for f in findings if f["severity"] == "warning")
            self.publish_status(
                experiment_id,
                f"Quality audit: {n_block} blocking issue(s), {n_warn} warning(s).",
                0.08,
            )

        return {"status": status, "findings": findings}

    # ── Check 1: Replicate correlation ────────────────────────────────────

    def _check_replicate_correlation(self, files: list, design: dict) -> list:
        """
        Loads the count matrix, computes log1p(CPM) Pearson correlations, and
        compares within-group vs cross-group similarity for every sample.

        Flags:
          blocking — any sample correlates more with another group than its own
          warning  — any sample has within-group r < 0.85
        """
        if not files:
            return []

        try:
            import numpy as np
            import pandas as pd
        except ImportError:
            return []

        counts = _load_count_matrix(files)
        if counts is None or counts.empty:
            return []

        group_labels = _infer_group_labels(counts, design)
        if not group_labels or len(set(group_labels.values())) < 2:
            return []

        # Filter to expressed genes and normalize to log1p(CPM)
        expressed = counts.loc[counts.sum(axis=1) >= 10]
        if expressed.empty:
            return []

        cpm   = expressed.div(expressed.sum(axis=0) / 1e6)
        lcpm  = np.log1p(cpm)

        corr = lcpm.corr(method="pearson")

        findings = []
        samples  = list(group_labels.keys())

        for sample in samples:
            grp = group_labels[sample]
            same_grp   = [s for s in samples if group_labels[s] == grp and s != sample]
            other_grp  = [s for s in samples if group_labels[s] != grp]

            if not same_grp:
                continue

            within_r = corr.loc[sample, same_grp].mean()
            cross_r  = corr.loc[sample, other_grp].mean() if other_grp else 0.0

            if within_r < cross_r:
                # Sample correlates more with the other group → possible swap
                best_other = corr.loc[sample, other_grp].idxmax()
                findings.append({
                    "severity": "blocking",
                    "check":    "replicate_correlation",
                    "message":  (
                        f"Sample '{sample}' ({grp}) correlates more strongly with "
                        f"'{best_other}' ({group_labels[best_other]}) "
                        f"(r={corr.loc[sample, best_other]:.3f}) than with its own "
                        f"replicates (mean r={within_r:.3f}). Possible sample swap "
                        f"or mislabeling."
                    ),
                    "recommendation": (
                        "Verify sample identity. If confirmed swap, correct the design "
                        "labels before proceeding. If uncertain, exclude the sample."
                    ),
                })
            elif within_r < _CORR_WARN_THR:
                findings.append({
                    "severity": "warning",
                    "check":    "replicate_correlation",
                    "message":  (
                        f"Sample '{sample}' ({grp}) has low mean correlation with its "
                        f"replicates (r={within_r:.3f}, threshold {_CORR_WARN_THR}). "
                        f"Possible outlier or technical artifact."
                    ),
                    "recommendation": (
                        f"Inspect '{sample}' for library size anomalies or "
                        f"unexpected gene expression patterns. Consider excluding it "
                        f"if the outlier status cannot be explained biologically."
                    ),
                })

        return findings

    # ── Check 2: PCA batch dominance ─────────────────────────────────────

    def _check_pca_batch_dominance(self, files: list, design: dict) -> list:
        """
        Runs a quick PCA on log1p(CPM) counts and checks whether PC1 separates
        samples by batch more than by biological condition.

        Only runs when the design explicitly has a batch_factor.
        """
        if not files:
            return []

        try:
            import numpy as np
            import pandas as pd
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import LabelEncoder
        except ImportError:
            return []

        counts = _load_count_matrix(files)
        if counts is None or counts.empty:
            return []

        group_labels = _infer_group_labels(counts, design)
        if not group_labels:
            return []

        batch_labels = _infer_batch_labels(counts, design)
        if not batch_labels or len(set(batch_labels.values())) < 2:
            return []

        expressed = counts.loc[counts.sum(axis=1) >= 10]
        if expressed.shape[0] < 10:
            return []

        cpm  = expressed.div(expressed.sum(axis=0) / 1e6)
        lcpm = np.log1p(cpm).T.values  # shape: (n_samples, n_genes)

        pca = PCA(n_components=min(2, lcpm.shape[0] - 1))
        coords = pca.fit_transform(lcpm)
        pc1    = coords[:, 0]
        var1   = pca.explained_variance_ratio_[0]

        if var1 < _PC1_BATCH_THR:
            return []  # PC1 is not dominant enough to matter

        samples = [s for s in counts.columns if s in group_labels]
        cond_vals  = [group_labels.get(s, "?") for s in samples]
        batch_vals = [batch_labels.get(s, "?") for s in samples]

        # Variance of PC1 means per label (higher = better separation)
        def _between_var(labels: list, values) -> float:
            from collections import defaultdict
            groups: dict = defaultdict(list)
            for lbl, v in zip(labels, values):
                groups[lbl].append(v)
            grand_mean = np.mean(values)
            return float(np.mean([
                len(v) * (np.mean(v) - grand_mean) ** 2
                for v in groups.values()
            ]))

        cond_bvar  = _between_var(cond_vals,  pc1)
        batch_bvar = _between_var(batch_vals, pc1)

        if batch_bvar > cond_bvar:
            return [{
                "severity": "warning",
                "check":    "pca_batch_dominance",
                "message":  (
                    f"PC1 ({var1*100:.1f}% of variance) separates samples more strongly "
                    f"by batch than by experimental condition "
                    f"(batch spread={batch_bvar:.2f}, condition spread={cond_bvar:.2f}). "
                    f"Batch effect may dominate the analysis."
                ),
                "recommendation": (
                    "Consider including the batch variable as a covariate in the DESeq2 "
                    "design formula (e.g., ~ batch + condition). ARIA will propagate "
                    "this automatically if the batch factor is confirmed in the design."
                ),
            }]

        return []

    # ── Check 3: STAR alignment rate ─────────────────────────────────────

    def _check_star_alignment(self, files: list) -> list:
        """
        Searches for STAR Log.final.out files adjacent to the count matrix files.
        Parses 'Uniquely mapped reads %' and flags samples below threshold.
        """
        if not files:
            return []

        star_logs = _find_star_logs(files)
        if not star_logs:
            return []  # Pre-processed data — no logs available, skip silently

        findings = []
        for sample_name, log_path in star_logs.items():
            pct = _parse_star_unique_pct(log_path)
            if pct is None:
                continue
            if pct < _ALIGN_WARN_THR:
                findings.append({
                    "severity": "warning",
                    "check":    "star_alignment_rate",
                    "message":  (
                        f"Sample '{sample_name}' has a low unique mapping rate "
                        f"({pct:.1f}%, threshold {_ALIGN_WARN_THR}%). This may "
                        f"indicate the wrong reference genome, rRNA contamination, "
                        f"or poor library quality."
                    ),
                    "recommendation": (
                        "Verify that the reference genome matches the sample organism "
                        "and assembly. Check for rRNA or adapter contamination in the "
                        "raw reads. Consider excluding the sample if mapping cannot "
                        "be improved."
                    ),
                })

        return findings


# ── Module-level helpers (no bus, pure computation) ───────────────────────────

def _load_count_matrix(files: list):
    """Load and merge count files into a single (genes × samples) DataFrame."""
    try:
        import pandas as pd

        frames = []
        for f in files:
            if not f or not Path(f).exists():
                continue
            sep = "\t" if str(f).endswith((".tsv", ".txt")) else ","
            df  = pd.read_csv(f, sep=sep, index_col=0, comment="#")
            # Drop non-numeric columns (featureCounts metadata)
            df  = df.select_dtypes(include="number")
            frames.append(df)

        if not frames:
            return None

        if len(frames) == 1:
            return frames[0]

        # Multiple files: each file is one sample (one numeric column)
        combined = frames[0].copy()
        for frame in frames[1:]:
            for col in frame.columns:
                if col not in combined.columns:
                    combined[col] = frame[col]
        return combined

    except Exception as e:
        log.warning(f"_load_count_matrix failed: {e}")
        return None


def _infer_group_labels(counts, design: dict) -> dict:
    """
    Map column names → group labels using design["groups"].
    Returns {column: group} for columns that can be matched.
    """
    groups = design.get("groups", {})  # {group_name: [sample_stem, ...]}
    if not groups:
        return {}

    sample_to_group = {
        sample: grp
        for grp, samples in groups.items()
        for sample in samples
    }

    result = {}
    for col in counts.columns:
        for stem, grp in sample_to_group.items():
            if stem in col:
                result[col] = grp
                break

    return result


def _infer_batch_labels(counts, design: dict) -> dict:
    """
    Map column names → batch labels using design["batch_map"] or
    design["batches"] if present.
    """
    batch_map = design.get("batch_map", {})   # {sample_stem: batch_label}
    if not batch_map:
        return {}

    result = {}
    for col in counts.columns:
        for stem, batch in batch_map.items():
            if stem in col:
                result[col] = str(batch)
                break

    return result


def _find_star_logs(files: list) -> dict:
    """
    Walk directories of the provided files looking for STAR Log.final.out files.
    Returns {sample_name: Path}.
    """
    logs = {}
    checked_dirs = set()

    for f in files:
        p = Path(f)
        for search_dir in [p.parent, p.parent.parent]:
            if search_dir in checked_dirs:
                continue
            checked_dirs.add(search_dir)
            for log_file in search_dir.rglob("*Log.final.out"):
                # Use the parent directory name as sample name
                sample = log_file.parent.name
                logs[sample] = log_file

    return logs


def _parse_star_unique_pct(log_path: Path) -> Optional[float]:
    """Parse 'Uniquely mapped reads %' from a STAR Log.final.out file."""
    try:
        text = Path(log_path).read_text()
        m = re.search(r"Uniquely mapped reads %\s*\|\s*([\d.]+)%", text)
        return float(m.group(1)) if m else None
    except Exception:
        return None
