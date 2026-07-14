"""scRNAAgent pseudobulk condition DE mixin (A7 extraction from scrna_agent.py; bodies verbatim)."""
from __future__ import annotations

from aria.agents.scrna._base import *  # noqa: F401,F403


class PseudobulkMixin:
    PSEUDOBULK_KEYWORDS = (
        # Aging / time
        "aging", "age", "young", "old", "lifespan", "senescen",
        # Comparison / contrast verbs
        "compare", "comparison", "contrast", "between", "differential",
        "differentially expressed", "deg", "versus", " vs ", "vs.",
        # Disease vs healthy
        "disease", "disorder", "healthy", "control", "wildtype", "wild-type",
        # Treatment / perturbation
        "treatment", "treated", "untreated", "perturb", "stimulus",
        "knockout", "knock-out", "knockdown", "ko ", "wt ",
        # Common biological contrasts
        "tumor", "normal", "responder", "non-responder",
    )

    def _needs_pseudobulk(self, intent: dict, exp_ctx: dict) -> bool:
        """
        True when DesignAgent identified enough biological groups (the
        prerequisite for any between-condition DE) AND the user's intent
        carries comparison semantics. Conservative: returns False if either
        side is missing — we never run a between-condition test against a
        single group, and we don't surprise the user with extra compute
        when the question is purely descriptive (cell-type characterisation).

        v4.4 publication-readiness rule: n>=3 per group is recommended.
        n=2 remains supported with a low-power warning, but only when CP2
        selected optional supported analyses.
        """
        design = (exp_ctx or {}).get("design", {}) or {}
        groups = design.get("analysis_groups") or design.get("groups", {}) or {}
        pb_cfg = design.get("pseudobulk", {}) or {}
        has_obs_design = bool(
            pb_cfg.get("condition_col") and pb_cfg.get("replicate_col")
            and (len(groups) >= 2 or pb_cfg.get("comparisons"))
        )
        if len(groups) < 2:
            if not has_obs_design:
                return False
        if groups:
            if all(len(s or []) >= 3 for s in groups.values()):
                pass
            elif all(len(s or []) >= 2 for s in groups.values()):
                if not (exp_ctx or {}).get("run_optional_supported"):
                    return False
            else:
                return False
        # Honor the approved plan. When DesignIntelligence recommended
        # pseudobulk DE from an explicit obs design (condition + replicate +
        # groups/comparisons) and the user approved that plan at CP2, the
        # between-condition contrast is well-defined and must run regardless of
        # how the free-text question is phrased. Re-gating on question keywords
        # here previously OVERRODE the recommendation silently: an IFN-β
        # "response programs / signaling networks" question carried no keyword,
        # so DE was dropped without a logged decision or a reported skip even
        # though DesignIntelligence had recommended it (audit 2026-05-28
        # PBMC-blocker; F-ENG-E2E class of silent plan/dispatch disagreement).
        di = (exp_ctx or {}).get("design_intelligence", {}) or {}
        recommended = " ".join(di.get("recommended", []) or []).lower()
        if has_obs_design and "pseudobulk" in recommended:
            return True

        # Fallback for designs without an explicit DI recommendation (e.g.
        # filename-inferred groups): the question must carry comparison
        # semantics so we do not surprise the user with extra compute on a
        # purely descriptive cell-characterisation question.
        text = (
            (intent.get("summary", "") or "") + " "
            + " ".join(intent.get("biological_entities", []) or [])
            + " " + (exp_ctx.get("user_question", "") or "")
        ).lower()
        return any(kw in text for kw in self.PSEUDOBULK_KEYWORDS)

    def _inject_condition_obs(self, h5ad_in: str, design: dict,
                                output_path: str) -> dict:
        """
        Write a copy of `h5ad_in` with `obs[main_factor]` populated from
        the sample_id → group mapping in `design.groups`. The cell-level
        sample column is preferred from obs in this priority:
        `sample_id` (set by rna_concat) → `batch` → `orig.ident`.

        Returns:
            {"status": "success" | "skipped" | "error",
             "output_path": str | None,
             "condition_col": str | None,
             "replicate_col": str | None,
             "matched_cells": int, "unmatched_cells": int}
        """
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_inject_condition.py",
            params={
                "data_path":    h5ad_in,
                "groups":       design.get("groups", {}),
                "factor":       design.get("main_factor", "condition"),
                "batch_col":    design.get("batch_covariate"),
                "replicate_units": (
                    ((design.get("replicate_handling") or {}).get("sample_to_unit") or {})
                    if (design.get("replicate_handling") or {}).get("mode")
                    == "technical_aggregate"
                    else {}
                ),
                "output_path":  output_path,
            },
        )
        return result

    @staticmethod
    def _trusted_annotation_groupby(annotation: dict | None) -> str | None:
        annotation = annotation or {}
        trusted = annotation.get("trusted_groupby_for_inference")
        if trusted:
            return str(trusted)
        if annotation.get("annotation_source") == "input_obs":
            label_col = annotation.get("label_col")
            return str(label_col) if label_col else None
        celltypist = annotation.get("celltypist") or {}
        if celltypist.get("status") == "success":
            label_col = annotation.get("label_col") or celltypist.get("label_col")
            return str(label_col) if label_col else "cell_type_celltypist"
        return None

    @staticmethod
    def _annotation_is_report_only(annotation: dict | None) -> bool:
        annotation = annotation or {}
        if annotation.get("trusted_groupby_for_inference"):
            return False
        if annotation.get("annotation_source") in {
            "llm_marker_only",
            "unresolved_marker_fallback",
        }:
            return True
        celltypist = annotation.get("celltypist") or {}
        return bool(annotation.get("label_col")
                    and celltypist.get("status") != "success"
                    and annotation.get("annotation_for_report"))

    @staticmethod
    def _integration_qc_has_blocking(integration_qc: dict | None) -> bool:
        return any(
            isinstance(issue, dict) and issue.get("severity") == "blocking"
            for issue in ((integration_qc or {}).get("issues") or [])
        )

    @staticmethod
    def _integration_qc_blocking_reason(integration_qc: dict | None) -> str:
        for issue in ((integration_qc or {}).get("issues") or []):
            if isinstance(issue, dict) and issue.get("severity") == "blocking":
                msg = issue.get("message") or issue.get("check")
                rec = issue.get("recommendation")
                return " ".join(str(x) for x in (msg, rec) if x).strip()
        return "Integration QC marked the integrated embedding as blocking."

    def _run_pseudobulk(self, experiment_id: str,
                         current_h5ad: str,
                         exp_ctx: dict, intent: dict,
                         annotation: dict,
                         raw_counts_h5ad: str | None = None) -> dict:
        """
        Pseudobulk DE between condition groups identified by DesignAgent.

        Aggregates counts per (cell_type × replicate), fits a pyDESeq2
        model with design ~ condition [+ batch], and runs ORA per
        (cell_type × comparison) on the top-N DE genes.
        """
        from pathlib import Path
        design = (exp_ctx or {}).get("design", {}) or {}
        groups = design.get("groups", {}) or {}
        pb_cfg = design.get("pseudobulk", {}) or {}
        factor = pb_cfg.get("condition_col") or design.get("main_factor", "condition")
        batch_cov = design.get("batch_covariate")

        comparisons = self._normalise_pseudobulk_comparisons(
            pb_cfg.get("comparisons")
        )
        if not comparisons:
            suggestions = self._suggest_pseudobulk_comparisons(groups)
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        "Pseudobulk DE was not run because no explicit "
                        "test/reference comparison was confirmed."
                    ),
                    "suggested_comparisons": suggestions,
                },
                Confidence.INSUFFICIENT,
            )
            self.publish_escalation(
                experiment_id=experiment_id,
                checkpoint="scrna.pseudobulk.contrast",
                question=(
                    "Pseudobulk differential expression requires an explicit "
                    "comparison. Choose the test level and reference level "
                    "before ARIA runs DE."
                ),
                options=[
                    "Skip pseudobulk DE until a comparison is confirmed",
                    *[
                        f"Run {test} vs {ref} (reference={ref})"
                        for test, ref in suggestions[:6]
                    ],
                ],
                context={
                    "analysis_type": "scrna_pseudobulk_de",
                    "parameter_name": "comparison",
                    "suggested_comparisons": suggestions,
                },
            )
            return {
                "status": "skipped",
                "reason": "explicit_comparison_required",
                "suggested_comparisons": suggestions,
            }

        # 1. Prefer h5ad-native obs design when CP1 inferred one. Otherwise
        # inject condition obs from sample → group mapping.
        workspace = self._workspace(experiment_id, "pseudobulk")
        technical_units = (
            (design.get("replicate_handling") or {}).get("sample_to_unit") or {}
        )
        use_obs_design = bool(
            pb_cfg.get("from_obs") and factor and pb_cfg.get("replicate_col")
            and not technical_units
        )
        if use_obs_design:
            pb_input = current_h5ad
            replicate_col = pb_cfg.get("replicate_col")
        else:
            injected = workspace / "with_condition.h5ad"
            inj = self._inject_condition_obs(current_h5ad, design, str(injected))
            if inj.get("status") != "success":
                return {"status": "skipped",
                        "reason": f"inject_failed: {inj.get('reason', '?')}"}
            pb_input = str(injected)
            replicate_col = inj.get("replicate_col") or "sample_id"

        # 2. Choose groupby column. Report-only LLM/marker labels never define
        # donor-level pseudobulk inferential units.
        trusted_annotation_groupby = self._trusted_annotation_groupby(annotation)
        cell_type_col = pb_cfg.get("groupby_col")
        if (cell_type_col == "cell_type_celltypist"
                and (annotation or {}).get("celltypist", {}).get("status") != "success"
                and trusted_annotation_groupby):
            cell_type_col = trusted_annotation_groupby
        if not cell_type_col:
            cell_type_col = trusted_annotation_groupby or "leiden"

        if (cell_type_col == "leiden"
                and self._annotation_is_report_only(annotation)):
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        "Pseudobulk DE was not run because cell-type labels "
                        "were report-only LLM/marker annotations, not a "
                        "trusted inferential grouping."
                    ),
                    "annotation_source": (annotation or {}).get("annotation_source"),
                    "label_col": (annotation or {}).get("label_col"),
                },
                Confidence.INSUFFICIENT,
            )
            return {
                "status": "skipped",
                "reason": "llm_only_annotation_not_trusted_for_inference",
            }

        integration_qc = (exp_ctx or {}).get("integration_qc") or {}
        if self._integration_qc_has_blocking(integration_qc):
            reason = self._integration_qc_blocking_reason(integration_qc)
            self.publish_finding(
                experiment_id,
                {
                    "summary": (
                        "Pseudobulk DE and differential abundance were not run "
                        "because integration QC marked the embedding as "
                        "blocking for overcorrection."
                    ),
                    "reason": reason,
                    "integration_qc": integration_qc,
                    "groupby": cell_type_col,
                },
                Confidence.INSUFFICIENT,
            )
            return {
                "status": "skipped",
                "reason": "integration_overcorrection_blocking",
                "details": reason,
            }

        # Hard skip: when annotation failed for an unrecoverable matrix and
        # we would be left grouping by numeric Leiden IDs, pseudobulk between
        # conditions is not interpretable. Bail with a clear reason instead.
        if cell_type_col == "leiden" and (annotation or {}).get(
                "annotation_source") == "unrecoverable_matrix":
            return {"status": "skipped",
                    "reason": "no_celltype_groupby_after_annotation_failure"}
        covariates = list(pb_cfg.get("covariates") or [])
        if not covariates and batch_cov:
            covariates = [batch_cov]

        self._log_decision(
            experiment_id,
            checkpoint="scRNA",
            question="Pseudobulk design",
            decision=(
                f"groupby={cell_type_col}; condition={factor}; "
                f"replicate={replicate_col}; covariates="
                f"{', '.join(covariates) if covariates else 'none'}"
            ),
            rationale=(
                "Single-cell condition contrasts use donor/sample-level "
                "pseudobulk rather than treating cells as independent replicates."
            ),
            made_by="scrna_agent",
        )

        # 2b. (T1.1) Differential abundance BEFORE pseudobulk DE.
        # Two purposes: (a) report shifts in cell-type proportions as a
        # primary observation; (b) decide whether the pseudobulk design
        # should include a composition covariate to control for those
        # shifts when comparing within-cell-type expression.
        da_result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_diff_abundance.py",
            params={
                "data_path":          str(pb_input),
                "groupby":            cell_type_col,
                "condition_col":      factor,
                "replicate_col":      replicate_col,
                "comparisons":        comparisons,
                "covariates":         covariates,
                "significance_alpha": 0.10,
                "output_dir":         str(workspace),
            },
        )
        da_significant = bool(
            da_result.get("status") == "success"
            and da_result.get("any_significant")
        )
        if da_result.get("status") == "success":
            self._log_decision(
                experiment_id,
                checkpoint="scRNA",
                question="Compositional shift before pseudobulk DE",
                decision=(
                    "composition_covariate=ON"
                    if da_significant else "composition_covariate=OFF"
                ),
                rationale=(
                    f"rna_diff_abundance flagged at least one cell type as "
                    f"significantly shifting (alpha=0.10); adding "
                    f"log(cells_in_group/total_cells) as a covariate in the "
                    f"DESeq2 design so per-cell-type DE is not confounded by "
                    f"mixture changes."
                    if da_significant else
                    "No cell type shifted significantly (alpha=0.10); the "
                    "DESeq2 design stays free of a composition covariate."
                ),
                made_by="scrna_agent",
            )

        group_sizes = [
            len(samples or []) for samples in (groups or {}).values()
        ]
        low_power_optional = (
            bool(group_sizes)
            and min(group_sizes) == 2
            and bool((exp_ctx or {}).get("run_optional_supported"))
        )
        min_replicates = 2 if low_power_optional else 3

        # P0-7: propagate the user-confirmed CP3 thresholds instead of
        # hardcoding padj/lfc here. The orchestrator records them in
        # exp_ctx["global_padj"]/["global_lfc"]; AnalysisThresholds resolves them
        # (scRNA pseudobulk keeps lfc default 0.5 when none was confirmed).
        from aria.utils.thresholds import AnalysisThresholds
        thresholds = AnalysisThresholds.from_exp_context(
            exp_ctx, modality="scRNA",
            log2fc_default=0.5, min_cells=10, min_replicates=min_replicates,
        )

        # 3. Run pseudobulk DE
        pb = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_pseudobulk_de.py",
            params={
                "data_path":     str(pb_input),
                **({"counts_data_path": str(raw_counts_h5ad)}
                   if raw_counts_h5ad else {}),
                "groupby":       cell_type_col,
                "condition_col": factor,
                "replicate_col": replicate_col,
                "comparisons":   comparisons,
                "covariates":    covariates,
                "composition_covariate": da_significant,
                **thresholds.as_pseudobulk_params(),
                "top_n":         50,
                "auto_paired_donor_covariate": True,
                "output_dir":    str(workspace),
            },
        )
        if pb.get("status") != "success":
            return {"status": "error",
                    "reason": f"pseudobulk_de_failed: {pb.get('error_type', '?')}",
                    "details": pb.get("details", "")[:300]}

        # 4. Pathway ORA per (group × comparison) — reuse the per-cluster
        # ORA script with synthetic cluster IDs like "Astro::old_vs_young".
        de_for_ora: dict = {}
        # C2 (audit 2026-05-29): give each per-cluster ORA its own universe —
        # the genes tested in that cell type's pseudobulk — instead of the
        # whole-dataset expressed set, which inflates per-cluster enrichment.
        background_by_cluster: dict = {}
        for group, info in (pb.get("per_group") or {}).items():
            if info.get("status") == "skipped":
                continue
            grp_bg = info.get("background_genes") or []
            for comp_key, comp in (info.get("per_comparison") or {}).items():
                if comp.get("status") != "success" or not comp.get("all_sig"):
                    continue
                cluster_id = f"{group}::{comp_key}"
                de_for_ora[cluster_id] = [
                    {"gene": r["gene"], "log2fc": r["log2fc"], "padj": r["padj"]}
                    for r in comp["all_sig"]
                ]
                if grp_bg:
                    background_by_cluster[cluster_id] = grp_bg

        pw_findings = None
        if de_for_ora:
            pw = self.env.run_in_stack(
                stack="rna",
                script_path="aria/scripts/rna_pathway_per_cluster.py",
                params={
                    "de_genes_by_cluster":   de_for_ora,
                    "organism":              exp_ctx.get("organism",
                                                          "Homo sapiens"),
                    "top_genes_per_cluster": 200,
                    "background_genes":       pb.get("background_genes", []),
                    "background_genes_by_cluster": background_by_cluster,
                    "padj_db_max":           0.05,
                    "output_dir":            str(workspace),
                },
            )
            if pw.get("status") == "success":
                pw_findings = {
                    "organism":    pw.get("organism"),
                    "databases":   pw.get("databases", {}),
                    "background_size": pw.get("background_size"),
                    "background_source": pw.get("background_source"),
                    "params_sha256": pw.get("params_sha256"),
                    "per_cluster": pw.get("per_cluster", {}),
                }

        # 5. Build the finding payload + announce on the bus
        pb_payload = {
            "groupby":       pb.get("groupby"),
            "condition_col": pb.get("condition_col"),
            "replicate_col": pb.get("replicate_col"),
            "covariates":    pb.get("covariates", []),
            "from_obs":      use_obs_design,
            "thresholds":    pb.get("thresholds", {}),
            "multiple_testing": pb.get("multiple_testing", {}),
            "background_size": pb.get("background_size"),
            "background_source": pb.get("background_source"),
            "count_source": pb.get("count_source"),
            "count_source_data_path": pb.get("count_source_data_path"),
            "lognorm_recovered": pb.get("lognorm_recovered"),
            "params_sha256": pb.get("params_sha256"),
            "differential_abundance_params_sha256": da_result.get("params_sha256"),
            "paired_design": pb.get("paired_design"),
            "auto_paired_donor_covariate": pb.get("auto_paired_donor_covariate"),
            "n_groups":      pb.get("n_groups"),
            "per_group":     pb.get("per_group", {}),
        }
        # Summary line
        per_group = pb.get("per_group", {}) or {}
        n_with_de = sum(
            1 for g in per_group.values()
            for c in (g.get("per_comparison", {}) or {}).values()
            if c.get("status") == "success" and c.get("n_significant", 0) > 0
        )
        comp_str = ", ".join(f"{t}_vs_{r}" for t, r in comparisons)
        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Pseudobulk DE (DESeq2) on {pb.get('n_groups', 0)} cell "
                f"types across {comp_str}: {n_with_de} (group × comparison) "
                f"blocks with significant DE."),
             "comparisons": comparisons,
             "n_groups":    pb.get("n_groups"),
             "output_csv":  pb.get("output_csv")},
            Confidence.HIGH if n_with_de > 0 else Confidence.INSUFFICIENT,
        )

        da_payload = None
        if da_result.get("status") == "success":
            da_payload = {
                "method":                 da_result.get("method"),
                "groupby":                da_result.get("groupby"),
                "condition_col":          da_result.get("condition_col"),
                "replicate_col":          da_result.get("replicate_col"),
                "covariates":             da_result.get("covariates", []),
                "paired_design":          da_result.get("paired_design"),
                "significance_alpha":     da_result.get("significance_alpha"),
                "any_significant":        da_result.get("any_significant"),
                "n_replicates_per_group": da_result.get("n_replicates_per_group"),
                "per_comparison":         da_result.get("per_comparison", {}),
                "output_path":            da_result.get("output_path"),
                "params_sha256":          da_result.get("params_sha256"),
                "warnings":               da_result.get("warnings", []),
            }
        elif da_result.get("status"):
            da_payload = {
                "status":     da_result.get("status"),
                "error_type": da_result.get("error_type"),
                "details":    (da_result.get("details") or "")[:300],
            }

        return {
            "status":                 "success",
            "differential_abundance": da_payload,
            "pseudobulk_de":          pb_payload,
            "pseudobulk_pathways":    pw_findings,
        }

    @staticmethod
    def _normalise_pseudobulk_comparisons(raw) -> list[list[str]]:
        if not raw:
            return []
        comparisons = []
        for comp in raw:
            if isinstance(comp, dict):
                test = comp.get("test") or comp.get("case") or comp.get("contrast")
                ref = comp.get("ref") or comp.get("reference") or comp.get("control")
                if test and ref:
                    comparisons.append([str(test), str(ref)])
            elif isinstance(comp, (list, tuple)) and len(comp) >= 2:
                comparisons.append([str(comp[0]), str(comp[1])])
        return comparisons

    @staticmethod
    def _suggest_pseudobulk_comparisons(groups: dict) -> list[list[str]]:
        names = sorted(str(g) for g in (groups or {}).keys())
        if len(names) < 2:
            return []
        lower = {g.lower(): g for g in names}
        for ref_key in ("control", "ctrl", "wt", "wildtype", "healthy", "untreated", "baseline"):
            if ref_key in lower:
                ref = lower[ref_key]
                return [[g, ref] for g in names if g != ref]
        return [[test, ref] for i, ref in enumerate(names) for test in names[i + 1:]]

    # ── Trajectory ────────────────────────────────────────────────────────

