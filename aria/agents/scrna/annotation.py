"""scRNAAgent cell-type annotation mixin (A7 extraction from scrna_agent.py; bodies verbatim)."""
from __future__ import annotations

from aria.agents.scrna._base import *  # noqa: F401,F403


class AnnotationMixin:
    @staticmethod
    def _infer_tissue_hint(exp_ctx: dict, intent: dict) -> str | None:
        """Use only an explicit tissue/model hint; do not infer one from prose."""
        for source in (exp_ctx or {}, intent or {}):
            for key in ("celltypist_tissue_hint", "tissue_hint", "tissue"):
                value = source.get(key)
                if isinstance(value, str) and value.strip():
                    return value.strip().lower()
        return None

    @staticmethod
    def _allow_default_immune_model(exp_ctx: dict, intent: dict) -> bool:
        """Return True only when the user/config explicitly accepts fallback."""
        for source in (exp_ctx or {}, intent or {}):
            value = source.get("allow_default_immune_model")
            if isinstance(value, bool):
                return value
            if isinstance(value, str) and value.strip().lower() in {
                "1", "true", "yes", "on",
            }:
                return True
        return False

    def _predefined_celltype_col(self, h5ad_path: str,
                                  exp_ctx: dict) -> str | None:
        """
        Return a cell-type column name to reuse from obs, or None.

        DataAuditAgent populates exp_ctx['design']['pseudobulk']['groupby_col']
        when the input h5ad shipped with usable cell-type annotation (Seurat
        `subclass`, `cell_type_celltypist`, etc.). We trust that audit decision
        but also verify the column survives the current pipeline state (post-
        QC, post-concat) before promoting it to the canonical grouping.
        """
        design = (exp_ctx or {}).get("design", {}) or {}
        pb_cfg = design.get("pseudobulk", {}) or {}
        col = pb_cfg.get("groupby_col")
        if not col:
            inferred = (exp_ctx or {}).get("inferred_design", {}) or {}
            col = inferred.get("groupby_col")
        if not col or col == "leiden":
            return None
        # Verify the column survives in the working h5ad.
        try:
            import anndata as ad
            adata = ad.read_h5ad(h5ad_path, backed="r")
            try:
                if col not in adata.obs.columns:
                    return None
                vals = adata.obs[col].astype(str)
                levels = [v for v in vals.unique() if v and v.lower() != "nan"]
                if len(levels) < 2:
                    return None
            finally:
                backing = getattr(adata, "file", None)
                if backing is not None:
                    backing.close()
        except Exception as e:
            log.debug(f"Predefined celltype check failed for {col}: {e}")
            return None
        return col

    def _annotation_from_obs(self, experiment_id: str,
                              clustered_h5ad: str,
                              cell_type_col: str,
                              cluster_sizes: dict,
                              top_markers: dict) -> dict:
        """
        Build a findings.cell_types payload from a pre-existing obs column,
        bypassing CellTypist and the LLM reinterpretation. Each unique label
        becomes one entry whose `cluster_id` is the label itself (downstream
        groupby is the same column, so cluster IDs and labels coincide).
        """
        cell_types: dict = {}
        for label, n in cluster_sizes.items():
            cell_types[str(label)] = {
                "cell_type":         str(label),
                "celltypist_label":  None,
                "agrees_with_celltypist": None,
                "confidence":        "high",
                "rationale":         (
                    "Pre-existing annotation from input obs column "
                    f"'{cell_type_col}'. ARIA reused the user-supplied "
                    "labels and did not re-cluster or re-annotate."
                ),
                "key_markers":       top_markers.get(str(label), [])[:5],
                "n_cells":           int(n),
                "annotation_source": "input_obs",
            }

        labels_preview = [v["cell_type"] for v in list(cell_types.values())[:5]]
        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Reused {len(cell_types)} cell types from input "
                f"obs['{cell_type_col}']: {labels_preview}"
             ),
             "cell_types":   cell_types,
             "celltypist":   {"ran": False, "reason": "predefined_obs"},
             "markers_used": {k: v[:5] for k, v in (top_markers or {}).items()}},
            Confidence.HIGH,
        )
        return {
            "cell_types":        cell_types,
            "markers_used":      top_markers or {},
            "n_clusters":        len(cell_types),
            "celltypist":        {"ran": False, "reason": "predefined_obs"},
            "tissue_hint":       None,
            "label_col":         cell_type_col,
            "annotation_for_report": cell_type_col,
            "trusted_groupby_for_inference": cell_type_col,
            "annotated_h5ad":    clustered_h5ad,
            "annotation_source": "input_obs",
        }

    def _annotate_cell_types(self, experiment_id: str,
                              clustered_h5ad: str,
                              top_markers: dict,
                              exp_ctx: dict, intent: dict) -> dict:
        if not top_markers:
            return {"cell_types": {}, "markers_used": {}, "n_clusters": 0,
                    "celltypist": {"ran": False}}

        # ── Layer 1: CellTypist (database-backed, code-guarantee) ────────
        tissue_hint = self._infer_tissue_hint(exp_ctx, intent)
        allow_default_immune_model = self._allow_default_immune_model(
            exp_ctx, intent
        )
        celltypist_result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_celltypist.py",
            params={
                "data_path":    clustered_h5ad,
                "organism":     exp_ctx.get("organism", "Homo sapiens"),
                "tissue_hint":  tissue_hint,
                "allow_default_immune_model": allow_default_immune_model,
                "cluster_col":  "leiden",
                "majority_voting": True,
                "output_dir":   str(self._workspace(experiment_id, "annotation")),
            },
        )

        # Trim markers for the LLM prompt (independent of CellTypist status).
        markers_for_prompt = {
            str(k): [g for g in (v or [])[:20] if g and str(g) != "nan"]
            for k, v in top_markers.items()
        }

        # ── Hard stop: when CellTypist rejects the matrix because it is
        # scaled (Seurat-style export with no recoverable log1p-CPM and no
        # pre-existing cell-type labels in obs), do NOT fall back to LLM-only
        # annotation on noisy markers. Surface the failure as an explicit
        # warning so the report does not silently report `annotation_failed`.
        unrecoverable = (
            celltypist_result.get("status") == "error"
            and celltypist_result.get("error_type") in {
                "InvalidExpressionMatrix",
                "NoUsableCountsOrLogNorm",
            }
        )
        if unrecoverable:
            reason = celltypist_result.get("details", "")[:300]
            log.warning(
                f"CellTypist rejected matrix: {celltypist_result.get('error_type')}. "
                "ARIA will not synthesise cell-type labels from noisy markers."
            )
            self.publish_finding(
                experiment_id,
                {"summary": (
                    "Cell-type annotation skipped: input expression matrix is "
                    "scaled/normalised in a way CellTypist cannot validate, "
                    "and no recoverable log1p-CPM or counts were found. "
                    "Provide raw counts or include cell-type labels in obs."
                 ),
                 "cell_types":   {},
                 "celltypist":   {
                     "ran":         False,
                     "error_type":  celltypist_result.get("error_type"),
                     "details":     reason,
                     "tissue_hint": tissue_hint,
                 },
                 "markers_used": {k: v[:5] for k, v in markers_for_prompt.items()}},
                Confidence.LOW,
            )
            return {
                "cell_types":     {},
                "markers_used":   markers_for_prompt,
                "n_clusters":     len(markers_for_prompt),
                "celltypist":     celltypist_result,
                "tissue_hint":    tissue_hint,
                "label_col":      None,
                "annotation_for_report": None,
                "trusted_groupby_for_inference": None,
                "annotated_h5ad": None,
                "annotation_source": "unrecoverable_matrix",
                "warnings":       [
                    f"Cell-type annotation unavailable: "
                    f"{celltypist_result.get('error_type')}. {reason}"
                ],
            }

        # ── Layer 2: LLM reinterprets CellTypist results in biological
        #            context — does NOT invent labels from markers alone ──
        if celltypist_result.get("status") == "success":
            per_cluster = celltypist_result.get("per_cluster", {})
            celltypist_evidence = per_cluster
            prompt = build_untrusted_prompt(
                task="""Reinterpret database-backed CellTypist calls in biological
context. `frequency` is the fraction of raw per-cell calls matching the cluster
label; `mean_confidence` is the model's mean per-cell probability;
`alt_labels` are competing raw calls.

Rules:
- Default to the CellTypist label. Override only if markers contradict it and
  the rationale names the specific mismatch.
- HIGH: frequency >=0.85 and confidence null/>=0.7 with consistent markers.
- MEDIUM: frequency 0.5-0.85, confidence 0.5-0.7, or partial marker support.
- LOW: frequency <0.5, confidence <0.5, or marker contradiction.
- Never invent a label unsupported by CellTypist or the supplied markers.""",
                fields=[
                    PromptDataField("celltypist_model", celltypist_result.get(
                        "model_used"
                    ), "identifier", "celltypist"),
                    PromptDataField("organism", exp_ctx.get("organism", "?"),
                                    "metadata", "data_audit"),
                    PromptDataField("biological_question", intent.get(
                        "summary", exp_ctx.get("user_question", "?")
                    ), "user_text", "user"),
                    PromptDataField("tissue_hint", tissue_hint, "metadata",
                                    "user_or_data_audit"),
                    PromptDataField("celltypist_per_cluster", per_cluster,
                                    "structured_result", "celltypist"),
                    PromptDataField("marker_genes_by_cluster", markers_for_prompt,
                                    "structured_result", "rna_markers"),
                ],
                response_contract='''Return JSON only:
{
  "cluster_id": {
    "cell_type": "chosen final label",
    "celltypist_label": "database label",
    "agrees_with_celltypist": true,
    "confidence": "high|medium|low",
    "rationale": "one sentence",
    "key_markers": ["gene1", "gene2"]
  }
}''',
            )
        else:
            log.warning(
                f"CellTypist failed ({celltypist_result.get('error_type', '?')}); "
                f"falling back to LLM-only annotation."
            )
            celltypist_evidence = None
            prompt = build_untrusted_prompt(
                task=(
                    "CellTypist was unavailable. Annotate clusters from marker "
                    "genes alone, conservatively. If markers are ambiguous, say "
                    "so and list possible types. Maximum confidence is MEDIUM."
                ),
                fields=[
                    PromptDataField("organism", exp_ctx.get(
                        "organism", "unknown"
                    ), "metadata", "data_audit"),
                    PromptDataField("biological_question", intent.get(
                        "summary", "unknown"
                    ), "user_text", "user"),
                    PromptDataField("marker_genes_by_cluster", markers_for_prompt,
                                    "structured_result", "rna_markers"),
                ],
                response_contract=(
                    'Return JSON only: {"cluster_id": {"cell_type": "name", '
                    '"confidence": "medium|low", "key_markers": ["gene1"], '
                    '"rationale": "one sentence"}}'
                ),
            )

        cell_types: dict = {}
        try:
            raw = self.llm.complete(
                prompt=prompt,
                system=system_with_untrusted_boundary(SCRNA_SYSTEM),
                tier=TaskTier.HEAVY, max_tokens=1500,
            )
            cell_types = self._parse_annotation_json(raw)
        except Exception as e:
            log.warning(f"LLM annotation failed: {e}")

        if not cell_types:
            # Fall back to raw CellTypist labels if we have them; otherwise keep
            # clusters explicitly unresolved. Runtime marker panels would encode
            # species/tissue assumptions and violate ADR-011.
            if celltypist_result.get("status") == "success":
                cell_types = {
                    cl: {
                        "cell_type":               info["label"],
                        "celltypist_label":        info["label"],
                        "agrees_with_celltypist":  True,
                        "confidence":              "medium",
                        "rationale": (
                            f"LLM unavailable; using CellTypist label directly "
                            f"({info['frequency']*100:.0f}% of cluster)."
                        ),
                        "key_markers": markers_for_prompt.get(cl, [])[:5],
                    }
                    for cl, info in celltypist_result.get("per_cluster", {}).items()
                }
            else:
                cell_types = self._marker_based_annotation(markers_for_prompt)

        # Confidence summary for the bus message.
        agree_count = sum(
            1 for v in cell_types.values()
            if isinstance(v, dict) and v.get("agrees_with_celltypist") is True
        )
        labels_preview = [
            v.get("cell_type", "?") if isinstance(v, dict) else str(v)
            for v in list(cell_types.values())[:5]
        ]

        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Annotated {len(markers_for_prompt)} clusters "
                f"(CellTypist: {celltypist_result.get('model_used', 'N/A')}, "
                f"{agree_count}/{len(cell_types)} agree with LLM): "
                f"{labels_preview}"
             ),
             "cell_types":   cell_types,
             "celltypist":   {
                 "ran":         celltypist_result.get("status") == "success",
                 "model_used":  celltypist_result.get("model_used"),
                 "tissue_hint": tissue_hint,
                 "per_cluster": celltypist_result.get("per_cluster", {}),
             },
             "markers_used": {k: v[:5] for k, v in markers_for_prompt.items()}},
            Confidence.HIGH if (
                celltypist_result.get("status") == "success" and
                agree_count == len(cell_types)
            ) else Confidence.MEDIUM,
        )
        celltypist_success = celltypist_result.get("status") == "success"
        annotated_h5ad = celltypist_result.get("output_path")
        label_col = celltypist_result.get("label_col")
        trusted_groupby = label_col if celltypist_success else None
        annotation_source = "celltypist" if celltypist_success else "llm_marker_only"
        if (not celltypist_success and cell_types
                and all(isinstance(v, dict)
                        and v.get("annotation_source") == "unresolved_marker_fallback"
                        for v in cell_types.values())):
            annotation_source = "unresolved_marker_fallback"
        if not annotated_h5ad and cell_types and celltypist_success:
            applied = self.env.run_in_stack(
                stack="rna",
                script_path="aria/scripts/rna_apply_cluster_labels.py",
                params={
                    "data_path": clustered_h5ad,
                    "labels": cell_types,
                    "cluster_col": "leiden",
                    "label_col": "cell_type_marker",
                    "output_dir": str(self._workspace(experiment_id, "annotation")),
                },
            )
            if applied.get("status") == "success":
                annotated_h5ad = applied.get("output_path")
                label_col = applied.get("label_col", "cell_type_marker")
                trusted_groupby = label_col

        return {
            "cell_types":     cell_types,
            "markers_used":   markers_for_prompt,
            "n_clusters":     len(markers_for_prompt),
            "celltypist":     celltypist_result,
            "tissue_hint":    tissue_hint,
            "label_col":      label_col,
            "annotation_for_report": cell_types,
            "trusted_groupby_for_inference": trusted_groupby,
            "annotated_h5ad": annotated_h5ad,
            "annotation_source": annotation_source,
        }

    @staticmethod
    def _marker_based_annotation(markers_by_cluster: dict) -> dict:
        out = {}
        for cluster, markers in markers_by_cluster.items():
            out[str(cluster)] = {
                "cell_type": f"Unresolved cluster {cluster}",
                "confidence": "low",
                "key_markers": list(markers or [])[:5],
                "rationale": (
                    "LLM/CellTypist unavailable; ARIA does not infer cell "
                    "identity from hardcoded marker panels."
                ),
                "annotation_source": "unresolved_marker_fallback",
            }
        return out

    @staticmethod
    def _parse_annotation_json(raw: str) -> dict:
        """
        Robust JSON extraction: handles ```json fences, leading prose, and
        bare JSON. Returns {} if nothing parseable is found.
        """
        if not raw:
            return {}
        stripped = raw.strip()
        # Markdown fence
        fence = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", stripped, re.DOTALL)
        candidate = fence.group(1) if fence else None
        if candidate is None:
            # Greedy first balanced object
            m = re.search(r"\{.*\}", stripped, re.DOTALL)
            candidate = m.group(0) if m else None
        if candidate is None:
            return {}
        try:
            data = json.loads(candidate)
            return data if isinstance(data, dict) else {}
        except (json.JSONDecodeError, TypeError):
            return {}

    # ── Differential expression per cluster ──────────────────────────────

