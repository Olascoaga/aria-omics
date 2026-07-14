"""scRNAAgent differential expression / pathway mixin (A7 extraction from scrna_agent.py; bodies verbatim)."""
from __future__ import annotations

from aria.agents.scrna._base import *  # noqa: F401,F403


class DEMixin:
    def _differential_expression(self, experiment_id: str,
                                   clustered_h5ad: str,
                                   intent: dict, exp_ctx: dict,
                                   groupby: str = "leiden") -> dict:
        # Respect user-confirmed thresholds from CP3 (Strict/Standard/etc.).
        padj_max = float(exp_ctx.get("global_padj", 0.05))
        lfc_min  = float(exp_ctx.get("global_lfc", 0.5))

        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_de_per_cluster.py",
            params={
                "data_path": clustered_h5ad,
                "groupby":   groupby,
                "padj_max":  padj_max,
                "lfc_min":   lfc_min,
                "top_n":     20,
                "output_dir": str(self._workspace(experiment_id, "de_per_cluster")),
            },
        )

        if result.get("status") != "success":
            log.warning(f"DE per cluster failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")
            return result

        n_sig = result.get("n_significant_total", 0)
        conf  = (Confidence.HIGH   if n_sig > 50 else
                 Confidence.MEDIUM if n_sig > 10 else
                 Confidence.LOW    if n_sig > 0  else
                 Confidence.INSUFFICIENT)

        # LLM interpretation grounded in actual numbers.
        de_for_llm = {
            k: v[:5] for k, v in result.get("de_genes_by_cluster", {}).items()
        }
        try:
            interpretation = self.llm.complete(
                prompt=build_untrusted_prompt(
                    task=(
                        "Interpret the supplied per-cluster differential-expression "
                        "results in 3-4 sentences. Focus on biology, not statistics, "
                        "and do not introduce entities absent from the results."
                    ),
                    fields=[
                        PromptDataField("biological_question", intent.get(
                            "summary", ""
                        ), "user_text", "user"),
                        PromptDataField("organism", exp_ctx.get("organism", ""),
                                        "metadata", "data_audit"),
                        PromptDataField("top_de_by_cluster", de_for_llm,
                                        "structured_result", "rna_de_per_cluster"),
                        PromptDataField("padj_max", padj_max, "metadata",
                                        "analysis_thresholds"),
                        PromptDataField("lfc_min", lfc_min, "metadata",
                                        "analysis_thresholds"),
                        PromptDataField("n_significant_genes", n_sig, "metadata",
                                        "rna_de_per_cluster"),
                    ],
                    response_contract="Return plain scientific prose only.",
                ),
                system=system_with_untrusted_boundary(SCRNA_SYSTEM),
                tier=TaskTier.HEAVY, max_tokens=400,
            )
        except Exception:
            interpretation = f"{n_sig} significant DE genes across clusters."

        self.publish_finding(
            experiment_id,
            {"summary":             interpretation,
             "n_significant_genes": n_sig,
             "thresholds":          {"padj_max": padj_max, "lfc_min": lfc_min}},
            conf,
        )

        return {
            "status":               result["status"],
            "padj_max":             padj_max,
            "lfc_min":              lfc_min,
            "n_significant_genes":  n_sig,
            "de_genes_by_cluster":  result.get("de_genes_by_cluster", {}),
            "n_sig_by_cluster":     result.get("n_sig_by_cluster", {}),
            "interpretation":       interpretation,
            "output_csv":           result.get("output_csv"),
        }

    # ── Pathway enrichment per cluster (ORA) ─────────────────────────────

    def _run_pathway_per_cluster(self, experiment_id: str,
                                   de_result: dict, exp_ctx: dict,
                                   data_path: str | None = None) -> dict:
        """
        Run ORA per Leiden cluster against GO_BP / KEGG / Reactome via the
        rna_pathway_per_cluster.py subprocess. Anchors cell-type biology in
        actual pathway hits instead of "what does the LLM guess about this
        marker gene".
        """
        de_by_cluster = de_result.get("de_genes_by_cluster", {})
        if not de_by_cluster:
            return {"status": "skipped", "reason": "no DE genes to enrich"}
        background_genes = self._expressed_gene_background(data_path)

        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_pathway_per_cluster.py",
            params={
                "de_genes_by_cluster":   de_by_cluster,
                "organism":              exp_ctx.get("organism", "Homo sapiens"),
                "top_genes_per_cluster": 200,
                "background_genes":       background_genes,
                "padj_db_max":           0.05,
                "output_dir":            str(self._workspace(experiment_id, "pathways")),
            },
            # P1-7/W-PRIV: ORA runs locally (hypergeometric) by default and is
            # fast. The opt-in Enrichr fallback is rate-limited (~8s/call), so
            # keep a generous 30 min ceiling for dense datasets in that mode.
            timeout=1800,
        )

        if result.get("status") != "success":
            log.warning(
                f"Pathway per cluster failed: "
                f"{result.get('error_type', '?')} — "
                f"{result.get('details', '')[:200]}"
            )
            return result

        per_cluster = result.get("per_cluster", {})
        total_sig   = sum(c.get("n_significant", 0) for c in per_cluster.values())
        top_terms_preview = []
        for cl, info in list(per_cluster.items())[:4]:
            for db, terms in info.get("results", {}).items():
                if terms:
                    top_terms_preview.append(
                        f"cluster {cl}/{db}: {terms[0]['term']}"
                    )
                    break

        self.publish_finding(
            experiment_id,
            {"summary": (
                f"Per-cluster ORA: {total_sig} significant pathway hits "
                f"across {len(per_cluster)} clusters "
                f"({', '.join(result.get('databases', {}).keys())}). "
                f"Top: {' | '.join(top_terms_preview[:3])}."
             ),
             "n_significant": total_sig,
             "databases":     result.get("databases", {}),
             "output_csv":    result.get("output_csv")},
            Confidence.HIGH if total_sig > 20 else
            Confidence.MEDIUM if total_sig > 0 else
            Confidence.INSUFFICIENT,
        )
        return result

    @staticmethod
    def _expressed_gene_background(data_path: str | None) -> list[str]:
        """Return genes detected at least once in the retained h5ad."""
        if not data_path:
            return []
        try:
            import numpy as np
            import anndata as ad
            from scipy import sparse
            adata = ad.read_h5ad(data_path)
            if adata.raw is not None:
                mat = adata.raw.X
                genes = list(adata.raw.var_names)
            else:
                mat = adata.X
                genes = list(adata.var_names)
            if sparse.issparse(mat):
                mask = np.asarray((mat > 0).sum(axis=0)).ravel() > 0
            else:
                mask = (np.asarray(mat) > 0).sum(axis=0) > 0
            return [str(g) for g, keep in zip(genes, mask) if keep]
        except Exception as exc:
            log.debug(f"Could not compute expressed-gene background: {exc}")
            return []

    # ── Pseudobulk DE between conditions ─────────────────────────────────

