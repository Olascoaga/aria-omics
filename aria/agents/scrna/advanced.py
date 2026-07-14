"""scRNAAgent trajectory / cell communication mixin (A7 extraction from scrna_agent.py; bodies verbatim)."""
from __future__ import annotations

from aria.agents.scrna._base import *  # noqa: F401,F403


class AdvancedMixin:
    def _needs_trajectory(self, intent: dict) -> bool:
        keywords = ["differentiat", "develop", "pseudotime", "trajectory",
                    "progenitor", "stem", "lineage", "time course",
                    "progression", "maturation", "hematopoiesis"]
        text = (intent.get("summary", "") + " " +
                " ".join(intent.get("biological_entities", []))).lower()
        return any(kw in text for kw in keywords)

    def _run_trajectory(self, experiment_id: str,
                         clustered_h5ad: str,
                         annotation: dict, intent: dict) -> dict:
        # Prefer CellTypist labels for group naming if available — the
        # pseudotime-by-group output is far more interpretable with real
        # cell type names than with leiden numbers.
        cell_type_col = (
            "cell_type_celltypist"
            if annotation.get("celltypist", {}).get("status") == "success"
            else annotation.get("label_col") or "leiden"
        )
        self._log_decision(
            experiment_id,
            checkpoint="scRNA",
            question="Trajectory grouping",
            decision=f"PAGA/DPT grouped by {cell_type_col}",
            rationale=(
                "Trajectory analysis uses the same trusted cell grouping used "
                "for annotation and reporting."
            ),
            made_by="scrna_agent",
        )
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_trajectory.py",
            params={
                "data_path":      clustered_h5ad,
                "root_cell_type": intent.get("root_cell_type"),
                "cell_type_col":  cell_type_col,
                "output_dir":     str(self._workspace(experiment_id, "trajectory")),
            },
        )

        if result.get("status") == "success":
            paga       = result.get("paga", {})
            top_conn   = paga.get("top_connections", {})
            pseudotime = result.get("pseudotime", {})
            velocity   = result.get("velocity", {})

            self.publish_finding(
                experiment_id,
                {"summary": f"Trajectory: PAGA {paga.get('n_connections', 0)} transitions, "
                            f"DPT computed={pseudotime.get('computed', False)}, "
                            f"RNA velocity computed={velocity.get('computed', False)}.",
                 "paga_top":   top_conn,
                 "pseudotime": pseudotime},
                Confidence.MEDIUM,
            )
        else:
            log.warning(f"Trajectory failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")
        return result

    # ── Cell-cell communication ───────────────────────────────────────────

    def _needs_cell_communication(self, intent: dict) -> bool:
        keywords = ["signal", "interact", "ligand", "receptor", "crosstalk",
                    "communication", "niche", "paracrine", "secreted",
                    "co-culture", "coculture", "microenvironment"]
        text = (intent.get("summary", "") + " " +
                " ".join(intent.get("biological_entities", []))).lower()
        return any(kw in text for kw in keywords)

    def _run_cell_communication(self, experiment_id: str,
                                 clustered_h5ad: str,
                                 exp_ctx: dict,
                                 annotation: dict | None = None) -> dict:
        # Use CellTypist labels as cell type groups when available; falls
        # back to leiden otherwise. rna_cellcomm.py already auto-falls to
        # leiden if the requested column is missing.
        cell_type_col = (
            "cell_type_celltypist"
            if (annotation or {}).get("celltypist", {}).get("status") == "success"
            else (annotation or {}).get("label_col") or "leiden"
        )
        self._log_decision(
            experiment_id,
            checkpoint="scRNA",
            question="Cell-cell communication grouping",
            decision=f"LIANA grouped by {cell_type_col}; n_perms=1000",
            rationale=(
                "Ligand-receptor analysis requires annotated sender and "
                "receiver groups; ARIA reused the trusted grouping column."
            ),
            made_by="scrna_agent",
        )
        result = self.env.run_in_stack(
            stack="rna",
            script_path="aria/scripts/rna_cellcomm.py",
            params={
                "data_path":     clustered_h5ad,
                # P0-1: canonical IPC key is `groupby` (the contract requires it);
                # `cell_type_col` remains an accepted alias for legacy callers.
                "groupby":       cell_type_col,
                "organism":      exp_ctx.get("organism", "Homo sapiens"),
                "n_perms":       1000,
                "output_dir":    str(self._workspace(experiment_id, "cellcomm")),
            },
        )

        if result.get("status") == "success":
            n_ia      = result.get("n_interactions", 0)
            n_types   = result.get("n_cell_types", 0)
            method    = result.get("method", "?")
            top_pairs = result.get("top_pairs", [])
            self.publish_finding(
                experiment_id,
                {"summary": f"Cell-cell communication ({method}): "
                            f"{n_ia} interactions across {n_types} clusters. "
                            f"Top pairs: {', '.join(top_pairs[:3])}."},
                Confidence.MEDIUM,
            )
        elif result.get("status") == "skipped":
            log.info(f"Cell-comm skipped: {result.get('reason', '?')}")
        else:
            log.warning(f"Cell-comm failed: "
                        f"{result.get('error_type', '?')} — "
                        f"{result.get('details', '')[:200]}")
        return result

