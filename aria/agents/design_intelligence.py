"""
ARIA DesignIntelligence
-----------------------
Rules-first feasibility and opportunity assessment across modalities.

This layer does not execute analysis. It inspects the detected data, confirmed
design, and user intent, then records what ARIA should recommend, offer as
optional, or explicitly avoid.
"""

from __future__ import annotations

from pathlib import Path
import re


class DesignIntelligence:
    """Build cross-modality design profiles before computation starts."""

    def evaluate(self, exp_context: dict, intent: dict | None = None) -> dict:
        intent = intent or {}
        modalities = exp_context.get("modalities", {}) or {}
        profiles = []
        for modality in modalities:
            if modality == "scRNA":
                profiles.append(self._scrna_profile(exp_context, intent))
            elif modality in {"bulk_RNA", "bulk_RNA_raw"}:
                profiles.append(self._bulk_rna_profile(modality, exp_context, intent))
            elif modality in {"scATAC", "bulk_ATAC", "ChIP", "CUT_AND_RUN", "CUT_AND_TAG"}:
                profiles.append(self._chromatin_profile(modality, exp_context, intent))
            elif modality == "HiC":
                profiles.append(self._hic_profile(exp_context, intent))
            else:
                profiles.append(self._generic_profile(modality))

        supported = []
        optional = []
        unsupported = []
        warnings = []
        for profile in profiles:
            supported.extend(profile.get("recommended", []))
            optional.extend(profile.get("optional", []))
            unsupported.extend(profile.get("unsupported", []))
            warnings.extend(profile.get("warnings", []))

        integration = self._integration_assessment(modalities, exp_context)
        if integration.get("recommended"):
            optional.append(integration["recommended"])
        unsupported.extend(integration.get("unsupported", []))
        warnings.extend(integration.get("warnings", []))

        return {
            "status": "success",
            "profiles": profiles,
            "recommended": supported,
            "optional": optional,
            "unsupported": unsupported,
            "warnings": warnings,
            "summary": self._summary(supported, optional, unsupported, warnings),
        }

    def _scrna_profile(self, exp_context: dict, intent: dict) -> dict:
        design = exp_context.get("design", {}) or {}
        inferred = exp_context.get("inferred_design", {}) or {}
        pb_cfg = (design.get("pseudobulk") or inferred.get("pseudobulk") or {})
        condition = pb_cfg.get("condition_col") or design.get("main_factor")
        replicate = pb_cfg.get("replicate_col")
        groupby = pb_cfg.get("groupby_col") or inferred.get("groupby_col")
        covariates = pb_cfg.get("covariates") or inferred.get("covariates") or []
        groups = design.get("groups") or inferred.get("groups") or {}
        files = (exp_context.get("modalities", {}) or {}).get("scRNA", [])
        h5ad_meta = self._inspect_h5ad(files[:1])
        focus = self._infer_focus(intent, exp_context, h5ad_meta.get("group_values", {}).get(groupby, []))

        recommended = ["QC using available h5ad obs metrics or count data."]
        optional = []
        unsupported = []
        warnings = []

        if condition and replicate and self._has_replicates(groups, min_reps=2):
            recommended.append(
                f"Donor/sample-level pseudobulk DE: condition={condition}, "
                f"replicate={replicate}, groupby={groupby or 'cell group'}."
            )
            optional.append("Pathway/ORA enrichment on significant pseudobulk DE genes.")
        else:
            unsupported.append(
                "Between-condition scRNA DE is not supported without a condition "
                "column and at least two biological replicates per group."
            )

        if covariates:
            recommended.append(
                "Use candidate covariates if balanced and not collinear: "
                + ", ".join(map(str, covariates))
            )
        else:
            optional.append("No usable covariates detected; run unadjusted model unless user supplies metadata.")

        if focus:
            recommended.append(
                f"Focus computation on {groupby or 'cell group'} in "
                f"{{{', '.join(focus)}}} before downstream analysis."
            )

        n_focus = len(focus)
        if groupby and (not focus or n_focus >= 2):
            optional.append("LIANA non-autocrine ligand-receptor analysis between annotated cell groups.")
        else:
            unsupported.append(
                "LIANA non-autocrine crosstalk is not informative after focusing "
                "to one cell group; report as skipped unless autocrine signaling is requested."
            )

        if h5ad_meta.get("has_velocity_layers"):
            optional.append("RNA velocity with scVelo is possible because spliced/unspliced layers exist.")
        else:
            unsupported.append("RNA velocity is not supported: no spliced/unspliced layers detected.")

        if self._lineage_intent(intent) and (not focus or len(focus) >= 2):
            optional.append("PAGA + DPT trajectory as exploratory manifold ordering, not causal lineage proof.")
        else:
            unsupported.append(
                "PAGA/DPT is not recommended unless the question and selected cell groups define a plausible lineage."
            )

        if groupby:
            optional.append("Cell-type composition can be summarized from donor-level proportions before cell focusing.")
        else:
            unsupported.append("Cell-type composition requires a trusted obs grouping column.")

        if not h5ad_meta.get("raw_counts_likely", True):
            warnings.append("Processed/log-normalized h5ad input may require count recovery for pseudobulk.")

        return {
            "modality": "scRNA",
            "condition": condition,
            "replicate": replicate,
            "groupby": groupby,
            "covariates": covariates,
            "focus": focus,
            "recommended": recommended,
            "optional": optional,
            "unsupported": unsupported,
            "warnings": warnings,
        }

    def _bulk_rna_profile(self, modality: str, exp_context: dict, intent: dict) -> dict:
        design = exp_context.get("design", {}) or {}
        groups = design.get("groups", {}) or {}
        recommended = []
        unsupported = []
        optional = ["PCA/sample-distance QC and pathway enrichment after DE."]
        if self._has_replicates(groups, min_reps=2):
            recommended.append("Bulk RNA DESeq2 differential expression with biological replicates.")
        else:
            unsupported.append("Bulk RNA differential expression needs at least two replicates per condition.")
        if modality == "bulk_RNA_raw":
            recommended.append("FASTQ QC/quantification before DE.")
        return {
            "modality": modality,
            "recommended": recommended,
            "optional": optional,
            "unsupported": unsupported,
            "warnings": [],
        }

    def _chromatin_profile(self, modality: str, exp_context: dict, intent: dict) -> dict:
        recommended = ["Chromatin QC and peak/feature-level summaries when inputs are validated."]
        optional = ["Differential accessibility/occupancy if replicate metadata supports contrasts."]
        unsupported = []
        warnings = []
        if modality == "scATAC":
            optional.extend(["LSI clustering", "motif enrichment", "gene activity summaries"])
            warnings.append("scATAC remains roadmap/beta until v4.4 validation is completed.")
        else:
            optional.append("Peak calling and consensus peak analysis for bulk chromatin assays.")
        return {
            "modality": modality,
            "recommended": recommended,
            "optional": optional,
            "unsupported": unsupported,
            "warnings": warnings,
        }

    def _hic_profile(self, exp_context: dict, intent: dict) -> dict:
        return {
            "modality": "HiC",
            "recommended": ["Hi-C contact-map QC and resolution feasibility assessment."],
            "optional": ["Compartments, TADs, and loops when depth/resolution supports them."],
            "unsupported": ["Causal 3D regulatory claims require orthogonal expression/chromatin evidence."],
            "warnings": [],
        }

    @staticmethod
    def _generic_profile(modality: str) -> dict:
        return {
            "modality": modality,
            "recommended": [f"Inspect {modality} inputs before analysis."],
            "optional": [],
            "unsupported": [],
            "warnings": [],
        }

    @staticmethod
    def _integration_assessment(modalities: dict, exp_context: dict) -> dict:
        assay_modalities = [m for m in modalities if m not in {"unknown"}]
        if len(assay_modalities) < 2:
            return {"unsupported": ["Multimodal integration is not supported with a single detected modality."]}
        return {
            "recommended": "Consider integration only after each modality passes standalone QC.",
            "unsupported": [],
            "warnings": ["Integration requires matched cells, donors, or features; unmatched inputs should not be forced."],
        }

    @staticmethod
    def _has_replicates(groups: dict, min_reps: int = 2) -> bool:
        return bool(groups) and all(len(v or []) >= min_reps for v in groups.values())

    @staticmethod
    def _lineage_intent(intent: dict) -> bool:
        text = " ".join([
            str(intent.get("summary", "")),
            " ".join(intent.get("biological_entities", []) or []),
        ]).lower()
        return any(k in text for k in ("trajectory", "lineage", "differentiat", "pseudotime", "progenitor"))

    def _infer_focus(self, intent: dict, exp_context: dict, available: list[str]) -> list[str]:
        if not available:
            return []
        text = self._cell_focus_text(exp_context, intent)
        if not text:
            return []
        aliases = {
            "microglia": {"Microglia"},
            "microglía": {"Microglia"},
            "opc": {"OPC"},
            "opcs": {"OPC"},
            "oligo": {"Oligo"},
            "oligodendrocyte": {"Oligo"},
            "oligodendrocytes": {"Oligo"},
            "oligodendroglial": {"OPC", "Oligo"},
            "oligodendrocito": {"Oligo"},
            "oligodendrocitos": {"Oligo"},
            "astrocyte": {"Astro"},
            "astrocytes": {"Astro"},
            "astrocito": {"Astro"},
            "astrocitos": {"Astro"},
        }
        focus = set()
        for value in available:
            if re.search(rf"\b{re.escape(value.lower())}\b", text):
                focus.add(value)
        for token, values in aliases.items():
            if re.search(rf"\b{re.escape(token)}\b", text):
                focus.update(v for v in values if v in available)
        return sorted(focus) if 0 < len(focus) < len(available) else []

    @staticmethod
    def _cell_focus_text(exp_context: dict, intent: dict) -> str:
        raw = str((exp_context or {}).get("user_question", "") or "")
        if not raw:
            raw = str((intent or {}).get("user_question", "") or "")
        if not raw:
            raw = str((intent or {}).get("summary", "") or "")
        clauses = [
            c.strip() for c in re.split(r"[\n.;]+", raw)
            if c and c.strip()
        ]
        focus_markers = (
            "focus", "focused", "focusing", "restrict", "restricted",
            "subset", "only", "exclusively", "obs[", "==",
            "solo", "sólo", "unicamente", "únicamente", "enfoc",
            "centr", "limita", "limitar",
        )
        selected = [
            c for c in clauses
            if any(marker in c.lower() for marker in focus_markers)
        ]
        return " ".join(selected).lower()

    @staticmethod
    def _inspect_h5ad(files: list[str]) -> dict:
        meta = {
            "has_velocity_layers": False,
            "raw_counts_likely": True,
            "group_values": {},
        }
        if not files:
            return meta
        try:
            import anndata as ad
        except ImportError:
            return meta
        path = files[0]
        if not str(path).lower().endswith(".h5ad") or not Path(path).exists():
            return meta
        try:
            adata = ad.read_h5ad(path, backed="r")
            layers = set(map(str, adata.layers.keys()))
            meta["has_velocity_layers"] = {"spliced", "unspliced"} <= layers
            for col in ("subclass", "cell_type", "celltype", "cell_type_celltypist", "leiden"):
                if col in adata.obs:
                    vals = [str(v) for v in adata.obs[col].dropna().unique()]
                    meta["group_values"][col] = sorted(v for v in vals if v and v.lower() != "nan")
            if "nCount_RNA" in adata.obs and not layers:
                meta["raw_counts_likely"] = False
            backing_file = getattr(adata, "file", None)
            if backing_file is not None:
                backing_file.close()
        except Exception:
            return meta
        return meta

    @staticmethod
    def _summary(recommended: list[str], optional: list[str],
                 unsupported: list[str], warnings: list[str]) -> str:
        return (
            f"{len(recommended)} recommended, {len(optional)} optional, "
            f"{len(unsupported)} unsupported/not recommended, "
            f"{len(warnings)} warning(s)."
        )


def format_design_intelligence(di: dict, max_items: int = 6) -> str:
    """Human-readable block for checkpoints."""
    if not di:
        return ""
    lines = ["\nDesign Intelligence:"]
    if di.get("recommended"):
        lines.append("  Recommended:")
        for item in di["recommended"][:max_items]:
            lines.append(f"    - {item}")
    if di.get("optional"):
        lines.append("  Optional / supported:")
        for item in di["optional"][:max_items]:
            lines.append(f"    - {item}")
    if di.get("unsupported"):
        lines.append("  Not supported / not recommended:")
        for item in di["unsupported"][:max_items]:
            lines.append(f"    - {item}")
    if di.get("warnings"):
        lines.append("  Warnings:")
        for item in di["warnings"][:max_items]:
            lines.append(f"    - {item}")
    return "\n".join(lines)
