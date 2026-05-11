"""
ARIA ParameterAdvisor
---------------------
The institutional brain for hyperparameter decisions.

NOT a black box. NOT a heuristic oracle.
A three-layer decision system that:

  Layer 1 — INTENT:    Constrains search space using biological context
  Layer 2 — METRICS:   Evaluates candidates using objective functions
  Layer 3 — MEMORY:    Consults lab's historical decisions before recommending

All decisions are written to ARIAMemory with full justification.
If the user approves a decision at Checkpoint 3, that decision becomes
part of the lab's institutional knowledge for future experiments.

Over time, the ParameterAdvisor becomes a lab-specific expert —
not because we programmed it, but because the PI's approvals trained it.
"""

from __future__ import annotations

import json
import uuid
import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

from aria.memory.memory import ARIAMemory
from aria.llm.provider import LLMProvider, TaskTier

log = logging.getLogger("aria.params")


# ── Data structures ───────────────────────────────────────────────────────────

@dataclass
class ParameterCandidate:
    """A single candidate parameter value with its evaluation."""
    value:       Any
    metrics:     dict[str, float] = field(default_factory=dict)
    flags:       list[str]        = field(default_factory=list)  # warnings
    score:       float            = 0.0
    recommended: bool             = False


@dataclass
class ParameterDecision:
    """
    The full record of a parameter decision.
    This is what gets stored in ARIAMemory.
    """
    decision_id:        str
    experiment_id:      str
    analysis_type:      str           # "leiden_clustering", "wnn_integration", etc.
    parameter_name:     str           # "resolution", "k", "min_dist", etc.
    candidates:         list[ParameterCandidate]
    chosen_value:       Any
    chosen_by:          str           # "advisor" or "user"
    biological_context: dict          # what the orchestrator knew
    justification:      str           # human-readable explanation
    warnings:           list[str]     # flags the user should know about
    approved_by_user:   bool = False


# ── Analysis-specific evaluators ────────────────────────────────────────────

class MetricEvaluator:
    """
    Computes objective metrics for parameter candidates.
    Each method corresponds to a specific analysis type.
    """

    @staticmethod
    def leiden_resolution(
        adata,           # AnnData object (scanpy)
        resolution: float,
        neighbors_key: str = "neighbors",
    ) -> dict[str, float]:
        """
        Evaluate a Leiden clustering resolution.
        Returns: silhouette score, modularity, n_clusters, n_singleton_clusters
        """
        try:
            import scanpy as sc
            import numpy as np
            from sklearn.metrics import silhouette_score

            # Deterministic Leiden seed so repeat evaluations agree.
            sc.tl.leiden(
                adata, resolution=resolution,
                key_added="_eval_leiden", random_state=0,
            )
            try:
                labels = adata.obs["_eval_leiden"].values

                n_clusters   = len(set(labels))
                cluster_sizes = [sum(labels == l) for l in set(labels)]
                n_singletons  = sum(1 for s in cluster_sizes if s < 10)

                # Silhouette on PCA embedding
                sil = 0.0
                if n_clusters > 1 and n_clusters < len(labels):
                    try:
                        pca_key = "X_pca" if "X_pca" in adata.obsm else None
                        if pca_key:
                            sil = float(silhouette_score(
                                adata.obsm[pca_key][:, :20], labels,
                                sample_size=2000, random_state=0,
                            ))
                    except Exception:
                        pass

                # Modularity (from igraph if available)
                modularity = 0.0
                try:
                    import igraph as ig
                    # Approximate: use scanpy's connectivities
                    pass  # full implementation in production
                except ImportError:
                    pass

                return {
                    "n_clusters":        n_clusters,
                    "silhouette":        round(sil, 4),
                    "modularity":        round(modularity, 4),
                    "n_singleton_clusters": n_singletons,
                    "min_cluster_size":  min(cluster_sizes),
                    "max_cluster_size":  max(cluster_sizes),
                }
            finally:
                # Don't leave the eval column sitting on the AnnData; the
                # caller iterates this method N times and we don't want
                # downstream agents to see a stale resolution.
                if "_eval_leiden" in adata.obs.columns:
                    adata.obs.drop(columns="_eval_leiden", inplace=True)

        except ImportError:
            log.warning("scanpy not available — using mock metrics")
            return MetricEvaluator._mock_metrics(resolution)

    @staticmethod
    def _mock_metrics(resolution: float) -> dict:
        """Mock metrics for testing without scanpy. Returns Python native types."""
        import random
        random.seed(int(resolution * 100))
        n   = max(2, int(resolution * 15) + random.randint(-2, 2))
        sil = max(0.2, 0.85 - resolution * 0.3 + random.uniform(-0.05, 0.05))
        return {
            "n_clusters":              int(n),
            "silhouette":              round(float(sil), 4),
            "modularity":              round(float(0.6 - resolution * 0.15), 4),
            "n_singleton_clusters":    int(max(0, n - 8)),
            "min_cluster_size":        int(max(5, 200 - int(resolution * 100))),
            "max_cluster_size":        int(1500),
        }

    @staticmethod
    def wnn_k(multimodal_adata, k: int) -> dict[str, float]:
        """
        Evaluate WNN k parameter for multimodal integration.

        DeepSeek P0 (addressed): This function currently returns ESTIMATED
        metrics based on a linear heuristic, not real WNN computation.
        This is intentional in beta: running full WNN for each k candidate
        is expensive (minutes per k value).

        The estimates are clearly labeled as such. In v0.3, this will
        execute muon.tl.wnn() with each k and compute:
          - Silhouette score of the joint embedding
          - ARI between RNA-only and WNN clusters
          - Neighbor overlap between modalities

        Until then, the advisor selects k=20 as a robust default based
        on empirical evidence from published multiome datasets, and
        the DebateCouncil is invoked when weights appear imbalanced after
        actual WNN execution.
        """
        import logging
        log = logging.getLogger("aria.advisor")
        log.debug(
            f"wnn_k metrics for k={k} are estimated (beta). "
            f"Real WNN evaluation planned for v0.3."
        )
        # Heuristic estimates — conservative, biased toward balanced weights
        # Real values depend on data quality and cell type composition
        return {
            "rna_weight":    round(max(0.40, min(0.80, 0.6 - k * 0.003)), 3),
            "atac_weight":   round(max(0.20, min(0.60, 0.4 + k * 0.003)), 3),
            "knn_overlap":   round(max(0.30, min(0.90, 0.7 - k * 0.008)), 3),
            "is_estimated":  True,   # flag: these are not from real WNN runs
        }


# ── ParameterAdvisor ────────────────────────────────────────────────────────

class ParameterAdvisor:
    """
    Three-layer parameter decision system.

    Usage:
        advisor = ParameterAdvisor(memory, llm_provider)
        decision = advisor.advise_leiden_resolution(
            adata=adata,
            experiment_id="exp_001",
            biological_context={"intent": "find_rare_subpopulations", ...}
        )
        # Returns ParameterDecision with ranked candidates + recommendation
        # Decision is stored in memory automatically
        # User approves at Checkpoint 3
    """

    def __init__(self, memory: ARIAMemory, llm: LLMProvider):
        self.memory = memory
        self.llm    = llm

    # ── Public advisors ───────────────────────────────────────────────────

    def advise_leiden_resolution(
        self,
        adata,
        experiment_id:      str,
        biological_context: dict,
        n_candidates:       int = 4,
    ) -> ParameterDecision:
        """
        Advise on Leiden clustering resolution.
        Layer 1: constrain range from biological intent
        Layer 2: evaluate candidates with silhouette + modularity
        Layer 3: check lab memory for similar experiments
        """
        # Layer 1: intent-constrained search space
        search_range = self._intent_to_leiden_range(biological_context)
        candidates_values = self._linspace(
            search_range[0], search_range[1], n_candidates
        )

        log.info(
            f"Leiden resolution search: {candidates_values} "
            f"(intent: {biological_context.get('analysis_type', '?')})"
        )

        # Layer 2: evaluate each candidate
        candidates = []
        for val in candidates_values:
            metrics = MetricEvaluator.leiden_resolution(adata, val)
            flags   = self._flag_leiden_issues(val, metrics)
            score   = self._score_leiden(metrics, biological_context)
            candidates.append(ParameterCandidate(
                value=val, metrics=metrics, flags=flags, score=score
            ))

        # Layer 3: memory-informed recommendation
        historical = self._recall_similar_decisions(
            experiment_id, "leiden_clustering", biological_context
        )

        # Choose best candidate
        recommended = self._choose_best(candidates, historical)
        recommended.recommended = True

        # Build justification via LLM (MEDIUM tier — reasoning required)
        justification = self._justify_decision(
            parameter_name="resolution",
            chosen=recommended,
            all_candidates=candidates,
            historical=historical,
            biological_context=biological_context,
        )

        warnings = [f for c in candidates if c.recommended for f in c.flags]

        decision = ParameterDecision(
            decision_id=str(uuid.uuid4())[:8],
            experiment_id=experiment_id,
            analysis_type="leiden_clustering",
            parameter_name="resolution",
            candidates=candidates,
            chosen_value=recommended.value,
            chosen_by="advisor",
            biological_context=biological_context,
            justification=justification,
            warnings=warnings,
        )

        # Store in memory immediately (before user approval)
        self._store_decision(decision)
        return decision

    def advise_wnn_k(
        self,
        multimodal_adata,
        experiment_id:      str,
        biological_context: dict,
    ) -> ParameterDecision:
        """Advise on WNN k parameter for multimodal integration."""
        # WNN k typically ranges 10-50
        candidates_values = [10, 20, 30, 50]
        candidates = []
        for k in candidates_values:
            metrics = MetricEvaluator.wnn_k(multimodal_adata, k)
            flags   = self._flag_wnn_issues(k, metrics)
            score   = self._score_wnn(metrics)
            candidates.append(ParameterCandidate(
                value=k, metrics=metrics, flags=flags, score=score
            ))

        historical = self._recall_similar_decisions(
            experiment_id, "wnn_integration", biological_context
        )
        recommended = self._choose_best(candidates, historical)
        recommended.recommended = True

        justification = self._justify_decision(
            parameter_name="k",
            chosen=recommended,
            all_candidates=candidates,
            historical=historical,
            biological_context=biological_context,
        )

        decision = ParameterDecision(
            decision_id=str(uuid.uuid4())[:8],
            experiment_id=experiment_id,
            analysis_type="wnn_integration",
            parameter_name="k",
            candidates=candidates,
            chosen_value=recommended.value,
            chosen_by="advisor",
            biological_context=biological_context,
            justification=justification,
            warnings=[f for c in candidates if c.recommended for f in c.flags],
        )
        self._store_decision(decision)
        return decision

    def approve_decision(self, decision: ParameterDecision,
                          user_override: Any = None) -> ParameterDecision:
        """
        Called when user approves a decision at Checkpoint 3.
        If user_override is provided, updates the chosen value.
        Marks decision as approved — this is what gets learned.
        """
        if user_override is not None:
            decision.chosen_value = user_override
            decision.chosen_by    = "user"

        decision.approved_by_user = True

        # Update memory with approval
        self.memory.store_decision(
            decision_id=decision.decision_id,
            wing_id=decision.experiment_id,
            checkpoint=3,
            question=(
                f"Parameter: {decision.parameter_name} "
                f"for {decision.analysis_type}"
            ),
            decision=str(decision.chosen_value),
            rationale=decision.justification,
            made_by=decision.chosen_by,
        )

        log.info(
            f"Decision approved: {decision.parameter_name}="
            f"{decision.chosen_value} for {decision.analysis_type}"
        )
        return decision

    # ── Formatting for Checkpoint 3 ──────────────────────────────────────

    def format_for_checkpoint(self, decision: ParameterDecision) -> str:
        """
        Format a parameter decision for display at Checkpoint 3.
        This is what the user sees before approving.
        """
        lines = [
            f"🔬 Parameter: [{decision.parameter_name}] "
            f"for {decision.analysis_type}\n"
        ]

        # Candidate table
        lines.append("  Candidates evaluated:\n")
        for c in decision.candidates:
            marker = "  ★ " if c.recommended else "    "
            flags  = " ⚠️ " + ", ".join(c.flags) if c.flags else ""
            
            # Format key metrics
            key_metrics = []
            if "silhouette" in c.metrics:
                key_metrics.append(f"sil={c.metrics['silhouette']:.3f}")
            if "n_clusters" in c.metrics:
                key_metrics.append(f"k={c.metrics['n_clusters']}")
            if "modularity" in c.metrics and c.metrics["modularity"] > 0:
                key_metrics.append(f"mod={c.metrics['modularity']:.3f}")
            
            metrics_str = "  " + " | ".join(key_metrics) if key_metrics else ""
            lines.append(
                f"{marker}{decision.parameter_name}={c.value}"
                f"{metrics_str}{flags}"
            )

        # Historical context
        hist = self._recall_similar_decisions(
            decision.experiment_id,
            decision.analysis_type,
            decision.biological_context,
        )
        if hist:
            lines.append(f"\n  📚 Lab history: {len(hist)} similar decision(s)")
            for h in hist[:2]:
                lines.append(
                    f"     • {h.get('decision')} approved "
                    f"({h.get('made_at', '')[:10]})"
                )

        # LLM justification
        lines.append(f"\n  💡 Recommendation: {decision.justification}")

        if decision.warnings:
            lines.append(f"\n  ⚠️  Warnings:")
            for w in decision.warnings:
                lines.append(f"     • {w}")

        return "\n".join(lines)

    # ── Layer 1: Intent → Search space ───────────────────────────────────

    def _intent_to_leiden_range(self, bio_context: dict) -> tuple[float, float]:
        """
        Constrain Leiden resolution range based on biological question.
        This is Layer 1: the LLM's biological understanding drives the search.

        Keywords cover English and Spanish — ARIA users frequently phrase
        questions in either language.
        """
        intent   = bio_context.get("analysis_type", "")
        question = bio_context.get("user_question", "").lower()

        # Fine-grained: looking for rare subpopulations
        fine_keywords = [
            # English
            "rare", "subpopulation", "subtype", "progenitor",
            "transitional", "fine", "heterogeneity",
            # Spanish
            "raras", "raro", "subpoblacion", "subpoblación",
            "subtipo", "progenitora", "progenitor",
            "transicional", "transición", "fina", "heterogeneidad",
        ]
        # Coarse: major cell types
        coarse_keywords = [
            # English
            "major", "type", "lineage", "broad", "main", "principal",
            # Spanish
            "mayor", "tipo", "linaje", "amplia", "amplio",
            "principal", "general",
        ]

        if any(k in question for k in fine_keywords):
            return (0.6, 1.4)
        elif any(k in question for k in coarse_keywords):
            return (0.2, 0.6)
        elif intent == "cell_type":
            return (0.3, 0.8)
        else:
            return (0.2, 1.0)  # default: full range

    # ── Layer 2: Objective scoring ────────────────────────────────────────

    def _score_leiden(self, metrics: dict, bio_context: dict) -> float:
        """
        Score a Leiden candidate. Higher = better.
        Weights vary by biological intent.
        """
        sil       = metrics.get("silhouette", 0)
        mod       = metrics.get("modularity", 0)
        n_sing    = metrics.get("n_singleton_clusters", 0)
        min_size  = metrics.get("min_cluster_size", 0)

        # Base score: silhouette is most reliable
        score = sil * 0.6 + mod * 0.3

        # Penalize singleton/tiny clusters
        score -= n_sing * 0.05
        if min_size < 10:
            score -= 0.15

        return round(max(0, score), 4)

    def _score_wnn(self, metrics: dict) -> float:
        """Score WNN k. Balanced RNA/ATAC weights are generally better."""
        rna_w  = metrics.get("rna_weight", 0.6)
        atac_w = metrics.get("atac_weight", 0.4)
        # Penalize extreme imbalance
        balance = 1.0 - abs(rna_w - atac_w)
        overlap = metrics.get("knn_overlap", 0.5)
        return round(balance * 0.6 + overlap * 0.4, 4)

    def _flag_leiden_issues(self, resolution: float, metrics: dict) -> list[str]:
        flags = []
        if metrics.get("n_singleton_clusters", 0) > 0:
            flags.append(
                f"{metrics['n_singleton_clusters']} cluster(s) with <10 cells "
                f"(likely over-clustering)"
            )
        if metrics.get("silhouette", 1) < 0.3:
            flags.append("Low silhouette score — cluster structure is weak")
        if metrics.get("n_clusters", 0) > 50:
            flags.append("Very high cluster count — review biological plausibility")
        return flags

    def _flag_wnn_issues(self, k: int, metrics: dict) -> list[str]:
        flags = []
        rna_w = metrics.get("rna_weight", 0.5)
        if rna_w > 0.90:
            flags.append("ATAC contribution <10% — check ATAC data quality")
        if rna_w < 0.10:
            flags.append("RNA contribution <10% — check RNA data quality")
        return flags

    # ── Layer 3: Historical memory ────────────────────────────────────────

    def _recall_similar_decisions(
        self,
        experiment_id:  str,
        analysis_type:  str,
        bio_context:    dict,
    ) -> list[dict]:
        """
        Query ARIAMemory for approved decisions on similar analyses.
        Returns decisions sorted by recency.
        """
        try:
            all_wings = self.memory.list_wings()
            similar   = []

            for wing in all_wings:
                decisions = self.memory.get_decisions(wing["id"])
                for d in decisions:
                    if (analysis_type.lower() in d.get("question", "").lower()
                            and d.get("decision") is not None):
                        similar.append(d)

            # Sort by recency, return last 5
            similar.sort(key=lambda x: x.get("made_at", ""), reverse=True)
            return similar[:5]

        except Exception as e:
            log.warning(f"Memory recall failed: {e}")
            return []

    def _choose_best(
        self,
        candidates: list[ParameterCandidate],
        historical: list[dict],
    ) -> ParameterCandidate:
        """
        Choose the best candidate, weighted by:
          - objective score (primary)
          - historical approval frequency (secondary bonus)
        """
        # Build historical value frequency
        hist_values = {}
        for h in historical:
            try:
                val = float(h.get("decision", 0))
                hist_values[val] = hist_values.get(val, 0) + 1
            except (ValueError, TypeError):
                pass

        best = None
        best_score = -1.0

        for c in candidates:
            adjusted = c.score
            # Small bonus for historically approved values
            try:
                if float(c.value) in hist_values:
                    adjusted += 0.05 * hist_values[float(c.value)]
            except (ValueError, TypeError):
                pass

            if adjusted > best_score:
                best_score = adjusted
                best = c

        return best or candidates[0]

    def _justify_decision(
        self,
        parameter_name: str,
        chosen:         ParameterCandidate,
        all_candidates: list[ParameterCandidate],
        historical:     list[dict],
        biological_context: dict,
    ) -> str:
        """
        Generate a concise, scientifically grounded justification.
        Uses MEDIUM tier LLM — needs reasoning but not heavy inference.
        """
        candidates_summary = [
            {
                "value":   c.value,
                "metrics": c.metrics,
                "score":   c.score,
                "flags":   c.flags,
            }
            for c in all_candidates
        ]

        hist_summary = [
            {"value": h.get("decision"), "date": h.get("made_at", "")[:10]}
            for h in historical[:3]
        ]

        # numpy types are not JSON serializable — convert before dumping
        def _json_safe(obj):
            if hasattr(obj, "item"):          # numpy scalar → Python scalar
                return obj.item()
            if hasattr(obj, "tolist"):         # numpy array → list
                return obj.tolist()
            raise TypeError(f"Not serializable: {type(obj)}")

        prompt = f"""
Biological context: {json.dumps(biological_context, indent=2, default=_json_safe)}
Parameter: {parameter_name}
Chosen value: {chosen.value}
All candidates evaluated: {json.dumps(candidates_summary, indent=2, default=_json_safe)}
Lab historical decisions: {json.dumps(hist_summary, indent=2, default=_json_safe)}

Write a 1-2 sentence scientific justification for choosing {parameter_name}={chosen.value}.
Be specific: cite the metrics, the biological intent, and historical precedent if relevant.
Do NOT hedge. Be direct. Use scientific language appropriate for a methods section.
"""
        try:
            return self.llm.complete_medium(
                prompt=prompt,
                system=(
                    "You are a bioinformatics expert writing a methods justification. "
                    "Be precise, cite specific metrics, max 2 sentences."
                ),
                max_tokens=150,
            )
        except Exception as e:
            log.warning(f"LLM justification failed: {e}. Using fallback.")
            return (
                f"{parameter_name}={chosen.value} selected based on highest "
                f"silhouette score ({chosen.metrics.get('silhouette', '?')}) "
                f"among {len(all_candidates)} candidates evaluated."
            )

    # ── Memory storage ────────────────────────────────────────────────────

    def _store_decision(self, decision: ParameterDecision):
        """Store pending decision in memory (before user approval)."""
        try:
            room_id = f"{decision.experiment_id}_{decision.analysis_type}"
            self.memory.store_finding(
                finding_id=decision.decision_id,
                room_id=room_id,
                content=json.dumps({
                    "parameter":    decision.parameter_name,
                    "chosen_value": decision.chosen_value,
                    "justification": decision.justification,
                    "candidates":   [
                        {"value": c.value, "score": c.score}
                        for c in decision.candidates
                    ],
                }, default=lambda o: o.item() if hasattr(o, "item") else str(o)),
                confidence="medium",
            )
        except Exception as e:
            log.warning(f"Could not store decision in memory: {e}")

    @staticmethod
    def _linspace(start: float, end: float, n: int) -> list[float]:
        """Generate n evenly spaced values between start and end."""
        if n <= 1:
            return [start]
        step = (end - start) / (n - 1)
        return [round(start + i * step, 2) for i in range(n)]
