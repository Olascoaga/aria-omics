# Graph Report - /home/medusa/Samael/ARIA  (2026-05-31)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 2345 nodes · 5514 edges · 95 communities (86 shown, 9 thin omitted)
- Extraction: 87% EXTRACTED · 13% INFERRED · 0% AMBIGUOUS · INFERRED: 729 edges (avg confidence: 0.51)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `d9e11ebf`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 10|Community 10]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 13|Community 13]]
- [[_COMMUNITY_Community 14|Community 14]]
- [[_COMMUNITY_Community 15|Community 15]]
- [[_COMMUNITY_Community 16|Community 16]]
- [[_COMMUNITY_Community 17|Community 17]]
- [[_COMMUNITY_Community 18|Community 18]]
- [[_COMMUNITY_Community 19|Community 19]]
- [[_COMMUNITY_Community 20|Community 20]]
- [[_COMMUNITY_Community 21|Community 21]]
- [[_COMMUNITY_Community 22|Community 22]]
- [[_COMMUNITY_Community 23|Community 23]]
- [[_COMMUNITY_Community 24|Community 24]]
- [[_COMMUNITY_Community 25|Community 25]]
- [[_COMMUNITY_Community 26|Community 26]]
- [[_COMMUNITY_Community 27|Community 27]]
- [[_COMMUNITY_Community 28|Community 28]]
- [[_COMMUNITY_Community 29|Community 29]]
- [[_COMMUNITY_Community 30|Community 30]]
- [[_COMMUNITY_Community 31|Community 31]]
- [[_COMMUNITY_Community 32|Community 32]]
- [[_COMMUNITY_Community 33|Community 33]]
- [[_COMMUNITY_Community 34|Community 34]]
- [[_COMMUNITY_Community 35|Community 35]]
- [[_COMMUNITY_Community 36|Community 36]]
- [[_COMMUNITY_Community 37|Community 37]]
- [[_COMMUNITY_Community 38|Community 38]]
- [[_COMMUNITY_Community 39|Community 39]]
- [[_COMMUNITY_Community 40|Community 40]]
- [[_COMMUNITY_Community 41|Community 41]]
- [[_COMMUNITY_Community 42|Community 42]]
- [[_COMMUNITY_Community 43|Community 43]]
- [[_COMMUNITY_Community 44|Community 44]]
- [[_COMMUNITY_Community 45|Community 45]]
- [[_COMMUNITY_Community 46|Community 46]]
- [[_COMMUNITY_Community 47|Community 47]]
- [[_COMMUNITY_Community 48|Community 48]]
- [[_COMMUNITY_Community 49|Community 49]]
- [[_COMMUNITY_Community 50|Community 50]]
- [[_COMMUNITY_Community 51|Community 51]]
- [[_COMMUNITY_Community 52|Community 52]]
- [[_COMMUNITY_Community 53|Community 53]]
- [[_COMMUNITY_Community 54|Community 54]]
- [[_COMMUNITY_Community 55|Community 55]]
- [[_COMMUNITY_Community 56|Community 56]]
- [[_COMMUNITY_Community 57|Community 57]]
- [[_COMMUNITY_Community 58|Community 58]]
- [[_COMMUNITY_Community 59|Community 59]]
- [[_COMMUNITY_Community 60|Community 60]]
- [[_COMMUNITY_Community 61|Community 61]]
- [[_COMMUNITY_Community 62|Community 62]]
- [[_COMMUNITY_Community 63|Community 63]]
- [[_COMMUNITY_Community 64|Community 64]]
- [[_COMMUNITY_Community 65|Community 65]]
- [[_COMMUNITY_Community 66|Community 66]]
- [[_COMMUNITY_Community 67|Community 67]]
- [[_COMMUNITY_Community 68|Community 68]]
- [[_COMMUNITY_Community 69|Community 69]]
- [[_COMMUNITY_Community 70|Community 70]]
- [[_COMMUNITY_Community 71|Community 71]]
- [[_COMMUNITY_Community 72|Community 72]]
- [[_COMMUNITY_Community 73|Community 73]]
- [[_COMMUNITY_Community 74|Community 74]]
- [[_COMMUNITY_Community 75|Community 75]]
- [[_COMMUNITY_Community 76|Community 76]]
- [[_COMMUNITY_Community 77|Community 77]]
- [[_COMMUNITY_Community 78|Community 78]]
- [[_COMMUNITY_Community 79|Community 79]]
- [[_COMMUNITY_Community 80|Community 80]]
- [[_COMMUNITY_Community 81|Community 81]]
- [[_COMMUNITY_Community 82|Community 82]]
- [[_COMMUNITY_Community 83|Community 83]]
- [[_COMMUNITY_Community 84|Community 84]]
- [[_COMMUNITY_Community 85|Community 85]]
- [[_COMMUNITY_Community 86|Community 86]]
- [[_COMMUNITY_Community 87|Community 87]]
- [[_COMMUNITY_Community 88|Community 88]]
- [[_COMMUNITY_Community 89|Community 89]]
- [[_COMMUNITY_Community 90|Community 90]]
- [[_COMMUNITY_Community 91|Community 91]]
- [[_COMMUNITY_Community 92|Community 92]]
- [[_COMMUNITY_Community 94|Community 94]]

## God Nodes (most connected - your core abstractions)
1. `LLMProvider` - 128 edges
2. `Confidence` - 111 edges
3. `BaseAgent` - 107 edges
4. `TaskTier` - 100 edges
5. `NarrativeBlock` - 72 edges
6. `ParameterAdvisor` - 65 edges
7. `NarrativeAgent` - 64 edges
8. `MessageType` - 62 edges
9. `OrchestratorAgent` - 56 edges
10. `scRNAAgent` - 55 edges

## Surprising Connections (you probably didn't know these)
- `NarrativeAgent` --references--> `Claim Compiler`  [INFERRED]
  aria/agents/narrative_agent.py → README.md
- `NarrativeAgent` --references--> `DebateCouncil`  [INFERRED]
  aria/agents/narrative_agent.py → README.md
- `test_shrinkage_clause_in_narrative()` --calls--> `_lfc_shrinkage_clause()`  [EXTRACTED]
  tests/test_lfc_shrinkage.py → aria/agents/_narrative_scrna.py
- `CompletedProcess` --uses--> `DataAuditAgent`  [INFERRED]
  tests/test_pytest_smoke.py → aria/agents/data_audit_agent.py
- `bool` --uses--> `DataAuditAgent`  [INFERRED]
  tests/test_pytest_smoke.py → aria/agents/data_audit_agent.py

## Import Cycles
- None detected.

## Communities (95 total, 9 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.03
Nodes (38): object, _assert_legacy_script_passed(), T1.3 — CP2 recommended-only must not run n=2 pseudobulk DE.     Design Intellige, T1.3 — n=2 pseudobulk remains supported with warnings when CP2     selected reco, T1.1 — scrna_agent._run_pseudobulk must pass     composition_covariate=True to r, T1.4 — scRNA ORA receives the dataset-expressed gene universe, not     only the, Generalization regression — the public Kang et al. PBMC IFN-β     dataset (8 don, When the LLM hallucinates a sentence as a 'factor name', the CP     2.3 menu mus (+30 more)

### Community 1 - "Community 1"
Cohesion: 0.07
Nodes (66): ABC, BaseAgent, ARIA BaseAgent -------------- All ARIA agents inherit from this class. Handles:, Base class for all ARIA agents.      Uses LLMProvider for all LLM calls — comple, _infer_lfc_threshold(), _is_fastq(), ARIA BulkRNAAgent (v3.10 → v4.0) --------------------------------- Bulk RNA-seq, ARIA ChromatinAgent ------------------- Handles all chromatin accessibility and (+58 more)

### Community 2 - "Community 2"
Cohesion: 0.05
Nodes (33): Robust JSON extraction: handles ```json fences, leading prose, and         bare, Run ORA per Leiden cluster against GO_BP / KEGG / Reactome via the         rna_p, Return genes detected at least once in the retained h5ad., True when DesignAgent identified enough biological groups (the         prerequis, Write a copy of `h5ad_in` with `obs[main_factor]` populated from         the sam, Pseudobulk DE between condition groups identified by DesignAgent.          Aggre, If the user asks to focus on specific obs cell types, materialize a         focu, Derive a stable per-sample label from a 10x .h5 / MEX / .h5ad path.         Stri (+25 more)

### Community 3 - "Community 3"
Cohesion: 0.05
Nodes (29): NarrativeAgent, Honest accounting of conflicts and limitations., Construct (and mkdir) ~/.aria/reports/aria_<ts>_<slug>_<suffix>/         with fi, Orchestrate scRNA figure rendering (UMAPs + DE bar + pathway dotplots         +, Render the full HTML report.         Self-contained: CSS embedded, no external d, One-paragraph executive summary for the PI.         Uses deterministic structure, P-LEDGER: render the planned-vs-run manifest. Any analysis the plan         call, Build a short URL-safe slug from biological entities or the question.         Ex (+21 more)

### Community 4 - "Community 4"
Cohesion: 0.06
Nodes (66): Path, str, ask_biological_question(), _discard_queued_stdin_lines(), _drain_checkpoints(), _live_analysis_loop(), main(), print_agent_message() (+58 more)

### Community 5 - "Community 5"
Cohesion: 0.08
Nodes (62): _annotation_state(), build_scrna_html_section(), build_scrna_integrated_interpretation(), build_scrna_methods(), _concise_question(), _describe_abundance_de_relationship(), _describe_cellcomm_context(), _describe_pathway_support() (+54 more)

### Community 6 - "Community 6"
Cohesion: 0.05
Nodes (29): ChromatinAgent, Run chromatin analysis based on available modalities.          context must cont, Single-cell ATAC-seq pipeline.         Key: LSI dimensionality reduction, discar, Bulk ATAC-seq: QC, peak calling, differential accessibility., ChIP-seq: QC, peak calling with input control, IDR for replicates.         Separ, CUT&RUN / CUT&TAG pipeline.         Key differences from ChIP: very low backgrou, Use ParameterAdvisor to decide LSI parameters.          Key decision: how many S, TF motif enrichment analysis.         DebateCouncil is MANDATORY here — Tn5 bias (+21 more)

### Community 7 - "Community 7"
Cohesion: 0.08
Nodes (42): DataAuditAgent, ARIA DataAuditAgent ------------------- The gatekeeper. Always runs FIRST, befor, Main audit pipeline.          context must contain:           - data_dir: path t, Recursively scan directory for all files., Map files to their omics modality., Remove ARIA-generated intermediates from modality inputs. This prevents, Infer genome and organism from filenames, paths, and question., Infer organism from h5ad feature names when metadata is absent. (+34 more)

### Community 8 - "Community 8"
Cohesion: 0.06
Nodes (26): IntegrationAgent, Select and execute the appropriate integration strategy.          context must c, Select the optimal integration strategy.          Decision tree:           same-, WNN integration of scRNA-seq + scATAC-seq from the same cells.         Uses Para, Evaluate WNN modality weights and invoke DebateCouncil         when imbalanced w, Link ATAC peaks to nearby genes using correlation across cells.         Stores h, Store high-confidence peak-to-gene links as Tunnels in ARIAMemory.         Tunne, Invoke DebateCouncil for high-confidence peak-gene links.         Requires 2 ind (+18 more)

### Community 9 - "Community 9"
Cohesion: 0.08
Nodes (33): ARIA scRNAAgent (v4.3.1 — subprocess-only) -------------------------------------, ARIAMemory, LLMProvider, Any, float, int, str, MetricEvaluator (+25 more)

### Community 10 - "Community 10"
Cohesion: 0.06
Nodes (21): GenomeArchAgent, Inspect Hi-C files: format, available resolutions, estimated size.         Uses, Infer file info from extensions when cooler is unavailable., Advise on analysis resolution based on:           1. Biological question (compar, Run topology analysis: compartments + TADs + optional loops.         Uses out-of, Calibrate Insulation Score window_size using chr1 as proxy.          Runs 3 wind, Publish compartment A/B finding.         DebateCouncil review for any A-to-B or, str (+13 more)

### Community 11 - "Community 11"
Cohesion: 0.07
Nodes (44): bool, float, str, bool, str, _basic_chromatin_qc(), _bulk_chromatin_qc(), chromatin_qc() (+36 more)

### Community 12 - "Community 12"
Cohesion: 0.07
Nodes (37): Any, float, int, str, _check_env_file_permissions(), _check_synthetic_assets(), console_main(), main() (+29 more)

### Community 13 - "Community 13"
Cohesion: 0.09
Nodes (33): AuditAgent, _find_star_logs(), _infer_batch_labels(), _infer_group_labels(), _load_count_matrix(), _parse_star_unique_pct(), ARIA AuditAgent (v4.1) — Quality Linter ----------------------------------------, Loads the count matrix, computes log1p(CPM) Pearson correlations, and         co (+25 more)

### Community 14 - "Community 14"
Cohesion: 0.10
Nodes (14): BulkRNAAgent, _normalise_sample_token(), Use the DesignAgent-confirmed groups and factor to build         sample mapping,, P0-6: inferring the experimental design from file/column names is a         gues, Covariates confirmed at DesignAgent CHECKPOINT 2.4, forwarded to the         DES, Return only caller-confirmed contrasts with valid test/ref levels., Suggest candidate contrasts; suggestions do not authorize DE., Infer group labels from column names using the same patterns as rna_bulk_de.py. (+6 more)

### Community 15 - "Community 15"
Cohesion: 0.11
Nodes (22): ARIAMemory, LLMProvider, float, int, Path, str, ContextManager, ModelProfile (+14 more)

### Community 16 - "Community 16"
Cohesion: 0.12
Nodes (39): float, Path, str, adapt(), _entities_from_pb(), _first_mt_threshold(), _genome_for_organism(), _input_record() (+31 more)

### Community 17 - "Community 17"
Cohesion: 0.09
Nodes (20): Any, float, int, str, _LazyMessageBus, MessageBus, Turn on append-only durability at run start (R6). Safe to call on the         pr, Block until an escalation has a user_decision, or timeout expires. (+12 more)

### Community 18 - "Community 18"
Cohesion: 0.13
Nodes (35): Any, bool, float, NarrativeBlock, str, annotate_claim_tiers(), _block_evidence_category(), _block_subject() (+27 more)

### Community 19 - "Community 19"
Cohesion: 0.14
Nodes (32): Raw input ingestion agent for canonical ARIA workspace artifacts., Any, bool, int, Path, str, hash_params(), Order-invariant SHA-256 for JSON-serializable parameter dicts. (+24 more)

### Community 20 - "Community 20"
Cohesion: 0.10
Nodes (33): Any, bool, NarrativeBlock, Any, bool, str, _block_composition_ok(), _block_low_power() (+25 more)

### Community 21 - "Community 21"
Cohesion: 0.11
Nodes (16): DebateCouncil, DebateResult, DebateRound, DebateVerdict, ARIA DebateCouncil ------------------ Internal peer review system for high-stake, Two-agent internal peer review for high-stakes biological interpretations., Run a full debate on a biological interpretation.          Args:             top, Write the final user-facing consensus in normal prose. (+8 more)

### Community 22 - "Community 22"
Cohesion: 0.13
Nodes (9): DesignAgent, Process the user's answer for the current step and advance the         state mac, Parse a free-text manual group assignment from the user. Supports         common, Map user-supplied sample tokens (e.g. 'hc1153') back to the canonical         pa, Reduce a free-text LLM reply to a single short identifier suitable for         a, Parse sample names, merging paired-end read pairs (R1/R2)., Not used in state-machine mode; use start_design() and         handle_user_respo, Kick off the design process. Publishes the first checkpoint (groups)         and (+1 more)

### Community 23 - "Community 23"
Cohesion: 0.09
Nodes (26): bool, Any, bool, Path, str, Apply per-cluster biological labels to an AnnData obs column.  Used by scRNAAgen, rna_apply_cluster_labels(), ARIA RNA Cell-Cell Communication Script ---------------------------------------- (+18 more)

### Community 24 - "Community 24"
Cohesion: 0.14
Nodes (16): Provisions everything ARIA needs before analysis starts.     Runs as Checkpoint, Return path to conda/mamba, installing Miniforge if needed., Download and install Miniforge silently., List installed conda environment names., Which ARIA environments does this experiment need?, Create aria conda env from bundled YAML if not present., Infer genome key from organism name or biological question., Download FASTA + GTF from Ensembl if not in ~/.aria/genomes/.         Decompress (+8 more)

### Community 25 - "Community 25"
Cohesion: 0.15
Nodes (8): OrchestratorAgent, Run AuditAgent synchronously. If blocking issues found, publish CP 3.5         a, CP 3.5 — user decided whether to proceed despite blocking audit issues., A dispatch-gated modality may carry an explicit experimental opt-in         env, Modalities that are dispatch-gated by default but allowed to run only         be, Generate a minimal plan if LLM fails., Message, str

### Community 26 - "Community 26"
Cohesion: 0.11
Nodes (30): float, Path, str, _auto_contrasts(), _generate_plots(), _gseapy_organism(), _gtf_to_symbol_map(), _load_gene_annotation() (+22 more)

### Community 27 - "Community 27"
Cohesion: 0.15
Nodes (28): bool, NarrativeBlock, Path, str, _apply_causal_guard(), _apply_low_confidence_warning(), _apply_trajectory_guard(), collect_named_entities() (+20 more)

### Community 28 - "Community 28"
Cohesion: 0.19
Nodes (13): bool, EvidenceItem, NarrativeBlock, Path, str, Caveat, _design_issues(), _evidence() (+5 more)

### Community 29 - "Community 29"
Cohesion: 0.10
Nodes (25): Any, bool, int, str, Any, float, str, X8 integration QC red-flags + X9 annotation-coherence checks. (+17 more)

### Community 30 - "Community 30"
Cohesion: 0.10
Nodes (18): bool, NarrativeBlock, Path, str, NarrativeBlock, ModalityNarrator, ModalityNarrator, Protocols for modality-specific narrative plugins. (+10 more)

### Community 31 - "Community 31"
Cohesion: 0.14
Nodes (13): DesignIntelligence, format_design_intelligence(), ARIA DesignIntelligence ----------------------- Rules-first feasibility and oppo, Build cross-modality design profiles before computation starts., Human-readable block for checkpoints., bool, int, str (+5 more)

### Community 32 - "Community 32"
Cohesion: 0.19
Nodes (26): bool, int, str, _calibrate_insulation(), _compute_compartments(), _compute_loops(), _compute_tads(), _cooltools_dots() (+18 more)

### Community 33 - "Community 33"
Cohesion: 0.13
Nodes (18): Any, bool, int, Path, str, ContractIssue, test_failed_run_archive_redacts_input_json(), EnvironmentManager (+10 more)

### Community 34 - "Community 34"
Cohesion: 0.12
Nodes (14): bool, int, str, Build the final message list that fits within the model's context.         Appli, Public token counter for external use., 4-step cascade to fit history into token budget.         Each step is logged so, Drop oldest messages until history fits budget.         Always keep the most rec, Apply CavemanULTRA compression to all assistant messages in history.         Use (+6 more)

### Community 35 - "Community 35"
Cohesion: 0.18
Nodes (13): bool, EvidenceItem, NarrativeBlock, Path, str, BulkRnaNarrator, _evidence(), _parse_gsea_row() (+5 more)

### Community 36 - "Community 36"
Cohesion: 0.16
Nodes (20): collect_version_metadata(), _fallback_source_hash(), _git_bytes(), _git_text(), Any, Path, str, Single source of truth for ARIA package version and build stamp. (+12 more)

### Community 37 - "Community 37"
Cohesion: 0.13
Nodes (19): Any, str, ModelConfig, Return the first available model for a tier (for display)., Configuration for a single model endpoint., build_robustness_multiverse(), Deterministic robustness summaries for methodology provenance., Summarize available multiverse checks without hidden reruns.      P-MULTIVERSE o (+11 more)

### Community 38 - "Community 38"
Cohesion: 0.11
Nodes (15): Path, str, float, _json_serializer(), ARIA Analysis Script Base Contract ------------------------------------ Every sc, Write a structured error response to output JSON., Handle numpy and other non-JSON-serializable types.     Called by json.dump as t, Standard entry point for all ARIA analysis scripts.      Reads params from input (+7 more)

### Community 39 - "Community 39"
Cohesion: 0.15
Nodes (10): Publish a user checkpoint for the UI/headless runner to resolve., Publish an in-dispatch checkpoint and wait for the user's decision., Called by MessageBus when this agent receives a message., Main execution method for the agent.         Must return a dict with at minimum:, Core LLM call. Provider-agnostic via LLMProvider.         Caveman compression in, LLM call returning structured JSON. Always ULTRA caveman., float, int (+2 more)

### Community 40 - "Community 40"
Cohesion: 0.38
Nodes (19): NarrativeBlock, str, compose_block_prose(), _compose_cellcomm(), _compose_composition(), _compose_de(), _compose_gsea(), _compose_non_success() (+11 more)

### Community 41 - "Community 41"
Cohesion: 0.12
Nodes (16): bool, _get_gene_sets(), _infer_groups(), _load_or_infer_metadata(), _mock_pathways(), Convert Ensembl IDs to HGNC symbols. Genes without a symbol are     dropped (Enr, Run pathway enrichment via gseapy.     Tests GO Biological Process, KEGG, and Re, Return gene set databases appropriate for organism. (+8 more)

### Community 42 - "Community 42"
Cohesion: 0.17
Nodes (19): bool, float, int, Path, str, _categorical_colors(), _compute_umap(), _embedding_label() (+11 more)

### Community 43 - "Community 43"
Cohesion: 0.19
Nodes (18): bool, int, Path, str, _align_sample(), _build_star_index(), _index_exists(), _mock_alignment() (+10 more)

### Community 44 - "Community 44"
Cohesion: 0.17
Nodes (10): str, Narrative kernel primitives.  This package defines the structured contract used, EvidenceItem, NarrativeBlock, Serializable narrative block schema used by modality narrators., test_render_blocks_shows_claim_evidence_caveats_and_validates_files(), test_block_round_trips_to_dict_for_methodology_json(), test_success_block_requires_claim_and_evidence() (+2 more)

### Community 45 - "Community 45"
Cohesion: 0.16
Nodes (16): float, int, str, export_de_table(), export_pathways_table(), make_gsea_running_sums(), make_ora_dotplot(), _ranked_signature_frame() (+8 more)

### Community 46 - "Community 46"
Cohesion: 0.14
Nodes (17): _abs_corr(), Absolute Pearson correlation, or None when either side has no variance., rna_pseudobulk_de(), _make_collinear_h5ad(), Stage 4 C3: the pseudobulk composition covariate is the cell type's own log-prop, test_abs_corr_detects_collinearity(), test_composition_covariate_dropped_when_collinear(), T1.3 — pseudobulk DE with exactly 2 replicates per condition must     surface lo (+9 more)

### Community 47 - "Community 47"
Cohesion: 0.16
Nodes (16): bool, bool, _cache_matches(), _cache_params(), ARIA RNA Clustering Script -------------------------- Normalization, dimensional, rna_clustering(), _cache_matches(), _cache_params() (+8 more)

### Community 48 - "Community 48"
Cohesion: 0.31
Nodes (14): Any, bool, str, BaseModel, FieldType, ContractField, ContractIssue, _f() (+6 more)

### Community 49 - "Community 49"
Cohesion: 0.16
Nodes (16): bool, int, Path, str, str, _categorical_colors(), _draw_paga(), make_paga_figures() (+8 more)

### Community 50 - "Community 50"
Cohesion: 0.19
Nodes (16): bool, int, Path, str, _detect_samples(), _fastp_outputs_valid(), _mock_fastp_result(), ARIA RNA-seq FASTQ QC Script ------------------------------ Runs fastp (trimming (+8 more)

### Community 51 - "Community 51"
Cohesion: 0.20
Nodes (16): bool, int, Path, str, _build_ensembl_to_symbol_map(), _clean_counts_matrix(), _counts_outputs_valid(), _detect_strandedness() (+8 more)

### Community 52 - "Community 52"
Cohesion: 0.16
Nodes (15): _build_design_formula(), _mock_de_result(), Build a DESeq2 design with covariates first and the factor of interest     last, Keep only covariates usable in this contrast subset: present as a column,     di, Run DESeq2 via pydeseq2 with correct design factor.     Returns (result_dict, wa, Mock DE result for environments without pydeseq2., _resolve_covariates(), _run_deseq2() (+7 more)

### Community 53 - "Community 53"
Cohesion: 0.18
Nodes (16): find_causal_language(), Return the first causal pattern found in ``text``, or None.      Reusable so any, _causal_block(), _ora_block_with_regulatory_term_name(), Causal-language guard hardening (audit 2026-05-28, F-SCI-CAUSAL).  Covers the br, A real causal CLAIM is still downgraded even when evidence has term     names —, Honest ORA block whose enriched term name embeds a causal verb., test_causal_claim_still_caught_despite_term_name_evidence() (+8 more)

### Community 54 - "Community 54"
Cohesion: 0.16
Nodes (16): str, Path, Resolve the grouping obs column from the IPC params.      Canonical key is `grou, _resolve_groupby(), _agent(), P0-1 regression: LIANA must run through the agent path.  The scRNA agent dispatc, Agent-level (not just script): the params the agent dispatches for LIANA     mus, Backward compatibility: a legacy caller sending only `cell_type_col`     (the pr (+8 more)

### Community 55 - "Community 55"
Cohesion: 0.25
Nodes (15): bool, int, str, _balance_cooler(), hic_qc_and_balance(), _mock_hic_qc(), _process_cooler(), _process_hic() (+7 more)

### Community 56 - "Community 56"
Cohesion: 0.37
Nodes (15): NarrativeBlock, Path, str, _claim_tier_badge(), group_blocks_by_prefix(), _group_key(), _image_uri(), HTML composer for structured narrative blocks. (+7 more)

### Community 57 - "Community 57"
Cohesion: 0.19
Nodes (14): AnswerPolicy, drain_pending_checkpoints(), _find_report_on_disk(), HeadlessResult, bool, float, int, str (+6 more)

### Community 58 - "Community 58"
Cohesion: 0.21
Nodes (14): str, Exception, _get_gene_coordinates(), _get_peak_coordinates(), integration_peak2gene(), _load_atac_matrix(), MissingGTFError, _mock_peak2gene() (+6 more)

### Community 59 - "Community 59"
Cohesion: 0.15
Nodes (11): bulk_rna_de(), _compute_tpm(), _contrast_overlap(), _format_top_genes(), _prune_outliers_for_design(), Compute TPM (Transcripts Per Million) from raw counts.      TPM = (reads_per_gen, Keep QC from destroying the statistical design.      Sample-level QC can be nois, Format top genes from a DE result. Adds symbol if mapping available. (+3 more)

### Community 60 - "Community 60"
Cohesion: 0.19
Nodes (13): Any, int, Path, str, test_hash_file_stable_under_chunk_size(), test_provenance_block_contains_required_fields(), collect_llm_usage(), collect_provenance() (+5 more)

### Community 61 - "Community 61"
Cohesion: 0.27
Nodes (12): Tests for the shared raw-count classifier (audit 2026-05-29, P-RAWCLASS).  Cover, R7: a matrix whose first rows are all-zero (e.g. ordered by cell type)     must, _raw_counts(), test_classification_is_order_independent_R7(), test_integer_dtype_small_max_is_raw(), test_lognormalized_matrix_is_not_raw(), test_raw_integer_counts_classified_as_raw(), test_sample_rows_densifies_only_the_slice() (+4 more)

### Community 62 - "Community 62"
Cohesion: 0.23
Nodes (11): bool, int, str, mocks_allowed(), Return True only when simulated analysis output is explicitly enabled.      Prod, integration_mofa(), _load_modality(), _mock_mofa() (+3 more)

### Community 63 - "Community 63"
Cohesion: 0.23
Nodes (11): bool, int, ndarray, test_sampling_is_deterministic_given_seed(), _is_integer_like(), Shared raw-count classifier for the DE entry points.  Audit 2026-05-29, P-RAWCLA, Return up to ``n_rows`` sorted, randomly chosen row indices in [0, n_total)., Densify and return up to ``n_rows`` randomly chosen rows of ``mat``.      Only t (+3 more)

### Community 64 - "Community 64"
Cohesion: 0.18
Nodes (11): int, _plot_pca_mds(), Select the top-N most variable genes from a VST-transformed matrix.      If biot, Generate PCA + MDS plots from a VST-transformed, variable-filtered matrix., Helper — single-axis DR plot for manuscript figures., Sample-level QC for bulk RNA-seq.      Two complementary checks:       (a) Pairw, Variance-Stabilizing Transformation via pydeseq2.      Returns a DataFrame (gene, _run_vst() (+3 more)

### Community 65 - "Community 65"
Cohesion: 0.29
Nodes (10): F-SCI-LOGNORM (audit 2026-05-28): when DESeq2 counts were     reverse-engineered, F-SCI-FDR (audit 2026-05-28): when the pseudobulk script declared the     per-cl, X8/X9: integration-overcorrection and annotation-coherence flags must     surfac, _synthetic_scrna_findings(), test_scrna_narrator_emits_data_quality_block_for_qc_flags(), test_scrna_narrator_flags_lognorm_recovered_counts(), test_scrna_narrator_generates_blocks_for_all_synthetic_results(), test_scrna_narrator_methods_reuse_legacy_methods() (+2 more)

### Community 66 - "Community 66"
Cohesion: 0.22
Nodes (10): float, str, _effective_alpha_from_significant(), _power_disclosure_for_strategy(), ARIA Pseudobulk Differential Expression ----------------------------------------, Find the apeGLM LFC coefficient column for test-vs-ref (C4 shrinkage).      pyde, Return the largest raw p-value among rows passing the applied rule., _shrink_coeff() (+2 more)

### Community 67 - "Community 67"
Cohesion: 0.22
Nodes (10): Any, bool, str, test_privacy_redacts_paths_and_secrets(), air_gapped_enabled(), Privacy helpers for local IPC archives and LLM/cache controls., Return True when ARIA must avoid cloud/network LLM calls., Recursively redact paths/secrets from diagnostic JSON.      Runtime IPC still ne (+2 more)

### Community 69 - "Community 69"
Cohesion: 0.27
Nodes (9): default_answer_policy(), Accept ARIA's high-confidence inference / recommended path.      The choice for, _drive_design(), End-to-end coverage for the design checkpoint state machine.  This closes the F-, Rejecting the final design must cancel, not silently proceed., Run start_design then resolve each published checkpoint via the policy.      Ret, Full CP2.1->2.6 walk on a high-confidence inferred design.      Batch (2.4) and, test_design_state_machine_cancels_on_reject() (+1 more)

### Community 70 - "Community 70"
Cohesion: 0.22
Nodes (9): int, Path, str, FileHandler, attach_experiment_log(), detach_experiment_log(), ARIA per-experiment file logging.  Every `logging.getLogger("aria.X")` is a chil, Attach a FileHandler for this experiment_id to the "aria" root logger.     Retur (+1 more)

### Community 71 - "Community 71"
Cohesion: 0.20
Nodes (10): ARIAMemory, AuditAgent, Claim Compiler, DataAuditAgent, DebateCouncil, DesignAgent, DesignIntelligence, MessageBus (+2 more)

### Community 72 - "Community 72"
Cohesion: 0.29
Nodes (9): float, T1.5 — power should increase with replicate count when dispersion     and effect, test_power_estimate_monotone_in_n(), _as_pair(), bulk_power_estimate(), pseudobulk_power_estimate(), Approximate RNA-seq power estimates for ARIA reports.  These helpers intentional, Approximate two-sided Wald power for a bulk RNA-seq contrast. (+1 more)

### Community 73 - "Community 73"
Cohesion: 0.28
Nodes (8): bool, Path, str, load_aria_env(), _parse_env_line(), ARIA environment loader.  Loads private runtime settings from ~/.aria/.env witho, Load KEY=VALUE pairs from ~/.aria/.env into os.environ.      Args:         path:, Parse one .env line, accepting optional leading 'export'.

### Community 75 - "Community 75"
Cohesion: 0.46
Nodes (7): _load_counts(), Load counts matrix from various formats.     Returns (DataFrame, warnings_list,, Bulk DE raw-count guard (audit 2026-05-29, B10).  `_load_counts` must hard-refus, test_lognorm_matrix_is_hard_refused_by_default(), test_nonraw_matrix_coerced_only_when_allowed(), test_raw_counts_load_and_are_tagged_raw(), _write_matrix()

### Community 76 - "Community 76"
Cohesion: 0.25
Nodes (7): _global_bh(), Return BH-adjusted p-values for one global family of tests., T1.2 — global BH across all block-gene tests can only reduce the     number of c, test_global_fdr_is_more_conservative_than_local_family(), bh_correct(), Small statistical helpers shared by ARIA scripts., Benjamini-Hochberg correction without requiring statsmodels.

### Community 77 - "Community 77"
Cohesion: 0.32
Nodes (7): bool, _cache_matches(), _cache_params(), ARIA RNA QC Script ------------------ Runs quality control on scRNA-seq or bulk, rna_qc(), test_rna_qc_empty_h5ad_returns_structured_error(), test_rna_qc_uses_existing_h5ad_obs_metrics_for_processed_input()

### Community 78 - "Community 78"
Cohesion: 0.29
Nodes (4): test_global_bus_and_env_manager_are_lazy_accessors(), _LazyEnvironmentManager, ARIA EnvironmentManager ----------------------- Isolates bioinformatics tool exe, Create the global EnvironmentManager only when it is first used.

### Community 79 - "Community 79"
Cohesion: 0.29
Nodes (3): P0-3 regression: Hi-C must not dispatch by default.  Hi-C was a `scaffold` modal, The Hi-C opt-in must not unblock unrelated scaffold modalities., test_other_scaffolds_ignore_the_hic_flag()

### Community 80 - "Community 80"
Cohesion: 0.52
Nodes (6): error(), info(), step(), success(), warn(), install.sh script

### Community 81 - "Community 81"
Cohesion: 0.43
Nodes (6): _make_h5ad(), Stage 4 C4: apeGLM LFC shrinkage in pseudobulk DE.  Raw MLE log2 fold changes ov, _run(), test_shrinkage_applied_and_raw_preserved(), test_shrinkage_can_be_disabled(), test_shrinkage_clause_in_narrative()

### Community 82 - "Community 82"
Cohesion: 0.38
Nodes (6): _build_sample_to_group(), inject(), _pick_sample_col(), ARIA scRNA — inject condition obs column from DesignAgent group mapping.  Used b, Invert {group: [samples]} → {sample: group}., Return (sample_col, sample_values_array). The sample column is the     one whose

### Community 83 - "Community 83"
Cohesion: 0.29
Nodes (6): ARIA Differential Abundance Script ----------------------------------- Tests whe, rna_diff_abundance(), T1.1 — when a cell type's abundance doubles in one condition, the     Poisson-of, T1.1 — flat data must yield no significant shifts and the agent's     decision l, test_diff_abundance_detects_2x_shift(), test_diff_abundance_no_signal_returns_none_significant()

### Community 84 - "Community 84"
Cohesion: 0.40
Nodes (4): aria-rna-env, BulkRNAAgent, EnvironmentManager, scRNAAgent

### Community 86 - "Community 86"
Cohesion: 0.60
Nodes (5): _manager(), test_environment_manager_attaches_contract_metadata_on_success(), test_environment_manager_validates_script_output_contract(), test_script_contract_rejects_missing_required_input(), test_script_contract_rejects_version_mismatch()

### Community 87 - "Community 87"
Cohesion: 0.83
Nodes (3): env_installed(), snapshot_one(), generate_locks.sh script

## Knowledge Gaps
- **64 isolated node(s):** `APPLY_fixes.sh script`, `str`, `Any`, `bool`, `Any` (+59 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **9 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Confidence` connect `Community 1` to `Community 0`, `Community 2`, `Community 3`, `Community 68`, `Community 6`, `Community 7`, `Community 8`, `Community 9`, `Community 10`, `Community 39`, `Community 13`, `Community 14`, `Community 17`, `Community 21`, `Community 22`, `Community 24`, `Community 25`?**
  _High betweenness centrality (0.123) - this node is a cross-community bridge._
- **Why does `LLMProvider` connect `Community 15` to `Community 0`, `Community 1`, `Community 2`, `Community 3`, `Community 68`, `Community 37`, `Community 6`, `Community 39`, `Community 8`, `Community 9`, `Community 10`, `Community 74`, `Community 7`, `Community 14`, `Community 21`, `Community 22`, `Community 24`, `Community 25`?**
  _High betweenness centrality (0.122) - this node is a cross-community bridge._
- **Why does `TaskTier` connect `Community 1` to `Community 0`, `Community 2`, `Community 3`, `Community 68`, `Community 37`, `Community 6`, `Community 39`, `Community 8`, `Community 9`, `Community 10`, `Community 74`, `Community 14`, `Community 15`, `Community 16`, `Community 21`, `Community 22`?**
  _High betweenness centrality (0.106) - this node is a cross-community bridge._
- **Are the 89 inferred relationships involving `LLMProvider` (e.g. with `BaseAgent` and `BulkRNAAgent`) actually correct?**
  _`LLMProvider` has 89 INFERRED edges - model-reasoned connections that need verification._
- **Are the 89 inferred relationships involving `Confidence` (e.g. with `AuditAgent` and `BaseAgent`) actually correct?**
  _`Confidence` has 89 INFERRED edges - model-reasoned connections that need verification._
- **Are the 78 inferred relationships involving `BaseAgent` (e.g. with `AuditAgent` and `CavemanMode`) actually correct?**
  _`BaseAgent` has 78 INFERRED edges - model-reasoned connections that need verification._
- **Are the 77 inferred relationships involving `TaskTier` (e.g. with `BaseAgent` and `BulkRNAAgent`) actually correct?**
  _`TaskTier` has 77 INFERRED edges - model-reasoned connections that need verification._