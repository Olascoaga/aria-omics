# Graph Report - /home/medusa/Samael/ARIA  (2026-05-31)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 2446 nodes · 5712 edges · 109 communities (95 shown, 14 thin omitted)
- Extraction: 87% EXTRACTED · 13% INFERRED · 0% AMBIGUOUS · INFERRED: 738 edges (avg confidence: 0.51)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `b91dd265`
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
- [[_COMMUNITY_Community 93|Community 93]]
- [[_COMMUNITY_Community 94|Community 94]]
- [[_COMMUNITY_Community 95|Community 95]]
- [[_COMMUNITY_Community 96|Community 96]]
- [[_COMMUNITY_Community 97|Community 97]]
- [[_COMMUNITY_Community 98|Community 98]]
- [[_COMMUNITY_Community 99|Community 99]]
- [[_COMMUNITY_Community 100|Community 100]]
- [[_COMMUNITY_Community 102|Community 102]]
- [[_COMMUNITY_Community 103|Community 103]]
- [[_COMMUNITY_Community 104|Community 104]]
- [[_COMMUNITY_Community 105|Community 105]]
- [[_COMMUNITY_Community 106|Community 106]]
- [[_COMMUNITY_Community 107|Community 107]]

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
- `EnvironmentManager` --references--> `aria-chromatin-env`  [INFERRED]
  aria/utils/environment_manager.py → envs/aria-chromatin-env.yml
- `EnvironmentManager` --references--> `aria-hic-env`  [INFERRED]
  aria/utils/environment_manager.py → envs/aria-hic-env.yml
- `EnvironmentManager` --references--> `aria-integration-env`  [INFERRED]
  aria/utils/environment_manager.py → envs/aria-integration-env.yml
- `EnvironmentManager` --references--> `aria-rna-env`  [INFERRED]
  aria/utils/environment_manager.py → envs/aria-rna-env.yml
- `test_shrinkage_clause_in_narrative()` --calls--> `_lfc_shrinkage_clause()`  [EXTRACTED]
  tests/test_lfc_shrinkage.py → aria/agents/_narrative_scrna.py

## Import Cycles
- None detected.

## Communities (109 total, 14 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.05
Nodes (35): Robust JSON extraction: handles ```json fences, leading prose, and         bare, Run ORA per Leiden cluster against GO_BP / KEGG / Reactome via the         rna_p, Return genes detected at least once in the retained h5ad., True when DesignAgent identified enough biological groups (the         prerequis, Write a copy of `h5ad_in` with `obs[main_factor]` populated from         the sam, Pseudobulk DE between condition groups identified by DesignAgent.          Aggre, If the user asks to focus on specific obs cell types, materialize a         focu, Derive a stable per-sample label from a 10x .h5 / MEX / .h5ad path.         Stri (+27 more)

### Community 1 - "Community 1"
Cohesion: 0.03
Nodes (38): object, _assert_legacy_script_passed(), T1.3 — CP2 recommended-only must not run n=2 pseudobulk DE.     Design Intellige, T1.3 — n=2 pseudobulk remains supported with warnings when CP2     selected reco, T1.1 — scrna_agent._run_pseudobulk must pass     composition_covariate=True to r, T1.4 — scRNA ORA receives the dataset-expressed gene universe, not     only the, Generalization regression — the public Kang et al. PBMC IFN-β     dataset (8 don, When the LLM hallucinates a sentence as a 'factor name', the CP     2.3 menu mus (+30 more)

### Community 2 - "Community 2"
Cohesion: 0.05
Nodes (29): NarrativeAgent, Honest accounting of conflicts and limitations., Construct (and mkdir) ~/.aria/reports/aria_<ts>_<slug>_<suffix>/         with fi, Orchestrate scRNA figure rendering (UMAPs + DE bar + pathway dotplots         +, Render the full HTML report.         Self-contained: CSS embedded, no external d, One-paragraph executive summary for the PI.         Uses deterministic structure, P-LEDGER: render the planned-vs-run manifest. Any analysis the plan         call, Build a short URL-safe slug from biological entities or the question.         Ex (+21 more)

### Community 3 - "Community 3"
Cohesion: 0.08
Nodes (59): ABC, BaseAgent, ARIA BaseAgent -------------- All ARIA agents inherit from this class. Handles:, Publish a user checkpoint for the UI/headless runner to resolve., Publish an in-dispatch checkpoint and wait for the user's decision., Called by MessageBus when this agent receives a message., Main execution method for the agent.         Must return a dict with at minimum:, Base class for all ARIA agents.      Uses LLMProvider for all LLM calls — comple (+51 more)

### Community 4 - "Community 4"
Cohesion: 0.06
Nodes (41): bool, int, str, float, int, Path, str, ContextManager (+33 more)

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
Nodes (47): Any, bool, float, int, str, _check_env_file_permissions(), _check_synthetic_assets(), console_main() (+39 more)

### Community 9 - "Community 9"
Cohesion: 0.07
Nodes (44): bool, float, str, bool, str, _basic_chromatin_qc(), _bulk_chromatin_qc(), chromatin_qc() (+36 more)

### Community 10 - "Community 10"
Cohesion: 0.09
Nodes (28): ARIAMemory, LLMProvider, float, int, str, MetricEvaluator, ParameterAdvisor, ParameterCandidate (+20 more)

### Community 11 - "Community 11"
Cohesion: 0.06
Nodes (27): str, str, Exception, _get_gene_coordinates(), _get_peak_coordinates(), integration_peak2gene(), _load_atac_matrix(), MissingGTFError (+19 more)

### Community 12 - "Community 12"
Cohesion: 0.10
Nodes (14): BulkRNAAgent, _normalise_sample_token(), Use the DesignAgent-confirmed groups and factor to build         sample mapping,, P0-6: inferring the experimental design from file/column names is a         gues, Covariates confirmed at DesignAgent CHECKPOINT 2.4, forwarded to the         DES, Return only caller-confirmed contrasts with valid test/ref levels., Suggest candidate contrasts; suggestions do not authorize DE., Infer group labels from column names using the same patterns as rna_bulk_de.py. (+6 more)

### Community 13 - "Community 13"
Cohesion: 0.09
Nodes (18): DesignAgent, Process the user's answer for the current step and advance the         state mac, Parse a free-text manual group assignment from the user. Supports         common, Map user-supplied sample tokens (e.g. 'hc1153') back to the canonical         pa, Reduce a free-text LLM reply to a single short identifier suitable for         a, Parse sample names, merging paired-end read pairs (R1/R2)., Not used in state-machine mode; use start_design() and         handle_user_respo, Kick off the design process. Publishes the first checkpoint (groups)         and (+10 more)

### Community 14 - "Community 14"
Cohesion: 0.13
Nodes (23): bool, EvidenceItem, NarrativeBlock, Path, str, Caveat, _design_issues(), _evidence() (+15 more)

### Community 15 - "Community 15"
Cohesion: 0.12
Nodes (39): float, Path, str, adapt(), _entities_from_pb(), _first_mt_threshold(), _genome_for_organism(), _input_record() (+31 more)

### Community 16 - "Community 16"
Cohesion: 0.09
Nodes (30): AuditAgent, _find_star_logs(), _infer_batch_labels(), _infer_group_labels(), _load_count_matrix(), ARIA AuditAgent (v4.1) — Quality Linter ----------------------------------------, Loads the count matrix, computes log1p(CPM) Pearson correlations, and         co, Runs a quick PCA on log1p(CPM) counts and checks whether PC1 separates         s (+22 more)

### Community 17 - "Community 17"
Cohesion: 0.08
Nodes (21): Any, float, int, str, _LazyMessageBus, MessageBus, Turn on append-only durability at run start (R6). Safe to call on the         pr, Block until an escalation has a user_decision, or timeout expires. (+13 more)

### Community 18 - "Community 18"
Cohesion: 0.13
Nodes (35): Any, bool, float, NarrativeBlock, str, annotate_claim_tiers(), _block_evidence_category(), _block_subject() (+27 more)

### Community 19 - "Community 19"
Cohesion: 0.15
Nodes (32): Raw input ingestion agent for canonical ARIA workspace artifacts., Any, bool, int, Path, str, hash_params(), Order-invariant SHA-256 for JSON-serializable parameter dicts. (+24 more)

### Community 20 - "Community 20"
Cohesion: 0.08
Nodes (29): bool, bool, Any, bool, Path, str, Apply per-cluster biological labels to an AnnData obs column.  Used by scRNAAgen, rna_apply_cluster_labels() (+21 more)

### Community 21 - "Community 21"
Cohesion: 0.11
Nodes (16): DebateCouncil, DebateResult, DebateRound, DebateVerdict, ARIA DebateCouncil ------------------ Internal peer review system for high-stake, Two-agent internal peer review for high-stakes biological interpretations., Run a full debate on a biological interpretation.          Args:             top, Write the final user-facing consensus in normal prose. (+8 more)

### Community 22 - "Community 22"
Cohesion: 0.10
Nodes (34): float, Path, str, _auto_contrasts(), _generate_plots(), _gtf_to_symbol_map(), _load_gene_annotation(), _load_symbol_map() (+26 more)

### Community 23 - "Community 23"
Cohesion: 0.10
Nodes (19): bool, NarrativeBlock, Path, str, NarrativeBlock, ModalityNarrator, ModalityNarrator, Protocols for modality-specific narrative plugins. (+11 more)

### Community 24 - "Community 24"
Cohesion: 0.08
Nodes (13): str, hic_inspect(), _inspect_cool(), _inspect_file(), _inspect_hic(), _inspect_mcool(), ARIA Hi-C File Inspector ------------------------- Reads Hi-C file metadata WITH, Inspect a single Hi-C file — metadata only. (+5 more)

### Community 25 - "Community 25"
Cohesion: 0.14
Nodes (16): Provisions everything ARIA needs before analysis starts.     Runs as Checkpoint, Return path to conda/mamba, installing Miniforge if needed., Download and install Miniforge silently., List installed conda environment names., Which ARIA environments does this experiment need?, Create aria conda env from bundled YAML if not present., Infer genome key from organism name or biological question., Download FASTA + GTF from Ensembl if not in ~/.aria/genomes/.         Decompress (+8 more)

### Community 26 - "Community 26"
Cohesion: 0.12
Nodes (16): IntegrationAgent, Select and execute the appropriate integration strategy.          context must c, Select the optimal integration strategy.          Decision tree:           same-, WNN integration of scRNA-seq + scATAC-seq from the same cells.         Uses Para, Evaluate WNN modality weights and invoke DebateCouncil         when imbalanced w, Link ATAC peaks to nearby genes using correlation across cells.         Stores h, Store high-confidence peak-to-gene links as Tunnels in ARIAMemory.         Tunne, Invoke DebateCouncil for high-confidence peak-gene links.         Requires 2 ind (+8 more)

### Community 27 - "Community 27"
Cohesion: 0.15
Nodes (30): ask_biological_question(), _discard_queued_stdin_lines(), _drain_checkpoints(), _live_analysis_loop(), main(), print_agent_message(), print_banner(), print_checkpoint() (+22 more)

### Community 28 - "Community 28"
Cohesion: 0.15
Nodes (28): bool, NarrativeBlock, Path, str, _apply_causal_guard(), _apply_low_confidence_warning(), _apply_trajectory_guard(), collect_named_entities() (+20 more)

### Community 29 - "Community 29"
Cohesion: 0.15
Nodes (8): OrchestratorAgent, Run AuditAgent synchronously. If blocking issues found, publish CP 3.5         a, CP 3.5 — user decided whether to proceed despite blocking audit issues., A dispatch-gated modality may carry an explicit experimental opt-in         env, Modalities that are dispatch-gated by default but allowed to run only         be, Generate a minimal plan if LLM fails., Message, str

### Community 30 - "Community 30"
Cohesion: 0.11
Nodes (27): str, Any, bool, str, _get_gene_sets(), _gseapy_organism(), ARIA RNA Per-Cluster Pathway Enrichment ----------------------------------------, gseapy 1.1.13+ requires lowercase organism strings. (+19 more)

### Community 31 - "Community 31"
Cohesion: 0.10
Nodes (25): Any, bool, int, str, Any, float, str, X8 integration QC red-flags + X9 annotation-coherence checks. (+17 more)

### Community 32 - "Community 32"
Cohesion: 0.09
Nodes (20): _infer_lfc_threshold(), _is_fastq(), ARIA BulkRNAAgent (v3.10 → v4.0) --------------------------------- Bulk RNA-seq, DesignStep, ARIA DesignAgent (v4.0) ------------------------ Interactive experimental design, ARIA IntegrationAgent --------------------- Synthesizes findings across omics mo, ARIA scRNAAgent (v4.3.1 — subprocess-only) -------------------------------------, bool (+12 more)

### Community 33 - "Community 33"
Cohesion: 0.14
Nodes (12): GenomeArchAgent, ARIA GenomeArchAgent -------------------- 3D genome organization analysis: TADs,, Inspect Hi-C files: format, available resolutions, estimated size.         Uses, Infer file info from extensions when cooler is unavailable., Advise on analysis resolution based on:           1. Biological question (compar, Run topology analysis: compartments + TADs + optional loops.         Uses out-of, Calibrate Insulation Score window_size using chr1 as proxy.          Runs 3 wind, Publish compartment A/B finding.         DebateCouncil review for any A-to-B or (+4 more)

### Community 34 - "Community 34"
Cohesion: 0.16
Nodes (20): Path, str, _cb(), GEOConnector, _gse_prefix(), _infer_data_type(), _parse_soft_text(), ARIA GEO/SRA Connector (v4.3) ------------------------------ Downloads and prepr (+12 more)

### Community 35 - "Community 35"
Cohesion: 0.13
Nodes (18): Any, bool, int, Path, str, ContractIssue, test_failed_run_archive_redacts_input_json(), EnvironmentManager (+10 more)

### Community 36 - "Community 36"
Cohesion: 0.19
Nodes (26): bool, int, str, _calibrate_insulation(), _compute_compartments(), _compute_loops(), _compute_tads(), _cooltools_dots() (+18 more)

### Community 37 - "Community 37"
Cohesion: 0.18
Nodes (13): bool, EvidenceItem, NarrativeBlock, Path, str, BulkRnaNarrator, _evidence(), _parse_gsea_row() (+5 more)

### Community 38 - "Community 38"
Cohesion: 0.10
Nodes (17): Any, str, Narrative kernel primitives.  This package defines the structured contract used, build_robustness_multiverse(), Deterministic robustness summaries for methodology provenance., Summarize available multiverse checks without hidden reruns.      P-MULTIVERSE o, EvidenceItem, Serializable narrative block schema used by modality narrators. (+9 more)

### Community 39 - "Community 39"
Cohesion: 0.16
Nodes (20): collect_version_metadata(), _fallback_source_hash(), _git_bytes(), _git_text(), Any, Path, str, Single source of truth for ARIA package version and build stamp. (+12 more)

### Community 40 - "Community 40"
Cohesion: 0.15
Nodes (23): bool, int, ndarray, Tests for the shared raw-count classifier (audit 2026-05-29, P-RAWCLASS).  Cover, R7: a matrix whose first rows are all-zero (e.g. ordered by cell type)     must, _raw_counts(), test_classification_is_order_independent_R7(), test_integer_dtype_small_max_is_raw() (+15 more)

### Community 41 - "Community 41"
Cohesion: 0.10
Nodes (18): Path, str, bool, _json_serializer(), ARIA Analysis Script Base Contract ------------------------------------ Every sc, Write a structured error response to output JSON., Handle numpy and other non-JSON-serializable types.     Called by json.dump as t, Standard entry point for all ARIA analysis scripts.      Reads params from input (+10 more)

### Community 42 - "Community 42"
Cohesion: 0.09
Nodes (20): bool, _get_gene_sets(), _gseapy_organism(), _infer_groups(), _load_or_infer_metadata(), _mock_pathways(), Convert Ensembl IDs to HGNC symbols. Genes without a symbol are     dropped (Enr, Run pathway enrichment via gseapy.     Tests GO Biological Process, KEGG, and Re (+12 more)

### Community 43 - "Community 43"
Cohesion: 0.14
Nodes (22): Any, bool, str, build_devils_advocate(), Annotate each challengeable block in place with a devil's-advocate caveat     an, build_run_ledger(), _is_planned(), _plan_phrases() (+14 more)

### Community 44 - "Community 44"
Cohesion: 0.18
Nodes (9): DesignIntelligence, Build cross-modality design profiles before computation starts., bool, str, T1.3 — Design Intelligence must move pseudobulk DE from     'recommended' to 'op, T1.3 — same downgrade rule applies to bulk RNA DE., test_design_intelligence_downgrades_bulk_at_n2(), test_design_intelligence_downgrades_pseudobulk_at_n2() (+1 more)

### Community 45 - "Community 45"
Cohesion: 0.13
Nodes (20): _build_design_formula(), _mock_de_result(), Build a DESeq2 design with covariates first and the factor of interest     last, Keep only covariates usable in this contrast subset: present as a column,     di, Run DESeq2 via pydeseq2 with correct design factor.     Returns (result_dict, wa, Mock DE result for environments without pydeseq2., _resolve_covariates(), _run_deseq2() (+12 more)

### Community 46 - "Community 46"
Cohesion: 0.38
Nodes (19): NarrativeBlock, str, compose_block_prose(), _compose_cellcomm(), _compose_composition(), _compose_de(), _compose_gsea(), _compose_non_success() (+11 more)

### Community 47 - "Community 47"
Cohesion: 0.17
Nodes (19): bool, float, int, Path, str, _categorical_colors(), _compute_umap(), _embedding_label() (+11 more)

### Community 48 - "Community 48"
Cohesion: 0.19
Nodes (18): bool, int, Path, str, _align_sample(), _build_star_index(), _index_exists(), _mock_alignment() (+10 more)

### Community 49 - "Community 49"
Cohesion: 0.16
Nodes (16): float, int, str, export_de_table(), export_pathways_table(), make_gsea_running_sums(), make_ora_dotplot(), _ranked_signature_frame() (+8 more)

### Community 50 - "Community 50"
Cohesion: 0.18
Nodes (16): float, P1-1 (c): FDR across the bulk contrast family.  Bulk DE controlled FDR within ea, test_bulk_de_global_contrast_family_is_recorded(), test_contrast_family_can_skip_lfc_gate_when_wald_test_used_lfc_null(), test_contrast_family_significance_pools_and_gates_on_lfc(), test_pooled_bh_is_more_conservative_than_per_group(), test_preregister_contrast_family_normalizes_and_declares(), bh_correct() (+8 more)

### Community 51 - "Community 51"
Cohesion: 0.31
Nodes (14): Any, bool, str, BaseModel, FieldType, ContractField, ContractIssue, _f() (+6 more)

### Community 52 - "Community 52"
Cohesion: 0.18
Nodes (17): AST, FunctionDef, _call_name(), _module_gates_mocks(), bool, Call, int, str (+9 more)

### Community 53 - "Community 53"
Cohesion: 0.14
Nodes (17): _abs_corr(), Absolute Pearson correlation, or None when either side has no variance., rna_pseudobulk_de(), _make_collinear_h5ad(), Stage 4 C3: the pseudobulk composition covariate is the cell type's own log-prop, test_abs_corr_detects_collinearity(), test_composition_covariate_dropped_when_collinear(), T1.3 — pseudobulk DE with exactly 2 replicates per condition must     surface lo (+9 more)

### Community 54 - "Community 54"
Cohesion: 0.16
Nodes (16): str, Path, Resolve the grouping obs column from the IPC params.      Canonical key is `grou, _resolve_groupby(), _agent(), P0-1 regression: LIANA must run through the agent path.  The scRNA agent dispatc, Agent-level (not just script): the params the agent dispatches for LIANA     mus, Backward compatibility: a legacy caller sending only `cell_type_col`     (the pr (+8 more)

### Community 55 - "Community 55"
Cohesion: 0.18
Nodes (16): find_causal_language(), Return the first causal pattern found in ``text``, or None.      Reusable so any, _causal_block(), _ora_block_with_regulatory_term_name(), Causal-language guard hardening (audit 2026-05-28, F-SCI-CAUSAL).  Covers the br, A real causal CLAIM is still downgraded even when evidence has term     names —, Honest ORA block whose enriched term name embeds a causal verb., test_causal_claim_still_caught_despite_term_name_evidence() (+8 more)

### Community 56 - "Community 56"
Cohesion: 0.19
Nodes (16): bool, int, Path, str, _detect_samples(), _fastp_outputs_valid(), _mock_fastp_result(), ARIA RNA-seq FASTQ QC Script ------------------------------ Runs fastp (trimming (+8 more)

### Community 57 - "Community 57"
Cohesion: 0.20
Nodes (16): bool, int, Path, str, _build_ensembl_to_symbol_map(), _clean_counts_matrix(), _counts_outputs_valid(), _detect_strandedness() (+8 more)

### Community 58 - "Community 58"
Cohesion: 0.37
Nodes (15): NarrativeBlock, Path, str, _claim_tier_badge(), group_blocks_by_prefix(), _group_key(), _image_uri(), HTML composer for structured narrative blocks. (+7 more)

### Community 59 - "Community 59"
Cohesion: 0.19
Nodes (14): AnswerPolicy, drain_pending_checkpoints(), _find_report_on_disk(), HeadlessResult, bool, float, int, str (+6 more)

### Community 60 - "Community 60"
Cohesion: 0.27
Nodes (14): bool, int, str, _balance_cooler(), _mock_hic_qc(), _process_cooler(), _process_hic(), _process_single_file() (+6 more)

### Community 61 - "Community 61"
Cohesion: 0.19
Nodes (14): str, _make_h5ad(), P1-2: the multiple-testing family is pre-registered, not chosen post-hoc.  Five, _run(), test_post_hoc_switch_is_rejected(), test_preregister_normalizes_and_declares(), test_primary_column_depends_only_on_strategy(), test_pseudobulk_records_preregistered_fdr_family() (+6 more)

### Community 62 - "Community 62"
Cohesion: 0.19
Nodes (13): Any, int, Path, str, test_hash_file_stable_under_chunk_size(), test_provenance_block_contains_required_fields(), collect_llm_usage(), collect_provenance() (+5 more)

### Community 63 - "Community 63"
Cohesion: 0.18
Nodes (10): BulkRNAAgent, aria-chromatin-env, aria-hic-env, aria-integration-env, aria-rna-env, EnvironmentManager, ParameterAdvisor, rna_cellcomm.py (LIANA) (+2 more)

### Community 64 - "Community 64"
Cohesion: 0.31
Nodes (12): Any, bool, NarrativeBlock, str, _block_composition_ok(), _block_low_power(), _challenges_for(), _contains() (+4 more)

### Community 65 - "Community 65"
Cohesion: 0.21
Nodes (12): bool, int, str, mocks_allowed(), Return True only when simulated analysis output is explicitly enabled.      Prod, hic_qc_and_balance(), integration_mofa(), _load_modality() (+4 more)

### Community 66 - "Community 66"
Cohesion: 0.18
Nodes (9): bulk_rna_de(), _compute_tpm(), _contrast_overlap(), _format_top_genes(), Compute TPM (Transcripts Per Million) from raw counts.      TPM = (reads_per_gen, Format top genes from a DE result. Adds symbol if mapping available., Compute DE gene overlap between contrasts.     Uses the FULL list of significant, P0-5 regression: DE must not choose reference levels implicitly.  Bulk and pseud (+1 more)

### Community 67 - "Community 67"
Cohesion: 0.24
Nodes (11): bool, _cache_params(), _cache_matches(), _cache_params(), ARIA RNA Concatenation Script ------------------------------ Concatenates a list, rna_concat(), test_clustering_cache_requires_matching_parameters(), test_clustering_skip_leiden_when_cluster_col_provided() (+3 more)

### Community 68 - "Community 68"
Cohesion: 0.33
Nodes (8): _organism_from_gene_symbols(), Peek at gene IDs in a count matrix (first column, first 200 rows) to     infer o, Deferred WIP closeout: GEO multi-organism (spike-in) handling.  Two pure helpers, test_organism_from_ensembl_prefixes(), test_organism_from_human_hgnc_symbols(), test_organism_from_mouse_mgi_symbols(), test_organism_inference_is_uncertain_for_mixed_ids(), _write_counts()

### Community 69 - "Community 69"
Cohesion: 0.22
Nodes (10): float, str, _effective_alpha_from_significant(), _power_disclosure_for_strategy(), ARIA Pseudobulk Differential Expression ----------------------------------------, Find the apeGLM LFC coefficient column for test-vs-ref (C4 shrinkage).      pyde, Return the largest raw p-value among rows passing the applied rule., _shrink_coeff() (+2 more)

### Community 70 - "Community 70"
Cohesion: 0.27
Nodes (10): bool, int, Path, str, _categorical_colors(), _draw_paga(), make_paga_figures(), ARIA scRNA PAGA + DPT figure generator.  Reads a trajectory.h5ad (from rna_traje (+2 more)

### Community 72 - "Community 72"
Cohesion: 0.22
Nodes (9): int, Path, str, FileHandler, attach_experiment_log(), detach_experiment_log(), ARIA per-experiment file logging.  Every `logging.getLogger("aria.X")` is a chil, Attach a FileHandler for this experiment_id to the "aria" root logger.     Retur (+1 more)

### Community 73 - "Community 73"
Cohesion: 0.29
Nodes (9): float, T1.5 — power should increase with replicate count when dispersion     and effect, test_power_estimate_monotone_in_n(), _as_pair(), bulk_power_estimate(), pseudobulk_power_estimate(), Approximate RNA-seq power estimates for ARIA reports.  These helpers intentional, Approximate two-sided Wald power for a bulk RNA-seq contrast. (+1 more)

### Community 74 - "Community 74"
Cohesion: 0.28
Nodes (8): bool, Path, str, load_aria_env(), _parse_env_line(), ARIA environment loader.  Loads private runtime settings from ~/.aria/.env witho, Load KEY=VALUE pairs from ~/.aria/.env into os.environ.      Args:         path:, Parse one .env line, accepting optional leading 'export'.

### Community 75 - "Community 75"
Cohesion: 0.31
Nodes (8): _discover_dispatched_scripts(), Call, str, P1-12: every dispatchable script must declare an IPC contract.  EnvironmentManag, Extract the script path from a run_in_stack(...) call.      Supports both the ke, _script_path_of(), test_contract_lookup_resolves_for_each_contract_key(), test_every_dispatched_script_on_disk_has_a_contract()

### Community 77 - "Community 77"
Cohesion: 0.46
Nodes (7): _load_counts(), Load counts matrix from various formats.     Returns (DataFrame, warnings_list,, Bulk DE raw-count guard (audit 2026-05-29, B10).  `_load_counts` must hard-refus, test_lognorm_matrix_is_hard_refused_by_default(), test_nonraw_matrix_coerced_only_when_allowed(), test_raw_counts_load_and_are_tagged_raw(), _write_matrix()

### Community 78 - "Community 78"
Cohesion: 0.29
Nodes (3): P0-3 regression: Hi-C must not dispatch by default.  Hi-C was a `scaffold` modal, The Hi-C opt-in must not unblock unrelated scaffold modalities., test_other_scaffolds_ignore_the_hic_flag()

### Community 79 - "Community 79"
Cohesion: 0.29
Nodes (7): int, Select the top-N most variable genes from a VST-transformed matrix.      If biot, Sample-level QC for bulk RNA-seq.      Two complementary checks:       (a) Pairw, Variance-Stabilizing Transformation via pydeseq2.      Returns a DataFrame (gene, _run_vst(), _sample_qc(), _select_variable_genes()

### Community 80 - "Community 80"
Cohesion: 0.33
Nodes (3): _LazyEnvironmentManager, ARIA EnvironmentManager ----------------------- Isolates bioinformatics tool exe, Create the global EnvironmentManager only when it is first used.

### Community 81 - "Community 81"
Cohesion: 0.43
Nodes (6): _make_h5ad(), Stage 4 C4: apeGLM LFC shrinkage in pseudobulk DE.  Raw MLE log2 fold changes ov, _run(), test_shrinkage_applied_and_raw_preserved(), test_shrinkage_can_be_disabled(), test_shrinkage_clause_in_narrative()

### Community 82 - "Community 82"
Cohesion: 0.52
Nodes (6): error(), info(), step(), success(), warn(), install.sh script

### Community 83 - "Community 83"
Cohesion: 0.38
Nodes (6): _build_sample_to_group(), inject(), _pick_sample_col(), ARIA scRNA — inject condition obs column from DesignAgent group mapping.  Used b, Invert {group: [samples]} → {sample: group}., Return (sample_col, sample_values_array). The sample column is the     one whose

### Community 84 - "Community 84"
Cohesion: 0.29
Nodes (6): ARIA Differential Abundance Script ----------------------------------- Tests whe, rna_diff_abundance(), T1.1 — when a cell type's abundance doubles in one condition, the     Poisson-of, T1.1 — flat data must yield no significant shifts and the agent's     decision l, test_diff_abundance_detects_2x_shift(), test_diff_abundance_no_signal_returns_none_significant()

### Community 85 - "Community 85"
Cohesion: 0.33
Nodes (4): format_design_intelligence(), ARIA DesignIntelligence ----------------------- Rules-first feasibility and oppo, Human-readable block for checkpoints., int

### Community 86 - "Community 86"
Cohesion: 0.60
Nodes (5): _manager(), test_environment_manager_attaches_contract_metadata_on_success(), test_environment_manager_validates_script_output_contract(), test_script_contract_rejects_missing_required_input(), test_script_contract_rejects_version_mismatch()

### Community 87 - "Community 87"
Cohesion: 0.40
Nodes (6): ARIAMemory, DataAuditAgent, DesignAgent, DesignIntelligence, MessageBus, OrchestratorAgent

### Community 89 - "Community 89"
Cohesion: 0.40
Nodes (5): float, _graph_modularity(), ARIA Leiden Resolution Advisor Script -------------------------------------- Eva, Compute Leiden cluster modularity from Scanpy's neighbor connectivities., rna_advise_resolution()

### Community 90 - "Community 90"
Cohesion: 0.40
Nodes (5): print_agent_progress(), print_finding(), Message, Display a STATUS message from an agent with clear identification., Display a FINDING from an agent with confidence badge.

### Community 91 - "Community 91"
Cohesion: 0.50
Nodes (3): _infer_design(), Infer experimental groups from sample characteristics.     Prefer experimental-d, test_geo_design_prefers_experimental_keys_and_uses_geo_ids()

### Community 92 - "Community 92"
Cohesion: 0.50
Nodes (4): _parse_star_unique_pct(), Parse 'Uniquely mapped reads %' from a STAR Log.final.out file., float, Path

### Community 93 - "Community 93"
Cohesion: 0.50
Nodes (4): Claim Compiler, DebateCouncil / Devil's Advocate, NarrativeAgent, Narrative Kernel

### Community 94 - "Community 94"
Cohesion: 0.83
Nodes (3): env_installed(), snapshot_one(), generate_locks.sh script

### Community 95 - "Community 95"
Cohesion: 0.50
Nodes (4): _global_bh(), Return BH-adjusted p-values for one global family of tests., T1.2 — global BH across all block-gene tests can only reduce the     number of c, test_global_fdr_is_more_conservative_than_local_family()

## Knowledge Gaps
- **81 isolated node(s):** `APPLY_fixes.sh script`, `str`, `Any`, `bool`, `Any` (+76 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **14 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Confidence` connect `Community 3` to `Community 0`, `Community 1`, `Community 2`, `Community 6`, `Community 7`, `Community 10`, `Community 11`, `Community 12`, `Community 13`, `Community 16`, `Community 17`, `Community 21`, `Community 25`, `Community 26`, `Community 29`, `Community 32`, `Community 33`, `Community 71`, `Community 92`?**
  _High betweenness centrality (0.139) - this node is a cross-community bridge._
- **Why does `LLMProvider` connect `Community 3` to `Community 32`, `Community 33`, `Community 2`, `Community 0`, `Community 4`, `Community 1`, `Community 6`, `Community 71`, `Community 7`, `Community 38`, `Community 10`, `Community 12`, `Community 13`, `Community 76`, `Community 21`, `Community 25`, `Community 26`, `Community 29`?**
  _High betweenness centrality (0.107) - this node is a cross-community bridge._
- **Why does `BaseAgent` connect `Community 3` to `Community 32`, `Community 33`, `Community 2`, `Community 0`, `Community 1`, `Community 6`, `Community 7`, `Community 8`, `Community 10`, `Community 12`, `Community 13`, `Community 16`, `Community 19`, `Community 25`, `Community 26`, `Community 92`, `Community 29`?**
  _High betweenness centrality (0.090) - this node is a cross-community bridge._
- **Are the 89 inferred relationships involving `LLMProvider` (e.g. with `BaseAgent` and `BulkRNAAgent`) actually correct?**
  _`LLMProvider` has 89 INFERRED edges - model-reasoned connections that need verification._
- **Are the 89 inferred relationships involving `Confidence` (e.g. with `AuditAgent` and `BaseAgent`) actually correct?**
  _`Confidence` has 89 INFERRED edges - model-reasoned connections that need verification._
- **Are the 78 inferred relationships involving `BaseAgent` (e.g. with `AuditAgent` and `CavemanMode`) actually correct?**
  _`BaseAgent` has 78 INFERRED edges - model-reasoned connections that need verification._
- **Are the 77 inferred relationships involving `TaskTier` (e.g. with `BaseAgent` and `BulkRNAAgent`) actually correct?**
  _`TaskTier` has 77 INFERRED edges - model-reasoned connections that need verification._