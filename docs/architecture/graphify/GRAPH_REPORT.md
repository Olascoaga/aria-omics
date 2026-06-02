# Graph Report - /home/medusa/Samael/ARIA  (2026-06-01)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 1995 nodes · 4693 edges · 89 communities (81 shown, 8 thin omitted)
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `41dfbea5`
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
- [[_COMMUNITY_Community 77|Community 77]]
- [[_COMMUNITY_Community 78|Community 78]]
- [[_COMMUNITY_Community 79|Community 79]]
- [[_COMMUNITY_Community 80|Community 80]]
- [[_COMMUNITY_Community 81|Community 81]]
- [[_COMMUNITY_Community 82|Community 82]]
- [[_COMMUNITY_Community 84|Community 84]]

## God Nodes (most connected - your core abstractions)
1. `scRNAAgent` - 45 edges
2. `BaseAgent` - 40 edges
3. `NarrativeAgent` - 39 edges
4. `LLMProvider` - 38 edges
5. `read_h5ad()` - 37 edges
6. `BulkRNAAgent` - 36 edges
7. `NarrativeBlock` - 35 edges
8. `ParameterAdvisor` - 34 edges
9. `DesignAgent` - 32 edges
10. `OrchestratorAgent` - 32 edges

## Surprising Connections (you probably didn't know these)
- `test_shrinkage_clause_in_narrative()` --calls--> `_lfc_shrinkage_clause()`  [EXTRACTED]
  tests/test_lfc_shrinkage.py → aria/agents/_narrative_scrna.py
- `test_bulk_methods_disclose_a_dropped_covariate()` --calls--> `BulkRnaNarrator`  [EXTRACTED]
  tests/test_bulk_covariates.py → aria/agents/narrative/narrators/bulk_rna.py
- `test_bulk_methods_report_the_fitted_covariate_formula()` --calls--> `BulkRnaNarrator`  [EXTRACTED]
  tests/test_bulk_covariates.py → aria/agents/narrative/narrators/bulk_rna.py
- `test_timeout_env_override()` --calls--> `LLMProvider`  [EXTRACTED]
  tests/test_llm_reliability.py → aria/llm/provider.py
- `test_golden_bulk_de_recovers_planted_genes()` --calls--> `bulk_rna_de()`  [EXTRACTED]
  tests/test_bulk_rna.py → aria/scripts/rna_bulk_de.py

## Import Cycles
- None detected.

## Communities (89 total, 8 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.03
Nodes (18): CompletedProcess, _assert_legacy_script_passed(), bool, Path, str, test_hash_file_stable_under_chunk_size(), test_hash_params_order_invariant(), test_legacy_script_wrapper_rejects_printed_failures() (+10 more)

### Community 1 - "Community 1"
Cohesion: 0.06
Nodes (57): Any, bool, float, int, str, _check_env_file_permissions(), _check_synthetic_assets(), console_main() (+49 more)

### Community 2 - "Community 2"
Cohesion: 0.06
Nodes (50): Path, str, str, Any, bool, str, _cb(), GEOConnector (+42 more)

### Community 3 - "Community 3"
Cohesion: 0.06
Nodes (42): Any, bool, int, Path, str, Any, bool, Path (+34 more)

### Community 4 - "Community 4"
Cohesion: 0.08
Nodes (34): bool, EvidenceItem, NarrativeBlock, Path, str, bool, EvidenceItem, NarrativeBlock (+26 more)

### Community 5 - "Community 5"
Cohesion: 0.07
Nodes (20): DesignIntelligence, format_design_intelligence(), OrchestratorAgent, bool, int, str, ARIAMemory, bool (+12 more)

### Community 6 - "Community 6"
Cohesion: 0.07
Nodes (26): Any, ARIAMemory, float, int, LLMProvider, str, float, int (+18 more)

### Community 7 - "Community 7"
Cohesion: 0.06
Nodes (54): float, str, bool, int, float, str, ndarray, _abs_corr() (+46 more)

### Community 8 - "Community 8"
Cohesion: 0.07
Nodes (18): NarrativeAgent, ARIAMemory, LLMProvider, Path, str, Any, bool, str (+10 more)

### Community 9 - "Community 9"
Cohesion: 0.09
Nodes (19): scRNAAgent, ARIAMemory, bool, LLMProvider, Path, str, _agent(), test_pseudobulk_falls_back_to_keywords_without_di_recommendation() (+11 more)

### Community 10 - "Community 10"
Cohesion: 0.07
Nodes (16): ChromatinAgent, ARIAMemory, bool, LLMProvider, str, float, int, str (+8 more)

### Community 11 - "Community 11"
Cohesion: 0.09
Nodes (14): DebateCouncil, DebateResult, DebateRound, DebateVerdict, IntegrationAgent, int, LLMProvider, str (+6 more)

### Community 12 - "Community 12"
Cohesion: 0.09
Nodes (26): str, bool, int, str, hic_inspect(), _inspect_cool(), _inspect_file(), _inspect_hic() (+18 more)

### Community 13 - "Community 13"
Cohesion: 0.08
Nodes (34): Path, str, float, str, bool, bool, Any, bool (+26 more)

### Community 14 - "Community 14"
Cohesion: 0.10
Nodes (18): DesignAgent, AnswerPolicy, ARIAMemory, float, LLMProvider, str, default_answer_policy(), drain_pending_checkpoints() (+10 more)

### Community 15 - "Community 15"
Cohesion: 0.11
Nodes (11): BulkRNAAgent, _is_fastq(), _normalise_sample_token(), ARIAMemory, bool, LLMProvider, str, _agent() (+3 more)

### Community 16 - "Community 16"
Cohesion: 0.17
Nodes (39): _annotation_state(), build_scrna_html_section(), build_scrna_integrated_interpretation(), build_scrna_methods(), _concise_question(), _describe_abundance_de_relationship(), _describe_cellcomm_context(), _describe_pathway_support() (+31 more)

### Community 17 - "Community 17"
Cohesion: 0.14
Nodes (16): DataAuditAgent, _usable_design_col(), _usable_groupby_col(), _usable_replicate_col(), ARIAMemory, bool, int, Path (+8 more)

### Community 18 - "Community 18"
Cohesion: 0.11
Nodes (17): Any, float, int, str, _LazyMessageBus, Message, MessageBus, _escalation() (+9 more)

### Community 19 - "Community 19"
Cohesion: 0.13
Nodes (28): bool, float, str, bool, str, _basic_chromatin_qc(), _bulk_chromatin_qc(), chromatin_qc() (+20 more)

### Community 20 - "Community 20"
Cohesion: 0.17
Nodes (29): RawIngestionAgent, ARIAMemory, str, Path, str, Any, bool, int (+21 more)

### Community 21 - "Community 21"
Cohesion: 0.14
Nodes (21): AuditAgent, _find_star_logs(), _infer_batch_labels(), _infer_group_labels(), _load_count_matrix(), _parse_star_unique_pct(), float, Path (+13 more)

### Community 22 - "Community 22"
Cohesion: 0.20
Nodes (27): Any, bool, float, NarrativeBlock, str, annotate_claim_tiers(), _block_evidence_category(), _block_subject() (+19 more)

### Community 23 - "Community 23"
Cohesion: 0.12
Nodes (23): _infer_lfc_threshold(), float, _infer_groups(), _load_or_infer_metadata(), _run_vst(), _sample_qc(), _select_variable_genes(), make_counts() (+15 more)

### Community 24 - "Community 24"
Cohesion: 0.21
Nodes (27): float, Path, str, adapt(), _entities_from_pb(), _first_mt_threshold(), _genome_for_organism(), _input_record() (+19 more)

### Community 25 - "Community 25"
Cohesion: 0.13
Nodes (20): Any, bool, int, str, Any, float, str, test_annotation_clean_when_all_labels_distinct() (+12 more)

### Community 26 - "Community 26"
Cohesion: 0.16
Nodes (11): ABC, BaseAgent, ARIAMemory, float, int, LLMProvider, Message, str (+3 more)

### Community 27 - "Community 27"
Cohesion: 0.22
Nodes (22): bool, NarrativeBlock, Path, str, _apply_causal_guard(), _apply_low_confidence_warning(), _apply_trajectory_guard(), collect_named_entities() (+14 more)

### Community 28 - "Community 28"
Cohesion: 0.26
Nodes (6): SetupAgent, ARIAMemory, bool, LLMProvider, Path, str

### Community 29 - "Community 29"
Cohesion: 0.18
Nodes (20): bool, float, int, Path, str, main(), int, Path (+12 more)

### Community 30 - "Community 30"
Cohesion: 0.21
Nodes (16): collect_version_metadata(), _fallback_source_hash(), _git_bytes(), _git_text(), Any, Path, str, _repo_root() (+8 more)

### Community 31 - "Community 31"
Cohesion: 0.17
Nodes (18): bool, float, int, _build_design_formula(), _mock_de_result(), _resolve_covariates(), _run_deseq2(), _run_outlier_sensitivity() (+10 more)

### Community 32 - "Community 32"
Cohesion: 0.23
Nodes (18): Any, bool, NarrativeBlock, str, build_evidence_card(), _claim_entities(), _claim_numbers(), EvidenceCard (+10 more)

### Community 33 - "Community 33"
Cohesion: 0.24
Nodes (5): bool, int, str, ContextManager, ModelProfile

### Community 34 - "Community 34"
Cohesion: 0.26
Nodes (6): float, int, Path, str, LLMProvider, test_llm_provider_loads_aria_env_file()

### Community 35 - "Community 35"
Cohesion: 0.16
Nodes (15): Any, str, ModelConfig, build_robustness_multiverse(), _fake_response(), _provider(), test_absent_tier_falls_back_to_medium(), test_heavy_only_config_does_not_keyerror_on_present_tier() (+7 more)

### Community 36 - "Community 36"
Cohesion: 0.23
Nodes (5): GenomeArchAgent, ARIAMemory, int, LLMProvider, str

### Community 37 - "Community 37"
Cohesion: 0.17
Nodes (15): float, bulk_rna_de(), _compute_tpm(), _contrast_overlap(), _format_top_genes(), test_bulk_de_global_contrast_family_is_recorded(), test_contrast_family_can_skip_lfc_gate_when_wald_test_used_lfc_null(), test_contrast_family_significance_pools_and_gates_on_lfc() (+7 more)

### Community 38 - "Community 38"
Cohesion: 0.36
Nodes (19): banner(), check_api_keys(), fail(), find_pbmc_data(), _find_pbmc_dataset(), info(), main(), ok() (+11 more)

### Community 39 - "Community 39"
Cohesion: 0.13
Nodes (10): bool, _json_serializer(), _cache_matches(), _cache_params(), rna_qc(), test_base_json_serializer_handles_numpy(), test_rna_qc_invalid_path_is_structured(), test_rna_qc_on_pbmc3k() (+2 more)

### Community 40 - "Community 40"
Cohesion: 0.23
Nodes (6): DesignStep, CavemanMode, Confidence, MessageType, Enum, TaskTier

### Community 41 - "Community 41"
Cohesion: 0.23
Nodes (17): Any, bool, NarrativeBlock, str, _block_composition_ok(), _block_low_power(), build_devils_advocate(), _challenges_for() (+9 more)

### Community 42 - "Community 42"
Cohesion: 0.19
Nodes (9): EvidenceItem, NarrativeBlock, test_render_blocks_fails_on_unsupported_claim_sentence(), test_render_blocks_shows_claim_evidence_caveats_and_validates_files(), test_render_blocks_stores_claim_verification_metadata(), test_block_round_trips_to_dict_for_methodology_json(), test_success_block_requires_claim_and_evidence(), _de_block() (+1 more)

### Community 43 - "Community 43"
Cohesion: 0.46
Nodes (17): NarrativeBlock, str, compose_block_prose(), _compose_cellcomm(), _compose_composition(), _compose_de(), _compose_gsea(), _compose_non_success() (+9 more)

### Community 44 - "Community 44"
Cohesion: 0.12
Nodes (4): _Decision, MockAdvisor, MockEnvManager, MockMemory

### Community 45 - "Community 45"
Cohesion: 0.26
Nodes (15): str, bool, str, _to_symbols(), _enrichr_enrichment(), _get_gene_sets(), _gseapy_organism(), _mock_pathways() (+7 more)

### Community 46 - "Community 46"
Cohesion: 0.22
Nodes (15): AST, FunctionDef, _call_name(), _module_gates_mocks(), bool, Call, int, str (+7 more)

### Community 47 - "Community 47"
Cohesion: 0.25
Nodes (14): float, str, _avg_pct_passed(), _build_lockfile_section(), _build_raw_ingestion_section(), _build_run_ledger_section(), _build_slug(), _collect_param_hashes() (+6 more)

### Community 48 - "Community 48"
Cohesion: 0.26
Nodes (14): bool, float, int, Path, str, _categorical_colors(), _compute_umap(), _embedding_label() (+6 more)

### Community 49 - "Community 49"
Cohesion: 0.46
Nodes (13): NarrativeBlock, Path, str, _claim_tier_badge(), group_blocks_by_prefix(), _group_key(), _image_uri(), _render_block() (+5 more)

### Community 50 - "Community 50"
Cohesion: 0.27
Nodes (12): find_causal_language(), _causal_block(), _ora_block_with_regulatory_term_name(), test_causal_claim_still_caught_despite_term_name_evidence(), test_find_causal_language_detects_broadened_terms(), test_find_causal_language_exclude_redacts_named_entities(), test_find_causal_language_passes_associative_text(), test_render_suppresses_warning_when_causal_evidence_declared() (+4 more)

### Community 51 - "Community 51"
Cohesion: 0.39
Nodes (11): bool, int, Path, str, _align_sample(), _build_star_index(), _index_exists(), _mock_alignment() (+3 more)

### Community 52 - "Community 52"
Cohesion: 0.27
Nodes (10): float, int, str, export_de_table(), export_pathways_table(), make_gsea_running_sums(), make_ora_dotplot(), _ranked_signature_frame() (+2 more)

### Community 53 - "Community 53"
Cohesion: 0.24
Nodes (6): bool, NarrativeBlock, Path, str, ModalityNarrator, Protocol

### Community 54 - "Community 54"
Cohesion: 0.18
Nodes (10): contrast, denominator, name, numerator, description, design_factor, max_false_up_among_null, min_recall (+2 more)

### Community 55 - "Community 55"
Cohesion: 0.36
Nodes (10): bool, int, Path, str, _detect_samples(), _fastp_outputs_valid(), _mock_fastp_result(), rna_fastq_qc() (+2 more)

### Community 56 - "Community 56"
Cohesion: 0.38
Nodes (10): bool, int, Path, str, _build_ensembl_to_symbol_map(), _clean_counts_matrix(), _counts_outputs_valid(), _detect_strandedness() (+2 more)

### Community 57 - "Community 57"
Cohesion: 0.27
Nodes (10): bool, _cache_params(), _cache_matches(), _cache_params(), rna_concat(), test_clustering_cache_requires_matching_parameters(), test_clustering_skip_leiden_when_cluster_col_provided(), test_concat_cache_requires_matching_manifest() (+2 more)

### Community 58 - "Community 58"
Cohesion: 0.53
Nodes (9): bool, int, str, _balance_cooler(), hic_qc_and_balance(), _mock_hic_qc(), _process_cooler(), _process_hic() (+1 more)

### Community 59 - "Community 59"
Cohesion: 0.38
Nodes (9): str, Exception, object, _get_gene_coordinates(), _get_peak_coordinates(), integration_peak2gene(), _load_atac_matrix(), MissingGTFError (+1 more)

### Community 60 - "Community 60"
Cohesion: 0.29
Nodes (5): NarrativeBlock, ModalityNarrator, NarrativeRegistry, registry_with(), test_registry_collects_first_accepting_narrator()

### Community 61 - "Community 61"
Cohesion: 0.36
Nodes (8): _build_clr_design(), _clr_transform_counts(), rna_diff_abundance(), test_clr_design_adds_donor_fixed_effect_for_paired_design(), test_clr_transform_rows_sum_to_zero(), test_diff_abundance_uses_paired_clr_model(), test_diff_abundance_detects_2x_shift(), test_diff_abundance_no_signal_returns_none_significant()

### Community 62 - "Community 62"
Cohesion: 0.36
Nodes (8): bool, int, Path, str, _categorical_colors(), _draw_paga(), make_paga_figures(), _spring_layout()

### Community 64 - "Community 64"
Cohesion: 0.39
Nodes (7): bool, int, str, mocks_allowed(), integration_mofa(), _load_modality(), _mock_mofa()

### Community 65 - "Community 65"
Cohesion: 0.43
Nodes (7): float, str, _generate_plots(), _plot_heatmap(), _plot_pca_mds(), _plot_sample_pca(), _save_single_dr_plot()

### Community 66 - "Community 66"
Cohesion: 0.57
Nodes (7): Path, str, _gtf_to_symbol_map(), _load_gene_annotation(), _load_symbol_map(), _locate_gtf(), _parse_gtf_biotype_and_length()

### Community 67 - "Community 67"
Cohesion: 0.52
Nodes (6): error(), info(), step(), success(), warn(), install.sh script

### Community 68 - "Community 68"
Cohesion: 0.33
Nodes (6): Any, int, test_provenance_block_contains_required_fields(), collect_llm_usage(), collect_provenance(), record_llm_usage()

### Community 69 - "Community 69"
Cohesion: 0.47
Nodes (5): bool, Path, str, load_aria_env(), _parse_env_line()

### Community 70 - "Community 70"
Cohesion: 0.33
Nodes (6): AuditAgent, DataAuditAgent, DesignAgent, MessageBus, OrchestratorAgent, SetupAgent

### Community 71 - "Community 71"
Cohesion: 0.73
Nodes (5): _load_counts(), test_lognorm_matrix_is_hard_refused_by_default(), test_nonraw_matrix_coerced_only_when_allowed(), test_raw_counts_load_and_are_tagged_raw(), _write_matrix()

### Community 74 - "Community 74"
Cohesion: 0.70
Nodes (4): str, integration_wnn(), _load_atac(), _load_rna()

### Community 77 - "Community 77"
Cohesion: 0.83
Nodes (3): env_installed(), snapshot_one(), generate_locks.sh script

### Community 78 - "Community 78"
Cohesion: 0.67
Nodes (3): filter_graph(), main(), int

### Community 80 - "Community 80"
Cohesion: 0.67
Nodes (3): Claim Compiler, Devil's Advocate, NarrativeAgent

## Knowledge Gaps
- **146 isolated node(s):** `APPLY_fixes.sh script`, `float`, `ARIAMemory`, `LLMProvider`, `CavemanMode` (+141 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **8 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `BaseAgent` connect `Community 26` to `Community 0`, `Community 1`, `Community 35`, `Community 36`, `Community 5`, `Community 6`, `Community 40`, `Community 8`, `Community 10`, `Community 11`, `Community 9`, `Community 14`, `Community 15`, `Community 17`, `Community 20`, `Community 21`, `Community 28`?**
  _High betweenness centrality (0.218) - this node is a cross-community bridge._
- **Why does `rna_pseudobulk_de()` connect `Community 7` to `Community 0`, `Community 1`, `Community 13`, `Community 21`?**
  _High betweenness centrality (0.102) - this node is a cross-community bridge._
- **Why does `scRNAAgent` connect `Community 9` to `Community 0`, `Community 3`, `Community 37`, `Community 6`, `Community 26`?**
  _High betweenness centrality (0.092) - this node is a cross-community bridge._
- **What connects `APPLY_fixes.sh script`, `float`, `ARIAMemory` to the rest of the system?**
  _146 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.030526315789473683 - nodes in this community are weakly interconnected._
- **Should `Community 1` be split into smaller, more focused modules?**
  _Cohesion score 0.0639269406392694 - nodes in this community are weakly interconnected._
- **Should `Community 2` be split into smaller, more focused modules?**
  _Cohesion score 0.05821917808219178 - nodes in this community are weakly interconnected._