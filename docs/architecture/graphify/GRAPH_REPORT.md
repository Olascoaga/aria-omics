# Graph Report - /home/medusa/Samael/ARIA  (2026-06-06)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 2399 nodes · 5777 edges · 105 communities (92 shown, 13 thin omitted)
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `1e9ac0f9`
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
- [[_COMMUNITY_Community 89|Community 89]]
- [[_COMMUNITY_Community 90|Community 90]]
- [[_COMMUNITY_Community 91|Community 91]]
- [[_COMMUNITY_Community 92|Community 92]]
- [[_COMMUNITY_Community 93|Community 93]]
- [[_COMMUNITY_Community 94|Community 94]]
- [[_COMMUNITY_Community 95|Community 95]]
- [[_COMMUNITY_Community 97|Community 97]]
- [[_COMMUNITY_Community 103|Community 103]]
- [[_COMMUNITY_Community 104|Community 104]]

## God Nodes (most connected - your core abstractions)
1. `scRNAAgent` - 45 edges
2. `NarrativeBlock` - 44 edges
3. `BaseAgent` - 40 edges
4. `OrchestratorAgent` - 39 edges
5. `LLMProvider` - 38 edges
6. `BulkRNAAgent` - 37 edges
7. `read_h5ad()` - 37 edges
8. `EvidenceItem` - 36 edges
9. `ParameterAdvisor` - 34 edges
10. `DataAuditAgent` - 33 edges

## Surprising Connections (you probably didn't know these)
- `test_shrinkage_clause_in_narrative()` --calls--> `_lfc_shrinkage_clause()`  [EXTRACTED]
  tests/test_lfc_shrinkage.py → aria/agents/_narrative_scrna.py
- `test_bulk_methods_disclose_a_dropped_covariate()` --calls--> `BulkRnaNarrator`  [EXTRACTED]
  tests/test_bulk_covariates.py → aria/agents/narrative/narrators/bulk_rna.py
- `test_bulk_methods_report_the_fitted_covariate_formula()` --calls--> `BulkRnaNarrator`  [EXTRACTED]
  tests/test_bulk_covariates.py → aria/agents/narrative/narrators/bulk_rna.py
- `test_orchestrator_does_not_dispatch_on_internal_cp3_resolution()` --calls--> `Message`  [EXTRACTED]
  tests/test_pytest_smoke.py → aria/bus/message_bus.py
- `test_geo_fetch_refuses_when_air_gapped()` --calls--> `GEOConnector`  [EXTRACTED]
  tests/test_egress_governance.py → aria/connectors/geo_connector.py

## Import Cycles
- None detected.

## Communities (105 total, 13 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.05
Nodes (90): BiologicalSynthesisAgent, bool, Any, bool, float, NarrativeBlock, str, Any (+82 more)

### Community 1 - "Community 1"
Cohesion: 0.07
Nodes (60): Any, bool, NarrativeBlock, str, Path, str, Any, bool (+52 more)

### Community 2 - "Community 2"
Cohesion: 0.06
Nodes (45): apply_metadata_corrections(), DataAuditAgent, DataAuditScanLimits, default_genome_for_organism(), _env_bool(), _env_float(), _env_int(), _env_optional_int() (+37 more)

### Community 3 - "Community 3"
Cohesion: 0.05
Nodes (45): str, Any, bool, int, Path, str, Any, bool (+37 more)

### Community 4 - "Community 4"
Cohesion: 0.06
Nodes (66): float, str, bool, float, int, str, float, str (+58 more)

### Community 5 - "Community 5"
Cohesion: 0.06
Nodes (24): scRNAAgent, ARIAMemory, bool, LLMProvider, Path, str, float, int (+16 more)

### Community 6 - "Community 6"
Cohesion: 0.03
Nodes (13): bool, CompletedProcess, _cache_params(), _cache_matches(), _assert_legacy_script_passed(), test_clustering_cache_requires_matching_parameters(), test_clustering_skip_leiden_when_cluster_col_provided(), test_concat_cache_requires_matching_manifest() (+5 more)

### Community 7 - "Community 7"
Cohesion: 0.07
Nodes (19): DebateCouncil, DebateResult, DebateRound, DebateVerdict, GenomeArchAgent, IntegrationAgent, int, LLMProvider (+11 more)

### Community 8 - "Community 8"
Cohesion: 0.06
Nodes (26): NarrativeAgent, ARIAMemory, LLMProvider, str, bool, NarrativeBlock, str, ChromatinNarrator (+18 more)

### Community 9 - "Community 9"
Cohesion: 0.07
Nodes (24): OrchestratorAgent, ARIAMemory, bool, float, int, LLMProvider, Message, str (+16 more)

### Community 10 - "Community 10"
Cohesion: 0.08
Nodes (44): AuditAgent, _find_star_logs(), _infer_batch_labels(), _infer_group_labels(), _load_count_matrix(), _parse_star_unique_pct(), _base_card(), _batch_condition_confounded() (+36 more)

### Community 11 - "Community 11"
Cohesion: 0.11
Nodes (42): float, int, Any, bool, float, int, str, ATACDACaller (+34 more)

### Community 12 - "Community 12"
Cohesion: 0.07
Nodes (16): ChromatinAgent, ARIAMemory, bool, LLMProvider, str, float, int, str (+8 more)

### Community 13 - "Community 13"
Cohesion: 0.12
Nodes (46): int, Path, str, ask_biological_question(), _collect_metadata_corrections(), _discard_queued_stdin_lines(), _drain_checkpoints(), _live_analysis_loop() (+38 more)

### Community 14 - "Community 14"
Cohesion: 0.10
Nodes (21): Any, ARIAMemory, float, int, LLMProvider, str, MetricEvaluator, ParameterAdvisor (+13 more)

### Community 15 - "Community 15"
Cohesion: 0.09
Nodes (26): str, bool, int, str, hic_inspect(), _inspect_cool(), _inspect_file(), _inspect_hic() (+18 more)

### Community 16 - "Community 16"
Cohesion: 0.08
Nodes (33): Path, str, float, bool, bool, Any, bool, Path (+25 more)

### Community 17 - "Community 17"
Cohesion: 0.12
Nodes (38): RawIngestionAgent, ARIAMemory, str, Any, int, Path, str, Any (+30 more)

### Community 18 - "Community 18"
Cohesion: 0.10
Nodes (18): DesignAgent, AnswerPolicy, ARIAMemory, float, LLMProvider, str, default_answer_policy(), drain_pending_checkpoints() (+10 more)

### Community 19 - "Community 19"
Cohesion: 0.11
Nodes (12): BulkRNAAgent, _is_fastq(), _normalise_sample_token(), ARIAMemory, bool, LLMProvider, Path, str (+4 more)

### Community 20 - "Community 20"
Cohesion: 0.17
Nodes (39): _annotation_state(), build_scrna_html_section(), build_scrna_integrated_interpretation(), build_scrna_methods(), _concise_question(), _describe_abundance_de_relationship(), _describe_cellcomm_context(), _describe_pathway_support() (+31 more)

### Community 21 - "Community 21"
Cohesion: 0.16
Nodes (29): bool, bytes, float, int, object, Path, str, test_data_audit_classifies_paired_rna_fastq_as_bulk_rna_raw() (+21 more)

### Community 22 - "Community 22"
Cohesion: 0.14
Nodes (20): bool, EvidenceItem, NarrativeBlock, Path, str, Caveat, BulkRnaNarrator, _evidence() (+12 more)

### Community 23 - "Community 23"
Cohesion: 0.13
Nodes (28): bool, float, str, bool, str, _basic_chromatin_qc(), _bulk_chromatin_qc(), chromatin_qc() (+20 more)

### Community 24 - "Community 24"
Cohesion: 0.15
Nodes (17): bool, EvidenceItem, NarrativeBlock, Path, str, _design_issues(), _evidence(), _first_design_issue() (+9 more)

### Community 25 - "Community 25"
Cohesion: 0.14
Nodes (14): Any, float, int, str, Message, MessageBus, _escalation(), _finding() (+6 more)

### Community 26 - "Community 26"
Cohesion: 0.17
Nodes (17): Path, str, _cb(), GEOConnector, _gse_prefix(), _infer_data_type(), _infer_design(), _organism_from_gene_symbols() (+9 more)

### Community 27 - "Community 27"
Cohesion: 0.12
Nodes (23): _infer_lfc_threshold(), float, _infer_groups(), _load_or_infer_metadata(), _run_vst(), _sample_qc(), _select_variable_genes(), make_counts() (+15 more)

### Community 28 - "Community 28"
Cohesion: 0.21
Nodes (27): float, Path, str, adapt(), _entities_from_pb(), _first_mt_threshold(), _genome_for_organism(), _input_record() (+19 more)

### Community 29 - "Community 29"
Cohesion: 0.17
Nodes (22): collect_image_metadata(), collect_version_metadata(), _fallback_source_hash(), _git_bytes(), _git_text(), Any, bytes, Path (+14 more)

### Community 30 - "Community 30"
Cohesion: 0.18
Nodes (7): ABC, DesignStep, CavemanMode, Confidence, MessageType, Enum, TaskTier

### Community 31 - "Community 31"
Cohesion: 0.13
Nodes (20): Any, bool, int, str, Any, float, str, test_annotation_clean_when_all_labels_distinct() (+12 more)

### Community 32 - "Community 32"
Cohesion: 0.22
Nodes (22): bool, NarrativeBlock, Path, str, _apply_causal_guard(), _apply_low_confidence_warning(), _apply_trajectory_guard(), collect_named_entities() (+14 more)

### Community 33 - "Community 33"
Cohesion: 0.19
Nodes (8): DesignIntelligence, format_design_intelligence(), bool, int, str, test_design_intelligence_downgrades_bulk_at_n2(), test_design_intelligence_downgrades_pseudobulk_at_n2(), test_design_intelligence_scrna_focused_group_feasibility()

### Community 34 - "Community 34"
Cohesion: 0.26
Nodes (6): SetupAgent, ARIAMemory, bool, LLMProvider, Path, str

### Community 35 - "Community 35"
Cohesion: 0.19
Nodes (11): EvidenceItem, NarrativeBlock, test_gsea_prose_nes_matches_evidence_after_rounding_fix(), test_render_blocks_does_not_duplicate_associative_badge_label(), test_render_blocks_fails_on_unsupported_claim_sentence(), test_render_blocks_shows_claim_evidence_caveats_and_validates_files(), test_render_blocks_stores_claim_verification_metadata(), test_render_blocks_strict_false_withholds_bad_block_without_aborting() (+3 more)

### Community 36 - "Community 36"
Cohesion: 0.18
Nodes (20): bool, float, int, Path, str, main(), int, Path (+12 more)

### Community 37 - "Community 37"
Cohesion: 0.18
Nodes (19): float, str, _avg_pct_passed(), _build_calibration_badge(), _build_lockfile_section(), _build_raw_ingestion_section(), _build_run_ledger_section(), _build_slug() (+11 more)

### Community 38 - "Community 38"
Cohesion: 0.19
Nodes (10): BaseAgent, ARIAMemory, float, int, LLMProvider, Message, str, CavemanMode (+2 more)

### Community 39 - "Community 39"
Cohesion: 0.15
Nodes (17): Any, str, ModelConfig, build_robustness_multiverse(), _fake_response(), _provider(), test_absent_tier_falls_back_to_medium(), test_heavy_only_config_does_not_keyerror_on_present_tier() (+9 more)

### Community 40 - "Community 40"
Cohesion: 0.17
Nodes (18): bool, float, int, _build_design_formula(), _mock_de_result(), _resolve_covariates(), _run_deseq2(), _run_outlier_sensitivity() (+10 more)

### Community 41 - "Community 41"
Cohesion: 0.23
Nodes (19): str, bool, str, _to_symbols(), _enrichr_enrichment(), _get_gene_sets(), _gseapy_organism(), _mock_pathways() (+11 more)

### Community 42 - "Community 42"
Cohesion: 0.12
Nodes (11): bool, _json_serializer(), _cache_matches(), _cache_params(), rna_qc(), _run_ambient_decontamination(), test_base_json_serializer_handles_numpy(), test_rna_qc_invalid_path_is_structured() (+3 more)

### Community 43 - "Community 43"
Cohesion: 0.24
Nodes (5): bool, int, str, ContextManager, ModelProfile

### Community 44 - "Community 44"
Cohesion: 0.25
Nodes (18): _check_env_file_permissions(), _check_llm(), _check_secrets(), _check_synthetic_assets(), console_main(), main(), _probe_latency(), int (+10 more)

### Community 45 - "Community 45"
Cohesion: 0.18
Nodes (17): Any, bool, str, RuntimeError, test_bulk_local_ora_still_runs_when_air_gapped(), test_egress_gate_tracks_air_gapped(), test_geo_fetch_refuses_when_air_gapped(), _write_human_gmts() (+9 more)

### Community 46 - "Community 46"
Cohesion: 0.39
Nodes (19): bool, NarrativeBlock, str, compose_block_prose(), _compose_cellcomm(), _compose_composition(), _compose_de(), _compose_gsea() (+11 more)

### Community 47 - "Community 47"
Cohesion: 0.36
Nodes (19): banner(), check_api_keys(), fail(), find_pbmc_data(), _find_pbmc_dataset(), info(), main(), ok() (+11 more)

### Community 48 - "Community 48"
Cohesion: 0.28
Nodes (5): float, int, Path, str, LLMProvider

### Community 49 - "Community 49"
Cohesion: 0.12
Nodes (4): _Decision, MockAdvisor, MockEnvManager, MockMemory

### Community 50 - "Community 50"
Cohesion: 0.39
Nodes (15): bool, NarrativeBlock, Path, str, _claim_tier_badge(), group_blocks_by_prefix(), _group_key(), _image_uri() (+7 more)

### Community 51 - "Community 51"
Cohesion: 0.19
Nodes (10): _import_bulk(), _ref_sf(), test_bulk_run_pathway_enrichment_defaults_to_local(), test_bulk_run_pathway_enrichment_skips_without_gmt_or_optin(), test_hypergeom_sf_matches_reference(), test_load_local_library_with_manifest(), test_local_ora_for_databases_reports_missing(), test_parse_gmt() (+2 more)

### Community 52 - "Community 52"
Cohesion: 0.22
Nodes (15): AST, FunctionDef, _call_name(), _module_gates_mocks(), bool, Call, int, str (+7 more)

### Community 53 - "Community 53"
Cohesion: 0.26
Nodes (14): bool, float, int, Path, str, _categorical_colors(), _compute_umap(), _embedding_label() (+6 more)

### Community 54 - "Community 54"
Cohesion: 0.23
Nodes (12): Any, float, int, str, test_clean_when_markers_are_cluster_specific(), test_empty_marker_lists_do_not_count_as_ubiquitous(), test_flags_ubiquitous_markers_across_clusters(), test_graceful_on_empty_markers() (+4 more)

### Community 55 - "Community 55"
Cohesion: 0.25
Nodes (12): Any, bool, str, test_candidate_column_corrected_by_integration_is_downgraded(), test_clean_when_no_candidate_batch_columns(), test_declared_batch_is_not_flagged(), test_flags_unmodeled_technical_column(), test_graceful_on_empty_inputs() (+4 more)

### Community 56 - "Community 56"
Cohesion: 0.19
Nodes (9): Exception, bool, Path, str, diagnose_llm_failure(), test_diagnose_points_to_the_api_key_the_user_actually_has(), test_diagnose_when_no_key_is_present(), load_aria_env() (+1 more)

### Community 57 - "Community 57"
Cohesion: 0.27
Nodes (12): find_causal_language(), _causal_block(), _ora_block_with_regulatory_term_name(), test_causal_claim_still_caught_despite_term_name_evidence(), test_find_causal_language_detects_broadened_terms(), test_find_causal_language_exclude_redacts_named_entities(), test_find_causal_language_passes_associative_text(), test_render_suppresses_warning_when_causal_evidence_declared() (+4 more)

### Community 58 - "Community 58"
Cohesion: 0.27
Nodes (10): float, int, str, export_de_table(), export_pathways_table(), make_gsea_running_sums(), make_ora_dotplot(), _ranked_signature_frame() (+2 more)

### Community 59 - "Community 59"
Cohesion: 0.32
Nodes (10): float, test_bulk_de_global_contrast_family_is_recorded(), test_contrast_family_can_skip_lfc_gate_when_wald_test_used_lfc_null(), test_contrast_family_significance_pools_and_gates_on_lfc(), test_pooled_bh_is_more_conservative_than_per_group(), test_preregister_contrast_family_normalizes_and_declares(), bh_correct(), contrast_family_significance() (+2 more)

### Community 60 - "Community 60"
Cohesion: 0.24
Nodes (6): bool, NarrativeBlock, Path, str, ModalityNarrator, Protocol

### Community 61 - "Community 61"
Cohesion: 0.20
Nodes (12): bool, Path, str, test_raw_ingestion_agent_updates_scrna_modalities(), test_raw_ingestion_converts_10x_triplet_with_reader_provenance(), test_raw_ingestion_detects_valid_10x_triplet(), test_raw_ingestion_fastq_plan_blocks_without_explicit_metadata(), test_raw_ingestion_kb_hash_errors_return_blocker() (+4 more)

### Community 62 - "Community 62"
Cohesion: 0.39
Nodes (11): bool, int, Path, str, _align_sample(), _build_star_index(), _index_exists(), _mock_alignment() (+3 more)

### Community 63 - "Community 63"
Cohesion: 0.18
Nodes (10): contrast, denominator, name, numerator, description, design_factor, max_false_up_among_null, min_recall (+2 more)

### Community 64 - "Community 64"
Cohesion: 0.36
Nodes (10): bool, int, Path, str, _detect_samples(), _fastp_outputs_valid(), _mock_fastp_result(), rna_fastq_qc() (+2 more)

### Community 65 - "Community 65"
Cohesion: 0.38
Nodes (10): bool, int, Path, str, _build_ensembl_to_symbol_map(), _clean_counts_matrix(), _counts_outputs_valid(), _detect_strandedness() (+2 more)

### Community 66 - "Community 66"
Cohesion: 0.29
Nodes (5): NarrativeBlock, ModalityNarrator, NarrativeRegistry, registry_with(), test_registry_collects_first_accepting_narrator()

### Community 67 - "Community 67"
Cohesion: 0.36
Nodes (8): _build_clr_design(), _clr_transform_counts(), rna_diff_abundance(), test_clr_design_adds_donor_fixed_effect_for_paired_design(), test_clr_transform_rows_sum_to_zero(), test_diff_abundance_uses_paired_clr_model(), test_diff_abundance_detects_2x_shift(), test_diff_abundance_no_signal_returns_none_significant()

### Community 68 - "Community 68"
Cohesion: 0.24
Nodes (5): Path, str, test_agent_registry_imports_and_script_contracts_are_valid(), check_registry_integrity(), _scripts()

### Community 69 - "Community 69"
Cohesion: 0.53
Nodes (9): bool, int, str, _balance_cooler(), hic_qc_and_balance(), _mock_hic_qc(), _process_cooler(), _process_hic() (+1 more)

### Community 70 - "Community 70"
Cohesion: 0.36
Nodes (8): bool, int, Path, str, _categorical_colors(), _draw_paga(), make_paga_figures(), _spring_layout()

### Community 71 - "Community 71"
Cohesion: 0.44
Nodes (8): object, str, _get_gene_coordinates(), _get_peak_coordinates(), integration_peak2gene(), _load_atac_matrix(), MissingGTFError, _mock_peak2gene()

### Community 73 - "Community 73"
Cohesion: 0.43
Nodes (7): float, str, _generate_plots(), _plot_heatmap(), _plot_pca_mds(), _plot_sample_pca(), _save_single_dr_plot()

### Community 74 - "Community 74"
Cohesion: 0.57
Nodes (7): Path, str, _gtf_to_symbol_map(), _load_gene_annotation(), _load_symbol_map(), _locate_gtf(), _parse_gtf_biotype_and_length()

### Community 76 - "Community 76"
Cohesion: 0.39
Nodes (7): bool, int, str, mocks_allowed(), integration_mofa(), _load_modality(), _mock_mofa()

### Community 77 - "Community 77"
Cohesion: 0.39
Nodes (5): _proj(), test_console_scripts_present(), test_core_dependencies_have_version_ceilings(), test_no_python_310_classifier(), test_requires_python_is_3_11()

### Community 78 - "Community 78"
Cohesion: 0.52
Nodes (6): error(), info(), step(), success(), warn(), install.sh script

### Community 79 - "Community 79"
Cohesion: 0.43
Nodes (6): str, _get_gene_sets(), _gseapy_organism(), rna_pathway_per_cluster(), test_per_cluster_ora_refuses_enrichr_when_air_gapped(), test_per_cluster_local_ora_success()

### Community 80 - "Community 80"
Cohesion: 0.73
Nodes (5): _load_counts(), test_lognorm_matrix_is_hard_refused_by_default(), test_nonraw_matrix_coerced_only_when_allowed(), test_raw_counts_load_and_are_tagged_raw(), _write_matrix()

### Community 84 - "Community 84"
Cohesion: 0.70
Nodes (4): test_disclosure_is_deterministic(), test_disclosure_keeps_bh_as_primary_and_never_claims_ihw(), test_disclosure_marks_ihw_and_svalues_unavailable_with_reasons(), fdr_advanced_methods_disclosure()

### Community 85 - "Community 85"
Cohesion: 0.80
Nodes (4): env_installed(), snapshot_one(), snapshot_requirements(), generate_locks.sh script

### Community 86 - "Community 86"
Cohesion: 0.70
Nodes (4): str, integration_wnn(), _load_atac(), _load_rna()

### Community 89 - "Community 89"
Cohesion: 0.50
Nodes (4): DataAuditAgent, DesignAgent, MessageBus, OrchestratorAgent

### Community 90 - "Community 90"
Cohesion: 0.67
Nodes (3): filter_graph(), main(), int

### Community 92 - "Community 92"
Cohesion: 1.00
Nodes (3): BiologicalSynthesisAgent, Claim Compiler, NarrativeAgent

## Knowledge Gaps
- **162 isolated node(s):** `APPLY_fixes.sh script`, `float`, `ARIAMemory`, `LLMProvider`, `CavemanMode` (+157 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **13 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `BaseAgent` connect `Community 38` to `Community 2`, `Community 34`, `Community 68`, `Community 5`, `Community 6`, `Community 7`, `Community 8`, `Community 9`, `Community 10`, `Community 12`, `Community 17`, `Community 18`, `Community 19`, `Community 30`?**
  _High betweenness centrality (0.211) - this node is a cross-community bridge._
- **Why does `rna_pseudobulk_de()` connect `Community 4` to `Community 6`, `Community 10`, `Community 11`, `Community 16`, `Community 84`?**
  _High betweenness centrality (0.105) - this node is a cross-community bridge._
- **Why does `OrchestratorAgent` connect `Community 9` to `Community 2`, `Community 38`, `Community 6`, `Community 10`, `Community 12`, `Community 13`, `Community 18`, `Community 30`?**
  _High betweenness centrality (0.080) - this node is a cross-community bridge._
- **What connects `APPLY_fixes.sh script`, `float`, `ARIAMemory` to the rest of the system?**
  _162 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.050314465408805034 - nodes in this community are weakly interconnected._
- **Should `Community 1` be split into smaller, more focused modules?**
  _Cohesion score 0.06893106893106893 - nodes in this community are weakly interconnected._
- **Should `Community 2` be split into smaller, more focused modules?**
  _Cohesion score 0.06220095693779904 - nodes in this community are weakly interconnected._