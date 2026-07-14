"""scRNA narrative helpers subpackage (A7 extraction).

The former monolithic ``aria/agents/_narrative_scrna.py`` (2.3k lines mixing
formatting, text, tables and figures) is split into cohesive modules here:

- ``_common``  — findings normalisation, formatting and selectors (leaf utils);
- ``text``     — descriptions, summary, integrated interpretation, methods;
- ``tables``   — HTML table bodies and supplementary TSV export;
- ``figures``  — matplotlib renders, HTML section, figure orchestration.

The full public surface is re-exported here (and again from the compatibility
facade ``aria/agents/_narrative_scrna.py``) so existing consumers
(``narrative_agent``, ``report_builder``, ``devils_advocate``, the scRNA
narrator) keep importing the same names. Behavior is unchanged and pinned by
``tests/test_narrative_scrna_contract.py``.
"""
from __future__ import annotations

from aria.agents.narrative.scrna._common import (
    _annotation_state,
    _fdr_primary_clause,
    _fmt_int,
    _fmt_stat,
    _find_pathway_block,
    _gene_brief,
    _gene_name,
    _group_label,
    _label_cell_type,
    _lfc_shrinkage_clause,
    _term_value,
    _top_de_blocks,
    _top_directional_genes,
    _top_pathway_blocks,
    _top_pathway_terms,
    unwrap_scrna_findings,
)
from aria.agents.narrative.scrna.text import (
    _concise_question,
    _describe_abundance_de_relationship,
    _describe_cellcomm_context,
    _describe_pathway_support,
    _describe_top_de_blocks,
    _describe_trajectory_context,
    build_scrna_integrated_interpretation,
    build_scrna_methods,
    summarize_scrna_text,
)
from aria.agents.narrative.scrna.tables import (
    _write_tsv,
    export_supplementary_tables,
    extract_cellcomm_table,
    extract_pseudobulk_de_table,
    extract_trajectory_tables,
)
from aria.agents.narrative.scrna.figures import (
    _embed_png,
    build_scrna_html_section,
    generate_figures,
    render_cellcomm_heatmap,
    render_cellcomm_top_pairs_bar,
    render_pathway_dotplots,
    render_per_celltype_de_bar,
)

__all__ = [
    # _common
    "unwrap_scrna_findings", "_fmt_int", "_fdr_primary_clause",
    "_lfc_shrinkage_clause", "_fmt_stat", "_group_label", "_label_cell_type",
    "_annotation_state", "_top_de_blocks", "_top_pathway_blocks", "_gene_name",
    "_gene_brief", "_top_directional_genes", "_find_pathway_block",
    "_top_pathway_terms",
    # text
    "_describe_top_de_blocks", "_describe_pathway_support",
    "_describe_abundance_de_relationship", "_describe_cellcomm_context",
    "_describe_trajectory_context", "summarize_scrna_text",
    "build_scrna_integrated_interpretation", "_concise_question",
    "build_scrna_methods",
    # tables
    "extract_pseudobulk_de_table", "extract_cellcomm_table",
    "extract_trajectory_tables", "_write_tsv", "_term_value",
    "export_supplementary_tables",
    # figures
    "render_pathway_dotplots", "render_cellcomm_heatmap",
    "render_cellcomm_top_pairs_bar", "render_per_celltype_de_bar", "_embed_png",
    "build_scrna_html_section", "generate_figures",
]
