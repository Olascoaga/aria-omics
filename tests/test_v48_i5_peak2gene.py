"""V48 Slice I5 — Paired Peak-to-Gene (`peak2gene_core` + `integration_peak2gene`) tests.

Covers:
  1. UCSC (`chr1..chr22,chrX,chrY,chrM`) <-> Ensembl (`1..22,X,Y,MT`) contig
     canonicalization, deterministic GTF duplicate `gene_name` resolution
     (`_PAR_Y` drop, `MANE_Select` / `protein_coding` priority, strand-aware TSS),
     and `ContigMismatch` fail-closed error when `<50%` peaks match GTF contigs.
  2. Pure-Python `.fai` byte-offset indexed FASTA per-peak GC fraction calculation.
  3. Rejection of global sequencing-depth confounding by the trans-chromosomal
     matched background null while recovering planted positive AND negative cis
     associations (`seed=0` bit-exact determinism + shared `bh_correct` FDR).
  4. Paired `.h5mu` end-to-end execution (`500kb` vs `250kb` window sensitivity,
     `gc_and_accessibility` background matching, `uns["aria_peak2gene"]` HDF5
     persistence, and `ScriptContract` beta output validation).
  5. Separate `.h5ad` deterministic barcode alignment under scrambled cell order,
     zero-link header-only CSV contract, structured `MissingGTF` / `MissingRNA`
     errors, `IntegrationAgent` associative tunnel storage (positive + negative),
     and AST no-fabrication guard across all Slice I5 modules.
"""

from __future__ import annotations

import csv
import inspect
import math
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _has_module(mod_name: str) -> bool:
    try:
        __import__(mod_name)
        return True
    except ImportError:
        return False


def _require_module(mod_name: str):
    try:
        return __import__(mod_name)
    except ImportError:
        try:
            import pytest

            pytest.skip(f"{mod_name} not installed in active environment")
        except ImportError as exc:
            raise RuntimeError(f"Required module {mod_name} missing") from exc


def _write_indexed_fasta(path: Path, records: list[tuple[str, str]], line_bases: int = 20) -> Path:
    """Write a multi-line FASTA and exact `.fai` index for unit testing."""
    fai_lines: list[str] = []
    offset = 0
    with open(path, "wb") as fh:
        for name, seq in records:
            header = f">{name}\n".encode("ascii")
            fh.write(header)
            offset += len(header)
            seq_start = offset
            for i in range(0, len(seq), line_bases):
                chunk = (seq[i : i + line_bases] + "\n").encode("ascii")
                fh.write(chunk)
                offset += len(chunk)
            fai_lines.append(
                f"{name}\t{len(seq)}\t{seq_start}\t{line_bases}\t{line_bases + 1}\n"
            )
    Path(f"{path}.fai").write_text("".join(fai_lines), encoding="utf-8")
    return path


def _build_synthetic_multiome_arrays(seed: int = 42):
    """Construct a 3-chromosome paired count dataset with depth confounding + cis truth.

    - `GENE_POS` (`chr1`, TSS=100,000, `+`): linked to promoter peak `chr1:99800-100300`
      and distal peak `chr1:150000-150500` (+50kb) + ultra-distal `chr1:450000-450500` (+350kb).
    - `GENE_NEG` (`chr2`, TSS=200,000, `-`): negatively correlated with `chr2:180000-180500`
      (20kb downstream of `-` TSS -> `tss_offset_bp > 0`).
    - `GENE_DEPTH` (`chr3`, TSS=300,000, `+`): housekeeping gene + cis peak `chr3:310000-310500`
      driven ONLY by global cell library depth (raw `r > 0.4`, must be rejected by
      trans-chromosomal background null!).
    """
    rng = np.random.default_rng(seed)
    n_cells = 140
    # Cell-specific global sequencing depth multiplier (confounder across all chromosomes).
    depth = rng.uniform(0.5, 3.5, size=n_cells)
    # Orthogonal biological cell-state programs.
    prog_a = np.repeat([0.0, 1.0], n_cells // 2) + rng.normal(0.0, 0.12, size=n_cells)
    prog_a = np.clip(prog_a, 0.05, 1.05)
    prog_b = np.tile([0.0, 1.0], n_cells // 2) + rng.normal(0.0, 0.12, size=n_cells)
    prog_b = np.clip(prog_b, 0.05, 1.05)

    rna_genes = ["GENE_POS", "GENE_NEG", "GENE_DEPTH", "HK1", "HK2", "HK3", "HK4", "HK5"]
    rna_counts = np.zeros((n_cells, len(rna_genes)), dtype=np.int32)
    # Scale biological programs by depth so CP10k normalization isolates prog_a / prog_b,
    # whereas raw ATAC counts still carry `depth` unless controlled by background null!
    rna_counts[:, 0] = rng.poisson(depth * (2.0 + 28.0 * prog_a))
    rna_counts[:, 1] = rng.poisson(depth * (2.0 + 28.0 * (1.15 - prog_b)))
    # Make GENE_DEPTH increase super-linearly with depth so even after CP10k it correlates
    # strongly with depth across cells (creating a severe raw cis correlation trap!).
    rna_counts[:, 2] = rng.poisson(depth * (2.0 + 12.0 * depth))
    for hk_i in range(3, len(rna_genes)):
        rna_counts[:, hk_i] = rng.poisson(depth * 25.0)

    peak_names: list[str] = []
    peak_cols: list[np.ndarray] = []

    # Planted cis peaks:
    peak_names.append("chr1:99800-100300")  # promoter of GENE_POS (dist=0)
    peak_cols.append(rng.poisson(depth * (0.4 + 4.5 * prog_a)))

    peak_names.append("chr1:150000-150500")  # +50kb distal cis of GENE_POS
    peak_cols.append(rng.poisson(depth * (0.4 + 4.2 * prog_a)))

    peak_names.append("chr1:450000-450500")  # +350kb distal cis of GENE_POS (in 500kb, out of 250kb)
    peak_cols.append(rng.poisson(depth * (0.4 + 4.0 * prog_a)))

    peak_names.append("chr2:180000-180500")  # negative cis peak near GENE_NEG (dist ~19.5kb)
    peak_cols.append(rng.poisson(depth * (0.4 + 4.5 * prog_b)))

    peak_names.append("chr3:310000-310500")  # depth-only cis peak near GENE_DEPTH (NO biological program!)
    peak_cols.append(rng.poisson(depth * 2.2))

    # 120 trans background peaks across chr1, chr2, chr3 (outside cis windows, >= 1.5 Mb away)
    # driven by the same depth process so the background null accurately captures depth correlation!
    for idx in range(40):
        pos = 1_500_000 + idx * 10_000
        for chrom in ("chr1", "chr2", "chr3"):
            peak_names.append(f"{chrom}:{pos}-{pos + 500}")
            rate = rng.uniform(1.2, 3.2)
            peak_cols.append(rng.poisson(depth * rate))

    atac_counts = np.column_stack(peak_cols).astype(np.int32)

    gtf_lines = [
        'chr1\ttest\tgene\t100000\t120000\t.\t+\t.\tgene_id "ENSG000001"; gene_name "GENE_POS"; gene_type "protein_coding"; tag "MANE_Select";\n',
        'chr2\ttest\tgene\t190000\t200000\t.\t-\t.\tgene_id "ENSG000002"; gene_name "GENE_NEG"; gene_type "protein_coding";\n',
        'chr3\ttest\tgene\t300000\t315000\t.\t+\t.\tgene_id "ENSG000003"; gene_name "GENE_DEPTH"; gene_type "protein_coding";\n',
    ]
    barcodes = [f"CELL_{i:04d}-1" for i in range(n_cells)]
    return barcodes, rna_genes, rna_counts, peak_names, atac_counts, "".join(gtf_lines)


def test_contig_canonicalization_and_gtf_duplicate_resolution(tmp_path):
    from aria.utils.peak2gene_core import (
        canonical_contig,
        compute_peak_to_gene_links,
        genes_to_peak_indices,
        parse_gtf_tss,
        parse_gtf_tss_detailed,
    )

    assert canonical_contig("1") == "chr1"
    assert canonical_contig("chr1") == "chr1"
    assert canonical_contig("X") == "chrX"
    assert canonical_contig("MT") == "chrM"
    assert canonical_contig("chrM") == "chrM"
    assert canonical_contig("GL000195.1") == "GL000195.1"

    gtf = tmp_path / "ensembl_dups.gtf"
    gtf.write_text(
        # PAR_Y copy should be dropped
        'Y\tensembl\tgene\t10000\t20000\t.\t+\t.\tgene_id "ENSG000100_PAR_Y"; gene_name "PLCXD1"; gene_biotype "protein_coding";\n'
        'X\tensembl\tgene\t10000\t20000\t.\t+\t.\tgene_id "ENSG000100"; gene_name "PLCXD1"; gene_biotype "protein_coding";\n'
        # Duplicate symbol: pseudogene vs MANE protein_coding on minus strand
        '1\tensembl\tgene\t500000\t501000\t.\t+\t.\tgene_id "ENSG000999"; gene_name "DUPGENE"; gene_biotype "processed_pseudogene";\n'
        '1\tensembl\tgene\t800000\t830000\t.\t-\t.\tgene_id "ENSG000200"; gene_name "DUPGENE"; gene_biotype "protein_coding"; tag "Ensembl_canonical";\n',
        encoding="utf-8",
    )

    detailed, stats = parse_gtf_tss_detailed(str(gtf))
    assert stats["n_par_y_dropped"] == 1
    assert stats["n_duplicate_symbols_resolved"] == 1
    assert detailed["PLCXD1"]["canon_chrom"] == "chrX"
    assert detailed["DUPGENE"]["gene_id"] == "ENSG000200"
    assert detailed["DUPGENE"]["tss"] == 830000  # minus strand end
    assert detailed["DUPGENE"]["canon_chrom"] == "chr1"

    # Ensembl GTF ("1", "X") seamlessly matches UCSC ATAC peaks ("chr1:...", "chrX:...")
    tss_map = parse_gtf_tss(str(gtf))
    hits = genes_to_peak_indices(
        ["chr1:825000-826000", "chrX:12000-12500", "chr2:100-200"],
        tss_map,
        window_bp=10_000,
    )
    assert hits == {"DUPGENE": [0], "PLCXD1": [1]}

    # ContigMismatch hard error when >=20 peaks belong to incompatible contigs
    bad_peaks = [f"scaffold_{i}:100-600" for i in range(25)]
    mismatch = compute_peak_to_gene_links(
        rna_matrix=np.ones((10, 2)),
        rna_var_names=["DUPGENE", "PLCXD1"],
        atac_matrix=np.ones((10, 25)),
        atac_var_names=bad_peaks,
        gene_records=detailed,
        min_cells_detected=1,
    )
    assert mismatch["status"] == "error"
    assert mismatch["error_type"] == "ContigMismatch"


def test_indexed_fasta_gc_computation(tmp_path):
    from aria.utils.peak2gene_core import compute_peak_gc_from_indexed_fasta

    seq_chr1 = "A" * 40 + "G" * 20 + "C" * 20 + "T" * 40  # positions 41..80 are 100% GC
    seq_chr2 = "ACGT" * 30  # exact 50% GC everywhere
    fa = _write_indexed_fasta(
        tmp_path / "mini.fa",
        [("1", seq_chr1), ("chr2", seq_chr2)],
        line_bases=17,
    )

    coords = [
        ("chr1", 1, 40),    # 0% GC via Ensembl "1" <-> UCSC "chr1" canonical lookup
        ("chr1", 41, 80),   # 100% GC across line wraps
        ("chr2", 11, 90),   # 50% GC
    ]
    gc, info = compute_peak_gc_from_indexed_fasta(str(fa), coords)
    assert info["available"] is True
    assert info["n_computed"] == 3
    assert math.isclose(float(gc[0]), 0.0, abs_tol=1e-6)
    assert math.isclose(float(gc[1]), 1.0, abs_tol=1e-6)
    assert math.isclose(float(gc[2]), 0.5, abs_tol=1e-6)


def test_background_null_rejects_depth_confounder_and_is_seed_deterministic(tmp_path):
    from aria.utils.peak2gene_core import (
        compute_peak_to_gene_links,
        parse_gtf_tss_detailed,
    )

    _, rna_genes, rna_counts, peak_names, atac_counts, gtf_txt = (
        _build_synthetic_multiome_arrays(seed=42)
    )
    gtf_p = tmp_path / "test.gtf"
    gtf_p.write_text(gtf_txt, encoding="utf-8")
    gene_records, _ = parse_gtf_tss_detailed(str(gtf_p))

    # 1. Without background null, the depth-confounded pair (GENE_DEPTH, chr3:310000-310500)
    # exhibits a strong spurious positive correlation (r > 0.35).
    naive = compute_peak_to_gene_links(
        rna_matrix=rna_counts,
        rna_var_names=rna_genes,
        atac_matrix=atac_counts,
        atac_var_names=peak_names,
        gene_records=gene_records,
        distance_bp=500_000,
        min_corr=0.30,
        rna_is_raw_counts=True,
        atac_is_raw_counts=True,
        preprocess="paired",
        background_null=False,
        seed=0,
    )
    naive_genes = {r["gene"] for r in naive["links"]}
    assert "GENE_DEPTH" in naive_genes

    # 2. With trans-chromosomal matched background null + BH FDR (max_padj=0.05),
    # GENE_DEPTH is REJECTED while true cis links (GENE_POS positive, GENE_NEG negative) survive!
    res0 = compute_peak_to_gene_links(
        rna_matrix=rna_counts,
        rna_var_names=rna_genes,
        atac_matrix=atac_counts,
        atac_var_names=peak_names,
        gene_records=gene_records,
        distance_bp=500_000,
        min_corr=0.30,
        max_padj=0.05,
        rna_is_raw_counts=True,
        atac_is_raw_counts=True,
        preprocess="paired",
        background_null=True,
        seed=0,
    )
    res0_repeat = compute_peak_to_gene_links(
        rna_matrix=rna_counts,
        rna_var_names=rna_genes,
        atac_matrix=atac_counts,
        atac_var_names=peak_names,
        gene_records=gene_records,
        distance_bp=500_000,
        min_corr=0.30,
        max_padj=0.05,
        rna_is_raw_counts=True,
        atac_is_raw_counts=True,
        preprocess="paired",
        background_null=True,
        seed=0,
    )

    assert res0["status"] == "success"
    assert res0["links"] == res0_repeat["links"]

    by_pair = {(r["gene"], r["peak"]): r for r in res0["links"]}
    assert ("GENE_POS", "chr1:99800-100300") in by_pair
    assert ("GENE_POS", "chr1:150000-150500") in by_pair
    assert ("GENE_NEG", "chr2:180000-180500") in by_pair
    assert all(r["gene"] != "GENE_DEPTH" for r in res0["links"])

    pos_prom = by_pair[("GENE_POS", "chr1:99800-100300")]
    assert pos_prom["direction"] == "positive"
    assert pos_prom["region_class"] == "promoter_proximal"
    assert pos_prom["distance_bp"] == 0
    assert pos_prom["background_scope"] == "trans"
    assert pos_prom["padj"] < 0.01

    neg_distal = by_pair[("GENE_NEG", "chr2:180000-180500")]
    assert neg_distal["direction"] == "negative"
    assert neg_distal["region_class"] == "distal_cis"
    # Minus-strand TSS=200,000, peak midpoint=180,250 -> downstream (+19,750 bp)
    assert neg_distal["tss_offset_bp"] == 19750
    assert neg_distal["padj"] < 0.01


def test_integration_peak2gene_paired_h5mu_windows_and_uns_persistence(tmp_path):
    ad = _require_module("anndata")
    md = _require_module("mudata")
    import scipy.sparse as sp

    from aria.scripts.integration_peak2gene import integration_peak2gene

    barcodes, rna_genes, rna_counts, peak_names, atac_counts, gtf_txt = (
        _build_synthetic_multiome_arrays(seed=7)
    )
    gtf_p = tmp_path / "annot.gtf"
    gtf_p.write_text(gtf_txt, encoding="utf-8")

    # Build synthetic indexed FASTA covering chr1..chr3 up to 2 Mb so GC matching activates.
    rng_seq = np.random.default_rng(11)
    bases = np.array(list("ACGT"))
    fa_records = [
        (chrom, "".join(bases[rng_seq.integers(0, 4, size=2_000_000)]))
        for chrom in ("chr1", "chr2", "chr3")
    ]
    fa_p = _write_indexed_fasta(tmp_path / "genome.fa", fa_records, line_bases=80)

    rna_ad = ad.AnnData(sp.csr_matrix(rna_counts))
    rna_ad.obs_names = barcodes
    rna_ad.var_names = rna_genes
    atac_ad = ad.AnnData(sp.csr_matrix(atac_counts))
    atac_ad.obs_names = barcodes
    atac_ad.var_names = peak_names

    mdata = md.MuData({"rna": rna_ad, "atac": atac_ad})
    h5mu_p = tmp_path / "paired.h5mu"
    mdata.write(h5mu_p)

    res_500 = integration_peak2gene(
        {
            "h5mu_path": str(h5mu_p),
            "gtf_path": str(gtf_p),
            "genome_fasta": str(fa_p),
            "distance_kb": 500,
            "min_corr": 0.30,
            "max_padj": 0.05,
            "write_h5mu": True,
            "output_dir": str(tmp_path / "out_500"),
            "seed": 0,
        }
    )
    assert res_500["status"] == "success"
    assert res_500["validation_level"] == "production"
    assert res_500["background_match"] == "gc_and_accessibility"
    assert res_500["gc_matching"]["available"] is True
    assert res_500["ctcf_validated"] is False
    assert res_500["hic_corroborated"] is False
    assert res_500["n_positive_correlations"] >= 2
    assert res_500["n_negative_correlations"] >= 1
    assert any("NEGATIVE correlation" in w for w in res_500["warnings"])

    if _has_module("pydantic"):
        from aria.utils.script_contracts import contract_for_script

        contract = contract_for_script("aria/scripts/integration_peak2gene.py")
        assert contract is not None
        assert contract.validation_level == "production"
        assert contract.validate_result(res_500) == []

    with open(res_500["output_path"], newline="", encoding="utf-8") as fh:
        rows_500 = list(csv.DictReader(fh))
    peaks_500 = {(r["gene"], r["peak"]) for r in rows_500}
    assert ("GENE_POS", "chr1:450000-450500") in peaks_500
    assert ("GENE_POS", "chr1:150000-150500") in peaks_500

    # Verify persisted .h5mu metadata
    assert res_500["h5mu_output_path"] is not None
    m_reloaded = md.read(res_500["h5mu_output_path"])
    assert "aria_peak2gene" in m_reloaded.uns
    assert int(m_reloaded.uns["aria_peak2gene"]["n_links"]) == res_500["n_links"]

    # 250kb window must drop the +350kb link while keeping the +50kb link!
    res_250 = integration_peak2gene(
        {
            "h5mu_path": str(h5mu_p),
            "gtf_path": str(gtf_p),
            "genome_fasta": str(fa_p),
            "distance_kb": 250,
            "min_corr": 0.30,
            "max_padj": 0.05,
            "output_dir": str(tmp_path / "out_250"),
            "seed": 0,
        }
    )
    assert res_250["status"] == "success"
    with open(res_250["output_path"], newline="", encoding="utf-8") as fh:
        rows_250 = list(csv.DictReader(fh))
    peaks_250 = {(r["gene"], r["peak"]) for r in rows_250}
    assert ("GENE_POS", "chr1:450000-450500") not in peaks_250
    assert ("GENE_POS", "chr1:150000-150500") in peaks_250


def test_scrambled_h5ad_barcode_alignment_zero_links_and_structured_errors(tmp_path):
    ad = _require_module("anndata")
    import scipy.sparse as sp

    from aria.scripts.integration_peak2gene import integration_peak2gene

    barcodes, rna_genes, rna_counts, peak_names, atac_counts, gtf_txt = (
        _build_synthetic_multiome_arrays(seed=19)
    )
    gtf_p = tmp_path / "annot.gtf"
    gtf_p.write_text(gtf_txt, encoding="utf-8")

    rna_ad = ad.AnnData(sp.csr_matrix(rna_counts))
    rna_ad.obs_names = barcodes
    rna_ad.var_names = rna_genes
    rna_p = tmp_path / "rna.h5ad"
    rna_ad.write_h5ad(rna_p)

    # Reverse ATAC cell order completely to test canonical barcode alignment!
    rev_idx = np.arange(len(barcodes))[::-1]
    atac_ad = ad.AnnData(sp.csr_matrix(atac_counts[rev_idx, :]))
    atac_ad.obs_names = [barcodes[i] for i in rev_idx]
    atac_ad.var_names = peak_names
    atac_p = tmp_path / "atac_scrambled.h5ad"
    atac_ad.write_h5ad(atac_p)

    res = integration_peak2gene(
        {
            "rna_files": [str(rna_p)],
            "atac_files": [str(atac_p)],
            "gtf_path": str(gtf_p),
            "genome": "synthetic_no_fasta",
            "distance_kb": 500,
            "min_corr": 0.30,
            "max_padj": 0.05,
            "output_dir": str(tmp_path / "out_aligned"),
            "seed": 0,
        }
    )
    assert res["status"] == "success"
    assert res["background_match"] == "accessibility_and_length"
    linked = {(r["gene"], r["peak"]): r["direction"] for r in res["top_links"]}
    assert linked.get(("GENE_POS", "chr1:150000-150500")) == "positive"
    assert linked.get(("GENE_NEG", "chr2:180000-180500")) == "negative"

    # Zero-link run still writes a valid header CSV and satisfies ScriptContract!
    res_zero = integration_peak2gene(
        {
            "rna_files": [str(rna_p)],
            "atac_files": [str(atac_p)],
            "gtf_path": str(gtf_p),
            "min_corr": 0.9999,
            "output_dir": str(tmp_path / "out_zero"),
        }
    )
    assert res_zero["status"] == "success"
    assert res_zero["n_links"] == 0
    assert Path(res_zero["output_path"]).is_file()
    if _has_module("pydantic"):
        from aria.utils.script_contracts import contract_for_script

        contract = contract_for_script("aria/scripts/integration_peak2gene.py")
        assert contract.validate_result(res_zero) == []

    # Structured error checks
    err_rna = integration_peak2gene({"rna_files": ["/no/such/rna.h5ad"], "atac_files": []})
    assert err_rna["status"] == "error" and err_rna["error_type"] == "MissingRNA"

    err_gtf = integration_peak2gene(
        {
            "rna_files": [str(rna_p)],
            "atac_files": [str(atac_p)],
            "gtf_path": str(tmp_path / "missing.gtf"),
        }
    )
    assert err_gtf["status"] == "error" and err_gtf["error_type"] == "MissingGTF"


def test_agent_associative_tunnels_orchestrator_gate_and_no_fabrication_guard(tmp_path):
    import ast
    import importlib.util

    repo_root = Path(__file__).resolve().parents[1]
    target_rels = (
        "aria/utils/peak2gene_core.py",
        "aria/scripts/integration_peak2gene.py",
        "aria/scripts/chromatin_regulatory.py",
    )

    for rel in target_rels:
        src = (repo_root / rel).read_text(encoding="utf-8")
        tree = ast.parse(src, filename=rel)
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                assert not node.name.startswith("_mock_"), f"Unexpected mock fn in {rel}: {node.name}"

    if _has_module("pytest"):
        spec = importlib.util.spec_from_file_location(
            "test_no_fabrication_guard",
            repo_root / "tests" / "test_no_fabrication_guard.py",
        )
        assert spec is not None and spec.loader is not None
        fab_mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fab_mod)
        for rel in target_rels:
            src = (repo_root / rel).read_text(encoding="utf-8")
            assert fab_mod.scan_source(rel, src) == [], f"Fabrication guard flagged {rel}"

    if not _has_module("pydantic"):
        return

    from aria.agents.integration_agent import IntegrationAgent
    from aria.agents.orchestrator_agent import INTEGRATION_VALIDATION
    from aria.utils.script_contracts import contract_for_script

    # Promoted to paired same-cell production (ADR-079)!
    assert INTEGRATION_VALIDATION["level"] == "production"
    assert INTEGRATION_VALIDATION["dispatch_enabled"] is True
    assert INTEGRATION_VALIDATION["requires_ack"] is False
    assert INTEGRATION_VALIDATION["production_promoted"] is True

    # ScriptContract is promoted to production and validates canonical success payload
    contract = contract_for_script("aria/scripts/integration_peak2gene.py")
    assert contract is not None
    assert contract.validation_level == "production"
    dummy_csv = tmp_path / "peak_gene_links.csv"
    dummy_csv.write_text("gene,peak,correlation\n", encoding="utf-8")
    assert (
        contract.validate_result(
            {
                "status": "success",
                "output_path": str(dummy_csv),
                "n_links": 0,
                "n_positive_correlations": 0,
                "n_negative_correlations": 0,
                "top_links": [],
                "ctcf_validated": False,
                "hic_corroborated": False,
                "warnings": [],
            }
        )
        == []
    )

    # IntegrationAgent._store_peak_gene_tunnels stores both strong positive and negative links
    tunnels: list[dict] = []

    class _Mem:
        def create_tunnel(self, **kw):
            tunnels.append(kw)

    agent = IntegrationAgent.__new__(IntegrationAgent)
    agent.memory = _Mem()
    agent._store_peak_gene_tunnels(
        "exp_i5",
        {
            "top_links": [
                {"gene": "G1", "peak": "chr1:10-20", "correlation": 0.65, "distance_kb": 15.0},
                {"gene": "G2", "peak": "chr2:30-40", "correlation": -0.48, "distance_kb": 42.5},
                {"gene": "G3", "peak": "chr3:50-60", "correlation": 0.31, "distance_kb": 10.0},
            ]
        },
    )
    assert len(tunnels) == 2
    assert {t["entity"] for t in tunnels} == {"G1", "G2"}
    assert all("associative" in t["description"] for t in tunnels)


if __name__ == "__main__":
    test_fns = [
        obj
        for name, obj in list(globals().items())
        if name.startswith("test_") and callable(obj)
    ]
    n_pass = 0
    for fn in test_fns:
        sig = inspect.signature(fn)
        if "tmp_path" in sig.parameters:
            with tempfile.TemporaryDirectory(prefix="test_v48_i5_") as td:
                fn(Path(td))
        else:
            fn()
        n_pass += 1
        print(f"  PASS {fn.__name__}")
    print(f"\nAll {n_pass}/{len(test_fns)} V48 Slice I5 tests passed.")
