#!/usr/bin/env python3
"""Deterministic Public Dataset Provisioner for the ARIA 11-Modality Showcase Atlas.

Stages 100% public, peer-reviewed omics datasets (NCBI GEO, EBI ENA, ENCODE
Portal, Bioconductor ExperimentHub, and NCBI RefSeq) into an ext4 workspace
(default: ``~/.aria/public_atlas_datasets``), recording SHA-256 checksums,
public accessions, DOIs, source URLs, and experimental metadata in per-modality
``dataset_manifest.json`` receipts while keeping analyzable payload files in an
isolated ``inputs/`` subdirectory.

All 11 production analytical workflows in ARIA v4.8+ are covered:
  01_bulk_rna_matrix         — SEQC/MAQC-III (GEO GSE47792, BGI UHRR vs HBRR, 4v4)
  02_bulk_rna_raw_fastq      — Human Airway Dex (GEO GSE52778 / ENA SRP033351, 4v4)
  03_scrna_multidonor        — Human Lupus PBMC IFN-beta (GEO GSE96583 / EH2259, 8v8 donors)
  04_scatac_multidonor       — Human Hippocampus Lifespan scATAC (GEO GSE278576, 3v3 donors)
  05_bulk_atac_tobias        — ENCODE GM12878 vs K562 Bulk ATAC + TOBIAS (3v3)
  06_chip_seq                — ENCODE GM12878 vs K562 H3K27ac ChIP-seq + Input (2v2 + 2 Input)
  07_cut_and_run_spikein     — U2OS SCAF1-AID Auxin Degron CUT&RUN + E. coli spike-in (GEO GSE273126, 3v3 + ctrl)
  08_cut_and_tag             — Kaya-Okur et al. 2019 H1 vs K562 H3K27me3 CUT&Tag (GEO GSE124557, 2v2 + 2 IgG)
  09_hic_3d_genome           — Human RPE-1 Diploid vs Tetraploid Hi-C .mcool (GEO GSE345433, 3v3)
  10_paired_multiome         — Human Hippocampus Same-Cell 10x Multiome (GEO GSE278576, 3v3 donors)
  11_cross_modal_concordance — ENCODE Matched GM12878 vs K562 3v3 Bulk RNA x 3v3 Bulk ATAC (N=6)
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import zlib

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ATLAS_ROOT = Path.home() / ".aria" / "public_atlas_datasets"

# Six gene-dense, biologically canonical GRCh38 chromosomes (~752 Mb, >10,000 genes)
# spanning chr1 (general/housekeeping), chr6 (MHC/immune), chr7 (HOXA/Polycomb),
# chr11 (beta-globin/SPI1/erythroid), chr17 (HOXB/TP53/BRCA1), and chr19 (lymphoid/ZNF).
ATLAS_CHROMS = ("chr1", "chr6", "chr7", "chr11", "chr17", "chr19")
BULK_CHROMS = ("chr6", "chr7", "chr11", "chr17", "chr19")

ALL_MODALITY_IDS = (
    "01_bulk_rna_matrix",
    "02_bulk_rna_raw_fastq",
    "03_scrna_multidonor",
    "04_scatac_multidonor",
    "05_bulk_atac_tobias",
    "06_chip_seq",
    "07_cut_and_run_spikein",
    "08_cut_and_tag",
    "09_hic_3d_genome",
    "10_paired_multiome",
    "11_cross_modal_concordance",
)

# Public NCBI GEO GSM mapping for GSE278576 human hippocampal Multiome donors
GSE278576_DONOR_GSMS: dict[str, dict[str, str]] = {
    "hc29": {"rna_gsm": "GSM8549398", "atac_gsm": "GSM8549450"},
    "hc77": {"rna_gsm": "GSM8549419", "atac_gsm": "GSM8549471"},
    "hc76": {"rna_gsm": "GSM8549418", "atac_gsm": "GSM8549470"},
    "hc11": {"rna_gsm": "GSM8549386", "atac_gsm": "GSM8549438"},
    "hc12": {"rna_gsm": "GSM8549392", "atac_gsm": "GSM8549444"},
    "hc19": {"rna_gsm": "GSM8549395", "atac_gsm": "GSM8549447"},
}


def _log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] [fetch_atlas] {msg}", flush=True)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _conda_run(env_name: str, cmd: list[str], *, cwd: Path | None = None) -> None:
    child_env = dict(os.environ)
    child_env["PYTHONPATH"] = (
        str(ROOT) + os.pathsep + child_env.get("PYTHONPATH", "")
    ).rstrip(os.pathsep)
    subprocess.run(
        ["conda", "run", "--no-capture-output", "-n", env_name, *cmd],
        check=True,
        cwd=str(cwd or ROOT),
        env=child_env,
    )


def _download_file(url: str, dest: Path, *, retries: int = 4, timeout: int = 180) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            _log(f"Downloading {url} -> {dest.name} (attempt {attempt}/{retries})")
            req = urllib.request.Request(
                url, headers={"User-Agent": "ARIA-Omics-Public-Atlas/4.8"}
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp, tmp.open("wb") as out:
                shutil.copyfileobj(resp, out, length=1024 * 1024)
            if tmp.stat().st_size == 0:
                raise RuntimeError(f"Empty download from {url}")
            tmp.replace(dest)
            return dest
        except Exception as exc:
            last_err = exc
            if tmp.exists():
                tmp.unlink(missing_ok=True)
            time.sleep(min(12, 3 * attempt))
    raise RuntimeError(f"Failed downloading {url}: {last_err}") from last_err


def _stream_fastq_head_gz(
    url: str, dest: Path, *, n_reads: int, retries: int = 4, timeout: int = 180
) -> Path:
    """Stream the exact first ``n_reads`` FASTQ records from an ENA .fastq.gz URL.

    ENA read ordering (SRRxxxx.1 .. SRRxxxx.N) is immutable, making prefix
    extraction byte-reproducible across machines without downloading 10+ GB SRA
    archives.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file() and dest.stat().st_size > 1024:
        return dest
    target_lines = int(n_reads) * 4
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            _log(
                f"Streaming first {n_reads:,} reads from {url} -> {dest.name} "
                f"(attempt {attempt}/{retries})"
            )
            req = urllib.request.Request(
                url, headers={"User-Agent": "ARIA-Omics-Public-Atlas/4.8"}
            )
            lines_written = 0
            with urllib.request.urlopen(req, timeout=timeout) as resp, gzip.GzipFile(
                fileobj=resp, mode="rb"
            ) as gz_in, gzip.open(
                tmp, "wb", compresslevel=6
            ) as gz_out:
                while lines_written < target_lines:
                    line = gz_in.readline()
                    if not line:
                        break
                    if not line.endswith(b"\n"):
                        line = line + b"\n"
                    gz_out.write(line)
                    lines_written += 1
            if lines_written < 4000 or lines_written % 4 != 0:
                raise RuntimeError(
                    f"Incomplete FASTQ stream ({lines_written} lines) from {url}"
                )
            tmp.replace(dest)
            return dest
        except Exception as exc:
            last_err = exc
            tmp.unlink(missing_ok=True)
            time.sleep(min(12, 3 * attempt))
    raise RuntimeError(f"Failed streaming FASTQ {url}: {last_err}") from last_err


def _inputs_dir(mod_dir: Path) -> Path:
    d = mod_dir / "inputs"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_ready_manifest(mod_dir: Path, payload: dict) -> Path:
    in_dir = _inputs_dir(mod_dir)
    files_meta = []
    for fname in sorted(payload.get("payload_files", [])):
        p = in_dir / fname
        files_meta.append(
            {
                "file": f"inputs/{fname}",
                "filename": fname,
                "size_bytes": p.stat().st_size,
                "sha256": _sha256(p),
            }
        )
    out = {
        **payload,
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "verified_files": files_meta,
    }
    manifest_path = mod_dir / "dataset_manifest.json"
    manifest_path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def _is_ready(mod_dir: Path) -> bool:
    mp = mod_dir / "dataset_manifest.json"
    if not mp.is_file():
        return False
    try:
        meta = json.loads(mp.read_text(encoding="utf-8"))
        for row in meta.get("verified_files", []):
            fp = mod_dir / row["file"]
            if not fp.is_file() or fp.stat().st_size != int(row["size_bytes"]):
                return False
        return bool(meta.get("verified_files"))
    except Exception:
        return False


# ── Shared Reference Builder (UCSC chr-normalized GRCh38 + E. coli MG1655) ────


def _normalize_to_ucsc_chr(chrom: str) -> str | None:
    c = chrom.strip()
    if c.startswith("chr"):
        base = c[3:]
    else:
        base = c
    if base in {str(i) for i in range(1, 23)} | {"X", "Y"}:
        return f"chr{base}"
    if base in {"M", "MT"}:
        return "chrM"
    return None


def ensure_references(atlas_root: Path) -> dict[str, Path]:
    """Build UCSC ``chr``-prefixed GRCh38 6-chromosome FASTA/GTF/bwa-mem2 index,
    full nuclear UCSC GTF (for genome-wide Hi-C), and E. coli K-12 MG1655 spike-in."""
    ref_root = atlas_root / "_references"
    hg38_sub = ref_root / "hg38_6chrom_ucsc"
    hg38_full_ucsc = ref_root / "hg38_full_ucsc"
    ecoli_dir = ref_root / "ecoli_mg1655"
    for d in (hg38_sub, hg38_full_ucsc, ecoli_dir):
        d.mkdir(parents=True, exist_ok=True)

    src_hg38_fa = Path.home() / ".aria" / "genomes" / "hg38" / "genome.fa"
    src_hg38_gtf = Path.home() / ".aria" / "genomes" / "hg38" / "annotation.gtf"
    if not src_hg38_fa.is_file() or not src_hg38_gtf.is_file():
        raise FileNotFoundError(
            "Expected local GRCh38 reference in ~/.aria/genomes/hg38/{genome.fa,annotation.gtf}. "
            "Run ARIA SetupAgent first to stage GRCh38."
        )

    fa_out = hg38_sub / "genome.fa"
    fai_out = hg38_sub / "genome.fa.fai"
    gtf_out = hg38_sub / "annotation.gtf"
    bwa_marker = hg38_sub / "genome.fa.bwt.2bit.64"
    full_gtf_out = hg38_full_ucsc / "annotation.gtf"

    keep_ucsc = set(ATLAS_CHROMS)

    if not fa_out.is_file() or not fai_out.is_file():
        _log(f"Building UCSC chr-normalized FASTA for {ATLAS_CHROMS} ...")
        tmp_fa = fa_out.with_suffix(".tmp.fa")
        writing = False
        seen: set[str] = set()
        with src_hg38_fa.open("r", encoding="utf-8") as fin, tmp_fa.open(
            "w", encoding="utf-8"
        ) as fout:
            for line in fin:
                if line.startswith(">"):
                    raw_name = line[1:].strip().split()[0]
                    ucsc = _normalize_to_ucsc_chr(raw_name)
                    if ucsc in keep_ucsc and ucsc not in seen:
                        writing = True
                        seen.add(ucsc)
                        fout.write(f">{ucsc}\n")
                    else:
                        writing = False
                elif writing:
                    fout.write(line)
        if seen != keep_ucsc:
            tmp_fa.unlink(missing_ok=True)
            raise RuntimeError(
                f"Expected {sorted(keep_ucsc)} in {src_hg38_fa}, found {sorted(seen)}"
            )
        tmp_fa.replace(fa_out)
        _conda_run("aria-atacseq-env", ["samtools", "faidx", str(fa_out)])

    if not gtf_out.is_file() or not full_gtf_out.is_file():
        _log("Building UCSC chr-normalized GRCh38 GTF annotations ...")
        tmp_sub = gtf_out.with_suffix(".tmp.gtf")
        tmp_full = full_gtf_out.with_suffix(".tmp.gtf")
        with src_hg38_gtf.open("r", encoding="utf-8") as fin, tmp_sub.open(
            "w", encoding="utf-8"
        ) as f_sub, tmp_full.open("w", encoding="utf-8") as f_full:
            for line in fin:
                if line.startswith("#"):
                    f_sub.write(line)
                    f_full.write(line)
                    continue
                chrom, rest = line.split("\t", 1)
                ucsc = _normalize_to_ucsc_chr(chrom)
                if not ucsc:
                    continue
                norm_line = f"{ucsc}\t{rest}"
                f_full.write(norm_line)
                if ucsc in keep_ucsc:
                    f_sub.write(norm_line)
        tmp_sub.replace(gtf_out)
        tmp_full.replace(full_gtf_out)

    if not bwa_marker.is_file():
        _log("Building bwa-mem2 index for hg38_6chrom_ucsc (~45s) ...")
        _conda_run("aria-atacseq-env", ["bwa-mem2", "index", str(fa_out)])

    ecoli_fa = ecoli_dir / "ecoli_mg1655.fa"
    ecoli_bwa = ecoli_dir / "ecoli_mg1655.fa.bwt.2bit.64"
    if not ecoli_fa.is_file():
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            "?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text"
        )
        _download_file(url, ecoli_fa)
    if not ecoli_bwa.is_file():
        _log("Building bwa-mem2 index for E. coli K-12 MG1655 spike-in ...")
        _conda_run(
            "aria-atacseq-env",
            [
                "bash",
                "-c",
                f"samtools faidx '{ecoli_fa}' && bwa-mem2 index '{ecoli_fa}'",
            ],
        )

    return {
        "hg38_chr_fa": fa_out,
        "hg38_chr_fai": fai_out,
        "hg38_chr_gtf": gtf_out,
        "hg38_full_ucsc_gtf": full_gtf_out,
        "hg38_ensembl_fa": src_hg38_fa,
        "hg38_ensembl_gtf": src_hg38_gtf,
        "ecoli_fa": ecoli_fa,
    }


def _slice_encode_bam(
    accession: str,
    dest_bam: Path,
    *,
    local_candidates: list[Path] | None = None,
    chroms: tuple[str, ...] = BULK_CHROMS,
    raw_cache_dir: Path,
) -> Path:
    """Produce a coordinate-sorted, indexed GRCh38 BAM slice preserving @PG headers."""
    bai_dest = Path(str(dest_bam) + ".bai")
    if dest_bam.is_file() and bai_dest.is_file():
        return dest_bam
    dest_bam.parent.mkdir(parents=True, exist_ok=True)
    raw_cache_dir.mkdir(parents=True, exist_ok=True)

    source_bam: Path | None = None
    for cand in local_candidates or []:
        if cand.is_file() and cand.stat().st_size > 1_000_000:
            source_bam = cand
            break

    if source_bam is None:
        dl_bam = raw_cache_dir / f"{accession}.bam"
        url = f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.bam"
        _download_file(url, dl_bam, timeout=360)
        source_bam = dl_bam

    src_bai = Path(str(source_bam) + ".bai")
    if not src_bai.is_file():
        _log(f"Indexing {source_bam.name} ...")
        _conda_run("aria-atacseq-env", ["samtools", "index", "-@", "8", str(source_bam)])

    tmp_bam = dest_bam.with_suffix(".tmp.bam")
    chrom_str = " ".join(chroms)
    _log(f"Slicing {accession} ({chrom_str}) -> {dest_bam.name}")
    _conda_run(
        "aria-atacseq-env",
        [
            "bash",
            "-c",
            f"samtools view -@ 8 -b '{source_bam}' {chrom_str} > '{tmp_bam}' "
            f"&& samtools index -@ 8 '{tmp_bam}'",
        ],
    )
    tmp_bam.replace(dest_bam)
    Path(str(tmp_bam) + ".bai").replace(bai_dest)
    return dest_bam


# ── Modality 01: bulk_RNA (Count Matrix — SEQC GSE47792) ──────────────────────


def stage_01_bulk_rna_matrix(atlas_root: Path) -> Path:
    mod_dir = atlas_root / "01_bulk_rna_matrix"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)

    seqc_src = Path.home() / ".aria" / "benchmarks" / "seqc_maqc_BGI" / "counts.tsv"
    if not seqc_src.is_file():
        _log("Downloading SEQC/MAQC-III GSE47792 BGI count matrix from NCBI GEO ...")
        gz_path = atlas_root / "_raw_cache" / "GSE47792_SEQC_ILM_BGI.txt.gz"
        _download_file(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47792/suppl/"
            "GSE47792_SEQC_ILM_BGI_GE_N_20120619.txt.gz",
            gz_path,
            timeout=240,
        )
        seqc_src = gz_path

    keep_cols = ["A_1", "A_2", "A_3", "A_4", "B_1", "B_2", "B_3", "B_4"]
    rename_map = {
        "A_1": "UHRR_rep1",
        "A_2": "UHRR_rep2",
        "A_3": "UHRR_rep3",
        "A_4": "UHRR_rep4",
        "B_1": "HBRR_rep1",
        "B_2": "HBRR_rep2",
        "B_3": "HBRR_rep3",
        "B_4": "HBRR_rep4",
    }
    out_counts = in_dir / "seqc_uhrr_vs_hbrr_counts.tsv"
    opener = gzip.open if seqc_src.suffix == ".gz" else open
    agg_counts: dict[str, list[int]] = {}
    with opener(seqc_src, "rt", encoding="utf-8") as fin:
        header = fin.readline().rstrip("\n").split("\t")
        col_idx = [header.index(c) for c in keep_cols]
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            gene = parts[0].strip()
            if not gene or gene.startswith("ERCC-"):
                continue
            vals = [int(round(float(parts[i]))) for i in col_idx]
            cur = agg_counts.get(gene)
            if cur is None:
                agg_counts[gene] = vals
            else:
                for k in range(len(vals)):
                    cur[k] += vals[k]

    with out_counts.open("w", encoding="utf-8") as fout:
        fout.write("gene\t" + "\t".join(rename_map[c] for c in keep_cols) + "\n")
        for gene in sorted(agg_counts):
            if sum(agg_counts[gene]) == 0:
                continue
            fout.write(gene + "\t" + "\t".join(str(v) for v in agg_counts[gene]) + "\n")

    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "01_bulk_rna_matrix",
            "aria_modality": "bulk_RNA",
            "title": "SEQC/MAQC-III Human Reference Transcriptomes (UHRR vs HBRR, BGI Center)",
            "accession": "GEO GSE47792",
            "doi": "10.1038/nbt.2957",
            "citation": "SEQC/MAQC-III Consortium, Nature Biotechnology 32, 903-914 (2014)",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "4 Universal Human Reference RNA (UHRR) vs 4 Human Brain Reference RNA (HBRR) replicates (ERCC synthetic spike-ins excluded)",
            "source_urls": [
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47792"
            ],
            "payload_files": ["seqc_uhrr_vs_hbrr_counts.tsv"],
        },
    )
    return mod_dir


# ── Modality 02: bulk_RNA_raw (Paired-End FASTQ — Airway GSE52778) ────────────


def stage_02_bulk_rna_raw_fastq(atlas_root: Path) -> Path:
    mod_dir = atlas_root / "02_bulk_rna_raw_fastq"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)

    runs = [
        ("untreated_rep1", "SRR1039508", "008"),
        ("untreated_rep2", "SRR1039512", "002"),
        ("untreated_rep3", "SRR1039516", "006"),
        ("untreated_rep4", "SRR1039520", "000"),
        ("dexamethasone_rep1", "SRR1039509", "009"),
        ("dexamethasone_rep2", "SRR1039513", "003"),
        ("dexamethasone_rep3", "SRR1039517", "007"),
        ("dexamethasone_rep4", "SRR1039521", "001"),
    ]
    payload_files = []
    urls = []
    for sid, srr, sub in runs:
        for mate in (1, 2):
            url = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/{sub}/{srr}/{srr}_{mate}.fastq.gz"
            fname = f"{sid}_R{mate}.fastq.gz"
            _stream_fastq_head_gz(url, in_dir / fname, n_reads=250_000)
            payload_files.append(fname)
            urls.append(url)

    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "02_bulk_rna_raw_fastq",
            "aria_modality": "bulk_RNA_raw",
            "title": "Human Airway Smooth Muscle Glucocorticoid (Dexamethasone) Response",
            "accession": "GEO GSE52778 / ENA SRP033351",
            "doi": "10.1371/journal.pone.0099625",
            "citation": "Himes BE et al., PLOS ONE 9(6): e99625 (2014)",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "4 untreated vs 4 dexamethasone-treated primary human airway smooth muscle donor lines (first 250,000 paired-end reads per ENA run)",
            "source_urls": urls,
            "payload_files": payload_files,
        },
    )
    return mod_dir


# ── Modality 03: scRNA (Multi-Donor PBMC IFN-beta — Kang 2018 GSE96583) ───────


def stage_03_scrna_multidonor(atlas_root: Path) -> Path:
    mod_dir = atlas_root / "03_scrna_multidonor"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)

    rda_path = Path.home() / ".aria" / "benchmarks" / "kang18_8vs8.rda"
    if not rda_path.is_file():
        rda_path = atlas_root / "_raw_cache" / "kang18_8vs8.rda"
        _download_file(
            "https://experimenthub.bioconductor.org/fetch/2292",
            rda_path,
            timeout=240,
        )

    out_h5ad = in_dir / "kang2018_pbmc_8donors_ifnb.h5ad"
    if not out_h5ad.is_file():
        tmp_export = mod_dir / "_tmp_r_export"
        tmp_export.mkdir(parents=True, exist_ok=True)
        r_code = f"""
        suppressPackageStartupMessages({{
            library(SingleCellExperiment)
            library(Matrix)
        }})
        env <- new.env()
        load("{rda_path}", envir = env)
        obj <- env[[ls(env)[1]]]
        cd <- as.data.frame(colData(obj))
        if ("multiplets" %in% colnames(cd)) {{
            keep_s <- cd$multiplets == "singlet" & !is.na(cd$cell)
            obj <- obj[, keep_s]
            cd <- as.data.frame(colData(obj))
        }}
        set.seed(42)
        cd$grp <- paste0(cd$ind, "_", cd$stim)
        idx_list <- split(seq_len(nrow(cd)), cd$grp)
        sel <- unlist(lapply(idx_list, function(ix) sample(ix, min(length(ix), 400L))), use.names = FALSE)
        sel <- sort(sel)
        sub_obj <- obj[, sel]
        cts <- counts(sub_obj)
        Matrix::writeMM(cts, file = "{tmp_export}/matrix.mtx")
        writeLines(rownames(cts), con = "{tmp_export}/genes.txt")
        sub_cd <- as.data.frame(colData(sub_obj))
        sub_cd$barcode <- colnames(cts)
        write.table(sub_cd, file = "{tmp_export}/obs.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)
        """
        _conda_run("aria-bench-env", ["Rscript", "-e", r_code])

        py_code = f"""
import anndata as ad
import pandas as pd
from scipy import io
from pathlib import Path

exp = Path("{tmp_export}")
X = io.mmread(str(exp / "matrix.mtx")).T.tocsr().astype("int32")
genes = (exp / "genes.txt").read_text().splitlines()
obs = pd.read_csv(exp / "obs.tsv", sep="\\t")
obs["donor_id"] = "donor_" + obs["ind"].astype(str)
obs["condition"] = obs["stim"].map({{"ctrl": "control", "stim": "ifnb_stimulated"}}).fillna(obs["stim"].astype(str))
obs["sample_id"] = obs["donor_id"] + "_" + obs["condition"]
obs["author_cell_type"] = obs["cell"].astype(str)
obs.index = [f"{{r.sample_id}}_cell{{i}}" for i, r in enumerate(obs.itertuples())]
var = pd.DataFrame(index=pd.Index(genes, name="gene_symbol"))
adata = ad.AnnData(X=X, obs=obs[["sample_id", "donor_id", "condition", "author_cell_type"]], var=var)
adata.var_names_make_unique()
adata.write_h5ad("{out_h5ad}")
"""
        _conda_run("aria-rna-env", ["python", "-c", py_code])
        shutil.rmtree(tmp_export, ignore_errors=True)

    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "03_scrna_multidonor",
            "aria_modality": "scRNA",
            "title": "Human Lupus PBMC IFN-beta Stimulation Across 8 Donors",
            "accession": "GEO GSE96583 / Bioconductor EH2259",
            "doi": "10.1038/nbt.4042",
            "citation": "Kang HM et al., Nature Biotechnology 36, 89-94 (2018)",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "8 matched human donors x 2 conditions (control vs ifnb_stimulated, 16 donor-condition libraries, 6,400 singlet PBMCs)",
            "source_urls": [
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583",
                "https://experimenthub.bioconductor.org/fetch/2292",
            ],
            "payload_files": ["kang2018_pbmc_8donors_ifnb.h5ad"],
        },
    )
    return mod_dir


# ── Helper for Public GEO GSE278576 Hippocampus Multiome / scATAC ─────────────


def _fetch_gsm_suppl_urls(gsm_id: str) -> list[str]:
    prefix = gsm_id[:7] + "nnn"
    soft_url = (
        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        f"?acc={gsm_id}&targ=self&form=text&view=brief"
    )
    req = urllib.request.Request(
        soft_url, headers={"User-Agent": "ARIA-Omics-Public-Atlas/4.8"}
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        txt = resp.read().decode("utf-8", "replace")
    urls = []
    for line in txt.splitlines():
        if line.startswith("!Sample_supplementary_file"):
            u = line.split("=", 1)[1].strip()
            if u.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
                u = "https://" + u[len("ftp://") :]
            if u and u.lower() != "none":
                urls.append(u)
    return urls


def _ensure_gse278576_donor_h5mu(donor_id: str, raw_cache: Path) -> Path:
    local_oligos = (
        Path.home()
        / "Samael"
        / "Oligos"
        / "data"
        / "interim"
        / "paired"
        / f"{donor_id}_paired.h5mu"
    )
    if local_oligos.is_file():
        return local_oligos

    cached = raw_cache / "gse278576" / f"{donor_id}_paired.h5mu"
    if cached.is_file():
        return cached

    gsm_info = GSE278576_DONOR_GSMS.get(donor_id)
    if not gsm_info:
        raise KeyError(f"Unknown GSE278576 donor_id: {donor_id}")

    dl_dir = raw_cache / "gse278576" / donor_id
    dl_dir.mkdir(parents=True, exist_ok=True)
    for gsm in (gsm_info["rna_gsm"], gsm_info["atac_gsm"]):
        for u in _fetch_gsm_suppl_urls(gsm):
            fname = u.rstrip("/").split("/")[-1]
            _download_file(u, dl_dir / fname, timeout=360)

    h5_files = sorted(dl_dir.glob("*.h5"))
    if not h5_files:
        raise RuntimeError(f"No .h5 files found in GEO download for {donor_id}")

    assemble_params = {
        "mode": "arc_h5",
        "arc_h5_path": str(h5_files[0]),
        "sample_id": donor_id,
        "output_path": str(cached),
    }
    p_json = dl_dir / "assemble_params.json"
    p_json.write_text(json.dumps(assemble_params), encoding="utf-8")
    _conda_run(
        "aria-integration-env",
        ["python", "aria/scripts/integration_assemble_multiome.py", str(p_json)],
    )
    return cached


# ── Modality 04: scATAC (Multi-Donor Human Hippocampus — GEO GSE278576) ───────


def stage_04_scatac_multidonor(atlas_root: Path) -> Path:
    mod_dir = atlas_root / "04_scatac_multidonor"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)
    raw_cache = atlas_root / "_raw_cache"

    donors = [
        ("hc29", "young_20_39", "M", 28),
        ("hc77", "young_20_39", "M", 20),
        ("hc76", "young_20_39", "F", 26),
        ("hc11", "aged_80_100", "M", 86),
        ("hc12", "aged_80_100", "M", 83),
        ("hc19", "aged_80_100", "F", 87),
    ]
    ref_h5ad = (
        Path.home()
        / "Samael"
        / "Oligos"
        / "data"
        / "raw"
        / "reference"
        / "GSE278576_hippocampus_RNA.h5ad"
    )

    src_map = {d[0]: str(_ensure_gse278576_donor_h5mu(d[0], raw_cache)) for d in donors}
    py_script = f"""
import json
from pathlib import Path
import numpy as np
import mudata as md
import anndata as ad

donors = {json.dumps(donors)}
src_map = {json.dumps(src_map)}
keep_chroms = set({json.dumps(list(BULK_CHROMS))})
in_dir = Path("{in_dir}")
ref_path = Path("{ref_h5ad}")

label_lookup = {{}}
if ref_path.is_file():
    ref_ad = ad.read_h5ad(str(ref_path), backed="r")
    obs_df = ref_ad.obs[["orig.ident", "subclass"]].copy()
    for idx_val, row in obs_df.iterrows():
        bc = str(idx_val).split("_")[-1].split("#")[-1]
        label_lookup[(str(row["orig.ident"]), bc)] = str(row["subclass"])

rng = np.random.default_rng(42)
for donor_id, age_group, sex, age in donors:
    out_p = in_dir / f"{{donor_id}}_scatac.h5mu"
    if out_p.is_file():
        continue
    mdata = md.read_h5mu(src_map[donor_id])
    atac = mdata.mod["atac"].copy()
    peak_chroms = [str(p).split(":")[0].split("-")[0] for p in atac.var_names]
    peak_mask = np.array([c in keep_chroms for c in peak_chroms], dtype=bool)
    atac = atac[:, peak_mask].copy()
    totals = np.asarray(atac.X.sum(axis=1)).ravel()
    valid = np.where(totals >= 50)[0]
    chosen = np.sort(rng.choice(valid, size=min(500, len(valid)), replace=False))
    atac = atac[chosen, :].copy()
    subclasses = []
    for bc in atac.obs_names:
        raw_bc = str(bc).split("#")[-1]
        subclasses.append(label_lookup.get((donor_id, raw_bc), "Glia_Neuron"))
    atac.obs["sample_id"] = donor_id
    atac.obs["donor_id"] = donor_id
    atac.obs["age_group"] = age_group
    atac.obs["condition"] = age_group
    atac.obs["replicate"] = donor_id
    atac.obs["sex"] = sex
    atac.obs["age"] = int(age)
    atac.obs["subclass"] = subclasses
    atac.obs_names = [f"{{donor_id}}#{{str(b).split('#')[-1]}}" for b in atac.obs_names]
    out_md = md.MuData({{"atac": atac}})
    for col in ("sample_id", "donor_id", "age_group", "condition", "replicate", "sex", "subclass"):
        out_md.obs[col] = atac.obs[col].values
    out_md.write_h5mu(str(out_p))
"""
    _conda_run("aria-chromatin-env", ["python", "-c", py_script])

    payload_files = [f"{d[0]}_scatac.h5mu" for d in donors]
    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "04_scatac_multidonor",
            "aria_modality": "scATAC",
            "title": "Human Hippocampus Lifespan Single-Cell Chromatin Accessibility",
            "accession": "GEO GSE278576",
            "doi": "10.1126/science.adt8307",
            "citation": "Zemke NR et al., Science (2026) doi:10.1126/science.adt8307",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "6 post-mortem human hippocampal donors (3 young_20_39 [hc29, hc77, hc76] vs 3 aged_80_100 [hc11, hc12, hc19], 500 nuclei/donor = 3,000 single-nucleus ATAC profiles across chr6/7/11/17/19)",
            "source_urls": [
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278576"
            ],
            "payload_files": payload_files,
        },
    )
    return mod_dir


# ── Modality 05: bulk_ATAC + TOBIAS (ENCODE GM12878 vs K562 3v3) ──────────────


def stage_05_bulk_atac_tobias(atlas_root: Path) -> Path:
    mod_dir = atlas_root / "05_bulk_atac_tobias"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)
    raw_cache = atlas_root / "_raw_cache" / "encode_atac"

    local_d = Path("/mnt/d/Samael/datasets/atac_bulk/bams")
    samples = [
        ("GM12878_rep1_atac.bam", "ENCFF981FXV", [local_d / "GM12878_rep1.ENCFF981FXV.bam"]),
        ("GM12878_rep2_atac.bam", "ENCFF962FMH", [local_d / "GM12878_rep2.ENCFF962FMH.bam"]),
        ("GM12878_rep3_atac.bam", "ENCFF440GRZ", []),
        ("K562_rep1_atac.bam", "ENCFF512VEZ", [local_d / "K562_rep1.ENCFF512VEZ.bam"]),
        ("K562_rep2_atac.bam", "ENCFF987XOV", [local_d / "K562_rep2.ENCFF987XOV.bam"]),
        ("K562_rep3_atac.bam", "ENCFF557RSA", []),
    ]
    payload_files = []
    urls = []
    for fname, encff, cands in samples:
        _slice_encode_bam(
            encff,
            in_dir / fname,
            local_candidates=cands,
            chroms=BULK_CHROMS,
            raw_cache_dir=raw_cache,
        )
        payload_files.append(fname)
        urls.append(f"https://www.encodeproject.org/files/{encff}/@@download/{encff}.bam")

    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "05_bulk_atac_tobias",
            "aria_modality": "bulk_ATAC",
            "title": "ENCODE Human Hematopoietic Chromatin Accessibility & TF Footprinting (GM12878 vs K562)",
            "accession": "ENCODE ENCSR637XSC, ENCSR483RKN, ENCSR859USB",
            "doi": "10.1038/s41586-020-2493-4",
            "citation": "ENCODE Project Consortium, Nature 583, 699-710 (2020)",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "3 GM12878 B-lymphoblastoid replicates vs 3 K562 erythroleukemia replicates (GRCh38 chr6, chr7, chr11, chr17, chr19) enabling inferential pyDESeq2 DA + TOBIAS differential footprinting",
            "source_urls": urls,
            "payload_files": payload_files,
        },
    )
    return mod_dir


# ── Modality 06: ChIP-seq (ENCODE Broad H3K27ac GM12878 vs K562 + Input) ──────


def stage_06_chip_seq(atlas_root: Path) -> Path:
    mod_dir = atlas_root / "06_chip_seq"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)
    raw_cache = atlas_root / "_raw_cache" / "encode_chip"

    samples = [
        ("GM12878_H3K27ac_rep1_chip.bam", "ENCFF804NCH"),
        ("GM12878_H3K27ac_rep2_chip.bam", "ENCFF948GTC"),
        ("K562_H3K27ac_rep1_chip.bam", "ENCFF301TVL"),
        ("K562_H3K27ac_rep2_chip.bam", "ENCFF879BWC"),
        ("GM12878_Input_ctrl_chip.bam", "ENCFF873ZWP"),
        ("K562_Input_ctrl_chip.bam", "ENCFF829HZY"),
    ]
    payload_files = []
    urls = []
    for fname, encff in samples:
        _slice_encode_bam(
            encff,
            in_dir / fname,
            chroms=BULK_CHROMS,
            raw_cache_dir=raw_cache,
        )
        payload_files.append(fname)
        urls.append(f"https://www.encodeproject.org/files/{encff}/@@download/{encff}.bam")

    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "06_chip_seq",
            "aria_modality": "ChIP",
            "title": "ENCODE Broad Histone H3K27ac ChIP-seq with Matched Input Controls (GM12878 vs K562)",
            "accession": "ENCODE ENCSR000AKC, ENCSR000AKP, ENCSR000AKJ, ENCSR000AKY",
            "doi": "10.1038/nature11247",
            "citation": "ENCODE Project Consortium / Bernstein Lab, Nature 489, 57-74 (2012)",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "2 GM12878 H3K27ac IP + Input vs 2 K562 H3K27ac IP + Input biological replicates (GRCh38 chr6, chr7, chr11, chr17, chr19)",
            "source_urls": urls,
            "payload_files": payload_files,
        },
    )
    return mod_dir


# ── Modality 07: CUT&RUN + E. coli Spike-In (U2OS SCAF1-AID Auxin GSE273126) ──


def stage_07_cut_and_run_spikein(atlas_root: Path) -> Path:
    mod_dir = atlas_root / "07_cut_and_run_spikein"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)

    runs = [
        ("Control_rep1_cutnrun", "SRR29997464", "064"),
        ("Control_rep2_cutnrun", "SRR29997462", "062"),
        ("Control_rep3_cutnrun", "SRR29997460", "060"),
        ("Auxin_rep1_cutnrun", "SRR29997463", "063"),
        ("Auxin_rep2_cutnrun", "SRR29997461", "061"),
        ("Auxin_rep3_cutnrun", "SRR29997459", "059"),
        ("IgG_Untagged_ctrl_cutnrun", "SRR29997458", "058"),
    ]
    fq_dump = shutil.which("fastq-dump") or "/home/medusa/anaconda3/envs/rnaseq_env/bin/fastq-dump"
    payload_files = []
    urls = []
    for sid, srr, sub in runs:
        r1_url = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/{sub}/{srr}/{srr}_1.fastq.gz"
        r2_url = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/{sub}/{srr}/{srr}_2.fastq.gz"
        r1_dest = in_dir / f"{sid}_R1.fastq.gz"
        r2_dest = in_dir / f"{sid}_R2.fastq.gz"
        if not (r1_dest.is_file() and r2_dest.is_file() and r1_dest.stat().st_size > 1024 and r2_dest.stat().st_size > 1024):
            try:
                _stream_fastq_head_gz(r1_url, r1_dest, n_reads=350_000, retries=2)
                _stream_fastq_head_gz(r2_url, r2_dest, n_reads=350_000, retries=2)
            except Exception as exc:
                _log(f"EBI stream failed for {srr} ({exc}); falling back to NCBI SRA fastq-dump ...")
                r1_dest.unlink(missing_ok=True)
                r2_dest.unlink(missing_ok=True)
                tmp_srr = in_dir / f"_tmp_{srr}"
                shutil.rmtree(tmp_srr, ignore_errors=True)
                tmp_srr.mkdir(parents=True, exist_ok=True)
                subprocess.run(
                    [fq_dump, "-X", "350000", "--split-files", "--gzip", "-O", str(tmp_srr), srr],
                    check=True,
                )
                (tmp_srr / f"{srr}_1.fastq.gz").replace(r1_dest)
                (tmp_srr / f"{srr}_2.fastq.gz").replace(r2_dest)
                shutil.rmtree(tmp_srr, ignore_errors=True)
        payload_files.extend([r1_dest.name, r2_dest.name])
        urls.extend([r1_url, r2_url])

    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "07_cut_and_run_spikein",
            "aria_modality": "CUT_AND_RUN",
            "title": "Human U2OS SCAF1-AID Acute Auxin Degron CUT&RUN with Heterologous E. coli Spike-In",
            "accession": "GEO GSE273126 / ENA PRJNA1140214",
            "doi": "10.1016/j.molcel.2024.09.018",
            "citation": "GSE273126 — SCAF1-AID Auxin Degron HA CUT&RUN in Human U2OS Cells (hg38 + E. coli K-12 MG1655 spike-in)",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "3 Control vs 3 Auxin (1 uM, 6h) biological replicates + 1 Untagged negative control (first 350,000 paired-end reads/library preserving human + heterologous E. coli MG1655 spike-in reads)",
            "source_urls": urls,
            "payload_files": payload_files,
        },
    )
    return mod_dir


# ── Modality 08: CUT&Tag (Kaya-Okur et al. 2019 GSE124557 H1 vs K562) ─────────


def stage_08_cut_and_tag(atlas_root: Path) -> Path:
    mod_dir = atlas_root / "08_cut_and_tag"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)

    runs = [
        ("H1_H3K27me3_rep1_cuttag", "SRR838", "002", "SRR8383482"),
        ("H1_H3K27me3_rep2_cuttag", "SRR838", "003", "SRR8383483"),
        ("K562_H3K27me3_rep1_cuttag", "SRR838", "000", "SRR8383510"),
        ("K562_H3K27me3_rep2_cuttag", "SRR838", "001", "SRR8383511"),
        ("H1_IgG_ctrl_cuttag", "SRR843", "008", "SRR8435038"),
        ("K562_IgG_ctrl_cuttag", "SRR843", "001", "SRR8435051"),
    ]
    payload_files = []
    urls = []
    for sid, prefix, sub, srr in runs:
        for mate in (1, 2):
            url = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{prefix}/{sub}/{srr}/{srr}_{mate}.fastq.gz"
            fname = f"{sid}_R{mate}.fastq.gz"
            _download_file(url, in_dir / fname, timeout=240)
            payload_files.append(fname)
            urls.append(url)

    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "08_cut_and_tag",
            "aria_modality": "CUT_AND_TAG",
            "title": "Founding Human CUT&Tag Polycomb H3K27me3 Profiling (H1 hESC vs K562 + IgG)",
            "accession": "GEO GSE124557 / ENA PRJNA512492",
            "doi": "10.1038/s41467-019-09982-5",
            "citation": "Kaya-Okur HS, Wu SJ, Codomo CA, Pledger ES, Bryson TD, Henikoff JG, Ahmad K, Henikoff S, Nature Communications 10, 1930 (2019)",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "Full unmodified ENA paired-end FASTQ libraries: 2 H1 hESC vs 2 K562 H3K27me3 biological replicates + matched IgG negative controls",
            "source_urls": urls,
            "payload_files": payload_files,
        },
    )
    return mod_dir


# ── Modality 09: Hi-C 3D Genome (RPE-1 Diploid vs Tetraploid GSE345433) ───────


def stage_09_hic_3d_genome(atlas_root: Path) -> Path:
    mod_dir = atlas_root / "09_hic_3d_genome"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)

    base = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM10003nnn"
    samples = [
        ("Diploid_rep1.mcool", f"{base}/GSM10003420/suppl/GSM10003420_RPE_2N-MS02_hic_1.mcool"),
        ("Diploid_rep2.mcool", f"{base}/GSM10003421/suppl/GSM10003421_RPE_2N-MS02_hic_2.mcool"),
        ("Diploid_rep3.mcool", f"{base}/GSM10003422/suppl/GSM10003422_RPE_2N-MS1_hic_1.mcool"),
        ("Tetraploid_rep1.mcool", f"{base}/GSM10003424/suppl/GSM10003424_RPE_4N-MS02_hic_1.mcool"),
        ("Tetraploid_rep2.mcool", f"{base}/GSM10003425/suppl/GSM10003425_RPE_4N-MS02_hic_2.mcool"),
        ("Tetraploid_rep3.mcool", f"{base}/GSM10003426/suppl/GSM10003426_RPE_4N-MS1_hic_1.mcool"),
    ]
    payload_files = []
    urls = []
    for fname, url in samples:
        _download_file(url, in_dir / fname, timeout=420)
        payload_files.append(fname)
        urls.append(url)

    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "09_hic_3d_genome",
            "aria_modality": "HiC",
            "title": "Human RPE-1 Whole-Genome Duplication 3D Chromatin Organization (Diploid vs Tetraploid G1)",
            "accession": "GEO GSE345433",
            "doi": "NCBI GEO GSE345433",
            "citation": "GSE345433 — Whole genome duplication through mitotic slippage 3D Hi-C in human RPE-1 cells (GRCh38 .mcool)",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "3 Diploid_G1 vs 3 Tetraploid_G1 biological replicates (multi-resolution GRCh38 .mcool matrices)",
            "source_urls": urls,
            "payload_files": payload_files,
        },
    )
    return mod_dir


# ── Modality 10: Paired Same-Cell Multiome (GEO GSE278576 6-Donor Cohort) ─────


def stage_10_paired_multiome(atlas_root: Path) -> Path:
    mod_dir = atlas_root / "10_paired_multiome"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)
    raw_cache = atlas_root / "_raw_cache"

    donors = [
        ("hc29", "young_20_39", "M", 28),
        ("hc77", "young_20_39", "M", 20),
        ("hc76", "young_20_39", "F", 26),
        ("hc11", "aged_80_100", "M", 86),
        ("hc12", "aged_80_100", "M", 83),
        ("hc19", "aged_80_100", "F", 87),
    ]
    ref_h5ad = (
        Path.home()
        / "Samael"
        / "Oligos"
        / "data"
        / "raw"
        / "reference"
        / "GSE278576_hippocampus_RNA.h5ad"
    )
    src_map = {d[0]: str(_ensure_gse278576_donor_h5mu(d[0], raw_cache)) for d in donors}

    py_script = f"""
import json
from pathlib import Path
import numpy as np
import mudata as md
import anndata as ad

donors = {json.dumps(donors)}
src_map = {json.dumps(src_map)}
keep_chroms = set({json.dumps(list(ATLAS_CHROMS))})
in_dir = Path("{in_dir}")
ref_path = Path("{ref_h5ad}")

label_lookup = {{}}
if ref_path.is_file():
    ref_ad = ad.read_h5ad(str(ref_path), backed="r")
    obs_df = ref_ad.obs[["orig.ident", "subclass"]].copy()
    for idx_val, row in obs_df.iterrows():
        bc = str(idx_val).split("_")[-1].split("#")[-1]
        label_lookup[(str(row["orig.ident"]), bc)] = str(row["subclass"])

rng = np.random.default_rng(2026)
for donor_id, age_group, sex, age in donors:
    out_p = in_dir / f"{{donor_id}}_paired_multiome.h5mu"
    if out_p.is_file():
        continue
    mdata = md.read_h5mu(src_map[donor_id])
    rna = mdata.mod["rna"].copy()
    atac = mdata.mod["atac"].copy()
    shared = rna.obs_names.intersection(atac.obs_names)
    rna = rna[shared, :].copy()
    atac = atac[shared, :].copy()
    peak_chroms = [str(p).split(":")[0].split("-")[0] for p in atac.var_names]
    peak_mask = np.array([c in keep_chroms for c in peak_chroms], dtype=bool)
    atac = atac[:, peak_mask].copy()
    r_tot = np.asarray(rna.X.sum(axis=1)).ravel()
    a_tot = np.asarray(atac.X.sum(axis=1)).ravel()
    valid = np.where((r_tot >= 200) & (a_tot >= 100))[0]
    chosen = np.sort(rng.choice(valid, size=min(400, len(valid)), replace=False))
    rna = rna[chosen, :].copy()
    atac = atac[chosen, :].copy()
    new_ids = [f"{{donor_id}}#{{str(b).split('#')[-1]}}" for b in rna.obs_names]
    subclasses = [
        label_lookup.get((donor_id, str(b).split("#")[-1]), "Hippocampal_Cell")
        for b in rna.obs_names
    ]
    rna.obs_names = new_ids
    atac.obs_names = new_ids
    for ad_mod in (rna, atac):
        ad_mod.obs["sample_id"] = donor_id
        ad_mod.obs["donor_id"] = donor_id
        ad_mod.obs["age_group"] = age_group
        ad_mod.obs["condition"] = age_group
        ad_mod.obs["replicate"] = donor_id
        ad_mod.obs["sex"] = sex
        ad_mod.obs["age"] = int(age)
        ad_mod.obs["subclass"] = subclasses
    out_md = md.MuData({{"rna": rna, "atac": atac}})
    for col in ("sample_id", "donor_id", "age_group", "condition", "replicate", "sex", "subclass"):
        out_md.obs[col] = rna.obs[col].values
    out_md.write_h5mu(str(out_p))
"""
    _conda_run("aria-chromatin-env", ["python", "-c", py_script])

    payload_files = [f"{d[0]}_paired_multiome.h5mu" for d in donors]
    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "10_paired_multiome",
            "aria_modality": "Paired_Multiome",
            "title": "Human Hippocampus Same-Cell 10x Multiome RNA+ATAC Cohort (3v3 Donors)",
            "accession": "GEO GSE278576",
            "doi": "10.1126/science.adt8307",
            "citation": "Zemke NR et al., Science (2026) doi:10.1126/science.adt8307",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "6 post-mortem human hippocampal donors (3 young_20_39 [hc29, hc77, hc76] vs 3 aged_80_100 [hc11, hc12, hc19], 400 paired same-cell nuclei/donor = 2,400 paired nuclei across chr1/6/7/11/17/19)",
            "source_urls": [
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278576"
            ],
            "payload_files": payload_files,
        },
    )
    return mod_dir


# ── Modality 11: Cross-Modal RNA x ATAC (ENCODE GM12878 vs K562 3v3 Matched) ──


def stage_11_cross_modal_concordance(atlas_root: Path, refs: dict[str, Path]) -> Path:
    mod_dir = atlas_root / "11_cross_modal_concordance"
    if _is_ready(mod_dir):
        return mod_dir
    in_dir = _inputs_dir(mod_dir)

    atac_dir = _inputs_dir(stage_05_bulk_atac_tobias(atlas_root))
    sample_order = [
        "GM12878_rep1",
        "GM12878_rep2",
        "GM12878_rep3",
        "K562_rep1",
        "K562_rep2",
        "K562_rep3",
    ]
    atac_bams = []
    for sid in sample_order:
        bname = f"{sid}.bam"
        atac_bams.append(bname)
        for ext in ("", ".bai"):
            src = atac_dir / f"{sid}_atac.bam{ext}"
            dst = in_dir / f"{bname}{ext}"
            if not dst.exists():
                try:
                    os.link(src, dst)
                except OSError:
                    shutil.copy2(src, dst)

    raw_rna = atlas_root / "_raw_cache" / "encode_rna"
    rna_tsvs = [
        ("GM12878_rep1", "ENCFF700JWD"),
        ("GM12878_rep2", "ENCFF935TRK"),
        ("GM12878_rep3", "ENCFF480HNH"),
        ("K562_rep1", "ENCFF643FJJ"),
        ("K562_rep2", "ENCFF939LDL"),
        ("K562_rep3", "ENCFF241QPH"),
    ]
    tsv_paths: dict[str, Path] = {}
    urls = []
    for sid, encff in rna_tsvs:
        url = f"https://www.encodeproject.org/files/{encff}/@@download/{encff}.tsv"
        urls.append(url)
        tsv_paths[sid] = _download_file(url, raw_rna / f"{encff}.tsv", timeout=180)

    # Filter RNA counts to BULK_CHROMS via the UCSC GRCh38 GTF so RNA and ATAC
    # share the exact 5-chromosome regulatory universe (chr6, chr7, chr11, chr17, chr19).
    gtf_path = refs["hg38_chr_gtf"]
    keep_bulk = set(BULK_CHROMS)
    ensg_to_sym: dict[str, str] = {}
    gene_id_re = re.compile(r'gene_id "([^"]+)"')
    gene_name_re = re.compile(r'gene_name "([^"]+)"')
    with gtf_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "gene" or parts[0] not in keep_bulk:
                continue
            m_id = gene_id_re.search(parts[8])
            m_nm = gene_name_re.search(parts[8])
            if m_id and m_nm:
                ensg_to_sym[m_id.group(1).split(".")[0]] = m_nm.group(1)

    gene_counts: dict[str, list[int]] = {}
    for s_idx, sid in enumerate(sample_order):
        with tsv_paths[sid].open("r", encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            gid_col = header.index("gene_id")
            cnt_col = header.index("expected_count")
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                raw_gid = parts[gid_col].split(".")[0]
                if not raw_gid.startswith("ENSG"):
                    continue
                sym = ensg_to_sym.get(raw_gid)
                if not sym:
                    continue
                cnt = int(round(float(parts[cnt_col])))
                vec = gene_counts.setdefault(sym, [0] * len(sample_order))
                vec[s_idx] += cnt

    counts_out = in_dir / "encode_gm12878_vs_k562_counts.tsv"
    with counts_out.open("w", encoding="utf-8") as out:
        out.write("gene\t" + "\t".join(sample_order) + "\n")
        for sym in sorted(gene_counts):
            if sum(gene_counts[sym]) == 0:
                continue
            out.write(sym + "\t" + "\t".join(str(v) for v in gene_counts[sym]) + "\n")

    _write_ready_manifest(
        mod_dir,
        {
            "modality_id": "11_cross_modal_concordance",
            "aria_modality": "Cross_Modal_RNA_x_ATAC",
            "title": "ENCODE Matched Human Hematopoietic polyA+ RNA-seq x Bulk ATAC-seq (GM12878 vs K562)",
            "accession": "ENCODE ENCSR843RJV, ENCSR889TRN, ENCSR594NJP, ENCSR596ACL, ENCSR637XSC, ENCSR483RKN, ENCSR859USB",
            "doi": "10.1038/s41586-020-2493-4",
            "citation": "ENCODE Project Consortium, Nature 583, 699-710 (2020)",
            "organism": "Homo sapiens",
            "genome": "hg38",
            "design_summary": "N=6 matched biological units (GM12878_rep1..3 vs K562_rep1..3) profiled across both bulk polyA+ RNA-seq and bulk ATAC-seq on GRCh38 chr6, chr7, chr11, chr17, chr19",
            "source_urls": urls,
            "payload_files": ["encode_gm12878_vs_k562_counts.tsv", *atac_bams],
        },
    )
    return mod_dir


def provision_modality(mod_id: str, atlas_root: Path, refs: dict[str, Path] | None = None) -> Path:
    if mod_id == "01_bulk_rna_matrix":
        return stage_01_bulk_rna_matrix(atlas_root)
    if mod_id == "02_bulk_rna_raw_fastq":
        return stage_02_bulk_rna_raw_fastq(atlas_root)
    if mod_id == "03_scrna_multidonor":
        return stage_03_scrna_multidonor(atlas_root)
    if mod_id == "04_scatac_multidonor":
        return stage_04_scatac_multidonor(atlas_root)
    if mod_id == "05_bulk_atac_tobias":
        return stage_05_bulk_atac_tobias(atlas_root)
    if mod_id == "06_chip_seq":
        return stage_06_chip_seq(atlas_root)
    if mod_id == "07_cut_and_run_spikein":
        return stage_07_cut_and_run_spikein(atlas_root)
    if mod_id == "08_cut_and_tag":
        return stage_08_cut_and_tag(atlas_root)
    if mod_id == "09_hic_3d_genome":
        return stage_09_hic_3d_genome(atlas_root)
    if mod_id == "10_paired_multiome":
        return stage_10_paired_multiome(atlas_root)
    if mod_id == "11_cross_modal_concordance":
        if refs is None:
            refs = ensure_references(atlas_root)
        return stage_11_cross_modal_concordance(atlas_root, refs)
    raise ValueError(f"Unknown modality_id: {mod_id}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ATLAS_ROOT,
        help="Target ext4 directory for public atlas datasets.",
    )
    ap.add_argument(
        "--modality",
        default="all",
        choices=["all", "references", *ALL_MODALITY_IDS],
        help="Specific modality ID to stage, 'references', or 'all'.",
    )
    args = ap.parse_args()

    atlas_root = args.root.expanduser().resolve()
    atlas_root.mkdir(parents=True, exist_ok=True)

    refs = ensure_references(atlas_root)
    if args.modality == "references":
        _log("Shared references ready.")
        return 0

    targets = ALL_MODALITY_IDS if args.modality == "all" else (args.modality,)
    for mid in targets:
        t0 = time.time()
        _log(f"=== Provisioning {mid} ===")
        p = provision_modality(mid, atlas_root, refs=refs)
        _log(f"=== Ready: {mid} at {p} ({time.time() - t0:.1f}s) ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
