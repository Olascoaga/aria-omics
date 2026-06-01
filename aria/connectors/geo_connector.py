"""
ARIA GEO/SRA Connector (v4.3)
------------------------------
Downloads and preprocesses public omics datasets from NCBI GEO and SRA.

Supported accession formats:
  GSExxxxxx  — GEO Series (most common: processed count matrices)
  SRPxxxxxx  — SRA Study (raw FASTQs via pysradb)
  PRJNAxxxxxx — BioProject (mapped to SRA)

Strategy:
  1. Fetch SOFT file for sample metadata (organism, groups, characteristics)
  2. Download supplementary processed files if available (count matrices, h5ad, etc.)
  3. Fall back to SRA raw FASTQs if no processed data found (requires pysradb)
  4. Infer experimental design from sample characteristics

Usage:
    from aria.connectors.geo_connector import GEOConnector
    gc = GEOConnector()
    result = gc.fetch("GSE183948", status_callback=print)
    # result["local_dir"]       → Path with downloaded files
    # result["inferred_design"] → groups, organism, factor
"""

from __future__ import annotations

import gzip
import io
import logging
import re
import shutil
import urllib.request
from pathlib import Path
from typing import Callable, Optional

log = logging.getLogger("aria.geo")

# GEO FTP base URL
_GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

# Organism → (canonical name, genome assembly)
_GENOME_MAP: dict[str, tuple[str, str]] = {
    "homo sapiens":                ("Homo sapiens",           "hg38"),
    "mus musculus":                ("Mus musculus",           "mm39"),
    "rattus norvegicus":           ("Rattus norvegicus",       "rn7"),
    "danio rerio":                 ("Danio rerio",             "danRer11"),
    "drosophila melanogaster":     ("Drosophila melanogaster", "dm6"),
    "caenorhabditis elegans":      ("Caenorhabditis elegans",  "ce11"),
    "saccharomyces cerevisiae":    ("Saccharomyces cerevisiae","sacCer3"),
    "arabidopsis thaliana":        ("Arabidopsis thaliana",    "TAIR10"),
    "macaca mulatta":              ("Macaca mulatta",          "rheMac10"),
    "sus scrofa":                  ("Sus scrofa",              "susScr11"),
}

# Supplementary file extensions that indicate processed data
_COUNT_PATTERNS  = re.compile(
    r"(count|raw|matrix|expression|tpm|fpkm|rpkm)", re.IGNORECASE
)
_SKIP_PATTERNS   = re.compile(
    r"(readme|supplementary|supp_note|figures?|table|filelist"
    r"|deseq2?|edger|limma|_significant|_sig\.|fold.?change|normalized|norm\.)",
    re.IGNORECASE
)


class GEOConnector:
    """
    Downloads and parses GEO/SRA datasets into a local cache directory.
    """

    def __init__(self, cache_dir: Optional[str] = None):
        self.cache_dir = Path(cache_dir or Path.home() / ".aria" / "geo_cache")
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    # ── Public API ────────────────────────────────────────────────────────

    def fetch(self, accession: str,
              status_callback: Optional[Callable[[str], None]] = None) -> dict:
        """
        Main entry point.

        Returns:
            {
              "accession":       str,
              "title":           str,
              "organism":        str,
              "genome":          str,
              "data_type":       "bulk_RNA" | "scRNA" | "unknown",
              "n_samples":       int,
              "local_dir":       Path,
              "files":           {"counts": [], "h5ad": [], "h5": [], "mtx": []},
              "inferred_design": {
                  "groups":  {group: [sample_names]},
                  "factor":  str,
                  "organism": str,
                  "genome":   str,
              },
              "geo_metadata":    dict,
            }
        """
        # W-PRIV (P1-7/P1-8): fetching a public accession is network egress to
        # NCBI; refuse it under air-gapped mode instead of leaking the request.
        from aria.utils.privacy import assert_egress_allowed
        assert_egress_allowed("GEO/SRA")

        acc = accession.strip().upper()
        if acc.startswith("GSE"):
            return self._fetch_gse(acc, status_callback)
        elif acc.startswith(("SRP", "PRJNA", "ERP", "DRP")):
            return self._fetch_sra(acc, status_callback)
        else:
            raise ValueError(
                f"Unrecognised accession '{accession}'. "
                "Expected GSExxxxxx, SRPxxxxxx, or PRJNAxxxxxx."
            )

    # ── GEO path ──────────────────────────────────────────────────────────

    def _fetch_gse(self, gse_id: str,
                   status_cb: Optional[Callable]) -> dict:
        local_dir = self.cache_dir / gse_id
        local_dir.mkdir(parents=True, exist_ok=True)

        _cb(status_cb, f"[GEO] Fetching metadata for {gse_id}...")
        metadata = self._parse_soft(gse_id, local_dir, status_cb)

        _cb(status_cb, f"[GEO] Downloading supplementary files...")
        files = self._download_supplementary(gse_id, local_dir, status_cb)

        organism, genome = _resolve_organism(metadata.get("organism", ""))

        # When multiple organisms appear in the SOFT (e.g. spike-in experiments),
        # use gene symbol style in the count matrix to identify the experimental organism.
        sample_organisms = {
            s.get("organism", "").lower()
            for s in metadata.get("samples", [])
            if s.get("organism")
        }
        if len(sample_organisms) > 1:
            count_files = files.get("counts", []) or files.get("h5ad", [])
            if count_files:
                sym_org = _organism_from_gene_symbols(count_files[0])
                if sym_org:
                    sym_resolved, sym_genome = _resolve_organism(sym_org)
                    if sym_resolved.lower() != organism.lower():
                        _cb(status_cb,
                            f"[GEO] Multiple organisms in metadata ({', '.join(sample_organisms)}); "
                            f"gene symbols in matrix suggest {sym_resolved} — using that.")
                        organism, genome = sym_resolved, sym_genome

        design = _infer_design(metadata)
        design["organism"] = organism
        design["genome"]   = genome
        data_type = _infer_data_type(metadata, files)

        return {
            "accession":       gse_id,
            "title":           metadata.get("title", ""),
            "organism":        organism,
            "genome":          genome,
            "data_type":       data_type,
            "n_samples":       len(metadata.get("samples", [])),
            "local_dir":       str(local_dir),
            "files":           files,
            "inferred_design": design,
            "geo_metadata":    metadata,
        }

    def _parse_soft(self, gse_id: str, local_dir: Path,
                    status_cb: Optional[Callable]) -> dict:
        """Parse GEO SOFT file. Tries GEOparse first, falls back to manual."""
        try:
            import GEOparse
            _cb(status_cb, f"[GEO] Parsing SOFT via GEOparse...")
            gse = GEOparse.get_GEO(geo=gse_id, destdir=str(local_dir),
                                    silent=True)
            return _soft_from_geoparse(gse)
        except ImportError:
            pass
        except Exception as e:
            log.warning(f"GEOparse failed ({e}), using manual parser.")

        return self._parse_soft_manual(gse_id, local_dir, status_cb)

    def _parse_soft_manual(self, gse_id: str, local_dir: Path,
                            status_cb: Optional[Callable]) -> dict:
        """Manual SOFT file parser (no GEOparse dependency)."""
        prefix   = _gse_prefix(gse_id)
        soft_url = (f"{_GEO_FTP}/{prefix}/{gse_id}/soft/"
                    f"{gse_id}_family.soft.gz")
        soft_path = local_dir / f"{gse_id}_family.soft.gz"

        if not soft_path.exists():
            _cb(status_cb, f"[GEO] Downloading SOFT file...")
            try:
                urllib.request.urlretrieve(soft_url, str(soft_path))
            except Exception as e:
                log.warning(f"SOFT download failed: {e}")
                return {"title": gse_id, "organism": "", "samples": [],
                        "suppl_files": []}

        try:
            with gzip.open(str(soft_path), "rt", encoding="utf-8",
                           errors="replace") as fh:
                text = fh.read()
        except Exception as e:
            log.warning(f"SOFT read failed: {e}")
            return {"title": gse_id, "organism": "", "samples": [],
                    "suppl_files": []}

        return _parse_soft_text(text)

    def _download_supplementary(self, gse_id: str, local_dir: Path,
                                  status_cb: Optional[Callable]) -> dict:
        """
        List and download supplementary files from GEO FTP.
        Classifies files into counts, h5ad, h5, mtx buckets.
        """
        files: dict = {"counts": [], "h5ad": [], "h5": [], "mtx": []}
        prefix  = _gse_prefix(gse_id)
        ftp_url = f"{_GEO_FTP}/{prefix}/{gse_id}/suppl/"

        try:
            with urllib.request.urlopen(ftp_url, timeout=30) as r:
                html = r.read().decode("utf-8", errors="replace")
        except Exception as e:
            log.warning(f"Could not list supplementary files: {e}")
            return files

        # Extract bare filenames from FTP HTML listing.
        # The NCBI FTP via HTTP returns an Apache-style index; href values are
        # plain filenames (no slashes, no scheme).  Filter out anything that
        # looks like a relative path, absolute URL, or navigation link.
        fnames = re.findall(r'href="([^"]+)"', html)
        fnames = [
            f for f in fnames
            if "." in f
            and "/" not in f
            and "://" not in f
            and not f.startswith(("?", ".."))
        ]

        for fname in fnames:
            if _SKIP_PATTERNS.search(fname):
                continue

            fpath = local_dir / fname
            if not fpath.exists():
                _cb(status_cb, f"[GEO] Downloading {fname}...")
                try:
                    urllib.request.urlretrieve(ftp_url + fname, str(fpath))
                except Exception as e:
                    log.warning(f"Download failed for {fname}: {e}")
                    continue

            fl = fname.lower()
            if fl.endswith(".h5ad"):
                files["h5ad"].append(str(fpath))
            elif fl.endswith(".h5"):
                files["h5"].append(str(fpath))
            elif "matrix.mtx" in fl or fl.endswith(".mtx.gz"):
                files["mtx"].append(str(fpath))
            elif re.search(r"\.(txt|tsv|csv)(\.gz)?$", fl):
                files["counts"].append(str(fpath))

        if not any(files.values()):
            log.warning(f"No processable supplementary files found for {gse_id}.")

        return files

    # ── SRA path ──────────────────────────────────────────────────────────

    def _fetch_sra(self, sra_id: str,
                   status_cb: Optional[Callable]) -> dict:
        """
        Fetch metadata and optionally download FASTQs via pysradb.
        """
        try:
            from pysradb.sraweb import SRAweb
        except ImportError:
            raise ImportError(
                "pysradb is required for SRA downloads. "
                "Install with: pip install pysradb"
            )

        _cb(status_cb, f"[SRA] Fetching metadata for {sra_id}...")
        db = SRAweb()

        try:
            df = db.sra_metadata(sra_id, detailed=True)
        except Exception as e:
            raise RuntimeError(f"SRA metadata fetch failed: {e}") from e

        if df is None or df.empty:
            raise RuntimeError(f"No metadata found for {sra_id}")

        local_dir = self.cache_dir / sra_id
        local_dir.mkdir(parents=True, exist_ok=True)

        organism = df["organism_name"].iloc[0] if "organism_name" in df.columns else ""
        org_canonical, genome = _resolve_organism(organism)

        # Build metadata dict compatible with the rest of the pipeline
        samples = []
        for _, row in df.iterrows():
            samples.append({
                "id":    str(row.get("run_accession", row.get("experiment_accession", ""))),
                "title": str(row.get("sample_title", row.get("experiment_title", ""))),
                "organism": organism,
                "characteristics": {
                    "source": str(row.get("source_name", "")),
                    "treatment": str(row.get("treatment", row.get("condition", ""))),
                },
                "library_strategy": str(row.get("library_strategy", "")),
            })

        metadata = {
            "title":    sra_id,
            "organism": organism,
            "samples":  samples,
            "suppl_files": [],
            "library_strategy": df["library_strategy"].iloc[0]
                                 if "library_strategy" in df.columns else "",
        }

        design    = _infer_design(metadata)
        design.update({"organism": org_canonical, "genome": genome})
        data_type = _infer_data_type(metadata, {})

        # Save metadata CSV
        df.to_csv(str(local_dir / f"{sra_id}_metadata.csv"), index=False)

        _cb(status_cb,
            f"[SRA] Found {len(samples)} runs. "
            f"To download FASTQs, run: pysradb download --srp {sra_id}")

        return {
            "accession":       sra_id,
            "title":           metadata["title"],
            "organism":        org_canonical,
            "genome":          genome,
            "data_type":       data_type,
            "n_samples":       len(samples),
            "local_dir":       str(local_dir),
            "files":           {"counts": [], "h5ad": [], "h5": [], "mtx": [],
                                "fastq_pending": True},
            "inferred_design": design,
            "geo_metadata":    metadata,
            "sra_metadata_csv": str(local_dir / f"{sra_id}_metadata.csv"),
        }


# ── Module-level helpers ──────────────────────────────────────────────────────

def _cb(fn: Optional[Callable], msg: str):
    if fn:
        fn(msg)
    log.info(msg)


def _gse_prefix(gse_id: str) -> str:
    """GSE183948 → GSE183nnn"""
    return gse_id[:-3] + "nnn"


def _resolve_organism(raw: str) -> tuple[str, str]:
    key = raw.strip().lower()
    for pattern, (name, asm) in _GENOME_MAP.items():
        if pattern in key:
            return name, asm
    return raw.strip() or "Unknown", "unknown"


def _parse_soft_text(text: str) -> dict:
    """Parse raw SOFT file text into a metadata dict."""
    metadata: dict = {
        "title": "", "organism": "",
        "samples": [], "suppl_files": [],
        "library_strategy": "",
    }
    current: Optional[dict] = None

    for line in text.splitlines():
        line = line.strip()
        val  = line.split(" = ", 1)[-1].strip() if " = " in line else ""

        if line.startswith("^SAMPLE"):
            if current:
                metadata["samples"].append(current)
            current = {"id": val, "title": "", "organism": "",
                       "characteristics": {}, "library_strategy": ""}

        elif line.startswith("!Series_title"):
            metadata["title"] = val
        elif line.startswith("!Series_supplementary_file") and val not in ("", "NONE"):
            metadata["suppl_files"].append(val)

        elif current is not None:
            if line.startswith("!Sample_title"):
                current["title"] = val
            elif line.startswith("!Sample_organism_ch1"):
                current["organism"] = val
                if not metadata["organism"]:
                    metadata["organism"] = val
            elif line.startswith("!Sample_characteristics_ch1") and ":" in val:
                k, v = val.split(":", 1)
                current["characteristics"][k.strip().lower()] = v.strip()
            elif line.startswith("!Sample_library_strategy"):
                current["library_strategy"] = val
                if not metadata["library_strategy"]:
                    metadata["library_strategy"] = val

    if current:
        metadata["samples"].append(current)

    return metadata


def _soft_from_geoparse(gse) -> dict:
    """Convert GEOparse GSE object to our metadata dict."""
    meta = gse.metadata
    title = (meta.get("title", [""])[0]
             if isinstance(meta.get("title"), list)
             else meta.get("title", ""))
    suppl = (meta.get("supplementary_file", [])
             if isinstance(meta.get("supplementary_file"), list)
             else [])

    samples = []
    for gsm_id, gsm in gse.gsms.items():
        gm = gsm.metadata
        organism = (gm.get("organism_ch1", [""])[0]
                    if isinstance(gm.get("organism_ch1"), list)
                    else gm.get("organism_ch1", ""))
        chars = {}
        for ch in gm.get("characteristics_ch1", []):
            if ":" in ch:
                k, v = ch.split(":", 1)
                chars[k.strip().lower()] = v.strip()
        lib_strat = (gm.get("library_strategy", [""])[0]
                     if isinstance(gm.get("library_strategy"), list)
                     else gm.get("library_strategy", ""))
        title_s = (gm.get("title", [""])[0]
                   if isinstance(gm.get("title"), list)
                   else gm.get("title", gsm_id))
        samples.append({
            "id": gsm_id, "title": title_s,
            "organism": organism,
            "characteristics": chars,
            "library_strategy": lib_strat,
        })

    organism = samples[0]["organism"] if samples else ""
    lib_strat = samples[0]["library_strategy"] if samples else ""

    return {
        "title": title,
        "organism": organism,
        "library_strategy": lib_strat,
        "samples": samples,
        "suppl_files": suppl,
    }


def _infer_design(metadata: dict) -> dict:
    """
    Infer experimental groups from sample characteristics.
    Prefer experimental-design keys over incidental sample attributes and use
    GEO/SRA accessions as the canonical sample IDs.  Titles are kept as aliases
    so downstream count matrices can match either GSM/SRR columns or readable
    sample labels.
    """
    samples = metadata.get("samples", [])
    if not samples:
        return {"groups": {}, "factor": "condition", "main_factor": "condition",
                "condition_col": "condition", "n_groups": 0}

    all_keys: set = set()
    for s in samples:
        all_keys.update(s.get("characteristics", {}).keys())

    best_key: Optional[str] = None
    best_score = -1.0
    n = len(samples)
    priority = {
        "condition": 100,
        "treatment": 95,
        "group": 90,
        "experimental group": 90,
        "genotype": 85,
        "perturbation": 80,
        "knockout": 80,
        "stim": 78,
        "stimulation": 78,
        "timepoint": 70,
        "time point": 70,
        "dose": 65,
        "disease": 60,
        "phenotype": 55,
        "source": 25,
        "source_name": 25,
    }

    for key in all_keys:
        vals = [
            str(s.get("characteristics", {}).get(key, "")).strip()
            for s in samples
        ]
        vals = [v for v in vals if v]
        unique = len(set(vals))
        if unique <= 1 or unique == n:
            continue
        # Prefer a balanced number of groups, but make semantic design keys
        # outrank short incidental keys such as sex.
        balance = 1.0 - abs(unique - 2) / max(2, n)
        semantic = priority.get(key.lower(), 0)
        if semantic == 0:
            for token, weight in priority.items():
                if token in key.lower():
                    semantic = max(semantic, weight)
        score = semantic + (10 * balance) - (0.05 * len(key))
        if score > best_score:
            best_score = score
            best_key = key

    if not best_key:
        # Fall back: any key with > 1 unique value
        for key in sorted(all_keys):
            vals = [
                str(s.get("characteristics", {}).get(key, "")).strip()
                for s in samples
            ]
            unique = len(set(v for v in vals if v))
            if 1 < unique < n:
                best_key = key
                break

    groups: dict = {}
    sample_aliases: dict = {}
    if best_key:
        for s in samples:
            val = s.get("characteristics", {}).get(best_key, "unknown")
            # Sanitise group name
            val = re.sub(r"[^a-zA-Z0-9_\-]", "_", val).strip("_") or "group"
            sample_id = str(s.get("id") or s.get("title") or "sample")
            title = str(s.get("title") or "")
            groups.setdefault(val, []).append(sample_id)
            aliases = [sample_id]
            if title and title != sample_id:
                aliases.append(title)
            sample_aliases[sample_id] = aliases
    else:
        for s in samples:
            sample_id = str(s.get("id") or s.get("title") or "sample")
            title = str(s.get("title") or "")
            groups.setdefault("all_samples", []).append(sample_id)
            sample_aliases[sample_id] = [sample_id] + (
                [title] if title and title != sample_id else []
            )

    return {
        "groups":   groups,
        "factor":   best_key or "condition",
        "main_factor": best_key or "condition",
        "condition_col": best_key or "condition",
        "n_groups": len(groups),
        "sample_aliases": sample_aliases,
        "source": "GEO/SRA metadata",
        "confidence": "high" if best_key and len(groups) > 1 else "low",
        "reasoning": (
            f"Groups inferred from GEO/SRA sample characteristic '{best_key}'."
            if best_key else
            "No multi-level GEO/SRA characteristic was suitable for grouping."
        ),
    }


def _organism_from_gene_symbols(count_path: str) -> str:
    """
    Peek at gene IDs in a count matrix (first column, first 200 rows) to
    infer organism by symbol style.  Returns a lowercase string suitable
    for _resolve_organism(), or "" if uncertain.

    Patterns:
      Human (HGNC): ALL-CAPS 2-10 chars — CXCL1, IL6, TP53, GAPDH
      Mouse  (MGI):  Title-case 2-10 chars — Cxcl1, Il6, Trp53
      Fly (FlyBase): CG numbers — CG1234, CG9870
      Ensembl:       ENSG / ENSMUSG / FBgn prefixes
    """
    try:
        path = Path(count_path)
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            fh.readline()          # skip header
            genes = []
            for _ in range(200):
                line = fh.readline()
                if not line:
                    break
                sep  = "\t" if "\t" in line else ","
                gene = line.split(sep)[0].strip().strip('"').strip("'")
                if gene:
                    genes.append(gene)
    except Exception:
        return ""

    if not genes:
        return ""

    n = len(genes)

    # Ensembl / FlyBase IDs are unambiguous
    if sum(1 for g in genes if g.startswith("ENSG"))   / n > 0.3: return "homo sapiens"
    if sum(1 for g in genes if g.startswith("ENSMUS")) / n > 0.3: return "mus musculus"
    if sum(1 for g in genes if g.startswith("FBgn"))   / n > 0.3: return "drosophila melanogaster"
    if sum(1 for g in genes if re.match(r"^CG\d+$", g))/ n > 0.2: return "drosophila melanogaster"

    # Symbol-style heuristics
    human_like = sum(1 for g in genes if re.match(r"^[A-Z][A-Z0-9][A-Z0-9\-\.]{0,8}$", g))
    mouse_like = sum(1 for g in genes if re.match(r"^[A-Z][a-z][a-z0-9\-]{0,8}$",      g))

    if human_like / n > 0.4: return "homo sapiens"
    if mouse_like / n > 0.4: return "mus musculus"

    return ""


def _infer_data_type(metadata: dict, files: dict) -> str:
    """Infer bulk_RNA, scRNA, or unknown from metadata and file names."""
    lib = metadata.get("library_strategy", "").lower()
    title_desc = (metadata.get("title", "") + " " +
                  " ".join(s.get("title", "") for s in metadata.get("samples", []))
                  ).lower()

    # scRNA indicators
    if any(kw in title_desc or kw in lib
           for kw in ("scrna", "single cell", "single-cell", "10x", "droplet",
                      "smart-seq", "smartseq", "drop-seq")):
        return "scRNA"
    if files.get("h5ad") or files.get("h5") or files.get("mtx"):
        return "scRNA"

    # Bulk RNA indicators
    if "rna-seq" in lib or "rna seq" in title_desc:
        return "bulk_RNA"
    if files.get("counts"):
        return "bulk_RNA"

    return "unknown"
