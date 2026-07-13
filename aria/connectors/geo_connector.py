"""
ARIA GEO/SRA Connector
----------------------
Atomically retrieves and validates public omics datasets from NCBI GEO/SRA.

Supported accession formats:
  GSExxxxxx  — GEO Series (most common: processed count matrices)
  SRP/ERP/DRPxxxxxx — SRA/ENA/DRA Study (raw FASTQs via SRA Toolkit)
  PRJNAxxxxxx — BioProject (mapped to SRA)

Strategy:
  1. Fetch SOFT file for sample metadata (organism, groups, characteristics)
  2. Download supplementary processed files if available (count matrices, h5ad, etc.)
  3. Safely extract tar/tgz/zip payloads and validate their content
  4. Fall back to official SRA RunInfo + prefetch/vdb-validate/fasterq-dump
  5. Publish only a complete, hash-manifested accession generation
  6. Infer experimental design from sample characteristics

Usage:
    from aria.connectors.geo_connector import GEOConnector
    gc = GEOConnector()
    result = gc.fetch("GSE183948", status_callback=print)
    # result["local_dir"]       → Path with downloaded files
    # result["inferred_design"] → groups, organism, factor
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import logging
import os
import re
import signal
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path
from typing import Callable, Optional
from urllib.parse import quote, unquote, urlencode

from aria.utils.atomic_retrieval import (
    MANIFEST_NAME,
    RetrievalError,
    decompress_gzip_atomic,
    download_atomic,
    open_url_with_retry,
    publish_directory_atomic,
    safe_extract_archive,
    validate_payload,
    validate_retrieval_manifest,
    write_retrieval_manifest,
)

log = logging.getLogger("aria.geo")

# GEO FTP base URL
_GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
_NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_RESULT_NAME = "retrieval_result.json"

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
              "data_type":       "bulk_RNA" | "scRNA" | "bulk_ATAC"
                                 | "scATAC" | "unknown",
              "n_samples":       int,
              "local_dir":       Path,
              "files":           {"counts": [], "h5ad": [], "h5": [], "mtx": [],
                                  "fragments": [], "peaks": [], "bam": []},
              "inferred_design": {
                  "groups":  {group: [sample_names]},
                  "factor":  str,
                  "organism": str,
                  "genome":   str,
              },
              "geo_metadata":    dict,
            }
        """
        acc = accession.strip().upper()
        if not re.fullmatch(
            r"(?:GSE|SRP|ERP|DRP|PRJNA|PRJEB|PRJDB)\d+", acc
        ):
            raise ValueError(
                f"Unrecognised accession '{accession}'. "
                "Expected GSExxxxxx, SRP/ERP/DRPxxxxxx, or PRJN/EB/DBxxxxxx."
            )

        # A fully hash-validated local generation performs no egress and is safe
        # to reuse even after the experiment switches to air-gapped mode.
        cached = self._load_cached_result(acc)
        if cached is not None:
            _cb(status_callback, f"[GEO/SRA] Reusing validated cache for {acc}.")
            return cached

        # W-PRIV (P1-7/P1-8): an uncached accession requires network egress.
        from aria.utils.privacy import assert_egress_allowed
        assert_egress_allowed("GEO/SRA")

        if acc.startswith("GSE"):
            return self._fetch_gse(acc, status_callback)
        if acc.startswith(("SRP", "PRJ", "ERP", "DRP")):
            return self._fetch_sra(acc, status_callback)
        raise AssertionError(f"validated accession was not routed: {acc}")

    # ── GEO path ──────────────────────────────────────────────────────────

    def _fetch_gse(self, gse_id: str,
                   status_cb: Optional[Callable]) -> dict:
        target = self.cache_dir / gse_id
        staging = self._new_staging_dir(gse_id)
        try:
            _cb(status_cb, f"[GEO] Fetching metadata for {gse_id}...")
            metadata = self._parse_soft(gse_id, staging, status_cb)

            _cb(status_cb, "[GEO] Downloading supplementary files...")
            files = self._download_supplementary(gse_id, staging, status_cb)
            file_modalities: dict[str, str] = {}
            source_accessions = [gse_id]
            fallback_bundle = None
            if not _has_payload(files):
                related = _preferred_sra_accessions(
                    metadata.get("sra_accessions") or []
                )
                if not related:
                    raise RetrievalError(
                        f"{gse_id} has no validated analyzable supplementary "
                        "payload and no related SRA/BioProject accession"
                    )
                _cb(
                    status_cb,
                    "[GEO] Falling back to SRA accession(s) "
                    + ", ".join(related)
                    + "...",
                )
                fallback_bundle = self._retrieve_sra_bundle(
                    related, staging, status_cb
                )
                files = fallback_bundle["files"]
                file_modalities = fallback_bundle["file_modalities"]
                source_accessions.extend(fallback_bundle["source_accessions"])

            organism, genome = _resolve_organism(metadata.get("organism", ""))
            sample_organisms = {
                s.get("organism", "").lower()
                for s in metadata.get("samples", []) if s.get("organism")
            }
            if len(sample_organisms) > 1:
                count_files = files.get("counts", []) or files.get("h5ad", [])
                if count_files:
                    sym_org = _organism_from_gene_symbols(count_files[0])
                    if sym_org:
                        sym_resolved, sym_genome = _resolve_organism(sym_org)
                        if sym_resolved.lower() != organism.lower():
                            _cb(
                                status_cb,
                                f"[GEO] Multiple organisms in metadata "
                                f"({', '.join(sample_organisms)}); gene symbols "
                                f"suggest {sym_resolved} — using that.",
                            )
                            organism, genome = sym_resolved, sym_genome

            if fallback_bundle and organism == "Unknown":
                organism = fallback_bundle["organism"]
                genome = fallback_bundle["genome"]
            design = _infer_design(metadata)
            if fallback_bundle and not design.get("sample_sheet"):
                design = fallback_bundle["inferred_design"]
            design.update({"organism": organism, "genome": genome})
            data_type = (
                fallback_bundle["data_type"] if fallback_bundle
                else _infer_data_type(metadata, files)
            )
            result = {
                "accession": gse_id,
                "source_accessions": source_accessions,
                "title": metadata.get("title", ""),
                "organism": organism,
                "genome": genome,
                "data_type": data_type,
                "data_types": (
                    fallback_bundle.get("data_types", []) if fallback_bundle
                    else [data_type]
                ),
                "n_samples": (
                    len(metadata.get("samples", []))
                    or (fallback_bundle.get("n_samples", 0) if fallback_bundle else 0)
                ),
                "local_dir": str(target),
                "files": _relocate_paths(files, staging, target),
                "file_modalities": _relocate_paths(file_modalities, staging, target),
                "inferred_design": design,
                "geo_metadata": metadata,
                "retrieval_status": "complete",
                "retrieval_manifest": str(target / MANIFEST_NAME),
            }
            return self._publish_result(
                staging, target, result, source_accessions=source_accessions
            )
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise

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
                download_atomic(soft_url, soft_path)
            except Exception as e:
                log.warning(f"SOFT download failed: {e}")
                return {"title": gse_id, "organism": "", "samples": [],
                        "suppl_files": [], "sra_accessions": []}

        try:
            with gzip.open(str(soft_path), "rt", encoding="utf-8",
                           errors="replace") as fh:
                text = fh.read()
        except Exception as e:
            log.warning(f"SOFT read failed: {e}")
            return {"title": gse_id, "organism": "", "samples": [],
                    "suppl_files": [], "sra_accessions": []}

        return _parse_soft_text(text)

    def _download_supplementary(self, gse_id: str, local_dir: Path,
                                  status_cb: Optional[Callable]) -> dict:
        """
        List and download supplementary files from GEO FTP.
        Classifies files into counts, h5ad, h5, mtx (RNA) and fragments, peaks,
        bam (ATAC) buckets.
        """
        files = _empty_file_buckets()
        prefix  = _gse_prefix(gse_id)
        ftp_url = f"{_GEO_FTP}/{prefix}/{gse_id}/suppl/"

        try:
            with open_url_with_retry(ftp_url, timeout=30) as r:
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
            unquote(f) for f in fnames
            if "." in f
            and "/" not in f
            and "://" not in f
            and not f.startswith(("?", ".."))
        ]

        download_records = []
        for fname in sorted(set(fnames)):
            if _SKIP_PATTERNS.search(fname):
                continue
            bucket = _classify_suppl_file(fname)
            archive = _is_archive_name(fname)
            if bucket is None and not archive:
                continue

            fpath = local_dir / "downloads" / fname
            _cb(status_cb, f"[GEO] Downloading {fname}...")
            download_records.append(
                download_atomic(ftp_url + quote(fname), fpath)
            )
            if archive:
                extract_dir = local_dir / "extracted" / _archive_label(fname)
                extracted = safe_extract_archive(fpath, extract_dir)
                for extracted_path in extracted:
                    extracted_bucket = _classify_suppl_file(extracted_path.name)
                    if extracted_bucket is None:
                        continue
                    analyzable = _validated_analyzable_payload(
                        extracted_path, extracted_bucket
                    )
                    files[extracted_bucket].append(str(analyzable))
            elif bucket is not None:
                analyzable = _validated_analyzable_payload(fpath, bucket)
                files[bucket].append(str(analyzable))

        if download_records:
            (local_dir / "geo_downloads.json").write_text(
                json.dumps(download_records, indent=2, sort_keys=True),
                encoding="utf-8",
            )

        if not any(files.values()):
            log.warning(f"No processable supplementary files found for {gse_id}.")

        return files

    # ── SRA path ──────────────────────────────────────────────────────────

    def _fetch_sra(self, sra_id: str,
                   status_cb: Optional[Callable]) -> dict:
        """Fetch RunInfo and materialize every run through the SRA Toolkit."""
        target = self.cache_dir / sra_id
        staging = self._new_staging_dir(sra_id)
        try:
            bundle = self._retrieve_sra_bundle(sra_id, staging, status_cb)
            result = {
                "accession": sra_id,
                "source_accessions": [sra_id],
                "title": bundle["title"],
                "organism": bundle["organism"],
                "genome": bundle["genome"],
                "data_type": bundle["data_type"],
                "data_types": bundle["data_types"],
                "n_samples": bundle["n_samples"],
                "local_dir": str(target),
                "files": _relocate_paths(bundle["files"], staging, target),
                "file_modalities": _relocate_paths(
                    bundle["file_modalities"], staging, target
                ),
                "inferred_design": bundle["inferred_design"],
                "geo_metadata": bundle["geo_metadata"],
                "sra_metadata_csv": str(
                    target / Path(bundle["sra_metadata_csv"]).relative_to(staging)
                ),
                "retrieval_status": "complete",
                "retrieval_manifest": str(target / MANIFEST_NAME),
            }
            return self._publish_result(
                staging, target, result, source_accessions=[sra_id]
            )
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise

    def _fetch_sra_runinfo(self, accession: str) -> list[dict[str, str]]:
        """Resolve an accession with ESearch history, then fetch official RunInfo."""
        try:
            search_url = f"{_NCBI_EUTILS}/esearch.fcgi?" + urlencode({
                "db": "sra",
                "term": accession,
                "retmax": "0",
                "usehistory": "y",
                "retmode": "json",
            })
            with open_url_with_retry(search_url, timeout=120) as response:
                search = json.loads(response.read().decode("utf-8"))["esearchresult"]
            count = int(search.get("count") or 0)
            if count <= 0:
                raise RetrievalError(f"No public SRA records found for {accession}")
            max_runs = int(os.environ.get("ARIA_SRA_MAX_RUNS", "10000"))
            if count > max_runs:
                raise RetrievalError(
                    f"{accession} resolves to {count} SRA records; configured limit "
                    f"ARIA_SRA_MAX_RUNS={max_runs}"
                )
            fetch_url = f"{_NCBI_EUTILS}/efetch.fcgi?" + urlencode({
                "db": "sra",
                "query_key": search["querykey"],
                "WebEnv": search["webenv"],
                "rettype": "runinfo",
                "retmode": "text",
                "retmax": str(count),
            })
            with open_url_with_retry(fetch_url, timeout=120) as response:
                text = response.read().decode("utf-8-sig", errors="strict")
        except RetrievalError:
            raise
        except Exception as exc:
            raise RetrievalError(
                f"SRA RunInfo fetch failed for {accession}: {exc}"
            ) from exc
        rows = [
            {str(key): str(value or "") for key, value in row.items() if key}
            for row in csv.DictReader(io.StringIO(text))
        ]
        rows = [row for row in rows if _sra_run_accession(row)]
        if re.fullmatch(r"[SED]RR\d+", accession.upper()):
            rows = [
                row for row in rows
                if _sra_run_accession(row) == accession.upper()
            ]
        if not rows:
            raise RetrievalError(f"No public SRA runs found for {accession}")
        return rows

    def _retrieve_sra_bundle(
        self,
        accession: str | list[str],
        staging: Path,
        status_cb: Optional[Callable],
    ) -> dict:
        accessions = list(dict.fromkeys(
            [accession] if isinstance(accession, str) else accession
        ))
        rows_by_run: dict[str, dict[str, str]] = {}
        for source_accession in accessions:
            _cb(status_cb, f"[SRA] Fetching RunInfo for {source_accession}...")
            for row in self._fetch_sra_runinfo(source_accession):
                run = _sra_run_accession(row)
                if run:
                    rows_by_run.setdefault(run, row)
        rows = list(rows_by_run.values())
        if not rows:
            raise RetrievalError(
                "No public SRA runs found for " + ", ".join(accessions)
            )
        files = _empty_file_buckets()
        file_modalities: dict[str, str] = {}
        samples = []
        observed_types: list[str] = []

        for row in rows:
            run = _sra_run_accession(row)
            _cb(status_cb, f"[SRA] Retrieving {run} ({len(samples) + 1}/{len(rows)})...")
            outputs = self._retrieve_sra_run(
                row, staging / "sra_runs" / run, status_cb
            )
            if not outputs:
                raise RetrievalError(f"SRA Toolkit produced no FASTQ for {run}")
            modality = _sra_row_modality(row)
            observed_types.append(_sra_row_data_type(row))
            for output in outputs:
                validate_payload(output, "fastq")
                files["fastq"].append(str(output))
                file_modalities[str(output)] = modality
            samples.append(_sra_sample_metadata(row))

        organism_raw = _row_value(rows[0], "ScientificName", "organism_name")
        organism, genome = _resolve_organism(organism_raw)
        library_strategies = [
            _row_value(row, "LibraryStrategy", "library_strategy") for row in rows
        ]
        metadata = {
            "title": accessions[0],
            "organism": organism_raw,
            "samples": samples,
            "suppl_files": [],
            "library_strategy": (
                library_strategies[0]
                if len(set(library_strategies)) == 1 else "mixed"
            ),
            "sra_accessions": accessions,
        }
        design = _infer_design(metadata)
        design.update({"organism": organism, "genome": genome})
        data_types = sorted(set(observed_types))
        data_type = data_types[0] if len(data_types) == 1 else "mixed"

        metadata_path = staging / f"{accessions[0]}_runinfo.csv"
        fieldnames = sorted({key for row in rows for key in row})
        with open(metadata_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        return {
            "title": accessions[0],
            "source_accessions": accessions,
            "organism": organism,
            "genome": genome,
            "data_type": data_type,
            "data_types": data_types,
            "n_samples": len(samples),
            "files": files,
            "file_modalities": file_modalities,
            "inferred_design": design,
            "geo_metadata": metadata,
            "sra_metadata_csv": str(metadata_path),
        }

    def _retrieve_sra_run(
        self,
        row: dict,
        destination: Path,
        status_cb: Optional[Callable] = None,
    ) -> list[Path]:
        """Prefetch, validate, and convert one Run; any failure aborts the batch."""
        run = _sra_run_accession(row)
        tools = {
            name: shutil.which(name)
            for name in ("prefetch", "vdb-validate", "fasterq-dump")
        }
        missing = [name for name, path in tools.items() if not path]
        if missing:
            raise RetrievalError(
                "SRA Toolkit is required for raw retrieval; missing: "
                + ", ".join(missing)
            )
        destination.mkdir(parents=True, exist_ok=True)
        sra_cache = destination / "sra"
        fastq_dir = destination / "fastq"
        sra_cache.mkdir()
        fastq_dir.mkdir()
        timeout = float(os.environ.get("ARIA_SRA_RUN_TIMEOUT", "86400"))
        threads = max(1, int(os.environ.get("ARIA_SRA_THREADS", "8")))

        commands = []
        commands.append(self._run_sra_command(
            [tools["prefetch"], "--output-directory", str(sra_cache), run],
            timeout=timeout,
        ))
        candidates = sorted(sra_cache.rglob(f"{run}.sra"))
        validation_target = candidates[0] if candidates else sra_cache / run
        if not validation_target.exists():
            raise RetrievalError(f"prefetch did not materialize {run}")
        commands.append(self._run_sra_command(
            [tools["vdb-validate"], str(validation_target)], timeout=timeout
        ))
        commands.append(self._run_sra_command(
            [
                tools["fasterq-dump"], "--split-files", "--threads", str(threads),
                "--outdir", str(fastq_dir), str(validation_target),
            ],
            timeout=timeout,
        ))
        outputs = sorted(fastq_dir.glob(f"{run}*.fastq"))
        layout = _row_value(row, "LibraryLayout", "library_layout").upper()
        if layout == "PAIRED":
            names = {path.name for path in outputs}
            if f"{run}_1.fastq" not in names or f"{run}_2.fastq" not in names:
                raise RetrievalError(
                    f"paired SRA run {run} did not produce both mate FASTQs"
                )
        provenance = {
            "run_accession": run,
            "library_layout": layout or "unknown",
            "tools": {
                name: {"path": path, "version": _tool_version(path)}
                for name, path in tools.items()
            },
            "commands": commands,
        }
        (destination / "sra_retrieval.json").write_text(
            json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8"
        )
        return outputs

    @staticmethod
    def _run_sra_command(command: list[str], *, timeout: float) -> dict:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            try:
                os.killpg(process.pid, signal.SIGTERM)
                process.communicate(timeout=5)
            except (ProcessLookupError, subprocess.TimeoutExpired):
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.communicate()
            raise RetrievalError(
                f"SRA command timed out: {Path(command[0]).name}"
            ) from exc
        if process.returncode != 0:
            detail = (stderr or stdout or "unknown failure").strip()[-1000:]
            raise RetrievalError(
                f"SRA command failed ({Path(command[0]).name}): {detail}"
            )
        return {
            "argv": command,
            "returncode": process.returncode,
            "stdout_tail": (stdout or "")[-1000:],
            "stderr_tail": (stderr or "")[-1000:],
        }

    def _new_staging_dir(self, accession: str) -> Path:
        return Path(tempfile.mkdtemp(
            prefix=f".{accession}.", suffix=".staging", dir=str(self.cache_dir)
        ))

    def _publish_result(
        self,
        staging: Path,
        target: Path,
        result: dict,
        *,
        source_accessions: list[str],
    ) -> dict:
        result_path = staging / _RESULT_NAME
        result_path.write_text(
            json.dumps(result, indent=2, sort_keys=True), encoding="utf-8"
        )
        payloads = [
            path for path in staging.rglob("*")
            if path.is_file() and path.name != MANIFEST_NAME
        ]
        write_retrieval_manifest(
            staging,
            accession=result["accession"],
            source_accessions=source_accessions,
            payloads=payloads,
        )
        staged_validation = validate_retrieval_manifest(staging)
        if staged_validation.get("status") != "valid":
            raise RetrievalError(
                "staged retrieval failed manifest validation: "
                + "; ".join(staged_validation.get("errors") or [])
            )
        publish_directory_atomic(staging, target)
        published = self._load_cached_result(result["accession"])
        if published is None:
            raise RetrievalError(
                f"published retrieval cannot be revalidated: {result['accession']}"
            )
        return published

    def _load_cached_result(self, accession: str) -> dict | None:
        target = self.cache_dir / accession
        if not target.is_dir():
            return None
        validation = validate_retrieval_manifest(
            target, expected_accession=accession
        )
        if validation.get("status") != "valid":
            return None
        result_path = target / _RESULT_NAME
        try:
            result = json.loads(result_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        if result.get("retrieval_status") != "complete":
            return None
        if str(result.get("accession") or "").upper() != accession.upper():
            return None
        return result


# ── Module-level helpers ──────────────────────────────────────────────────────

def _empty_file_buckets() -> dict[str, list[str]]:
    return {
        "counts": [], "h5ad": [], "h5": [], "mtx": [], "fragments": [],
        "peaks": [], "bam": [], "fastq": [],
    }


def _has_payload(files: dict) -> bool:
    return any(bool(value) for value in (files or {}).values() if isinstance(value, list))


def _relocate_paths(value, old_root: Path, new_root: Path):
    """Replace staging-root paths recursively, including mapping keys."""
    old_prefix = str(old_root) + os.sep

    def relocate(item):
        if isinstance(item, str) and item.startswith(old_prefix):
            return str(new_root / Path(item).relative_to(old_root))
        return item

    if isinstance(value, dict):
        return {
            relocate(key): _relocate_paths(item, old_root, new_root)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_relocate_paths(item, old_root, new_root) for item in value]
    return relocate(value)


def _is_archive_name(name: str) -> bool:
    lower = name.lower()
    return lower.endswith((".tar", ".tar.gz", ".tgz", ".zip"))


def _archive_label(name: str) -> str:
    label = re.sub(r"(?i)(\.tar\.gz|\.tgz|\.tar|\.zip)$", "", Path(name).name)
    return re.sub(r"[^A-Za-z0-9_.-]", "_", label) or "archive"


def _validated_analyzable_payload(path: Path, bucket: str) -> Path:
    """Validate a supplement and materialize containers downstream cannot read."""
    validate_payload(path, bucket)
    if bucket in {"h5", "h5ad"} and path.name.lower().endswith(".gz"):
        destination = path.with_suffix("")
        decompress_gzip_atomic(path, destination)
        validate_payload(destination, bucket)
        return destination
    return path


def _row_value(row: dict, *keys: str) -> str:
    lowered = {str(key).lower(): value for key, value in row.items()}
    for key in keys:
        value = row.get(key)
        if value in (None, ""):
            value = lowered.get(key.lower())
        if value not in (None, "", "nan"):
            return str(value).strip()
    return ""


def _sra_run_accession(row: dict) -> str:
    run = _row_value(row, "Run", "run_accession")
    return run if re.fullmatch(r"[SED]RR\d+", run) else ""


def _sra_row_modality(row: dict) -> str:
    inferred = _sra_row_data_type(row)
    if inferred == "bulk_RNA":
        return "bulk_RNA_raw"
    return inferred


def _sra_row_data_type(row: dict) -> str:
    strategy = _row_value(row, "LibraryStrategy", "library_strategy")
    title = _row_value(
        row, "SampleName", "sample_title", "LibraryName", "experiment_title"
    )
    metadata = {
        "title": title,
        "samples": [{"title": title}],
        "library_strategy": strategy,
    }
    return _infer_data_type(metadata, {})


def _tool_version(path: str) -> str:
    try:
        result = subprocess.run(
            [path, "--version"], capture_output=True, text=True, timeout=10
        )
        text = (result.stdout or result.stderr or "").strip()
        return text.splitlines()[0][:300] if text else "unknown"
    except Exception:
        return "unknown"


def _sra_sample_metadata(row: dict) -> dict:
    run = _sra_run_accession(row)
    title = _row_value(
        row, "SampleName", "sample_title", "LibraryName", "experiment_title"
    ) or run
    characteristics = {
        "source": _row_value(row, "source", "source_name"),
        "treatment": _row_value(row, "treatment", "condition"),
        "biosample": _row_value(row, "BioSample", "biosample"),
    }
    return {
        "id": run,
        "title": title,
        "organism": _row_value(row, "ScientificName", "organism_name"),
        "characteristics": {
            key: value for key, value in characteristics.items() if value
        },
        "library_strategy": _row_value(
            row, "LibraryStrategy", "library_strategy"
        ),
    }


def _extract_public_sequence_accessions(text: str) -> list[str]:
    return re.findall(
        r"\b(?:SRP|ERP|DRP|PRJNA|PRJEB|PRJDB|SRX|ERX|DRX|SRR|ERR|DRR)\d+\b",
        str(text or "").upper(),
    )


def _preferred_sra_accessions(accessions) -> list[str]:
    """Prefer study/project relations; otherwise union sample/run relations."""
    unique = list(dict.fromkeys(str(value).upper() for value in accessions if value))
    broad = [
        value for value in unique
        if re.fullmatch(r"(?:SRP|ERP|DRP|PRJNA|PRJEB|PRJDB)\d+", value)
    ]
    return broad or unique


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
        "library_strategy": "", "sra_accessions": [],
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
        elif line.startswith("!Series_relation"):
            metadata["sra_accessions"].extend(
                _extract_public_sequence_accessions(val)
            )

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
            elif line.startswith("!Sample_relation"):
                metadata["sra_accessions"].extend(
                    _extract_public_sequence_accessions(val)
                )

    if current:
        metadata["samples"].append(current)

    metadata["sra_accessions"] = list(dict.fromkeys(metadata["sra_accessions"]))
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

    sra_accessions = []
    relations = meta.get("relation", [])
    if not isinstance(relations, list):
        relations = [relations]
    for relation in relations:
        sra_accessions.extend(_extract_public_sequence_accessions(relation))

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
        sample_relations = gm.get("relation", [])
        if not isinstance(sample_relations, list):
            sample_relations = [sample_relations]
        for relation in sample_relations:
            sra_accessions.extend(_extract_public_sequence_accessions(relation))
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
        "sra_accessions": list(dict.fromkeys(sra_accessions)),
    }


def _infer_design(metadata: dict) -> dict:
    """
    Infer experimental groups from sample characteristics.
    Prefer experimental-design keys over incidental sample attributes and use
    GEO/SRA accessions as the canonical sample IDs.  Titles are kept as aliases
    so downstream count matrices can match either GSM/SRR columns or readable
    sample labels.
    """
    from aria.utils.design_matrix import factors_confounded_with_condition

    samples = metadata.get("samples", [])
    if not samples:
        return {"groups": {}, "factor": "condition", "main_factor": "condition",
                "condition_col": "condition", "n_groups": 0,
                "sample_sheet": [], "covariates": [], "donor_col": None,
                "confounded_covariates": [], "unresolved_confounding": False}

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

    condition_col = best_key or "condition"

    # Preserve the full multifactorial structure instead of collapsing the study
    # to a single "best" characteristic.  Every sample keeps all of its
    # characteristics in an explicit sample sheet, secondary multi-level factors
    # are surfaced as covariates, and any secondary factor perfectly aliased with
    # the chosen condition is flagged so the design phase blocks (B2).
    donor_tokens = ("donor", "subject", "patient", "individual", "replicate")
    sample_sheet: list[dict] = []
    donor_col: Optional[str] = None
    for s in samples:
        sample_id = str(s.get("id") or s.get("title") or "sample")
        row = {"sample": sample_id}
        row.update(s.get("characteristics", {}) or {})
        sample_sheet.append(row)

    covariates: list[str] = []
    for key in sorted(all_keys):
        if key == best_key:
            continue
        vals = [
            str(s.get("characteristics", {}).get(key, "")).strip()
            for s in samples
        ]
        vals = [v for v in vals if v]
        unique = len(set(vals))
        if unique < 2:
            continue
        donor_like = any(tok in key.lower() for tok in donor_tokens)
        # A field unique per sample is only design-relevant when it names a
        # donor/subject/replicate block; otherwise it is an incidental id.
        if unique == n and not donor_like:
            continue
        covariates.append(key)
        if donor_like and donor_col is None:
            donor_col = key

    confounded_covariates = (
        factors_confounded_with_condition(sample_sheet, condition_col, covariates)
        if best_key else []
    )

    return {
        "groups":   groups,
        "factor":   condition_col,
        "main_factor": condition_col,
        "condition_col": condition_col,
        "n_groups": len(groups),
        "sample_aliases": sample_aliases,
        "sample_sheet": sample_sheet,
        "covariates": covariates,
        "donor_col": donor_col,
        "confounded_covariates": confounded_covariates,
        "unresolved_confounding": bool(confounded_covariates),
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


def _classify_suppl_file(fname: str) -> str | None:
    """Map a supplementary filename to a file bucket, or None if not processable.

    ATAC artifacts are matched BEFORE the generic mtx/tsv buckets because a
    scATAC fragments file is `*fragments.tsv.gz` and a peak file `*.bed.gz`
    would otherwise be mis-bucketed as a count table.
    """
    fl = fname.lower()
    if fl.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return "fastq"
    if fl.endswith((".h5ad", ".h5ad.gz")):
        return "h5ad"
    if fl.endswith((".h5", ".h5.gz")):
        return "h5"
    if "fragments" in fl and re.search(r"\.(tsv|bed)(\.gz)?$", fl):
        return "fragments"
    if (re.search(r"\.(narrowpeak|broadpeak|gappedpeak)(\.gz)?$", fl)
            or ("peak" in fl and re.search(r"\.bed(\.gz)?$", fl))):
        return "peaks"
    if fl.endswith((".bam", ".cram")):
        return "bam"
    if "matrix.mtx" in fl or fl.endswith(".mtx.gz"):
        return "mtx"
    if re.search(r"\.(txt|tsv|csv)(\.gz)?$", fl):
        return "counts"
    return None


def _infer_data_type(metadata: dict, files: dict) -> str:
    """Infer bulk_RNA, scRNA, bulk_ATAC, scATAC, or unknown from metadata and
    file names. ATAC is resolved BEFORE RNA because single-cell ATAC studies
    also carry the 10x/single-cell keywords that the RNA branch keys on."""
    lib = metadata.get("library_strategy", "").lower()
    title_desc = (metadata.get("title", "") + " " +
                  " ".join(s.get("title", "") for s in metadata.get("samples", []))
                  ).lower()

    single_cell = any(
        kw in title_desc or kw in lib
        for kw in ("scrna", "scatac", "snatac", "single cell", "single-cell",
                   "single nucleus", "single-nucleus", "single-nuclei",
                   "10x", "droplet", "smart-seq", "smartseq", "drop-seq")
    )

    # ── ATAC indicators (checked first) ───────────────────────────────────
    atac = (
        "atac" in lib
        or any(kw in title_desc
               for kw in ("atac-seq", "atac seq", "atacseq", "scatac", "snatac",
                          "chromatin accessibility", "transposase-accessible",
                          "transposase accessible", "transposase",
                          "transposition of native chromatin"))
        or bool(files.get("fragments"))
    )
    if atac:
        # A fragments file is the canonical single-cell ATAC artifact.
        if files.get("fragments") or single_cell:
            return "scATAC"
        return "bulk_ATAC"

    # ── scRNA indicators ──────────────────────────────────────────────────
    if single_cell:
        return "scRNA"
    if files.get("h5ad") or files.get("h5") or files.get("mtx"):
        return "scRNA"

    # ── Bulk RNA indicators ───────────────────────────────────────────────
    if "rna-seq" in lib or "rna seq" in title_desc:
        return "bulk_RNA"
    if files.get("counts"):
        return "bulk_RNA"

    return "unknown"
