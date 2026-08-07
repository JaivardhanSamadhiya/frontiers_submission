"""
fetch_data.py - Input verification and online data acquisition
==============================================================
Step 0 of run_pipeline.py (skipped with --skip-download for network only).

MANUAL (must exist locally — verified, never downloaded):
  - data/raw/VirusHostInter.csv
  - data/raw/phage_ncbi_refseq_def_info.txt
  - data/fastas/hosts/*.fasta

ONLINE (dedicated fetch functions below):
  - fetch_virushostdb()        → data/raw/virushostdb_raw.tsv
  - fetch_phagesdb_cache()     → data/raw/phagesdb_records_cache.json
  - fetch_phage_fastas()       → data/phages/*.fasta  (NCBI Entrez)

DERIVED from manual VHI (no network):
  - ensure_validation_csvs()   → saureus_*_interactions.csv

OPTIONAL manual:
  - data/raw/inphared/*.tsv    (INPHARED; skipped if absent)
"""
from __future__ import annotations

import gzip
import json
import logging
import re
import shutil
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd

log = logging.getLogger(__name__)

PHAGESDB_API = "https://phagesdb.org/api/phages/"
VHDB_URL = "https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv"
USER_AGENT = "PrecisionPhage-Pipeline/1.0 (research)"


class MissingInputError(FileNotFoundError):
    """Raised when a required manual input file is absent."""


# ---------------------------------------------------------------------------
# Network helper
# ---------------------------------------------------------------------------
def _fetch(url: str, timeout: int = 120, retries: int = 3) -> bytes:
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (URLError, HTTPError, TimeoutError, OSError) as e:
            last_err = e
            wait = min(2 ** attempt, 30)
            log.warning("[fetch] %s attempt %d/%d failed (%s); sleeping %ds",
                        url, attempt, retries, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"Failed to fetch {url}: {last_err}")


def _clean_phage(name: str) -> str:
    if not isinstance(name, str):
        return ""
    return re.sub(r"\s+", " ", name.strip().lower()).strip()


def _clean_host(name: str) -> str:
    if not isinstance(name, str):
        return ""
    n = name.strip().lower().replace("_", " ")
    n = re.sub(r"[^a-z0-9 ]+", " ", n)
    return re.sub(r"\s+", " ", n).strip()


# ---------------------------------------------------------------------------
# Manual input verification
# ---------------------------------------------------------------------------
def verify_manual_inputs(cfg: dict) -> None:
    """Fail fast if required local files are missing."""
    missing: list[str] = []

    vhi = Path(cfg["paths"]["vhi_csv"])
    if not vhi.exists():
        missing.append(str(vhi))

    refseq_meta = Path(cfg["paths"]["raw_dir"]) / "phage_ncbi_refseq_def_info.txt"
    if not refseq_meta.exists():
        missing.append(str(refseq_meta))

    host_dir = Path(cfg["paths"]["host_fasta_dir"])
    host_fastas = list(host_dir.glob("*.fasta")) if host_dir.exists() else []
    if len(host_fastas) < 10:
        missing.append(f"{host_dir} (need host *.fasta files, found {len(host_fastas)})")

    if missing:
        raise MissingInputError(
            "Required manual inputs missing:\n  - "
            + "\n  - ".join(missing)
            + "\nSee DATA.md for what to place locally."
        )

    log.info("[verify] manual inputs OK:")
    log.info("[verify]   %s", vhi.name)
    log.info("[verify]   %s", refseq_meta.name)
    log.info("[verify]   %d host FASTAs in %s", len(host_fastas), host_dir)


# ---------------------------------------------------------------------------
# Online source 1: Virus-Host DB (KEGG FTP)
# ---------------------------------------------------------------------------
def fetch_virushostdb(raw_dir: Path) -> Path:
    """Download Virus-Host DB TSV if missing."""
    out = raw_dir / "virushostdb_raw.tsv"
    if out.exists() and out.stat().st_size > 10_000:
        log.info("[fetch:virushostdb] using cached %s (%d bytes)",
                 out.name, out.stat().st_size)
        return out

    log.info("[fetch:virushostdb] downloading from KEGG FTP...")
    raw_dir.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tsv.tmp")
    tmp.write_bytes(_fetch(VHDB_URL))
    tmp.replace(out)
    log.info("[fetch:virushostdb] saved %s (%d bytes)", out.name, out.stat().st_size)
    return out


# ---------------------------------------------------------------------------
# Online source 2: PhagesDB REST API
# ---------------------------------------------------------------------------
def fetch_phagesdb_cache(cache_path: Path, min_records: int = 1000) -> Path:
    """Paginate PhagesDB REST API into a local JSON cache."""
    if cache_path.exists():
        try:
            records = json.loads(cache_path.read_text(encoding="utf-8"))
            if isinstance(records, list) and len(records) >= min_records:
                log.info("[fetch:phagesdb] using cache (%d records)", len(records))
                return cache_path
        except Exception:
            pass

    log.info("[fetch:phagesdb] downloading records (paginated API)...")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    all_records: list[dict] = []
    url: str | None = PHAGESDB_API
    page = 0
    while url:
        page += 1
        payload = json.loads(_fetch(url, timeout=180))
        batch = payload.get("results", [])
        if not isinstance(batch, list):
            break
        all_records.extend(batch)
        url = payload.get("next")
        if page % 10 == 0:
            log.info("[fetch:phagesdb]   page %d (%d records)", page, len(all_records))
        time.sleep(0.15)

    tmp = cache_path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(all_records, indent=2), encoding="utf-8")
    tmp.replace(cache_path)
    log.info("[fetch:phagesdb] saved %d records -> %s", len(all_records), cache_path.name)
    return cache_path


# ---------------------------------------------------------------------------
# Online source 3: NCBI RefSeq phage FASTAs (Entrez)
# ---------------------------------------------------------------------------
def _parse_refseq_accessions(refseq_file: Path) -> list[str]:
    accessions: set[str] = set()
    with refseq_file.open(encoding="utf-8", errors="ignore") as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            acc = parts[1].strip()
            desc = parts[2].lower() if len(parts) > 2 else ""
            if "phage" not in desc:
                continue
            if acc.startswith(("NC_", "NZ_", "CP_", "NW_")):
                accessions.add(acc)
    return sorted(accessions)


def fetch_phage_fastas(cfg: dict) -> int:
    """Download bacteriophage RefSeq FASTAs listed in phage_ncbi_refseq_def_info.txt."""
    try:
        from Bio import Entrez
    except ImportError as e:
        raise ImportError(
            "biopython is required for NCBI phage FASTA download. "
            "Install with: pip install biopython"
        ) from e

    raw_dir = Path(cfg["paths"]["raw_dir"])
    phage_dir = Path(cfg["paths"]["phage_fasta_dir"])
    phage_dir.mkdir(parents=True, exist_ok=True)

    refseq_file = raw_dir / "phage_ncbi_refseq_def_info.txt"
    accessions = _parse_refseq_accessions(refseq_file)
    if not accessions:
        log.warning("[fetch:ncbi-phages] no accessions in %s", refseq_file.name)
        return 0

    existing = sum(1 for acc in accessions
                   if (phage_dir / f"{acc}.fasta").exists()
                   and (phage_dir / f"{acc}.fasta").stat().st_size > 200)
    log.info("[fetch:ncbi-phages] %d accessions (%d already on disk)",
             len(accessions), existing)

    Entrez.email = cfg.get("fetch", {}).get("ncbi_email", "precisionphage@example.com")
    n_new = 0
    for i, acc in enumerate(accessions):
        out = phage_dir / f"{acc}.fasta"
        if out.exists() and out.stat().st_size > 200:
            continue
        try:
            with Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text") as handle:
                seq = handle.read()
            if seq and len(seq) > 100:
                out.write_text(seq, encoding="utf-8")
                n_new += 1
        except Exception as e:
            log.warning("[fetch:ncbi-phages] efetch failed for %s: %s", acc, e)
        if (i + 1) % 100 == 0:
            log.info("[fetch:ncbi-phages]   %d/%d processed (%d new)",
                     i + 1, len(accessions), n_new)
        time.sleep(0.34)

    total = sum(1 for acc in accessions
                if (phage_dir / f"{acc}.fasta").exists())
    log.info("[fetch:ncbi-phages] done: %d new, %d/%d total in %s",
             n_new, total, len(accessions), phage_dir)
    return n_new


# ---------------------------------------------------------------------------
# Derived: S. aureus validation CSVs from VirusHostInter
# ---------------------------------------------------------------------------
def ensure_validation_csvs(cfg: dict) -> None:
    """Build held-out S. aureus validation CSVs from VirusHostInter if missing."""
    vhi_path = Path(cfg["paths"]["vhi_csv"])
    exp_out = Path(cfg["paths"]["saureus_experimental"])
    phage_out = Path(cfg["paths"]["saureus_phage"])

    if exp_out.exists() and phage_out.exists():
        log.info("[derive:saureus] validation CSVs already present")
        return

    df = pd.read_csv(vhi_path)
    host_col = next((c for c in df.columns if c.lower() in
                     ("hostname", "host", "host_name", "bacterium")), None)
    phage_col = next((c for c in df.columns if c.lower() in
                      ("phagename", "phage", "phage_name", "virus")), None)
    inf_col = next((c for c in df.columns if c.lower() in
                    ("infection", "label", "outcome", "interacts")), None)
    if host_col is None or phage_col is None or inf_col is None:
        log.warning("[derive:saureus] VirusHostInter columns not recognized; skipping")
        return

    mask = (df[host_col].astype(str).str.lower()
            .str.replace("_", " ", regex=False)
            .str.contains("staphylococcus aureus", na=False))
    sa = df.loc[mask].copy()
    if len(sa) == 0:
        log.warning("[derive:saureus] no S. aureus rows in VirusHostInter")
        return

    exp_out.parent.mkdir(parents=True, exist_ok=True)
    if not exp_out.exists():
        sa.to_csv(exp_out, index=False)
        log.info("[derive:saureus] wrote %s (%d rows)", exp_out.name, len(sa))

    if not phage_out.exists():
        pos = sa[sa[inf_col].astype(str).str.lower().isin(
            ["inf", "infection", "1", "true", "infects", "yes", "positive"])]
        rows = [{
            "phage_id": _clean_phage(str(r[phage_col])),
            "host": _clean_host(str(r[host_col])),
            "flag": "curated_positive",
        } for _, r in pos.iterrows()]
        pd.DataFrame(rows).drop_duplicates(["phage_id", "host"]).to_csv(phage_out, index=False)
        log.info("[derive:saureus] wrote %s (%d rows)", phage_out.name, len(rows))


# ---------------------------------------------------------------------------
# Optional online: INPHARED metadata TSV
# ---------------------------------------------------------------------------
def fetch_inphared_tsv(cfg: dict) -> Path | None:
    """Download optional INPHARED TSV if URL configured and file absent."""
    url = cfg.get("fetch", {}).get("inphared_tsv_url")
    if not url:
        return None

    inphared_dir = Path(cfg["paths"]["inphared_dir"])
    inphared_dir.mkdir(parents=True, exist_ok=True)
    existing = [f for f in inphared_dir.glob("*.tsv")
                if f.stat().st_size > 500
                and not f.read_bytes()[:20].startswith((b"<?xml", b"<!DOCTYPE"))]
    if existing:
        log.info("[fetch:inphared] using %s", existing[0].name)
        return existing[0]

    log.info("[fetch:inphared] downloading from configured URL...")
    raw = _fetch(url, timeout=300)
    out = inphared_dir / Path(url.split("?")[0]).name
    if out.suffix == ".gz" or raw[:2] == b"\x1f\x8b":
        tmp_gz = out if out.suffix == ".gz" else out.with_suffix(out.suffix + ".gz")
        tmp_gz.write_bytes(raw)
        out = tmp_gz.with_suffix("")
        with gzip.open(tmp_gz, "rb") as fin, out.open("wb") as fout:
            shutil.copyfileobj(fin, fout)
        tmp_gz.unlink(missing_ok=True)
    else:
        out.write_bytes(raw)
    log.info("[fetch:inphared] saved %s (%d bytes)", out.name, out.stat().st_size)
    return out


# ---------------------------------------------------------------------------
# Orchestrator (called by run_pipeline.py step 0)
# ---------------------------------------------------------------------------
def fetch_all(cfg: dict) -> None:
    """Fetch all online sources (manual inputs must already be verified)."""
    raw_dir = Path(cfg["paths"]["raw_dir"])
    cache_path = Path(cfg["paths"]["phagesdb_cache"])

    fetch_virushostdb(raw_dir)
    fetch_phagesdb_cache(cache_path)
    fetch_inphared_tsv(cfg)
    ensure_validation_csvs(cfg)
    fetch_phage_fastas(cfg)
