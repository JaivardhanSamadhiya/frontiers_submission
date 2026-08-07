"""
data_collection.py - Assemble positive phage-host interaction pairs
=====================================================================
Sources, in order:
  1. VirusHostInter.csv (local)
  2. Virus-Host DB (optional network)
  3. INPHARED TSV (local, if present)
  4. NCBI nuccore (optional network)
  5. PhagesDB cache (local JSON, optionally refreshed)

This module NEVER constructs negatives.  All returned pairs are positives
(infection observed in the wet lab).  The build_labels module is responsible
for constructing negatives.

The held-out S. aureus experimental validation files are loaded by
load_experimental_validation() ONLY - they are never used during training,
feature computation, hyperparameter search, or negative sampling.
"""
from __future__ import annotations

import hashlib
import json
import logging
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# String helpers
# ---------------------------------------------------------------------------
_STRIP_PATTERNS = [
    r"\s+strain\b.*$",
    r"\s+subsp\.?\b.*$",
    r"\s+ssp\.?\b.*$",
    r"\s+serovar\b.*$",
    r"\s+sv\b.*$",
    r"\s+pv\.?\b.*$",
    r"\s+biovar\b.*$",
    r"\s+bv\b.*$",
    r"\s+sp\.?\s*$",
    r"\s+\(.*\)\s*$",
    r"\bphage\b",
    r"\bvirus\b",
    r"\bbacteriophage\b",
]


def _clean_host(name: str) -> str:
    if not isinstance(name, str):
        return ""
    n = name.strip().lower()
    n = re.sub(r"_", " ", n)
    for pat in _STRIP_PATTERNS:
        n = re.sub(pat, "", n, flags=re.IGNORECASE)
    n = re.sub(r"[^a-z0-9 ]+", " ", n)
    n = re.sub(r"\s+", " ", n).strip()
    return n


def _clean_phage(name: str) -> str:
    if not isinstance(name, str):
        return ""
    n = name.strip().lower()
    n = re.sub(r"_", " ", n)
    n = re.sub(r"\s+", " ", n).strip()
    return n


def _genus_of(host: str) -> str:
    h = _clean_host(host)
    return h.split(" ", 1)[0] if h else ""


# ---------------------------------------------------------------------------
# PhagesDB cache helpers
# ---------------------------------------------------------------------------
def _find_phagesdb_cache(cfg: dict) -> Path | None:
    candidates: list[Path] = [cfg["paths"]["phagesdb_cache"]]
    if "fallback_paths" in cfg:
        candidates.extend(cfg["fallback_paths"].get("phagesdb_cache_candidates", []))
    for c in candidates:
        try:
            if c.exists():
                return c
        except Exception:
            pass
    return None


def _load_phagesdb_cache(cache_path: Path) -> list:
    if cache_path is None:
        return []
    try:
        if cache_path.exists():
            return json.loads(cache_path.read_text(encoding="utf-8"))
    except Exception as e:
        log.warning(f"[phagesdb] failed to load cache {cache_path}: {e}")
    return []


def _save_phagesdb_cache(records: list, cache_path: Path) -> None:
    tmp = cache_path.with_suffix(cache_path.suffix + ".tmp")
    tmp.write_text(json.dumps(records, indent=2), encoding="utf-8")
    tmp.replace(cache_path)


# ---------------------------------------------------------------------------
# Source 1: VirusHostInter.csv
# ---------------------------------------------------------------------------
def _from_vhi(cfg: dict) -> pd.DataFrame:
    p = cfg["paths"]["vhi_csv"]
    if not p.exists():
        log.warning(f"[vhi] file not found: {p}")
        return pd.DataFrame(columns=["phage", "host", "source"])
    log.info(f"[vhi] reading {p}")
    df = pd.read_csv(p)
    # detect columns
    host_col = next((c for c in df.columns if c.lower() in
                     ("hostname", "host", "host_name", "bacterium")), None)
    phage_col = next((c for c in df.columns if c.lower() in
                      ("phagename", "phage", "phage_name", "virus")), None)
    inf_col = next((c for c in df.columns if c.lower() in
                    ("infection", "label", "outcome", "interacts")), None)
    if host_col is None or phage_col is None:
        log.error(f"[vhi] could not detect host/phage columns: {df.columns.tolist()}")
        return pd.DataFrame(columns=["phage", "host", "source"])
    out = pd.DataFrame({
        "phage": df[phage_col].astype(str).map(_clean_phage),
        "host":  df[host_col].astype(str).map(_clean_host),
    })
    out["source"] = "VHI"
    if inf_col is not None:
        mask = df[inf_col].astype(str).str.lower().isin(
            ["inf", "infection", "1", "true", "infects", "yes", "positive"])
        out = out.loc[mask].reset_index(drop=True)
    out = out[(out["phage"] != "") & (out["host"] != "")]
    log.info(f"[vhi] kept {len(out)} positive pairs")
    return out.drop_duplicates(["phage", "host"])


# ---------------------------------------------------------------------------
# Source 2: Virus-Host DB (local cache from fetch_data.fetch_virushostdb)
# ---------------------------------------------------------------------------
def _pick_column(columns: list[str], exact: list[str], contains: list[str]) -> str | None:
    lower = {c: c.lower() for c in columns}
    for want in exact:
        for col, lc in lower.items():
            if lc == want:
                return col
    for needle in contains:
        for col, lc in lower.items():
            if needle in lc:
                return col
    return None


def _from_virushostdb(cfg: dict, skip_download: bool) -> pd.DataFrame:
    local = cfg["paths"]["raw_dir"] / "virushostdb_raw.tsv"
    if not local.exists():
        if skip_download:
            log.info("[vhdb] virushostdb_raw.tsv missing (--skip-download)")
        else:
            log.warning("[vhdb] virushostdb_raw.tsv missing; run without --skip-download")
        return pd.DataFrame(columns=["phage", "host", "source"])
    log.info(f"[vhdb] reading {local}")
    try:
        df = pd.read_csv(local, sep="\t", on_bad_lines="skip", low_memory=False)
    except Exception as e:
        log.warning(f"[vhdb] failed to read: {e}")
        return pd.DataFrame(columns=["phage", "host", "source"])

    host_col = _pick_column(df.columns.tolist(),
                            exact=["host name"],
                            contains=["host name", "host lineage"])
    phage_col = _pick_column(df.columns.tolist(),
                             exact=["virus name"],
                             contains=["virus name"])
    if host_col is None or phage_col is None:
        log.warning(f"[vhdb] columns not detected; got {df.columns.tolist()[:8]}")
        return pd.DataFrame(columns=["phage", "host", "source"])
    out = pd.DataFrame({
        "phage": df[phage_col].astype(str).map(_clean_phage),
        "host":  df[host_col].astype(str).map(
            lambda x: _clean_host(x.split(";")[-1] if isinstance(x, str) else "")),
    })
    out["source"] = "Virus-Host-DB"
    out = out[(out["phage"] != "") & (out["host"] != "")]
    out = out[out["phage"].str.contains("phage|virus|bacteriop", case=False, na=False)]
    log.info(f"[vhdb] kept {len(out)} pairs")
    return out.drop_duplicates(["phage", "host"])


# ---------------------------------------------------------------------------
# Source 3: INPHARED
# ---------------------------------------------------------------------------
def _from_inphared(cfg: dict) -> pd.DataFrame:
    candidates = [cfg["paths"]["inphared_dir"]]
    if "fallback_paths" in cfg:
        candidates.extend(cfg["fallback_paths"].get("inphared_dirs", []))
    base = None
    for c in candidates:
        try:
            if c.exists() and c.is_dir():
                base = c
                break
        except Exception:
            pass
    if base is None:
        log.info("[inphared] no INPHARED directory found - skipping")
        return pd.DataFrame(columns=["phage", "host", "source"])

    tsvs = list(base.glob("*data_excluding*.tsv")) + list(base.glob("*INPHARED*.tsv"))
    if not tsvs:
        tsvs = list(base.glob("*.tsv"))
    tsvs = [f for f in tsvs if f.stat().st_size > 500
            and not f.read_bytes()[:20].startswith((b"<?xml", b"<!DOCTYPE", b"<html"))]
    if not tsvs:
        log.info(f"[inphared] no valid TSV files in {base} (optional source)")
        return pd.DataFrame(columns=["phage", "host", "source"])

    frames = []
    for f in tsvs[:2]:
        try:
            df = pd.read_csv(f, sep="\t", on_bad_lines="skip", low_memory=False)
            if len(df) < 10 or len(df.columns) < 3:
                log.warning(f"[inphared] skipping {f.name}: too few rows/columns")
                continue
            log.info(f"[inphared] {f.name}: {len(df)} rows")
        except Exception as e:
            log.warning(f"[inphared] failed to read {f}: {e}")
            continue
        host_col = next((c for c in df.columns if "host" in c.lower()), None)
        name_col = next((c for c in df.columns if c.lower() in
                         ("description", "name", "phage", "accession")), None)
        if host_col is None or name_col is None:
            continue
        sub = pd.DataFrame({
            "phage": df[name_col].astype(str).map(_clean_phage),
            "host":  df[host_col].astype(str).map(_clean_host),
        })
        sub["source"] = "INPHARED"
        frames.append(sub)
    if not frames:
        return pd.DataFrame(columns=["phage", "host", "source"])
    out = pd.concat(frames, ignore_index=True)
    out = out[(out["phage"] != "") & (out["host"] != "") & (out["host"] != "unspecified")]
    log.info(f"[inphared] kept {len(out)} pairs")
    return out.drop_duplicates(["phage", "host"])


# ---------------------------------------------------------------------------
# Source 4: NCBI nuccore (lightweight - phage_ncbi_refseq_def_info.txt if present)
# ---------------------------------------------------------------------------
def _from_ncbi(cfg: dict, skip_download: bool) -> pd.DataFrame:
    local = cfg["paths"]["raw_dir"] / "phage_ncbi_refseq_def_info.txt"
    if not local.exists():
        log.info("[ncbi] no local refseq def info; skipping (cannot safely query "
                 "NCBI from pipeline without rate-limit cooperation)")
        return pd.DataFrame(columns=["phage", "host", "source"])
    rows = []
    try:
        text = local.read_text(encoding="utf-8", errors="ignore")
    except Exception as e:
        log.warning(f"[ncbi] could not read {local}: {e}")
        return pd.DataFrame(columns=["phage", "host", "source"])
    for line in text.splitlines():
        if not line.strip():
            continue
        m_host = re.search(r"host[=:\s]+([A-Z][a-z]+\s[a-z]+)", line)
        m_name = re.search(r"\b([A-Za-z0-9_\-]+)\s+phage\s+([A-Za-z0-9_\-\.]+)", line)
        if m_host and m_name:
            rows.append({"phage": _clean_phage(m_name.group(2)),
                         "host": _clean_host(m_host.group(1)),
                         "source": "NCBI"})
    out = pd.DataFrame(rows)
    if len(out) == 0:
        return pd.DataFrame(columns=["phage", "host", "source"])
    out = out[(out["phage"] != "") & (out["host"] != "")]
    log.info(f"[ncbi] kept {len(out)} pairs")
    return out.drop_duplicates(["phage", "host"])


# ---------------------------------------------------------------------------
# Source 5: PhagesDB cache (no network unless missing)
# ---------------------------------------------------------------------------
def _from_phagesdb(cfg: dict, skip_download: bool) -> pd.DataFrame:
    cache_path = _find_phagesdb_cache(cfg)
    cached = _load_phagesdb_cache(cache_path) if cache_path is not None else []
    if not cached:
        if skip_download:
            log.info("[phagesdb] cache missing (--skip-download)")
        else:
            log.warning("[phagesdb] cache missing; run without --skip-download to fetch from API")
        return pd.DataFrame(columns=["phage", "host", "source"])
    log.info(f"[phagesdb] loaded {len(cached)} cached records (cache: {cache_path})")

    rows = []
    for rec in cached:
        if not isinstance(rec, dict):
            continue
        name = rec.get("name") or rec.get("name_original") or ""
        host_g = rec.get("host_genus", "") or ""
        host_s = rec.get("host_species", "") or ""
        host = f"{host_g} {host_s}".strip()
        if not name or not host:
            continue
        rows.append({"phage": _clean_phage(str(name)),
                     "host":  _clean_host(host),
                     "source": "PhagesDB"})
    out = pd.DataFrame(rows)
    if len(out) == 0:
        return out
    out = out[(out["phage"] != "") & (out["host"] != "")]
    log.info(f"[phagesdb] kept {len(out)} pairs")
    # Atomic re-save if path was found (idempotent)
    if cache_path is not None:
        try:
            _save_phagesdb_cache(cached, cache_path)
        except Exception as e:
            log.warning(f"[phagesdb] atomic re-save failed: {e}")
    return out.drop_duplicates(["phage", "host"])


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def assemble_dataset(cfg: dict, skip_download: bool = False) -> pd.DataFrame:
    """Return positive (phage, host) pairs only.  Never generates negatives."""
    frames: list[pd.DataFrame] = []
    for src_fn, name in [
        (_from_vhi, "VHI"),
        (_from_virushostdb, "Virus-Host-DB"),
        (_from_inphared, "INPHARED"),
        (_from_ncbi, "NCBI"),
        (_from_phagesdb, "PhagesDB"),
    ]:
        try:
            if name in ("Virus-Host-DB",):
                df = src_fn(cfg, skip_download)
            elif name in ("NCBI", "PhagesDB"):
                df = src_fn(cfg, skip_download)
            else:
                df = src_fn(cfg)
            frames.append(df)
        except Exception as e:
            log.exception(f"[{name}] failed: {e}")
            frames.append(pd.DataFrame(columns=["phage", "host", "source"]))

    if not any(len(f) for f in frames):
        log.error("All data sources empty - cannot continue")
        raise RuntimeError("assemble_dataset produced 0 positive pairs")

    pos = pd.concat(frames, ignore_index=True)
    pos = pos[(pos["phage"] != "") & (pos["host"] != "")]
    pos["host"] = pos["host"].map(_clean_host)
    pos["phage"] = pos["phage"].map(_clean_phage)
    pos["host_genus"] = pos["host"].map(_genus_of)
    pos["label"] = 1
    pos["is_constructed_negative"] = False

    # de-dup preferring earlier sources (VHI first)
    pos = pos.drop_duplicates(["phage", "host"]).reset_index(drop=True)

    log.info(f"[assemble] total unique positive pairs: {len(pos)}")
    by_source = pos.groupby("source").size().to_dict()
    for src, n in by_source.items():
        log.info(f"  {src}: {n}")

    # Save a snapshot in raw/
    snap = cfg["paths"]["results_dir"] / "positives_snapshot.csv"
    pos.to_csv(snap, index=False)
    log.info(f"[assemble] positives snapshot saved to {snap}")
    return pos


def load_experimental_validation(cfg: dict) -> pd.DataFrame | None:
    """Load held-out S. aureus experimental validation pairs.

    These are NEVER used for training.  Returns a DataFrame with columns:
        phage, host, label  (0 = no infection, 1 = infection)
    Returns None if neither file can be parsed.
    """
    exp_path = cfg["paths"]["saureus_experimental"]
    sec_path = cfg["paths"]["saureus_phage"]
    rows: list[dict] = []
    n_pos_e = n_neg_e = 0
    n_pos_s = n_neg_s = 0

    # File 1: saureus_experimental_interactions.csv (VHI-like schema)
    if exp_path.exists():
        try:
            df = pd.read_csv(exp_path)
            host_col = next((c for c in df.columns if c.lower() in
                             ("hostname", "host", "host_name")), None)
            phage_col = next((c for c in df.columns if c.lower() in
                              ("phagename", "phage", "phage_name")), None)
            inf_col = next((c for c in df.columns if c.lower() in
                            ("infection", "label", "outcome")), None)
            if host_col is not None and phage_col is not None and inf_col is not None:
                for _, r in df.iterrows():
                    lab_raw = str(r[inf_col]).lower()
                    label = 1 if lab_raw in ("inf", "1", "true", "yes", "infects",
                                             "infection", "positive") else 0
                    rows.append({
                        "phage": _clean_phage(str(r[phage_col])),
                        "host":  _clean_host(str(r[host_col])),
                        "label": int(label),
                        "source": "saureus_experimental",
                    })
                    if label == 1:
                        n_pos_e += 1
                    else:
                        n_neg_e += 1
                log.info(f"  saureus_experimental: {n_pos_e} pos, {n_neg_e} neg "
                         f"({exp_path.name})")
            else:
                log.warning(f"  could not detect columns in {exp_path.name}")
        except Exception as e:
            log.warning(f"  failed reading {exp_path.name}: {e}")
    else:
        log.warning(f"  {exp_path} not found")

    # File 2: saureus_phage_interactions.csv (phage_id / host / flag schema)
    if sec_path.exists():
        try:
            df2 = pd.read_csv(sec_path, low_memory=False)
            host_col = next((c for c in df2.columns if c.lower() in
                             ("host", "host_species", "host_name")), None)
            phage_col = next((c for c in df2.columns if c.lower() in
                              ("phage_id", "phagename", "phage", "phage_def")), None)
            flag_col = next((c for c in df2.columns if c.lower() in
                             ("flag", "label", "split", "type")), None)
            if host_col is not None and phage_col is not None:
                for _, r in df2.iterrows():
                    flag = str(r.get(flag_col, "")).lower() if flag_col else ""
                    # treat any annotated row in this file as a positive
                    # unless its flag explicitly encodes a negative
                    if "neg" in flag or "no_inf" in flag:
                        label = 0
                    else:
                        label = 1
                    rows.append({
                        "phage": _clean_phage(str(r[phage_col])),
                        "host":  _clean_host(str(r[host_col])),
                        "label": int(label),
                        "source": "saureus_phage_curated",
                    })
                    if label == 1:
                        n_pos_s += 1
                    else:
                        n_neg_s += 1
                log.info(f"  saureus_phage_curated: {n_pos_s} pos, {n_neg_s} neg "
                         f"({sec_path.name})")
            else:
                log.warning(f"  could not detect columns in {sec_path.name}")
        except Exception as e:
            log.warning(f"  failed reading {sec_path.name}: {e}")
    else:
        log.warning(f"  {sec_path} not found")

    if not rows:
        return None
    out = pd.DataFrame(rows)
    out = out[(out["phage"] != "") & (out["host"] != "")]
    out = out.drop_duplicates(["phage", "host"]).reset_index(drop=True)
    out["host_genus"] = out["host"].map(_genus_of)
    return out


def freeze_dataset(df: pd.DataFrame, version: str, results_dir: Path) -> str:
    """Atomically write the frozen dataset and its SHA-256 hash."""
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / f"frozen_dataset_v{version}.csv"
    tmp = out_path.with_suffix(".csv.tmp")
    df.to_csv(tmp, index=False)
    tmp.replace(out_path)
    sha256 = hashlib.sha256(out_path.read_bytes()).hexdigest()
    hash_path = results_dir / f"frozen_dataset_v{version}.sha256"
    hash_path.write_text(sha256, encoding="utf-8")
    log.info(f"Dataset frozen -> {out_path.name}")
    log.info(f"SHA-256: {sha256}")
    log.info("Archive this file on Zenodo/Figshare before submission.")
    return sha256
