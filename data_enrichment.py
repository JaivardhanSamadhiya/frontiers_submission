"""
data_enrichment.py - Compute genomic features for phages and hosts
==================================================================
- Loads FASTA sequences from local directories (with fallback paths).
- Computes di/tri/tetra-nucleotide composition, codon-usage bias, GC content,
  genome length, name-character SVD embeddings, and pair-level features
  (tetra-correlation, codon distance, GC-match, length ratio).
- Reduces tetra (256-dim) and codon (64-dim) to small PCA representations.
- Sanitizes phage names against host-genus tokens for the leakage ablation.

All randomness is local via numpy default_rng(seed).  No global seeds set.
Parallel feature computation uses a strict-ordering helper so that the
output array order matches the input phage/host list order.
"""
from __future__ import annotations

import json
import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import CountVectorizer

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_NUCS = "ACGT"
_VALID_BASES = set("ACGTN")


def _normalize_seq(seq: str) -> str:
    if not isinstance(seq, str):
        return ""
    s = seq.upper().replace("U", "T")
    return "".join(c for c in s if c in _VALID_BASES)


def _read_fasta(path: Path) -> str | None:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None
    parts = []
    for line in text.splitlines():
        if line.startswith(">") or not line.strip():
            continue
        parts.append(line.strip())
    seq = _normalize_seq("".join(parts))
    if len(seq) < 200:
        return None
    return seq


def _find_dirs(cfg: dict, primary_key: str, fallback_key: str) -> list[Path]:
    primary = cfg["paths"][primary_key]
    dirs: list[Path] = []
    if primary.exists() and primary.is_dir():
        dirs.append(primary)
    if "fallback_paths" in cfg:
        for d in cfg["fallback_paths"].get(fallback_key, []):
            try:
                if d.exists() and d.is_dir() and d not in dirs:
                    dirs.append(d)
            except Exception:
                pass
    return dirs


def _slug(name: str) -> str:
    n = str(name).lower().strip()
    n = re.sub(r"[^a-z0-9]+", "_", n)
    return n.strip("_")


def _find_fasta_for(name: str, dirs: list[Path]) -> Path | None:
    if not name:
        return None
    slug = _slug(name)
    for d in dirs:
        # exact
        for ext in (".fasta", ".fa", ".fna", ".fasta.gz"):
            p = d / f"{slug}{ext}"
            if p.exists():
                return p
        # contains
        try:
            for p in d.glob("*.fasta"):
                if slug in _slug(p.stem):
                    return p
            for p in d.glob("*.fa"):
                if slug in _slug(p.stem):
                    return p
        except Exception:
            pass
    return None


def _parallel_compute(fn: Callable, items: list, n_workers: int, desc: str = ""):
    results: list = [None] * len(items)
    with ThreadPoolExecutor(max_workers=max(1, n_workers)) as pool:
        futures = {pool.submit(fn, item): i for i, item in enumerate(items)}
        for future in as_completed(futures):
            i = futures[future]
            try:
                results[i] = future.result()
            except Exception as e:
                log.warning(f"{desc} index {i} failed: {e}")
                results[i] = None
    bad = [i for i, r in enumerate(results) if r is None]
    if bad:
        log.warning(f"{desc}: {len(bad)} item(s) produced None - filling with empty defaults")
        for i in bad:
            results[i] = {}
    return results


# ---------------------------------------------------------------------------
# k-mer composition
# ---------------------------------------------------------------------------
def _kmer_counts(seq: str, k: int) -> np.ndarray:
    if not seq or len(seq) < k:
        return np.zeros(4 ** k, dtype=np.float32)
    # build kmer index
    indices = []
    for i in range(len(seq) - k + 1):
        sub = seq[i:i + k]
        if "N" in sub:
            continue
        idx = 0
        ok = True
        for c in sub:
            v = _NUCS.find(c)
            if v < 0:
                ok = False
                break
            idx = idx * 4 + v
        if ok:
            indices.append(idx)
    if not indices:
        return np.zeros(4 ** k, dtype=np.float32)
    counts = np.bincount(np.array(indices, dtype=np.int64), minlength=4 ** k).astype(np.float32)
    s = counts.sum()
    return counts / s if s > 0 else counts


def _gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for c in seq if c in "GC")
    n = sum(1 for c in seq if c in "ACGT")
    return float(gc) / float(n) if n else 0.0


# codon usage: 64 codons
def _codon_usage(seq: str) -> np.ndarray:
    if not seq or len(seq) < 3:
        return np.zeros(64, dtype=np.float32)
    counts = np.zeros(64, dtype=np.float32)
    for frame in (0, 1, 2):
        s = seq[frame:]
        for i in range(0, len(s) - 2, 3):
            tri = s[i:i + 3]
            if "N" in tri:
                continue
            idx = 0
            ok = True
            for c in tri:
                v = _NUCS.find(c)
                if v < 0:
                    ok = False
                    break
                idx = idx * 4 + v
            if ok:
                counts[idx] += 1
    s = counts.sum()
    return counts / s if s > 0 else counts


# ---------------------------------------------------------------------------
# Per-phage / per-host feature workers
# ---------------------------------------------------------------------------
def _compute_phage_feats(args) -> dict:
    name, seq = args
    if not seq:
        return {
            "name": name,
            "p_di": np.zeros(16, dtype=np.float32),
            "p_tri": np.zeros(64, dtype=np.float32),
            "p_tet": np.zeros(256, dtype=np.float32),
            "p_cub": np.zeros(64, dtype=np.float32),
            "p_gc": 0.0,
            "p_len": 0,
            "has_seq": False,
        }
    return {
        "name": name,
        "p_di": _kmer_counts(seq, 2),
        "p_tri": _kmer_counts(seq, 3),
        "p_tet": _kmer_counts(seq, 4),
        "p_cub": _codon_usage(seq),
        "p_gc": _gc_content(seq),
        "p_len": len(seq),
        "has_seq": True,
    }


def _compute_host_feats(args) -> dict:
    name, seq = args
    if not seq:
        return {
            "name": name,
            "h_di": np.zeros(16, dtype=np.float32),
            "h_tri": np.zeros(64, dtype=np.float32),
            "h_tet": np.zeros(256, dtype=np.float32),
            "h_cub": np.zeros(64, dtype=np.float32),
            "h_gc": 0.0,
            "h_len": 0,
            "has_seq": False,
        }
    return {
        "name": name,
        "h_di": _kmer_counts(seq, 2),
        "h_tri": _kmer_counts(seq, 3),
        "h_tet": _kmer_counts(seq, 4),
        "h_cub": _codon_usage(seq),
        "h_gc": _gc_content(seq),
        "h_len": len(seq),
        "has_seq": True,
    }


# ---------------------------------------------------------------------------
# Sequence loading
# ---------------------------------------------------------------------------
def _load_all_sequences(names: list[str], dirs: list[Path], desc: str) -> list[str | None]:
    seqs: list[str | None] = []
    found = 0
    for n in names:
        p = _find_fasta_for(n, dirs)
        s = _read_fasta(p) if p is not None else None
        if s:
            found += 1
        seqs.append(s)
    log.info(f"[{desc}] sequences found: {found}/{len(names)}")
    return seqs


def _ncbi_backfill(names: list[str], seqs: list[str | None], cfg: dict) -> list[str | None]:
    # No-op stub: real NCBI backfill is intentionally disabled to keep the
    # pipeline deterministic and offline-safe.  The placeholder keeps the
    # phase structure in place for future extension.
    return seqs


# ---------------------------------------------------------------------------
# Name SVD embeddings
# ---------------------------------------------------------------------------
def _svd_embed(strings: list[str], dim: int, max_features: int, seed: int) -> np.ndarray:
    if not strings:
        return np.zeros((0, dim), dtype=np.float32)
    vec = CountVectorizer(analyzer="char_wb", ngram_range=(2, 4),
                          max_features=max_features, lowercase=True)
    try:
        X = vec.fit_transform(strings)
    except ValueError:
        return np.zeros((len(strings), dim), dtype=np.float32)
    n_features = X.shape[1]
    target = min(dim, max(2, n_features - 1), max(2, len(strings) - 1))
    if target < 2:
        return np.zeros((len(strings), dim), dtype=np.float32)
    svd = TruncatedSVD(n_components=target, random_state=seed)
    Z = svd.fit_transform(X).astype(np.float32)
    if Z.shape[1] < dim:
        pad = np.zeros((Z.shape[0], dim - Z.shape[1]), dtype=np.float32)
        Z = np.hstack([Z, pad])
    return Z


def sanitize_phage_name(name: str, host_genera: Iterable[str]) -> str:
    """Replace any host genus appearing in a phage name with HOSTGENUS token."""
    sanitized = str(name)
    for genus in host_genera:
        if not genus:
            continue
        pattern = re.compile(re.escape(genus), re.IGNORECASE)
        sanitized = pattern.sub("HOSTGENUS", sanitized)
    return sanitized


# ---------------------------------------------------------------------------
# Pair-level features
# ---------------------------------------------------------------------------
def _pair_features(phage_feats: dict, host_feats: dict) -> dict:
    """Compute pair-level features and a per-feature computability flag."""
    p_has = phage_feats.get("has_seq", False)
    h_has = host_feats.get("has_seq", False)

    out = {"tetra_corr": 0.0, "cub_dist": 0.0, "gc_match": 0.0, "len_ratio": 0.0,
           "computable": bool(p_has and h_has)}

    if not (p_has and h_has):
        return out

    p_tet = phage_feats["p_tet"]
    h_tet = host_feats["h_tet"]
    try:
        # Pearson correlation, robust to constant vectors
        a = p_tet - p_tet.mean()
        b = h_tet - h_tet.mean()
        denom = (np.linalg.norm(a) * np.linalg.norm(b)) + 1e-12
        out["tetra_corr"] = float(np.dot(a, b) / denom)
    except Exception:
        out["tetra_corr"] = 0.0
    try:
        diff = phage_feats["p_cub"] - host_feats["h_cub"]
        out["cub_dist"] = float(np.linalg.norm(diff))
    except Exception:
        out["cub_dist"] = 0.0
    out["gc_match"] = float(1.0 - abs(phage_feats["p_gc"] - host_feats["h_gc"]))
    plen = phage_feats["p_len"] or 1
    hlen = host_feats["h_len"] or 1
    out["len_ratio"] = float(plen) / float(hlen)
    return out


# ---------------------------------------------------------------------------
# Feature manifest
# ---------------------------------------------------------------------------
def save_feature_manifest(all_feat_cols: list[str], pair_feat_cols: list[str],
                          results_dir: Path) -> dict:
    SEQ_PREFIXES = ("p_di", "p_tri", "p_tet", "p_cub", "h_di", "h_tri",
                    "h_tet", "h_cub", "p_gc", "p_len", "h_gc", "h_len",
                    "tetra", "cub_", "gc_", "len_")
    manifest = {
        "generated_at": datetime.now().isoformat(),
        "total_features": len(all_feat_cols),
        "pair_features": pair_feat_cols,
        "features": [
            {
                "name": c,
                "type": ("sequence-derived"
                         if any(c.startswith(p) for p in SEQ_PREFIXES)
                         else "interaction-derived"),
                "missing_strategy": ("0.0 - sequence unavailable"
                                     if any(c.startswith(p) for p in
                                            ("tetra", "cub_", "gc_", "len_"))
                                     else "q75 imputation"),
                "computable_column": (f"{c}_computable"
                                      if c in pair_feat_cols else None),
            }
            for c in all_feat_cols
        ],
    }
    results_dir.mkdir(parents=True, exist_ok=True)
    out = results_dir / "feature_manifest.json"
    tmp = out.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    tmp.replace(out)
    return manifest


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def compute_all_features(dataset: pd.DataFrame, cfg: dict, seed: int):
    """Compute all features and merge into the dataset.

    Returns
    -------
    dataset_aug : pd.DataFrame  - input + pair-level feature columns + computability flags
    phage_feat_df : pd.DataFrame  - per-phage feature table indexed by phage name
    host_feat_df : pd.DataFrame  - per-host feature table indexed by host name
    feature_meta : dict
    """
    n_workers = int(cfg["features"]["n_feature_workers"])
    n_tet = int(cfg["features"]["tet_pca_dim"])
    n_cub = int(cfg["features"]["cub_pca_dim"])
    svd_dim = int(cfg["features"]["svd_dim"])
    char_max = int(cfg["features"]["char_ngram_max"])

    phage_dirs = _find_dirs(cfg, "phage_fasta_dir", "phage_fasta_dirs")
    host_dirs = _find_dirs(cfg, "host_fasta_dir", "host_fasta_dirs")
    log.info(f"[features] phage FASTA dirs: {[str(p) for p in phage_dirs]}")
    log.info(f"[features] host  FASTA dirs: {[str(p) for p in host_dirs]}")

    phage_names = sorted(dataset["phage"].unique().tolist())
    host_names = sorted(dataset["host"].unique().tolist())

    # Phase 1: load all sequences
    phage_seq_list = _load_all_sequences(phage_names, phage_dirs, "phage")
    host_seq_list = _load_all_sequences(host_names, host_dirs, "host")

    # Phase 2: NCBI backfill (no-op by default)
    phage_seq_list = _ncbi_backfill(phage_names, phage_seq_list, cfg)
    host_seq_list = _ncbi_backfill(host_names, host_seq_list, cfg)

    # Phase 3: verify completion
    n_p = sum(1 for s in phage_seq_list if s and len(s) >= 500)
    n_h = sum(1 for s in host_seq_list if s and len(s) >= 500)
    log.info(f"[features] usable phage seqs: {n_p}/{len(phage_names)}")
    log.info(f"[features] usable host  seqs: {n_h}/{len(host_names)}")

    # Phase 4: parallel feature computation
    ph_results = _parallel_compute(_compute_phage_feats,
                                   list(zip(phage_names, phage_seq_list)),
                                   n_workers, desc="phage_feats")
    ho_results = _parallel_compute(_compute_host_feats,
                                   list(zip(host_names, host_seq_list)),
                                   n_workers, desc="host_feats")

    # Phase 5: stack arrays AFTER all features computed
    p_tet_arr = np.vstack([r.get("p_tet", np.zeros(256, dtype=np.float32))
                           for r in ph_results])
    p_cub_arr = np.vstack([r.get("p_cub", np.zeros(64, dtype=np.float32))
                           for r in ph_results])
    h_tet_arr = np.vstack([r.get("h_tet", np.zeros(256, dtype=np.float32))
                           for r in ho_results])
    h_cub_arr = np.vstack([r.get("h_cub", np.zeros(64, dtype=np.float32))
                           for r in ho_results])
    p_di_arr = np.vstack([r.get("p_di", np.zeros(16, dtype=np.float32))
                          for r in ph_results])
    h_di_arr = np.vstack([r.get("h_di", np.zeros(16, dtype=np.float32))
                          for r in ho_results])

    # Phase 6: PCA after vstack
    def _pca_safe(X: np.ndarray, n: int, label: str) -> np.ndarray:
        if X.shape[0] < 2:
            return np.zeros((X.shape[0], n), dtype=np.float32)
        n_use = min(n, X.shape[0] - 1, X.shape[1])
        if n_use < 1:
            return np.zeros((X.shape[0], n), dtype=np.float32)
        try:
            Z = PCA(n_components=n_use, random_state=seed).fit_transform(X).astype(np.float32)
        except Exception as e:
            log.warning(f"[features] PCA failed for {label}: {e}")
            return np.zeros((X.shape[0], n), dtype=np.float32)
        if Z.shape[1] < n:
            pad = np.zeros((Z.shape[0], n - Z.shape[1]), dtype=np.float32)
            Z = np.hstack([Z, pad])
        return Z

    p_tet_red = _pca_safe(p_tet_arr, n_tet, "p_tet")
    p_cub_red = _pca_safe(p_cub_arr, n_cub, "p_cub")
    h_tet_red = _pca_safe(h_tet_arr, n_tet, "h_tet")
    h_cub_red = _pca_safe(h_cub_arr, n_cub, "h_cub")

    # Build per-phage / per-host feature tables
    phage_records = []
    for i, name in enumerate(phage_names):
        rec = {"phage": name,
               "p_gc": float(ph_results[i].get("p_gc", 0.0)),
               "p_len": float(ph_results[i].get("p_len", 0)),
               "phage_has_seq": bool(ph_results[i].get("has_seq", False))}
        for k in range(p_di_arr.shape[1]):
            rec[f"p_di_{k}"] = float(p_di_arr[i, k])
        for k in range(p_tet_red.shape[1]):
            rec[f"p_tet_{k}"] = float(p_tet_red[i, k])
        for k in range(p_cub_red.shape[1]):
            rec[f"p_cub_{k}"] = float(p_cub_red[i, k])
        phage_records.append(rec)
    phage_feat_df = pd.DataFrame(phage_records).set_index("phage")

    host_records = []
    for i, name in enumerate(host_names):
        rec = {"host": name,
               "h_gc": float(ho_results[i].get("h_gc", 0.0)),
               "h_len": float(ho_results[i].get("h_len", 0)),
               "host_has_seq": bool(ho_results[i].get("has_seq", False))}
        for k in range(h_di_arr.shape[1]):
            rec[f"h_di_{k}"] = float(h_di_arr[i, k])
        for k in range(h_tet_red.shape[1]):
            rec[f"h_tet_{k}"] = float(h_tet_red[i, k])
        for k in range(h_cub_red.shape[1]):
            rec[f"h_cub_{k}"] = float(h_cub_red[i, k])
        host_records.append(rec)
    host_feat_df = pd.DataFrame(host_records).set_index("host")

    # Pair-level features
    log.info("[features] computing pair-level features...")
    pair_cols = ["tetra_corr", "cub_dist", "gc_match", "len_ratio"]
    phage_idx = {n: i for i, n in enumerate(phage_names)}
    host_idx = {n: i for i, n in enumerate(host_names)}
    pair_arr = np.zeros((len(dataset), len(pair_cols)), dtype=np.float32)
    computable = np.zeros(len(dataset), dtype=bool)
    for r_i, (_, r) in enumerate(dataset.iterrows()):
        pi = phage_idx.get(r["phage"])
        hi = host_idx.get(r["host"])
        if pi is None or hi is None:
            continue
        pf = ph_results[pi]
        hf = ho_results[hi]
        feats = _pair_features(pf, hf)
        for c_i, c in enumerate(pair_cols):
            pair_arr[r_i, c_i] = float(feats[c])
        computable[r_i] = bool(feats["computable"])

    dataset_aug = dataset.copy().reset_index(drop=True)
    for c_i, c in enumerate(pair_cols):
        dataset_aug[c] = pair_arr[:, c_i]
    phages_with_seq = {n for i, n in enumerate(phage_names) if ph_results[i].get("has_seq")}
    hosts_with_seq = {n for i, n in enumerate(host_names) if ho_results[i].get("has_seq")}
    for c in pair_cols:
        dataset_aug[f"{c}_computable"] = (
            dataset_aug["phage"].isin(phages_with_seq) &
            dataset_aug["host"].isin(hosts_with_seq)
        )

    # SVD on names (original + sanitized)
    host_genera = sorted({n.split(" ", 1)[0] for n in host_names if n})
    sanitized_names = [sanitize_phage_name(n, host_genera) for n in phage_names]
    log.info(f"[features] computing SVD embeddings (dim={svd_dim})...")
    name_emb = _svd_embed(phage_names, svd_dim, char_max, seed)
    sanit_emb = _svd_embed(sanitized_names, svd_dim, char_max, seed)
    for k in range(name_emb.shape[1]):
        phage_feat_df[f"p_name_{k}"] = name_emb[:, k]
    for k in range(sanit_emb.shape[1]):
        phage_feat_df[f"p_name_sanit_{k}"] = sanit_emb[:, k]

    # Host name SVD (smaller scale)
    host_emb = _svd_embed(host_names, min(32, svd_dim), char_max, seed)
    for k in range(host_emb.shape[1]):
        host_feat_df[f"h_name_{k}"] = host_emb[:, k]

    # Feature metadata
    feature_meta = {
        "pair_features": pair_cols,
        "phage_dim": int(phage_feat_df.shape[1]),
        "host_dim": int(host_feat_df.shape[1]),
        "n_phages_with_seq": int(len(phages_with_seq)),
        "n_hosts_with_seq": int(len(hosts_with_seq)),
        "svd_dim": int(svd_dim),
        "sanitized_name_cols": [c for c in phage_feat_df.columns if c.startswith("p_name_sanit_")],
        "name_cols": [c for c in phage_feat_df.columns if c.startswith("p_name_") and not c.startswith("p_name_sanit_")],
    }

    # Manifest
    pair_feat_cols = pair_cols
    all_feat_cols = list(phage_feat_df.columns) + list(host_feat_df.columns) + pair_feat_cols
    try:
        save_feature_manifest(all_feat_cols, pair_feat_cols, cfg["paths"]["results_dir"])
    except Exception as e:
        log.warning(f"[features] could not save manifest: {e}")

    # Cache features to disk
    try:
        fdir = cfg["paths"]["features_dir"]
        fdir.mkdir(parents=True, exist_ok=True)
        phage_feat_df.to_csv(fdir / "phage_features.csv")
        host_feat_df.to_csv(fdir / "host_features.csv")
    except Exception as e:
        log.warning(f"[features] could not cache features: {e}")

    return dataset_aug, phage_feat_df, host_feat_df, feature_meta
