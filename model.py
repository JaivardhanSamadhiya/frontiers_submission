"""
model.py - Feature matrix assembly and device helpers
=====================================================
Stores RAW (unscaled) feature matrices in a dict consumed by every
downstream training loop.  Per-fold scaling (fit on training only) is
mandatory and lives inside each fold's run_fold function.

Provided functions:
    build_feature_matrices(dataset, phage_feat_df, host_feat_df,
                            feature_meta, cfg, seed) -> dict
    get_device() -> torch.device | None
    scale_train_test(X_train, X_test) -> (X_train_scaled, X_test_scaled, scaler)
        Fits StandardScaler on X_train ONLY then transforms both.
    edge_feature_columns(dataset) -> list[str]
        Canonical ordering of edge-level numeric columns.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

log = logging.getLogger(__name__)


# Numeric columns that historically come from VirusHostInter.csv at the pair
# level (each row is a phage-host pair).
_VHI_NUMERIC = ["k3dist", "k6dist", "GCdiff", "Homology"]


def edge_feature_columns(dataset: pd.DataFrame) -> list[str]:
    """Return canonical edge-level column ordering present in the dataset."""
    pair_feats = ["tetra_corr", "cub_dist", "gc_match", "len_ratio"]
    cols = []
    for c in _VHI_NUMERIC + pair_feats:
        if c in dataset.columns:
            cols.append(c)
    return cols


def get_device():
    """Return a torch.device or None if torch is unavailable."""
    try:
        import torch
        if torch.cuda.is_available():
            return torch.device("cuda")
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    except ImportError:
        return None


def scale_train_test(X_train: np.ndarray, X_test: np.ndarray):
    """Fit a StandardScaler on X_train, then transform X_train and X_test."""
    scaler = StandardScaler()
    if X_train.size == 0:
        return X_train, X_test, scaler
    X_train_s = scaler.fit_transform(X_train).astype(np.float32)
    X_test_s = scaler.transform(X_test).astype(np.float32) if X_test.size else X_test
    # guard against NaN/Inf
    X_train_s = np.where(np.isfinite(X_train_s), X_train_s, 0.0).astype(np.float32)
    if X_test.size:
        X_test_s = np.where(np.isfinite(X_test_s), X_test_s, 0.0).astype(np.float32)
    return X_train_s, X_test_s, scaler


def build_feature_matrices(dataset: pd.DataFrame,
                           phage_feat_df: pd.DataFrame,
                           host_feat_df: pd.DataFrame,
                           feature_meta: dict,
                           cfg: dict,
                           seed: int) -> dict:
    """Assemble RAW feature matrices and an index of phage/host nodes."""
    log.info("[matrices] assembling raw feature matrices")

    phage_list = list(phage_feat_df.index.astype(str))
    host_list = list(host_feat_df.index.astype(str))
    phage2idx = {p: i for i, p in enumerate(phage_list)}
    host2idx = {h: i for i, h in enumerate(host_list)}

    # phage embedding columns: split into "original-name" and "sanitized-name"
    name_cols = feature_meta.get("name_cols", [])
    sanit_cols = feature_meta.get("sanitized_name_cols", [])

    # base phage features WITHOUT either name embedding
    other_phage_cols = [c for c in phage_feat_df.columns
                        if c not in name_cols and c not in sanit_cols
                        and c != "phage_has_seq"]
    other_host_cols = [c for c in host_feat_df.columns if c != "host_has_seq"]

    # Build matrices
    # PHAGE_BASE = [other_phage] + [name_cols]
    # PHAGE_BASE_SANITIZED = [other_phage] + [sanit_cols]
    phage_other = phage_feat_df[other_phage_cols].to_numpy(dtype=np.float32, copy=True)
    phage_name = (phage_feat_df[name_cols].to_numpy(dtype=np.float32, copy=True)
                  if name_cols else np.zeros((len(phage_list), 0), dtype=np.float32))
    phage_sanit = (phage_feat_df[sanit_cols].to_numpy(dtype=np.float32, copy=True)
                   if sanit_cols else np.zeros((len(phage_list), 0), dtype=np.float32))

    PHAGE_BASE = np.hstack([phage_other, phage_name]).astype(np.float32, copy=False)
    PHAGE_BASE_SANITIZED = np.hstack([phage_other, phage_sanit]).astype(np.float32, copy=False)
    HOST_BASE = host_feat_df[other_host_cols].to_numpy(dtype=np.float32, copy=True)

    # Replace NaNs with zeros (safe default)
    PHAGE_BASE = np.where(np.isfinite(PHAGE_BASE), PHAGE_BASE, 0.0).astype(np.float32)
    PHAGE_BASE_SANITIZED = np.where(np.isfinite(PHAGE_BASE_SANITIZED),
                                    PHAGE_BASE_SANITIZED, 0.0).astype(np.float32)
    HOST_BASE = np.where(np.isfinite(HOST_BASE), HOST_BASE, 0.0).astype(np.float32)

    # Edge features (per pair) - RAW
    edge_cols = edge_feature_columns(dataset)
    if edge_cols:
        EDGE_FEATS_RAW = dataset[edge_cols].to_numpy(dtype=np.float32, copy=True)
        EDGE_FEATS_RAW = np.where(np.isfinite(EDGE_FEATS_RAW),
                                  EDGE_FEATS_RAW, 0.0).astype(np.float32)
    else:
        EDGE_FEATS_RAW = np.zeros((len(dataset), 0), dtype=np.float32)

    # Map dataset rows to node indices (in-place add to dataset)
    if "phage_idx" not in dataset.columns:
        dataset["phage_idx"] = dataset["phage"].map(phage2idx).astype("Int64")
    if "host_idx" not in dataset.columns:
        dataset["host_idx"] = dataset["host"].map(host2idx).astype("Int64")

    name_dim = phage_name.shape[1]
    sanit_dim = phage_sanit.shape[1]
    name_start = phage_other.shape[1]  # column offset of name dims inside PHAGE_BASE

    matrices = {
        "PHAGE_BASE": PHAGE_BASE,
        "PHAGE_BASE_SANITIZED": PHAGE_BASE_SANITIZED,
        "HOST_BASE": HOST_BASE,
        "EDGE_FEATS_RAW": EDGE_FEATS_RAW,
        "phage_list": phage_list,
        "host_list": host_list,
        "phage2idx": phage2idx,
        "host2idx": host2idx,
        "n_phages": len(phage_list),
        "n_hosts": len(host_list),
        "all_feat_cols": list(other_phage_cols) + list(name_cols) + list(other_host_cols),
        "ALL_EDGE_FEATS": list(edge_cols),
        "name_dim": int(name_dim),
        "sanit_dim": int(sanit_dim),
        "name_offset": int(name_start),
    }
    log.info(f"[matrices] PHAGE_BASE shape={PHAGE_BASE.shape} "
             f"HOST_BASE shape={HOST_BASE.shape} EDGE shape={EDGE_FEATS_RAW.shape}")
    return matrices


# ---------------------------------------------------------------------------
# Per-row helpers for classical ML and GNN concatenation
# ---------------------------------------------------------------------------
def pair_feature_matrix(dataset: pd.DataFrame, matrices: dict,
                        use_sanitized_names: bool = False) -> np.ndarray:
    """Concatenate per-row [phage_features | host_features | edge_features]."""
    P = matrices["PHAGE_BASE_SANITIZED"] if use_sanitized_names else matrices["PHAGE_BASE"]
    H = matrices["HOST_BASE"]
    E = matrices["EDGE_FEATS_RAW"]
    p_idx = dataset["phage_idx"].astype("Int64").to_numpy()
    h_idx = dataset["host_idx"].astype("Int64").to_numpy()
    # safety: any missing index -> zeros vector
    n = len(dataset)
    p_dim = P.shape[1]
    h_dim = H.shape[1]
    e_dim = E.shape[1]
    out = np.zeros((n, p_dim + h_dim + e_dim), dtype=np.float32)
    for i in range(n):
        pi = p_idx[i]
        hi = h_idx[i]
        if not pd.isna(pi):
            out[i, :p_dim] = P[int(pi)]
        if not pd.isna(hi):
            out[i, p_dim:p_dim + h_dim] = H[int(hi)]
        if e_dim:
            out[i, p_dim + h_dim:] = E[i]
    return out
