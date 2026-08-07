"""
ablation_study.py - Feature and name-leakage ablation studies
==============================================================
- run_ablation: classical (RF) ablation over feature subsets across LOSO folds
  with per-fold scaler fit on training only.
- run_name_leakage_ablation: critical for publication peer review.

Both functions ALWAYS run, even when GNN deps are absent (they internally
use a classical fallback RandomForest).
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from build_labels import DataLeakageError, validate_no_leakage
from model import pair_feature_matrix

log = logging.getLogger(__name__)


# Feature index helpers
def _feature_groups(matrices: dict, dataset: pd.DataFrame) -> dict[str, np.ndarray]:
    """Return boolean column masks within the concatenated pair-feature matrix.

    Pair matrix layout (see model.pair_feature_matrix):
        [phage_features (p_dim) | host_features (h_dim) | edge_features (e_dim)]
    """
    p_dim = matrices["PHAGE_BASE"].shape[1]
    h_dim = matrices["HOST_BASE"].shape[1]
    e_dim = matrices["EDGE_FEATS_RAW"].shape[1]
    total = p_dim + h_dim + e_dim

    name_offset = matrices.get("name_offset", p_dim)
    name_dim = matrices.get("name_dim", 0)

    masks = {
        "phage_genomic": np.zeros(total, dtype=bool),  # p_di/p_tri/p_tet/p_cub etc.
        "phage_name": np.zeros(total, dtype=bool),
        "host_block": np.zeros(total, dtype=bool),
        "edge_block": np.zeros(total, dtype=bool),
    }
    masks["phage_genomic"][:name_offset] = True
    masks["phage_name"][name_offset:name_offset + name_dim] = True
    masks["host_block"][p_dim:p_dim + h_dim] = True
    masks["edge_block"][p_dim + h_dim:] = True

    # Sub-edge masks
    edge_cols = matrices.get("ALL_EDGE_FEATS", [])
    baseline_cols = ["k3dist", "k6dist", "GCdiff", "Homology"]
    pair_cols = ["tetra_corr", "cub_dist", "gc_match", "len_ratio"]
    masks["edge_baseline"] = np.zeros(total, dtype=bool)
    masks["edge_pair"] = np.zeros(total, dtype=bool)
    for k, c in enumerate(edge_cols):
        idx = p_dim + h_dim + k
        if c in baseline_cols:
            masks["edge_baseline"][idx] = True
        elif c in pair_cols:
            masks["edge_pair"][idx] = True
    return masks


def _rf_loso(dataset: pd.DataFrame, X_raw: np.ndarray,
             valid_species: list[str], seed: int, label: str) -> pd.DataFrame:
    """Run LOSO-CV with a fresh RandomForest at each fold."""
    rows = []
    y = dataset["label"].astype(float).values
    for sp in valid_species:
        test_mask = (dataset["host"] == sp).values
        train_mask = ~test_mask
        if test_mask.sum() == 0 or train_mask.sum() == 0:
            continue
        train_df = dataset[train_mask].reset_index(drop=True)
        test_df = dataset[test_mask].reset_index(drop=True)
        try:
            validate_no_leakage(train_df, test_df, context=f"ablation {label} {sp}")
        except DataLeakageError as e:
            log.error(f"  Skipping ablation fold {sp}: {e}")
            continue
        scaler = StandardScaler()
        Xtr = scaler.fit_transform(X_raw[train_mask]).astype(np.float32)
        Xte = scaler.transform(X_raw[test_mask]).astype(np.float32)
        Xtr = np.where(np.isfinite(Xtr), Xtr, 0.0)
        Xte = np.where(np.isfinite(Xte), Xte, 0.0)
        clf = RandomForestClassifier(n_estimators=300, n_jobs=-1, random_state=seed,
                                     class_weight="balanced_subsample")
        try:
            clf.fit(Xtr, y[train_mask])
            probs = clf.predict_proba(Xte)[:, 1]
        except Exception as e:
            log.warning(f"  ablation {label} on {sp} failed: {e}")
            continue
        try:
            auc = float(roc_auc_score(y[test_mask], probs)) if len(np.unique(y[test_mask])) > 1 else 0.5
        except Exception:
            auc = 0.5
        rows.append({"species": sp, "auc": auc,
                     "n_test": int(test_mask.sum()),
                     "ablation": label})
    return pd.DataFrame(rows)


def run_ablation(dataset: pd.DataFrame, matrices: dict,
                 valid_species: list[str], cfg: dict, seed: int,
                 gnn_available: bool = False) -> dict:
    """Run classical RF ablation across canonical feature subsets."""
    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / "ablation_results.csv"

    X_full = pair_feature_matrix(dataset, matrices, use_sanitized_names=False)
    masks = _feature_groups(matrices, dataset)
    total_cols = X_full.shape[1]

    conditions: dict[str, np.ndarray] = {}
    # Full: keep everything
    conditions["full"] = np.ones(total_cols, dtype=bool)
    # No genomic (drop phage_genomic prefix subset i.e. di/tri/tet/cub on phage)
    cond = np.ones(total_cols, dtype=bool)
    cond &= ~masks["phage_genomic"]
    conditions["no_phage_genomic"] = cond
    # No pair features
    cond = np.ones(total_cols, dtype=bool)
    cond &= ~masks["edge_pair"]
    conditions["no_pair_features"] = cond
    # Baseline only
    cond = np.zeros(total_cols, dtype=bool)
    cond |= masks["edge_baseline"]
    conditions["baseline_only"] = cond if cond.any() else np.ones(total_cols, dtype=bool)

    # Hard negatives only: filter dataset to within-genus constructed negatives + positives
    if "source" in dataset.columns:
        hard_mask = ((dataset["label"] == 1) |
                     (dataset["source"] == "constructed_within_genus")).values
    else:
        hard_mask = np.ones(len(dataset), dtype=bool)

    all_rows: list[pd.DataFrame] = []
    for cond_name, cmask in conditions.items():
        if cmask.sum() == 0:
            continue
        df = _rf_loso(dataset, X_full[:, cmask], valid_species, seed, cond_name)
        all_rows.append(df)

    # Hard-negatives condition reuses the full feature set
    if hard_mask.sum() > 0 and hard_mask.sum() < len(dataset):
        subset = dataset[hard_mask].reset_index(drop=True)
        X_sub = X_full[hard_mask]
        df = _rf_loso(subset, X_sub, valid_species, seed, "hard_negatives_only")
        all_rows.append(df)

    if all_rows:
        out_df = pd.concat(all_rows, ignore_index=True)
    else:
        out_df = pd.DataFrame(columns=["species", "auc", "n_test", "ablation"])
    out_df.to_csv(out_path, index=False)

    summary = (out_df.groupby("ablation")["auc"]
               .agg(["mean", "std", "count"]).reset_index()
               .rename(columns={"mean": "mean_auc", "std": "std_auc"}))
    summary.to_csv(results_dir / "ablation_summary.csv", index=False)
    for _, r in summary.iterrows():
        log.info(f"  ablation {r['ablation']}: mean AUC={r['mean_auc']:.4f} "
                 f"(n={int(r['count'])})")

    return {"per_fold": out_df, "summary": summary}


def run_name_leakage_ablation(dataset: pd.DataFrame, matrices: dict,
                              valid_species: list[str], cfg: dict, seed: int,
                              gnn_available: bool = False) -> dict:
    """Compare AUC with original names, no names, and sanitized names."""
    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / "name_leakage_ablation.csv"

    X_orig = pair_feature_matrix(dataset, matrices, use_sanitized_names=False)
    X_sanit = pair_feature_matrix(dataset, matrices, use_sanitized_names=True)

    p_dim = matrices["PHAGE_BASE"].shape[1]
    name_offset = matrices.get("name_offset", p_dim)
    name_dim = matrices.get("name_dim", 0)

    # Without-names: zero out the name slice in X_orig
    X_no_name = X_orig.copy()
    if name_dim > 0:
        X_no_name[:, name_offset:name_offset + name_dim] = 0.0

    df_with = _rf_loso(dataset, X_orig, valid_species, seed, "with_names")
    df_without = _rf_loso(dataset, X_no_name, valid_species, seed, "without_names")
    df_sanit = _rf_loso(dataset, X_sanit, valid_species, seed, "sanitized_names")
    df = pd.concat([df_with, df_without, df_sanit], ignore_index=True)
    df.to_csv(out_path, index=False)

    aucs = {label: float(sub["auc"].mean()) if len(sub) else 0.5
            for label, sub in df.groupby("ablation")}
    a_w = aucs.get("with_names", 0.5)
    a_s = aucs.get("sanitized_names", 0.5)
    a_n = aucs.get("without_names", 0.5)
    drop = a_w - a_s

    # Wilcoxon: paired by species
    try:
        merged = df_with.merge(df_sanit, on="species", suffixes=("_with", "_sanit"))
        if len(merged) >= 5:
            stat, p_val = wilcoxon(merged["auc_with"].values,
                                   merged["auc_sanit"].values)
            p_val = float(p_val)
        else:
            stat, p_val = (np.nan, np.nan)
    except Exception:
        stat, p_val = (np.nan, np.nan)

    if drop < 0.02:
        conclusion = ("Name embeddings do not introduce meaningful label leakage "
                      f"(drop={drop:.4f} < 0.02).")
    else:
        conclusion = ("WARNING: Name embeddings may inflate performance. "
                      f"AUC drop from with_names to sanitized_names = {drop:.4f} "
                      ">= 0.02.  Re-evaluate with sanitized names for publication.")
        log.warning(conclusion)

    out = {
        "with_names": a_w,
        "without_names": a_n,
        "sanitized_names": a_s,
        "auc_drop_with_minus_sanit": drop,
        "wilcoxon_stat": float(stat) if np.isfinite(stat) else None,
        "wilcoxon_p": p_val if np.isfinite(p_val) else None,
        "conclusion": conclusion,
        "per_fold": df,
    }
    return out
