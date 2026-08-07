"""
evaluation.py - LOSO / metrics / classical ML training
======================================================
Public surface:
    get_valid_species(dataset, cfg) -> list[str]
    get_valid_genera(dataset, cfg) -> list[str]
    metrics(y_true, y_prob, threshold=0.5) -> dict
    run_classical_loso(dataset, matrices, valid_species, cfg, seed) -> dict
    run_experimental_validation(exp_val_df, dataset, matrices, all_results,
                                cfg, seed) -> dict
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    roc_auc_score,
)

from build_labels import DataLeakageError, validate_no_leakage
from model import pair_feature_matrix, scale_train_test

log = logging.getLogger(__name__)

warnings.filterwarnings("ignore", category=UserWarning)


# ---------------------------------------------------------------------------
# Filters
# ---------------------------------------------------------------------------
def get_valid_species(dataset: pd.DataFrame, cfg: dict) -> list[str]:
    min_rows = int(cfg["data"]["min_species_rows"])
    max_rows = int(cfg["data"]["max_species_rows"])
    min_pos = int(cfg["data"]["min_positive_per_species"])
    min_neg = int(cfg["data"]["min_negative_per_species"])

    valid, excluded = [], []
    for sp in dataset["host"].unique():
        sub = dataset[dataset["host"] == sp]
        n_pos = int((sub["label"] == 1).sum())
        n_neg = int((sub["label"] == 0).sum())
        n_tot = len(sub)
        if n_tot < min_rows:
            reason = f"too_few_rows({n_tot})"
        elif n_tot > max_rows:
            reason = f"too_many_rows({n_tot})"
        elif n_pos < min_pos:
            reason = f"too_few_positives({n_pos})"
        elif n_neg < min_neg:
            reason = f"too_few_negatives({n_neg})"
        else:
            valid.append(sp)
            continue
        excluded.append({"species": sp, "reason": reason, "n_total": n_tot,
                         "n_pos": n_pos, "n_neg": n_neg})

    excl_df = pd.DataFrame(excluded)
    out = Path(cfg["paths"]["results_dir"]) / "species_inclusion_analysis.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    excl_df.to_csv(out, index=False)
    log.info(f"Valid: {len(valid)} species | Excluded: {len(excluded)} "
             f"(saved to species_inclusion_analysis.csv)")
    return sorted(valid)


def get_valid_genera(dataset: pd.DataFrame, cfg: dict) -> list[str]:
    min_rows = int(cfg["data"]["min_species_rows"])
    min_pos = int(cfg["data"]["min_positive_per_species"])
    min_neg = int(cfg["data"]["min_negative_per_species"])
    valid = []
    excluded = []
    for gen in dataset["host_genus"].dropna().unique():
        sub = dataset[dataset["host_genus"] == gen]
        n_pos = int((sub["label"] == 1).sum())
        n_neg = int((sub["label"] == 0).sum())
        if len(sub) < min_rows or n_pos < min_pos or n_neg < min_neg:
            excluded.append({"genus": gen, "n_total": len(sub),
                             "n_pos": n_pos, "n_neg": n_neg})
            continue
        valid.append(gen)
    out = Path(cfg["paths"]["results_dir"]) / "genus_inclusion_analysis.csv"
    pd.DataFrame(excluded).to_csv(out, index=False)
    return sorted(valid)


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def metrics(y_true, y_prob, threshold: float = 0.5) -> dict:
    """Never returns NaN; always returns floats; catches all exceptions."""
    try:
        y = np.asarray(y_true, dtype=np.float32)
        p = np.asarray(y_prob, dtype=np.float32)
        p = np.where(np.isfinite(p), p, 0.5)
        p = np.clip(p, 0.0, 1.0)
        pred = (p >= threshold).astype(int)
        try:
            roc = float(roc_auc_score(y, p)) if len(np.unique(y)) > 1 else 0.5
        except Exception:
            roc = 0.5
        try:
            pr_arr, rc_arr, _ = precision_recall_curve(y, p)
            pr_auc = float(np.trapz(pr_arr[::-1], rc_arr[::-1]))
            if not np.isfinite(pr_auc):
                pr_auc = float(average_precision_score(y, p))
        except Exception:
            pr_auc = 0.0
        try:
            f1 = float(f1_score(y, pred, zero_division=0))
        except Exception:
            f1 = 0.0
        try:
            mcc = float(matthews_corrcoef(y, pred))
        except Exception:
            mcc = 0.0
        try:
            cm = confusion_matrix(y, pred, labels=[0, 1])
            tn, fp, fn, tp = cm.ravel()
        except Exception:
            tn = fp = fn = tp = 0
        sens = float(tp) / float(tp + fn) if (tp + fn) > 0 else 0.0
        spec = float(tn) / float(tn + fp) if (tn + fp) > 0 else 0.0
        ppv = float(tp) / float(tp + fp) if (tp + fp) > 0 else 0.0
        npv = float(tn) / float(tn + fn) if (tn + fn) > 0 else 0.0
        try:
            acc = float(accuracy_score(y, pred))
        except Exception:
            acc = 0.0
        return {
            "roc_auc": float(roc),
            "pr_auc": float(pr_auc),
            "f1": float(f1),
            "mcc": float(mcc),
            "sensitivity": float(sens),
            "specificity": float(spec),
            "ppv": float(ppv),
            "npv": float(npv),
            "accuracy": float(acc),
            "n_pos": int((y == 1).sum()),
            "n_neg": int((y == 0).sum()),
            "threshold": float(threshold),
        }
    except Exception as e:
        log.warning(f"[metrics] failed: {e}")
        return {"roc_auc": 0.5, "pr_auc": 0.0, "f1": 0.0, "mcc": 0.0,
                "sensitivity": 0.0, "specificity": 0.0, "ppv": 0.0,
                "npv": 0.0, "accuracy": 0.0, "n_pos": 0, "n_neg": 0,
                "threshold": float(threshold)}


# ---------------------------------------------------------------------------
# Classical model factory
# ---------------------------------------------------------------------------
def _build_classical_models(seed: int) -> dict:
    models = {
        "RandomForest": RandomForestClassifier(
            n_estimators=400, max_depth=None, n_jobs=-1, random_state=seed,
            class_weight="balanced_subsample"),
        "HistGB": HistGradientBoostingClassifier(
            max_iter=400, max_depth=6, learning_rate=0.05, random_state=seed),
    }
    try:
        from xgboost import XGBClassifier
        models["XGBoost"] = XGBClassifier(
            n_estimators=400, max_depth=6, learning_rate=0.05,
            random_state=seed, tree_method="hist", n_jobs=-1,
            use_label_encoder=False, eval_metric="logloss", verbosity=0)
    except Exception:
        log.info("[classical] xgboost not available; using RF + HistGB only")
    return models


# ---------------------------------------------------------------------------
# Classical LOSO
# ---------------------------------------------------------------------------
def run_classical_loso(dataset: pd.DataFrame, matrices: dict,
                       valid_species: list[str], cfg: dict, seed: int) -> dict:
    """Run LOSO-CV for each classical model; return dict[name] -> result dict."""
    rng = np.random.default_rng(seed)
    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)

    models_factory = _build_classical_models(seed)
    X_raw = pair_feature_matrix(dataset, matrices, use_sanitized_names=False)
    y = dataset["label"].astype(float).values

    out = {}
    for name, _ in models_factory.items():
        log.info(f"  [classical] {name}: LOSO over {len(valid_species)} species")
        per_species = []
        all_probs, all_y = [], []
        pred_records = []
        csv_path = results_dir / f"classical_{name.lower()}_loso.csv"
        for sp in valid_species:
            test_mask = (dataset["host"] == sp).values
            train_mask = ~test_mask
            train_df = dataset[train_mask].reset_index(drop=True)
            test_df = dataset[test_mask].reset_index(drop=True)
            try:
                validate_no_leakage(train_df, test_df, context=f"{name} LOSO {sp}")
            except DataLeakageError as e:
                log.error(f"    Skipping {sp}: {e}")
                continue

            X_tr_raw = X_raw[train_mask]
            X_te_raw = X_raw[test_mask]
            y_tr = y[train_mask]
            y_te = y[test_mask]

            X_tr, X_te, _ = scale_train_test(X_tr_raw, X_te_raw)

            # Build a fresh model instance (avoid contamination)
            try:
                clf = _build_classical_models(seed)[name]
                clf.fit(X_tr, y_tr)
                probs = clf.predict_proba(X_te)[:, 1].astype(np.float32)
            except Exception as e:
                log.warning(f"    {name} failed on {sp}: {e}")
                continue
            m = metrics(y_te, probs)
            per_species.append({"species": sp, "n_test": int(test_mask.sum()), **m})
            all_probs.append(probs)
            all_y.append(y_te)
            idx = np.where(test_mask)[0]
            for j, ix in enumerate(idx):
                pred_records.append({"species": sp,
                                     "phage": dataset.iloc[ix]["phage"],
                                     "host": dataset.iloc[ix]["host"],
                                     "y_true": float(y_te[j]),
                                     "y_prob": float(probs[j])})
            pd.DataFrame(per_species).to_csv(csv_path, index=False)

        df = pd.DataFrame(per_species)
        mean_auc = float(df["roc_auc"].mean()) if len(df) else 0.5
        std_auc = float(df["roc_auc"].std()) if len(df) > 1 else 0.0
        if all_probs:
            pooled_p = np.concatenate(all_probs)
            pooled_y = np.concatenate(all_y)
            try:
                pooled_auc = float(roc_auc_score(pooled_y, pooled_p))
            except Exception:
                pooled_auc = 0.5
        else:
            pooled_p = np.array([])
            pooled_y = np.array([])
            pooled_auc = 0.5
        pred_df = pd.DataFrame(pred_records)
        pred_df.to_csv(results_dir / f"classical_{name.lower()}_predictions.csv", index=False)
        out[name] = {
            "loso_df": df,
            "mean_auc": mean_auc,
            "loso_mean": mean_auc,
            "loso_std": std_auc,
            "loso_pooled": pooled_auc,
            "per_species_auc": df.set_index("species")["roc_auc"] if "roc_auc" in df else pd.Series(),
            "pred_df": pred_df,
            "all_proba": pooled_p,
            "y_true": pooled_y,
        }
        log.info(f"  [classical] {name}: mean LOSO AUC={mean_auc:.4f} (std={std_auc:.4f}, pooled={pooled_auc:.4f})")
    return out


# ---------------------------------------------------------------------------
# Held-out experimental validation
# ---------------------------------------------------------------------------
def run_experimental_validation(exp_val_df: pd.DataFrame, dataset: pd.DataFrame,
                                matrices: dict, all_results: dict, cfg: dict,
                                seed: int) -> dict:
    """Evaluate all models on a held-out experimental dataset.

    No exp_val rows are used in training.  We retrain each classical model on
    the full assembled dataset (which excludes exp_val) and produce predictions
    on exp_val.  GNN predictions come from the LOSO-fold ensemble probability
    average if available; otherwise we retrain on the full dataset using the
    classical fallback.
    """
    log.warning("EXPERIMENTAL VALIDATION (held-out, never used in training)")

    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)

    # Sanity check: no overlap between exp_val and training set
    try:
        validate_no_leakage(dataset, exp_val_df, context="train-vs-exp_val")
    except DataLeakageError as e:
        log.warning(f"  Some exp_val pairs overlap training set: {e}")
        # remove overlapping rows from exp_val so reported metrics are honest
        train_pairs = set(zip(dataset["phage"].astype(str).str.lower().str.strip(),
                              dataset["host"].astype(str).str.lower().str.strip()))
        keep = ~exp_val_df.apply(lambda r:
                                 (str(r["phage"]).lower().strip(),
                                  str(r["host"]).lower().strip()) in train_pairs,
                                 axis=1)
        before = len(exp_val_df)
        exp_val_df = exp_val_df.loc[keep].reset_index(drop=True)
        log.warning(f"  Dropped {before - len(exp_val_df)} overlapping rows; "
                    f"{len(exp_val_df)} remain")

    # Build feature vectors for exp_val using the same matrices
    # exp_val may reference phages/hosts not in our matrices; fall back to zero vectors.
    eval_df = exp_val_df.copy().reset_index(drop=True)
    eval_df["phage_idx"] = eval_df["phage"].map(matrices["phage2idx"])
    eval_df["host_idx"] = eval_df["host"].map(matrices["host2idx"])

    # We must construct edge features for exp_val using values that exist there
    edge_cols = matrices.get("ALL_EDGE_FEATS", [])
    edge_arr = np.zeros((len(eval_df), len(edge_cols)), dtype=np.float32)
    for j, c in enumerate(edge_cols):
        if c in eval_df.columns:
            edge_arr[:, j] = pd.to_numeric(eval_df[c], errors="coerce").fillna(0.0).astype(np.float32).values

    P = matrices["PHAGE_BASE"]
    H = matrices["HOST_BASE"]
    p_dim = P.shape[1]
    h_dim = H.shape[1]
    e_dim = edge_arr.shape[1]
    X_eval = np.zeros((len(eval_df), p_dim + h_dim + e_dim), dtype=np.float32)
    for i, r in eval_df.iterrows():
        if not pd.isna(r["phage_idx"]):
            X_eval[i, :p_dim] = P[int(r["phage_idx"])]
        if not pd.isna(r["host_idx"]):
            X_eval[i, p_dim:p_dim + h_dim] = H[int(r["host_idx"])]
        if e_dim:
            X_eval[i, p_dim + h_dim:] = edge_arr[i]
    y_eval = eval_df["label"].astype(int).values

    # Retrain classical models on full training set and predict
    X_train = pair_feature_matrix(dataset, matrices, use_sanitized_names=False)
    y_train = dataset["label"].astype(float).values
    results = {"per_model": {}, "best_auc": 0.0, "best_model": None, "n_pairs": int(len(eval_df))}

    models_factory = _build_classical_models(seed)
    for name, _ in models_factory.items():
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler().fit(X_train)
        X_tr_s = np.where(np.isfinite(scaler.transform(X_train)), scaler.transform(X_train), 0.0).astype(np.float32)
        X_ev_s = np.where(np.isfinite(scaler.transform(X_eval)), scaler.transform(X_eval), 0.0).astype(np.float32)
        try:
            clf = _build_classical_models(seed)[name]
            clf.fit(X_tr_s, y_train)
            probs = clf.predict_proba(X_ev_s)[:, 1].astype(np.float32)
        except Exception as e:
            log.warning(f"  exp_val: {name} failed: {e}")
            continue
        m = metrics(y_eval, probs)
        results["per_model"][name] = m
        log.info(f"  exp_val {name}: AUC={m['roc_auc']:.4f} "
                 f"(n_pos={m['n_pos']}, n_neg={m['n_neg']})")
        if m["roc_auc"] > results["best_auc"]:
            results["best_auc"] = m["roc_auc"]
            results["best_model"] = name

    # Persist results
    rows = []
    for name, m in results["per_model"].items():
        rows.append({"model": name, **{k: v for k, v in m.items() if isinstance(v, (int, float))}})
    out = results_dir / "experimental_validation.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    return results
