"""
logo_validation.py - Leave-One-Genus-Out cross-validation
=========================================================
Used by gnn.run_gnn_pipeline and by the ablation/ensemble modules.
Per-fold scaler is fit on the training fold only; validate_no_leakage
is invoked at the start of every fold.

Public:
    run_logo(arch, dataset, matrices, valid_genera, cfg, seed, fold_fn) -> pd.DataFrame
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from build_labels import DataLeakageError, validate_no_leakage
from evaluation import metrics as compute_metrics

log = logging.getLogger(__name__)


def run_logo(arch: str, dataset: pd.DataFrame, matrices: dict,
             valid_genera: list[str], cfg: dict, seed: int,
             fold_fn) -> pd.DataFrame:
    """Run LOGO-CV using a caller-provided fold function.

    fold_fn(train_mask, test_mask) -> (probs, y_true, meta)
    """
    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / f"gnn_{arch.lower()}_logo.csv"
    rows = []
    for gen in valid_genera:
        test_mask = (dataset["host_genus"] == gen).values
        train_mask = ~test_mask
        if test_mask.sum() == 0 or train_mask.sum() == 0:
            continue
        train_df = dataset[train_mask].reset_index(drop=True)
        test_df = dataset[test_mask].reset_index(drop=True)
        try:
            validate_no_leakage(train_df, test_df, context=f"{arch} LOGO {gen}")
        except DataLeakageError as e:
            log.error(f"  Skipping LOGO genus {gen}: {e}")
            continue
        try:
            probs, y_true, meta = fold_fn(train_mask, test_mask)
        except Exception as e:
            log.warning(f"  LOGO fold failed for {gen}: {e}")
            continue
        if meta.get("skipped") or len(probs) == 0:
            continue
        m = compute_metrics(y_true, probs)
        rows.append({"genus": gen, "n_test": int(test_mask.sum()), **m})
        pd.DataFrame(rows).to_csv(out_path, index=False)
    return pd.DataFrame(rows)
