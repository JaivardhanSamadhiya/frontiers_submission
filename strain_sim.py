"""
strain_sim.py - Monte-Carlo unseen-strain simulation
====================================================
Each MC round uses np.random.default_rng(seed + round_number) so rounds are
reproducible and independent.  Per-fold leakage check is enforced.

Public:
    run_unseen_strain(arch, dataset, matrices, cfg, seed, fold_fn) -> pd.DataFrame
        Returns a DataFrame with columns: round, auc, n_unseen, n_train
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from build_labels import DataLeakageError, validate_no_leakage
from evaluation import metrics as compute_metrics

log = logging.getLogger(__name__)


def run_unseen_strain(arch: str, dataset: pd.DataFrame, matrices: dict,
                      cfg: dict, seed: int, fold_fn) -> pd.DataFrame:
    n_rounds = int(cfg["evaluation"]["n_mc_rounds"])
    unseen_frac = float(cfg["evaluation"]["unseen_frac"])
    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / f"gnn_{arch.lower()}_unseen.csv"

    all_phages = sorted(dataset["phage"].unique().tolist())
    rows = []
    for r in range(n_rounds):
        rng = np.random.default_rng(seed + r)
        k = max(1, int(len(all_phages) * unseen_frac))
        unseen = set(rng.choice(all_phages, size=k, replace=False).tolist())
        test_mask = dataset["phage"].isin(unseen).values
        train_mask = ~test_mask
        if test_mask.sum() == 0 or train_mask.sum() == 0:
            continue
        train_df = dataset[train_mask].reset_index(drop=True)
        test_df = dataset[test_mask].reset_index(drop=True)
        try:
            validate_no_leakage(train_df, test_df, context=f"{arch} unseen-strain round {r}")
        except DataLeakageError as e:
            log.error(f"  Skipping MC round {r}: {e}")
            continue
        try:
            probs, y_true, meta = fold_fn(train_mask, test_mask)
        except Exception as e:
            log.warning(f"  MC round {r} failed: {e}")
            continue
        if meta.get("skipped") or len(probs) == 0:
            continue
        m = compute_metrics(y_true, probs)
        rows.append({"round": int(r), "auc": float(m["roc_auc"]),
                     "n_unseen": int(test_mask.sum()),
                     "n_train": int(train_mask.sum())})
        pd.DataFrame(rows).to_csv(out_path, index=False)
    return pd.DataFrame(rows)
