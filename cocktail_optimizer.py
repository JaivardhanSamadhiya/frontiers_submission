"""
cocktail_optimizer.py - Greedy phage cocktail selection
=======================================================
Uses LOSO-fold predictions (already computed) rather than re-running models.
All randomness uses np.random.default_rng(seed).

Public:
    run_cocktail(arch, dataset, pred_df, cfg, seed) -> pd.DataFrame
    strain_coverage(per_strain_probs, k, seed) -> float
    greedy_cocktail(per_strain_probs, k, seed) -> list[str]
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def greedy_cocktail(prob_matrix: pd.DataFrame, k: int, seed: int,
                    threshold: float = 0.5) -> list[str]:
    """Pick k phages that maximise the count of strains with >=1 infecting hit.

    prob_matrix: index=phages, columns=hosts (strains), values=P(infection)
    """
    rng = np.random.default_rng(seed)
    chosen: list[str] = []
    if prob_matrix is None or len(prob_matrix) == 0:
        return chosen
    covered = pd.Series(False, index=prob_matrix.columns)
    candidates = list(prob_matrix.index)
    while len(chosen) < k and candidates:
        best_phage = None
        best_gain = -1
        for ph in candidates:
            row = prob_matrix.loc[ph] >= threshold
            gain = int((row & ~covered).sum())
            if gain > best_gain:
                best_gain = gain
                best_phage = ph
        if best_phage is None:
            # tie-break randomly
            best_phage = rng.choice(candidates)
        chosen.append(best_phage)
        covered |= prob_matrix.loc[best_phage] >= threshold
        candidates.remove(best_phage)
        if covered.all():
            break
    return chosen


def strain_coverage(prob_matrix: pd.DataFrame, k: int, seed: int,
                    threshold: float = 0.5) -> dict:
    chosen = greedy_cocktail(prob_matrix, k=k, seed=seed, threshold=threshold)
    if not chosen:
        return {"cocktail": [], "coverage": 0.0, "n_strains": int(prob_matrix.shape[1])
                if hasattr(prob_matrix, "shape") else 0}
    covered = (prob_matrix.loc[chosen] >= threshold).any(axis=0)
    cov = float(covered.mean())
    return {"cocktail": chosen, "coverage": cov,
            "n_strains": int(prob_matrix.shape[1])}


def run_cocktail(arch: str, dataset: pd.DataFrame, pred_df: pd.DataFrame,
                 cfg: dict, seed: int) -> pd.DataFrame:
    """Run cocktail optimization per host-genus using LOSO predictions.

    pred_df columns: species, phage, host, y_true, y_prob
    """
    if pred_df is None or len(pred_df) == 0:
        return pd.DataFrame()

    k = int(cfg["evaluation"]["k_cocktail"])
    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / f"gnn_{arch.lower()}_cocktail.csv"

    # Join genus
    host_to_genus = dataset.drop_duplicates("host").set_index("host")["host_genus"].to_dict()
    pred_df = pred_df.copy()
    pred_df["host_genus"] = pred_df["host"].map(host_to_genus).fillna("unknown")

    rows = []
    for gen, sub in pred_df.groupby("host_genus"):
        if sub["host"].nunique() < 2 or sub["phage"].nunique() < 2:
            continue
        pivot = sub.pivot_table(index="phage", columns="host",
                                values="y_prob", aggfunc="max")
        pivot = pivot.fillna(0.0)
        res = strain_coverage(pivot, k=k, seed=seed)
        rows.append({"genus": gen, "k": k,
                     "n_phages": int(pivot.shape[0]),
                     "n_strains": int(pivot.shape[1]),
                     "mean_coverage": float(res["coverage"]),
                     "coverage_geq_75": int(res["coverage"] >= 0.75),
                     "cocktail": ",".join(res["cocktail"])})
    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    log.info(f"  [cocktail] {arch}: {len(df)} genera evaluated; "
             f"mean coverage={df['mean_coverage'].mean() if len(df) else 0:.3f}")
    return df
