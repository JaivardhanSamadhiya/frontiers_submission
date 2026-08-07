#!/usr/bin/env python3
"""Step 8: phage cocktail optimisation on the predicted susceptibility matrix.

Pipeline:
  1. Leakage-free out-of-fold GBM predictions for every covered (phage, host)
     pair via StratifiedGroupKFold grouped by host genome-cluster.
  2. Build predicted coverage A (pred >= operating threshold) and TRUE coverage
     T (observed label == 1) over tested pairs only.
  3. Optimise cocktails (greedy + exact ILP) and score them on TRUE coverage:
       - coverage-vs-size curve (greedy vs random vs truth-oracle),
       - minimum cocktail to cover all coverable hosts (greedy vs ILP vs oracle),
       - robustness: k-redundant cocktails (k = 1,2,3).

Run (from /tmp, then cd in):
  env -u PYTHONHOME -u PYTHONPATH .venv/bin/python -u experiments/08_cocktail.py
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from sklearn.metrics import f1_score  # noqa: E402
from sklearn.model_selection import StratifiedGroupKFold  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

from precisionphage.cocktail import (  # noqa: E402
    greedy_cover, ilp_max_cover, ilp_min_cover, true_coverage,
    true_coverage_curve,
)
from precisionphage.features.assembly import build_covered_dataset  # noqa: E402
from precisionphage.models import fit_predict_gbm  # noqa: E402
from precisionphage.splits import build_clusters, sketch_entities  # noqa: E402
from precisionphage.utils import (  # noqa: E402
    ensure_dirs, get_logger, limit_threads, load_config, set_determinism,
)

log = get_logger("cocktail")


def _assign_host_clusters(cfg, data):
    sp = cfg["splits"]
    cache = (cfg["paths"]["cache_dir"]
             / f"clusters_k{sp['mash_k']}_d{sp['mash_max_distance']}.json")
    if cache.exists():
        obj = json.loads(cache.read_text())
        if len(obj.get("host", {})) == len(data.hosts):
            return obj["host"]
    h_sk = sketch_entities(data.hosts, data.host_index, sp["mash_k"],
                           sp["minhash_num"], cfg["features"]["n_workers"])
    return build_clusters(data.hosts, h_sk, sp["mash_max_distance"],
                          sp["mash_k"], sp["minhash_num"])


def _oof_predictions(X, y, groups, seed, n_splits=5):
    oof = np.full(len(y), np.nan, dtype=np.float32)
    sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for tr, te in sgkf.split(X, y, groups):
        sc = StandardScaler().fit(X[tr])
        Xtr = np.nan_to_num(sc.transform(X[tr])).astype(np.float32)
        Xte = np.nan_to_num(sc.transform(X[te])).astype(np.float32)
        oof[te] = fit_predict_gbm(Xtr, y[tr], Xte, seed)
    return oof


def main() -> None:
    cfg = load_config(ROOT / "configs" / "default.yaml")
    seed = cfg["seed"]
    set_determinism(seed)
    ensure_dirs(cfg)
    limit_threads(1)
    rd = cfg["paths"]["results_dir"]

    data = build_covered_dataset(cfg)
    cov = data.df.reset_index(drop=True)
    X = data.X_flat()
    y = cov["label"].to_numpy().astype(int)
    hc = _assign_host_clusters(cfg, data)
    groups = cov["host"].map(hc).to_numpy()

    log.info("Generating leakage-free OOF GBM predictions (grouped by host cluster)")
    oof = _oof_predictions(X, y, groups, seed)
    oof = np.nan_to_num(oof, nan=0.0)

    # operating threshold = F1-optimal on OOF
    grid = np.round(np.arange(0.10, 0.91, 0.02), 3)
    f1s = [f1_score(y, (oof >= t).astype(int), zero_division=0) for t in grid]
    thr = float(grid[int(np.argmax(f1s))])
    log.info("OOF F1-optimal threshold = %.3f (F1=%.3f)", thr, max(f1s))

    # build phage x host coverage matrices over TESTED pairs only
    phages = sorted(cov["phage"].unique())
    hosts = sorted(cov["host"].unique())
    pix = {p: i for i, p in enumerate(phages)}
    hix = {h: i for i, h in enumerate(hosts)}
    n_p, n_h = len(phages), len(hosts)
    A = np.zeros((n_p, n_h), dtype=bool)            # predicted coverage
    T = np.zeros((n_p, n_h), dtype=bool)            # true coverage
    for r in range(len(cov)):
        pi = pix[cov.at[r, "phage"]]
        hi = hix[cov.at[r, "host"]]
        if oof[r] >= thr:
            A[pi, hi] = True
        if y[r] == 1:
            T[pi, hi] = True

    targets = np.where(T.sum(0) >= 1)[0]            # coverable target hosts
    log.info("panel: %d phages x %d hosts; coverable targets=%d", n_p, n_h,
             len(targets))

    # --- coverage-vs-size curves ---
    g_order = greedy_cover(A, targets, k=1)
    g_curve = true_coverage_curve(g_order, T, targets, k=1)
    oracle_order = greedy_cover(T, targets, k=1)
    o_curve = true_coverage_curve(oracle_order, T, targets, k=1)
    rng = np.random.default_rng(seed)
    rand_curves = []
    cand = np.where(A.any(1))[0]
    for _ in range(50):
        perm = rng.permutation(cand)
        rand_curves.append(true_coverage_curve(perm[:len(g_order) + 5], T, targets, 1))
    L = min(len(c) for c in rand_curves)
    rand_curve = np.mean([c[:L] for c in rand_curves], axis=0)

    def size_to_reach(curve, frac):
        idx = np.where(curve >= frac)[0]
        return int(idx[0] + 1) if len(idx) else None

    log.info("coverage-vs-size (TRUE coverage of %d targets):", len(targets))
    for frac in (0.5, 0.8, 0.9, 0.95):
        log.info("  reach %.0f%%: greedy=%s  oracle=%s  random(mean)=%s",
                 frac * 100, size_to_reach(g_curve, frac),
                 size_to_reach(o_curve, frac), size_to_reach(rand_curve, frac))

    # --- minimum full-cover cocktail (greedy vs ILP vs oracle) ---
    g_full = greedy_cover(A, targets, k=1)
    ilp_full = ilp_min_cover(A, targets, k=1)
    ilp_oracle = ilp_min_cover(T, targets, k=1)
    log.info("MIN COCKTAIL for full predicted-cover of coverable targets:")
    log.info("  greedy : size=%d  true-coverage=%.3f", len(g_full),
             true_coverage(g_full, T, targets, 1))
    if ilp_full is not None:
        log.info("  ILP    : size=%d  true-coverage=%.3f", len(ilp_full),
                 true_coverage(ilp_full, T, targets, 1))
    if ilp_oracle is not None:
        log.info("  ILP(oracle on truth): size=%d (theoretical minimum)",
                 len(ilp_oracle))

    # --- robustness: k-redundant cocktails ---
    rob_rows = []
    for k in (1, 2, 3):
        sel = greedy_cover(A, targets, k=k)
        tc1 = true_coverage(sel, T, targets, 1)
        tck = true_coverage(sel, T, targets, k)
        rob_rows.append({"k": k, "cocktail_size": len(sel),
                         "true_cover_>=1": round(tc1, 3),
                         "true_cover_>=k": round(tck, 3)})
    rob = pd.DataFrame(rob_rows)
    log.info("ROBUSTNESS (k-redundant greedy cocktails):\n%s",
             rob.to_string(index=False))

    # --- max-coverage under budget (exact ILP) ---
    budget_rows = []
    for B in (5, 10, 20, 30):
        sel, _ = ilp_max_cover(A, targets, budget=B, k=1)
        budget_rows.append({"budget": B,
                            "ilp_true_coverage": round(true_coverage(sel, T, targets, 1), 3),
                            "greedy_true_coverage": round(
                                float(g_curve[min(B, len(g_curve)) - 1]), 3)})
    budget = pd.DataFrame(budget_rows)
    log.info("MAX-COVERAGE under budget (ILP vs greedy, TRUE coverage):\n%s",
             budget.to_string(index=False))

    # --- persist ---
    curve_df = pd.DataFrame({
        "size": np.arange(1, len(g_curve) + 1),
        "greedy": g_curve,
        "oracle": np.pad(o_curve, (0, max(0, len(g_curve) - len(o_curve))),
                         constant_values=o_curve[-1] if len(o_curve) else 0)[:len(g_curve)],
        "random_mean": np.pad(rand_curve, (0, max(0, len(g_curve) - len(rand_curve))),
                              constant_values=rand_curve[-1] if len(rand_curve) else 0)[:len(g_curve)],
    })
    curve_df.to_csv(rd / "cocktail_curve.csv", index=False)
    rob.to_csv(rd / "cocktail_robustness.csv", index=False)
    budget.to_csv(rd / "cocktail_budget.csv", index=False)
    (rd / "cocktail_summary.json").write_text(json.dumps({
        "threshold": thr, "n_phages": n_p, "n_hosts": n_h,
        "n_targets": int(len(targets)),
        "greedy_min_size": len(g_full),
        "ilp_min_size": int(len(ilp_full)) if ilp_full is not None else None,
        "ilp_oracle_min_size": int(len(ilp_oracle)) if ilp_oracle is not None else None,
        "robustness": rob_rows, "budget": budget_rows,
    }, indent=2, default=float))

    # --- figure ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 5))
        s = curve_df["size"]
        ax.plot(s, curve_df["oracle"], label="Oracle (truth-greedy)", lw=2, color="#2ca02c")
        ax.plot(s, curve_df["greedy"], label="Model-driven greedy", lw=2, color="#1f77b4")
        ax.plot(s, curve_df["random_mean"], label="Random (mean of 50)", lw=2,
                ls="--", color="#888888")
        if ilp_full is not None:
            ax.axvline(len(ilp_full), color="#d62728", ls=":",
                       label=f"ILP min cocktail (n={len(ilp_full)})")
        ax.set_xlabel("Cocktail size (number of phages)")
        ax.set_ylabel("Fraction of target hosts truly covered")
        ax.set_title("Phage cocktail coverage vs size")
        ax.set_ylim(0, 1.02); ax.legend(loc="lower right"); ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(rd / "cocktail_coverage.png", dpi=150)
        log.info("Wrote figure %s", rd / "cocktail_coverage.png")
    except Exception as e:
        log.warning("figure skipped (%s)", e)

    log.info("Wrote cocktail_curve.csv + cocktail_*.csv + cocktail_summary.json")


if __name__ == "__main__":
    main()
