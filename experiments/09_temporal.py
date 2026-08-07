#!/usr/bin/env python3
"""Step 9: eco-evolutionary temporal simulation of phage therapy strategies.

Using leakage-free OOF GBM susceptibility predictions, we pick a multi-strain
"infection" panel (the host genome-cluster with the most coverable strains) and
simulate four strategies forward in time with resistance evolution:
  * no phage (control),
  * best single phage (monophage),
  * model-designed cocktail (greedy, k=1),
  * robust cocktail (greedy, k=2 redundancy).

Resistance to the cocktail requires independent resistance to every targeting
phage (mu_eff = mu ** n_targeting), so cocktails suppress resistance emergence.
Outputs trajectories, a summary metrics table, and a figure.

Run (from /tmp, then cd in):
  env -u PYTHONHOME -u PYTHONPATH .venv/bin/python -u experiments/09_temporal.py
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

from precisionphage.cocktail import greedy_cover  # noqa: E402
from precisionphage.features.assembly import build_covered_dataset  # noqa: E402
from precisionphage.models import fit_predict_gbm  # noqa: E402
from precisionphage.splits import build_clusters, sketch_entities  # noqa: E402
from precisionphage.temporal import TherapyParams, simulate, therapy_metrics  # noqa: E402
from precisionphage.utils import (  # noqa: E402
    ensure_dirs, get_logger, limit_threads, load_config, set_determinism,
)

log = get_logger("temporal")


def _host_clusters(cfg, data):
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


def _oof(X, y, groups, seed, n_splits=5):
    oof = np.full(len(y), np.nan, dtype=np.float32)
    sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for tr, te in sgkf.split(X, y, groups):
        sc = StandardScaler().fit(X[tr])
        Xtr = np.nan_to_num(sc.transform(X[tr])).astype(np.float32)
        Xte = np.nan_to_num(sc.transform(X[te])).astype(np.float32)
        oof[te] = fit_predict_gbm(Xtr, y[tr], Xte, seed)
    return np.nan_to_num(oof, nan=0.0)


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
    hc = _host_clusters(cfg, data)
    groups = cov["host"].map(hc).to_numpy()
    oof = _oof(X, y, groups, seed)

    grid = np.round(np.arange(0.10, 0.91, 0.02), 3)
    thr = float(grid[int(np.argmax([f1_score(y, (oof >= t).astype(int),
                                              zero_division=0) for t in grid]))])

    phages = sorted(cov["phage"].unique())
    hosts = sorted(cov["host"].unique())
    pix = {p: i for i, p in enumerate(phages)}
    hix = {h: i for i, h in enumerate(hosts)}
    n_p, n_h = len(phages), len(hosts)
    Aprob = np.zeros((n_p, n_h), dtype=float)
    T = np.zeros((n_p, n_h), dtype=bool)
    for r in range(len(cov)):
        pi, hi = pix[cov.at[r, "phage"]], hix[cov.at[r, "host"]]
        Aprob[pi, hi] = max(Aprob[pi, hi], float(oof[r]))
        if y[r] == 1:
            T[pi, hi] = True
    A = Aprob >= thr

    # infection panel = the best-covered host strains (each with multiple
    # candidate phages) so cocktails can provide per-host redundancy and the
    # resistance-prevention benefit of k>1 cocktails is observable.
    coverable = T.sum(0) >= 1
    cand_count = A.sum(0)                              # predicted phages per host
    elig = np.where(coverable & (cand_count >= 3))[0]
    if len(elig) == 0:                                 # fallback: most-covered
        elig = np.where(coverable)[0]
    panel = elig[np.argsort(cand_count[elig])[::-1][:5]]
    log.info("infection panel: %d strains (candidate phages/strain=%s): %s",
             len(panel), cand_count[panel].tolist(),
             [hosts[i][:24] for i in panel])

    # candidate phages with >=1 predicted edge to the panel
    cand = np.where(A[:, panel].any(1))[0]
    Apanel = A[np.ix_(cand, panel)]
    Pprob = Aprob[np.ix_(cand, panel)]
    panel_local = np.arange(len(panel))

    # strategies -> selected phage rows (indices into cand) and susceptibility
    mono = [int(np.argmax(Apanel.sum(1)))]
    greedy1 = greedy_cover(Apanel, panel_local, k=1)
    greedy2 = greedy_cover(Apanel, panel_local, k=2)
    strategies = {
        "control": np.zeros((1, len(panel))),
        "monophage": Pprob[mono, :],
        "cocktail_k1": Pprob[greedy1, :] if greedy1 else Pprob[mono, :],
        "robust_k2": Pprob[greedy2, :] if greedy2 else Pprob[mono, :],
    }
    sizes = {"control": 0, "monophage": 1, "cocktail_k1": len(greedy1),
             "robust_k2": len(greedy2)}

    S0 = np.full(len(panel), 1e6)
    rows, traj = [], {}
    for name, A_sel in strategies.items():
        pp = TherapyParams(dose=0.0 if name == "control" else 1e8)
        out = simulate(A_sel, S0, pp)
        m = therapy_metrics(out, S0.sum())
        traj[name] = {"t": out["t"], "total": out["total"],
                      "R": out["R"].sum(0)}
        rows.append({"strategy": name, "n_phages": sizes[name],
                     "end_load": m["end_load"], "nadir": m["nadir"],
                     "log10_drop": round(m["log10_drop"], 2),
                     "resistant_frac_end": round(m["resistant_fraction_end"], 3),
                     "rebound": m["rebound"]})
    tbl = pd.DataFrame(rows)
    log.info("THERAPY OUTCOMES (panel of %d strains):\n%s", len(panel),
             tbl.to_string(index=False))

    tbl.to_csv(rd / "temporal_outcomes.csv", index=False)
    (rd / "temporal_summary.json").write_text(json.dumps(
        {"threshold": thr, "panel_strains": [hosts[i] for i in panel],
         "panel_size": int(len(panel)), "outcomes": rows}, indent=2, default=float))

    # Persist trajectories so figures can be regenerated without rerunning simulation
    traj_save = {k: v for k, v in traj.items()}
    np.savez(rd / "temporal_trajectory.npz", **traj_save)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        colors = {"control": "#7f7f7f", "monophage": "#d62728",
                  "cocktail_k1": "#1f77b4", "robust_k2": "#2ca02c"}
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for name, d in traj.items():
            lab = f"{name} (n={sizes[name]})"
            ax[0].plot(d["t"], np.maximum(d["total"], 1), label=lab,
                       color=colors[name], lw=2)
            ax[1].plot(d["t"], np.maximum(d["R"], 1), label=lab,
                       color=colors[name], lw=2)
        for a, ttl, yl in ((ax[0], "Total bacterial load", "CFU/mL (total)"),
                           (ax[1], "Resistant subpopulation", "CFU/mL (resistant)")):
            a.set_yscale("log"); a.set_xlabel("time (h)"); a.set_ylabel(yl)
            a.set_title(ttl); a.grid(alpha=0.3); a.legend(loc="lower right")
        fig.suptitle("Eco-evolutionary phage therapy simulation")
        fig.tight_layout()
        fig.savefig(rd / "temporal_dynamics.png", dpi=150)
        log.info("Wrote figure %s", rd / "temporal_dynamics.png")
    except Exception as e:
        log.warning("figure skipped (%s)", e)

    log.info("Wrote temporal_outcomes.csv + temporal_summary.json")


if __name__ == "__main__":
    main()
