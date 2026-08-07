"""
statistical_tests.py - Publication-grade statistical analysis
=============================================================
Implements:
    - pairwise_wilcoxon         (paired species-level AUC tests)
    - bh_fdr_correction         (no external deps)
    - delong_auc_test           (DeLong et al. 1988 method)
    - mcnemar_test              (binary prediction comparison)
    - bootstrap_ci              (standard + cluster bootstrap)

run_all_tests(all_results, dataset, cfg, seed) -> dict
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm, wilcoxon

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pairwise Wilcoxon
# ---------------------------------------------------------------------------
def pairwise_wilcoxon(all_loso_aucs: dict) -> list[dict]:
    rows = []
    names = list(all_loso_aucs.keys())
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            m1, m2 = names[i], names[j]
            s1 = all_loso_aucs[m1]
            s2 = all_loso_aucs[m2]
            shared = s1.index.intersection(s2.index)
            if len(shared) < 5:
                continue
            v1 = s1.loc[shared].astype(float).values
            v2 = s2.loc[shared].astype(float).values
            try:
                stat, p = wilcoxon(v1, v2, zero_method="wilcox")
            except Exception:
                stat, p = np.nan, np.nan
            diff = v1 - v2
            denom = np.std(diff, ddof=1) + 1e-9
            d = float(np.mean(diff) / denom)
            rows.append({
                "model_1": m1, "model_2": m2,
                "n_shared": int(len(shared)),
                "mean_auc_1": round(float(v1.mean()), 4),
                "mean_auc_2": round(float(v2.mean()), 4),
                "wilcoxon_stat": round(float(stat), 4) if np.isfinite(stat) else None,
                "p_raw": round(float(p), 6) if np.isfinite(p) else None,
                "effect_size_d": round(d, 4),
            })
    return rows


# ---------------------------------------------------------------------------
# BH-FDR correction
# ---------------------------------------------------------------------------
def bh_fdr_correction(p_values: list[float], alpha: float = 0.05):
    n = len(p_values)
    if n == 0:
        return [], []
    arr = np.array(p_values, dtype=float)
    finite_mask = np.isfinite(arr)
    order = np.argsort(arr)
    sorted_p = arr[order]
    bh_critical = (np.arange(1, n + 1) / n) * alpha
    reject_sorted = sorted_p <= bh_critical
    if reject_sorted.any():
        last = int(np.where(reject_sorted)[0].max())
        reject_sorted[:last + 1] = True
    q = np.minimum(1.0, sorted_p * n / np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    reject_out = np.zeros(n, dtype=bool)
    q_out = np.ones(n)
    reject_out[order] = reject_sorted
    q_out[order] = q
    reject_out[~finite_mask] = False
    q_out[~finite_mask] = 1.0
    return reject_out.tolist(), q_out.tolist()


# ---------------------------------------------------------------------------
# DeLong AUC test
# ---------------------------------------------------------------------------
def delong_auc_test(y_true: np.ndarray, prob1: np.ndarray,
                    prob2: np.ndarray) -> dict:
    y = np.asarray(y_true)
    p1 = np.asarray(prob1, dtype=float)
    p2 = np.asarray(prob2, dtype=float)
    pos_mask = y == 1
    neg_mask = y == 0

    def structural_components(pos_scores, neg_scores):
        m = len(pos_scores)
        n = len(neg_scores)
        if m == 0 or n == 0:
            return 0.5, 1e-6, np.zeros(max(m, 1)), np.zeros(max(n, 1))
        V10 = np.array([
            (pos_scores[i] > neg_scores).mean() +
            0.5 * (pos_scores[i] == neg_scores).mean()
            for i in range(m)])
        V01 = np.array([
            (neg_scores[j] < pos_scores).mean() +
            0.5 * (neg_scores[j] == pos_scores).mean()
            for j in range(n)])
        auc = float(V10.mean())
        var = (np.var(V10, ddof=1) / m + np.var(V01, ddof=1) / n
               if m > 1 and n > 1 else 1e-6)
        return auc, var, V10, V01

    auc1, var1, V10_1, V01_1 = structural_components(p1[pos_mask], p1[neg_mask])
    auc2, var2, V10_2, V01_2 = structural_components(p2[pos_mask], p2[neg_mask])
    m = int(pos_mask.sum())
    n = int(neg_mask.sum())
    cov = ((np.cov(V10_1, V10_2)[0, 1] / m if m > 1 and len(V10_1) == len(V10_2) else 0) +
           (np.cov(V01_1, V01_2)[0, 1] / n if n > 1 and len(V01_1) == len(V01_2) else 0))
    var_diff = max(var1 + var2 - 2 * cov, 1e-12)
    z = (auc1 - auc2) / np.sqrt(var_diff)
    p_val = float(2 * (1 - norm.cdf(abs(z))))
    ci_lo = (auc1 - auc2) - 1.96 * np.sqrt(var_diff)
    ci_hi = (auc1 - auc2) + 1.96 * np.sqrt(var_diff)
    return {"auc1": round(auc1, 4), "auc2": round(auc2, 4),
            "diff": round(auc1 - auc2, 4), "z": round(float(z), 4),
            "p_value": round(p_val, 6),
            "ci_lo": round(float(ci_lo), 4), "ci_hi": round(float(ci_hi), 4)}


# ---------------------------------------------------------------------------
# McNemar's test
# ---------------------------------------------------------------------------
def mcnemar_test(y_true: np.ndarray, pred1: np.ndarray,
                 pred2: np.ndarray) -> dict:
    y = np.asarray(y_true).astype(int)
    p1 = np.asarray(pred1).astype(int)
    p2 = np.asarray(pred2).astype(int)
    b = int(((p1 == 1) & (p2 == 0) & (y == 1)).sum() +
            ((p1 == 0) & (p2 == 1) & (y == 0)).sum())
    c = int(((p1 == 0) & (p2 == 1) & (y == 1)).sum() +
            ((p1 == 1) & (p2 == 0) & (y == 0)).sum())
    n = b + c
    if n == 0:
        return {"statistic": 0.0, "p_value": 1.0, "b": 0, "c": 0}
    stat = float((abs(b - c) - 1) ** 2 / n)
    p = float(1 - chi2.cdf(stat, df=1))
    return {"statistic": round(stat, 4), "p_value": round(p, 6), "b": b, "c": c}


# ---------------------------------------------------------------------------
# Cluster bootstrap CI
# ---------------------------------------------------------------------------
def bootstrap_ci(values: np.ndarray, cluster_labels=None,
                 n_boot: int = 10000, ci: float = 0.95,
                 seed: int = 42) -> dict:
    rng = np.random.default_rng(seed)
    vals = np.asarray(values, dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    n = len(vals)
    if n == 0:
        return {"mean": 0.0, "std_ci_lo": 0.0, "std_ci_hi": 0.0,
                "cluster_ci_lo": 0.0, "cluster_ci_hi": 0.0, "n": 0,
                "note": "no data"}
    alpha = (1 - ci) / 2

    boot_means = vals[rng.integers(0, n, (n_boot, n))].mean(axis=1)
    std_lo = float(np.percentile(boot_means, alpha * 100))
    std_hi = float(np.percentile(boot_means, (1 - alpha) * 100))

    clust_lo, clust_hi = std_lo, std_hi
    if cluster_labels is not None and len(cluster_labels) == n:
        cluster_labels = np.asarray(cluster_labels)
        clusters = np.unique(cluster_labels)
        if len(clusters) >= 2:
            c_means = []
            for _ in range(n_boot):
                chosen = rng.choice(clusters, size=len(clusters), replace=True)
                sampled = np.concatenate(
                    [vals[cluster_labels == c] for c in chosen
                     if (cluster_labels == c).any()])
                if len(sampled) > 0:
                    c_means.append(sampled.mean())
            if c_means:
                clust_lo = float(np.percentile(c_means, alpha * 100))
                clust_hi = float(np.percentile(c_means, (1 - alpha) * 100))

    return {
        "mean": round(float(vals.mean()), 4),
        "std_ci_lo": round(std_lo, 4),
        "std_ci_hi": round(std_hi, 4),
        "cluster_ci_lo": round(clust_lo, 4),
        "cluster_ci_hi": round(clust_hi, 4),
        "n": int(n),
        "note": ("Cluster bootstrap is preferred for publication - "
                 "more conservative, accounts for taxonomic non-independence"),
    }


# ---------------------------------------------------------------------------
# Pipeline driver
# ---------------------------------------------------------------------------
def _genus_of(host: str) -> str:
    return str(host).split(" ", 1)[0] if isinstance(host, str) else ""


def run_all_tests(all_results: dict, dataset: pd.DataFrame, cfg: dict,
                  seed: int) -> dict:
    """Run all statistical tests and write the supplementary tables."""
    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)
    n_boot = int(cfg["evaluation"]["bootstrap_n"])
    ci = float(cfg["evaluation"]["bootstrap_ci"])

    # Build per-species AUC series for each model
    per_species: dict[str, pd.Series] = {}
    for name, res in all_results.items():
        s = res.get("per_species_auc")
        if isinstance(s, pd.Series) and len(s):
            per_species[name] = s
            continue
        df = res.get("loso_df")
        if isinstance(df, pd.DataFrame) and "species" in df.columns and "roc_auc" in df.columns:
            per_species[name] = df.set_index("species")["roc_auc"]

    pairwise = pairwise_wilcoxon(per_species) if per_species else []

    # BH-FDR over pairwise p-values
    p_vals = [r["p_raw"] for r in pairwise if r["p_raw"] is not None]
    if p_vals:
        reject, q = bh_fdr_correction(p_vals, alpha=0.05)
        # weave back
        it = iter(zip(reject, q))
        for r in pairwise:
            if r["p_raw"] is None:
                r["reject_005"] = False
                r["q_bh"] = None
            else:
                rej, qv = next(it)
                r["reject_005"] = bool(rej)
                r["q_bh"] = round(float(qv), 6)
                r["reject_001"] = bool(qv < 0.01)
    else:
        for r in pairwise:
            r["reject_005"] = False
            r["q_bh"] = None
            r["reject_001"] = False

    # DeLong tests on pooled probabilities
    delong = []
    names = list(all_results.keys())
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            m1, m2 = names[i], names[j]
            r1 = all_results[m1]
            r2 = all_results[m2]
            p1 = r1.get("all_proba")
            p2 = r2.get("all_proba")
            y1 = r1.get("y_true")
            y2 = r2.get("y_true")
            if (p1 is None or p2 is None or y1 is None or y2 is None or
                    len(p1) == 0 or len(p2) == 0 or len(p1) != len(p2)):
                continue
            if not np.array_equal(np.asarray(y1, dtype=int),
                                  np.asarray(y2, dtype=int)):
                continue
            try:
                d = delong_auc_test(y1, p1, p2)
                d.update({"model_1": m1, "model_2": m2})
                delong.append(d)
            except Exception as e:
                log.warning(f"  DeLong failed for {m1} vs {m2}: {e}")

    # McNemar tests
    mcn = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            m1, m2 = names[i], names[j]
            r1 = all_results[m1]
            r2 = all_results[m2]
            p1 = r1.get("all_proba")
            p2 = r2.get("all_proba")
            y1 = r1.get("y_true")
            if (p1 is None or p2 is None or y1 is None or len(p1) == 0 or
                    len(p2) == 0 or len(p1) != len(p2)):
                continue
            pred1 = (np.asarray(p1) >= 0.5).astype(int)
            pred2 = (np.asarray(p2) >= 0.5).astype(int)
            try:
                m = mcnemar_test(y1, pred1, pred2)
                m.update({"model_1": m1, "model_2": m2})
                mcn.append(m)
            except Exception as e:
                log.warning(f"  McNemar failed for {m1} vs {m2}: {e}")

    # Bootstrap CI per model on per-species AUC, clustered by genus
    boot_rows = []
    species_to_genus = (dataset.dropna(subset=["host_genus"])
                        .drop_duplicates("host")
                        .set_index("host")["host_genus"].to_dict())
    for name, s in per_species.items():
        vals = s.values
        clusters = np.array([species_to_genus.get(idx, "unknown") for idx in s.index])
        out = bootstrap_ci(vals, cluster_labels=clusters, n_boot=n_boot, ci=ci, seed=seed)
        boot_rows.append({"model": name, "metric": "roc_auc", **out})

    # Save tables
    test_rows: list[dict] = []
    for r in pairwise:
        test_rows.append({
            "test_name": "Wilcoxon",
            "model_1": r["model_1"], "model_2": r["model_2"],
            "statistic": r.get("wilcoxon_stat"),
            "p_raw": r.get("p_raw"),
            "p_adjusted_bh": r.get("q_bh"),
            "q_bh": r.get("q_bh"),
            "reject_005": r.get("reject_005"),
            "reject_001": r.get("reject_001"),
            "effect_size": r.get("effect_size_d"),
            "n_samples": r.get("n_shared"),
            "notes": "paired species-level LOSO AUC",
        })
    for d in delong:
        test_rows.append({
            "test_name": "DeLong",
            "model_1": d["model_1"], "model_2": d["model_2"],
            "statistic": d.get("z"),
            "p_raw": d.get("p_value"),
            "p_adjusted_bh": None,
            "q_bh": None,
            "reject_005": (d.get("p_value", 1.0) or 1.0) < 0.05,
            "reject_001": (d.get("p_value", 1.0) or 1.0) < 0.01,
            "effect_size": d.get("diff"),
            "n_samples": None,
            "notes": "pooled ROC-AUC comparison",
        })
    for m in mcn:
        test_rows.append({
            "test_name": "McNemar",
            "model_1": m["model_1"], "model_2": m["model_2"],
            "statistic": m.get("statistic"),
            "p_raw": m.get("p_value"),
            "p_adjusted_bh": None,
            "q_bh": None,
            "reject_005": (m.get("p_value", 1.0) or 1.0) < 0.05,
            "reject_001": (m.get("p_value", 1.0) or 1.0) < 0.01,
            "effect_size": None,
            "n_samples": m.get("b", 0) + m.get("c", 0),
            "notes": f"b={m.get('b', 0)}, c={m.get('c', 0)}",
        })

    pd.DataFrame(test_rows).to_csv(results_dir / "statistical_tests.csv", index=False)
    pd.DataFrame(boot_rows).to_csv(results_dir / "confidence_intervals.csv", index=False)

    log.info(f"  Saved statistical_tests.csv ({len(test_rows)} rows)")
    log.info(f"  Saved confidence_intervals.csv ({len(boot_rows)} rows)")

    if boot_rows:
        for r in boot_rows:
            log.info(f"    {r['model']}: mean AUC={r['mean']:.4f}, "
                     f"cluster CI=[{r['cluster_ci_lo']:.4f}, {r['cluster_ci_hi']:.4f}]")

    return {
        "pairwise": pairwise,
        "delong": delong,
        "mcnemar": mcn,
        "bootstrap": boot_rows,
    }
